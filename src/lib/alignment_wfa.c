/*
 alignment_wfa.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// The score of the best global alignment by the wavefront algorithm (WFA:
// Marco-Sola et al. 2021, Bioinformatics 37:456-463), written from the
// paper, for sizing a band (see alignment_fill_banded()); and, keeping
// every wavefront, the alignment itself (alignment_wfa_align(), for
// needleman_wunsch --wfa). Its time grows with the alignment's penalty,
// so it is quick for similar sequences and gives up past a budget.
// WFA minimises penalties with no match bonus. With match bonus l, and a
// mismatch, gap open and gap extension costing x, o and e (a gap of k
// letters o + k e, as here), every alignment of a (m letters) and b (n)
// has 2M + 2X + Gl = m + n (M matches, X mismatches, Gl gap letters), so
// twice its score is l(m+n) less (2l+2x) X + 2o G + (2e+l) Gl, G the gaps
// opened (as Eizenga & Paten, bioRxiv 10.1101/2022.01.12.476087, convert
// it): the best score is (l(m+n) - P)/2 for the least penalty P with
// mismatch 2l+2x, gap open 2o and gap extension 2e+l

#include <stdlib.h>
#include <stdint.h>
#include <stdatomic.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include "alignment.h"
#include "alignment_alloc.h"
#include "alignment_simd.h"

// (not in the scalar-only build, -DALN_NO_SIMD)
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && \
    !defined(ALN_TSAN) && !defined(ALN_NO_SIMD)
#  include <immintrin.h>
#  define WFA_HAS_X86_AVX2 1
#endif

// Diagonal k = h - v of text b (h) and pattern a (v); offsets are h
#define WFA_NONE (INT32_MIN / 2)

typedef struct {
  int32_t lo, hi; // active diagonals held; none if lo > hi
  size_t cap;
  int32_t *buf; // allocation base (what free() receives)
  int32_t *off; // points into buf; off[k - lo] is diagonal k
} wavefront_t;

// The wavefronts of one score: match (M), insertion (I) and deletion (D)
typedef struct { wavefront_t m, i, d; } wf_slot_t;

// Reusable all-WFA_NONE buffer so wf_compute_interior can run vectorised
// even when one or more source wavefronts are empty. Shared by every step
// of one alignment, it also records whether any of them ran out of memory
// (the steps otherwise just report failure, as they do past the budget)
typedef struct {
  int32_t *buf;
  size_t cap;
  bool nomem;
} wf_null_t;

static bool wf_null_ensure(wf_null_t *nl, size_t need)
{
  if(need <= nl->cap) return true;
  size_t cap = nl->cap ? nl->cap : 256;
  while(cap < need) cap *= 2;
  int32_t *b = aln_realloc(nl->buf, cap * sizeof(int32_t));
  if(b == NULL) return false;
  for(size_t i = nl->cap; i < cap; i++) b[i] = WFA_NONE;
  nl->buf = b;
  nl->cap = cap;
  return true;
}

static inline int32_t wf_get(const wavefront_t *w, int32_t k)
{
  return k < w->lo || k > w->hi ? WFA_NONE : w->off[k - w->lo];
}

static bool wf_reserve(wavefront_t *w, int32_t lo, int32_t hi)
{
  w->lo = lo;
  w->hi = hi;
  size_t need = (size_t)(hi - lo + 1);
  if(need > w->cap) {
    int32_t *buf = aln_realloc(w->buf, need * sizeof(int32_t));
    if(buf == NULL) return false;
    w->buf = buf;
    w->cap = need;
  }
  w->off = w->buf;
  return true;
}

static inline void wf_trim(wavefront_t *w)
{
  int32_t lo = w->lo, hi = w->hi;
  if(lo > hi) return;
  int32_t *off = w->off;
  int32_t len = hi - lo + 1;
  int32_t l = 0, r = len - 1;
  while(l <= r && off[l] == WFA_NONE) l++;
  if(l > r) {
    w->lo = 1;
    w->hi = 0;
    return;
  }
  while(r > l && off[r] == WFA_NONE) r--;
  w->lo = lo + l;
  w->hi = lo + r;
  w->off = off + l;
}

static long long wfa_gcd(long long a, long long b)
{
  if(a < 0) a = -a;
  if(b < 0) b = -b;
  while(b != 0) { long long t = a % b; a = b; b = t; }
  return a;
}

// Matching letters from h, v on: their number. Letters are class numbers
// with 8 different padding bytes after each sequence, so words of 8 can be
// compared without bounds checks
static inline int32_t wf_extend(const uint8_t *a, const uint8_t *b,
                                int32_t v, int32_t h)
{
  int32_t start = h;
  for(;;) {
    uint64_t x, y;
    memcpy(&x, a + v, 8);
    memcpy(&y, b + h, 8);
    uint64_t d = x ^ y;
    if(d != 0) {
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
      int32_t same = __builtin_ctzll(d) >> 3;
#else
      int32_t same = 0;
      while(a[v + same] == b[h + same]) same++;
#endif
      return h + same - start;
    }
    v += 8;
    h += 8;
  }
}

// A pair ready for the wavefronts: letters as class numbers (a, b, with
// padding), the match score l, and the penalties of a mismatch (px), a gap
// open (po) and a gap letter (pe), doubled as above
typedef struct {
  uint8_t *a, *b;
  int32_t m, nb;
  long long l, px, po, pe, scale;
} wfa_pair_t;

// Set pair up for seq_a and seq_b, or return false where the scoring is not
// one match and one mismatch score with gap penalties, or out of memory
// (*nomem set then)
static bool wfa_setup(const sub_table_t *table, const scoring_t *scoring,
                      const char *seq_a, size_t len_a,
                      const char *seq_b, size_t len_b, wfa_pair_t *pair,
                      bool *nomem)
{
  const size_t n = table->n;
  size_t i, j;
  if(len_a == 0 || len_b == 0 || len_a >= INT32_MAX / 2 ||
     len_b >= INT32_MAX / 2 ||
     scoring->no_start_gap_penalty || scoring->no_end_gap_penalty ||
     scoring->no_gaps_in_a || scoring->no_gaps_in_b ||
     scoring->no_gaps_in_longest || scoring->no_mismatches) return false;

  // Letters become class numbers under 0xFE, one per letter (case aside if
  // !case_sensitive). The pairs the fill reads (a letter of a against one
  // of b) must score one match score within a class and one mismatch score
  // across classes: no wildcards or other scores
  bool in_a[256] = {false}, in_b[256] = {false};
  uint8_t cls[256], fold_cls[256];
  size_t nclasses = 0;
  memset(fold_cls, 0xFF, sizeof(fold_cls));
  for(i = 0; i < len_a; i++) in_a[(unsigned char)seq_a[i]] = true;
  for(i = 0; i < len_b; i++) in_b[(unsigned char)seq_b[i]] = true;
  for(i = 0; i < 256; i++) {
    if(!in_a[i] && !in_b[i]) continue;
    const size_t f = scoring->case_sensitive ? i : (size_t)tolower((int)i);
    if(fold_cls[f] == 0xFF) {
      if(nclasses >= 0xFE) return false;
      fold_cls[f] = (uint8_t)nclasses++;
    }
    cls[table->code[i]] = fold_cls[f];
  }
  long long l = 0, mis = 0;
  bool have_l = false, have_mis = false;
  for(i = 0; i < 256; i++) {
    if(!in_a[i]) continue;
    for(j = 0; j < 256; j++) {
      if(!in_b[j]) continue;
      const size_t ci = table->code[i], cj = table->code[j];
      const long long sc = table->scores[cj*n + ci];
      if(cls[ci] == cls[cj]) {
        if(have_l && sc != l) return false;
        l = sc; have_l = true;
      }
      else {
        if(have_mis && sc != mis) return false;
        mis = sc; have_mis = true;
      }
    }
  }
  if(!have_l) return false; // no letter in both

  const long long x = have_mis ? l - mis : 1; // mismatch less than a match
  const long long o = -(long long)scoring->gap_open;
  const long long e = -(long long)scoring->gap_extend;
  pair->l = l;
  pair->px = 2 * x;
  pair->po = 2 * o;
  pair->pe = 2 * e + l;
  if(pair->px <= 0 || pair->po < 0 || pair->pe <= 0 || pair->px > 1000000 ||
     pair->po + pair->pe > 1000000) return false;
  const long long g = wfa_gcd(pair->px, wfa_gcd(pair->po, pair->pe));
  pair->scale = g > 0 ? g : 1;
  pair->px /= pair->scale;
  pair->po /= pair->scale;
  pair->pe /= pair->scale;

  // Class numbers of the letters, then 32 padding bytes each (different),
  // so both 8-byte scalar and 32-byte AVX2 extensions need no bounds checks
  pair->m = (int32_t)len_a;
  pair->nb = (int32_t)len_b;
  pair->a = aln_malloc(len_a + len_b + 64);
  if(pair->a == NULL) { *nomem = true; return false; }
  pair->b = pair->a + len_a + 32;
  for(i = 0; i < len_a; i++)
    pair->a[i] = cls[table->code[(unsigned char)seq_a[i]]];
  for(i = 0; i < len_b; i++)
    pair->b[i] = cls[table->code[(unsigned char)seq_b[i]]];
  memset(pair->a + len_a, 0xFE, 32);
  memset(pair->b + len_b, 0xFF, 32);
  return true;
}

// The wavefronts of a pair, score by score, until one reaches the last
// cell. With keep, every score's wavefronts are kept (slots[s], for the
// traceback); otherwise only the last px or po + pe scores', in a ring.
// Returns the least penalty, or -1 past budget diagonal updates (cells
// kept, with keep) or if out of memory (*nomem set then). *slots and
// *nslots hold what to free (wfa_free())
static long long wfa_fill(const wfa_pair_t *pair, size_t budget, bool keep,
                          wf_slot_t **slots_out, size_t *nslots_out,
                          bool *nomem)
{
  const int32_t m = pair->m, nb = pair->nb, kend = nb - m;
  const long long px = pair->px, po = pair->po, pe = pair->pe;
  const size_t ring = (size_t)(px > po + pe ? px : po + pe) + 1;
  size_t nslots = keep ? 0 : ring, i;
  wf_slot_t *slots = keep ? NULL : aln_calloc(ring, sizeof(wf_slot_t));
  const wavefront_t none = {1, 0, 0, NULL, NULL};
  long long s, result = -1;
  size_t work = 0, empty_streak = 0;
  if(!keep && slots == NULL) {
    *slots_out = NULL;
    *nslots_out = 0;
    *nomem = true;
    return -1;
  }
  *slots_out = slots;
  *nslots_out = nslots;
  for(i = 0; i < nslots; i++) {
    slots[i].m = slots[i].i = slots[i].d = none;
  }
  #define WF_SLOT(t) (keep ? &slots[(t)] : &slots[(size_t)(t) % ring])

  for(s = 0; ; s++) {
    if(keep) {
      // One more slot, kept (doubling the room as needed)
      if((size_t)s >= nslots) {
        size_t grow = nslots ? 2 * nslots : 64;
        wf_slot_t *g = aln_realloc(slots, grow * sizeof(wf_slot_t));
        if(g == NULL) { *nomem = true; break; }
        slots = g;
        for(i = nslots; i < grow; i++)
          slots[i].m = slots[i].i = slots[i].d = none;
        nslots = grow;
        *slots_out = slots;
        *nslots_out = nslots;
      }
    }
    wf_slot_t *cur = WF_SLOT(s);
    const wavefront_t *mx = s >= px ? &WF_SLOT(s - px)->m : &none;
    const wavefront_t *mo = s >= po + pe ? &WF_SLOT(s - po - pe)->m : &none;
    const wavefront_t *ie = s >= pe ? &WF_SLOT(s - pe)->i : &none;
    const wavefront_t *de = s >= pe ? &WF_SLOT(s - pe)->d : &none;
    int32_t lo, hi, k;

    if(s == 0) {
      if(!wf_reserve(&cur->m, 0, 0)) { *nomem = true; break; }
      cur->m.off[0] = wf_extend(pair->a, pair->b, 0, 0);
      cur->i.lo = cur->d.lo = 1;
      cur->i.hi = cur->d.hi = 0;
    }
    else {
      // The diagonals some source reaches
      lo = INT32_MAX; hi = INT32_MIN;
      if(mx->lo <= mx->hi) { lo = mx->lo; hi = mx->hi; }
      if(mo->lo <= mo->hi) {
        lo = mo->lo - 1 < lo ? mo->lo - 1 : lo;
        hi = mo->hi + 1 > hi ? mo->hi + 1 : hi;
      }
      if(ie->lo <= ie->hi) {
        lo = ie->lo + 1 < lo ? ie->lo + 1 : lo;
        hi = ie->hi + 1 > hi ? ie->hi + 1 : hi;
      }
      if(de->lo <= de->hi) {
        lo = de->lo - 1 < lo ? de->lo - 1 : lo;
        hi = de->hi - 1 > hi ? de->hi - 1 : hi;
      }
      if(lo > hi) {
        cur->m.lo = cur->i.lo = cur->d.lo = 1;
        cur->m.hi = cur->i.hi = cur->d.hi = 0;
        if(++empty_streak >= ring) break;
        continue;
      }
      empty_streak = 0;
      if(lo < -m) lo = -m;
      if(hi > nb) hi = nb;
      if(lo > hi) {
        cur->m.lo = cur->i.lo = cur->d.lo = 1;
        cur->m.hi = cur->i.hi = cur->d.hi = 0;
        if(++empty_streak >= ring) break;
        continue;
      }
      empty_streak = 0;
      work += (size_t)(hi - lo + 1);
      if(work > budget) break;
      if(!wf_reserve(&cur->m, lo, hi) || !wf_reserve(&cur->i, lo, hi) ||
         !wf_reserve(&cur->d, lo, hi)) { *nomem = true; break; }
      for(k = lo; k <= hi; k++) {
        // Insertion (a letter of b): from k-1, one along b
        int32_t vi = wf_get(mo, k - 1), ei = wf_get(ie, k - 1);
        int32_t ins = (vi > ei ? vi : ei) + 1;
        if(ins > nb || ins - k > m || ins < 0) ins = WFA_NONE;
        // Deletion (a letter of a): from k+1, at the same h
        int32_t vd = wf_get(mo, k + 1), ed = wf_get(de, k + 1);
        int32_t del = vd > ed ? vd : ed;
        if(del - k > m || del < 0) del = WFA_NONE;
        // Mismatch: along k
        int32_t sub = wf_get(mx, k) + 1;
        if(sub > nb || sub - k > m || sub < 0) sub = WFA_NONE;
        int32_t best = sub > ins ? sub : ins;
        best = best > del ? best : del;
        cur->i.off[k - lo] = ins;
        cur->d.off[k - lo] = del;
        if(best >= 0 && best - k >= 0) {
          int32_t run = wf_extend(pair->a, pair->b, best - k, best);
          // (kept wavefronts count as cells, not as letters compared)
          if(!keep) work += (size_t)run / 8;
          best += run;
        }
        else best = WFA_NONE;
        cur->m.off[k - lo] = best;
      }
      wf_trim(&cur->m);
      wf_trim(&cur->i);
      wf_trim(&cur->d);
    }
    if(wf_get(&cur->m, kend) >= nb) { result = s; break; }
    if(s > (long long)INT32_MAX / 2) break;
  }
  #undef WF_SLOT
  return result;
}

static void wfa_free(wf_slot_t *slots, size_t nslots)
{
  if(slots == NULL) return;
  size_t i;
  for(i = 0; i < nslots; i++) {
    free(slots[i].m.buf);
    free(slots[i].i.buf);
    free(slots[i].d.buf);
  }
  free(slots);
}

// The best score, (l(m+n) - P)/2, if it fits a score_t
static bool wfa_best(const wfa_pair_t *pair, long long penalty,
                     score_t *score)
{
  const long long best = (pair->l * ((long long)pair->m + pair->nb) -
                          penalty * pair->scale) / 2;
  if(best < INT_MIN || best > INT_MAX) return false;
  *score = (score_t)best;
  return true;
}

bool alignment_wfa_score(const sub_table_t *table, const scoring_t *scoring,
                         const char *seq_a, size_t len_a,
                         const char *seq_b, size_t len_b,
                         size_t budget, score_t *score)
{
  wfa_pair_t pair;
  bool nomem = false; // scored or not either way: the caller falls back
  if(!wfa_setup(table, scoring, seq_a, len_a, seq_b, len_b, &pair, &nomem))
    return false;
  wf_slot_t *slots;
  size_t nslots;
  const long long penalty = wfa_fill(&pair, budget, false, &slots, &nslots,
                                     &nomem);
  wfa_free(slots, nslots);
  free(pair.a);
  return penalty >= 0 && wfa_best(&pair, penalty, score);
}

// The alignment from the kept wavefronts, back from the last cell: out_a
// and out_b get it reversed, and the number of columns is returned (0 if
// the wavefronts are not consistent, which would be a bug). At each step
// the predecessor that gives the offset is taken, in the order insertion,
// deletion, mismatch from a match cell and extension before opening in a
// gap: one of the best alignments, though not always the one the
// dynamic programming traceback picks among equal-scoring ones
static size_t wfa_trace(const wfa_pair_t *pair, const wf_slot_t *slots,
                        long long s, const char *seq_a, const char *seq_b,
                        char *out_a, char *out_b)
{
  const int32_t m = pair->m, nb = pair->nb;
  const long long px = pair->px, po = pair->po, pe = pair->pe;
  enum { WM, WI, WD } state = WM;
  int32_t k = nb - m, h = nb;
  size_t n = 0;

  for(;;) {
    const wf_slot_t *cur = &slots[s];
    if(state == WM) {
      // The offset before this score's run of matches
      int32_t src = 0, ins = WFA_NONE, del = WFA_NONE, sub = WFA_NONE;
      if(s > 0) {
        ins = wf_get(&cur->i, k);
        del = wf_get(&cur->d, k);
        if(s >= px) {
          sub = wf_get(&slots[s - px].m, k) + 1;
          if(sub > nb || sub - k > m || sub < 0) sub = WFA_NONE;
        }
        src = sub > ins ? sub : ins;
        src = src > del ? src : del;
        if(src < 0) return 0;
      }
      for(; h > src; h--, n++) {
        out_a[n] = seq_a[h - k - 1];
        out_b[n] = seq_b[h - 1];
      }
      if(s == 0) return h == 0 && k == 0 ? n : 0;
      if(src == ins) state = WI;
      else if(src == del) state = WD;
      else {
        // A mismatch
        out_a[n] = seq_a[h - k - 1];
        out_b[n] = seq_b[h - 1];
        n++;
        h--;
        s -= px;
      }
    }
    else if(state == WI) {
      // A letter of b against a gap, from k-1: extended or opened
      out_a[n] = '-';
      out_b[n] = seq_b[h - 1];
      n++;
      h--;
      k--;
      if(s >= pe && wf_get(&slots[s - pe].i, k) == h) s -= pe;
      else if(s >= po + pe && wf_get(&slots[s - po - pe].m, k) == h) {
        s -= po + pe;
        state = WM;
      }
      else return 0;
    }
    else {
      // A letter of a against a gap, from k+1 at the same h
      out_a[n] = seq_a[h - k - 1];
      out_b[n] = '-';
      n++;
      k++;
      if(s >= pe && wf_get(&slots[s - pe].d, k) == h) s -= pe;
      else if(s >= po + pe && wf_get(&slots[s - po - pe].m, k) == h) {
        s -= po + pe;
        state = WM;
      }
      else return 0;
    }
  }
}

// Bounded 8-byte match extension for subproblems without trailing sentinels
static inline int32_t wf_extend_lim(const uint8_t *a, const uint8_t *b,
                                    int32_t v, int32_t h,
                                    int32_t m, int32_t nb)
{
  int32_t start = h;
  while(v + 8 <= m && h + 8 <= nb) {
    uint64_t x, y;
    memcpy(&x, a + v, 8);
    memcpy(&y, b + h, 8);
    uint64_t d = x ^ y;
    if(d != 0) {
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
      int32_t same = __builtin_ctzll(d) >> 3;
#else
      int32_t same = 0;
      while(a[v + same] == b[h + same]) same++;
#endif
      return h + same - start;
    }
    v += 8;
    h += 8;
  }
  while(v < m && h < nb && a[v] == b[h]) { v++; h++; }
  return h - start;
}

#ifdef WFA_HAS_X86_AVX2
static inline bool wfa_cpu_avx2(void)
{
  static _Atomic int cached = -1;
  int v = atomic_load_explicit(&cached, memory_order_relaxed);
  if(v < 0) {
    v = __builtin_cpu_supports("avx2") ? 1 : 0;
    atomic_store_explicit(&cached, v, memory_order_relaxed);
  }
  return v != 0;
}

__attribute__((target("avx2")))
static int32_t wf_extend_loop_avx2(const uint8_t *a, const uint8_t *b,
                                   int32_t lo, int32_t hi, int32_t *moff,
                                   int32_t mh)
{
  for(int32_t k = lo; k <= hi; k++) {
    int32_t best = moff[k - lo];
    if(best >= 0) {
      int32_t v = best - k, h = best;
      for(;;) {
        __m256i va = _mm256_loadu_si256((const __m256i *)(a + v));
        __m256i vb = _mm256_loadu_si256((const __m256i *)(b + h));
        uint32_t mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(va, vb));
        if(mask != 0xFFFFFFFFu) {
          h += __builtin_ctz(~mask);
          break;
        }
        v += 32;
        h += 32;
      }
      moff[k - lo] = h;
      if(h > mh) mh = h;
    }
  }
  return mh;
}

__attribute__((target("avx2")))
static int32_t wf_extend_lim_loop_avx2(const uint8_t *a, const uint8_t *b,
                                       int32_t m, int32_t nb,
                                       int32_t lo, int32_t hi, int32_t *moff,
                                       int32_t mh)
{
  for(int32_t k = lo; k <= hi; k++) {
    int32_t best = moff[k - lo];
    if(best >= 0) {
      int32_t v = best - k, h = best;
      while(v + 32 <= m && h + 32 <= nb) {
        __m256i va = _mm256_loadu_si256((const __m256i *)(a + v));
        __m256i vb = _mm256_loadu_si256((const __m256i *)(b + h));
        uint32_t mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(va, vb));
        if(mask != 0xFFFFFFFFu) {
          h += __builtin_ctz(~mask);
          goto done_k;
        }
        v += 32;
        h += 32;
      }
      h += wf_extend_lim(a, b, v, h, m, nb);
    done_k:
      moff[k - lo] = h;
      if(h > mh) mh = h;
    }
  }
  return mh;
}
#endif

typedef enum { WFA_M = 0, WFA_I = 1, WFA_D = 2 } wfa_state_t;

// One direction of a pair for wf_step(): the class strings a and b (padded:
// sentinel bytes after each, so extensions need no bounds), their lengths,
// the penalties, the state the alignment starts in (s0), and the store of
// wavefronts: every score kept (keep) in slots, or a ring of `ring` scores
typedef struct {
  const uint8_t *a, *b;
  int32_t m, nb;
  long long px, po, pe;
  bool padded;
  wfa_state_t s0;
  wf_slot_t *slots;
  size_t ring;
  bool keep;
  wf_null_t *null_wf;
} wf_dir_t;

// What every sub-problem of the bidirectional alignment shares: the
// penalties, the null wavefront buffer (and the out-of-memory flag) and the
// diagonal updates left in the budget
typedef struct {
  long long px, po, pe;
  wf_null_t *null_wf;
  size_t *budget_left;
} biwfa_ctx_t;

// A sub-problem: a[0..m) against b[0..nb) as class strings forward (a, b)
// and reversed (ra, rb: the reversed whole strings' suffixes) and as letters
// for the output (sa, sb), entered in state s0 and left in state s1.
// fwd_pad / rev_pad: the forward / reversed strings have their sentinel
// bytes after them
typedef struct {
  const uint8_t *a, *ra;
  const char *sa;
  int32_t m;
  const uint8_t *b, *rb;
  const char *sb;
  int32_t nb;
  wfa_state_t s0, s1;
  bool fwd_pad, rev_pad;
} biwfa_sub_t;

// Where the forward and reverse wavefronts met: the least penalty and the
// cell (v, h) and state of the meeting
typedef struct {
  long long s;
  int32_t v, h;
  wfa_state_t state;
} biwfa_meet_t;

// Vectorisable hot interior compute loop (multi-versioned for AVX2/SSE4.1/SSE2).
// All pointers are strictly in-bounds for [0 .. len-1], and each candidate
// (ins, del, sub) is range-checked before taking the maximum, matching the
// scalar edge loops identically.
ALN_TARGET_CLONES
static void wf_compute_interior(int32_t k0, int32_t len, int32_t m, int32_t nb,
                                const int32_t *restrict mxp,
                                const int32_t *restrict molp,
                                const int32_t *restrict morp,
                                const int32_t *restrict iep,
                                const int32_t *restrict dep,
                                int32_t *restrict moff,
                                int32_t *restrict ioff,
                                int32_t *restrict doff)
{
  int32_t idx;
  for(idx = 0; idx < len; idx++) {
    int32_t k = k0 + idx;
    int32_t vi = molp[idx], ei = iep[idx];
    int32_t ins = (vi > ei ? vi : ei) + 1;
    ins = (ins <= nb && ins - k <= m && ins >= 0) ? ins : WFA_NONE;
    int32_t vd = morp[idx], ed = dep[idx];
    int32_t del = vd > ed ? vd : ed;
    del = (del - k <= m && del >= 0) ? del : WFA_NONE;
    int32_t sub = mxp[idx] + 1;
    sub = (sub <= nb && sub - k <= m && sub >= 0) ? sub : WFA_NONE;
    int32_t best = sub > ins ? sub : ins;
    best = best > del ? best : del;
    ioff[idx] = ins;
    doff[idx] = del;
    moff[idx] = (best >= 0 && best - k >= 0) ? best : WFA_NONE;
  }
}

// Advance one wavefront score s in ring or kept storage. Returns diagonal
// updates performed, or SIZE_MAX on OOM. Writes max offset reached in *max_h.
// When padded is true, a[m..m+7] and b[nb..nb+7] hold distinct sentinel bytes
// so unchecked 8-byte wf_extend() is safe.
static size_t wf_step(const wf_dir_t *d, long long s, int32_t *max_h)
{
  const uint8_t *const a = d->a, *const b = d->b;
  const int32_t m = d->m, nb = d->nb;
  const long long px = d->px, po = d->po, pe = d->pe;
  const wfa_state_t s0 = d->s0;
  wf_slot_t *const slots = d->slots;
  const size_t ring = d->ring;
  const bool keep = d->keep, padded = d->padded;
  wf_null_t *const null_wf = d->null_wf;
  const wavefront_t none = {1, 0, 0, NULL, NULL};
  #define SLOT(t) (keep ? &slots[(t)] : &slots[(size_t)(t) % ring])
  wf_slot_t *cur = SLOT(s);
  if(s == 0) {
    if(!wf_reserve(&cur->m, 0, 0)) return SIZE_MAX;
    int32_t h0 = padded ? wf_extend(a, b, 0, 0) : wf_extend_lim(a, b, 0, 0, m, nb);
    cur->m.off[0] = h0;
    if(max_h && h0 > *max_h) *max_h = h0;
    if(s0 == WFA_I) {
      if(!wf_reserve(&cur->i, 0, 0)) return SIZE_MAX;
      cur->i.off[0] = 0;
      cur->d.lo = 1; cur->d.hi = 0;
    }
    else if(s0 == WFA_D) {
      if(!wf_reserve(&cur->d, 0, 0)) return SIZE_MAX;
      cur->d.off[0] = 0;
      cur->i.lo = 1; cur->i.hi = 0;
    }
    else {
      cur->i.lo = cur->d.lo = 1;
      cur->i.hi = cur->d.hi = 0;
    }
    return 1;
  }
  const wavefront_t *mx = s >= px ? &SLOT(s - px)->m : &none;
  const wavefront_t *mo = s >= po + pe ? &SLOT(s - po - pe)->m : &none;
  const wavefront_t *ie = s >= pe ? &SLOT(s - pe)->i : &none;
  const wavefront_t *de = s >= pe ? &SLOT(s - pe)->d : &none;
  int32_t lo = INT32_MAX, hi = INT32_MIN, k;
  if(mx->lo <= mx->hi) { lo = mx->lo; hi = mx->hi; }
  if(mo->lo <= mo->hi) {
    lo = mo->lo - 1 < lo ? mo->lo - 1 : lo;
    hi = mo->hi + 1 > hi ? mo->hi + 1 : hi;
  }
  if(ie->lo <= ie->hi) {
    lo = ie->lo + 1 < lo ? ie->lo + 1 : lo;
    hi = ie->hi + 1 > hi ? ie->hi + 1 : hi;
  }
  if(de->lo <= de->hi) {
    lo = de->lo - 1 < lo ? de->lo - 1 : lo;
    hi = de->hi - 1 > hi ? de->hi - 1 : hi;
  }
  if(lo < -m) lo = -m;
  if(hi > nb) hi = nb;
  if(lo > hi) {
    cur->m.lo = cur->i.lo = cur->d.lo = 1;
    cur->m.hi = cur->i.hi = cur->d.hi = 0;
    return 0;
  }
  if(!wf_reserve(&cur->m, lo, hi) || !wf_reserve(&cur->i, lo, hi) ||
     !wf_reserve(&cur->d, lo, hi)) return SIZE_MAX;
  int32_t *moff = cur->m.off;
  int32_t *ioff = cur->i.off;
  int32_t *doff = cur->d.off;
  const bool has_mx = mx->lo <= mx->hi;
  const bool has_mo = mo->lo <= mo->hi;
  const bool has_ie = ie->lo <= ie->hi;
  const bool has_de = de->lo <= de->hi;
  // Interior range where every non-empty source wavefront is in-bounds
  // (empty sources read from null_wf->buf filled with WFA_NONE)
  int32_t in_lo = lo, in_hi = hi;
  if(has_mx) {
    if(mx->lo > in_lo) in_lo = mx->lo;
    if(mx->hi < in_hi) in_hi = mx->hi;
  }
  if(has_mo) {
    if(mo->lo + 1 > in_lo) in_lo = mo->lo + 1;
    if(mo->hi - 1 < in_hi) in_hi = mo->hi - 1;
  }
  if(has_ie) {
    if(ie->lo + 1 > in_lo) in_lo = ie->lo + 1;
    if(ie->hi + 1 < in_hi) in_hi = ie->hi + 1;
  }
  if(has_de) {
    if(de->lo - 1 > in_lo) in_lo = de->lo - 1;
    if(de->hi - 1 < in_hi) in_hi = de->hi - 1;
  }
  if(in_lo <= in_hi && (!has_mx || !has_mo || !has_ie || !has_de)) {
    if(!wf_null_ensure(null_wf, (size_t)(in_hi - in_lo + 1))) return SIZE_MAX;
  }

  // Loop 1a: boundary diagonals before in_lo
  const int32_t b_end = in_lo <= in_hi ? in_lo - 1 : hi;
  for(k = lo; k <= b_end; k++) {
    int32_t vi = wf_get(mo, k - 1), ei = wf_get(ie, k - 1);
    int32_t ins = (vi > ei ? vi : ei) + 1;
    if(ins > nb || ins - k > m || ins < 0) ins = WFA_NONE;
    int32_t vd = wf_get(mo, k + 1), ed = wf_get(de, k + 1);
    int32_t del = vd > ed ? vd : ed;
    if(del - k > m || del < 0) del = WFA_NONE;
    int32_t sub = wf_get(mx, k) + 1;
    if(sub > nb || sub - k > m || sub < 0) sub = WFA_NONE;
    int32_t best = sub > ins ? sub : ins;
    best = best > del ? best : del;
    ioff[k - lo] = ins;
    doff[k - lo] = del;
    moff[k - lo] = (best >= 0 && best - k >= 0) ? best : WFA_NONE;
  }
  // Loop 1b: hot interior (zero bounds checks, AVX2/SSE4.1/SSE2 multi-versioned)
  if(in_lo <= in_hi) {
    const int32_t in_len = in_hi - in_lo + 1;
    const int32_t off_d = in_lo - lo;
    const int32_t *mxp  = has_mx ? mx->off + (in_lo - mx->lo)     : null_wf->buf;
    const int32_t *molp = has_mo ? mo->off + (in_lo - 1 - mo->lo) : null_wf->buf;
    const int32_t *morp = has_mo ? mo->off + (in_lo + 1 - mo->lo) : null_wf->buf;
    const int32_t *iep  = has_ie ? ie->off + (in_lo - 1 - ie->lo) : null_wf->buf;
    const int32_t *dep  = has_de ? de->off + (in_lo + 1 - de->lo) : null_wf->buf;
    wf_compute_interior(in_lo, in_len, m, nb, mxp, molp, morp, iep, dep,
                        moff + off_d, ioff + off_d, doff + off_d);
    // Loop 1c: boundary diagonals after in_hi
    for(k = in_hi + 1; k <= hi; k++) {
      int32_t vi = wf_get(mo, k - 1), ei = wf_get(ie, k - 1);
      int32_t ins = (vi > ei ? vi : ei) + 1;
      if(ins > nb || ins - k > m || ins < 0) ins = WFA_NONE;
      int32_t vd = wf_get(mo, k + 1), ed = wf_get(de, k + 1);
      int32_t del = vd > ed ? vd : ed;
      if(del - k > m || del < 0) del = WFA_NONE;
      int32_t sub = wf_get(mx, k) + 1;
      if(sub > nb || sub - k > m || sub < 0) sub = WFA_NONE;
      int32_t best = sub > ins ? sub : ins;
      best = best > del ? best : del;
      ioff[k - lo] = ins;
      doff[k - lo] = del;
      moff[k - lo] = (best >= 0 && best - k >= 0) ? best : WFA_NONE;
    }
  }
  // Loop 2: match extension (32-byte AVX2 where supported, else 8-byte scalar)
  int32_t mh = max_h ? *max_h : 0;
#ifdef WFA_HAS_X86_AVX2
  if(wfa_cpu_avx2()) {
    mh = padded ? wf_extend_loop_avx2(a, b, lo, hi, moff, mh)
                : wf_extend_lim_loop_avx2(a, b, m, nb, lo, hi, moff, mh);
  } else
#endif
  if(padded) {
    for(k = lo; k <= hi; k++) {
      int32_t best = moff[k - lo];
      if(best >= 0) {
        best += wf_extend(a, b, best - k, best);
        moff[k - lo] = best;
        if(best > mh) mh = best;
      }
    }
  }
  else {
    for(k = lo; k <= hi; k++) {
      int32_t best = moff[k - lo];
      if(best >= 0) {
        best += wf_extend_lim(a, b, best - k, best, m, nb);
        moff[k - lo] = best;
        if(best > mh) mh = best;
      }
    }
  }
  if(max_h) *max_h = mh;
  wf_trim(&cur->m);
  wf_trim(&cur->i);
  wf_trim(&cur->d);
  #undef SLOT
  return (size_t)(hi - lo + 1);
}

// Base-case kept-wavefront fill & traceback with general start/end states.
// Writes forward-order alignment into out_a/out_b and returns length (0 on fail)
static size_t biwfa_base(const biwfa_ctx_t *c, const biwfa_sub_t *q,
                         char *out_a, char *out_b)
{
  const uint8_t *const a = q->a, *const b = q->b;
  const char *const sa = q->sa, *const sb = q->sb;
  const int32_t m = q->m, nb = q->nb;
  const long long px = c->px, po = c->po, pe = c->pe;
  const wfa_state_t s0 = q->s0, s1 = q->s1;
  wf_null_t *const null_wf = c->null_wf;
  size_t *const budget_left = c->budget_left;
  const int32_t kend = nb - m;
  size_t nslots = 0, i;
  wf_slot_t *slots = NULL;
  const wavefront_t none = {1, 0, 0, NULL, NULL};
  long long s, penalty = -1;

  for(s = 0; ; s++) {
    if((size_t)s >= nslots) {
      size_t grow = nslots ? 2 * nslots : 64;
      wf_slot_t *g = aln_realloc(slots, grow * sizeof(wf_slot_t));
      if(g == NULL) { null_wf->nomem = true; break; }
      slots = g;
      for(i = nslots; i < grow; i++)
        slots[i].m = slots[i].i = slots[i].d = none;
      nslots = grow;
    }
    const wf_dir_t dir = {a, b, m, nb, px, po, pe, q->fwd_pad, s0,
                          slots, 0, true, null_wf};
    size_t w = wf_step(&dir, s, NULL);
    if(w == SIZE_MAX) { null_wf->nomem = true; break; }
    if(w > *budget_left) break;
    *budget_left -= w;
    const wf_slot_t *cur = &slots[s];
    int32_t end_h = s1 == WFA_M ? wf_get(&cur->m, kend) :
                    s1 == WFA_I ? wf_get(&cur->i, kend) :
                                  wf_get(&cur->d, kend);
    if(end_h >= nb) { penalty = s; break; }
    if(s > (long long)INT32_MAX / 2) break;
  }
  if(penalty < 0) { wfa_free(slots, nslots); return 0; }

  wfa_state_t state = s1;
  int32_t k = kend, h = nb;
  size_t n = 0;
  s = penalty;
  for(;;) {
    const wf_slot_t *cur = &slots[s];
    if(state == WFA_M) {
      int32_t src = 0, ins = WFA_NONE, del = WFA_NONE, sub = WFA_NONE;
      if(s > 0) {
        ins = wf_get(&cur->i, k);
        del = wf_get(&cur->d, k);
        if(s >= px) {
          sub = wf_get(&slots[s - px].m, k) + 1;
          if(sub > nb || sub - k > m || sub < 0) sub = WFA_NONE;
        }
        src = sub > ins ? sub : ins;
        src = src > del ? src : del;
        if(src < 0) { n = 0; break; }
      }
      for(; h > src; h--, n++) {
        out_a[n] = sa[h - k - 1];
        out_b[n] = sb[h - 1];
      }
      if(s == 0) { if(h != 0 || k != 0) n = 0; break; }
      if(src == ins) state = WFA_I;
      else if(src == del) state = WFA_D;
      else {
        out_a[n] = sa[h - k - 1];
        out_b[n] = sb[h - 1];
        n++; h--; s -= px;
      }
    }
    else if(state == WFA_I) {
      if(s == 0 && s0 == WFA_I && h == 0 && k == 0) break;
      out_a[n] = '-';
      out_b[n] = sb[h - 1];
      n++; h--; k--;
      if(s >= pe && wf_get(&slots[s - pe].i, k) == h) s -= pe;
      else if(s >= po + pe && wf_get(&slots[s - po - pe].m, k) == h) {
        s -= po + pe; state = WFA_M;
      }
      else { n = 0; break; }
    }
    else {
      if(s == 0 && s0 == WFA_D && h == 0 && k == 0) break;
      out_a[n] = sa[h - k - 1];
      out_b[n] = '-';
      n++; k++;
      if(s >= pe && wf_get(&slots[s - pe].d, k) == h) s -= pe;
      else if(s >= po + pe && wf_get(&slots[s - po - pe].m, k) == h) {
        s -= po + pe; state = WFA_M;
      }
      else { n = 0; break; }
    }
  }
  wfa_free(slots, nslots);
  for(i = 0; i < n / 2; i++) {
    char t = out_a[i]; out_a[i] = out_a[n - 1 - i]; out_a[n - 1 - i] = t;
    t = out_b[i]; out_b[i] = out_b[n - 1 - i]; out_b[n - 1 - i] = t;
  }
  return n;
}

// Check overlap between forward slot at sf and reverse slot at sr
static void biwfa_check(const wf_slot_t *fs, long long sf,
                        const wf_slot_t *rs, long long sr,
                        int32_t m, int32_t nb, long long po,
                        biwfa_meet_t *best)
{
  const int32_t kend = nb - m;
  int32_t lo, hi, k;
  if(sf + sr < best->s && fs->m.lo <= fs->m.hi && rs->m.lo <= rs->m.hi) {
    lo = fs->m.lo > kend - rs->m.hi ? fs->m.lo : kend - rs->m.hi;
    hi = fs->m.hi < kend - rs->m.lo ? fs->m.hi : kend - rs->m.lo;
    for(k = lo; k <= hi; k++) {
      int32_t hf = fs->m.off[k - fs->m.lo];
      int32_t hr = rs->m.off[(kend - k) - rs->m.lo];
      if(hf >= 0 && hr >= 0 && hf + hr >= nb) {
        best->s = sf + sr;
        best->h = hf;
        best->v = hf - k;
        best->state = WFA_M;
        break;
      }
    }
  }
  if(sf + sr - po < best->s && fs->i.lo <= fs->i.hi && rs->i.lo <= rs->i.hi) {
    lo = fs->i.lo > kend - rs->i.hi ? fs->i.lo : kend - rs->i.hi;
    hi = fs->i.hi < kend - rs->i.lo ? fs->i.hi : kend - rs->i.lo;
    for(k = lo; k <= hi; k++) {
      int32_t hf = fs->i.off[k - fs->i.lo];
      int32_t hr = rs->i.off[(kend - k) - rs->i.lo];
      if(hf >= 1 && hr >= 1 && hf + hr >= nb) {
        best->s = sf + sr - po;
        best->h = hf;
        best->v = hf - k;
        best->state = WFA_I;
        break;
      }
    }
  }
  if(sf + sr - po < best->s && fs->d.lo <= fs->d.hi && rs->d.lo <= rs->d.hi) {
    lo = fs->d.lo > kend - rs->d.hi ? fs->d.lo : kend - rs->d.hi;
    hi = fs->d.hi < kend - rs->d.lo ? fs->d.hi : kend - rs->d.lo;
    for(k = lo; k <= hi; k++) {
      int32_t hf = fs->d.off[k - fs->d.lo];
      int32_t hr = rs->d.off[(kend - k) - rs->d.lo];
      if(hf >= 0 && hr >= 0 && hf + hr >= nb && (hf - k) >= 1 && (hr - (kend - k)) >= 1) {
        best->s = sf + sr - po;
        best->h = hf;
        best->v = hf - k;
        best->state = WFA_D;
        break;
      }
    }
  }
}

// O(s)-space Bidirectional WFA (Marco-Sola et al., 2023)
static size_t biwfa_rec(const biwfa_ctx_t *c, const biwfa_sub_t *q,
                        char *out_a, char *out_b, long long *penalty_out)
{
  const uint8_t *const a = q->a, *const ra = q->ra;
  const uint8_t *const b = q->b, *const rb = q->rb;
  const char *const sa = q->sa, *const sb = q->sb;
  const int32_t m = q->m, nb = q->nb;
  const long long px = c->px, po = c->po, pe = c->pe;
  const wfa_state_t s0 = q->s0, s1 = q->s1;
  wf_null_t *const null_wf = c->null_wf;
  size_t *const budget_left = c->budget_left;
  if(m == 0 && nb == 0) { if(penalty_out) *penalty_out = 0; return 0; }
  if(m == 0) {
    int32_t j;
    for(j = 0; j < nb; j++) { out_a[j] = '-'; out_b[j] = sb[j]; }
    if(penalty_out) *penalty_out = (s0 == WFA_I ? 0 : po) + (long long)nb * pe;
    return (size_t)nb;
  }
  if(nb == 0) {
    int32_t i;
    for(i = 0; i < m; i++) { out_a[i] = sa[i]; out_b[i] = '-'; }
    if(penalty_out) *penalty_out = (s0 == WFA_D ? 0 : po) + (long long)m * pe;
    return (size_t)m;
  }
  if((size_t)m + (size_t)nb <= 2048) return biwfa_base(c, q, out_a, out_b);

  const size_t ring = (size_t)(px > po + pe ? px : po + pe) + 1;
  wf_slot_t *fwd = aln_calloc(ring, sizeof(wf_slot_t));
  wf_slot_t *rev = aln_calloc(ring, sizeof(wf_slot_t));
  if(fwd == NULL || rev == NULL) {
    free(fwd);
    free(rev);
    null_wf->nomem = true;
    return 0;
  }
  const wavefront_t none = {1, 0, 0, NULL, NULL};
  size_t i;
  for(i = 0; i < ring; i++) {
    fwd[i].m = fwd[i].i = fwd[i].d = none;
    rev[i].m = rev[i].i = rev[i].d = none;
  }

  biwfa_meet_t meet = {LLONG_MAX, 0, 0, WFA_M};
  const wf_dir_t df = {a, b, m, nb, px, po, pe, q->fwd_pad, s0, fwd, ring,
                       false, null_wf};
  const wf_dir_t dr = {ra, rb, m, nb, px, po, pe, q->rev_pad, s1, rev, ring,
                       false, null_wf};
  long long sf = -1, sr = -1;
  int32_t max_hf = 0, max_hr = 0;
  size_t empty_f = 0, empty_r = 0;
  bool ok = true;

  for(;;) {
    sf++;
    size_t w = wf_step(&df, sf, &max_hf);
    if(w == SIZE_MAX) null_wf->nomem = true;
    if(w == SIZE_MAX || w > *budget_left) { ok = false; break; }
    *budget_left -= w;
    if(w == 0) { if(++empty_f >= ring) { ok = false; break; } }
    else empty_f = 0;
    if(sr >= 0 && w > 0 && max_hf + max_hr >= nb) {
      long long t0 = sr - (long long)ring + 1;
      if(t0 < 0) t0 = 0;
      for(long long t = t0; t <= sr; t++)
        biwfa_check(&fwd[(size_t)sf % ring], sf, &rev[(size_t)t % ring], t,
                    m, nb, po, &meet);
    }
    if(meet.s < LLONG_MAX && sf + sr >= meet.s + po + (long long)ring) break;

    sr++;
    w = wf_step(&dr, sr, &max_hr);
    if(w == SIZE_MAX) null_wf->nomem = true;
    if(w == SIZE_MAX || w > *budget_left) { ok = false; break; }
    *budget_left -= w;
    if(w == 0) { if(++empty_r >= ring) { ok = false; break; } }
    else empty_r = 0;
    if(w > 0 && max_hf + max_hr >= nb) {
      long long t0 = sf - (long long)ring + 1;
      if(t0 < 0) t0 = 0;
      for(long long t = t0; t <= sf; t++)
        biwfa_check(&fwd[(size_t)t % ring], t, &rev[(size_t)sr % ring], sr,
                    m, nb, po, &meet);
    }
    if(meet.s < LLONG_MAX && sf + sr >= meet.s + po + (long long)ring) break;
    if(sf > (long long)INT32_MAX / 2) { ok = false; break; }
  }

  wfa_free(fwd, ring);
  wfa_free(rev, ring);
  if(!ok || meet.s == LLONG_MAX) return 0;
  if(penalty_out) *penalty_out = meet.s;

  const int32_t bp_v = meet.v, bp_h = meet.h;
  if((bp_v <= 0 && bp_h <= 0) || (bp_v >= m && bp_h >= nb))
    return biwfa_base(c, q, out_a, out_b);

  // The two halves: up to the meeting cell (its reversed strings start
  // there, and the forward ones' padding is now out of reach), and from it
  biwfa_sub_t left = *q, right = *q;
  left.ra = ra + (m - bp_v); left.m = bp_v;
  left.rb = rb + (nb - bp_h); left.nb = bp_h;
  left.s1 = meet.state; left.fwd_pad = false;
  right.a = a + bp_v; right.sa = sa + bp_v; right.m = m - bp_v;
  right.b = b + bp_h; right.sb = sb + bp_h; right.nb = nb - bp_h;
  right.s0 = meet.state; right.rev_pad = false;
  size_t nl = biwfa_rec(c, &left, out_a, out_b, NULL);
  if(nl == 0 && (bp_v > 0 || bp_h > 0)) return 0;
  size_t nr = biwfa_rec(c, &right, out_a + nl, out_b + nl, NULL);
  if(nr == 0 && (m - bp_v > 0 || nb - bp_h > 0)) return 0;
  return nl + nr;
}

aln_status_t alignment_wfa_align(const sub_table_t *table,
                                 const scoring_t *scoring,
                                 const char *seq_a, size_t len_a,
                                 const char *seq_b, size_t len_b,
                                 size_t budget, alignment_t *result)
{
  wfa_pair_t pair;
  bool nomem = false;
  if(!wfa_setup(table, scoring, seq_a, len_a, seq_b, len_b, &pair, &nomem))
    return nomem ? ALN_ERR_NOMEM : ALN_ERR_UNSUPPORTED;

  if(len_a + len_b <= 2048) {
    wf_slot_t *slots;
    size_t nslots;
    aln_status_t status = ALN_ERR_UNSUPPORTED;
    score_t score;
    const long long penalty = wfa_fill(&pair, budget, true, &slots, &nslots,
                                       &nomem);
    if(nomem) status = ALN_ERR_NOMEM;
    else if(penalty >= 0 && wfa_best(&pair, penalty, &score)) {
      if(!alignment_reserve(result, len_a + len_b)) status = ALN_ERR_NOMEM;
      else {
        size_t n = wfa_trace(&pair, slots, penalty, seq_a, seq_b,
                             result->result_a, result->result_b);
        if(n > 0) {
          size_t i;
          for(i = 0; i < n / 2; i++) {
            char t = result->result_a[i];
            result->result_a[i] = result->result_a[n - 1 - i];
            result->result_a[n - 1 - i] = t;
            t = result->result_b[i];
            result->result_b[i] = result->result_b[n - 1 - i];
            result->result_b[n - 1 - i] = t;
          }
          result->result_a[n] = result->result_b[n] = '\0';
          result->length = n;
          result->score = score;
          result->pos_a = result->pos_b = 0;
          result->len_a = len_a;
          result->len_b = len_b;
          status = ALN_OK;
        }
      }
    }
    wfa_free(slots, nslots);
    free(pair.a);
    return status;
  }

  if(!alignment_reserve(result, len_a + len_b)) {
    free(pair.a);
    return ALN_ERR_NOMEM;
  }
  uint8_t *ra = aln_malloc(len_a + len_b + 64);
  if(ra == NULL) { free(pair.a); return ALN_ERR_NOMEM; }
  uint8_t *rb = ra + len_a + 32;
  size_t i;
  for(i = 0; i < len_a; i++) ra[i] = pair.a[len_a - 1 - i];
  for(i = 0; i < len_b; i++) rb[i] = pair.b[len_b - 1 - i];
  memset(ra + len_a, 0xFE, 32);
  memset(rb + len_b, 0xFF, 32);

  wf_null_t null_wf = {NULL, 0, false};
  size_t budget_left = budget;
  long long penalty = -1;
  const biwfa_ctx_t ctx = {pair.px, pair.po, pair.pe, &null_wf, &budget_left};
  const biwfa_sub_t whole = {pair.a, ra, seq_a, pair.m, pair.b, rb, seq_b,
                             pair.nb, WFA_M, WFA_M, true, true};
  size_t n = biwfa_rec(&ctx, &whole, result->result_a, result->result_b,
                       &penalty);
  free(null_wf.buf);
  free(ra);
  aln_status_t status = null_wf.nomem ? ALN_ERR_NOMEM : ALN_ERR_UNSUPPORTED;
  score_t score;
  if(!null_wf.nomem && n > 0 && penalty >= 0 &&
     wfa_best(&pair, penalty, &score)) {
    result->result_a[n] = result->result_b[n] = '\0';
    result->length = n;
    result->score = score;
    result->pos_a = result->pos_b = 0;
    result->len_a = len_a;
    result->len_b = len_b;
    status = ALN_OK;
  }
  free(pair.a);
  return status;
}
