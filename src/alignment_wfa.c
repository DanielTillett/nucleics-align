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
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include "alignment.h"
#include "alignment_alloc.h"

// Diagonal k = h - v of text b (h) and pattern a (v); offsets are h
#define WFA_NONE (INT32_MIN / 2)

typedef struct {
  int32_t lo, hi; // diagonals held; none if lo > hi
  size_t cap;
  int32_t *off; // off[k - lo]
} wavefront_t;

// The wavefronts of one score: match (M), insertion (I) and deletion (D)
typedef struct { wavefront_t m, i, d; } wf_slot_t;

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
    int32_t *off = aln_realloc(w->off, need * sizeof(int32_t));
    if(off == NULL) return false;
    w->off = off;
    w->cap = need;
  }
  return true;
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
  long long l, px, po, pe;
} wfa_pair_t;

// Set pair up for seq_a and seq_b, or return false where the scoring is not
// one match and one mismatch score with gap penalties (or out of memory)
static bool wfa_setup(const sub_table_t *table, const scoring_t *scoring,
                      const char *seq_a, size_t len_a,
                      const char *seq_b, size_t len_b, wfa_pair_t *pair)
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

  // Class numbers of the letters, then 8 padding bytes each, different
  pair->m = (int32_t)len_a;
  pair->nb = (int32_t)len_b;
  pair->a = aln_malloc(len_a + len_b + 16);
  if(pair->a == NULL) return false;
  pair->b = pair->a + len_a + 8;
  for(i = 0; i < len_a; i++)
    pair->a[i] = cls[table->code[(unsigned char)seq_a[i]]];
  for(i = 0; i < len_b; i++)
    pair->b[i] = cls[table->code[(unsigned char)seq_b[i]]];
  memset(pair->a + len_a, 0xFE, 8);
  memset(pair->b + len_b, 0xFF, 8);
  return true;
}

// The wavefronts of a pair, score by score, until one reaches the last
// cell. With keep, every score's wavefronts are kept (slots[s], for the
// traceback); otherwise only the last px or po + pe scores', in a ring.
// Returns the least penalty, or -1 past budget diagonal updates (cells
// kept, with keep) or if out of memory. *slots and *nslots hold what to
// free (wfa_free())
static long long wfa_fill(const wfa_pair_t *pair, size_t budget, bool keep,
                          wf_slot_t **slots_out, size_t *nslots_out)
{
  const int32_t m = pair->m, nb = pair->nb, kend = nb - m;
  const long long px = pair->px, po = pair->po, pe = pair->pe;
  const size_t ring = (size_t)(px > po + pe ? px : po + pe) + 1;
  size_t nslots = keep ? 0 : ring, i;
  wf_slot_t *slots = keep ? NULL : aln_calloc(ring, sizeof(wf_slot_t));
  const wavefront_t none = {1, 0, 0, NULL};
  long long s, result = -1;
  size_t work = 0;
  *slots_out = slots;
  *nslots_out = nslots;
  if(!keep && slots == NULL) return -1;
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
        if(g == NULL) break;
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
      if(!wf_reserve(&cur->m, 0, 0)) break;
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
        continue;
      }
      if(lo < -m) lo = -m;
      if(hi > nb) hi = nb;
      work += (size_t)(hi - lo + 1);
      if(work > budget) break;
      if(!wf_reserve(&cur->m, lo, hi) || !wf_reserve(&cur->i, lo, hi) ||
         !wf_reserve(&cur->d, lo, hi)) break;
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
    free(slots[i].m.off);
    free(slots[i].i.off);
    free(slots[i].d.off);
  }
  free(slots);
}

// The best score, (l(m+n) - P)/2, if it fits a score_t
static bool wfa_best(const wfa_pair_t *pair, long long penalty,
                     score_t *score)
{
  const long long best = (pair->l * ((long long)pair->m + pair->nb) -
                          penalty) / 2;
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
  if(!wfa_setup(table, scoring, seq_a, len_a, seq_b, len_b, &pair))
    return false;
  wf_slot_t *slots;
  size_t nslots;
  const long long penalty = wfa_fill(&pair, budget, false, &slots, &nslots);
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

aln_status_t alignment_wfa_align(const sub_table_t *table,
                                 const scoring_t *scoring,
                                 const char *seq_a, size_t len_a,
                                 const char *seq_b, size_t len_b,
                                 size_t budget, alignment_t *result)
{
  wfa_pair_t pair;
  if(!wfa_setup(table, scoring, seq_a, len_a, seq_b, len_b, &pair))
    return ALN_ERR_UNSUPPORTED;
  wf_slot_t *slots;
  size_t nslots;
  aln_status_t status = ALN_ERR_UNSUPPORTED;
  score_t score;
  const long long penalty = wfa_fill(&pair, budget, true, &slots, &nslots);
  if(penalty >= 0 && wfa_best(&pair, penalty, &score)) {
    if(!alignment_reserve(result, len_a + len_b)) status = ALN_ERR_NOMEM;
    else {
      size_t n = wfa_trace(&pair, slots, penalty, seq_a, seq_b,
                           result->result_a, result->result_b);
      if(n > 0) {
        // Reversed: put it the right way round
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
