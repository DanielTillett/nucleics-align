/*
 alignment_wfa.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// The score of the best global alignment by the wavefront algorithm (WFA:
// Marco-Sola et al. 2021, Bioinformatics 37:456-463), written from the
// paper, for sizing a band (see alignment_fill_banded()). Its time grows
// with the alignment's penalty, so it is quick for similar sequences and
// gives up past a budget.
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

bool alignment_wfa_score(const sub_table_t *table, const scoring_t *scoring,
                         const char *seq_a, size_t len_a,
                         const char *seq_b, size_t len_b,
                         size_t budget, score_t *score)
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
  const long long px = 2 * x, po = 2 * o, pe = 2 * e + l;
  if(px <= 0 || po < 0 || pe <= 0 || px > 1000000 || po + pe > 1000000)
    return false;

  // Class numbers of the letters, then 8 padding bytes each, different
  const int32_t m = (int32_t)len_a, nb = (int32_t)len_b;
  uint8_t *a = aln_malloc(len_a + len_b + 16);
  if(a == NULL) return false;
  uint8_t *b = a + len_a + 8;
  for(i = 0; i < len_a; i++)
    a[i] = cls[table->code[(unsigned char)seq_a[i]]];
  for(i = 0; i < len_b; i++)
    b[i] = cls[table->code[(unsigned char)seq_b[i]]];
  memset(a + len_a, 0xFE, 8);
  memset(b + len_b, 0xFF, 8);

  // Wavefronts of the last px or po + pe scores
  const size_t ring = (size_t)(px > po + pe ? px : po + pe) + 1;
  wf_slot_t *slots = aln_calloc(ring, sizeof(wf_slot_t));
  bool found = false;
  long long s = 0;
  size_t work = 0;
  const int32_t kend = nb - m;
  if(slots == NULL) goto done;
  for(i = 0; i < ring; i++)
    slots[i].m.lo = slots[i].i.lo = slots[i].d.lo = 1; // none
  for(i = 0; i < ring; i++)
    slots[i].m.hi = slots[i].i.hi = slots[i].d.hi = 0;

  for(s = 0; ; s++) {
    wf_slot_t *cur = &slots[s % ring];
    const wavefront_t none = {1, 0, 0, NULL};
    const wavefront_t *mx = s >= px ? &slots[(s - px) % ring].m : &none;
    const wavefront_t *mo = s >= po + pe ? &slots[(s - po - pe) % ring].m
                                         : &none;
    const wavefront_t *ie = s >= pe ? &slots[(s - pe) % ring].i : &none;
    const wavefront_t *de = s >= pe ? &slots[(s - pe) % ring].d : &none;
    int32_t lo, hi, k;

    if(s == 0) {
      if(!wf_reserve(&cur->m, 0, 0)) goto done;
      cur->m.off[0] = wf_extend(a, b, 0, 0);
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
      if(work > budget) goto done;
      if(!wf_reserve(&cur->m, lo, hi) || !wf_reserve(&cur->i, lo, hi) ||
         !wf_reserve(&cur->d, lo, hi)) goto done;
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
          int32_t run = wf_extend(a, b, best - k, best);
          work += (size_t)run / 8;
          best += run;
        }
        else best = WFA_NONE;
        cur->m.off[k - lo] = best;
      }
    }
    if(wf_get(&cur->m, kend) >= nb) { found = true; break; }
    if(s > (long long)INT32_MAX / 2) goto done;
  }

done:
  if(slots != NULL) {
    for(i = 0; i < ring; i++) {
      free(slots[i].m.off);
      free(slots[i].i.off);
      free(slots[i].d.off);
    }
    free(slots);
  }
  free(a);
  if(!found) return false;
  const long long best = (l * ((long long)len_a + (long long)len_b) - s) / 2;
  if(best < INT_MIN || best > INT_MAX) return false;
  *score = (score_t)best;
  return true;
}
