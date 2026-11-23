/*
 alignment_myers.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// Global unit-cost edit distance (match 0, mismatch or gap -1, gap open 0)
// by Myers' bit-vector algorithm (Myers 1999, J ACM 46:395-415), with
// Ukkonen exponential diagonal banding for --scoreonly and full-row output
// for Hirschberg --linear.

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "alignment.h"
#include "alignment_alloc.h"

// Whether the scoring is unit cost for these sequences' letters: gap open
// 0, gap extension -1, and every pair the table holds scoring 0 or -1
static bool myers_unit_cost(const sub_table_t *table, const scoring_t *scoring)
{
  if(scoring->gap_open != 0 || scoring->gap_extend != -1 ||
     scoring->no_start_gap_penalty || scoring->no_end_gap_penalty ||
     scoring->no_gaps_in_a || scoring->no_gaps_in_b ||
     scoring->no_gaps_in_longest || scoring->no_mismatches) return false;
  size_t i;
  for(i = 0; i < table->n * table->n; i++)
    if(table->scores[i] != 0 && table->scores[i] != -1) return false;
  return true;
}

static inline void myers_build_eq(const sub_table_t *table,
                                  const char *seq_a, size_t len_a,
                                  bool rev, uint64_t *eq, size_t nw)
{
  const size_t n = table->n;
  size_t i, c;
  for(i = 0; i < len_a; i++) {
    const unsigned char ch = (unsigned char)seq_a[rev ? (len_a - 1 - i) : i];
    const size_t ca = table->code[ch];
    const size_t w = i / 64;
    const uint64_t bit = (uint64_t)1 << (i % 64);
    for(c = 0; c < n; c++)
      if(table->scores[c * n + ca] == 0)
        eq[c * nw + w] |= bit;
  }
}

// Run one banded Myers pass with edit-distance threshold K >= |m - n|.
// Returns true and writes exact edit distance to *dist_out iff dist <= K.
static bool myers_banded_pass(const sub_table_t *table,
                              const uint64_t *eq, size_t nw,
                              size_t len_a, const char *seq_b, size_t len_b,
                              long long K, uint64_t *P, uint64_t *M,
                              long long *bscore, long long *dist_out)
{
  const long long m = (long long)len_a, nb = (long long)len_b;
  const long long d = m - nb;
  const long long ad = d < 0 ? -d : d;
  if(K < ad) return false;
  const long long slack = (K - ad) / 2 + 1;
  const long long lo_off = (d < 0 ? d : 0) - slack;
  const long long hi_off = (d > 0 ? d : 0) + slack;
  const uint64_t last = (uint64_t)1 << ((len_a - 1) % 64);
  const long long last_rows = (long long)((len_a - 1) % 64 + 1);

  // Column 0: initialize words up to the initial band bottom
  long long r0_max = hi_off;
  size_t prev_hi = r0_max >= 1 ? (size_t)((r0_max - 1) / 64) : 0;
  if(prev_hi >= nw) prev_hi = nw - 1;
  size_t w, j;
  for(w = 0; w <= prev_hi; w++) {
    P[w] = ~(uint64_t)0;
    M[w] = 0;
    bscore[w] = (w + 1 == nw) ? m : (long long)((w + 1) * 64);
  }

  for(j = 0; j < len_b; j++) {
    const long long col = (long long)j + 1;
    const long long r_min = col + lo_off;
    const long long r_max = col + hi_off;
    size_t w_lo = r_min > 1 ? (size_t)((r_min - 1) / 64) : 0;
    size_t w_hi = r_max >= 1 ? (size_t)((r_max - 1) / 64) : 0;
    if(w_hi >= nw) w_hi = nw - 1;
    if(w_lo > w_hi) return false;

    // Initialize newly entered words from column j-1's bottom score
    if(w_hi > prev_hi) {
      for(w = prev_hi + 1; w <= w_hi; w++) {
        P[w] = ~(uint64_t)0;
        M[w] = 0;
        bscore[w] = bscore[w - 1] + (w + 1 == nw ? last_rows : 64);
      }
      prev_hi = w_hi;
    }

    const uint64_t *e = eq + table->code[(unsigned char)seq_b[j]] * nw;
    int hin = 1; // +1 at row 0 (gap along b) and safe upper bound at w_lo > 0
    for(w = w_lo; w <= w_hi; w++) {
      uint64_t Pv = P[w], Mv = M[w], Eq = e[w];
      const uint64_t Xv = Eq | Mv;
      if(hin < 0) Eq |= 1;
      const uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
      uint64_t Ph = Mv | ~(Xh | Pv);
      uint64_t Mh = Pv & Xh;
      const uint64_t top = w + 1 == nw ? last : (uint64_t)1 << 63;
      const int hout = (Ph & top) ? 1 : (Mh & top) ? -1 : 0;
      Ph <<= 1;
      Mh <<= 1;
      if(hin < 0) Mh |= 1;
      else if(hin > 0) Ph |= 1;
      P[w] = Mh | ~(Xv | Ph);
      M[w] = Ph & Xv;
      bscore[w] += hout;
      hin = hout;
    }
  }

  if(prev_hi < nw - 1) return false;
  if(bscore[nw - 1] <= K) {
    *dist_out = bscore[nw - 1];
    return true;
  }
  return false;
}

static bool myers_score_impl(const sub_table_t *table, const scoring_t *scoring,
                             const char *seq_a, size_t len_a,
                             const char *seq_b, size_t len_b,
                             bool banded_only, score_t *score)
{
  if(!myers_unit_cost(table, scoring) || len_a == 0 ||
     len_a + len_b > (size_t)SCORE_ABS_MAX) return false;

  if(len_b == 0) {
    *score = -(score_t)len_a;
    return true;
  }

  const size_t n = table->n, nw = (len_a + 63) / 64;
  if(nw > SIZE_MAX / sizeof(uint64_t) / (n + 2)) return false;
  uint64_t *eq = aln_calloc((n + 2) * nw, sizeof(uint64_t));
  if(eq == NULL) return false;
  uint64_t *P = eq + n * nw, *M = P + nw;
  myers_build_eq(table, seq_a, len_a, false, eq, nw);

  // Ukkonen exponential banding when len_a spans enough words to pay off
  if(nw > 4) {
    long long *bscore = aln_malloc(nw * sizeof(long long));
    if(bscore != NULL) {
      const long long d = (long long)len_a - (long long)len_b;
      const long long ad = d < 0 ? -d : d;
      long long K = ad > 64 ? ad : 64;
      const long long limit = (long long)(len_a < len_b ? len_a : len_b) / 4 + ad;
      long long dist = 0;
      while(K < limit) {
        if(myers_banded_pass(table, eq, nw, len_a, seq_b, len_b,
                             K, P, M, bscore, &dist)) {
          free(bscore);
          free(eq);
          *score = (score_t)-dist;
          return true;
        }
        K = K > (long long)SCORE_ABS_MAX / 2 ? (long long)SCORE_ABS_MAX : K * 2;
      }
      free(bscore);
    }
    if(banded_only) {
      free(eq);
      return false;
    }
  }

  // Full-width Myers pass
  size_t w, j;
  for(w = 0; w < nw; w++) { P[w] = ~(uint64_t)0; M[w] = 0; }
  const uint64_t last = (uint64_t)1 << ((len_a - 1) % 64);
  long long dist = (long long)len_a;

  for(j = 0; j < len_b; j++) {
    const uint64_t *e = eq + table->code[(unsigned char)seq_b[j]] * nw;
    int hin = 1;
    for(w = 0; w < nw; w++) {
      uint64_t Pv = P[w], Mv = M[w], Eq = e[w];
      const uint64_t Xv = Eq | Mv;
      if(hin < 0) Eq |= 1;
      const uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
      uint64_t Ph = Mv | ~(Xh | Pv);
      uint64_t Mh = Pv & Xh;
      const uint64_t top = w + 1 == nw ? last : (uint64_t)1 << 63;
      const int hout = (Ph & top) ? 1 : (Mh & top) ? -1 : 0;
      Ph <<= 1;
      Mh <<= 1;
      if(hin < 0) Mh |= 1;
      else if(hin > 0) Ph |= 1;
      P[w] = Mh | ~(Xv | Ph);
      M[w] = Ph & Xv;
      hin = hout;
    }
    dist += hin;
  }

  free(eq);
  *score = (score_t)-dist;
  return true;
}

bool alignment_myers_score(const sub_table_t *table, const scoring_t *scoring,
                           const char *seq_a, size_t len_a,
                           const char *seq_b, size_t len_b, score_t *score)
{
  return myers_score_impl(table, scoring, seq_a, len_a, seq_b, len_b,
                          false, score);
}

bool alignment_myers_score_banded(const sub_table_t *table,
                                  const scoring_t *scoring,
                                  const char *seq_a, size_t len_a,
                                  const char *seq_b, size_t len_b,
                                  score_t *score)
{
  return myers_score_impl(table, scoring, seq_a, len_a, seq_b, len_b,
                          true, score);
}

bool alignment_myers_row(const sub_table_t *table, const scoring_t *scoring,
                         const char *seq_a, size_t len_a,
                         const char *seq_b, size_t len_b,
                         bool rev, score_t *cc)
{
  if(!myers_unit_cost(table, scoring) || len_a == 0 ||
     len_a + len_b > (size_t)SCORE_ABS_MAX) return false;

  const size_t n = table->n, nw = (len_a + 63) / 64;
  if(nw > SIZE_MAX / sizeof(uint64_t) / (n + 2)) return false;
  uint64_t *eq = aln_calloc((n + 2) * nw, sizeof(uint64_t));
  if(eq == NULL) return false;
  uint64_t *P = eq + n * nw, *M = P + nw;
  myers_build_eq(table, seq_a, len_a, rev, eq, nw);

  size_t w, j;
  for(w = 0; w < nw; w++) P[w] = ~(uint64_t)0;
  const uint64_t last = (uint64_t)1 << ((len_a - 1) % 64);
  long long dist = (long long)len_a;
  if(rev) cc[len_b] = (score_t)-dist;
  else cc[0] = (score_t)-dist;

  for(j = 0; j < len_b; j++) {
    const unsigned char ch = (unsigned char)seq_b[rev ? (len_b - 1 - j) : j];
    const uint64_t *e = eq + table->code[ch] * nw;
    int hin = 1;
    for(w = 0; w < nw; w++) {
      uint64_t Pv = P[w], Mv = M[w], Eq = e[w];
      const uint64_t Xv = Eq | Mv;
      if(hin < 0) Eq |= 1;
      const uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
      uint64_t Ph = Mv | ~(Xh | Pv);
      uint64_t Mh = Pv & Xh;
      const uint64_t top = w + 1 == nw ? last : (uint64_t)1 << 63;
      const int hout = (Ph & top) ? 1 : (Mh & top) ? -1 : 0;
      Ph <<= 1;
      Mh <<= 1;
      if(hin < 0) Mh |= 1;
      else if(hin > 0) Ph |= 1;
      P[w] = Mh | ~(Xv | Ph);
      M[w] = Ph & Xv;
      hin = hout;
    }
    dist += hin;
    if(rev) cc[len_b - 1 - j] = (score_t)-dist;
    else cc[j + 1] = (score_t)-dist;
  }

  free(eq);
  return true;
}
