/*
 alignment_myers.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// Global unit-cost edit distance (match 0, mismatch or gap -1, gap open 0)
// by Myers' bit-vector algorithm (Myers 1999, J ACM 46:395-415), with
// Ukkonen exponential diagonal banding, active-block frontier trimming, and
// branch-free 64-bit word updates for --scoreonly, plus full-row output for
// Hirschberg --linear.

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

// Branch-free 64-bit Myers word update for a full 64-row block.
static inline int myers_step_word(uint64_t *P, uint64_t *M,
                                  uint64_t Eq, int hin)
{
  const uint64_t Pv = *P, Mv = *M;
  const uint64_t Xv = Eq | Mv;
  const uint64_t hneg = (uint64_t)(hin < 0);
  const uint64_t hpos = (uint64_t)(hin > 0);
  Eq |= hneg;
  const uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
  uint64_t Ph = Mv | ~(Xh | Pv);
  uint64_t Mh = Pv & Xh;
  const int hout = (int)(Ph >> 63) - (int)(Mh >> 63);
  Ph = (Ph << 1) | hpos;
  Mh = (Mh << 1) | hneg;
  *P = Mh | ~(Xv | Ph);
  *M = Ph & Xv;
  return hout;
}

// Branch-free 64-bit Myers word update for the final (possibly partial) block,
// where shift = (len_a - 1) % 64.
static inline int myers_step_last(uint64_t *P, uint64_t *M,
                                  uint64_t Eq, int hin, unsigned shift)
{
  const uint64_t Pv = *P, Mv = *M;
  const uint64_t Xv = Eq | Mv;
  const uint64_t hneg = (uint64_t)(hin < 0);
  const uint64_t hpos = (uint64_t)(hin > 0);
  Eq |= hneg;
  const uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
  uint64_t Ph = Mv | ~(Xh | Pv);
  uint64_t Mh = Pv & Xh;
  const int hout = (int)((Ph >> shift) & 1) - (int)((Mh >> shift) & 1);
  Ph = (Ph << 1) | hpos;
  Mh = (Mh << 1) | hneg;
  *P = Mh | ~(Xv | Ph);
  *M = Ph & Xv;
  return hout;
}

// Exact lower bound on the total edit distance of any global path passing
// through block w (rows [r_lo, r_hi]) at column `col`, given its bottom cell
// distance `bs = D(r_hi, col)`. Because D(r, col) >= bs - (r_hi - r) for all
// r in [r_lo, r_hi], min_{r in [r_lo, r_hi]} (D(r, col) + |r - r_tgt|) is
// minimised at r_best = clamp(r_tgt, r_lo, r_hi), where r_tgt = m - (nb - col).
static inline long long myers_block_lb(size_t w, size_t nw, long long m,
                                       long long r_tgt, long long bs)
{
  const long long r_lo = (long long)(w * 64) + 1;
  const long long r_hi = (w + 1 == nw) ? m : (long long)((w + 1) * 64);
  const long long r_best = r_tgt < r_lo ? r_lo : (r_tgt > r_hi ? r_hi : r_tgt);
  const long long dg = r_best >= r_tgt ? (r_best - r_tgt) : (r_tgt - r_best);
  return bs - (r_hi - r_best) + dg;
}

// Run one banded Myers pass with edit-distance threshold K >= |m - n|,
// trimming inactive blocks at both ends and returning false immediately if
// the feasible frontier empties.
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
  const unsigned last_shift = (unsigned)((len_a - 1) % 64);
  const long long last_rows = (long long)(last_shift + 1);

  // Active block range [f_lo, f_hi] at column 0: rows up to (K + d)/2 + 1
  // (plus safety margin for block granularity)
  const long long r0_max = (K + d) / 2 + 1;
  size_t f_lo = 0;
  size_t f_hi = r0_max >= 1 ? (size_t)((r0_max - 1) / 64) : 0;
  if(f_hi >= nw) f_hi = nw - 1;
  size_t w, j;
  for(w = 0; w <= f_hi; w++) {
    P[w] = ~(uint64_t)0;
    M[w] = 0;
    bscore[w] = (w + 1 == nw) ? m : (long long)((w + 1) * 64);
  }

  for(j = 0; j < len_b; j++) {
    const long long col = (long long)j + 1;
    const long long r_tgt = m - (nb - col);
    const uint64_t *e = eq + table->code[(unsigned char)seq_b[j]] * nw;
    int hin = 1; // +1 at row 0 (gap along b) and safe upper bound at f_lo > 0

    // Run full words in [f_lo, f_hi], splitting the partial last word out
    const size_t full_end = (f_hi + 1 == nw) ? (nw - 1) : (f_hi + 1);
    for(w = f_lo; w < full_end; w++) {
      const int hout = myers_step_word(&P[w], &M[w], e[w], hin);
      bscore[w] += hout;
      hin = hout;
    }
    if(f_hi + 1 == nw) {
      const int hout = myers_step_last(&P[nw - 1], &M[nw - 1], e[nw - 1],
                                       hin, last_shift);
      bscore[nw - 1] += hout;
      hin = hout;
    }

    // If the bottom active block can still reach the next block within K,
    // extend f_hi by 1 block and compute it in this same column.
    if(f_hi + 1 < nw) {
      const long long prev_bs = bscore[f_hi];
      // A path entering block f_hi + 1 from cell (r_hi(f_hi), col) has score
      // at least prev_bs - 1 (or prev_bs + 1 from above); check if the top
      // cell of block f_hi + 1 can still complete within K.
      const long long next_r_lo = (long long)((f_hi + 1) * 64) + 1;
      const long long dg = next_r_lo >= r_tgt ? (next_r_lo - r_tgt)
                                              : (r_tgt - next_r_lo);
      if(prev_bs - 1 + dg <= K) {
        f_hi++;
        P[f_hi] = ~(uint64_t)0;
        M[f_hi] = 0;
        const long long bh = (f_hi + 1 == nw) ? last_rows : 64;
        // Bottom cell of block f_hi at column col-1: initialized from
        // (bscore[f_hi-1] - hin) + bh, then updated by this column's hout.
        bscore[f_hi] = (prev_bs - hin) + bh;
        const int hout = (f_hi + 1 == nw)
                         ? myers_step_last(&P[f_hi], &M[f_hi], e[f_hi],
                                           hin, last_shift)
                         : myers_step_word(&P[f_hi], &M[f_hi], e[f_hi], hin);
        bscore[f_hi] += hout;
      }
    }

    // Row 0 has exact distance `col` and remaining imbalance `|r_tgt|`; while
    // row 0 can still complete within K, block 0 must stay active so row 0
    // can feed row 1 on later columns (e.g. a = A^257, b = C^128 A^257).
    const bool row0_viable = (col + (r_tgt >= 0 ? r_tgt : -r_tgt) <= K);

    // Trim inactive blocks from the bottom
    while(myers_block_lb(f_hi, nw, m, r_tgt, bscore[f_hi]) > K) {
      if(f_hi == f_lo) {
        if(f_lo == 0 && row0_viable) break;
        return false;
      }
      f_hi--;
    }
    // Trim inactive blocks from the top (paths never move up, so trimmed
    // top blocks are never re-entered)
    while(myers_block_lb(f_lo, nw, m, r_tgt, bscore[f_lo]) > K) {
      if(f_lo == 0 && row0_viable) break;
      if(f_lo == f_hi) return false;
      f_lo++;
    }
  }

  if(f_hi < nw - 1) return false;
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

  // Ukkonen exponential banding with active-block trimming
  if(nw > 4) {
    long long *bscore = aln_malloc(nw * sizeof(long long));
    if(bscore != NULL) {
      const long long d = (long long)len_a - (long long)len_b;
      const long long ad = d < 0 ? -d : d;
      long long K = ad > 64 ? ad : 64;
      const long long limit = (long long)(len_a < len_b ? len_a : len_b) / 2 + ad;
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

  // Full-width branch-free Myers pass
  size_t w, j;
  for(w = 0; w < nw; w++) { P[w] = ~(uint64_t)0; M[w] = 0; }
  const unsigned last_shift = (unsigned)((len_a - 1) % 64);
  const size_t full_nw = nw - 1;
  long long dist = (long long)len_a;

  for(j = 0; j < len_b; j++) {
    const uint64_t *e = eq + table->code[(unsigned char)seq_b[j]] * nw;
    int hin = 1;
    for(w = 0; w < full_nw; w++)
      hin = myers_step_word(&P[w], &M[w], e[w], hin);
    hin = myers_step_last(&P[full_nw], &M[full_nw], e[full_nw],
                          hin, last_shift);
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
  const unsigned last_shift = (unsigned)((len_a - 1) % 64);
  const size_t full_nw = nw - 1;
  long long dist = (long long)len_a;
  if(rev) cc[len_b] = (score_t)-dist;
  else cc[0] = (score_t)-dist;

  for(j = 0; j < len_b; j++) {
    const unsigned char ch = (unsigned char)seq_b[rev ? (len_b - 1 - j) : j];
    const uint64_t *e = eq + table->code[ch] * nw;
    int hin = 1;
    for(w = 0; w < full_nw; w++)
      hin = myers_step_word(&P[w], &M[w], e[w], hin);
    hin = myers_step_last(&P[full_nw], &M[full_nw], e[full_nw],
                          hin, last_shift);
    dist += hin;
    if(rev) cc[len_b - 1 - j] = (score_t)-dist;
    else cc[j + 1] = (score_t)-dist;
  }

  free(eq);
  return true;
}
