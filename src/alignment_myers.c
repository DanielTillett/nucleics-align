/*
 alignment_myers.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// The score of the best global alignment under unit costs (a match 0, a
// mismatch or a gap letter -1, no gap open penalty), which is minus the
// edit distance, by Myers' bit-vector algorithm (Myers 1999, J ACM
// 46:395-415, doi:10.1145/316542.316550), written from the paper. A column
// of the dynamic programming matrix is held as two bit-vectors, the rows
// where a cell is one more (P) or one less (M) than the cell above, and
// is advanced a 64-row word at a time with a handful of word operations,
// the words of a column chained by the horizontal difference each passes
// down. That beats the striped passes on one thread for divergent pairs
// (100 kb unrelated: 1.44 -> 0.61 s), which is where it is used (see
// hb_score()); banded passes and several threads beat it otherwise.

#include <stdlib.h>
#include <stdint.h>

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

bool alignment_myers_score(const sub_table_t *table, const scoring_t *scoring,
                           const char *seq_a, size_t len_a,
                           const char *seq_b, size_t len_b, score_t *score)
{
  if(!myers_unit_cost(table, scoring) || len_a == 0 ||
     len_a + len_b > (size_t)SCORE_ABS_MAX) return false;

  // Rows run down seq_a, 64 to a word; columns along seq_b. eq[c*nw + w]
  // has a bit for each row of word w whose letter scores 0 against a
  // letter of seq_b with code c
  const size_t n = table->n, nw = (len_a + 63) / 64;
  if(nw > SIZE_MAX / sizeof(uint64_t) / (n + 2)) return false;
  uint64_t *eq = aln_calloc((n + 2) * nw, sizeof(uint64_t));
  if(eq == NULL) return false;
  uint64_t *P = eq + n * nw, *M = P + nw;
  size_t i, j, c, w;

  for(i = 0; i < len_a; i++) {
    const size_t ca = table->code[(unsigned char)seq_a[i]];
    for(c = 0; c < n; c++)
      if(table->scores[c*n + ca] == 0)
        eq[c*nw + i/64] |= (uint64_t)1 << (i % 64);
  }

  // Column 0: each row one more than the one above (a gap down seq_a)
  for(w = 0; w < nw; w++) P[w] = ~(uint64_t)0;
  const uint64_t last = (uint64_t)1 << ((len_a - 1) % 64);
  long long dist = (long long)len_a;

  for(j = 0; j < len_b; j++) {
    const uint64_t *e = eq + table->code[(unsigned char)seq_b[j]] * nw;
    int hin = 1; // row 0 grows by one per column (a gap along seq_b)
    for(w = 0; w < nw; w++) {
      uint64_t Pv = P[w], Mv = M[w], Eq = e[w];
      const uint64_t Xv = Eq | Mv;
      if(hin < 0) Eq |= 1;
      const uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
      uint64_t Ph = Mv | ~(Xh | Pv);
      uint64_t Mh = Pv & Xh;
      // The difference leaving the word's last row (the last row of
      // seq_a in the last word) goes down to the next word
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
    dist += hin; // the bottom row's step to this column
  }

  free(eq);
  *score = (score_t)-dist;
  return true;
}
