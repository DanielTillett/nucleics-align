/*
 alignment_fill.c
 Full score-matrix fill (scalar & striped) and compact traceback fill.
 */

#define _DEFAULT_SOURCE // madvise

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdatomic.h>
#ifdef __linux__
#  include <sys/mman.h>
#endif

#include "alignment_fill_internal.h"
#include "alignment_macros.h"
#include "alignment_alloc.h"
#include "alignment_testing.h"

// Lowest score a cell may take (cells are clamped to it)
score_t alignment_min_score(const scoring_t *scoring, char is_sw)
{
  int gap_open_penalty = scoring->gap_extend + scoring->gap_open;
  int gap_extend_penalty = scoring->gap_extend;

  // Penalties may have been changed since scoring_init() (e.g. by the command
  // line parser), so take the most negative penalty from the current values
  int worst_penalty = MIN3(scoring->min_penalty, gap_open_penalty, gap_extend_penalty);
  worst_penalty = MIN3(worst_penalty, scoring->match, scoring->mismatch);

  return is_sw ? 0 : SCORE_MIN + abs(worst_penalty);
}

// Fill in traceback matrix. Returns false if out of memory
bool alignment_fill_matrices(aligner_t *aligner, char is_sw)
{
  score_t *match_scores = aligner->match_scores;
  score_t *gap_a_scores = aligner->gap_a_scores;
  score_t *gap_b_scores = aligner->gap_b_scores;
  const scoring_t *scoring = aligner->scoring;
  size_t score_width = aligner->score_width;
  size_t score_height = aligner->score_height;
  size_t i, j;

  int gap_open_penalty = scoring->gap_extend + scoring->gap_open;
  int gap_extend_penalty = scoring->gap_extend;

  const score_t min = alignment_min_score(scoring, is_sw);

  size_t seq_i, seq_j, len_i = score_width-1, len_j = score_height-1;

  bool no_gaps_in_a, no_gaps_in_b;
  alignment_no_gaps(scoring, len_i, len_j, &no_gaps_in_a, &no_gaps_in_b);
  size_t index, index_left, index_up, index_upleft;

  sub_table_t table;
  if(!sub_table_build(&table, scoring, aligner->seq_a, len_i,
                      aligner->seq_b, len_j)) return false;

  // [0][0]
  match_scores[0] = 0;
  gap_a_scores[0] = 0;
  gap_b_scores[0] = 0;

  if(is_sw)
  {
    for(i = 1; i < score_width; i++)
      match_scores[i] = gap_a_scores[i] = gap_b_scores[i] = 0;
    for(j = 1, index = score_width; j < score_height; j++, index += score_width)
      match_scores[index] = gap_a_scores[index] = gap_b_scores[index] = min;
  }
  else
  {
    // The first row's gaps start the alignment, and end it too where seq_b
    // is empty (the first row is then the last); the first column's the
    // same with seq_a
    const bool free_row = scoring->no_start_gap_penalty ||
                          (scoring->no_end_gap_penalty && len_j == 0);
    const bool free_col = scoring->no_start_gap_penalty ||
                          (scoring->no_end_gap_penalty && len_i == 0);

    // work along first row -> [i][0]
    for(i = 1; i < score_width; i++)
    {
      match_scores[i] = min;

      // Think carefully about which way round these are
      gap_a_scores[i] = min;
      gap_b_scores[i] = free_row ? 0
                        : scoring->gap_open + (int)i * scoring->gap_extend;
    }

    // work down first column -> [0][j]
    for(j = 1, index = score_width; j < score_height; j++, index += score_width)
    {
      match_scores[index] = min;

      // Think carefully about which way round these are
      gap_a_scores[index] = free_col ? 0
                            : scoring->gap_open + (int)j * scoring->gap_extend;
      gap_b_scores[index] = min;
    }
  }

  #ifdef ALN_STRIPED
  // The striped fill handles the uniform recurrence only; options that
  // change it for particular cells use the scalar loop below (an empty
  // sequence has no cells past the edges just set: nothing to fill)
  if(aligner->opts.striped && len_i > 0 && len_j > 0 &&
     !scoring->no_end_gap_penalty && !scoring->no_mismatches &&
     !no_gaps_in_a && !no_gaps_in_b)
  {
    alignment_block_t tmpl = {
      .aligner = aligner, .table = &table, .min = min,
      .r0 = {match_scores, gap_a_scores, gap_b_scores},
      .lc = {match_scores, gap_a_scores, gap_b_scores},
      .lstride = score_width, .width = len_i, .y0 = 1, .y1 = len_j};
    bool filled = alignment_fill_blocks(&tmpl);
    sub_table_free(&table);
    return filled;
  }
  #endif

  // start at position [1][1]
  index_upleft = 0;
  index_up = 1;
  index_left = score_width;
  index = score_width+1;

  for(seq_j = 0; seq_j < len_j; seq_j++)
  {
    size_t row = table.code[(unsigned char)aligner->seq_b[seq_j]] * table.n;
    const int *sub_row = table.scores + row;
    const bool *match_row = table.is_match + row;

    for(seq_i = 0; seq_i < len_i; seq_i++)
    {
      // Update match_scores[i][j] with position [i-1][j-1]
      // substitution penalty
      unsigned char code_a = table.code[(unsigned char)aligner->seq_a[seq_i]];
      bool is_match = match_row[code_a];
      int substitution_penalty = sub_row[code_a];

      if(scoring->no_mismatches && !is_match)
      {
        match_scores[index] = min;
      }
      else
      {
        // substitution
        // 1) continue alignment
        // 2) close gap in seq_a
        // 3) close gap in seq_b
        // max(x+s, y+s, z+s) == max(x,y,z)+s, with one add instead of three
        match_scores[index]
          = MAX2(MAX3(match_scores[index_upleft],
                      gap_a_scores[index_upleft],
                      gap_b_scores[index_upleft]) + substitution_penalty,
                 min);
      }

      // (Cells hold at least min, which is SCORE_MIN plus the worst
      // penalty, so adding a penalty to one cannot overflow)

      // Update gap_a_scores[i][j] from position [i][j-1]
      if(seq_i == len_i-1 && scoring->no_end_gap_penalty)
      {
        gap_a_scores[index] = MAX3(match_scores[index_up],
                                   gap_a_scores[index_up],
                                   gap_b_scores[index_up]);
      }
      else if(!no_gaps_in_a || seq_i == len_i-1)
      {
        gap_a_scores[index]
          = MAX3(MAX2(match_scores[index_up],
                      gap_b_scores[index_up]) + gap_open_penalty,
                 gap_a_scores[index_up] + gap_extend_penalty,
                 min);
      }
      else
        gap_a_scores[index] = min;

      // Update gap_b_scores[i][j] from position [i-1][j]
      if(seq_j == len_j-1 && scoring->no_end_gap_penalty)
      {
        gap_b_scores[index] = MAX3(match_scores[index_left],
                                   gap_a_scores[index_left],
                                   gap_b_scores[index_left]);
      }
      else if(!no_gaps_in_b || seq_j == len_j-1)
      {
        gap_b_scores[index]
          = MAX3(MAX2(match_scores[index_left],
                      gap_a_scores[index_left]) + gap_open_penalty,
                 gap_b_scores[index_left] + gap_extend_penalty,
                 min);
      }
      else
        gap_b_scores[index] = min;

      index++;
      index_left++;
      index_up++;
      index_upleft++;
    }

    index++;
    index_left++;
    index_up++;
    index_upleft++;
  }

  sub_table_free(&table);
  return true;
}

// Allocate memory for one score matrix. On Linux, large matrices are 2 MB
// aligned and marked for transparent huge pages: with 4 KB pages a 5 kb x
// 5 kb alignment takes ~72,000 page faults and many TLB misses, with huge
// pages ~320 (whole program 0.14 -> 0.08 s; repeated alignments ~30% faster)
void* alignment_alloc_matrix(size_t mem)
{
#ifdef __linux__
  const size_t huge = (size_t)2 << 20;
  if(mem >= huge) {
    mem = (mem + huge - 1) / huge * huge; // aligned_alloc needs a multiple
    void *ptr = aln_aligned_alloc(huge, mem);
    if(ptr != NULL) madvise(ptr, mem, MADV_HUGEPAGE); // only a hint
    return ptr;
  }
#endif
  return aln_malloc(mem);
}

#ifdef ALN_STRIPED
// Room for the traceback bytes of `cells` cells. Returns false if out of
// memory
size_t alignment_grow_capacity(size_t cells)
{
  // Doubling a multi-gigabyte buffer (e.g. 2.5 GiB -> 4 GiB) would fail
  // where a little more room fits: a quarter extra past 1M cells
  return cells <= ((size_t)1 << 20) ? ROUNDUP2POW(cells)
         : (cells <= SIZE_MAX - cells / 4 ? cells + cells / 4 : cells);
}

bool alignment_trace_capacity(aligner_t *aligner, size_t cells)
{
  if(aligner->trace_capacity < cells) {
    free(aligner->trace);
    size_t cap = alignment_grow_capacity(cells);
    aligner->trace_capacity = cap;
    aligner->trace = alignment_alloc_matrix(aligner->trace_capacity);
    if(aligner->trace == NULL) {
      aligner->trace_capacity = 0;
      return false;
    }
  }
  return true;
}
#endif /* ALN_STRIPED */

// Compact mode: run the striped fill keeping one traceback byte per cell
// instead of the score matrices. Returns false (nothing done) where the
// striped fill does not apply, with free start gaps (the traceback changes
// penalties along the first row and column), or where the first row and
// column fall below min: working out the traceback bytes adds a penalty to
// each predecessor, which could then overflow (the fill adds it only to the
// best one). FILL_NOMEM: out of memory, with no trace
fill_result_t alignment_fill_trace(aligner_t *aligner, char is_sw)
{
#ifdef ALN_STRIPED
  const scoring_t *scoring = aligner->scoring;
  const size_t w = aligner->score_width, h = aligner->score_height;
  const size_t len_i = w-1, len_j = h-1;
  size_t i, m;

  bool no_gaps_in_a, no_gaps_in_b;
  alignment_no_gaps(scoring, len_i, len_j, &no_gaps_in_a, &no_gaps_in_b);
  if(!aligner->opts.striped || len_i == 0 || len_j == 0 ||
     len_i >= INT_MAX || len_j >= INT_MAX ||
     w + h > SIZE_MAX / (3 * sizeof(score_t)) ||
     scoring->no_start_gap_penalty ||
     scoring->no_end_gap_penalty || scoring->no_mismatches ||
     no_gaps_in_a || no_gaps_in_b) return FILL_SKIPPED;

  const score_t min = alignment_min_score(scoring, is_sw);
  if(!is_sw) {
    // Gap scores along the first row and column run from these two
    long long edge1 = (long long)scoring->gap_open + scoring->gap_extend;
    long long edgen = (long long)scoring->gap_open +
                      (long long)MAX2(len_i, len_j) * scoring->gap_extend;
    if(MIN2(edge1, edgen) < min || MAX2(edge1, edgen) > INT_MAX)
      return FILL_SKIPPED;
  }

  // Row 0 and column 0, as alignment_fill_matrices() sets them
  score_t *edges = aln_malloc(3 * (w + h) * sizeof(score_t));
  if(edges == NULL) return FILL_NOMEM;
  score_t *row0[3], *col0[3];
  for(m = 0; m < 3; m++) {
    row0[m] = edges + m * w;
    col0[m] = edges + 3 * w + m * h;
    row0[m][0] = col0[m][0] = 0;
  }
  for(i = 1; i < w; i++) {
    row0[MATCH][i] = row0[GAP_A][i] = is_sw ? 0 : min;
    row0[GAP_B][i] = is_sw ? 0 : scoring->gap_open + (int)i * scoring->gap_extend;
  }
  for(i = 1; i < h; i++) {
    col0[MATCH][i] = col0[GAP_B][i] = min;
    col0[GAP_A][i] = is_sw ? min : scoring->gap_open + (int)i * scoring->gap_extend;
  }

  sub_table_t table;
  if(!sub_table_build(&table, scoring, aligner->seq_a, len_i,
                      aligner->seq_b, len_j)) {
    free(edges);
    return FILL_NOMEM;
  }
  aligner->trace_banded = false;
  aligner->trace_striped = false;
  trace_strips_free(aligner->trace_strips);
  aligner->trace_strips = NULL;
  fill_result_t done = FILL_SKIPPED;
  if(!is_sw) done = alignment_fill_banded(aligner, &table, min, row0, col0);
  if(done == FILL_SKIPPED)
  {
    score_t dfloor = 0, dceiling = 0;
    int diff_bits = (!is_sw && aligner->opts.diff)
                    ? alignment_diff_bits(&table, scoring, &dfloor, &dceiling)
                    : 0;
    size_t trace_cells = w * h;
    size_t trace_stride = 0;
    if(diff_bits != 0) {
      const size_t vl = (diff_bits == 8) ? (size_t)(4 * ALN_VL)
                                         : (size_t)(2 * ALN_VL);
      const size_t nb = alignment_fill_threads(&aligner->opts, len_i, len_j, true);
      const size_t bw = (len_i / nb) / ALN_VL * ALN_VL;
      const size_t last_w = len_i - (nb - 1) * bw;
      const size_t seg0 = (bw + vl - 1) / vl;
      const size_t seg_last = (last_w + vl - 1) / vl;
      trace_stride = (nb - 1) * seg0 * vl + seg_last * vl;
      if(len_j > SIZE_MAX / trace_stride) done = FILL_NOMEM;
      else {
        trace_cells = len_j * trace_stride;
        aligner->trace_striped = true;
        aligner->trace_nblocks = nb;
        aligner->trace_bwidth = bw;
        aligner->trace_seg0 = seg0;
        aligner->trace_seg_last = seg_last;
        aligner->trace_vl = vl;
        aligner->trace_stride = trace_stride;
      }
    }
    if(done == FILL_SKIPPED && !alignment_trace_capacity(aligner, trace_cells)) {
      aligner->trace_striped = false;
      done = FILL_NOMEM;
    }
    else if(done == FILL_SKIPPED) {
      alignment_block_t tmpl = {
        .aligner = aligner, .table = &table, .min = min,
        .r0 = {row0[MATCH], row0[GAP_A], row0[GAP_B]},
        .lc = {col0[MATCH], col0[GAP_A], col0[GAP_B]}, .lstride = 1,
        .trace = aligner->trace, .track_best = is_sw,
        .width = len_i, .y0 = 1, .y1 = len_j,
        .diff = diff_bits, .floor = dfloor, .ceiling = dceiling,
        .trace_stride = trace_stride};
      if(is_sw) alignment_lanes(&tmpl, &aligner->opts, scoring, len_i);
      done = alignment_fill_blocks(&tmpl) ? FILL_DONE : FILL_NOMEM;
      if(done != FILL_DONE) aligner->trace_striped = false;
    }
  }
  sub_table_free(&table);
  free(edges);
  return done;
#else
  (void)aligner; (void)is_sw;
  return FILL_SKIPPED;
#endif
}
