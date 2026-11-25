/*
 alignment_fill_internal.h
 Internal cross-module seam for alignment.c, alignment_kernel.c,
 alignment_fill.c, alignment_band.c, and alignment_debug.c.
 */

#ifndef ALIGNMENT_FILL_INTERNAL_HEADER_SEEN
#define ALIGNMENT_FILL_INTERNAL_HEADER_SEEN

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <limits.h>
#include <pthread.h>
#include <stdatomic.h>

#include "alignment.h"
#include "alignment_simd.h"

// Which sequences may not contain gaps once no_gaps_in_longest is resolved
static inline void alignment_no_gaps(const scoring_t *scoring,
                                     size_t len_a, size_t len_b,
                                     bool *no_gaps_in_a, bool *no_gaps_in_b)
{
  *no_gaps_in_a = scoring->no_gaps_in_a ||
                  (scoring->no_gaps_in_longest && len_a >= len_b);
  *no_gaps_in_b = scoring->no_gaps_in_b ||
                  (scoring->no_gaps_in_longest && len_a < len_b);
}

// Whether the striped kernels may fill len_i x len_j cells under scoring:
// they handle the uniform recurrence only, so the options that change it
// for particular cells (no gaps or mismatches, free end gaps) take the
// scalar fill. edges: the kernels are given the first row and column
// (traceback bytes and scans), so free start gaps, which change those, are
// refused too, and the lengths must fit the kernels' counters and the edge
// buffers; the matrix fill sets its own edges
static inline bool alignment_striped_applies(const alignment_options_t *o,
                                             const scoring_t *scoring,
                                             size_t len_i, size_t len_j,
                                             bool edges)
{
  bool no_gaps_in_a, no_gaps_in_b;
  alignment_no_gaps(scoring, len_i, len_j, &no_gaps_in_a, &no_gaps_in_b);
  if(!o->striped || len_i == 0 || len_j == 0 ||
     scoring->no_end_gap_penalty || scoring->no_mismatches ||
     no_gaps_in_a || no_gaps_in_b) return false;
  return !edges ||
         (!scoring->no_start_gap_penalty &&
          len_i < INT_MAX && len_j < INT_MAX &&
          len_i + len_j + 2 <= SIZE_MAX / (3 * sizeof(score_t)));
}

// What a compact fill did
typedef enum { FILL_SKIPPED, FILL_DONE, FILL_NOMEM } fill_result_t;

// Internal diagnostic counters (defined in alignment.c)
extern atomic_size_t aln_stat_narrow_fills, aln_stat_narrow_retries;
extern atomic_size_t aln_stat_byte_fills, aln_stat_byte_retries;
extern atomic_size_t aln_stat_diff_fills;
extern atomic_size_t aln_stat_band_tries, aln_stat_band_wfa_sized;

#ifdef ALN_STRIPED

typedef struct {
  _Alignas(64) atomic_size_t v;
  atomic_uint waiters;
  pthread_mutex_t mu;
  pthread_cond_t cv;
} aln_done_t;

// One column block of the striped fill (see alignment_kernel.c): rows
// y0 .. y1 of columns c0+1 .. c0+width. Which fields are set says what it
// produces (the kernels in alignment_kernel.inc read them so):
//   full matrices  trace == NULL, !score_only: the three score matrices,
//                  r0/lc the matrices themselves (lstride the row width)
//   traceback      trace != NULL: a byte per cell (trace_rows: a band's
//                  strip, by row; else trace_boff into the striped layout
//                  when diff != 0), the last column in rc (from row rc_y0)
//                  for the block on the right, and with track_best the
//                  first best cell; lastrow for the strip below
//   score scan     score_only: rc and the best cell (track_best), or the
//                  cells scoring target (target > 0: tgt_x, tgt_y)
// bits (0, 16 or 8) picks the lane width, narrow lanes stopping every
// block (*stop) past narrow_limit; diff != 0 the difference-encoded
// kernel with its floor/ceiling. given_up, best_*, tgt_* and nomem come
// back out
typedef struct {
  aligner_t *aligner; const sub_table_t *table; score_t min;
  size_t c0, width; aln_done_t *left_done, *done;
  const score_t *r0[3], *lc[3]; size_t lstride;
  uint8_t *trace; score_t *rc[3];
  bool track_best; score_t best_score; size_t best_x, best_y;
  size_t y0, y1, ly0;
  size_t rc_y0; // the row rc[c][0] holds (rc runs from row y0-1 to y1)
  const ptrdiff_t *trace_rows;
  score_t *lastrow[3];
  bool score_only;
  score_t target;
  size_t tgt_x, tgt_y;
  double giveup_rate;
  atomic_bool *gave_up;
  bool band;
  bool given_up;
  score_t *ckpt;
  size_t ckpt_every;
  int bits;
  score_t narrow_limit;
  int diff;
  score_t floor, ceiling;
  size_t trace_stride, trace_boff;
  atomic_bool *stop, *nomem;
} alignment_block_t;

// Defined in alignment_kernel.c
int alignment_diff_bits(const sub_table_t *table,
                        const scoring_t *scoring, score_t *floor,
                        score_t *ceiling);
bool alignment_narrow_limit(const sub_table_t *table,
                            const scoring_t *scoring,
                            score_t min_score, size_t len_a,
                            int bits, score_t *limit);
int alignment_lane_bits(const sub_table_t *table,
                        const alignment_options_t *opts,
                        const scoring_t *scoring, score_t min, size_t len_i,
                        score_t need, bool allow8, score_t *limit);
void alignment_lanes(alignment_block_t *tmpl,
                     const alignment_options_t *o,
                     const scoring_t *scoring, size_t len_a);
// whole: a fill across the whole width (false: one strip of a band)
size_t alignment_fill_threads(const alignment_options_t *o,
                              size_t width, size_t rows, bool whole);
bool alignment_fill_blocks(alignment_block_t *tmpl);

// Defined in alignment_fill.c
bool alignment_trace_capacity(aligner_t *aligner, size_t cells);

// Defined in alignment_band.c
void trace_strips_free(void *ts_ptr);
fill_result_t alignment_fill_banded(aligner_t *aligner,
                                    const sub_table_t *table,
                                    score_t min,
                                    score_t *const r0[3],
                                    score_t *const c0[3]);

#endif /* ALN_STRIPED */

// Always declared (stubbed when ALN_STRIPED is not defined)
void trace_strips_free(void *ts_ptr);

// Defined in alignment_fill.c
score_t alignment_min_score(const scoring_t *scoring, char is_sw);
void* alignment_alloc_matrix(size_t mem);
// Room for cells with some to spare: the next power of two up to 1M cells,
// then a quarter more (doubling a multi-gigabyte buffer would fail where
// a little more fits). Used by every build: alignment.c grows the score
// matrices with it
size_t alignment_grow_capacity(size_t cells);
bool alignment_fill_matrices(aligner_t *aligner, char is_sw);
fill_result_t alignment_fill_trace(aligner_t *aligner, char is_sw);

#endif /* ALIGNMENT_FILL_INTERNAL_HEADER_SEEN */
