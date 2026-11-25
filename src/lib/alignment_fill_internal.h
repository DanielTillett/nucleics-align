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

// One column block of the striped fill (see alignment_kernel.c)
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
void alignment_lanes(alignment_block_t *tmpl,
                     const alignment_options_t *o,
                     const scoring_t *scoring, size_t len_a);
size_t alignment_fill_threads(const alignment_options_t *o,
                              size_t len_i, size_t len_j, bool band);
bool alignment_fill_blocks(alignment_block_t *tmpl);

// Defined in alignment_fill.c
bool alignment_trace_capacity(aligner_t *aligner, size_t cells);
// Room for cells with some to spare: the next power of two up to 1M cells,
// then a quarter more (doubling a multi-gigabyte buffer would fail where
// a little more fits)
size_t alignment_grow_capacity(size_t cells);

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
bool alignment_fill_matrices(aligner_t *aligner, char is_sw);
fill_result_t alignment_fill_trace(aligner_t *aligner, char is_sw);

#endif /* ALIGNMENT_FILL_INTERNAL_HEADER_SEEN */
