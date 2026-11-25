/*
 smith_waterman_internal.h
 Internal cross-module seam for smith_waterman.c and smith_waterman_rect.c.
 */

#ifndef SMITH_WATERMAN_INTERNAL_HEADER_SEEN
#define SMITH_WATERMAN_INTERNAL_HEADER_SEEN

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "smith_waterman.h"
#include "hirschberg.h"
#include "alignment_macros.h"
#include "alignment_alloc.h"
#include "alignment_testing.h"

typedef struct {
  uint32_t *b; size_t l, s; // l is bits, s in uint32_t
} BitSet;

static inline BitSet* bitset_alloc(BitSet *bs, size_t l) {
  bs->s = l/32 + (l%32 != 0);
  bs->l = l;
  bs->b = aln_calloc(bs->s, sizeof(bs->b[0]));
  return bs->b ? bs : NULL;
}

static inline void bitset_dealloc(BitSet *bs) {
  free(bs->b);
  memset(bs, 0, sizeof(*bs));
}

static inline BitSet* bitset_set_length(BitSet *bs, size_t l) {
  size_t ss = l/32 + (l%32 != 0);
  if(ss > bs->s) {
    uint32_t *b = aln_realloc(bs->b, ss*sizeof(bs->b[0]));
    if(b == NULL) return NULL;
    memset(b+bs->s, 0, (ss-bs->s)*sizeof(b[0]));
    bs->b = b;
    bs->s = ss;
  }
  bs->l = l;
  return bs;
}

typedef struct {
  BitSet match_scores_mask;
  size_t *sorted_match_indices, hits_capacity, num_of_hits, next_hit;
  bool indexed, first_done;
} sw_history_t;

// How the last alignment was made, by the flags _sw_align() leaves set
// (_sw_fetch() reads them to find the first hit and, for a second, to fill
// the score matrices):
//   full matrices    !compact (or compact but the scoring needs the scalar
//                    fill): aligner.match_scores hold every cell
//   traceback bytes  compact, !rect: aligner.trace holds a byte per cell,
//                    aligner.traced; best cell at aligner.best_index
//   rectangle        rect, !strips, !fast_hit: the bytes of the rectangle
//                    at (off_a, off_b) holding the best hit, whose cell in
//                    the whole matrix is full_best; or no hit at all
//                    (aligner.best_score < min_score), nothing filled
//   strips           rect && strips: ckpt holds every strip_rows-th row of
//                    the rectangle; sw_follow_strips() refills strips
//   fast hit         rect && fast_hit: the first hit is already in
//                    fast_res (wavefront or linear-memory reconstruction)
// In every mode but the first, aligner.traced is set, whether or not
// aligner.trace was filled
struct sw_aligner_t {
  aligner_t aligner;
  sw_history_t history;
  sub_table_t table;
  bool fast_walk; // _sw_reverse_move() may take its quick step
  score_t min_score;
  bool compact;
  bool rect;
  size_t off_a, off_b, full_best;
  const char *seq_a, *seq_b;
  size_t len_a, len_b;
  bool strips;
  score_t *ckpt;
  size_t strip_rows;
  bool fast_hit;
  alignment_t fast_res;
};

// One step back along traceback bytes, shared by the walks over a whole
// traced matrix (_sw_reverse_move()) and over a strip (sw_follow_strips()):
// from cell (x, y) in matrix *m, whose byte says the step goes to matrix
// next, to the predecessor, taking the step's score off *score. Exits on
// ALN_TRACE_NONE (a byte with no predecessor: a program error), naming
// where
static inline void sw_trace_step(const sub_table_t *t,
                                 const scoring_t *scoring,
                                 const char *a, const char *b, int next,
                                 enum Matrix *m, size_t *x, size_t *y,
                                 score_t *score, const char *where)
{
  const score_t open = scoring->gap_open + scoring->gap_extend;
  const score_t ext = scoring->gap_extend;
  score_t step;
  if(next == ALN_TRACE_NONE) {
    fprintf(stderr, "Program error: traceback fail (%s)\n", where);
    exit(EXIT_FAILURE);
  }
  if(*m == MATCH)
    step = t->scores[t->code[(unsigned char)b[*y-1]] * t->n +
                     t->code[(unsigned char)a[*x-1]]];
  else if(*m == GAP_A) step = next == GAP_A ? ext : open;
  else step = next == GAP_B ? ext : open;
  if(*m != GAP_A) (*x)--;
  if(*m != GAP_B) (*y)--;
  *m = (enum Matrix)next;
  *score -= step;
}

// Describe the sequences as a fill of them would, for a result found
// without one: the best cell (best_index in a matrix of (len_a + 1)
// columns) scoring best, found in the rectangle at (off_a, off_b)
static inline void sw_describe(sw_aligner_t *sw, const scoring_t *scoring,
                               const char *a, const char *b,
                               size_t len_a, size_t len_b,
                               score_t best, size_t best_index,
                               size_t off_a, size_t off_b)
{
  aligner_t *aligner = &sw->aligner;
  aligner->scoring = scoring;
  aligner->seq_a = a;
  aligner->seq_b = b;
  aligner->score_width = len_a + 1;
  aligner->score_height = len_b + 1;
  aligner->traced = true;
  aligner->trace_banded = false;
  aligner->best_score = best;
  aligner->best_index = best_index;
  sw->rect = true;
  sw->off_a = off_a;
  sw->off_b = off_b;
}

// Defined in smith_waterman_rect.c
aln_status_t sw_strips_fill(sw_aligner_t *sw, size_t x0, size_t y0,
                            size_t rx, size_t ry, const scoring_t *scoring,
                            score_t best, bool *used);
int sw_follow_strips(sw_aligner_t *sw, alignment_t *result);
bool sw_align_fast_rect(sw_aligner_t *sw, const char *a, const char *b,
                        size_t len_a, size_t len_b, const scoring_t *scoring);
bool sw_align_rect(sw_aligner_t *sw, const char *a, const char *b,
                   size_t len_a, size_t len_b, const scoring_t *scoring);

#endif /* SMITH_WATERMAN_INTERNAL_HEADER_SEEN */
