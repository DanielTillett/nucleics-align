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
  bs->b = calloc(bs->s, sizeof(bs->b[0]));
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

struct sw_aligner_t {
  aligner_t aligner;
  sw_history_t history;
  sub_table_t table;
  bool fast_walk;
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
