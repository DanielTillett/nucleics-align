/*
 hirschberg.h
 url: https://github.com/DanielTillett/seq-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

#ifndef HIRSCHBERG_HEADER_SEEN
#define HIRSCHBERG_HEADER_SEEN

#include <stdbool.h>

#include "seq_align.h"
#include "alignment.h"

#ifdef __cplusplus
extern "C" {
#endif

// Global alignment in O(len_b) memory: Hirschberg's divide and conquer with
// the affine gap extension of Myers & Miller (1988) Comput Appl Biosci 4:11-17.
// Scores are the same as needleman_wunsch_align() but, where several
// alignments share the optimal score, the one returned may differ.
// Supports substitution tables, wildcards and affine gaps. Returns false
// (and leaves result untouched) if any of no_start_gap_penalty,
// no_end_gap_penalty, no_gaps_in_a/b/longest or no_mismatches is set.
bool hirschberg_align(const char *a, const char *b,
                      const scoring_t *scoring, alignment_t *result);

bool hirschberg_align2(const char *a, const char *b,
                       size_t len_a, size_t len_b,
                       const scoring_t *scoring, alignment_t *result);

// Score an alignment as needleman_wunsch_align() would
score_t hirschberg_score_alignment(const char *result_a, const char *result_b,
                                   const scoring_t *scoring);

#ifdef __cplusplus
}
#endif

#endif
