/*
 hirschberg.h
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

#ifndef HIRSCHBERG_HEADER_SEEN
#define HIRSCHBERG_HEADER_SEEN

#include <stdbool.h>

#include "nucleics_align.h"
#include "alignment.h"

#ifdef __cplusplus
extern "C" {
#endif

// Global alignment in O(len_b) memory: Hirschberg's divide and conquer with
// the affine gap extension of Myers & Miller (1988) Comput Appl Biosci 4:11-17.
// Scores are the same as needleman_wunsch_align() but, where several
// alignments share the optimal score, the one returned may differ.
// Supports substitution tables, wildcards and affine gaps. Returns false
// (and leaves result untouched) if gap_open > 0 or any of
// no_start_gap_penalty, no_end_gap_penalty, no_gaps_in_a/b/longest or
// no_mismatches is set, or if the penalties are so large that a score could
// overflow (roughly 3*|gap_open| + (len_a+len_b)*|gap_extend| above 1e9).
bool hirschberg_align(const char *a, const char *b,
                      const scoring_t *scoring, alignment_t *result);

bool hirschberg_align2(const char *a, const char *b,
                       size_t len_a, size_t len_b,
                       const scoring_t *scoring, alignment_t *result);

// With settings of its own (opts; NULL: the defaults, see alignment_options_t)
bool hirschberg_align_opts(const char *a, const char *b,
                           size_t len_a, size_t len_b,
                           const scoring_t *scoring,
                           const alignment_options_t *opts,
                           alignment_t *result);

// Tuning (and testing) of column-block threading in large passes: a pass of
// at least hirschberg_block_min_cells cells is split into blocks of at least
// hirschberg_block_min_cols columns, one thread each, up to
// hirschberg_max_threads in all (0: half the CPUs); with two or more, big
// forward and reverse passes run at the same time with half each. Set
// before aligning
extern size_t hirschberg_block_min_cells, hirschberg_block_min_cols;
extern int hirschberg_max_threads;

// Banded passes: cells that cannot be on a best path are skipped, with the
// same result. hirschberg_band_enabled (default true) turns them off;
// hirschberg_band_strip_rows (0: automatic) and hirschberg_band_min_width
// (first band tried, default 64) are for tuning and tests
extern bool hirschberg_band_enabled;
extern size_t hirschberg_band_strip_rows, hirschberg_band_min_width;

// Rows reused down the recursion: a pass over all columns also saves the
// rows the next hirschberg_reuse_depth (default 2; 0: none) sub-problems on
// its side will need, if they have at least hirschberg_reuse_min_rows rows
// (default 32), so they skip that pass. Each level can hold up to four more
// rows of scores (the passes use four). The result is the same
extern size_t hirschberg_reuse_depth, hirschberg_reuse_min_rows;

// Whether this build and CPU run the SIMD passes (and so banding) under
// these options (NULL: the defaults)
bool hirschberg_simd_available(const alignment_options_t *options);

// Score of the best global alignment, without the alignment itself: one
// pass in O(len_b) memory. Same score as needleman_wunsch_align(), and the
// same unsupported options as hirschberg_align() (returns false)
bool hirschberg_score(const char *a, const char *b,
                      const scoring_t *scoring, score_t *score);

bool hirschberg_score2(const char *a, const char *b,
                       size_t len_a, size_t len_b,
                       const scoring_t *scoring, score_t *score);

bool hirschberg_score_opts(const char *a, const char *b,
                           size_t len_a, size_t len_b,
                           const scoring_t *scoring,
                           const alignment_options_t *opts, score_t *score);

// As the *_opts() functions, but an error (see aln_status_t: in the input,
// ALN_ERR_UNSUPPORTED where they return false, or ALN_ERR_NOMEM, running out
// of memory) is returned, with err set if not NULL, instead of ending the
// program
aln_status_t hirschberg_align_checked(const char *a, const char *b,
                                      size_t len_a, size_t len_b,
                                      const scoring_t *scoring,
                                      const alignment_options_t *opts,
                                      alignment_t *result, aln_error_t *err);

aln_status_t hirschberg_score_checked(const char *a, const char *b,
                                      size_t len_a, size_t len_b,
                                      const scoring_t *scoring,
                                      const alignment_options_t *opts,
                                      score_t *score, aln_error_t *err);

// Score an alignment as needleman_wunsch_align() would without the free
// start/end gap options (every gap is charged; those options are not
// applied). A '-' in either string is read as a gap and the strings end at
// '\0': hirschberg_align() scores its alignments from the steps taken
// instead
score_t hirschberg_score_alignment(const char *result_a, const char *result_b,
                                   const scoring_t *scoring);

#ifdef __cplusplus
}
#endif

#endif
