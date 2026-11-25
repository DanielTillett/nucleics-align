/*
 smith_waterman.h
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */

#ifndef SMITH_WATERMAN_HEADER_SEEN
#define SMITH_WATERMAN_HEADER_SEEN

#include "nucleics_align.h"
#include "alignment.h"

typedef struct sw_aligner_t sw_aligner_t;

#ifdef __cplusplus
extern "C" {
#endif

sw_aligner_t *smith_waterman_new(void);
sw_aligner_t *smith_waterman_try_new(void);
void smith_waterman_free(sw_aligner_t *sw_aligner);

aligner_t* smith_waterman_get_aligner(sw_aligner_t *sw);

/*
 Do not alter seq_a, seq_b or scoring whilst calling this method
 or between calls to smith_waterman_get_hit
*/
void smith_waterman_align(const char *seq_a, const char *seq_b,
                          const scoring_t *scoring, sw_aligner_t *sw);

void smith_waterman_align2(const char *seq_a, const char *seq_b,
                           size_t len_a, size_t len_b,
                           const scoring_t *scoring, sw_aligner_t *sw);

// As smith_waterman_align2(), but an error (see aln_status_t: in the input,
// or running out of memory) is returned, with err set if not NULL, instead of
// ending the program. After an error there are no hits to fetch
aln_status_t smith_waterman_align_checked(const char *seq_a, const char *seq_b,
                                          size_t len_a, size_t len_b,
                                          const scoring_t *scoring,
                                          sw_aligner_t *sw, aln_error_t *err);

// Only consider local alignments scoring at least min_score (default 1) in
// later calls to smith_waterman_align(). Hits come out best first, so the
// hits returned are exactly those of a full search up to the first one
// scoring below min_score; lower-scoring cells are just never indexed
void smith_waterman_set_min_score(sw_aligner_t *sw, score_t min_score);

// For callers that usually want just the best hit (default false): later
// calls to smith_waterman_align() keep one traceback byte per cell instead
// of the three score matrices, where the scoring options allow (see
// aligner_t.compact). The hits are the same. Fetching a second hit fills the
// matrices then; until it does, the aligner from smith_waterman_get_aligner()
// has no score matrices if its traced flag is set
void smith_waterman_set_compact(sw_aligner_t *sw, bool compact);

// Settings for later alignments (opts must stay valid while they run; NULL:
// the defaults, see alignment_options_t). For Needleman-Wunsch set the
// aligner's options field
void smith_waterman_set_options(sw_aligner_t *sw,
                                const alignment_options_t *opts);

// Compact alignments of at least this many cells (default 2^24, 16 MB of
// traceback bytes) find the best hit with O(len_b)-memory scans and then
// fill just the rectangle of the best hits: much less memory and time when
// the hit is short against the sequences (50 kb pairs sharing 2 kb: 0.55 ->
// 0.09 s, 2.9 GB -> 10 MB; 20 kb unrelated: 0.08 -> 0.01 s, 450 -> 5 MB).
// The first scan gives up once the hit is seen to span the sequences (see
// alignment_sw_scan_short()), so those cost about the same as before. Same
// hits. Not used where a gap scores above 0 (a hit can then start with
// one). Until a second hit is fetched, the aligner describes the rectangle,
// or the whole sequences if no hit reaches the minimum score. For tuning and
// tests
extern size_t smith_waterman_rect_min_cells;

// Compact best hits whose region (the rectangle above, or the whole matrix
// where there is none) has at least this many cells (default 2^28, 256 MB
// of traceback bytes) keep no traceback bytes for it: the scan for the best
// cell also keeps every k-th row of the matrices (k about sqrt(12 x rows)),
// and the hit is followed back strip by strip, each strip's bytes filled
// again from the row kept above it, only up to the hit's column. The bytes,
// and so the hit, are those of the whole fill, in O(columns x sqrt(rows))
// memory (50 kb pairs sharing the whole sequence: 2.9 GB -> MB). A second
// hit still fills the score matrices. For tuning and tests
extern size_t smith_waterman_strip_min_cells;

// Opt-in fast best-hit reconstruction (default false): after forward and
// reverse O(len_b)-memory scans locate an optimal local start and end pair,
// reconstruct the sub-rectangle with BiWFA or linear-space global alignment
// and verify exact optimal local score and strictly positive prefixes
// (falling back to the exact local fill on any mismatch). Preserves optimal
// local score; among equal-scoring paths, may choose another alignment
extern bool smith_waterman_fast_rect;

// An alignment to read from, and a pointer to memory to store the result
// returns 1 if an alignment was read, 0 otherwise
int smith_waterman_fetch(sw_aligner_t *sw, alignment_t *result);

// As smith_waterman_fetch(), but if out of memory returns 0 with
// err->status ALN_ERR_NOMEM (else ALN_OK) instead of ending the program.
// Fetching a second hit fills the score matrices if a compact alignment
// has not. After an error there are no more hits to fetch
int smith_waterman_fetch_checked(sw_aligner_t *sw, alignment_t *result,
                                 aln_error_t *err);

#ifdef __cplusplus
}
#endif

#endif
