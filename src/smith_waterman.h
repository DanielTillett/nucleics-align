/*
 smith_waterman.h
 url: https://github.com/DanielTillett/seq-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */

#ifndef SMITH_WATERMAN_HEADER_SEEN
#define SMITH_WATERMAN_HEADER_SEEN

#include "seq_align.h"
#include "alignment.h"

typedef struct sw_aligner_t sw_aligner_t;

#ifdef __cplusplus
extern "C" {
#endif

sw_aligner_t *smith_waterman_new(void);
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

// Compact alignments of at least this many cells (default 2^30, 1 GB of
// traceback bytes) find the best hit with O(len_b)-memory scans and then
// fill just the rectangle of the best hits: much less memory and time when
// the hit is short against the sequences (50 kb pairs sharing 2 kb: 0.55 ->
// 0.09 s, 2.9 GB -> 10 MB), but about 50% more time when it spans them.
// Same hits. Not used where a gap scores above 0 (a hit can then start with
// one). Until a second hit is fetched, the aligner describes the rectangle,
// or the whole sequences if no hit reaches the minimum score. For tuning and
// tests
extern size_t smith_waterman_rect_min_cells;

// An alignment to read from, and a pointer to memory to store the result
// returns 1 if an alignment was read, 0 otherwise
int smith_waterman_fetch(sw_aligner_t *sw, alignment_t *result);

#ifdef __cplusplus
}
#endif

#endif
