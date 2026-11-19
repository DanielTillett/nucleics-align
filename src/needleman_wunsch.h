/*
 needleman_wunsch.h
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */


#ifndef NEEDLEMAN_WUNSCH_HEADER_SEEN
#define NEEDLEMAN_WUNSCH_HEADER_SEEN

#include "seq_align.h"
#include "alignment.h"

typedef aligner_t nw_aligner_t;

#ifdef __cplusplus
extern "C" {
#endif

nw_aligner_t* needleman_wunsch_new(void);
void needleman_wunsch_free(nw_aligner_t *nw);

void needleman_wunsch_align(const char *a, const char *b,
                            const scoring_t *scoring,
                            nw_aligner_t *nw, alignment_t *result);

void needleman_wunsch_align2(const char *a, const char *b,
                             size_t len_a, size_t len_b,
                             const scoring_t *scoring,
                             nw_aligner_t *nw, alignment_t *result);

// As needleman_wunsch_align2(), but an error in the input (see aln_status_t)
// is returned, with err set if not NULL, instead of ending the program
aln_status_t needleman_wunsch_align_checked(const char *a, const char *b,
                                            size_t len_a, size_t len_b,
                                            const scoring_t *scoring,
                                            nw_aligner_t *nw,
                                            alignment_t *result,
                                            aln_error_t *err);

#ifdef __cplusplus
}
#endif

#endif
