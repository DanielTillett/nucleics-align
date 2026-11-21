/*
 alignment_scoring_load.h
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */

#ifndef ALIGNMENT_SCORING_LOAD_HEADER_SEEN
#define ALIGNMENT_SCORING_LOAD_HEADER_SEEN

#include <zlib.h>
#include "alignment_scoring.h"

#ifdef __cplusplus
extern "C" {
#endif

void align_scoring_load_matrix(gzFile file, const char* file_path,
                               scoring_t* scoring, char case_sensitive);

void align_scoring_load_pairwise(gzFile file, const char* file_path,
                                 scoring_t* scoring, char case_sensitive);

#ifdef __cplusplus
}
#endif

#endif
