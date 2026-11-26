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

// Add the substitution matrix, or pairs, in file to scoring (whose
// case_sensitive setting decides how letters are read), or print the
// error with file_path and exit. See scoring_parse_matrix()
void align_scoring_load_matrix(gzFile file, const char* file_path,
                               scoring_t* scoring);

void align_scoring_load_pairwise(gzFile file, const char* file_path,
                                 scoring_t* scoring);

#ifdef __cplusplus
}
#endif

#endif
