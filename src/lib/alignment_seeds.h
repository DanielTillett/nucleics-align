/*
 alignment_seeds.h
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// Private to the library: a bound on the score the rest of an alignment
// can add from a cell, from seed words of a and their copies in b (see
// alignment_seeds.c), used to narrow band strips (alignment_band.c)

#ifndef ALIGNMENT_SEEDS_HEADER_SEEN
#define ALIGNMENT_SEEDS_HEADER_SEEN

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include "alignment.h"

typedef struct seeds_t seeds_t;

// Seeds of a with their copies in b, for a scoring where letters match
// exactly when they score smax, the table's best score (NULL otherwise, or
// when seeds cannot help: the scoring has no mismatches, or out of memory).
// k letters per seed (0: by the length of b)
seeds_t *seeds_build(const sub_table_t *table, const scoring_t *scoring,
                     const char *a, size_t len_a,
                     const char *b, size_t len_b, score_t smax, size_t k);
void seeds_free(seeds_t *sd);

// No alignment of a[x..] with b[y..] scores more than this, whatever gap
// the path is in at (x, y)
int64_t seeds_bound(const seeds_t *sd, size_t x, size_t y);

// Narrow the diagonals *lo .. *hi (x - y) of a band's strip of rows y1+1
// .. y2 to those a path scoring lb or more can use, from row y1's scores
// (row: columns c0+1 .. c1 by matrix; col0: column 0 by row, read when
// c0 == 0; min elsewhere)
void seeds_narrow(seeds_t *sd, score_t *const row[3],
                  score_t *const col0[3], size_t y1, size_t c0, size_t c1,
                  size_t y2, score_t smax, int gap_extend, score_t lb,
                  ptrdiff_t *lo, ptrdiff_t *hi);

#endif /* ALIGNMENT_SEEDS_HEADER_SEEN */
