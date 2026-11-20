/*
 alignment_testing.h
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// Private: hooks for the unit tests (src/tools/tests.c), not part of the
// library's interface. They may change or go without notice

#ifndef ALIGNMENT_TESTING_HEADER_SEEN
#define ALIGNMENT_TESTING_HEADER_SEEN

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Fills started with 16-bit lanes, and those started again with 32-bit ones,
// since the program began; and the same for 8-bit lanes (started again with
// 16- or 32-bit ones)
void alignment_narrow_counts(size_t *fills, size_t *retries);
void alignment_byte_counts(size_t *fills, size_t *retries);
// Difference-encoded fills since the program began
size_t alignment_diff_count(void);

// When true, each SIMD pass of hirschberg.c is redone without SIMD and
// hirschberg_striped_mismatches counts any row that differs. While it is
// set, hirschberg_banded_passes counts banded passes (each try of a band)
extern bool hirschberg_check_striped;
extern size_t hirschberg_striped_mismatches, hirschberg_banded_passes;

// Rows per strip of a strip alignment (see smith_waterman_strip_min_cells;
// 0: automatic)
extern size_t smith_waterman_strip_rows;
// Strip alignments made so far (not safe with alignments on other threads)
extern size_t smith_waterman_strip_alignments;

// Set to n > 0: the nth allocation an alignment makes from then on fails,
// as if out of memory (on any thread), and the count returns to 0. Creating
// aligners and results (e.g. alignment_create()) does not count
void alignment_alloc_fail(size_t n);
// Allocations still to go before the one set to fail (0 once it has)
size_t alignment_alloc_fail_pending(void);

#ifdef __cplusplus
}
#endif

#endif /* ALIGNMENT_TESTING_HEADER_SEEN */
