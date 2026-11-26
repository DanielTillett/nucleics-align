/*
 tests/test_buffers.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 Alignment buffers, matrix growth, printing and the constructors
 */

#include "test_common.h"


void alignment_test_create_zero()
{
  alignment_t *result = alignment_create(0);
  ASSERT(result->capacity >= 1 && result->result_a[0] == '\0');
  alignment_ensure_capacity(result, 100);
  ASSERT(result->capacity >= 101);
  alignment_free(result);
  alignment_free(NULL); // as free(NULL): nothing

  // A length that cannot be held (SIZE_MAX letters need SIZE_MAX + 1
  // bytes) is refused, not wrapped to a no-op
  result = alignment_create(0);
  ASSERT(!alignment_reserve(result, SIZE_MAX));
  alignment_free(result);
}

// The full score matrices grow by a quarter past 1M cells, as the traceback
// buffer does, not to the next power of two (26 GB for a 13 GB pair)
void alignment_test_matrix_growth()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(16);
  static scoring_t scoring;
  static char a[1051], b[1051];
  scoring_system_default(&scoring);
  test_rand_seq(a, 1050, "acgt");
  test_rand_seq(b, 1050, "acgt");
  nw->compact = false;
  needleman_wunsch_align(a, b, &scoring, nw, aln);
  const size_t cells = 1051 * 1051;
  ASSERT(nw->match_scores != NULL && !nw->traced);
  ASSERT(nw->capacity >= cells && nw->capacity <= cells + cells / 4);
  needleman_wunsch_free(nw);
  alignment_free(aln);
}

// Colours: one per column (red a gap, green a mismatch), switched only
// where the next column's differs; a gap right after a mismatch stays red
void alignment_test_colour()
{
  char *text = NULL;
  size_t len = 0;
  FILE *out = open_memstream(&text, &len);
  alignment_colour_fprint_against(out, "CAAA", "T--A", 0);
  fclose(out);
  ASSERT(strcmp(text, "\033[92mC\033[0m\033[91mAA\033[0mA") == 0);
  free(text);
}

/* Rounding sizes up to a power of two never wraps to a smaller size */
void alignment_test_roundup()
{
  ASSERT(ROUNDUP2POW(0) == 1);
  ASSERT(ROUNDUP2POW(5) == 8);
  ASSERT(ROUNDUP2POW(((size_t)1 << 20) + 1) == (size_t)1 << 21);
  size_t top = (SIZE_MAX >> 1) + 1; // the highest power of two
  ASSERT(ROUNDUP2POW(top) == top);
  ASSERT(ROUNDUP2POW(top + 1) == top + 1);
  ASSERT(ROUNDUP2POW(SIZE_MAX) == SIZE_MAX);
}

/* The constructors that return NULL instead of exiting give working
 * aligners and results; a result that cannot be held is NULL */
void alignment_test_try_new()
{
  nw_aligner_t *nw = needleman_wunsch_try_new();
  sw_aligner_t *sw = smith_waterman_try_new();
  alignment_t *aln = alignment_try_create(8);
  ASSERT(nw != NULL && sw != NULL && aln != NULL);
  scoring_t scoring;
  scoring_system_default(&scoring);
  needleman_wunsch_align("ACGT", "AGT", &scoring, nw, aln);
  ASSERT(aln->length == 4 && test_ungapped_equals(aln->result_a, "ACGT") &&
         test_ungapped_equals(aln->result_b, "AGT"));
  smith_waterman_set_min_score(sw, 1);
  smith_waterman_align("ACGT", "AGT", &scoring, sw);
  ASSERT(smith_waterman_fetch(sw, aln) && aln->score == 2); // GT
  // (A capacity that cannot be allocated gives NULL too, but a request of
  // that size aborts under AddressSanitizer, so it is not asked for here;
  // alignment_test_create_zero() covers the length that cannot be held)
  needleman_wunsch_free(nw);
  smith_waterman_free(sw);
  alignment_free(aln);
}

void test_buffers(void)
{
  SUITE_START("Alignment buffers");
  RUN(alignment_test_create_zero());
  RUN(alignment_test_try_new());
  RUN(alignment_test_matrix_growth());
  RUN(alignment_test_colour());
  RUN(alignment_test_roundup());
  SUITE_END();
}
