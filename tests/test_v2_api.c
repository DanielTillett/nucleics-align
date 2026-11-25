/*
 tests/test_v2_api.c
 Pure public-API test suite: includes ONLY include/nucleics_align.h.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nucleics_align.h"

#define CHECK(cond) do { \
  if(!(cond)) { \
    fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
    exit(EXIT_FAILURE); \
  } \
} while(0)

static void test_version_and_errors(void)
{
  CHECK(strcmp(na_version(), NA_VERSION_STR) == 0);
  CHECK(na_strerror(NA_OK) != NULL);
  CHECK(na_strerror(NA_ERR_NOMEM) != NULL);
  CHECK(na_scoring_preset_names() != NULL);

  na_error_t err = { NA_ERR_SCORING_PAIR, 'X', 'Y', -1, "" };
  char buf[128];
  size_t n = na_error_format(&err, buf, sizeof(buf));
  CHECK(n > 0 && strstr(buf, "'X'") != NULL && strstr(buf, "'Y'") != NULL);

  na_error_t perr = { NA_ERR_PARSE, 0, 0, 7, "bad value" };
  n = na_error_format(&perr, buf, sizeof(buf));
  CHECK(n > 0 && strstr(buf, "7") != NULL && strstr(buf, "bad value") != NULL);
  printf("  [PASS] test_version_and_errors\n");
}

static void test_scoring_and_parsers(void)
{
  na_scoring_t *b62 = na_scoring_new_preset("BLOSUM62");
  CHECK(b62 != NULL);
  int s = 0; bool m = false;
  CHECK(na_scoring_lookup(b62, 'W', 'W', &s, &m) == NA_OK && s == 11 && m);
  na_scoring_free(b62);

  CHECK(na_scoring_new_preset("NO_SUCH_PRESET") == NULL);

  na_scoring_t *sc = na_scoring_new(2, -1, -4, -1);
  CHECK(sc != NULL);
  CHECK(na_scoring_set_match_mismatch(sc, false, 0, 0) == NA_OK);
  CHECK(na_scoring_set_flags(sc, NA_CASE_SENSITIVE) == NA_OK);
  /* Clearing NA_CASE_SENSITIVE after setting it is rejected */
  CHECK(na_scoring_set_flags(sc, 0) == NA_ERR_INVALID_ARG);
  /* Invalid flag bit is rejected */
  CHECK(na_scoring_set_flags(sc, 1u << 20) == NA_ERR_INVALID_ARG);

  const char *mat = "# comment\n  A C\nA 5 -3\nC -3 5\n";
  na_error_t err = {0};
  CHECK(na_scoring_parse_matrix(sc, mat, strlen(mat), &err) == NA_OK);
  CHECK(na_scoring_lookup(sc, 'A', 'A', &s, &m) == NA_OK && s == 5);
  CHECK(na_scoring_lookup(sc, 'A', 'C', &s, &m) == NA_OK && s == -3);
  CHECK(na_scoring_lookup(sc, 'A', 'G', &s, &m) == NA_ERR_SCORING_PAIR);

  /* Header-only matrix rejected */
  const char *hdr_only = "  A C\n";
  CHECK(na_scoring_parse_matrix(sc, hdr_only, strlen(hdr_only), &err) == NA_ERR_PARSE);
  /* Embedded NUL in numeric token rejected */
  const char nul_mat[] = "  A C\nA 5\0junk -3\n";
  CHECK(na_scoring_parse_matrix(sc, nul_mat, sizeof(nul_mat) - 1, &err) == NA_ERR_PARSE);
  /* NULL text rejected */
  CHECK(na_scoring_parse_matrix(sc, NULL, 0, &err) == NA_ERR_INVALID_ARG);

  /* Transactional failure: bad matrix leaves sc unchanged */
  const char *bad_mat = "  A C\nA 5 bad\n";
  CHECK(na_scoring_parse_matrix(sc, bad_mat, strlen(bad_mat), &err) == NA_ERR_PARSE);
  CHECK(err.line == 2);
  CHECK(na_scoring_lookup(sc, 'A', 'C', &s, &m) == NA_OK && s == -3);

  /* Pairs parser */
  const char *pairs = "G G 4\nG,A,-2\n";
  CHECK(na_scoring_parse_pairs(sc, pairs, strlen(pairs), &err) == NA_OK);
  CHECK(na_scoring_lookup(sc, 'G', 'G', &s, &m) == NA_OK && s == 4);
  CHECK(na_scoring_lookup(sc, 'G', 'A', &s, &m) == NA_OK && s == -2);

  na_scoring_free(sc);
  printf("  [PASS] test_scoring_and_parsers\n");
}

static void test_global_methods_and_edge_cases(void)
{
  na_scoring_t *sc = na_scoring_new(2, -4, -4, -2);
  na_aligner_t *al = na_aligner_new();
  na_alignment_t *res = na_alignment_new();
  CHECK(sc && al && res);

  const char *a = "ACGTACGTACGTACGTACGT";
  const char *b = "ACGTTCGTACGTACGGACGT";
  size_t la = strlen(a), lb = strlen(b);

  na_score_t ref_score = 0;
  na_method_t methods[] = { NA_METHOD_AUTO, NA_METHOD_FULL, NA_METHOD_LINEAR, NA_METHOD_WFA };
  for(size_t i = 0; i < 4; i++) {
    na_aligner_set_method(al, methods[i]);
    na_aligner_set_threads(al, 1);
    na_aligner_set_wfa_budget(al, 4096);
    CHECK(na_global_align(al, a, la, b, lb, sc, res) == NA_OK);
    na_score_t checked = 0;
    CHECK(na_score_alignment(na_alignment_a(res), na_alignment_b(res),
                             na_alignment_length(res), sc, &checked) == NA_OK);
    CHECK(checked == na_alignment_score(res));
    if(i == 0) ref_score = checked;
    else CHECK(checked == ref_score);
  }

  /* Matrix dump to non-stdout FILE* works for FULL and is safely suppressed after AUTO */
  na_aligner_set_method(al, NA_METHOD_FULL);
  na_aligner_set_linear_above(al, 1); /* FULL must still build full matrices */
  CHECK(na_global_align(al, "ACG", 3, "ATG", 3, sc, res) == NA_OK);
  FILE *tf = tmpfile();
  CHECK(tf != NULL);
  na_debug_fprint_matrices(al, tf);
  CHECK(ftell(tf) > 50);
  fclose(tf);

  na_aligner_set_method(al, NA_METHOD_AUTO);
  na_aligner_set_linear_above(al, 0);
  CHECK(na_global_align(al, a, la, b, lb, sc, res) == NA_OK);
  tf = tmpfile();
  CHECK(tf != NULL);
  na_debug_fprint_matrices(al, tf); /* must not read stale 4x4 matrix with 21x21 dims */
  CHECK(ftell(tf) == 0);
  fclose(tf);

  /* Free start/end gaps parity in na_score_alignment() */
  na_scoring_t *fg = na_scoring_new(2, -3, -4, -1);
  CHECK(na_scoring_set_flags(fg, NA_FREE_START_GAP | NA_FREE_END_GAP) == NA_OK);
  CHECK(na_global_align(al, "ACGT", 4, "TTTTACGTTTTT", 12, fg, res) == NA_OK);
  na_score_t fg_checked = 0;
  CHECK(na_score_alignment(na_alignment_a(res), na_alignment_b(res),
                           na_alignment_length(res), fg, &fg_checked) == NA_OK);
  CHECK(fg_checked == na_alignment_score(res));
  na_scoring_free(fg);

  /* Free gaps as the DP applies them: a gap in a is free at the start while
   * no letter of a is consumed, and at the end once all of a is; likewise
   * for b. Not "every leading column with a gap in either sequence" */
  na_scoring_t *fe = na_scoring_new(1, -100, -4, -1);
  CHECK(na_scoring_set_flags(fe, NA_FREE_END_GAP) == NA_OK);
  CHECK(na_global_align(al, "ACGT", 4, "", 0, fe, res) == NA_OK);
  CHECK(na_alignment_score(res) == 0);
  CHECK(na_score_alignment(na_alignment_a(res), na_alignment_b(res),
                           na_alignment_length(res), fe, &fg_checked) == NA_OK);
  CHECK(fg_checked == 0);
  na_scoring_free(fe);
  na_scoring_t *fs = na_scoring_new(1, -100, -4, -1);
  CHECK(na_scoring_set_flags(fs, NA_FREE_START_GAP) == NA_OK);
  CHECK(na_global_align(al, "AC", 2, "GT", 2, fs, res) == NA_OK);
  CHECK(na_alignment_score(res) == -6); /* one free gap, one charged */
  CHECK(na_score_alignment(na_alignment_a(res), na_alignment_b(res),
                           na_alignment_length(res), fs, &fg_checked) == NA_OK);
  CHECK(fg_checked == -6);
  CHECK(na_score_alignment("AC--", "--GT", 4, fs, &fg_checked) == NA_OK);
  CHECK(fg_checked == -6);
  CHECK(na_score_alignment("--AC", "GT--", 4, fs, &fg_checked) == NA_OK);
  CHECK(fg_checked == -6);
  na_scoring_free(fs);

  /* A rejected na_scoring_set_flags() (it would clear NA_CASE_SENSITIVE)
   * leaves the scoring as it was: "T" against "ACGT" under a free end gap
   * only is "T----"/"-ACGT", one charged gap of length 1 then free trailing
   * gaps (-5); with a free start gap instead it would score +1 */
  na_scoring_t *cs = na_scoring_new(1, -100, -4, -1);
  CHECK(na_scoring_set_flags(cs, NA_CASE_SENSITIVE | NA_FREE_END_GAP) == NA_OK);
  CHECK(na_global_align(al, "T", 1, "ACGT", 4, cs, res) == NA_OK);
  CHECK(na_alignment_score(res) == -5);
  CHECK(na_scoring_set_flags(cs, NA_FREE_START_GAP) == NA_ERR_INVALID_ARG);
  CHECK(na_global_align(al, "T", 1, "ACGT", 4, cs, res) == NA_OK);
  CHECK(na_alignment_score(res) == -5);
  na_scoring_free(cs);

  /* A mismatched column is an error under NA_NO_MISMATCHES, not score 0 */
  na_scoring_t *nm = na_scoring_new(1, -1, -4, -1);
  CHECK(na_scoring_set_flags(nm, NA_NO_MISMATCHES) == NA_OK);
  CHECK(na_score_alignment("A", "C", 1, nm, &fg_checked) == NA_ERR_SCORING_PAIR);
  CHECK(na_score_alignment("A", "A", 1, nm, &fg_checked) == NA_OK && fg_checked == 1);
  na_scoring_free(nm);

  /* Unsupported flag combination rejected instead of fatal traceback exit */
  na_scoring_t *badf = na_scoring_new(1, -1, -2, -1);
  CHECK(na_scoring_set_flags(badf, NA_NO_MISMATCHES | NA_NO_GAPS_A) == NA_OK);
  CHECK(na_global_align(al, "AA", 2, "CC", 2, badf, res) == NA_ERR_UNSUPPORTED);
  na_scoring_free(badf);

  na_score_t sc_only = 0;
  CHECK(na_global_score(al, a, la, b, lb, sc, &sc_only) == NA_OK);
  CHECK(sc_only == ref_score);

  /* Unit-cost Myers via na_global_score */
  na_scoring_t *ed = na_scoring_new(0, -1, 0, -1);
  CHECK(na_global_score(al, a, la, b, lb, ed, &sc_only) == NA_OK);
  CHECK(sc_only == -2);
  na_scoring_free(ed);

  /* Error propagation clears out */
  CHECK(na_global_align(al, NULL, 5, b, lb, sc, res) == NA_ERR_INVALID_ARG);
  CHECK(na_aligner_last_error(al)->status == NA_ERR_INVALID_ARG);
  CHECK(na_alignment_length(res) == 0 && na_alignment_score(res) == 0);

  na_alignment_free(res);
  na_aligner_free(al);
  na_scoring_free(sc);
  printf("  [PASS] test_global_methods_and_edge_cases\n");
}

static void test_local_iteration_and_invalidation(void)
{
  na_scoring_t *sc = na_scoring_new(2, -3, -4, -1);
  na_aligner_t *al = na_aligner_new();
  na_alignment_t *res = na_alignment_new();
  CHECK(sc && al && res);

  const char *a = "TTTTACGTACGTTTTTGGGGGGGGTTTT";
  const char *b = "CCCCACGTACGTCCCCGGGGGGGGCCCC";
  na_aligner_set_min_local_score(al, 8);
  na_aligner_set_fast_local(al, true);
  CHECK(na_local_align(al, a, strlen(a), b, strlen(b), sc) == NA_OK);

  int hits = 0;
  na_score_t prev = 1000000;
  na_status_t st;
  while((st = na_local_next(al, res)) == NA_OK) {
    CHECK(na_alignment_score(res) >= 8);
    CHECK(na_alignment_score(res) <= prev);
    prev = na_alignment_score(res);
    size_t pa = 0, la = 0, pb = 0, lb = 0;
    na_alignment_span(res, &pa, &la, &pb, &lb);
    CHECK(la > 0 && lb > 0);
    hits++;
  }
  CHECK(st == NA_DONE);
  CHECK(hits >= 2);

  /* Local iterator is invalidated by a subsequent na_global_align() */
  CHECK(na_local_align(al, a, strlen(a), b, strlen(b), sc) == NA_OK);
  CHECK(na_global_align(al, "ACGT", 4, "ACGT", 4, sc, res) == NA_OK);
  CHECK(na_local_next(al, res) == NA_DONE);

  /* NULL-safe free */
  na_alignment_free(NULL);
  na_aligner_free(NULL);
  na_scoring_free(NULL);

  na_alignment_free(res);
  na_aligner_free(al);
  na_scoring_free(sc);
  printf("  [PASS] test_local_iteration_and_invalidation\n");
}

int main(void)
{
  printf("Running v2.0 public API tests...\n");
  test_version_and_errors();
  test_scoring_and_parsers();
  test_global_methods_and_edge_cases();
  test_local_iteration_and_invalidation();
  printf("All v2.0 public API tests passed!\n");
  return 0;
}
