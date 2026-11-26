/*
 tests/test_regressions.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 Regression tests for bugs fixed in earlier versions
 */

#include "test_common.h"


/* Regression tests for overflow bounds, Myers flags, and case-insensitive
   scoring_add_mutation() (the command line parsers' ranges are pinned by
   the fixtures in tests/cli/cases) */
void test_bug_regressions()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(256);
  scoring_t scoring;
  aln_error_t err;
  char seq_a[2149], seq_c[2149];
  memset(seq_a, 'A', 2148); seq_a[2148] = '\0';
  memset(seq_c, 'C', 2148); seq_c[2148] = '\0';

  // Edge scores in [INT_MIN, min) must report ALN_ERR_OVERFLOW
  scoring_init(&scoring, 1, -2, -1000000, -1000000,
               false, false, false, false, false, true);
  ASSERT(needleman_wunsch_align_checked(seq_a, "A", 2146, 1, &scoring,
                                        nw, aln, &err) == ALN_ERR_OVERFLOW);
  ASSERT(needleman_wunsch_align_checked("A", seq_c, 1, 2146, &scoring,
                                        nw, aln, &err) == ALN_ERR_OVERFLOW);
  ASSERT(needleman_wunsch_align_checked(seq_a, "C", 2146, 1, &scoring,
                                        nw, aln, &err) == ALN_ERR_OVERFLOW);

  // no_mismatches all-indel downward overflow must report ALN_ERR_OVERFLOW
  scoring_init(&scoring, 1, -2, -1000000, -1000000,
               false, false, false, false, true, true);
  ASSERT(needleman_wunsch_align_checked(seq_a, seq_c, 1100, 1100, &scoring,
                                        nw, aln, &err) == ALN_ERR_OVERFLOW);
  ASSERT(needleman_wunsch_align_checked(seq_a, seq_c, 1072, 1072, &scoring,
                                        nw, aln, &err) == ALN_ERR_OVERFLOW);

  // Free start gaps on an empty-sequence alignment must not be rejected
  scoring_init(&scoring, 1, -2, -1000000, -1000000,
               true, false, false, false, false, true);
  ASSERT(needleman_wunsch_align_checked(seq_a, "", 2148, 0, &scoring,
                                        nw, aln, &err) == ALN_OK &&
         aln->score == 0);

  // alignment_myers_score() declines unsupported constraint flags
  sub_table_t table;
  score_t s = 0;
  scoring_init(&scoring, 0, -1, 0, -1, false, false, false, false, false, true);
  ASSERT(sub_table_build(&table, &scoring, "AA", 2, "A", 1));
  ASSERT(alignment_myers_score(&table, &scoring, "AA", 2, "A", 1, &s) && s == -1);
  scoring.no_start_gap_penalty = true;
  ASSERT(!alignment_myers_score(&table, &scoring, "AA", 2, "A", 1, &s));
  scoring.no_start_gap_penalty = false;
  scoring.no_end_gap_penalty = true;
  ASSERT(!alignment_myers_score(&table, &scoring, "AA", 2, "A", 1, &s));
  scoring.no_end_gap_penalty = false;
  scoring.no_gaps_in_a = true;
  ASSERT(!alignment_myers_score(&table, &scoring, "AA", 2, "A", 1, &s));
  scoring.no_gaps_in_a = false;
  scoring.no_gaps_in_b = true;
  ASSERT(!alignment_myers_score(&table, &scoring, "AA", 2, "A", 1, &s));
  scoring.no_gaps_in_b = false;
  scoring.no_gaps_in_longest = true;
  ASSERT(!alignment_myers_score(&table, &scoring, "AA", 2, "A", 1, &s));
  scoring.no_gaps_in_longest = false;
  scoring.no_mismatches = true;
  ASSERT(!alignment_myers_score(&table, &scoring, "AA", 2, "A", 1, &s));
  sub_table_free(&table);

  // scoring_add_mutation() normalizes case when !case_sensitive
  scoring_init(&scoring, 1, -2, -4, -1, false, false, false, false, false, false);
  scoring_add_mutation(&scoring, 'A', 'G', 5);
  int sc = 0;
  bool is_match = false;
  scoring_lookup(&scoring, 'a', 'g', &sc, &is_match);
  ASSERT(sc == 5);
  scoring_lookup(&scoring, 'A', 'G', &sc, &is_match);
  ASSERT(sc == 5);
  scoring_lookup(&scoring, 'A', 'g', &sc, &is_match);
  ASSERT(sc == 5);

  needleman_wunsch_free(nw);
  alignment_free(aln);
}

void test_myers_boundaries(void)
{
  scoring_t sc;
  scoring_init(&sc, 0, -1, 0, -1, false, false, false, false, false, true);
  const size_t lens[] = {1, 63, 64, 65, 127, 128, 129, 512, 2048};
  const unsigned int rates[] = {0, 5, 64, 150, 500};
  alignment_options_t opts;
  alignment_options_default(&opts);
  opts.hb_max_threads = 1;
  alignment_t *aln = alignment_create(64);

  for(size_t li = 0; li < sizeof(lens)/sizeof(lens[0]); li++) {
    for(size_t ri = 0; ri < sizeof(rates)/sizeof(rates[0]); ri++) {
      size_t la = lens[li];
      char *a = malloc(la + 1), *b = malloc(la * 2 + 64);
      ASSERT(a && b);
      test_rand_seq(a, la, "ACGT");
      test_mutate(b, a, rates[ri], "ACGT");
      size_t lb = strlen(b);

      sub_table_t tbl;
      ASSERT(sub_table_build(&tbl, &sc, a, la, b, lb));
      score_t ms = 0, bs = 0, hs = 0;
      ASSERT(alignment_myers_score(&tbl, &sc, a, la, b, lb, &ms));
      if(alignment_myers_score_banded(&tbl, &sc, a, la, b, lb, &bs))
        ASSERT(bs == ms);
      ASSERT(hirschberg_align_checked(a, b, la, lb, &sc, &opts, aln, NULL) == ALN_OK);
      hs = aln->score;
      ASSERT(ms == hs);
      ASSERT(hirschberg_score_alignment(aln->result_a, aln->result_b, &sc) == ms);
      ASSERT(test_ungapped_equals(aln->result_a, a) && test_ungapped_equals(aln->result_b, b));
      sub_table_free(&tbl);
      free(a);
      free(b);
    }
  }

  // Codex row-zero frontier regression: leading gap in b of length K = 128 / 160
  // where block 0's lower bound exceeds K at column 1 while row 0 is still viable.
  {
    const size_t cases[][2] = {{257, 128}, {320, 160}};
    for(size_t ci = 0; ci < 2; ci++) {
      size_t ma = cases[ci][0], pad = cases[ci][1], nb = ma + pad;
      char *a = malloc(ma + 1), *b = malloc(nb + 1);
      ASSERT(a && b);
      memset(a, 'A', ma); a[ma] = '\0';
      memset(b, 'C', pad); memset(b + pad, 'A', ma); b[nb] = '\0';
      sub_table_t tbl;
      ASSERT(sub_table_build(&tbl, &sc, a, ma, b, nb));
      score_t ms = 0, bs = 0;
      ASSERT(alignment_myers_score(&tbl, &sc, a, ma, b, nb, &ms));
      ASSERT(ms == -(score_t)pad);
      ASSERT(alignment_myers_score_banded(&tbl, &sc, a, ma, b, nb, &bs));
      ASSERT(bs == ms);
      sub_table_free(&tbl);
      free(a);
      free(b);
    }
  }
  alignment_free(aln);
}

void test_strip_col0_and_striped_trace(void)
{
  // 1. Codex counterexample: TAA+32*C vs ATA+32*C (+5/-4/-2/-1) with strip_rows=1.
  // The 32-C tail makes width-1 banding eligible (asserting trace_banded == true)
  // while preserving the exact leading-gap tie (-TAA / ATA-) at column 0.
  scoring_t sc;
  scoring_init(&sc, 5, -4, -2, -1, false, false, false, false, false, true);
  alignment_options_t o_scalar, o_band;
  alignment_options_default(&o_scalar);
  o_scalar.striped = false; o_scalar.band = false; o_scalar.diff = false; o_scalar.max_threads = 1;
  alignment_options_default(&o_band);
  o_band.striped = true; o_band.band = true; o_band.diff = false;
  o_band.band_strip_rows = 1; o_band.band_min_width = 1; o_band.max_threads = 1;

  char ta[36], tb[36];
  memcpy(ta, "TAA", 3); memset(ta + 3, 'C', 32); ta[35] = '\0';
  memcpy(tb, "ATA", 3); memset(tb + 3, 'C', 32); tb[35] = '\0';

  nw_aligner_t *nw1 = needleman_wunsch_new(), *nw2 = needleman_wunsch_new();
  alignment_t *r1 = alignment_create(128), *r2 = alignment_create(128);
  nw1->compact = false; nw1->options = &o_scalar;
  nw2->compact = true;  nw2->options = &o_band;
  ASSERT(needleman_wunsch_align_checked(ta, tb, 35, 35, &sc, nw1, r1, NULL) == ALN_OK);
  ASSERT(needleman_wunsch_align_checked(ta, tb, 35, 35, &sc, nw2, r2, NULL) == ALN_OK);
  ASSERT_STRIPED(nw2->trace_banded); // bands need the striped fill
  ASSERT(r1->score == r2->score);
  ASSERT(strcmp(r1->result_a, r2->result_a) == 0);
  ASSERT(strcmp(r1->result_b, r2->result_b) == 0);

  // 2. Striped vector-order trace (diff8 DNA + diff16 BLOSUM62, 1t & 4t) byte-identical to scalar
  scoring_t sc_dna, sc_b62;
  scoring_init(&sc_dna, 1, -2, -4, -1, false, false, false, false, false, true);
  scoring_system_BLOSUM62(&sc_b62);
  alignment_options_t o_diff1, o_diff4;
  alignment_options_default(&o_diff1);
  o_diff1.striped = true; o_diff1.band = false; o_diff1.diff = true; o_diff1.max_threads = 1;
  alignment_options_default(&o_diff4);
  o_diff4.striped = true; o_diff4.band = false; o_diff4.diff = true;
  o_diff4.max_threads = 4; o_diff4.block_min_cells = 1; o_diff4.block_min_cols = 32;

  static const char aa20[] = "ACDEFGHIKLMNPQRSTVWY";
  size_t test_lens[] = {17, 63, 64, 65, 257, 513};
  for(size_t li = 0; li < sizeof(test_lens)/sizeof(test_lens[0]); li++) {
    size_t la = test_lens[li], lb = test_lens[li] + (li % 3) * 7;
    char *a = malloc(la + 1), *b = malloc(lb + 1);
    ASSERT(a && b);
    test_rand_seq(a, la, "ACGT");
    test_rand_seq(b, lb, "ACGT");
    // DNA diff8
    nw1->options = &o_scalar;
    ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &sc_dna, nw1, r1, NULL) == ALN_OK);
    nw2->options = &o_diff1;
    ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &sc_dna, nw2, r2, NULL) == ALN_OK);
    ASSERT_STRIPED(nw2->trace_striped);
    size_t dna_vl = nw2->trace_vl;
    ASSERT(r1->score == r2->score && strcmp(r1->result_a, r2->result_a) == 0 && strcmp(r1->result_b, r2->result_b) == 0);
    nw2->options = &o_diff4;
    ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &sc_dna, nw2, r2, NULL) == ALN_OK);
    ASSERT_STRIPED(nw2->trace_striped && nw2->trace_vl == dna_vl);
    ASSERT(r1->score == r2->score && strcmp(r1->result_a, r2->result_a) == 0 && strcmp(r1->result_b, r2->result_b) == 0);
    // Protein BLOSUM62 diff16 (ensure 'W' is present so score range forces 16-bit lanes)
    for(size_t i = 0; i < la; i++) a[i] = aa20[test_rand(20)];
    for(size_t i = 0; i < lb; i++) b[i] = aa20[test_rand(20)];
    a[0] = 'W'; b[0] = 'W';
    nw1->options = &o_scalar;
    ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &sc_b62, nw1, r1, NULL) == ALN_OK);
    nw2->options = &o_diff1;
    ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &sc_b62, nw2, r2, NULL) == ALN_OK);
    ASSERT_STRIPED(nw2->trace_striped && nw2->trace_vl * 2 == dna_vl);
    ASSERT(r1->score == r2->score && strcmp(r1->result_a, r2->result_a) == 0 && strcmp(r1->result_b, r2->result_b) == 0);
    nw2->options = &o_diff4;
    ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &sc_b62, nw2, r2, NULL) == ALN_OK);
    ASSERT_STRIPED(nw2->trace_striped && nw2->trace_vl * 2 == dna_vl);
    ASSERT(r1->score == r2->score && strcmp(r1->result_a, r2->result_a) == 0 && strcmp(r1->result_b, r2->result_b) == 0);
    free(a); free(b);
  }
  alignment_free(r1); alignment_free(r2);
  needleman_wunsch_free(nw1); needleman_wunsch_free(nw2);
}

void test_regressions(void)
{
  SUITE_START("Regressions");
  RUN(test_bug_regressions());
  RUN(test_myers_boundaries());
  RUN(test_strip_col0_and_striped_trace());
  SUITE_END();
}
