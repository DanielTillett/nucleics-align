/*
 tests/test_wfa.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 Wavefront alignment and scores against the matrices and linear space
 */

#include "test_common.h"


/* With wfa_align, an alignment by the wavefront algorithm: the matrices'
 * score, a real alignment of the two sequences scoring that, and the
 * matrices where the scoring does not suit it */
void test_wfa_align()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(16), *mat = alignment_create(16);
  alignment_options_t opts;
  scoring_t scoring;
  static char seqa[900], seqb[1800];
  const unsigned int rates[] = {0, 5, 30, 150, 400};
  size_t i, used = 0;
  alignment_options_default(&opts);
  nw->options = &opts;

  for(i = 0; i < 400; i++)
  {
    scoring_init(&scoring, 1 + (int)test_rand(4), -(int)test_rand(8),
                 -(int)test_rand(10), -1 - (int)test_rand(4), false, false,
                 false, false, false, i % 2 == 0);
    test_rand_seq(seqa, 1 + test_rand(800), i % 2 ? "acgt" : "ACgt");
    if(i % 6 == 5) test_rand_seq(seqb, 1 + test_rand(800), "acgt");
    else test_mutate(seqb, seqa, rates[i % 5], "acgt");
    if(seqb[0] == '\0') strcpy(seqb, "a");

    opts.wfa_align = false;
    needleman_wunsch_align(seqa, seqb, &scoring, nw, mat);
    opts.wfa_align = true;
    ASSERT(needleman_wunsch_align_checked(seqa, seqb, strlen(seqa),
                                          strlen(seqb), &scoring, nw, aln,
                                          NULL) == ALN_OK);
    ASSERT(aln->score == mat->score);
    ASSERT(hirschberg_score_alignment(aln->result_a, aln->result_b,
                                      &scoring) == aln->score);
    ASSERT(test_ungapped_equals(aln->result_a, seqa) &&
           test_ungapped_equals(aln->result_b, seqb));
    used += !nw->traced && nw->match_scores == NULL;
  }

  // BLOSUM62 is not one match and one mismatch score: the matrices
  scoring_system_BLOSUM62(&scoring);
  opts.wfa_align = false;
  needleman_wunsch_align("MKTAYIAKQR", "MKTAHIAKQR", &scoring, nw, mat);
  opts.wfa_align = true;
  needleman_wunsch_align("MKTAYIAKQR", "MKTAHIAKQR", &scoring, nw, aln);
  ASSERT(aln->score == mat->score && strcmp(aln->result_a, mat->result_a) == 0);
  (void)used;

  // Running out of memory is reported as that, never as "unsupported"
  // (which would send the pair to the matrices), on both the short-pair
  // and the bidirectional paths; and through needleman_wunsch_align_checked()
  static char big_a[3001], big_b[6001];
  scoring_system_default(&scoring);
  sub_table_t table;
  size_t k, nomem;
  const size_t lens[] = {700, 3000};
  for(i = 0; i < 2; i++) {
    test_rand_seq(big_a, lens[i], "acgt");
    test_mutate(big_b, big_a, 20, "acgt");
    const size_t la = strlen(big_a), lb = strlen(big_b);
    ASSERT(sub_table_build(&table, &scoring, big_a, la, big_b, lb));
    for(k = 1, nomem = 0; k <= 60; k++) {
      alignment_alloc_fail(k);
      aln_status_t st = alignment_wfa_align(&table, &scoring, big_a, la,
                                            big_b, lb, 50000000, aln);
      alignment_alloc_fail(0);
      ASSERT(st == ALN_OK || st == ALN_ERR_NOMEM);
      nomem += st == ALN_ERR_NOMEM;
    }
    ASSERT(nomem > 0);
    sub_table_free(&table);
    for(k = 1, nomem = 0; k <= 60; k++) {
      aln_error_t e;
      alignment_alloc_fail(k);
      aln_status_t st = needleman_wunsch_align_checked(big_a, big_b, la, lb,
                                                       &scoring, nw, aln, &e);
      alignment_alloc_fail(0);
      ASSERT(st == e.status && (st == ALN_OK || st == ALN_ERR_NOMEM));
      nomem += st == ALN_ERR_NOMEM;
    }
    ASSERT(nomem > 0);
  }

  needleman_wunsch_free(nw);
  alignment_free(aln);
  alignment_free(mat);
}

/* The wavefront score is the full fill's best global score (Hirschberg's
 * here, which matches it), or it declines */
static void wfa_check(const char *a, const char *b, const scoring_t *scoring,
                      size_t budget, size_t *scored)
{
  const size_t la = strlen(a), lb = strlen(b);
  sub_table_t table;
  score_t wfa, full;
  ASSERT(sub_table_build(&table, scoring, a, la, b, lb));
  if(alignment_wfa_score(&table, scoring, a, la, b, lb, budget, &wfa)) {
    ASSERT(hirschberg_score(a, b, scoring, &full));
    ASSERT(wfa == full);
    if(wfa != full) printf("  wfa %i full %i: %s %s\n", wfa, full, a, b);
    (*scored)++;
  }
  sub_table_free(&table);
}

void test_wfa_score()
{
  scoring_t scoring;
  char seqa[1024], seqb[1024], base[1024];
  size_t p, la, lb, xa, xb, k, scored = 0;
  const size_t big = (size_t)1 << 30;

  // Every pair over {a,c} up to length 5, under several penalties
  const int penalties[][4] = {{1, -1, -1, -1}, {1, -15, -10, -1},
                              {2, -3, 0, -2}, {1, -5, -10, -5},
                              {3, 0, -6, -1}, {1, -20, -2, -1},
                              {0, -1, -2, -1}, {5, -4, -10, -1}};
  for(p = 0; p < sizeof(penalties)/sizeof(penalties[0]); p++) {
    scoring_init(&scoring, penalties[p][0], penalties[p][1],
                 penalties[p][2], penalties[p][3],
                 false, false, false, false, false, true);
    for(la = 1; la <= 5; la++)
      for(xa = 0; xa < (1UL << la); xa++)
        for(lb = 1; lb <= 5; lb++)
          for(xb = 0; xb < (1UL << lb); xb++) {
            for(k = 0; k < la; k++) seqa[k] = (xa >> k) & 1 ? 'c' : 'a';
            for(k = 0; k < lb; k++) seqb[k] = (xb >> k) & 1 ? 'c' : 'a';
            seqa[la] = seqb[lb] = '\0';
            wfa_check(seqa, seqb, &scoring, big, &scored);
          }
  }
  ASSERT(scored == 8 * (62 * 62 - 50)); // not a's only against c's only

  // Random scorings and pairs, related and not, in either case
  for(p = 0; p < 3000; p++) {
    const int match = (int)test_rand(6), mismatch = match - 1 - (int)test_rand(8);
    const int gap_open = -(int)test_rand(13);
    const int gap_extend = -(int)test_rand(5) - (match == 0);
    scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
                 false, false, false, false, false, test_rand(2) == 0);
    la = 1 + test_rand(p % 10 == 0 ? 400 : 60);
    test_rand_seq(base, la, p % 3 == 0 ? "acgtACGT" : "acgt");
    if(test_rand(4) == 0) {
      test_rand_seq(seqb, 1 + test_rand(400), "acgtn");
      strcpy(seqa, base);
    }
    else {
      test_mutate(seqa, base, test_rand(300), "acgt");
      test_mutate(seqb, base, test_rand(300), "acgt");
    }
    if(seqa[0] == '\0' || seqb[0] == '\0') continue;
    wfa_check(seqa, seqb, &scoring, big, &scored);
    wfa_check(seqa, seqb, &scoring, 50, &scored); // gives up, or exact
  }
  ASSERT(scored > 8 * (62 * 62 - 50) + 3000); // nearly all, and some at 50

  // Declines what it does not handle
  sub_table_t table;
  score_t s;
  scoring_init(&scoring, 1, -2, -4, -1, false, false, false, false, false, true);
  ASSERT(sub_table_build(&table, &scoring, "acgt", 4, "acgt", 4));
  ASSERT(alignment_wfa_score(&table, &scoring, "acgt", 4, "acgt", 4, big, &s) &&
         s == 4);
  ASSERT(!alignment_wfa_score(&table, &scoring, "acgt", 4, "", 0, big, &s));
  scoring.no_end_gap_penalty = true;
  ASSERT(!alignment_wfa_score(&table, &scoring, "acgt", 4, "acgt", 4, big, &s));
  sub_table_free(&table);
  scoring_init(&scoring, 1, -2, -4, -1, false, false, false, false, false, true);
  scoring_add_wildcard(&scoring, 'n', 0);
  ASSERT(sub_table_build(&table, &scoring, "acgn", 4, "acgt", 4));
  ASSERT(!alignment_wfa_score(&table, &scoring, "acgn", 4, "acgt", 4, big, &s));
  sub_table_free(&table);
  ASSERT(sub_table_build(&table, &scoring, "aaaa", 4, "cccc", 4));
  ASSERT(!alignment_wfa_score(&table, &scoring, "aaaa", 4, "cccc", 4, big, &s));
  sub_table_free(&table);
}

void test_biwfa_deep(void)
{
  const int params[][4] = {
    {1, -2, -4, -1}, // odd gcd = 1
    {2, -2, -4, -2}, // even gcd = 2
    {3, -5, -7, -2}, // high coprime penalties
  };
  char *a = malloc(4096), *b = malloc(6001);
  ASSERT(a && b);
  alignment_t *aln = alignment_create(8192);

  for(size_t pi = 0; pi < sizeof(params)/sizeof(params[0]); pi++) {
    scoring_t sc;
    scoring_init(&sc, params[pi][0], params[pi][1], params[pi][2], params[pi][3],
                 false, false, false, false, false, true);
    test_rand_seq(a, 3000, "ACGT");
    test_mutate(b, a, 20, "ACGT");
    if(b[0] == '\0') strcpy(b, "A");
    size_t la = strlen(a), lb = strlen(b);

    sub_table_t tbl;
    ASSERT(sub_table_build(&tbl, &sc, a, la, b, lb));
    score_t ref = 0;
    ASSERT(hirschberg_score2(a, b, la, lb, &sc, &ref));

    // Tiny budget must decline cleanly
    ASSERT(alignment_wfa_align(&tbl, &sc, a, la, b, lb, 50, aln) ==
           ALN_ERR_UNSUPPORTED);
    // Sufficient budget must succeed via BiWFA (> 2048 cutoff) and rescore identically
    ASSERT(alignment_wfa_align(&tbl, &sc, a, la, b, lb, 5000000, aln) == ALN_OK);
    ASSERT(aln->score == ref);
    ASSERT(hirschberg_score_alignment(aln->result_a, aln->result_b, &sc) == ref);
    ASSERT(test_ungapped_equals(aln->result_a, a) && test_ungapped_equals(aln->result_b, b));
    sub_table_free(&tbl);
  }
  alignment_free(aln);
  free(a);
  free(b);
}

void test_wfa(void)
{
  SUITE_START("Wavefront alignment");
  RUN(test_wfa_align());
  RUN(test_wfa_score());
  RUN(test_biwfa_deep());
  SUITE_END();
}
