/*
 tests/test_options.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 Per-aligner options, threads, linear-space fallback and linear_above
 */

#include "test_common.h"


/* Linear space scores literal '-' and '\0' characters as Needleman-Wunsch
 * does, and both global aligners describe the whole sequences in a result
 * that last held a Smith-Waterman hit */
void test_global_inputs()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *nw_aln = alignment_create(8), *hb_aln = alignment_create(8);
  scoring_t scoring;
  scoring_system_default(&scoring);
  const char *pairs[][2] = {{"-", "-"}, {"A-C", "AC"}, {"--", "A"},
                            {"AC-GT", "AC-GT"}};
  size_t i;

  for(i = 0; i < sizeof(pairs)/sizeof(pairs[0]); i++) {
    needleman_wunsch_align(pairs[i][0], pairs[i][1], &scoring, nw, nw_aln);
    ASSERT(hirschberg_align(pairs[i][0], pairs[i][1], &scoring, hb_aln));
    ASSERT(nw_aln->score == hb_aln->score);
  }

  // Embedded '\0' (lengths given)
  const char a[] = {'A', '\0', 'T'};
  needleman_wunsch_align2(a, a, 3, 3, &scoring, nw, nw_aln);
  ASSERT(hirschberg_align2(a, a, 3, 3, &scoring, hb_aln));
  ASSERT(nw_aln->score == 3 && hb_aln->score == 3);

  // Reused after a local hit
  smith_waterman_set_min_score(sw, 1);
  smith_waterman_align("TTTTACGT", "ACGT", &scoring, sw);
  ASSERT(smith_waterman_fetch(sw, nw_aln) && nw_aln->pos_a == 4);
  needleman_wunsch_align("TTTTACGT", "ACGT", &scoring, nw, nw_aln);
  ASSERT(nw_aln->pos_a == 0 && nw_aln->pos_b == 0);
  ASSERT(nw_aln->len_a == 8 && nw_aln->len_b == 4);
  hb_aln->pos_a = hb_aln->pos_b = 7;
  ASSERT(hirschberg_align("TTTTACGT", "ACGT", &scoring, hb_aln));
  ASSERT(hb_aln->pos_a == 0 && hb_aln->pos_b == 0);
  ASSERT(hb_aln->len_a == 8 && hb_aln->len_b == 4);

  needleman_wunsch_free(nw);
  smith_waterman_free(sw);
  alignment_free(nw_aln);
  alignment_free(hb_aln);
}

/* Aligners with options of their own give the same results as the defaults,
 * also on two threads at once with different options (no shared settings) */
#define OPT_PAIRS 40
typedef struct {
  const alignment_options_t *opts;
  char (*a)[401], (*b)[801];
  score_t nw[OPT_PAIRS], sw[OPT_PAIRS], hb[OPT_PAIRS];
} opt_run_t;

static void* opt_run(void *ptr)
{
  opt_run_t *r = ptr;
  nw_aligner_t *nw = needleman_wunsch_new();
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *aln = alignment_create(256);
  scoring_t scoring;
  size_t i;
  scoring_system_default(&scoring);
  nw->options = r->opts;
  nw->compact = true;
  smith_waterman_set_options(sw, r->opts);
  smith_waterman_set_compact(sw, true);
  smith_waterman_set_min_score(sw, 1);
  for(i = 0; i < OPT_PAIRS; i++) {
    needleman_wunsch_align(r->a[i], r->b[i], &scoring, nw, aln);
    r->nw[i] = aln->score;
    smith_waterman_align(r->a[i], r->b[i], &scoring, sw);
    r->sw[i] = smith_waterman_fetch(sw, aln) ? aln->score : -1;
    r->hb[i] = hirschberg_align_opts(r->a[i], r->b[i], strlen(r->a[i]),
                                     strlen(r->b[i]), &scoring, r->opts, aln)
               ? aln->score : -1;
  }
  needleman_wunsch_free(nw);
  smith_waterman_free(sw);
  alignment_free(aln);
  return NULL;
}

void options_test_per_aligner()
{
  static char a[OPT_PAIRS][401], b[OPT_PAIRS][801];
  alignment_options_t many, plain;
  size_t i;
  for(i = 0; i < OPT_PAIRS; i++) {
    test_rand_seq(a[i], 50 + test_rand(350), "acgt");
    test_mutate(b[i], a[i], i % 2 ? 50 : 300, "acgt");
  }
  // Every switch the other way, and blocks on several threads
  alignment_options_default(&many);
  many.max_threads = many.hb_max_threads = 4;
  many.block_min_cells = many.hb_block_min_cells = 1;
  many.block_min_cols = many.hb_block_min_cols = 8;
  many.band = many.hb_band = !many.band;
  many.narrow = !many.narrow;
  many.byte = !many.byte;
  many.sw_rect_min_cells = many.sw_strip_min_cells = 0;
  alignment_options_default(&plain);

  // The aligner's own options are the ones used: 8-bit lanes off there
  // start no 8-bit fill, while the defaults do
  {
    sw_aligner_t *sw = smith_waterman_new();
    alignment_options_t off;
    scoring_t scoring;
    size_t fills0, fills1, fills2, retries;
    scoring_system_default(&scoring);
    smith_waterman_set_compact(sw, true);
    alignment_options_default(&off);
    off.byte = false;
    alignment_byte_counts(&fills0, &retries);
    smith_waterman_set_options(sw, &off);
    smith_waterman_align(a[0], b[0], &scoring, sw);
    alignment_byte_counts(&fills1, &retries);
    smith_waterman_set_options(sw, NULL);
    smith_waterman_align(a[0], b[0], &scoring, sw);
    alignment_byte_counts(&fills2, &retries);
    ASSERT(fills1 == fills0);
#ifdef ALN_STRIPED
    ASSERT(fills2 > fills1);
#endif
    smith_waterman_free(sw);
  }

  opt_run_t base = {NULL, a, b, {0}, {0}, {0}};
  opt_run_t r1 = {&many, a, b, {0}, {0}, {0}};
  opt_run_t r2 = {&plain, a, b, {0}, {0}, {0}};
  opt_run(&base);
  pthread_t t1, t2;
  ASSERT(pthread_create(&t1, NULL, opt_run, &r1) == 0);
  ASSERT(pthread_create(&t2, NULL, opt_run, &r2) == 0);
  pthread_join(t1, NULL);
  pthread_join(t2, NULL);
  for(i = 0; i < OPT_PAIRS; i++) {
    ASSERT(r1.nw[i] == base.nw[i] && r2.nw[i] == base.nw[i]);
    ASSERT(r1.sw[i] == base.sw[i] && r2.sw[i] == base.sw[i]);
    ASSERT(r1.hb[i] == base.hb[i] && r2.hb[i] == base.hb[i]);
    ASSERT(base.hb[i] == base.nw[i]);
  }
}

/* Thousands of column blocks must not need a stack to match: the per-block
 * state is on the heap, so a fill on a thread with a small stack (as on
 * macOS, where threads get 512 KB) survives a high thread count */

typedef struct {
  const char *a, *b;
  alignment_options_t *opts;
  score_t nw_score, hb_score;
  bool ok;
} stack_case_t;

static void *stack_case_run(void *ptr)
{
  stack_case_t *c = ptr;
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(16);
  static scoring_t scoring; // 260 KB of tables: not on this thread's stack
  scoring_system_default(&scoring);
  nw->options = c->opts;
  nw->compact = true;
  const size_t la = strlen(c->a), lb = strlen(c->b);
  c->ok = needleman_wunsch_align_checked(c->a, c->b, la, lb, &scoring, nw,
                                         aln, NULL) == ALN_OK;
  c->nw_score = aln->score;
  c->ok = c->ok && hirschberg_align_checked(c->a, c->b, la, lb, &scoring,
                                            c->opts, aln, NULL) == ALN_OK;
  c->hb_score = aln->score;
  needleman_wunsch_free(nw);
  alignment_free(aln);
  return NULL;
}

void test_small_stack_many_blocks()
{
  static char a[32769], b[201];
  test_rand_seq(a, 32768, "acgt");
  test_rand_seq(b, 200, "acgt");
  alignment_options_t one, many;
  alignment_options_default(&one);
  one.max_threads = one.hb_max_threads = 1;
  many = one;
  many.max_threads = many.hb_max_threads = 4096;
  many.block_min_cells = many.hb_block_min_cells = 1;
  many.block_min_cols = many.hb_block_min_cols = 8; // 4096 blocks
  stack_case_t ref = {a, b, &one, 0, 0, false};
  stack_case_t small = {a, b, &many, 0, 0, false};
  stack_case_run(&ref);
  ASSERT(ref.ok);

  pthread_attr_t attr;
  pthread_t thread;
  ASSERT(pthread_attr_init(&attr) == 0);
  ASSERT(pthread_attr_setstacksize(&attr, 384 * 1024) == 0); // under macOS's
  ASSERT(pthread_create(&thread, &attr, stack_case_run, &small) == 0);
  pthread_join(thread, NULL);
  pthread_attr_destroy(&attr);
  ASSERT(small.ok && small.nw_score == ref.nw_score &&
         small.hb_score == ref.hb_score);
}

/* With linear_fallback, a Needleman-Wunsch alignment too big for the
 * matrices (here: out of memory) is made in linear space, the same as
 * hirschberg_align_opts() makes it; without it the error comes back */
void test_linear_fallback()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(16), *lin = alignment_create(16);
  alignment_options_t opts;
  scoring_t scoring;
  aln_error_t e;
  const char *a = "ACGTTGCAACGTAGCTAGCTAGGATCC", *b = "ACGTGCAACGTAGTAGCTAGGTTCC";
  const size_t la = strlen(a), lb = strlen(b);
  alignment_options_default(&opts);
  scoring_system_default(&scoring);
  nw->options = &opts;
  ASSERT(hirschberg_align_checked(a, b, la, lb, &scoring, &opts, lin,
                                  NULL) == ALN_OK);

  // The first allocation is a score matrix
  opts.linear_fallback = false;
  alignment_alloc_fail(1);
  ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &scoring, nw, aln,
                                        &e) == ALN_ERR_NOMEM);
  alignment_alloc_fail(0);
  ASSERT(e.status == ALN_ERR_NOMEM);

  opts.linear_fallback = true;
  alignment_alloc_fail(1);
  ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &scoring, nw, aln,
                                        &e) == ALN_OK);
  alignment_alloc_fail(0);
  ASSERT(e.status == ALN_OK);
  ASSERT(aln->score == lin->score && strcmp(aln->result_a, lin->result_a) == 0 &&
         strcmp(aln->result_b, lin->result_b) == 0);

  needleman_wunsch_free(nw);
  alignment_free(aln);
  alignment_free(lin);
}

/* With linear_above, a pair whose longer sequence is at least that long is
 * made in linear space (as hirschberg_align_opts() makes it), shorter pairs
 * with the matrices; scoring linear space refuses uses the matrices */
void test_linear_above()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(16), *lin = alignment_create(16);
  alignment_t *mat = alignment_create(16);
  alignment_options_t opts;
  scoring_t scoring;
  // Several alignments share the best score here, so the two aligners may
  // choose different ones: which made it shows in the strings
  const char *a = "AAAACCCCAAAAGGGGTTTT", *b = "AAAACCAAAAGGTTTTTT";
  const size_t la = strlen(a), lb = strlen(b);
  alignment_options_default(&opts);
  scoring_system_default(&scoring);
  nw->options = &opts;
  ASSERT(hirschberg_align_checked(a, b, la, lb, &scoring, &opts, lin,
                                  NULL) == ALN_OK);
  opts.linear_above = 0;
  ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &scoring, nw, mat,
                                        NULL) == ALN_OK);
  ASSERT(mat->score == lin->score);

  opts.linear_above = la; // the longer one: linear space
  needleman_wunsch_align(a, b, &scoring, nw, aln);
  ASSERT(aln->score == lin->score && strcmp(aln->result_a, lin->result_a) == 0 &&
         strcmp(aln->result_b, lin->result_b) == 0);
  opts.linear_above = la + 1; // too short: the matrices
  ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &scoring, nw, aln,
                                        NULL) == ALN_OK);
  ASSERT(strcmp(aln->result_a, mat->result_a) == 0 &&
         strcmp(aln->result_b, mat->result_b) == 0);

  // Out of memory in linear space is reported, not retried with the
  // matrices (which need far more)
  opts.linear_above = la;
  size_t k, nomem = 0;
  for(k = 1; k <= 30; k++) {
    aln_error_t e;
    alignment_alloc_fail(k);
    aln_status_t st = needleman_wunsch_align_checked(a, b, la, lb, &scoring,
                                                     nw, aln, &e);
    alignment_alloc_fail(0);
    ASSERT(st == e.status && (st == ALN_OK || st == ALN_ERR_NOMEM));
    if(st == ALN_OK)
      ASSERT(strcmp(aln->result_a, lin->result_a) == 0); // linear space did it
    nomem += st == ALN_ERR_NOMEM;
  }
  ASSERT(nomem > 0);

  // Free end gaps: linear space refuses them, so the matrices align it
  opts.linear_above = 1;
  scoring.no_end_gap_penalty = true;
  ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &scoring, nw, aln,
                                        NULL) == ALN_OK);
  opts.linear_above = 0;
  ASSERT(needleman_wunsch_align_checked(a, b, la, lb, &scoring, nw, mat,
                                        NULL) == ALN_OK);
  ASSERT(aln->score == mat->score && strcmp(aln->result_a, mat->result_a) == 0);

  needleman_wunsch_free(nw);
  alignment_free(aln);
  alignment_free(lin);
  alignment_free(mat);
}

void test_options(void)
{
  SUITE_START("Options and threads");
  RUN(test_global_inputs());
  RUN(options_test_per_aligner());
  RUN(test_small_stack_many_blocks());
  RUN(test_linear_fallback());
  RUN(test_linear_above());
  SUITE_END();
}
