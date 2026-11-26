/*
 tests/test_errors.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 The checked aligners report errors in the input and running out of memory
 */

#include "test_common.h"


/* The checked aligners return errors in the input, naming the pair
 * scoring_lookup() would stop at, rather than ending the program */
void test_checked()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *aln = alignment_create(16);
  scoring_t scoring;
  aln_error_t e;
  score_t score;

  // Pairs: only a/c defined (lower case), no match/mismatch
  scoring_init(&scoring, 1, -1, -2, -1, false, false, false, false, false,
               false);
  scoring.use_match_mismatch = false;
  scoring_add_mutation(&scoring, 'a', 'a', 1);
  scoring_add_mutation(&scoring, 'a', 'c', -1);
  scoring_add_mutation(&scoring, 'c', 'a', -1);
  scoring_add_mutation(&scoring, 'c', 'c', 1);

  ASSERT(needleman_wunsch_align_checked("ACCA", "ACA", 4, 3, &scoring, nw,
                                        aln, &e) == ALN_OK);
  ASSERT(aln->score == 0 && e.status == ALN_OK); // ACCA / AC-A
  // As sub_table_build() looks them up: by character of b, then of a
  ASSERT(needleman_wunsch_align_checked("AGT", "CA", 3, 2, &scoring, nw,
                                        aln, &e) == ALN_ERR_SCORING);
  ASSERT(e.a == 'g' && e.b == 'a');
  ASSERT(smith_waterman_align_checked("ACT", "AC", 3, 2, &scoring, sw,
                                      &e) == ALN_ERR_SCORING);
  ASSERT(e.a == 't' && e.b == 'a');
  // A failed check leaves no hits from the pair aligned before it
  ASSERT(smith_waterman_align_checked("ACCA", "ACA", 4, 3, &scoring, sw,
                                      &e) == ALN_OK);
  ASSERT(smith_waterman_align_checked("ACT", "AC", 3, 2, &scoring, sw,
                                      &e) == ALN_ERR_SCORING);
  ASSERT(!smith_waterman_fetch_checked(sw, aln, &e) && e.status == ALN_OK);
  ASSERT(hirschberg_align_checked("CA", "AT", 2, 2, &scoring, NULL, aln,
                                  &e) == ALN_ERR_SCORING);
  ASSERT(e.a == 'a' && e.b == 't');
  ASSERT(hirschberg_score_checked("CA", "AC", 2, 2, &scoring, NULL, &score,
                                  NULL) == ALN_OK && score == -2);

  // Scores that could overflow; options linear space refuses
  scoring_init(&scoring, 300000000, -1, -2, -1, false, false, false, false,
               false, false);
  ASSERT(needleman_wunsch_align_checked("ACGT", "ACGT", 4, 4, &scoring, nw,
                                        aln, &e) == ALN_ERR_OVERFLOW);
  ASSERT(smith_waterman_align_checked("ACGT", "ACGT", 4, 4, &scoring, sw,
                                      &e) == ALN_ERR_OVERFLOW);
  // The score-only scan refuses them too (a scan reports scores as found)
  score_t best;
  size_t bx, by;
  ASSERT(!alignment_sw_scan("AAAAA", 5, "AAAAA", 5, &scoring, 0, NULL, &best,
                            &bx, &by));
  scoring_init(&scoring, 1, -1, 2, -1, false, false, false, false, false,
               false);
  ASSERT(hirschberg_align_checked("AC", "AC", 2, 2, &scoring, NULL, aln,
                                  &e) == ALN_ERR_UNSUPPORTED);

  // Smith-Waterman: a gap that scores above 0 would grow out of the empty
  // (zero) cells, so such scores are refused; Needleman-Wunsch takes them
  scoring_init(&scoring, 2, -3, -2, 1, false, false, false, false, false,
               false);
  ASSERT(smith_waterman_align_checked("CA", "GA", 2, 2, &scoring, sw, &e)
         == ALN_ERR_UNSUPPORTED);
  ASSERT(needleman_wunsch_align_checked("CA", "GA", 2, 2, &scoring, nw, aln,
                                        &e) == ALN_OK);
  scoring_init(&scoring, 2, -3, 3, -1, false, false, false, false, false,
               false);
  ASSERT(smith_waterman_align_checked("CA", "GA", 2, 2, &scoring, sw, &e)
         == ALN_ERR_UNSUPPORTED && e.status == ALN_ERR_UNSUPPORTED);
  ASSERT(!smith_waterman_fetch_checked(sw, aln, &e));
  scoring_init(&scoring, 2, -3, -1, 0, false, false, false, false, false,
               false);
  ASSERT(smith_waterman_align_checked("CA", "GA", 2, 2, &scoring, sw, &e)
         == ALN_OK);

  // Printing the matrices of a compact alignment (there are none) says so
  // instead of reading through NULL
  scoring_init(&scoring, 1, -1, -2, -1, false, false, false, false, false,
               false);
  nw->compact = true;
  ASSERT(needleman_wunsch_align_checked("ACGT", "AGT", 4, 3, &scoring, nw,
                                        aln, &e) == ALN_OK);
  if(nw->traced) {
    char *text = NULL;
    size_t len = 0;
    FILE *out = open_memstream(&text, &len);
    alignment_fprint_matrices(out, nw);
    fclose(out);
    ASSERT(strstr(text, "no score matrices") != NULL);
    free(text);
  }
  nw->compact = false;

  // A score that would reach the matrices' floor exactly is refused as an
  // overflow: at the floor a cell cannot be told from an unreachable one
  scoring_init(&scoring, 1, -387483648, -530000000, -350000000, false, false,
               false, false, false, false);
  ASSERT(needleman_wunsch_align_checked("AA", "C", 2, 1, &scoring, nw, aln,
                                        &e) == ALN_ERR_OVERFLOW);

  // The messages, as the aligners print them before exiting
  char *text = NULL;
  size_t len = 0;
  FILE *out = open_memstream(&text, &len);
  e.status = ALN_ERR_SCORING; e.a = 'g'; e.b = 'a';
  alignment_fprint_error(out, &e);
  fclose(out);
  ASSERT(strcmp(text, "Error: Unknown character pair (g,a) and "
                      "match/mismatch have not been set\n") == 0);
  free(text);

  needleman_wunsch_free(nw);
  smith_waterman_free(sw);
  alignment_free(aln);
}

/* Out of memory: each allocation of an alignment fails in turn */

enum { OOM_NW, OOM_NW_COMPACT, OOM_SW, OOM_SW_COMPACT, OOM_HB, OOM_HB_SCORE,
       OOM_KINDS };

typedef struct {
  const char *a, *b;
  const alignment_options_t *opts;
  const scoring_t *scoring;
  nw_aligner_t *nw;
  sw_aligner_t *sw;
  alignment_t *aln;
} oom_case_t;

// Append an alignment's score and strings to sig
static void oom_sig_add(char *sig, size_t cap, const alignment_t *aln)
{
  size_t len = strlen(sig);
  snprintf(sig + len, cap - len, "%i %s %s;", aln->score, aln->result_a,
           aln->result_b);
}

// One alignment of the kind given (Smith-Waterman: its first few hits),
// described in sig
static aln_status_t oom_run(oom_case_t *c, int kind, char *sig, size_t cap)
{
  const size_t la = strlen(c->a), lb = strlen(c->b);
  aln_error_t e;
  aln_status_t status;
  score_t score = 0;
  sig[0] = '\0';
  switch(kind) {
    case OOM_NW:
    case OOM_NW_COMPACT:
      c->nw->options = c->opts;
      c->nw->compact = kind == OOM_NW_COMPACT;
      status = needleman_wunsch_align_checked(c->a, c->b, la, lb, c->scoring,
                                              c->nw, c->aln, &e);
      if(status == ALN_OK) oom_sig_add(sig, cap, c->aln);
      break;
    case OOM_SW:
    case OOM_SW_COMPACT:
      smith_waterman_set_options(c->sw, c->opts);
      smith_waterman_set_compact(c->sw, kind == OOM_SW_COMPACT);
      status = smith_waterman_align_checked(c->a, c->b, la, lb, c->scoring,
                                            c->sw, &e);
      int hits = 0;
      while(status == ALN_OK && hits++ < 4 &&
            smith_waterman_fetch_checked(c->sw, c->aln, &e))
        oom_sig_add(sig, cap, c->aln);
      if(status == ALN_OK) status = e.status;
      // After an error: no hits, and no reads of what was not filled
      if(status != ALN_OK) {
        aln_error_t e2;
        ASSERT(!smith_waterman_fetch_checked(c->sw, c->aln, &e2) &&
               e2.status == ALN_OK);
      }
      break;
    case OOM_HB:
      status = hirschberg_align_checked(c->a, c->b, la, lb, c->scoring,
                                        c->opts, c->aln, &e);
      if(status == ALN_OK) oom_sig_add(sig, cap, c->aln);
      break;
    default:
      status = hirschberg_score_checked(c->a, c->b, la, lb, c->scoring,
                                        c->opts, &score, &e);
      if(status == ALN_OK) snprintf(sig, cap, "%i", score);
  }
  ASSERT(e.status == status);
  return status;
}

void test_out_of_memory()
{
  enum { CAP = 1 << 14 };
  static char a[161], b[321], ref[CAP], sig[CAP];
  size_t i;
  // b: a with about one letter in 25 changed, and no gaps, so the narrow
  // bands forced below always hold the best path (a run deletion could
  // put it outside them, and the test then could not tell a band that
  // was not used from one whose allocations were made to fail)
  test_rand_seq(a, 160, "acgt");
  for(i = 0; i < 160; i++) {
    b[i] = a[i];
    if(test_rand(25) == 0) b[i] = "acgt"[(strchr("acgt", a[i]) - "acgt" + 1 + test_rand(3)) % 4];
  }
  b[160] = '\0';
  // Blocks on several threads, the Smith-Waterman rectangle, and bands
  // narrow enough to pay on these short sequences
  alignment_options_t many;
  alignment_options_default(&many);
  many.max_threads = many.hb_max_threads = 4;
  many.block_min_cells = many.hb_block_min_cells = 1;
  many.block_min_cols = many.hb_block_min_cols = 8;
  many.sw_rect_min_cells = many.sw_strip_min_cells = 0;
  many.band_strip_rows = many.hb_band_strip_rows = 8;
  many.band_min_width = many.hb_band_min_width = 4;
  many.hb_reuse_min_rows = 2; // rows reused down the recursion
  // Bands with absolute lanes (with difference-encoded ones the whole fill
  // is cheap enough that a band this wide would not pay at some lane
  // counts), and the whole fill with difference-encoded lanes on blocks
  many.diff = false;
  many.band_keep_max_cells = 1; // bands kept in strips
  alignment_options_t whole = many;
  whole.band = whole.hb_band = false;
  whole.diff = true;
  const alignment_options_t *sets[] = {NULL, &many, &whole};
  scoring_t scoring;
  scoring_system_default(&scoring);
  oom_case_t c = {a, b, NULL, &scoring, NULL, NULL, NULL};
  int kind, o;
  size_t n;

  for(kind = 0; kind < OOM_KINDS; kind++) {
    for(o = 0; o < 3; o++) {
      c.opts = sets[o];
      for(n = 0; ; n++) {
        // New aligners each time, so every allocation is made again
        c.nw = needleman_wunsch_new();
        c.sw = smith_waterman_new();
        c.aln = alignment_create(16);
        alignment_alloc_fail(n);
        size_t banded = hirschberg_banded_passes;
        hirschberg_check_striped = n == 0; // counts banded passes
        aln_status_t status = oom_run(&c, kind, n ? sig : ref, CAP);
        hirschberg_check_striped = false;
        if(n == 0 && o == 1) {
          // The forced bands are used, so their allocations fail too
          // (bands need the striped fill)
#ifdef ALN_STRIPED
          if(kind == OOM_NW_COMPACT) ASSERT(c.nw->trace_banded);
#endif
          if(kind == OOM_HB && hirschberg_simd_available(NULL))
            ASSERT(hirschberg_banded_passes > banded);
        }
        bool failed = n > 0 && alignment_alloc_fail_pending() == 0;
        alignment_alloc_fail(0);
        bool done = n > 0 && !failed;
        if(n == 0) {
          ASSERT(status == ALN_OK); // the result to match
        }
        else if(!failed) {
          // Fewer than n allocations: the whole alignment ran
          ASSERT(status == ALN_OK && strcmp(sig, ref) == 0);
        }
        else {
          // A failure that only stops a faster path gives the same result
          ASSERT(status == ALN_ERR_NOMEM ||
                 (status == ALN_OK && strcmp(sig, ref) == 0));
          // and the aligner works as before afterwards
          ASSERT(oom_run(&c, kind, sig, CAP) == ALN_OK &&
                 strcmp(sig, ref) == 0);
        }
        needleman_wunsch_free(c.nw);
        smith_waterman_free(c.sw);
        alignment_free(c.aln);
        if(done) break;
      }
      ASSERT(n > 1); // something was allocated
    }
  }

  // The message
  char *text = NULL;
  size_t len = 0;
  FILE *out = open_memstream(&text, &len);
  aln_error_t e = {ALN_ERR_NOMEM, 0, 0};
  alignment_fprint_error(out, &e);
  fclose(out);
  ASSERT(strcmp(text, "Error: Out of memory\n") == 0);
  free(text);
}

void test_errors(void)
{
  SUITE_START("Errors and out of memory");
  RUN(test_checked());
  RUN(test_out_of_memory());
  SUITE_END();
}
