/*
 tests/test_sw.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 Smith-Waterman: hits, reuse and the index of later hits
 */

#include "test_common.h"
#include "smith_waterman_internal.h" // sw_history_t, for the index-size test


void sw_test_no_gaps_smith_waterman()
{
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *result = alignment_create(256);

  const char* seq_a = "gacag";
  const char* seq_b = "tgaagt";

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;

  bool no_start_gap_penalty = false, no_end_gap_penalty = false;
  bool no_gaps_in_a = true, no_gaps_in_b = true;
  bool no_mismatches = false, case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b,
               no_mismatches, case_sensitive);

  smith_waterman_align(seq_a, seq_b, &scoring, sw);

  smith_waterman_fetch(sw, result);
  ASSERT(strcmp(result->result_a, "ga") == 0 &&
         strcmp(result->result_b, "ga") == 0);

  smith_waterman_fetch(sw, result);
  ASSERT(strcmp(result->result_a, "ag") == 0 &&
         strcmp(result->result_b, "ag") == 0);

  alignment_free(result);
  smith_waterman_free(sw);
}


// Fetch every hit, returning the number of hits and the sum of their scores
static size_t sw_fetch_all(sw_aligner_t *sw, alignment_t *result, long *score_sum)
{
  size_t n = 0;
  *score_sum = 0;
  while(smith_waterman_fetch(sw, result)) { n++; *score_sum += result->score; }
  return n;
}

/* Re-using an aligner must give the same hits as a fresh one */
void sw_test_reuse_aligner()
{
  sw_aligner_t *reused = smith_waterman_new(), *fresh = smith_waterman_new();
  alignment_t *r1 = alignment_create(256), *r2 = alignment_create(256);

  scoring_t scoring;
  scoring_system_default(&scoring);

  // Leave visited-cell marks behind in the reused aligner
  long sum;
  smith_waterman_align("aaaaaaaa", "aaaaaaaa", &scoring, reused);
  ASSERT(sw_fetch_all(reused, r1, &sum) > 0);

  smith_waterman_align("cccccccc", "cccccccc", &scoring, reused);
  smith_waterman_align("cccccccc", "cccccccc", &scoring, fresh);

  // Compare hit by hit
  size_t n = 0;
  int got1, got2;
  while(1) {
    got1 = smith_waterman_fetch(reused, r1);
    got2 = smith_waterman_fetch(fresh, r2);
    ASSERT(got1 == got2);
    if(!got1 || !got2) break;
    ASSERT(r1->score == r2->score && r1->pos_a == r2->pos_a &&
           r1->pos_b == r2->pos_b && strcmp(r1->result_a, r2->result_a) == 0);
    n++;
  }
  ASSERT(n > 0);

  alignment_free(r1);
  alignment_free(r2);
  smith_waterman_free(reused);
  smith_waterman_free(fresh);
}

// The index of later hits holds the cells scoring at least min_score, not
// every cell of the matrix (a 50 kb pair would need 34 GB for that)
void sw_test_index_size()
{
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *aln = alignment_create(256);
  scoring_t scoring;
  scoring_system_default(&scoring);
  static char a[401], b[401];
  test_rand_seq(a, 400, "acgt");
  test_rand_seq(b, 400, "acgt");
  memcpy(b + 100, a + 50, 20); // one clear hit, few other cells score >= 20
  smith_waterman_set_min_score(sw, 20);
  smith_waterman_align(a, b, &scoring, sw);
  ASSERT(smith_waterman_fetch(sw, aln));
  smith_waterman_fetch(sw, aln); // builds the index
  const sw_history_t *hist = &sw->history;
  ASSERT(hist->hits_capacity > 0);
  ASSERT(hist->hits_capacity < 401 * 401 / 8);
  alignment_free(aln);
  smith_waterman_free(sw);
}

// A strips fill of a region with no rows is nothing to do, not a division
// by zero (the strip height is clamped to the rows)
void sw_test_strips_empty()
{
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *aln = alignment_create(16);
  alignment_options_t opts;
  scoring_t scoring;
  scoring_system_default(&scoring);
  alignment_options_default(&opts);
  opts.sw_strip_min_cells = 0; // every compact fill through the strips
  opts.sw_rect_min_cells = 0;
  smith_waterman_set_options(sw, &opts);
  smith_waterman_set_compact(sw, true);
  aln_error_t e;
  ASSERT(smith_waterman_align_checked("ACGT", "", 4, 0, &scoring, sw, &e)
         == ALN_OK);
  ASSERT(!smith_waterman_fetch_checked(sw, aln, &e) && e.status == ALN_OK);
  ASSERT(smith_waterman_align_checked("", "ACGT", 0, 4, &scoring, sw, &e)
         == ALN_OK);
  ASSERT(!smith_waterman_fetch_checked(sw, aln, &e) && e.status == ALN_OK);
  alignment_free(aln);
  smith_waterman_free(sw);
}

void test_sw(void)
{
  SUITE_START("Smith-Waterman");
  RUN(sw_test_no_gaps_smith_waterman());
  RUN(sw_test_reuse_aligner());
  RUN(sw_test_index_size());
  RUN(sw_test_strips_empty());
  SUITE_END();
}
