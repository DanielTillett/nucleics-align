/*
 tests/test_hirschberg.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 Linear-space global alignment scores and aligns as Needleman-Wunsch
 */

#include "test_common.h"


/* Linear-space alignment must score the same as Needleman-Wunsch */

static void hb_check(const char *a, const char *b, const scoring_t *scoring,
                     nw_aligner_t *nw, alignment_t *nw_aln,
                     alignment_t *hb_aln)
{
  needleman_wunsch_align(a, b, scoring, nw, nw_aln);
  ASSERT(hirschberg_align(a, b, scoring, hb_aln));
  ASSERT(hb_aln->score == nw_aln->score);
  score_t score_only = 0;
  ASSERT(hirschberg_score(a, b, scoring, &score_only));
  ASSERT(score_only == nw_aln->score);
  ASSERT(strlen(hb_aln->result_a) == hb_aln->length &&
         strlen(hb_aln->result_b) == hb_aln->length);
  ASSERT(test_ungapped_equals(hb_aln->result_a, a));
  ASSERT(test_ungapped_equals(hb_aln->result_b, b));
  // The scorer agrees with Needleman-Wunsch on its own alignment
  ASSERT(hirschberg_score_alignment(nw_aln->result_a, nw_aln->result_b,
                                    scoring) == nw_aln->score);
}

void hb_test_rand_dna()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *nw_aln = alignment_create(256), *hb_aln = alignment_create(256);
  char seqa[61], seqb[61];
  scoring_t scoring;
  size_t i;

  for(i = 0; i < 3000; i++)
  {
    // Random affine penalties, including gap_open == 0 (linear gaps) and
    // mismatches costly enough that deleting both characters is better
    int match = 1 + (int)test_rand(5), mismatch = -(int)test_rand(16);
    int gap_open = -(int)test_rand(11), gap_extend = -1 - (int)test_rand(5);
    scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
                 false, false, false, false, false, false);
    if(test_rand(2)) scoring_add_wildcard(&scoring, 'n', 0);

    // Short sequences hit the base cases more often
    unsigned int max_len = i % 4 ? 61 : 9;
    test_rand_seq(seqa, test_rand(max_len), i % 3 ? "acgt" : "acgtn");
    test_rand_seq(seqb, test_rand(max_len), i % 2 ? "acgt" : "ACgtN");
    hb_check(seqa, seqb, &scoring, nw, nw_aln, hb_aln);
  }

  alignment_free(nw_aln);
  alignment_free(hb_aln);
  needleman_wunsch_free(nw);
}

void hb_test_rand_protein()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *nw_aln = alignment_create(256), *hb_aln = alignment_create(256);
  char seqa[201], seqb[201];
  scoring_t scoring;
  size_t i;

  scoring_system_BLOSUM62(&scoring);

  for(i = 0; i < 300; i++)
  {
    test_rand_seq(seqa, test_rand(201), "ARNDCQEGHILKMFPSTWYV");
    test_rand_seq(seqb, test_rand(201), "ARNDCQEGHILKMFPSTWYV");
    hb_check(seqa, seqb, &scoring, nw, nw_aln, hb_aln);
  }

  alignment_free(nw_aln);
  alignment_free(hb_aln);
  needleman_wunsch_free(nw);
}

void hb_test_edge_cases()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *nw_aln = alignment_create(256), *hb_aln = alignment_create(256);
  scoring_t scoring;
  scoring_init(&scoring, 1, -2, -4, -1, false, false, false, false,
               false, true);

  hb_check("a", "a", &scoring, nw, nw_aln, hb_aln);
  hb_check("a", "c", &scoring, nw, nw_aln, hb_aln);
  hb_check("aaaaacg", "acgt", &scoring, nw, nw_aln, hb_aln);
  hb_check("acgt", "aaaaacg", &scoring, nw, nw_aln, hb_aln);
  hb_check("aaaacccccgggg", "aaaattcccccc", &scoring, nw, nw_aln, hb_aln);
  // Long gaps crossing the split rows
  hb_check("acgtttttttttttttttttacgt", "acgtacgt", &scoring, nw, nw_aln, hb_aln);
  hb_check("acgtacgt", "acgtttttttttttttttttacgt", &scoring, nw, nw_aln, hb_aln);

  // Empty sequences
  ASSERT(hirschberg_align("", "acg", &scoring, hb_aln));
  ASSERT(strcmp(hb_aln->result_a, "---") == 0 &&
         strcmp(hb_aln->result_b, "acg") == 0 && hb_aln->score == -7);
  ASSERT(hirschberg_align("acg", "", &scoring, hb_aln));
  ASSERT(strcmp(hb_aln->result_a, "acg") == 0 &&
         strcmp(hb_aln->result_b, "---") == 0 && hb_aln->score == -7);
  ASSERT(hirschberg_align("", "", &scoring, hb_aln));
  ASSERT(hb_aln->length == 0 && hb_aln->score == 0);

  // Options that need the full matrices are refused
  scoring.no_end_gap_penalty = true;
  ASSERT(!hirschberg_align("acgt", "acg", &scoring, hb_aln));
  scoring.no_end_gap_penalty = false;
  scoring.no_gaps_in_longest = true;
  ASSERT(!hirschberg_align("acgt", "acg", &scoring, hb_aln));

  // Penalties large enough for a score to overflow are refused, not
  // answered wrongly (these used to give 1894967296 and 2094967296)
  score_t big_score;
  scoring_init(&scoring, 1, -1, -800000000, 0, false, false, false, false,
               false, false);
  ASSERT(!hirschberg_align("AA", "AA", &scoring, hb_aln));
  ASSERT(!hirschberg_score("AA", "AA", &scoring, &big_score));
  scoring_init(&scoring, 1, -1, -1100000000, 0, false, false, false, false,
               false, false);
  ASSERT(!hirschberg_align("A", "A", &scoring, hb_aln));
  ASSERT(!hirschberg_score("A", "A", &scoring, &big_score));

  alignment_free(nw_aln);
  alignment_free(hb_aln);
  needleman_wunsch_free(nw);
}

// Every pair of sequences over {a,c} up to length 6, under penalties that
// favour matches, long gaps, or deleting rather than mismatching
void hb_test_exhaustive_small()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *nw_aln = alignment_create(256), *hb_aln = alignment_create(256);
  const int penalties[][4] = {{1, -1, -1, -1}, {1, -15, -10, -1},
                              {2, -3, 0, -2}, {1, -5, -10, -5},
                              {3, 0, -6, -1}, {1, -20, -2, -1}};
  char seqa[7], seqb[7];
  scoring_t scoring;
  size_t p, la, lb, xa, xb, k;

  for(p = 0; p < sizeof(penalties)/sizeof(penalties[0]); p++)
  {
    scoring_init(&scoring, penalties[p][0], penalties[p][1],
                 penalties[p][2], penalties[p][3],
                 false, false, false, false, false, true);

    // xa, xb enumerate sequences as bit patterns (a=0, c=1)
    for(la = 0; la <= 6; la++)
      for(xa = 0; xa < (1UL << la); xa++)
        for(lb = 0; lb <= 6; lb++)
          for(xb = 0; xb < (1UL << lb); xb++)
          {
            for(k = 0; k < la; k++) seqa[k] = (xa >> k) & 1 ? 'c' : 'a';
            for(k = 0; k < lb; k++) seqb[k] = (xb >> k) & 1 ? 'c' : 'a';
            seqa[la] = seqb[lb] = '\0';
            hb_check(seqa, seqb, &scoring, nw, nw_aln, hb_aln);
          }
  }

  alignment_free(nw_aln);
  alignment_free(hb_aln);
  needleman_wunsch_free(nw);
}

void test_hirschberg(void)
{
  SUITE_START("Linear-space global alignment");
  RUN(hb_test_edge_cases());
  RUN(hb_test_exhaustive_small());
  RUN(hb_test_rand_dna());
  RUN(hb_test_rand_protein());
  SUITE_END();
}
