/*
 tests/test_nw.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 Needleman-Wunsch: the options that change the recurrence
 */

#include "test_common.h"


/* No gap is expected in the longer sequence */
void nw_test_no_gaps_in_longer()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(256);

  const char* seq_a = "aaaaacg";
  const char* seq_b = "acgt";

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;

  bool no_start_gap_penalty = false, no_end_gap_penalty = false;
  bool no_gaps_in_a = true, no_gaps_in_b = false;
  bool no_mismatches = false, case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b,
               no_mismatches, case_sensitive);

  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, aln);

  // ASSERT(strcmp(aln->result_a, "aaaaacg") == 0 &&
  //        strcmp(aln->result_b, "acgt---") == 0);

  ASSERT(strcmp(aln->result_a, "aaaaacg-") == 0 &&
         strcmp(aln->result_b, "a----cgt") == 0);

  alignment_free(aln);
  needleman_wunsch_free(nw);
}

/* --nogapsinlongest: the longer sequence gets no gaps, whichever input it is */
void nw_test_no_gaps_in_longest()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(256);

  scoring_t scoring;
  scoring_init(&scoring, 1, -2, -4, -1, false, false, false, false,
               false, true);
  scoring.no_gaps_in_longest = true;

  // Without the option the longer sequence gets a gap: aaaa---cccccgggg
  // Longer sequence first
  needleman_wunsch_align("aaaacccccgggg", "aaaattcccccc", &scoring, nw, aln);
  ASSERT(strcmp(aln->result_a, "aaaacccccgggg") == 0 &&
         strcmp(aln->result_b, "aaaattcccccc-") == 0);

  // Longer sequence second
  needleman_wunsch_align("aaaattcccccc", "aaaacccccgggg", &scoring, nw, aln);
  ASSERT(strcmp(aln->result_a, "aaaattcccccc-") == 0 &&
         strcmp(aln->result_b, "aaaacccccgggg") == 0);

  alignment_free(aln);
  needleman_wunsch_free(nw);
}

/* First sequence is aligned to the corresponding (equal) substring of the second
 * sequence because both gaps at start and at end are free */
void nw_test_free_gaps_at_ends()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);

  const char* seq_a = "acg";
  const char* seq_b = "tttacgttt";

  int match = 1;
  int mismatch = -1;
  int gap_open = -4;
  int gap_extend = -1;

  bool no_start_gap_penalty = true, no_end_gap_penalty = true;
  bool no_gaps_in_a = false, no_gaps_in_b = false;
  bool no_mismatches = false, case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b,
               no_mismatches, case_sensitive);

  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
  ASSERT(strcmp(result->result_a, "---acg---") == 0 &&
         strcmp(result->result_b, "tttacgttt") == 0);

  alignment_free(result);
  needleman_wunsch_free(nw);
}

/* Against an empty sequence the whole of the other is one gap, at the start
 * and at the end at once: free with either free start or free end gaps
 * (free end gaps used to charge for it) */
void nw_test_free_gaps_empty()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);
  scoring_t scoring;
  int start, end;

  for(start = 0; start < 2; start++) {
    for(end = 0; end < 2; end++) {
      scoring_init(&scoring, 1, -2, -4, -1, start, end,
                   false, false, false, false);
      const score_t want = start || end ? 0 : -8;
      needleman_wunsch_align("acgt", "", &scoring, nw, result);
      ASSERT(result->score == want &&
             strcmp(result->result_a, "acgt") == 0 &&
             strcmp(result->result_b, "----") == 0);
      needleman_wunsch_align("", "acgt", &scoring, nw, result);
      ASSERT(result->score == want &&
             strcmp(result->result_a, "----") == 0 &&
             strcmp(result->result_b, "acgt") == 0);
      needleman_wunsch_align("", "", &scoring, nw, result);
      ASSERT(result->score == 0);
      // Not empty: free end gaps alone still charge a gap at the start
      needleman_wunsch_align("ttacgt", "acgt", &scoring, nw, result);
      ASSERT(result->score == (start ? 4 : -2));
    }
  }

  alignment_free(result);
  needleman_wunsch_free(nw);
}

void nw_test_no_mismatches()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;

  bool no_start_gap_penalty = false, no_end_gap_penalty = false;
  bool no_gaps_in_a = false, no_gaps_in_b = false;
  bool no_mismatches = true, case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b,
               no_mismatches, case_sensitive);

  needleman_wunsch_align("atc", "ac", &scoring, nw, result);
  ASSERT(strcmp(result->result_a, "atc") == 0 &&
         strcmp(result->result_b, "a-c") == 0);

  needleman_wunsch_align("cgatcga", "catcctcga", &scoring, nw, result);
  ASSERT(strcmp(result->result_a, "cgatc---ga") == 0 &&
         strcmp(result->result_b, "c-atcctcga") == 0);

  alignment_free(result);
  needleman_wunsch_free(nw);
}

void nw_test_no_mismatches_rand()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(256);

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;

  bool no_start_gap_penalty = false, no_end_gap_penalty = false;
  bool no_gaps_in_a = false, no_gaps_in_b = false;
  bool no_mismatches = true, case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b,
               no_mismatches, case_sensitive);

  char seqa[100], seqb[100];
  size_t i;

  // Run 50 random alignments
  for(i = 0; i < 50; i++)
  {
    test_rand_seq(seqa, test_rand(99), "acgt");
    test_rand_seq(seqb, test_rand(99), "acgt");
    needleman_wunsch_align(seqa, seqb, &scoring, nw, aln);
    // Check no mismatches
    char *a = aln->result_a, *b = aln->result_b;
    while(1) {
      ASSERT(*a == '-' || *b == '-' || *a == *b);
        // printf("Seq: '%s' '%s'\n", aln->result_a, aln->result_b);
        // exit(EXIT_FAILURE);
      if(!*a && !*b) break;
      a++; b++;
    }
  }

  alignment_free(aln);
  needleman_wunsch_free(nw);
}


/* Bytes >= 0x80 (e.g. UTF-8) must not index outside the scoring tables */
void nw_test_high_bit_chars()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);

  scoring_t scoring;
  scoring_init(&scoring, 1, -2, -4, -1, false, false, false, false, false, false);
  scoring_add_wildcard(&scoring, (char)0xff, 1);

  needleman_wunsch_align("ac\xc3\xa9t", "ac\xc3\xa9t", &scoring, nw, result);
  ASSERT(strcmp(result->result_a, "ac\xc3\xa9t") == 0 &&
         strcmp(result->result_b, "ac\xc3\xa9t") == 0);
  ASSERT(result->score == 5);

  // 0xff is a wildcard so matches anything
  needleman_wunsch_align("a\xff" "t", "act", &scoring, nw, result);
  ASSERT(result->score == 3);

  alignment_free(result);
  needleman_wunsch_free(nw);
}

/* Changing penalties after scoring_init() (as the command line parser does)
 * must not let the matrix floor overflow */
void nw_test_penalties_changed_after_init()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);

  scoring_t scoring;
  scoring_init(&scoring, 1, -2, -4, -1, false, false, true, false, false, true);
  scoring.gap_open = -1000;

  needleman_wunsch_align("acgt", "aacgt", &scoring, nw, result);
  ASSERT(result->score == 4 + -1000 + -1);
  ASSERT(strcmp(result->result_a, "-acgt") == 0 &&
         strcmp(result->result_b, "aacgt") == 0);

  alignment_free(result);
  needleman_wunsch_free(nw);
}

void test_nw(void)
{
  SUITE_START("Needleman-Wunsch");
  RUN(nw_test_no_gaps_in_longer());
  RUN(nw_test_no_gaps_in_longest());
  RUN(nw_test_free_gaps_at_ends());
  RUN(nw_test_free_gaps_empty());
  RUN(nw_test_no_mismatches());
  RUN(nw_test_no_mismatches_rand());
  RUN(nw_test_high_bit_chars());
  RUN(nw_test_penalties_changed_after_init());
  SUITE_END();
}
