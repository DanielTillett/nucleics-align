/*
 tests/tests.c
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Feb 2015
 */

// request decent POSIX version (mkstemp, getpid)
#define _XOPEN_SOURCE 700

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h> // needed for rand()
#include <unistd.h>  // need for getpid() for getting setting rand number
#include <pthread.h>


#include "needleman_wunsch.h"
#include "smith_waterman.h"
#include "hirschberg.h"
#include "alignment_scoring_load.h"
#include "alignment_simd.h" // ALN_STRIPED

//
// Tests
//
const char *suite_name = NULL;
char suite_pass = 0;
int suites_run = 0, suites_failed = 0, suites_empty = 0;
int tests_in_suite = 0, tests_run = 0, tests_failed = 0;

#define QUOTE(str) #str
#define ASSERT(x) {tests_run++; tests_in_suite++; if(!(x)) \
  { fprintf(stderr, "failed assert [%s:%i] %s\n", __FILE__, __LINE__, QUOTE(x)); \
    suite_pass = 0; tests_failed++; }}

void SUITE_START(const char *name)
{
  suite_pass = 1;
  suite_name = name;
  suites_run++;
  tests_in_suite = 0;
}

void SUITE_END()
{
  printf("Testing %s ", suite_name);
  size_t suite_i;
  for(suite_i = strlen(suite_name); suite_i < 80-8-5; suite_i++) printf(".");
  printf("%s\n", suite_pass ? " pass" : " fail");
  if(!suite_pass) suites_failed++;
  if(!tests_in_suite) suites_empty++;
}

//
// Useful MACROs

#define die(fmt,...) do { \
  fprintf(stderr, "[%s:%i] Error: %s() "fmt"\n", __FILE__, __LINE__, __func__, ##__VA_ARGS__); \
  exit(EXIT_FAILURE); \
} while(0);

#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define MIN(x,y) ((x) <= (y) ? (x) : (y))

//

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

// size is buffer size including NUL byte
static char* make_rand_seq(char *seq, size_t size)
{
  if(size == 0) { return seq; }
  size_t i, len = rand() % (size-1);
  const char bases[] = "acgt";
  for(i = 0; i < len; i++) seq[i] = bases[rand() & 3];
  seq[len] = '\0';
  return seq;
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
    make_rand_seq(seqa, sizeof(seqa));
    make_rand_seq(seqb, sizeof(seqb));
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

/* Linear-space alignment must score the same as Needleman-Wunsch */

// Deterministic generator so failures can be reproduced
static unsigned long hb_rand_state = 1;
static unsigned int hb_rand(unsigned int n)
{
  hb_rand_state = hb_rand_state * 6364136223846793005UL + 1442695040888963407UL;
  return (unsigned int)((hb_rand_state >> 33) % n);
}

static void hb_rand_seq(char *seq, size_t len, const char *alphabet)
{
  size_t i, alen = strlen(alphabet);
  for(i = 0; i < len; i++) seq[i] = alphabet[hb_rand((unsigned int)alen)];
  seq[len] = '\0';
}

// Removing gaps from an alignment string gives back the sequence
static bool hb_ungapped_equals(const char *aln, const char *seq)
{
  for(; *aln != '\0'; aln++) {
    if(*aln == '-') continue;
    if(*aln != *seq) return false;
    seq++;
  }
  return *seq == '\0';
}

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
  ASSERT(hb_ungapped_equals(hb_aln->result_a, a));
  ASSERT(hb_ungapped_equals(hb_aln->result_b, b));
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
    int match = 1 + (int)hb_rand(5), mismatch = -(int)hb_rand(16);
    int gap_open = -(int)hb_rand(11), gap_extend = -1 - (int)hb_rand(5);
    scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
                 false, false, false, false, false, false);
    if(hb_rand(2)) scoring_add_wildcard(&scoring, 'n', 0);

    // Short sequences hit the base cases more often
    unsigned int max_len = i % 4 ? 61 : 9;
    hb_rand_seq(seqa, hb_rand(max_len), i % 3 ? "acgt" : "acgtn");
    hb_rand_seq(seqb, hb_rand(max_len), i % 2 ? "acgt" : "ACgtN");
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
    hb_rand_seq(seqa, hb_rand(201), "ARNDCQEGHILKMFPSTWYV");
    hb_rand_seq(seqb, hb_rand(201), "ARNDCQEGHILKMFPSTWYV");
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

void test_hirschberg()
{
  SUITE_START("Linear-space global alignment");

  hb_test_edge_cases();
  hb_test_exhaustive_small();
  hb_test_rand_dna();
  hb_test_rand_protein();

  SUITE_END();
}

/* The striped SIMD fill must give exactly the scalar matrices */
static void simd_check(const char *a, const char *b, const scoring_t *scoring,
                       char is_sw, aligner_t *scalar, aligner_t *striped)
{
  alignment_striped_enabled = false;
  aligner_align(scalar, a, b, strlen(a), strlen(b), scoring, is_sw);
  alignment_striped_enabled = true;
  aligner_align(striped, a, b, strlen(a), strlen(b), scoring, is_sw);

  size_t cells = scalar->score_width * scalar->score_height;
  ASSERT(memcmp(scalar->match_scores, striped->match_scores,
                cells * sizeof(score_t)) == 0);
  ASSERT(memcmp(scalar->gap_a_scores, striped->gap_a_scores,
                cells * sizeof(score_t)) == 0);
  ASSERT(memcmp(scalar->gap_b_scores, striped->gap_b_scores,
                cells * sizeof(score_t)) == 0);
}

void simd_test_matches_scalar()
{
  aligner_t scalar, striped;
  aligner_init(&scalar);
  aligner_init(&striped);
  char seqa[301], seqb[301];
  scoring_t scoring;
  size_t i;

  for(i = 0; i < 3000; i++)
  {
    // gap_extend == 0 and gap_open == 0 included; short seq_a gives
    // segments of one or two cells, long gaps exercise the lazy F loop
    int match = 1 + (int)hb_rand(5), mismatch = -(int)hb_rand(16);
    int gap_open = -(int)hb_rand(11), gap_extend = -(int)hb_rand(6);
    scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
                 hb_rand(2), false, false, false, false, false);
    if(hb_rand(2)) scoring_add_wildcard(&scoring, 'n', 0);

    unsigned int max_a = i % 4 == 0 ? 20 : 301, max_b = i % 3 == 0 ? 20 : 301;
    hb_rand_seq(seqa, hb_rand(max_a), "acgtn");
    hb_rand_seq(seqb, hb_rand(max_b), "ACgtN");
    simd_check(seqa, seqb, &scoring, (char)hb_rand(2), &scalar, &striped);
  }

  scoring_system_BLOSUM62(&scoring);
  for(i = 0; i < 300; i++)
  {
    hb_rand_seq(seqa, hb_rand(301), "ARNDCQEGHILKMFPSTWYV");
    hb_rand_seq(seqb, hb_rand(301), "ARNDCQEGHILKMFPSTWYV");
    simd_check(seqa, seqb, &scoring, (char)hb_rand(2), &scalar, &striped);
  }

  // A long gap in seq_b crossing every segment boundary
  scoring_init(&scoring, 1, -2, -4, -1, false, false, false, false,
               false, true);
  simd_check("acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt",
             "acgtacgt", &scoring, 0, &scalar, &striped);

  // Mismatches far worse than extending a gap: in the first rows the gap
  // from column 0 wins in every cell, so the lazy F loop has to carry it
  // through whole segments into later lanes
  char poly_a[201];
  memset(poly_a, 'a', 200);
  poly_a[200] = '\0';
  for(i = 1; i <= 200; i += 7)
  {
    poly_a[i] = '\0';
    scoring_init(&scoring, 1, -20, -2, -1, false, false, false, false,
                 false, true);
    simd_check(poly_a, "c", &scoring, 0, &scalar, &striped);
    simd_check(poly_a, "cc", &scoring, 0, &scalar, &striped);
    simd_check(poly_a, "cca", &scoring, 0, &scalar, &striped);
    scoring_init(&scoring, 1, -20, 0, -1, true, false, false, false,
                 false, true);
    simd_check(poly_a, "cgc", &scoring, 0, &scalar, &striped);
    poly_a[i] = 'a';
  }

  aligner_destroy(&scalar);
  aligner_destroy(&striped);
}

/* Positive or zero gap penalties are unusual but accepted */
void simd_test_odd_penalties()
{
  aligner_t scalar, striped;
  aligner_init(&scalar);
  aligner_init(&striped);
  char seqa[61], seqb[61];
  scoring_t scoring;
  size_t i;

  for(i = 0; i < 500; i++)
  {
    int gap_open = (int)hb_rand(11) - 5, gap_extend = (int)hb_rand(7) - 3;
    scoring_init(&scoring, 1 + (int)hb_rand(3), -(int)hb_rand(5),
                 gap_open, gap_extend, hb_rand(2), false, false, false,
                 false, false);
    hb_rand_seq(seqa, hb_rand(61), "acgt");
    hb_rand_seq(seqb, hb_rand(61), "acgt");
    simd_check(seqa, seqb, &scoring, (char)hb_rand(2), &scalar, &striped);
  }

  aligner_destroy(&scalar);
  aligner_destroy(&striped);
}

/* A substitution table that only defines pairs of a character from seq_a
 * with one from seq_b (no match/mismatch) must not ask for any other pair:
 * scoring_lookup() exits on an undefined pair */
void simd_test_directional_pairs()
{
  scoring_t scoring;
  scoring_init(&scoring, 0, 0, -4, -1, false, false, false, false,
               false, true);
  scoring.use_match_mismatch = 0;
  scoring_add_mutation(&scoring, 'a', 'x', 3);
  scoring_add_mutation(&scoring, 'a', 'y', -1);
  scoring_add_mutation(&scoring, 'c', 'x', -1);
  scoring_add_mutation(&scoring, 'c', 'y', 2);

  nw_aligner_t *nw = needleman_wunsch_new();
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *aln = alignment_create(256);

  needleman_wunsch_align("ac", "xy", &scoring, nw, aln);
  ASSERT(aln->score == 5);
  ASSERT(hirschberg_align("ac", "xy", &scoring, aln));
  ASSERT(aln->score == 5);
  smith_waterman_align("ac", "xy", &scoring, sw);
  ASSERT(smith_waterman_fetch(sw, aln) && aln->score == 5);

  alignment_free(aln);
  smith_waterman_free(sw);
  needleman_wunsch_free(nw);
}

/* Linear-space alignment must give the same alignment with or without SIMD */
static void simd_hb_check(const char *a, const char *b,
                          const scoring_t *scoring,
                          alignment_t *scalar, alignment_t *striped)
{
  if(scoring->gap_open > 0)
  {
    // Not supported in linear space: must be refused, not answered wrongly
    score_t score;
    ASSERT(!hirschberg_align(a, b, scoring, scalar));
    ASSERT(!hirschberg_score(a, b, scoring, &score));
    return;
  }

  alignment_striped_enabled = false;
  ASSERT(hirschberg_align(a, b, scoring, scalar));
  alignment_striped_enabled = true;
  hirschberg_check_striped = true;
  ASSERT(hirschberg_align(a, b, scoring, striped));
  hirschberg_check_striped = false;
  ASSERT(scalar->score == striped->score);
  score_t score_scalar = 0, score_striped = 0;
  alignment_striped_enabled = false;
  ASSERT(hirschberg_score(a, b, scoring, &score_scalar));
  alignment_striped_enabled = true;
  ASSERT(hirschberg_score(a, b, scoring, &score_striped));
  ASSERT(score_scalar == scalar->score && score_striped == scalar->score);
  ASSERT(strcmp(scalar->result_a, striped->result_a) == 0);
  ASSERT(strcmp(scalar->result_b, striped->result_b) == 0);
}

void simd_test_hirschberg()
{
  alignment_t *scalar = alignment_create(256), *striped = alignment_create(256);
  static char seqa[801], seqb[801];
  scoring_t scoring;
  size_t i;

  for(i = 0; i < 400; i++)
  {
    // Includes gap_open/gap_extend of 0, and positive values (scalar only)
    int gap_open = (int)hb_rand(13) - (i % 5 ? 11 : 1);
    int gap_extend = (int)hb_rand(6) - (i % 7 ? 5 : 1);
    scoring_init(&scoring, 1 + (int)hb_rand(5), -(int)hb_rand(16),
                 gap_open, gap_extend, false, false, false, false,
                 false, false);
    if(hb_rand(2)) scoring_add_wildcard(&scoring, 'n', 0);

    // Lengths either side of the thresholds for using SIMD
    unsigned int max_a = i % 3 ? 801 : 70, max_b = i % 4 ? 801 : 20;
    hb_rand_seq(seqa, hb_rand(max_a), "acgtn");
    hb_rand_seq(seqb, hb_rand(max_b), "ACgtN");
    simd_hb_check(seqa, seqb, &scoring, scalar, striped);
  }

  scoring_system_BLOSUM62(&scoring);
  for(i = 0; i < 100; i++)
  {
    hb_rand_seq(seqa, hb_rand(801), "ARNDCQEGHILKMFPSTWYV");
    hb_rand_seq(seqb, hb_rand(801), "ARNDCQEGHILKMFPSTWYV");
    simd_hb_check(seqa, seqb, &scoring, scalar, striped);
  }

  // Mismatches far worse than gaps: insertions carried through whole segments
  memset(seqa, 'a', 300); seqa[300] = 0;
  memset(seqb, 'c', 40); seqb[40] = 0;
  scoring_init(&scoring, 1, -20, -2, -1, false, false, false, false,
               false, true);
  simd_hb_check(seqa, seqb, &scoring, scalar, striped);
  simd_hb_check(seqb, seqa, &scoring, scalar, striped);

  // Huge gap_open with gap_extend 0: two gaps score below the SIMD
  // sentinel, which used to make the lazy F loop run forever, but nothing
  // overflows. Must fall back to the scalar passes
  memset(seqa, 'a', 32); seqa[32] = 0;
  memset(seqb, 'c', 32); seqb[32] = 0;
  scoring_init(&scoring, 0, 0, -300000000, 0, false, false, false, false,
               false, true);
  simd_hb_check(seqa, seqb, &scoring, scalar, striped);

  // Every SIMD row pass matched the scalar one exactly
  ASSERT(hirschberg_striped_mismatches == 0);

  alignment_free(scalar);
  alignment_free(striped);
}

/* Column-block threading of the linear-space passes: force small blocks so
 * short sequences are split across several threads */
void simd_test_hirschberg_blocks()
{
  size_t min_cells = hirschberg_block_min_cells;
  size_t min_cols = hirschberg_block_min_cols;
  int max_threads = hirschberg_max_threads;
  size_t mismatches = hirschberg_striped_mismatches;
  alignment_t *scalar = alignment_create(256), *striped = alignment_create(256);
  static char seqa[801], seqb[801];
  scoring_t scoring;
  const int threads[] = {2, 3, 5};
  size_t i, t;

  hirschberg_block_min_cells = 1;
  hirschberg_block_min_cols = 8;

  for(t = 0; t < sizeof(threads)/sizeof(threads[0]); t++)
  {
    hirschberg_max_threads = threads[t];
    for(i = 0; i < 80; i++)
    {
      scoring_init(&scoring, 1 + (int)hb_rand(5), -(int)hb_rand(16),
                   -(int)hb_rand(11), -(int)hb_rand(6), false, false, false,
                   false, false, false);
      hb_rand_seq(seqa, hb_rand(801), "acgt");
      hb_rand_seq(seqb, hb_rand(801), "acgt");
      simd_hb_check(seqa, seqb, &scoring, scalar, striped);
    }
    scoring_system_BLOSUM62(&scoring);
    for(i = 0; i < 20; i++)
    {
      hb_rand_seq(seqa, hb_rand(801), "ARNDCQEGHILKMFPSTWYV");
      hb_rand_seq(seqb, hb_rand(801), "ARNDCQEGHILKMFPSTWYV");
      simd_hb_check(seqa, seqb, &scoring, scalar, striped);
    }
    // Insertions carried through whole blocks
    memset(seqa, 'a', 300); seqa[300] = 0;
    memset(seqb, 'c', 200); seqb[200] = 0;
    scoring_init(&scoring, 1, -20, -2, -1, false, false, false, false,
                 false, true);
    simd_hb_check(seqa, seqb, &scoring, scalar, striped);
    simd_hb_check(seqb, seqa, &scoring, scalar, striped);
  }

  // Every block-threaded pass matched the scalar one exactly
  ASSERT(hirschberg_striped_mismatches == mismatches);

  hirschberg_block_min_cells = min_cells;
  hirschberg_block_min_cols = min_cols;
  hirschberg_max_threads = max_threads;
  alignment_free(scalar);
  alignment_free(striped);
}

/* A copy of in with about per_mille substitutions, a quarter as many
 * deletions and insertions, and now and then a deleted run of 20 */
static void hb_mutate(char *out, const char *in, unsigned int per_mille,
                      const char *alphabet)
{
  size_t i, o = 0, len = strlen(in), nalpha = strlen(alphabet);
  for(i = 0; i < len; i++)
  {
    unsigned int r = hb_rand(1000);
    if(r == 999) i += 20;
    else if(r < per_mille) out[o++] = alphabet[hb_rand(nalpha)];
    else if(r < per_mille + per_mille/4 + 1) continue;
    else if(r < per_mille + per_mille/2 + 2) {
      out[o++] = in[i];
      out[o++] = alphabet[hb_rand(nalpha)];
    }
    else out[o++] = in[i];
  }
  out[o] = '\0';
}

/* Banded linear-space passes: similar sequences, where the band pays, with
 * small strips and narrow first bands so that bands are often too narrow
 * and must be widened. Output must match the scalar code exactly */
void simd_test_hirschberg_band()
{
  size_t strip_rows = hirschberg_band_strip_rows;
  size_t min_width = hirschberg_band_min_width;
  size_t min_cells = hirschberg_block_min_cells;
  size_t min_cols = hirschberg_block_min_cols;
  int max_threads = hirschberg_max_threads;
  size_t banded = hirschberg_banded_passes;
  size_t mismatches = hirschberg_striped_mismatches;
  alignment_t *scalar = alignment_create(256), *striped = alignment_create(256);
  static char seqa[1201], seqb[2402];
  scoring_t scoring;
  const size_t strips[] = {0, 1, 5, 32};
  const unsigned int rates[] = {2, 10, 40, 120};
  size_t s, i;

  for(s = 0; s < sizeof(strips)/sizeof(strips[0]); s++)
  {
    hirschberg_band_strip_rows = strips[s];
    hirschberg_band_min_width = s % 2 ? 1 : 64;
    if(s == 3) {
      // Strips split into column blocks on several threads too
      hirschberg_block_min_cells = 1;
      hirschberg_block_min_cols = 8;
      hirschberg_max_threads = 3;
    }
    for(i = 0; i < 40; i++)
    {
      scoring_init(&scoring, 1 + (int)hb_rand(5), -(int)hb_rand(16),
                   -(int)hb_rand(11), -(int)hb_rand(6), false, false, false,
                   false, false, false);
      hb_rand_seq(seqa, 100 + hb_rand(1101), "acgt");
      hb_mutate(seqb, seqa, rates[i % 4], "acgt");
      simd_hb_check(seqa, seqb, &scoring, scalar, striped);
      simd_hb_check(seqb, seqa, &scoring, scalar, striped);
    }
    scoring_system_BLOSUM62(&scoring);
    for(i = 0; i < 8; i++)
    {
      hb_rand_seq(seqa, 100 + hb_rand(1101), "ARNDCQEGHILKMFPSTWYV");
      hb_mutate(seqb, seqa, rates[i % 4], "ARNDCQEGHILKMFPSTWYV");
      simd_hb_check(seqa, seqb, &scoring, scalar, striped);
    }
  }

  // The band was used (it needs the SIMD passes), and unbanded passes still
  // matched the scalar ones
  ASSERT(hirschberg_banded_passes > banded || !hirschberg_simd_available());
  ASSERT(hirschberg_striped_mismatches == mismatches);

  hirschberg_band_strip_rows = strip_rows;
  hirschberg_band_min_width = min_width;
  hirschberg_block_min_cells = min_cells;
  hirschberg_block_min_cols = min_cols;
  hirschberg_max_threads = max_threads;
  alignment_free(scalar);
  alignment_free(striped);
}

/* Column-block threading of the full-matrix fill: force small blocks so
 * short sequences are split across several threads */
void simd_test_fill_blocks()
{
  size_t min_cells = alignment_block_min_cells;
  size_t min_cols = alignment_block_min_cols;
  int max_threads = alignment_max_threads;
  aligner_t scalar, striped;
  aligner_init(&scalar);
  aligner_init(&striped);
  char seqa[401], seqb[401], poly_a[301], poly_c[201];
  scoring_t scoring;
  const int threads[] = {2, 3, 5};
  size_t i, t;

  alignment_block_min_cells = 1;
  alignment_block_min_cols = 8;

  for(t = 0; t < sizeof(threads)/sizeof(threads[0]); t++)
  {
    alignment_max_threads = threads[t];
    for(i = 0; i < 200; i++)
    {
      scoring_init(&scoring, 1 + (int)hb_rand(5), -(int)hb_rand(16),
                   -(int)hb_rand(11), -(int)hb_rand(6), hb_rand(2), false,
                   false, false, false, false);
      hb_rand_seq(seqa, hb_rand(401), "acgt");
      hb_rand_seq(seqb, hb_rand(401), "acgt");
      simd_check(seqa, seqb, &scoring, (char)hb_rand(2), &scalar, &striped);
    }
    scoring_system_BLOSUM62(&scoring);
    for(i = 0; i < 40; i++)
    {
      hb_rand_seq(seqa, hb_rand(401), "ARNDCQEGHILKMFPSTWYV");
      hb_rand_seq(seqb, hb_rand(401), "ARNDCQEGHILKMFPSTWYV");
      simd_check(seqa, seqb, &scoring, (char)hb_rand(2), &scalar, &striped);
    }
    // Gaps along seq_a carried through whole blocks
    memset(poly_a, 'a', 300); poly_a[300] = '\0';
    memset(poly_c, 'c', 200); poly_c[200] = '\0';
    scoring_init(&scoring, 1, -20, -2, -1, false, false, false, false,
                 false, true);
    simd_check(poly_a, poly_c, &scoring, 0, &scalar, &striped);
    simd_check(poly_c, poly_a, &scoring, 0, &scalar, &striped);
    simd_check(poly_a, "c", &scoring, 0, &scalar, &striped);
  }

  alignment_block_min_cells = min_cells;
  alignment_block_min_cols = min_cols;
  alignment_max_threads = max_threads;
  aligner_destroy(&scalar);
  aligner_destroy(&striped);
}

/* Block settings below one vector of columns (or zero) must not make empty
 * blocks: these used to corrupt the heap or give 37 for a score of 40 */
void simd_test_tiny_blocks()
{
  size_t a_cells = alignment_block_min_cells, a_cols = alignment_block_min_cols;
  int a_threads = alignment_max_threads;
  size_t h_cells = hirschberg_block_min_cells, h_cols = hirschberg_block_min_cols;
  int h_threads = hirschberg_max_threads;
  bool band = hirschberg_band_enabled;
  aligner_t scalar, striped;
  aligner_init(&scalar);
  aligner_init(&striped);
  alignment_t *hb_scalar = alignment_create(256);
  alignment_t *hb_striped = alignment_create(256);
  char forty[41];
  scoring_t scoring;
  const size_t cols[] = {0, 1};
  const int threads[] = {2, 40};
  size_t c, t;

  memset(forty, 'A', 40); forty[40] = '\0';
  scoring_init(&scoring, 1, -1, -3, -1, false, false, false, false,
               false, false);
  alignment_block_min_cells = hirschberg_block_min_cells = 1;
  hirschberg_band_enabled = false;

  for(c = 0; c < 2; c++)
    for(t = 0; t < 2; t++)
    {
      alignment_block_min_cols = hirschberg_block_min_cols = cols[c];
      alignment_max_threads = hirschberg_max_threads = threads[t];
      simd_check("AA", "A", &scoring, 0, &scalar, &striped);
      simd_check(forty, forty, &scoring, 0, &scalar, &striped);
      simd_hb_check(forty, forty, &scoring, hb_scalar, hb_striped);
      ASSERT(hb_striped->score == 40);
    }

  alignment_block_min_cells = a_cells; alignment_block_min_cols = a_cols;
  alignment_max_threads = a_threads;
  hirschberg_block_min_cells = h_cells; hirschberg_block_min_cols = h_cols;
  hirschberg_max_threads = h_threads;
  hirschberg_band_enabled = band;
  aligner_destroy(&scalar);
  aligner_destroy(&striped);
  alignment_free(hb_scalar);
  alignment_free(hb_striped);
}

static bool alignments_equal(const alignment_t *x, const alignment_t *y)
{
  return x->score == y->score && x->length == y->length &&
         strcmp(x->result_a, y->result_a) == 0 &&
         strcmp(x->result_b, y->result_b) == 0 &&
         x->pos_a == y->pos_a && x->pos_b == y->pos_b &&
         x->len_a == y->len_a && x->len_b == y->len_b;
}

/* Compact mode (traceback bytes instead of score matrices) gives the same
 * Needleman-Wunsch alignment and the same Smith-Waterman hits, in order:
 * first from the bytes, the rest after the matrices are filled */
void simd_test_compact()
{
  size_t min_cells = alignment_block_min_cells;
  size_t min_cols = alignment_block_min_cols;
  int max_threads = alignment_max_threads;
  nw_aligner_t *nw_full = needleman_wunsch_new(), *nw_compact = needleman_wunsch_new();
  sw_aligner_t *sw_full = smith_waterman_new(), *sw_compact = smith_waterman_new();
  alignment_t *r_full = alignment_create(256), *r_compact = alignment_create(256);
  char seqa[401], seqb[401];
  scoring_t scoring;
  const int threads[] = {1, 2, 5};
  const char *alphabets[] = {"acgt", "ac", "a", "acgtn", "ARNDCQEGHILKMFPSTWYV"};
  size_t i, t, traced_nw = 0, traced_sw = 0, hits = 0;

  nw_compact->compact = true;
  smith_waterman_set_compact(sw_compact, true);
  alignment_block_min_cells = 1;
  alignment_block_min_cols = 8;

  for(t = 0; t < sizeof(threads)/sizeof(threads[0]); t++)
  {
    alignment_max_threads = threads[t];
    for(i = 0; i < 400; i++)
    {
      const char *alpha = alphabets[hb_rand(5)];
      if(alpha[0] == 'A') scoring_system_BLOSUM62(&scoring);
      else {
        // Free start gaps (one in 8) are not traced: compared all the same
        scoring_init(&scoring, 1 + (int)hb_rand(5), -(int)hb_rand(16),
                     (int)hb_rand(13) - 10, (int)hb_rand(8) - 6,
                     hb_rand(8) == 0, false, false, false, false, false);
        if(hb_rand(2)) scoring_add_wildcard(&scoring, 'n', 0);
      }
      unsigned int max_len = i % 4 == 0 ? 12 : 401;
      hb_rand_seq(seqa, hb_rand(max_len), alpha);
      hb_rand_seq(seqb, hb_rand(max_len), alpha);

      needleman_wunsch_align(seqa, seqb, &scoring, nw_full, r_full);
      needleman_wunsch_align(seqa, seqb, &scoring, nw_compact, r_compact);
      ASSERT(!nw_full->traced);
      ASSERT(alignments_equal(r_full, r_compact));
      traced_nw += nw_compact->traced;

      score_t min_score = 1 + (score_t)hb_rand(4) * (score_t)hb_rand(10);
      smith_waterman_set_min_score(sw_full, min_score);
      smith_waterman_set_min_score(sw_compact, min_score);
      smith_waterman_align(seqa, seqb, &scoring, sw_full);
      smith_waterman_align(seqa, seqb, &scoring, sw_compact);
      traced_sw += smith_waterman_get_aligner(sw_compact)->traced;
      // Up to 6 hits, or all of them
      size_t h, max_hits = hb_rand(2) ? 6 : SIZE_MAX;
      for(h = 0; h < max_hits; h++) {
        int got_full = smith_waterman_fetch(sw_full, r_full);
        int got_compact = smith_waterman_fetch(sw_compact, r_compact);
        ASSERT(got_full == got_compact);
        if(!got_full || !got_compact) break;
        ASSERT(alignments_equal(r_full, r_compact));
        hits++;
      }
    }
  }
  // Traceback bytes come from the striped fill only
#ifdef ALN_STRIPED
  ASSERT(traced_nw > 600 && traced_sw > 600);
#else
  ASSERT(traced_nw == 0 && traced_sw == 0);
#endif
  ASSERT(hits > 1000);

  // A gap penalty so large that the first column falls below the clamp:
  // working out the bytes would overflow, so the matrices are used
  scoring_init(&scoring, 1, -1, 0, -1073741824, false, false, false, false,
               false, false);
  needleman_wunsch_align("AA", "A", &scoring, nw_full, r_full);
  needleman_wunsch_align("AA", "A", &scoring, nw_compact, r_compact);
  ASSERT(!nw_compact->traced);
  ASSERT(alignments_equal(r_full, r_compact));
  ASSERT(r_compact->score == -1073741823);

  alignment_block_min_cells = min_cells;
  alignment_block_min_cols = min_cols;
  alignment_max_threads = max_threads;
  needleman_wunsch_free(nw_full);
  needleman_wunsch_free(nw_compact);
  smith_waterman_free(sw_full);
  smith_waterman_free(sw_compact);
  alignment_free(r_full);
  alignment_free(r_compact);
}

/* Large compact Smith-Waterman alignments (forced here for small ones) find
 * the best hit from scans and a fill of just its rectangle: the same hits
 * as the full search, first and later ones, local hits in random flanks */
void simd_test_compact_rect()
{
  size_t min_cells = alignment_block_min_cells;
  size_t min_cols = alignment_block_min_cols;
  int max_threads = alignment_max_threads;
  size_t rect_cells = smith_waterman_rect_min_cells;
  sw_aligner_t *full = smith_waterman_new(), *rect = smith_waterman_new();
  alignment_t *r_full = alignment_create(256), *r_rect = alignment_create(256);
  char core[701], hit[901], seqa[2001], seqb[2001];
  scoring_t scoring;
  const unsigned int rates[] = {0, 10, 60, 200};
  size_t i, h, rects = 0, hits = 0;

  smith_waterman_set_compact(rect, true);
  smith_waterman_rect_min_cells = 0;
  alignment_block_min_cells = 1;
  alignment_block_min_cols = 8;

  for(i = 0; i < 600; i++)
  {
    alignment_max_threads = 1 + (int)(i % 3) * 2;
    const char *alpha = i % 2 ? "acgt" : "ac";
    if(i % 7 == 6) {
      scoring_system_BLOSUM62(&scoring);
      alpha = "ARNDCQEGHILKMFPSTWYV";
    }
    else if(i % 5 == 4) {
      // Gaps scoring 0 or above 0 (no rectangle then)
      int extend = (int)hb_rand(6) - 3;
      scoring_init(&scoring, 1 + (int)hb_rand(4), -(int)hb_rand(8),
                   (int)hb_rand(9) - 4 - extend, extend, false, false,
                   false, false, false, false);
    }
    else scoring_init(&scoring, 1 + (int)hb_rand(4), -(int)hb_rand(8),
                      -(int)hb_rand(8), -(int)hb_rand(4), false, false,
                      false, false, false, false);
    // A related piece somewhere in each of two longer sequences
    hb_rand_seq(core, 1 + hb_rand(700), alpha);
    hb_mutate(hit, core, rates[i % 4], alpha);
    size_t fa = hb_rand(600), fb = hb_rand(600), ta = hb_rand(600);
    hb_rand_seq(seqa, fa, alpha);
    strcat(seqa, core);
    hb_rand_seq(seqa + strlen(seqa), ta, alpha);
    hb_rand_seq(seqb, fb, alpha);
    strcat(seqb, hit);
    hb_rand_seq(seqb + strlen(seqb), hb_rand(600), alpha);

    score_t min_score = 1 + (score_t)hb_rand(3) * (score_t)hb_rand(20);
    smith_waterman_set_min_score(full, min_score);
    smith_waterman_set_min_score(rect, min_score);
    smith_waterman_align(seqa, seqb, &scoring, full);
    smith_waterman_align(seqa, seqb, &scoring, rect);
    const aligner_t *ra = smith_waterman_get_aligner(rect);
    rects += ra->traced && ra->score_width < strlen(seqa) + 1;
    if(ra->traced && ra->best_score < min_score) {
      // No hit: the aligner still describes the sequences
      ASSERT(ra->score_width == strlen(seqa) + 1);
      ASSERT(ra->score_height == strlen(seqb) + 1);
      ASSERT(ra->seq_a == seqa && ra->seq_b == seqb);
    }
    // The best hit, sometimes later ones too
    size_t max_hits = i % 3 == 0 ? 4 : 1;
    for(h = 0; h < max_hits; h++) {
      int got_full = smith_waterman_fetch(full, r_full);
      int got_rect = smith_waterman_fetch(rect, r_rect);
      ASSERT(got_full == got_rect);
      if(!got_full || !got_rect) break;
      ASSERT(alignments_equal(r_full, r_rect));
      hits++;
    }
  }
  // Rectangles smaller than the matrix were used (they need the striped
  // fill)
#ifdef ALN_STRIPED
  ASSERT(rects > 300);
#else
  ASSERT(rects == 0);
#endif
  ASSERT(hits > 500);

  smith_waterman_rect_min_cells = rect_cells;
  alignment_block_min_cells = min_cells;
  alignment_block_min_cols = min_cols;
  alignment_max_threads = max_threads;
  smith_waterman_free(full);
  smith_waterman_free(rect);
  alignment_free(r_full);
  alignment_free(r_rect);
}

/* Smith-Waterman fills and scans with 16-bit lanes give the same results as
 * 32-bit ones and the score matrices, including when scores grow past the
 * 16-bit range part-way (the fill starts again with 32-bit lanes), on one
 * thread and on column blocks */
void simd_test_narrow()
{
  size_t min_cells = alignment_block_min_cells;
  size_t min_cols = alignment_block_min_cols;
  int max_threads = alignment_max_threads;
  sw_aligner_t *full = smith_waterman_new(), *compact = smith_waterman_new();
  alignment_t *r_full = alignment_create(256), *r_compact = alignment_create(256);
  char seqa[1201], seqb[2001], core[1001];
  scoring_t scoring;
  const unsigned int rates[] = {0, 10, 60, 300};
  size_t i, past = 0, within = 0;
  size_t fills0, retries0, fills1, retries1;

  smith_waterman_set_compact(compact, true);
  alignment_block_min_cells = 1;
  alignment_block_min_cols = 8;
  alignment_narrow_counts(&fills0, &retries0);

  for(i = 0; i < 400; i++)
  {
    alignment_max_threads = 1 + (int)(i % 4);
    // Large scores in some, so the best passes 32767 - (highest score)
    int match = i % 2 ? 1 + (int)hb_rand(400) : 1 + (int)hb_rand(5);
    scoring_init(&scoring, match, -(int)hb_rand((unsigned int)match * 2 + 2),
                 -(int)hb_rand(300), -(int)hb_rand(100), false, false,
                 false, false, false, false);
    hb_rand_seq(core, hb_rand(1000), i % 3 ? "acgt" : "ac");
    hb_rand_seq(seqa, hb_rand(100), "acgt");
    strcat(seqa, core);
    hb_mutate(seqb, core, rates[i % 4], "acgt");
    hb_rand_seq(seqb + strlen(seqb), hb_rand(100), "acgt");

    smith_waterman_set_min_score(full, 1);
    smith_waterman_set_min_score(compact, 1);
    smith_waterman_align(seqa, seqb, &scoring, full);
    smith_waterman_align(seqa, seqb, &scoring, compact);
    int got_full = smith_waterman_fetch(full, r_full);
    int got_compact = smith_waterman_fetch(compact, r_compact);
    ASSERT(got_full == got_compact);
    if(got_full && got_compact) {
      ASSERT(alignments_equal(r_full, r_compact));
      if(r_full->score > INT16_MAX - match) past++;
      else within++;
    }

    // Scans: the best cell, then the highest cells with that score
    score_t best16, best32, t16, t32;
    size_t x16, y16, x32, y32, tx16, ty16, tx32, ty32;
    size_t la = strlen(seqa), lb = strlen(seqb);
    bool ok16 = alignment_sw_scan(seqa, la, seqb, lb, &scoring, 0, NULL,
                                  &best16, &x16, &y16);
    alignment_narrow_enabled = false;
    bool ok32 = alignment_sw_scan(seqa, la, seqb, lb, &scoring, 0, NULL,
                                  &best32, &x32, &y32);
#ifdef ALN_STRIPED
    ASSERT(ok16 == (la > 0 && lb > 0));
#endif
    ASSERT(ok16 == ok32);
    if(ok16 && ok32) {
      ASSERT(best16 == best32 && x16 == x32 && y16 == y32);
      if(best32 > 0) {
        ok32 = alignment_sw_scan(seqa, la, seqb, lb, &scoring, best32, NULL,
                                 &t32, &tx32, &ty32);
        alignment_narrow_enabled = true;
        ok16 = alignment_sw_scan(seqa, la, seqb, lb, &scoring, best32, NULL,
                                 &t16, &tx16, &ty16);
        ASSERT(ok16 == ok32);
        ASSERT(!ok16 || (t16 == t32 && tx16 == tx32 && ty16 == ty32));
      }
    }
    alignment_narrow_enabled = true;
  }
  // Both sides of the 16-bit range were tested
  ASSERT(past > 50);
  ASSERT(within > 50);
  // The 16-bit fills ran, and every best hit past the range was filled
  // again with 32-bit lanes (the scans may stop earlier)
  alignment_narrow_counts(&fills1, &retries1);
#ifdef ALN_STRIPED
  ASSERT(fills1 - fills0 > 600); // a fill and a scan for most pairs
  ASSERT(retries1 - retries0 >= past);
#else
  ASSERT(fills1 == fills0 && retries1 == retries0);
#endif

  alignment_block_min_cells = min_cells;
  alignment_block_min_cols = min_cols;
  alignment_max_threads = max_threads;
  smith_waterman_free(full);
  smith_waterman_free(compact);
  alignment_free(r_full);
  alignment_free(r_compact);
}

/* A gap scoring above 0 can open a local hit, which the rectangle's
 * reverse scan would miss: the compact aligner gives the same hit */
void simd_test_rect_gap_score()
{
  size_t rect_cells = smith_waterman_rect_min_cells;
  sw_aligner_t *full = smith_waterman_new(), *rect = smith_waterman_new();
  alignment_t *r_full = alignment_create(256), *r_rect = alignment_create(256);
  char seqa[400], seqb[400];
  scoring_t scoring;
  size_t pad;

  smith_waterman_set_compact(rect, true);
  smith_waterman_rect_min_cells = 0;
  scoring_init(&scoring, 2, -1000000, 4, -3, false, false,
               false, false, false, false);
  scoring_add_mutation(&scoring, 'a', 'c', -1);
  scoring_add_mutation(&scoring, 'a', 'g', -1);
  scoring_add_mutation(&scoring, 'c', 'g', -1);
  scoring_add_mutation(&scoring, 'x', 'g', 5);

  for(pad = 0; pad < 300; pad += 37)
  {
    strcpy(seqa, "ACCAAGG");
    memset(seqa + 7, 'Y', pad);
    strcpy(seqa + 7 + pad, "X");
    strcpy(seqb, "GACAAC");
    memset(seqb + 6, 'Z', pad + 2);
    seqb[8 + pad] = '\0';
    smith_waterman_set_min_score(full, 1);
    smith_waterman_set_min_score(rect, 1);
    smith_waterman_align(seqa, seqb, &scoring, full);
    smith_waterman_align(seqa, seqb, &scoring, rect);
    ASSERT(smith_waterman_fetch(full, r_full) == 1);
    ASSERT(smith_waterman_fetch(rect, r_rect) == 1);
    ASSERT(alignments_equal(r_full, r_rect));
  }

  smith_waterman_rect_min_cells = rect_cells;
  smith_waterman_free(full);
  smith_waterman_free(rect);
  alignment_free(r_full);
  alignment_free(r_rect);
}

/* Banded compact Needleman-Wunsch fills give the same alignment as the
 * score matrices: narrow first bands (widened when the score found needs
 * it), strips of a few rows, threaded strips, DNA and protein */
void simd_test_compact_band()
{
  size_t min_cells = alignment_block_min_cells;
  size_t min_cols = alignment_block_min_cols;
  int max_threads = alignment_max_threads;
  size_t strip_rows = alignment_band_strip_rows;
  size_t min_width = alignment_band_min_width;
  nw_aligner_t *full = needleman_wunsch_new(), *compact = needleman_wunsch_new();
  alignment_t *r_full = alignment_create(256), *r_compact = alignment_create(256);
  char seqa[1601], seqb[1801];
  scoring_t scoring;
  const size_t strips[] = {1, 5, 32, 0};
  const unsigned int rates[] = {0, 5, 30, 120};
  size_t s, i, banded = 0;

  compact->compact = true;
  alignment_block_min_cells = 1;
  alignment_block_min_cols = 8;

  for(s = 0; s < 4; s++)
  {
    alignment_band_strip_rows = strips[s];
    alignment_max_threads = 1 + (int)s % 3 * 2; // 1, 3, 5, 1
    for(i = 0; i < 150; i++)
    {
      alignment_band_min_width = i % 2 ? 1 : 1 + hb_rand(16);
      const char *alpha = "acgt";
      if(i % 5 == 4) {
        scoring_system_BLOSUM62(&scoring);
        alpha = "ARNDCQEGHILKMFPSTWYV";
      }
      else scoring_init(&scoring, 1 + (int)hb_rand(3), -(int)hb_rand(6),
                        -(int)hb_rand(8), -(int)hb_rand(3), false, false,
                        false, false, false, false);
      hb_rand_seq(seqa, 50 + hb_rand(1551), alpha);
      hb_mutate(seqb, seqa, rates[i % 4], alpha);

      needleman_wunsch_align(seqa, seqb, &scoring, full, r_full);
      needleman_wunsch_align(seqa, seqb, &scoring, compact, r_compact);
      ASSERT(alignments_equal(r_full, r_compact));
      banded += compact->traced && compact->trace_banded;
    }
  }
  // Tight bounds: no gap open, and b is a with one base moved, so the band
  // that holds is one diagonal either side and the best paths run along its
  // edges, where cells must still be exact
  scoring_init(&scoring, 1, -2, 0, -1, false, false, false, false,
               false, false);
  alignment_band_min_width = 1;
  for(i = 0; i < 400; i++)
  {
    alignment_band_strip_rows = 1 + i % 7;
    alignment_max_threads = 1 + (int)(i % 2) * 2;
    hb_rand_seq(seqa, 100 + hb_rand(600), "acgt");
    size_t len = strlen(seqa), from = hb_rand((unsigned int)len);
    size_t to = hb_rand((unsigned int)len), j, o = 0;
    for(j = 0; j < len; j++) {
      if(j == to) seqb[o++] = "acgt"[hb_rand(4)];
      if(j != from) seqb[o++] = seqa[j];
    }
    seqb[o] = '\0';
    needleman_wunsch_align(seqa, seqb, &scoring, full, r_full);
    needleman_wunsch_align(seqa, seqb, &scoring, compact, r_compact);
    ASSERT(alignments_equal(r_full, r_compact));
    banded += compact->traced && compact->trace_banded;
  }

  // Bands held for many (they need the striped fill)
#ifdef ALN_STRIPED
  ASSERT(banded > 400);
#else
  ASSERT(banded == 0);
#endif

  alignment_block_min_cells = min_cells;
  alignment_block_min_cols = min_cols;
  alignment_max_threads = max_threads;
  alignment_band_strip_rows = strip_rows;
  alignment_band_min_width = min_width;
  needleman_wunsch_free(full);
  needleman_wunsch_free(compact);
  alignment_free(r_full);
  alignment_free(r_compact);
}

void test_simd()
{
  SUITE_START("Striped SIMD fill");
  simd_test_compact_rect();
  simd_test_rect_gap_score();
  simd_test_narrow();
  simd_test_compact_band();
  simd_test_compact();
  simd_test_tiny_blocks();
  simd_test_fill_blocks();
  simd_test_hirschberg_band();
  simd_test_hirschberg_blocks();
  simd_test_hirschberg();
  simd_test_matches_scalar();
  simd_test_odd_penalties();
  simd_test_directional_pairs();
  SUITE_END();
}

void test_nw()
{
  SUITE_START("Needleman-Wunsch");

  nw_test_no_gaps_in_longer();
  nw_test_no_gaps_in_longest();
  nw_test_free_gaps_at_ends();
  nw_test_no_mismatches();
  nw_test_no_mismatches_rand();
  nw_test_high_bit_chars();
  nw_test_penalties_changed_after_init();

  SUITE_END();
}

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

void test_sw()
{
  SUITE_START("Smith-Waterman");

  sw_test_no_gaps_smith_waterman();
  sw_test_reuse_aligner();

  SUITE_END();
}

void alignment_test_create_zero()
{
  alignment_t *result = alignment_create(0);
  ASSERT(result->capacity >= 1 && result->result_a[0] == '\0');
  alignment_ensure_capacity(result, 100);
  ASSERT(result->capacity >= 101);
  alignment_free(result);
}

// Write text to a temporary file and load it as a substitution matrix
static void load_matrix_text(const char *text, scoring_t *scoring)
{
  char path[] = "/tmp/seq_align_testXXXXXX";
  int fd = mkstemp(path);
  if(fd < 0) die("Cannot create temp file");
  if(write(fd, text, strlen(text)) != (ssize_t)strlen(text)) die("write failed");
  close(fd);

  gzFile gz = gzopen(path, "r");
  if(gz == NULL) die("Cannot open %s", path);
  align_scoring_load_matrix(gz, path, scoring, false);
  gzclose(gz);
  unlink(path);
}

void scoring_test_load_matrix()
{
  scoring_t scoring;
  int score;
  bool is_match;

  // Whitespace separated, as in scoring/*.txt
  scoring_init(&scoring, 0, 0, -4, -1, false, false, false, false, false, false);
  load_matrix_text("# comment\n   A  C\nA  5 -3\nC -3  7\n", &scoring);
  scoring_lookup(&scoring, 'a', 'c', &score, &is_match);
  ASSERT(score == -3);
  scoring_lookup(&scoring, 'C', 'C', &score, &is_match);
  ASSERT(score == 7);

  // Single character separator
  scoring_init(&scoring, 0, 0, -4, -1, false, false, false, false, false, false);
  load_matrix_text(",A,C\nA,5,-3\nC,-3,7\n", &scoring);
  scoring_lookup(&scoring, 'a', 'a', &score, &is_match);
  ASSERT(score == 5);
  scoring_lookup(&scoring, 'c', 'a', &score, &is_match);
  ASSERT(score == -3);
  scoring_lookup(&scoring, 'c', 'c', &score, &is_match);
  ASSERT(score == 7);
}

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

/* scoring_make_case_sensitive() keeps scores already set, for upper-case
 * characters too, and tells a character from its other case */
void scoring_test_case_sensitive()
{
  scoring_t scoring;
  int s1, s2;
  bool m1, m2;

  scoring_system_BLOSUM62(&scoring);
  scoring_lookup(&scoring, 'w', 'r', &s1, &m1);
  scoring_make_case_sensitive(&scoring);
  scoring_lookup(&scoring, 'W', 'R', &s2, &m2);
  ASSERT(s1 == s2 && m1 == m2);
  scoring_lookup(&scoring, 'w', 'r', &s2, &m2);
  ASSERT(s1 == s2 && m1 == m2);
  scoring_lookup(&scoring, 'W', 'r', &s2, &m2);
  ASSERT(s1 == s2 && m1 == m2);
  scoring_lookup(&scoring, 'w', 'R', &s2, &m2);
  ASSERT(s1 == s2 && m1 == m2);
  scoring_lookup(&scoring, 'W', 'w', &s2, &m2); // match/mismatch
  ASSERT(s2 == scoring.mismatch && !m2);

  scoring_system_default(&scoring);
  scoring_add_wildcard(&scoring, 'n', 0);
  scoring_make_case_sensitive(&scoring);
  scoring_lookup(&scoring, 'A', 'a', &s1, &m1);
  ASSERT(s1 == scoring.mismatch && !m1);
  scoring_lookup(&scoring, 'N', 'C', &s1, &m1);
  ASSERT(s1 == 0 && m1);
}

/* Aligners with options of their own give the same results as the defaults,
 * also on two threads at once with different options (no shared settings) */
#define OPT_PAIRS 40
typedef struct {
  const alignment_options_t *opts;
  char (*a)[401], (*b)[701];
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

void test_options()
{
  static char a[OPT_PAIRS][401], b[OPT_PAIRS][701];
  alignment_options_t many, plain;
  size_t i;
  for(i = 0; i < OPT_PAIRS; i++) {
    hb_rand_seq(a[i], 50 + hb_rand(350), "acgt");
    hb_mutate(b[i], a[i], i % 2 ? 50 : 300, "acgt");
  }
  // Every switch the other way, and blocks on several threads
  alignment_options_default(&many);
  many.max_threads = many.hb_max_threads = 4;
  many.block_min_cells = many.hb_block_min_cells = 1;
  many.block_min_cols = many.hb_block_min_cols = 8;
  many.band = many.hb_band = !many.band;
  many.narrow = !many.narrow;
  many.sw_rect_min_cells = 0;
  alignment_options_default(&plain);

  // The aligner's own options are the ones used: 16-bit lanes off there
  // start no 16-bit fill, while the defaults do
  {
    sw_aligner_t *sw = smith_waterman_new();
    alignment_options_t off;
    scoring_t scoring;
    size_t fills0, fills1, fills2, retries;
    scoring_system_default(&scoring);
    smith_waterman_set_compact(sw, true);
    alignment_options_default(&off);
    off.narrow = false;
    alignment_narrow_counts(&fills0, &retries);
    smith_waterman_set_options(sw, &off);
    smith_waterman_align(a[0], b[0], &scoring, sw);
    alignment_narrow_counts(&fills1, &retries);
    smith_waterman_set_options(sw, NULL);
    smith_waterman_align(a[0], b[0], &scoring, sw);
    alignment_narrow_counts(&fills2, &retries);
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

void test_misc()
{
  SUITE_START("Alignment buffers and scoring files");

  alignment_test_create_zero();
  alignment_test_roundup();
  scoring_test_load_matrix();
  scoring_test_case_sensitive();
  test_global_inputs();
  test_options();
  test_checked();

  SUITE_END();
}

int main(int argc, char **argv)
{
  if(argc != 1)
  {
    printf("  Unused args '%s..'\n", argv[1]);
    printf("Usage: ./seq_align_tests\n");
    exit(EXIT_FAILURE);
  }

  // Initialise random number generator: SEQ_ALIGN_TEST_SEED repeats a run
  const char *seed_env = getenv("SEQ_ALIGN_TEST_SEED");
  unsigned int seed = seed_env != NULL
                      ? (unsigned int)strtoul(seed_env, NULL, 10)
                      : (unsigned int)time(NULL) + (unsigned int)getpid();
  srand(seed);
  printf("  Random seed: %u (set SEQ_ALIGN_TEST_SEED to repeat)\n", seed);

  printf("  Test Nucleics Align C library:\n\n");

  // Test suites go here
  test_nw();
  test_hirschberg();
  test_simd();
  test_sw();
  test_misc();

  printf("\n");
  printf(" %i / %i suites failed\n", suites_failed, suites_run);
  printf(" %i / %i suites empty\n", suites_empty, suites_run);
  printf(" %i / %i tests failed\n", tests_failed, tests_run);

  printf("\n THE END.\n");

  return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
