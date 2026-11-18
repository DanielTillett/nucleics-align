/*
 tests/tests.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
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


#include "needleman_wunsch.h"
#include "smith_waterman.h"
#include "alignment_scoring_load.h"

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

void test_misc()
{
  SUITE_START("Alignment buffers and scoring files");

  alignment_test_create_zero();
  scoring_test_load_matrix();

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

  // Initialise random number generator
  srand((unsigned int)time(NULL) + getpid());

  printf("  Test seq-align C library:\n\n");

  // Test suites go here
  test_nw();
  test_sw();
  test_misc();

  printf("\n");
  printf(" %i / %i suites failed\n", suites_failed, suites_run);
  printf(" %i / %i suites empty\n", suites_empty, suites_run);
  printf(" %i / %i tests failed\n", tests_failed, tests_run);

  printf("\n THE END.\n");

  return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
