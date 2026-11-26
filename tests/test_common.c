/*
 tests/test_common.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

#include "test_common.h"

//
// Suites and assertions
//

int tests_run = 0, tests_failed = 0;
int suites_run = 0, suites_failed = 0, suites_empty = 0;

static const char *suite_name = NULL;
static const char *test_name = "(no test)";
static bool suite_pass = true;
static int tests_in_suite = 0;
static unsigned int run_seed = 1;

void test_counted(void)
{
  tests_run++;
  tests_in_suite++;
}

void test_failed(const char *file, int line, const char *expr)
{
  fprintf(stderr, "failed assert [%s:%i] %s (in %s, seed %u)\n", file, line,
          expr, test_name, run_seed);
  suite_pass = false;
  tests_failed++;
}

void SUITE_START(const char *name)
{
  suite_pass = true;
  suite_name = name;
  suites_run++;
  tests_in_suite = 0;
}

void SUITE_END(void)
{
  printf("Testing %s ", suite_name);
  size_t suite_i;
  for(suite_i = strlen(suite_name); suite_i < 80-8-5; suite_i++) printf(".");
  printf("%s\n", suite_pass ? " pass" : " fail");
  fflush(stdout);
  if(!suite_pass) suites_failed++;
  if(!tests_in_suite) suites_empty++;
}

//
// The tuning globals: saved before each test, restored after it
//

#define TEST_GLOBALS(X) \
  X(bool, alignment_striped_enabled) \
  X(bool, alignment_band_enabled) \
  X(bool, alignment_seeds_enabled) \
  X(bool, alignment_narrow_enabled) \
  X(bool, alignment_byte_enabled) \
  X(bool, alignment_diff_enabled) \
  X(bool, alignment_linear_fallback) \
  X(bool, alignment_wfa_alignment) \
  X(bool, hirschberg_band_enabled) \
  X(bool, hirschberg_check_striped) \
  X(bool, hirschberg_myers_enabled) \
  X(bool, smith_waterman_fast_rect) \
  X(int, alignment_max_threads) \
  X(int, hirschberg_max_threads) \
  X(size_t, alignment_block_min_cells) \
  X(size_t, alignment_block_min_cols) \
  X(size_t, alignment_band_block_min_cells) \
  X(size_t, alignment_band_block_min_cols) \
  X(size_t, alignment_band_strip_rows) \
  X(size_t, alignment_band_min_width) \
  X(size_t, alignment_seed_length) \
  X(size_t, alignment_band_keep_max_cells) \
  X(size_t, alignment_band_wfa_budget) \
  X(size_t, alignment_linear_above) \
  X(size_t, alignment_wfa_budget) \
  X(size_t, hirschberg_block_min_cells) \
  X(size_t, hirschberg_block_min_cols) \
  X(size_t, hirschberg_band_strip_rows) \
  X(size_t, hirschberg_band_min_width) \
  X(size_t, hirschberg_reuse_depth) \
  X(size_t, hirschberg_reuse_min_rows) \
  X(size_t, hirschberg_part16_max_cols) \
  X(size_t, smith_waterman_rect_min_cells) \
  X(size_t, smith_waterman_strip_min_cells) \
  X(size_t, smith_waterman_strip_rows)

typedef struct {
#define X(type, name) type name;
  TEST_GLOBALS(X)
#undef X
} test_globals_t;

static test_globals_t saved_globals;

static void globals_save(test_globals_t *g)
{
#define X(type, name) g->name = name;
  TEST_GLOBALS(X)
#undef X
}

// Restores every global, naming those the test left changed
static void globals_restore(const test_globals_t *g, const char *test)
{
#define X(type, name) \
  if(name != g->name) { \
    fprintf(stderr, "note: %s left %s changed; restored\n", test, #name); \
    name = g->name; \
  }
  TEST_GLOBALS(X)
#undef X
}

//
// Random data
//

static unsigned long long rand_state = 1;

void test_set_seed(unsigned int seed)
{
  run_seed = seed;
}

unsigned int test_seed(void)
{
  return run_seed;
}

// FNV-1a, so each test's stream depends on its name and the run's seed
static unsigned long long name_hash(const char *name)
{
  unsigned long long h = 14695981039346656037ULL;
  for(; *name != '\0'; name++) {
    h ^= (unsigned char)*name;
    h *= 1099511628211ULL;
  }
  return h;
}

void test_begin(const char *name)
{
  test_name = name;
  rand_state = name_hash(name) ^ ((unsigned long long)run_seed << 32);
  globals_save(&saved_globals);
}

void test_end(const char *name)
{
  globals_restore(&saved_globals, name);
  // An allocation set to fail and not reached would fail in the next test
  if(alignment_alloc_fail_pending() != 0) {
    fprintf(stderr, "note: %s left an allocation set to fail; cleared\n", name);
    alignment_alloc_fail(0);
  }
  test_name = "(no test)";
}

unsigned int test_rand(unsigned int n)
{
  rand_state = rand_state * 6364136223846793005ULL + 1442695040888963407ULL;
  return (unsigned int)((rand_state >> 33) % n);
}

void test_rand_seq(char *seq, size_t len, const char *alphabet)
{
  size_t i, alen = strlen(alphabet);
  for(i = 0; i < len; i++) seq[i] = alphabet[test_rand((unsigned int)alen)];
  seq[len] = '\0';
}

void test_mutate(char *out, const char *in, unsigned int per_mille,
                 const char *alphabet)
{
  size_t i, o = 0, len = strlen(in), nalpha = strlen(alphabet);
  for(i = 0; i < len; i++)
  {
    unsigned int r = test_rand(1000);
    if(r == 999) i += 20;
    else if(r < per_mille) out[o++] = alphabet[test_rand((unsigned int)nalpha)];
    else if(r < per_mille + per_mille/4 + 1) continue;
    else if(r < per_mille + per_mille/2 + 2) {
      out[o++] = in[i];
      out[o++] = alphabet[test_rand((unsigned int)nalpha)];
    }
    else out[o++] = in[i];
  }
  out[o] = '\0';
}

//
// Comparisons
//

bool test_ungapped_equals(const char *aln, const char *seq)
{
  for(; *aln != '\0'; aln++) {
    if(*aln == '-') continue;
    if(*aln != *seq) return false;
    seq++;
  }
  return *seq == '\0';
}

bool alignments_equal(const alignment_t *x, const alignment_t *y)
{
  return x->score == y->score && x->length == y->length &&
         strcmp(x->result_a, y->result_a) == 0 &&
         strcmp(x->result_b, y->result_b) == 0 &&
         x->pos_a == y->pos_a && x->pos_b == y->pos_b &&
         x->len_a == y->len_a && x->len_b == y->len_b;
}

void simd_check(const char *a, const char *b, const scoring_t *scoring,
                char is_sw, aligner_t *scalar, aligner_t *striped)
{
  const bool was = alignment_striped_enabled;
  alignment_striped_enabled = false;
  aligner_align(scalar, a, b, strlen(a), strlen(b), scoring, is_sw);
  alignment_striped_enabled = true;
  aligner_align(striped, a, b, strlen(a), strlen(b), scoring, is_sw);
  alignment_striped_enabled = was;

  size_t cells = scalar->score_width * scalar->score_height;
  ASSERT(memcmp(scalar->match_scores, striped->match_scores,
                cells * sizeof(score_t)) == 0);
  ASSERT(memcmp(scalar->gap_a_scores, striped->gap_a_scores,
                cells * sizeof(score_t)) == 0);
  ASSERT(memcmp(scalar->gap_b_scores, striped->gap_b_scores,
                cells * sizeof(score_t)) == 0);
}

void simd_hb_check(const char *a, const char *b, const scoring_t *scoring,
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

  const bool was = alignment_striped_enabled;
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
  alignment_striped_enabled = was;
  ASSERT(score_scalar == scalar->score && score_striped == scalar->score);
  ASSERT(strcmp(scalar->result_a, striped->result_a) == 0);
  ASSERT(strcmp(scalar->result_b, striped->result_b) == 0);
}
