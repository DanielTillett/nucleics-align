/*
 tests/test_v130.c
 Targeted v1.3.0 regression tests:
  1. BiWFA forced > 2048 bp (odd & even GCD, gap-state joins, budget cap)
  2. Banded Myers & Myers-Hirschberg across 63/64/65/512/2000 word boundaries
  3. Smith-Waterman --fastlocal (ATATT/TAATT, ACAA/ACCCAA multi-hit marking, 3kb)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alignment.h"
#include "alignment_simd.h" // ALN_STRIPED
#include "needleman_wunsch.h"
#include "smith_waterman.h"
#include "hirschberg.h"

// Checks on how a fill was made that hold only where the striped fill
// exists (the scalar-only build, -DALN_NO_SIMD, still checks the results)
#ifdef ALN_STRIPED
#  define CHECK_STRIPED(x) CHECK(x)
#else
#  define CHECK_STRIPED(x) ((void)0)
#endif

// Always evaluated (assert() would vanish, calls and all, under NDEBUG)
#define CHECK(cond) do { \
  if(!(cond)) { \
    fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
    exit(EXIT_FAILURE); \
  } \
} while(0)

static unsigned int rng_state = 123456789u;
static unsigned int next_rand(void)
{
  rng_state = rng_state * 1664525u + 1013904223u;
  return rng_state >> 8;
}

static void rand_dna(char *buf, size_t len)
{
  static const char bases[] = "ACGT";
  for(size_t i = 0; i < len; i++) buf[i] = bases[next_rand() & 3];
  buf[len] = '\0';
}

static void mutate_dna(char *dst, const char *src, unsigned int permille)
{
  static const char bases[] = "ACGT";
  size_t o = 0;
  for(size_t i = 0; src[i] != '\0'; ) {
    if(next_rand() % 1000 < permille) {
      unsigned int r = next_rand() % 10;
      if(r < 6) {
        char c;
        do { c = bases[next_rand() & 3]; } while(c == src[i]);
        dst[o++] = c;
        i++;
      } else if(r < 8) {
        i++; // deletion
      } else {
        dst[o++] = bases[next_rand() & 3]; // insertion
      }
    } else {
      dst[o++] = src[i++];
    }
  }
  dst[o] = '\0';
}

static bool same_letters(const char *aln, const char *seq)
{
  for(; *aln != '\0'; aln++) {
    if(*aln == '-') continue;
    if(*aln != *seq) return false;
    seq++;
  }
  return *seq == '\0';
}

static score_t rescore_global(const char *a, const char *b, const scoring_t *sc)
{
  long long tot = 0;
  int st = 0;
  for(size_t i = 0; a[i] != '\0'; i++) {
    if(a[i] == '-') {
      tot += (st == 1 ? sc->gap_extend : sc->gap_open + sc->gap_extend);
      st = 1;
    } else if(b[i] == '-') {
      tot += (st == 2 ? sc->gap_extend : sc->gap_open + sc->gap_extend);
      st = 2;
    } else {
      int s = 0;
      bool m = false;
      scoring_lookup(sc, a[i], b[i], &s, &m);
      tot += s;
      st = 0;
    }
  }
  return (score_t)tot;
}

static void test_biwfa_deep(void)
{
  const int params[][4] = {
    {1, -2, -4, -1}, // odd gcd = 1
    {2, -2, -4, -2}, // even gcd = 2
    {3, -5, -7, -2}, // high coprime penalties
  };
  char *a = malloc(4096), *b = malloc(5000);
  CHECK(a && b);
  alignment_t *aln = alignment_create(8192);

  for(size_t pi = 0; pi < sizeof(params)/sizeof(params[0]); pi++) {
    scoring_t sc;
    scoring_init(&sc, params[pi][0], params[pi][1], params[pi][2], params[pi][3],
                 false, false, false, false, false, true);
    rand_dna(a, 3000);
    mutate_dna(b, a, 20);
    if(b[0] == '\0') strcpy(b, "A");
    size_t la = strlen(a), lb = strlen(b);

    sub_table_t tbl;
    CHECK(sub_table_build(&tbl, &sc, a, la, b, lb));
    score_t ref = 0;
    CHECK(hirschberg_score2(a, b, la, lb, &sc, &ref));

    // Tiny budget must decline cleanly
    CHECK(alignment_wfa_align(&tbl, &sc, a, la, b, lb, 50, aln) ==
           ALN_ERR_UNSUPPORTED);
    // Sufficient budget must succeed via BiWFA (> 2048 cutoff) and rescore identically
    CHECK(alignment_wfa_align(&tbl, &sc, a, la, b, lb, 5000000, aln) == ALN_OK);
    CHECK(aln->score == ref);
    CHECK(rescore_global(aln->result_a, aln->result_b, &sc) == ref);
    CHECK(same_letters(aln->result_a, a) && same_letters(aln->result_b, b));
    sub_table_free(&tbl);
  }
  alignment_free(aln);
  free(a);
  free(b);
  printf("  [PASS] test_biwfa_deep\n");
}

static void test_myers_boundaries(void)
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
      CHECK(a && b);
      rand_dna(a, la);
      mutate_dna(b, a, rates[ri]);
      size_t lb = strlen(b);

      sub_table_t tbl;
      CHECK(sub_table_build(&tbl, &sc, a, la, b, lb));
      score_t ms = 0, bs = 0, hs = 0;
      CHECK(alignment_myers_score(&tbl, &sc, a, la, b, lb, &ms));
      if(alignment_myers_score_banded(&tbl, &sc, a, la, b, lb, &bs))
        CHECK(bs == ms);
      CHECK(hirschberg_align_checked(a, b, la, lb, &sc, &opts, aln, NULL) == ALN_OK);
      hs = aln->score;
      CHECK(ms == hs);
      CHECK(rescore_global(aln->result_a, aln->result_b, &sc) == ms);
      CHECK(same_letters(aln->result_a, a) && same_letters(aln->result_b, b));
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
      CHECK(a && b);
      memset(a, 'A', ma); a[ma] = '\0';
      memset(b, 'C', pad); memset(b + pad, 'A', ma); b[nb] = '\0';
      sub_table_t tbl;
      CHECK(sub_table_build(&tbl, &sc, a, ma, b, nb));
      score_t ms = 0, bs = 0;
      CHECK(alignment_myers_score(&tbl, &sc, a, ma, b, nb, &ms));
      CHECK(ms == -(score_t)pad);
      CHECK(alignment_myers_score_banded(&tbl, &sc, a, ma, b, nb, &bs));
      CHECK(bs == ms);
      sub_table_free(&tbl);
      free(a);
      free(b);
    }
  }
  alignment_free(aln);
  printf("  [PASS] test_myers_boundaries\n");
}

static void test_sw_fastlocal(void)
{
  sw_aligner_t *sw = smith_waterman_new();
  smith_waterman_set_compact(sw, true);
  smith_waterman_set_min_score(sw, 1);
  alignment_t *r1 = alignment_create(64), *r2 = alignment_create(64), *rdef = alignment_create(64);

  // Force fast_rect even on tiny pairs for exactness regression tests
  size_t saved_min = smith_waterman_rect_min_cells;
  smith_waterman_rect_min_cells = 1;

  // 1. Codex counterexample: ATATT vs TAATT (+2/-1/-1/-1) -> local score 6
  scoring_t sc1;
  scoring_init(&sc1, 2, -1, -1, -1, false, false, false, false, false, true);
  smith_waterman_fast_rect = true;
  smith_waterman_align("ATATT", "TAATT", &sc1, sw);
  CHECK(smith_waterman_fetch(sw, r1) == 1);
  CHECK(r1->score == 6);

  // 2. Codex multi-hit counterexample: ACAA vs ACCCAA (+2/-1/0/-1)
  scoring_t sc2;
  scoring_init(&sc2, 2, -1, 0, -1, false, false, false, false, false, true);
  smith_waterman_align("ACAA", "ACCCAA", &sc2, sw);
  CHECK(smith_waterman_fetch(sw, r1) == 1);
  CHECK(r1->score == 6);
  int f2 = smith_waterman_fetch(sw, r2);
  CHECK(f2 == 0 || (f2 == 1 && r2->score <= r1->score));

  // 3. 3000 bp pair: exact optimal score match vs default Gotoh + 2nd hit
  scoring_t sc3;
  scoring_init(&sc3, 1, -2, -4, -1, false, false, false, false, false, true);
  char *a = malloc(3100), *b = malloc(3500);
  CHECK(a && b);
  rand_dna(a, 3000);
  mutate_dna(b, a, 50);

  smith_waterman_fast_rect = false;
  smith_waterman_align(a, b, &sc3, sw);
  CHECK(smith_waterman_fetch(sw, rdef) == 1);

  smith_waterman_fast_rect = true;
  smith_waterman_align(a, b, &sc3, sw);
  CHECK(smith_waterman_fetch(sw, r1) == 1);
  CHECK(r1->score == rdef->score);
  f2 = smith_waterman_fetch(sw, r2);
  CHECK(f2 == 0 || (f2 == 1 && r2->score <= r1->score));

  smith_waterman_fast_rect = false;
  smith_waterman_rect_min_cells = saved_min;
  alignment_free(r1);
  alignment_free(r2);
  alignment_free(rdef);
  smith_waterman_free(sw);
  free(a);
  free(b);
  printf("  [PASS] test_sw_fastlocal\n");
}

static void test_strip_col0_and_striped_trace(void)
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
  CHECK(needleman_wunsch_align_checked(ta, tb, 35, 35, &sc, nw1, r1, NULL) == ALN_OK);
  CHECK(needleman_wunsch_align_checked(ta, tb, 35, 35, &sc, nw2, r2, NULL) == ALN_OK);
  CHECK_STRIPED(nw2->trace_banded); // bands need the striped fill
  CHECK(r1->score == r2->score);
  CHECK(strcmp(r1->result_a, r2->result_a) == 0);
  CHECK(strcmp(r1->result_b, r2->result_b) == 0);

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
    CHECK(a && b);
    rand_dna(a, la);
    rand_dna(b, lb);
    // DNA diff8
    nw1->options = &o_scalar;
    CHECK(needleman_wunsch_align_checked(a, b, la, lb, &sc_dna, nw1, r1, NULL) == ALN_OK);
    nw2->options = &o_diff1;
    CHECK(needleman_wunsch_align_checked(a, b, la, lb, &sc_dna, nw2, r2, NULL) == ALN_OK);
    CHECK_STRIPED(nw2->trace_striped);
    size_t dna_vl = nw2->trace_vl;
    CHECK(r1->score == r2->score && strcmp(r1->result_a, r2->result_a) == 0 && strcmp(r1->result_b, r2->result_b) == 0);
    nw2->options = &o_diff4;
    CHECK(needleman_wunsch_align_checked(a, b, la, lb, &sc_dna, nw2, r2, NULL) == ALN_OK);
    CHECK_STRIPED(nw2->trace_striped && nw2->trace_vl == dna_vl);
    CHECK(r1->score == r2->score && strcmp(r1->result_a, r2->result_a) == 0 && strcmp(r1->result_b, r2->result_b) == 0);
    // Protein BLOSUM62 diff16 (ensure 'W' is present so score range forces 16-bit lanes)
    for(size_t i = 0; i < la; i++) a[i] = aa20[next_rand() % 20];
    for(size_t i = 0; i < lb; i++) b[i] = aa20[next_rand() % 20];
    a[0] = 'W'; b[0] = 'W';
    nw1->options = &o_scalar;
    CHECK(needleman_wunsch_align_checked(a, b, la, lb, &sc_b62, nw1, r1, NULL) == ALN_OK);
    nw2->options = &o_diff1;
    CHECK(needleman_wunsch_align_checked(a, b, la, lb, &sc_b62, nw2, r2, NULL) == ALN_OK);
    CHECK_STRIPED(nw2->trace_striped && nw2->trace_vl * 2 == dna_vl);
    CHECK(r1->score == r2->score && strcmp(r1->result_a, r2->result_a) == 0 && strcmp(r1->result_b, r2->result_b) == 0);
    nw2->options = &o_diff4;
    CHECK(needleman_wunsch_align_checked(a, b, la, lb, &sc_b62, nw2, r2, NULL) == ALN_OK);
    CHECK_STRIPED(nw2->trace_striped && nw2->trace_vl * 2 == dna_vl);
    CHECK(r1->score == r2->score && strcmp(r1->result_a, r2->result_a) == 0 && strcmp(r1->result_b, r2->result_b) == 0);
    free(a); free(b);
  }
  alignment_free(r1); alignment_free(r2);
  needleman_wunsch_free(nw1); needleman_wunsch_free(nw2);
  printf("  [PASS] test_strip_col0_and_striped_trace\n");
}

int main(void)
{
  printf("Running v1.3.0 targeted regression tests...\n");
  test_biwfa_deep();
  test_myers_boundaries();
  test_sw_fastlocal();
  test_strip_col0_and_striped_trace();
  printf("All v1.3.0 targeted tests passed!\n");
  return 0;
}
