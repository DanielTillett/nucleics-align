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
#include <assert.h>

#include "alignment.h"
#include "needleman_wunsch.h"
#include "smith_waterman.h"
#include "hirschberg.h"

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
  assert(a && b);
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
    assert(sub_table_build(&tbl, &sc, a, la, b, lb));
    score_t ref = 0;
    assert(hirschberg_score2(a, b, la, lb, &sc, &ref));

    // Tiny budget must decline cleanly
    assert(alignment_wfa_align(&tbl, &sc, a, la, b, lb, 50, aln) ==
           ALN_ERR_UNSUPPORTED);
    // Sufficient budget must succeed via BiWFA (> 2048 cutoff) and rescore identically
    assert(alignment_wfa_align(&tbl, &sc, a, la, b, lb, 5000000, aln) == ALN_OK);
    assert(aln->score == ref);
    assert(rescore_global(aln->result_a, aln->result_b, &sc) == ref);
    assert(same_letters(aln->result_a, a) && same_letters(aln->result_b, b));
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
      assert(a && b);
      rand_dna(a, la);
      mutate_dna(b, a, rates[ri]);
      size_t lb = strlen(b);

      sub_table_t tbl;
      assert(sub_table_build(&tbl, &sc, a, la, b, lb));
      score_t ms = 0, bs = 0, hs = 0;
      assert(alignment_myers_score(&tbl, &sc, a, la, b, lb, &ms));
      if(alignment_myers_score_banded(&tbl, &sc, a, la, b, lb, &bs))
        assert(bs == ms);
      assert(hirschberg_align_checked(a, b, la, lb, &sc, &opts, aln, NULL) == ALN_OK);
      hs = aln->score;
      assert(ms == hs);
      assert(rescore_global(aln->result_a, aln->result_b, &sc) == ms);
      assert(same_letters(aln->result_a, a) && same_letters(aln->result_b, b));
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
  assert(smith_waterman_fetch(sw, r1) == 1);
  assert(r1->score == 6);

  // 2. Codex multi-hit counterexample: ACAA vs ACCCAA (+2/-1/0/-1)
  scoring_t sc2;
  scoring_init(&sc2, 2, -1, 0, -1, false, false, false, false, false, true);
  smith_waterman_align("ACAA", "ACCCAA", &sc2, sw);
  assert(smith_waterman_fetch(sw, r1) == 1);
  assert(r1->score == 6);
  int f2 = smith_waterman_fetch(sw, r2);
  assert(f2 == 0 || (f2 == 1 && r2->score <= r1->score));

  // 3. 3000 bp pair: exact optimal score match vs default Gotoh + 2nd hit
  scoring_t sc3;
  scoring_init(&sc3, 1, -2, -4, -1, false, false, false, false, false, true);
  char *a = malloc(3100), *b = malloc(3500);
  assert(a && b);
  rand_dna(a, 3000);
  mutate_dna(b, a, 50);

  smith_waterman_fast_rect = false;
  smith_waterman_align(a, b, &sc3, sw);
  assert(smith_waterman_fetch(sw, rdef) == 1);

  smith_waterman_fast_rect = true;
  smith_waterman_align(a, b, &sc3, sw);
  assert(smith_waterman_fetch(sw, r1) == 1);
  assert(r1->score == rdef->score);
  f2 = smith_waterman_fetch(sw, r2);
  assert(f2 == 0 || (f2 == 1 && r2->score <= r1->score));

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

int main(void)
{
  printf("Running v1.3.0 targeted regression tests...\n");
  test_biwfa_deep();
  test_myers_boundaries();
  test_sw_fastlocal();
  printf("All v1.3.0 targeted tests passed!\n");
  return 0;
}
