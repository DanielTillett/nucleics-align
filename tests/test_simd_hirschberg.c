/*
 tests/test_simd_hirschberg.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 The SIMD linear-space passes give the scalar passes' alignment
 */

#include "test_common.h"


void simd_test_hirschberg()
{
  alignment_t *scalar = alignment_create(256), *striped = alignment_create(256);
  static char seqa[801], seqb[801];
  scoring_t scoring;
  size_t i;

  for(i = 0; i < 400; i++)
  {
    // Includes gap_open/gap_extend of 0, and positive values (scalar only)
    int gap_open = (int)test_rand(13) - (i % 5 ? 11 : 1);
    int gap_extend = (int)test_rand(6) - (i % 7 ? 5 : 1);
    scoring_init(&scoring, 1 + (int)test_rand(5), -(int)test_rand(16),
                 gap_open, gap_extend, false, false, false, false,
                 false, false);
    if(test_rand(2)) scoring_add_wildcard(&scoring, 'n', 0);

    // Lengths either side of the thresholds for using SIMD
    unsigned int max_a = i % 3 ? 801 : 70, max_b = i % 4 ? 801 : 20;
    test_rand_seq(seqa, test_rand(max_a), "acgtn");
    test_rand_seq(seqb, test_rand(max_b), "ACgtN");
    simd_hb_check(seqa, seqb, &scoring, scalar, striped);
  }

  scoring_system_BLOSUM62(&scoring);
  for(i = 0; i < 100; i++)
  {
    test_rand_seq(seqa, test_rand(801), "ARNDCQEGHILKMFPSTWYV");
    test_rand_seq(seqb, test_rand(801), "ARNDCQEGHILKMFPSTWYV");
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
      scoring_init(&scoring, 1 + (int)test_rand(5), -(int)test_rand(16),
                   -(int)test_rand(11), -(int)test_rand(6), false, false, false,
                   false, false, false);
      test_rand_seq(seqa, test_rand(801), "acgt");
      test_rand_seq(seqb, test_rand(801), "acgt");
      simd_hb_check(seqa, seqb, &scoring, scalar, striped);
    }
    scoring_system_BLOSUM62(&scoring);
    for(i = 0; i < 20; i++)
    {
      test_rand_seq(seqa, test_rand(801), "ARNDCQEGHILKMFPSTWYV");
      test_rand_seq(seqb, test_rand(801), "ARNDCQEGHILKMFPSTWYV");
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
      scoring_init(&scoring, 1 + (int)test_rand(5), -(int)test_rand(16),
                   -(int)test_rand(11), -(int)test_rand(6), false, false, false,
                   false, false, false);
      test_rand_seq(seqa, 100 + test_rand(1101), "acgt");
      test_mutate(seqb, seqa, rates[i % 4], "acgt");
      simd_hb_check(seqa, seqb, &scoring, scalar, striped);
      simd_hb_check(seqb, seqa, &scoring, scalar, striped);
    }
    scoring_system_BLOSUM62(&scoring);
    for(i = 0; i < 8; i++)
    {
      test_rand_seq(seqa, 100 + test_rand(1101), "ARNDCQEGHILKMFPSTWYV");
      test_mutate(seqb, seqa, rates[i % 4], "ARNDCQEGHILKMFPSTWYV");
      simd_hb_check(seqa, seqb, &scoring, scalar, striped);
    }
  }

  // The band was used (it needs the SIMD passes), and unbanded passes still
  // matched the scalar ones
  ASSERT(hirschberg_banded_passes > banded || !hirschberg_simd_available(NULL));
  ASSERT(hirschberg_striped_mismatches == mismatches);

  hirschberg_band_strip_rows = strip_rows;
  hirschberg_band_min_width = min_width;
  hirschberg_block_min_cells = min_cells;
  hirschberg_block_min_cols = min_cols;
  hirschberg_max_threads = max_threads;
  alignment_free(scalar);
  alignment_free(striped);
}

/* Rows reused down the linear-space recursion (hirschberg_reuse_depth):
 * every depth gives exactly the alignment without reuse, with and without
 * SIMD, bands, column blocks, and deletions across the middle row (which
 * shift the halves' rows so saved ones go unused) */
static void hb_reuse_check(const char *a, const char *b,
                           const scoring_t *scoring, size_t depth,
                           alignment_t *ref, alignment_t *aln)
{
  hirschberg_reuse_depth = 0;
  ASSERT(hirschberg_align(a, b, scoring, ref));
  hirschberg_reuse_depth = depth;
  hirschberg_check_striped = true; // reused rows must match the scalar ones
  ASSERT(hirschberg_align(a, b, scoring, aln));
  hirschberg_check_striped = false;
  ASSERT(aln->score == ref->score);
  ASSERT(strcmp(aln->result_a, ref->result_a) == 0);
  ASSERT(strcmp(aln->result_b, ref->result_b) == 0);
}

void simd_test_hirschberg_reuse()
{
  size_t depth = hirschberg_reuse_depth, min_rows = hirschberg_reuse_min_rows;
  size_t min_cells = hirschberg_block_min_cells;
  size_t min_cols = hirschberg_block_min_cols;
  int max_threads = hirschberg_max_threads;
  bool band = hirschberg_band_enabled, striped = alignment_striped_enabled;
  size_t reused = hirschberg_reused_rows;
  size_t mismatches = hirschberg_striped_mismatches;
  alignment_t *ref = alignment_create(256), *aln = alignment_create(256);
  static char seqa[701], seqb[1402];
  scoring_t scoring;
  const size_t depths[] = {1, 2, 3, 8, 50};
  const unsigned int rates[] = {0, 30, 150, 400};
  size_t i;

  for(i = 0; i < 480; i++)
  {
    hirschberg_reuse_min_rows = i % 3 ? 2 : 1 + test_rand(40);
    hirschberg_band_enabled = i % 4 == 0;
    alignment_striped_enabled = i % 5 != 0;
    if(i % 6 == 0) {
      // Passes split into column blocks on several threads
      hirschberg_block_min_cells = 1;
      hirschberg_block_min_cols = 8;
      hirschberg_max_threads = 2 + (int)test_rand(3);
    } else {
      hirschberg_block_min_cells = min_cells;
      hirschberg_block_min_cols = min_cols;
      hirschberg_max_threads = max_threads;
    }
    // Costly gap opens, so deletions often span the middle row
    scoring_init(&scoring, 1 + (int)test_rand(5), -(int)test_rand(16),
                 -(int)test_rand(i % 2 ? 11 : 40), -(int)test_rand(6),
                 false, false, false, false, false, false);
    test_rand_seq(seqa, test_rand(701), "acgt");
    if(i % 8 == 7) test_rand_seq(seqb, test_rand(701), "acgt");
    else test_mutate(seqb, seqa, rates[i % 4], "acgt");
    hb_reuse_check(seqa, seqb, &scoring, depths[i % 5], ref, aln);
    hb_reuse_check(seqb, seqa, &scoring, depths[(i+2) % 5], ref, aln);
  }

  alignment_striped_enabled = true;
  hirschberg_band_enabled = false;
  scoring_system_BLOSUM62(&scoring);
  for(i = 0; i < 40; i++)
  {
    hirschberg_reuse_min_rows = 2;
    test_rand_seq(seqa, test_rand(701), "ARNDCQEGHILKMFPSTWYV");
    test_mutate(seqb, seqa, rates[i % 4], "ARNDCQEGHILKMFPSTWYV");
    hb_reuse_check(seqa, seqb, &scoring, depths[i % 5], ref, aln);
  }

  // A deletion across the middle row of a long run
  memset(seqa, 'a', 600); seqa[600] = 0;
  memset(seqb, 'a', 300); seqb[300] = 0;
  scoring_init(&scoring, 1, -4, -30, -1, false, false, false, false,
               false, false);
  hb_reuse_check(seqa, seqb, &scoring, 2, ref, aln);
  hb_reuse_check(seqb, seqa, &scoring, 2, ref, aln);

  // Rows were reused, and every SIMD pass matched the scalar one
  ASSERT(hirschberg_reused_rows > reused + 1000);
  ASSERT(hirschberg_striped_mismatches == mismatches);

  hirschberg_reuse_depth = depth;
  hirschberg_reuse_min_rows = min_rows;
  hirschberg_block_min_cells = min_cells;
  hirschberg_block_min_cols = min_cols;
  hirschberg_max_threads = max_threads;
  hirschberg_band_enabled = band;
  alignment_striped_enabled = striped;
  alignment_free(ref);
  alignment_free(aln);
}

/* 16-bit lanes in the linear-space passes: the same alignment as 32-bit
 * lanes and rows matching the scalar passes, with blocks cut into many
 * narrow parts, column blocks on several threads, rows reused from the
 * passes above (their tops), penalties near the limit of what 16 bits hold
 * and penalties too big for them (32-bit lanes then) */
static void hb_narrow_check(const char *a, const char *b,
                            const scoring_t *scoring,
                            alignment_t *ref, alignment_t *aln)
{
  alignment_narrow_enabled = false;
  ASSERT(hirschberg_align(a, b, scoring, ref));
  alignment_narrow_enabled = true;
  hirschberg_check_striped = true;
  ASSERT(hirschberg_align(a, b, scoring, aln));
  hirschberg_check_striped = false;
  ASSERT(aln->score == ref->score);
  ASSERT(strcmp(aln->result_a, ref->result_a) == 0);
  ASSERT(strcmp(aln->result_b, ref->result_b) == 0);
  score_t s32 = 0, s16 = 0;
  alignment_narrow_enabled = false;
  ASSERT(hirschberg_score(a, b, scoring, &s32));
  alignment_narrow_enabled = true;
  ASSERT(hirschberg_score(a, b, scoring, &s16));
  ASSERT(s16 == s32 && s16 == ref->score);
}

void simd_test_hirschberg_narrow()
{
  size_t min_cells = hirschberg_block_min_cells;
  size_t min_cols = hirschberg_block_min_cols;
  int max_threads = hirschberg_max_threads;
  bool band = hirschberg_band_enabled, narrow = alignment_narrow_enabled;
  size_t part_cols = hirschberg_part16_max_cols, parts = hirschberg_parts16;
  size_t mismatches = hirschberg_striped_mismatches;
  alignment_t *ref = alignment_create(256), *aln = alignment_create(256);
  static char seqa[1401], seqb[1401];
  scoring_t scoring;
  const unsigned int rates[] = {0, 30, 150, 400};
  size_t i;

  hirschberg_band_enabled = false;
  for(i = 0; i < 300; i++)
  {
    // Parts of 4 vectors (the fewest allowed) up to unlimited
    hirschberg_part16_max_cols = i % 3 == 2 ? 0 : 1 + test_rand(400);
    if(i % 4 == 0) {
      hirschberg_block_min_cells = 1;
      hirschberg_block_min_cols = 8;
      hirschberg_max_threads = 2 + (int)test_rand(3);
    } else {
      hirschberg_block_min_cells = min_cells;
      hirschberg_block_min_cols = min_cols;
      hirschberg_max_threads = max_threads;
    }
    // Mostly small penalties; every tenth near or past the 16-bit limit
    int big = i % 10 == 9 ? (int)test_rand(4000) : 0;
    scoring_init(&scoring, 1 + (int)test_rand(5) + big / 4,
                 -(int)test_rand(16) - big, -(int)test_rand(12) - big,
                 -(int)test_rand(6) - big / 8,
                 false, false, false, false, false, false);
    if(test_rand(2)) scoring_add_wildcard(&scoring, 'n', 0);
    test_rand_seq(seqa, test_rand(1401), "acgtn");
    if(i % 8 == 7) test_rand_seq(seqb, test_rand(1401), "acgt");
    else test_mutate(seqb, seqa, rates[i % 4], "acgt");
    hb_narrow_check(seqa, seqb, &scoring, ref, aln);
  }

  scoring_system_BLOSUM62(&scoring);
  for(i = 0; i < 40; i++)
  {
    hirschberg_part16_max_cols = 1 + test_rand(300);
    test_rand_seq(seqa, test_rand(1401), "ARNDCQEGHILKMFPSTWYV");
    test_mutate(seqb, seqa, rates[i % 4], "ARNDCQEGHILKMFPSTWYV");
    hb_narrow_check(seqa, seqb, &scoring, ref, aln);
  }

  // Long insertions carried across many parts, and a mismatch-heavy pair
  // whose rows fall steeply across each part
  hirschberg_part16_max_cols = 1;
  memset(seqa, 'a', 1400); seqa[1400] = 0;
  memset(seqb, 'c', 60); seqb[60] = 0;
  scoring_init(&scoring, 1, -20, -2, -1, false, false, false, false,
               false, true);
  hb_narrow_check(seqa, seqb, &scoring, ref, aln);
  hb_narrow_check(seqb, seqa, &scoring, ref, aln);
  scoring_init(&scoring, 5, -40, -11, -9, false, false, false, false,
               false, true);
  hb_narrow_check(seqa, seqb, &scoring, ref, aln);
  hb_narrow_check(seqb, seqa, &scoring, ref, aln);

  // Blocks were cut into parts, and every pass matched the scalar one
  ASSERT(hirschberg_parts16 > parts + 100 || !hirschberg_simd_available(NULL));
  ASSERT(hirschberg_striped_mismatches == mismatches);

  hirschberg_part16_max_cols = part_cols;
  hirschberg_block_min_cells = min_cells;
  hirschberg_block_min_cols = min_cols;
  hirschberg_max_threads = max_threads;
  hirschberg_band_enabled = band;
  alignment_narrow_enabled = narrow;
  alignment_free(ref);
  alignment_free(aln);
}

/* Myers' bit-vector score (unit costs, one thread): the same score as the
 * linear-space passes and the full matrices, across word boundaries (63,
 * 64, 65 ... letters), with a wildcard, case-insensitive letters and
 * one-way substitution pairs; and not used for other scorings */
void simd_test_myers()
{
  int max_threads = hirschberg_max_threads;
  bool myers = hirschberg_myers_enabled;
  size_t uses = hirschberg_myers_scores;
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(256);
  static char seqa[700], seqb[1400];
  const unsigned int lens[] = {1, 2, 63, 64, 65, 127, 128, 129, 300, 640};
  scoring_t scoring;
  score_t s_myers, s_passes;
  size_t i;

  hirschberg_max_threads = 1;
  for(i = 0; i < 400; i++)
  {
    scoring_init(&scoring, 0, -1, 0, -1, false, false, false, false, false,
                 false);
    if(i % 3 == 1) scoring_add_wildcard(&scoring, 'n', i % 2 ? 0 : -1);
    if(i % 5 == 2) {
      // One way only: a against c scores 0, c against a stays -1
      scoring_add_mutation(&scoring, 'a', 'c', 0);
    }
    unsigned int la = i % 4 ? lens[test_rand(10)] : test_rand(640);
    unsigned int lb = i % 4 ? lens[test_rand(10)] : test_rand(640);
    test_rand_seq(seqa, la, i % 2 ? "acgtn" : "ACgtN");
    if(i % 3 == 0) test_mutate(seqb, seqa, 100, "acgt");
    else test_rand_seq(seqb, lb, "AcGtn");

    hirschberg_myers_enabled = true;
    ASSERT(hirschberg_score(seqa, seqb, &scoring, &s_myers));
    hirschberg_myers_enabled = false;
    ASSERT(hirschberg_score(seqa, seqb, &scoring, &s_passes));
    needleman_wunsch_align(seqa, seqb, &scoring, nw, aln);
    ASSERT(s_myers == s_passes && s_myers == aln->score);
  }
  ASSERT(hirschberg_myers_scores > uses + 200);

  // Not unit cost: the passes, not Myers
  hirschberg_myers_enabled = true;
  uses = hirschberg_myers_scores;
  scoring_init(&scoring, 1, -1, 0, -1, false, false, false, false, false,
               false);
  ASSERT(hirschberg_score("acgtacgt", "acgaacgt", &scoring, &s_myers));
  scoring_init(&scoring, 0, -1, -2, -1, false, false, false, false, false,
               false);
  ASSERT(hirschberg_score("acgtacgt", "acgaacgt", &scoring, &s_myers));
  ASSERT(hirschberg_myers_scores == uses);

  hirschberg_max_threads = max_threads;
  hirschberg_myers_enabled = myers;
  alignment_free(aln);
  needleman_wunsch_free(nw);
}

void test_simd_hirschberg(void)
{
  SUITE_START("SIMD linear-space passes");
  RUN(simd_test_hirschberg_band());
  RUN(simd_test_hirschberg_reuse());
  RUN(simd_test_hirschberg_narrow());
  RUN(simd_test_myers());
  RUN(simd_test_hirschberg_blocks());
  RUN(simd_test_hirschberg());
  SUITE_END();
}
