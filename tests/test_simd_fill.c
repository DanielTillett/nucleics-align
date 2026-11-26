/*
 tests/test_simd_fill.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 The striped SIMD fill gives the scalar matrices, traceback bytes and bands
 */

#include "test_common.h"
#include "alignment_seeds.h" // seeds_build(), for the bound test


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
    int match = 1 + (int)test_rand(5), mismatch = -(int)test_rand(16);
    int gap_open = -(int)test_rand(11), gap_extend = -(int)test_rand(6);
    scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
                 test_rand(2), false, false, false, false, false);
    if(test_rand(2)) scoring_add_wildcard(&scoring, 'n', 0);

    unsigned int max_a = i % 4 == 0 ? 20 : 301, max_b = i % 3 == 0 ? 20 : 301;
    test_rand_seq(seqa, test_rand(max_a), "acgtn");
    test_rand_seq(seqb, test_rand(max_b), "ACgtN");
    simd_check(seqa, seqb, &scoring, (char)test_rand(2), &scalar, &striped);
  }

  scoring_system_BLOSUM62(&scoring);
  for(i = 0; i < 300; i++)
  {
    test_rand_seq(seqa, test_rand(301), "ARNDCQEGHILKMFPSTWYV");
    test_rand_seq(seqb, test_rand(301), "ARNDCQEGHILKMFPSTWYV");
    simd_check(seqa, seqb, &scoring, (char)test_rand(2), &scalar, &striped);
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
    int gap_open = (int)test_rand(11) - 5, gap_extend = (int)test_rand(7) - 3;
    scoring_init(&scoring, 1 + (int)test_rand(3), -(int)test_rand(5),
                 gap_open, gap_extend, test_rand(2), false, false, false,
                 false, false);
    test_rand_seq(seqa, test_rand(61), "acgt");
    test_rand_seq(seqb, test_rand(61), "acgt");
    simd_check(seqa, seqb, &scoring, (char)test_rand(2), &scalar, &striped);
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

/* Column-block threading of the full-matrix fill: force small blocks so
 * short sequences are split across several threads */
void simd_test_fill_blocks()
{
  size_t min_cells = alignment_block_min_cells;
  size_t min_cols = alignment_block_min_cols;
  size_t band_cells = alignment_band_block_min_cells;
  size_t band_cols = alignment_band_block_min_cols;
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
  alignment_band_block_min_cells = 1;
  alignment_band_block_min_cols = 8;

  for(t = 0; t < sizeof(threads)/sizeof(threads[0]); t++)
  {
    alignment_max_threads = threads[t];
    for(i = 0; i < 200; i++)
    {
      scoring_init(&scoring, 1 + (int)test_rand(5), -(int)test_rand(16),
                   -(int)test_rand(11), -(int)test_rand(6), test_rand(2), false,
                   false, false, false, false);
      test_rand_seq(seqa, test_rand(401), "acgt");
      test_rand_seq(seqb, test_rand(401), "acgt");
      simd_check(seqa, seqb, &scoring, (char)test_rand(2), &scalar, &striped);
    }
    scoring_system_BLOSUM62(&scoring);
    for(i = 0; i < 40; i++)
    {
      test_rand_seq(seqa, test_rand(401), "ARNDCQEGHILKMFPSTWYV");
      test_rand_seq(seqb, test_rand(401), "ARNDCQEGHILKMFPSTWYV");
      simd_check(seqa, seqb, &scoring, (char)test_rand(2), &scalar, &striped);
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
  alignment_band_block_min_cells = band_cells;
  alignment_band_block_min_cols = band_cols;
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

/* Compact mode (traceback bytes instead of score matrices) gives the same
 * Needleman-Wunsch alignment and the same Smith-Waterman hits, in order:
 * first from the bytes, the rest after the matrices are filled */
void simd_test_compact()
{
  size_t min_cells = alignment_block_min_cells;
  size_t min_cols = alignment_block_min_cols;
  size_t band_cells = alignment_band_block_min_cells;
  size_t band_cols = alignment_band_block_min_cols;
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
  alignment_band_block_min_cells = 1;
  alignment_band_block_min_cols = 8;

  for(t = 0; t < sizeof(threads)/sizeof(threads[0]); t++)
  {
    alignment_max_threads = threads[t];
    for(i = 0; i < 400; i++)
    {
      const char *alpha = alphabets[test_rand(5)];
      if(alpha[0] == 'A') scoring_system_BLOSUM62(&scoring);
      else {
        // Free start gaps (one in 8) are not traced: compared all the same
        scoring_init(&scoring, 1 + (int)test_rand(5), -(int)test_rand(16),
                     (int)test_rand(13) - 10, (int)test_rand(8) - 6,
                     test_rand(8) == 0, false, false, false, false, false);
        if(test_rand(2)) scoring_add_wildcard(&scoring, 'n', 0);
      }
      unsigned int max_len = i % 4 == 0 ? 12 : 401;
      test_rand_seq(seqa, test_rand(max_len), alpha);
      test_rand_seq(seqb, test_rand(max_len), alpha);

      needleman_wunsch_align(seqa, seqb, &scoring, nw_full, r_full);
      needleman_wunsch_align(seqa, seqb, &scoring, nw_compact, r_compact);
      ASSERT(!nw_full->traced);
      ASSERT(alignments_equal(r_full, r_compact));
      traced_nw += nw_compact->traced;

      score_t min_score = 1 + (score_t)test_rand(4) * (score_t)test_rand(10);
      smith_waterman_set_min_score(sw_full, min_score);
      smith_waterman_set_min_score(sw_compact, min_score);
      smith_waterman_align(seqa, seqb, &scoring, sw_full);
      smith_waterman_align(seqa, seqb, &scoring, sw_compact);
      traced_sw += smith_waterman_get_aligner(sw_compact)->traced;
      // Up to 6 hits, or all of them
      size_t h, max_hits = test_rand(2) ? 6 : SIZE_MAX;
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

  // A gap penalty so large that the first row falls below the clamp:
  // refused as overflow rather than silently clamping or overflowing traceback
  scoring_init(&scoring, 1, -1, 0, -1073741824, false, false, false, false,
               false, false);
  ASSERT(needleman_wunsch_align_checked("AA", "A", 2, 1, &scoring, nw_full,
                                        r_full, NULL) == ALN_ERR_OVERFLOW);
  ASSERT(needleman_wunsch_align_checked("AA", "A", 2, 1, &scoring, nw_compact,
                                        r_compact, NULL) == ALN_ERR_OVERFLOW);

  alignment_block_min_cells = min_cells;
  alignment_block_min_cols = min_cols;
  alignment_band_block_min_cells = band_cells;
  alignment_band_block_min_cols = band_cols;
  alignment_max_threads = max_threads;
  needleman_wunsch_free(nw_full);
  needleman_wunsch_free(nw_compact);
  smith_waterman_free(sw_full);
  smith_waterman_free(sw_compact);
  alignment_free(r_full);
  alignment_free(r_compact);
}

/* Needleman-Wunsch traceback bytes from difference-encoded lanes are those
 * of the absolute 32-bit fill, cell for cell, with the same final scores:
 * scorings within the 8-bit range, at its edge, in the 16-bit range and
 * beyond it (absolute lanes then), on one thread and on column blocks */
void simd_test_diff()
{
  alignment_options_t on, off;
  alignment_options_default(&on);
  on.band = false; // the whole fill
  on.block_min_cells = 1;
  on.block_min_cols = 8;
  off = on;
  off.diff = false;
  aligner_t al_on, al_off;
  aligner_init(&al_on);
  aligner_init(&al_off);
  al_on.compact = al_off.compact = true;
  al_on.options = &on;
  al_off.options = &off;
  static char seqa[1201], seqb[2001], core[1001];
  scoring_t scoring;
  const unsigned int rates[] = {0, 20, 100, 400};
  size_t i, x, y, fills = alignment_diff_count(), traced = 0;

  for(i = 0; i < 300; i++)
  {
    on.max_threads = off.max_threads = 1 + (int)(i % 4);
    const char *alpha = i % 2 ? "acgt" : "ac";
    switch(i % 5) {
      case 0: // DNA-like, well within 8 bits
        scoring_init(&scoring, 1 + (int)test_rand(2), -(int)test_rand(4),
                     -(int)test_rand(8), -(int)test_rand(3), false, false, false,
                     false, false, false);
        break;
      case 1: // near the 8-bit limits, either side
        scoring_init(&scoring, 1 + (int)test_rand(12), -(int)test_rand(20),
                     -(int)test_rand(16), -(int)test_rand(6), false, false, false,
                     false, false, false);
        break;
      case 2: // 16 bits
        scoring_init(&scoring, 1 + (int)test_rand(300), -(int)test_rand(600),
                     -(int)test_rand(900), -(int)test_rand(300), false, false,
                     false, false, false, false);
        break;
      case 3: // a gap scoring above 0: absolute lanes
        scoring_init(&scoring, 1 + (int)test_rand(3), -(int)test_rand(4),
                     (int)test_rand(3), (int)test_rand(3) - 1, false, false,
                     false, false, false, false);
        break;
      default:
        scoring_system_BLOSUM62(&scoring);
        alpha = "ARNDCQEGHILKMFPSTWYV";
    }
    test_rand_seq(core, 1 + test_rand(800), alpha);
    test_rand_seq(seqa, test_rand(100), alpha);
    strcat(seqa, core);
    test_mutate(seqb, core, rates[i % 4], alpha);
    test_rand_seq(seqb + strlen(seqb), test_rand(100), alpha);
    const size_t la = strlen(seqa), lb = strlen(seqb);

    aligner_align(&al_on, seqa, seqb, la, lb, &scoring, 0);
    aligner_align(&al_off, seqa, seqb, la, lb, &scoring, 0);
    ASSERT(al_on.traced == al_off.traced);
    if(!al_on.traced || !al_off.traced) continue;
    traced++;
    bool same = true;
    for(y = 1; y <= lb; y++)
      for(x = 1; x <= la; x++)
        same &= al_on.trace[alignment_trace_index(&al_on, x, y)] ==
                al_off.trace[alignment_trace_index(&al_off, x, y)];
    ASSERT(same);
    ASSERT(memcmp(al_on.final_scores, al_off.final_scores,
                  sizeof(al_on.final_scores)) == 0);
  }
#ifdef ALN_STRIPED
  ASSERT(traced > 250); // the traceback bytes come from the striped fill
  ASSERT(alignment_diff_count() - fills > 150); // the 8- and 16-bit fills
#else
  ASSERT(alignment_diff_count() == fills);
#endif
  aligner_destroy(&al_on);
  aligner_destroy(&al_off);
}

/* Banded compact Needleman-Wunsch fills give the same alignment as the
 * score matrices: narrow first bands (widened when the score found needs
 * it), strips of a few rows, threaded strips, DNA and protein. With keep
 * false, every band keeps one strip's bytes at a time, filled again as the
 * traceback reaches it */
void simd_test_compact_band(bool keep)
{
  size_t keep_cells = alignment_band_keep_max_cells;
  size_t wfa_budget = alignment_band_wfa_budget;
  size_t min_cells = alignment_block_min_cells;
  size_t min_cols = alignment_block_min_cols;
  size_t band_cells = alignment_band_block_min_cells;
  size_t band_cols = alignment_band_block_min_cols;
  int max_threads = alignment_max_threads;
  size_t strip_rows = alignment_band_strip_rows;
  size_t min_width = alignment_band_min_width;
  size_t seed_length = alignment_seed_length;
  bool seeds_on = alignment_seeds_enabled;
  nw_aligner_t *full = needleman_wunsch_new(), *compact = needleman_wunsch_new();
  alignment_t *r_full = alignment_create(256), *r_compact = alignment_create(256);
  char seqa[1601], seqb[1801];
  scoring_t scoring;
  const size_t strips[] = {1, 5, 32, 0};
  const unsigned int rates[] = {0, 5, 30, 120};
  size_t s, i, banded = 0, kept = 0, sized = 0, tries, wfa, tries0, wfa0;

  compact->compact = true;
  alignment_block_min_cells = 1;
  alignment_block_min_cols = 8;
  alignment_band_block_min_cells = 1;
  alignment_band_block_min_cols = 8;
  if(!keep) alignment_band_keep_max_cells = 1;

  for(s = 0; s < 4; s++)
  {
    alignment_band_strip_rows = strips[s];
    alignment_max_threads = 1 + (int)s % 3 * 2; // 1, 3, 5, 1
    for(i = 0; i < 150; i++)
    {
      alignment_band_min_width = i % 2 ? 1 : 1 + test_rand(16);
      alignment_band_wfa_budget = i % 3 ? 0 : (size_t)1 << 20;
      // Short seeds, so they narrow strips of these short pairs
      alignment_seed_length = i % 3 ? 2 + test_rand(5) : 0;
      alignment_seeds_enabled = i % 7 != 6;
      const char *alpha = "acgt";
      if(i % 5 == 4) {
        scoring_system_BLOSUM62(&scoring);
        alpha = "ARNDCQEGHILKMFPSTWYV";
      }
      else scoring_init(&scoring, 1 + (int)test_rand(3), -(int)test_rand(6),
                        -(int)test_rand(8), -(int)test_rand(3), false, false,
                        false, false, false, false);
      test_rand_seq(seqa, 50 + test_rand(1551), alpha);
      test_mutate(seqb, seqa, rates[i % 4], alpha);

      needleman_wunsch_align(seqa, seqb, &scoring, full, r_full);
      alignment_band_counts(&tries0, &wfa0);
      needleman_wunsch_align(seqa, seqb, &scoring, compact, r_compact);
      alignment_band_counts(&tries, &wfa);
      // A band sized by the wavefront score passes on its first fill
      if(wfa > wfa0 && compact->trace_banded) {
        ASSERT(tries == tries0 + 1);
        sized++;
      }
      ASSERT(alignments_equal(r_full, r_compact));
      banded += compact->traced && compact->trace_banded;
      kept += compact->traced && compact->trace_strips != NULL;
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
    alignment_band_wfa_budget = i % 3 ? 0 : (size_t)1 << 20;
    alignment_seed_length = i % 3 ? 2 + test_rand(5) : 0;
    alignment_seeds_enabled = i % 7 != 6;
    test_rand_seq(seqa, 100 + test_rand(600), "acgt");
    size_t len = strlen(seqa), from = test_rand((unsigned int)len);
    size_t to = test_rand((unsigned int)len), j, o = 0;
    for(j = 0; j < len; j++) {
      if(j == to) seqb[o++] = "acgt"[test_rand(4)];
      if(j != from) seqb[o++] = seqa[j];
    }
    seqb[o] = '\0';
    needleman_wunsch_align(seqa, seqb, &scoring, full, r_full);
    alignment_band_counts(&tries0, &wfa0);
    needleman_wunsch_align(seqa, seqb, &scoring, compact, r_compact);
    alignment_band_counts(&tries, &wfa);
    if(wfa > wfa0 && compact->trace_banded) {
      ASSERT(tries == tries0 + 1);
      sized++;
    }
    ASSERT(alignments_equal(r_full, r_compact));
    banded += compact->traced && compact->trace_banded;
    kept += compact->traced && compact->trace_strips != NULL;
  }

  // One long gap, in a or in b: the best path crosses many diagonals in a
  // strip, which narrowed strips (see alignment_band_narrow()) must hold
  scoring_init(&scoring, 1, -3, -5, -1, false, false, false, false,
               false, false);
  for(i = 0; i < 300; i++)
  {
    alignment_band_strip_rows = i % 4 ? 1 + i % 9 : 0;
    alignment_band_min_width = i % 2 ? 1 : 64;
    alignment_max_threads = 1 + (int)(i % 2) * 2;
    alignment_band_wfa_budget = i / 3 % 2 ? 0 : (size_t)1 << 20;
    alignment_seed_length = i % 3 ? 2 + test_rand(5) : 0;
    alignment_seeds_enabled = i % 7 != 6;
    test_rand_seq(seqa, 800 + test_rand(600), "acgt");
    size_t len = strlen(seqa), at = test_rand((unsigned int)len);
    size_t gap = 1 + test_rand((unsigned int)len / 10), j, o = 0;
    char *in = i / 2 % 2 ? seqa : seqb, *out = i / 2 % 2 ? seqb : seqa;
    if(in == seqb) memcpy(seqb, seqa, len + 1);
    for(j = 0; j < len; j++) {
      if(j == at) while(o < at + gap) out[o++] = "acgt"[test_rand(4)];
      out[o++] = in[j];
    }
    out[o] = '\0';
    needleman_wunsch_align(seqa, seqb, &scoring, full, r_full);
    alignment_band_counts(&tries0, &wfa0);
    needleman_wunsch_align(seqa, seqb, &scoring, compact, r_compact);
    alignment_band_counts(&tries, &wfa);
    // Narrowed by the exact score, the band still passes first time
    if(wfa > wfa0 && compact->trace_banded) {
      ASSERT(tries == tries0 + 1);
      sized++;
    }
    ASSERT(alignments_equal(r_full, r_compact));
    banded += compact->traced && compact->trace_banded;
    kept += compact->traced && compact->trace_strips != NULL;
  }

  // Bands held for many (they need the striped fill), kept in strips only
  // where forced
#ifdef ALN_STRIPED
  ASSERT(banded > 500);
  ASSERT(sized > 100);
#else
  ASSERT(banded == 0 && sized == 0);
#endif
  ASSERT(keep ? kept == 0 : kept == banded);
  alignment_band_keep_max_cells = keep_cells;
  alignment_band_wfa_budget = wfa_budget;

  alignment_block_min_cells = min_cells;
  alignment_block_min_cols = min_cols;
  alignment_band_block_min_cells = band_cells;
  alignment_band_block_min_cols = band_cols;
  alignment_max_threads = max_threads;
  alignment_band_strip_rows = strip_rows;
  alignment_band_min_width = min_width;
  alignment_seed_length = seed_length;
  alignment_seeds_enabled = seeds_on;
  needleman_wunsch_free(full);
  needleman_wunsch_free(compact);
  alignment_free(r_full);
  alignment_free(r_compact);
}

// The seed bound (alignment_seeds.c) is never under the best score the
// rest of an alignment can add from a cell, whatever gap the path is in
// there: checked at every cell of small pairs against a reverse fill that
// keeps the type of the step into the cell (a gap of the same type goes
// on without opening), with short seeds of various scorings
void simd_test_seeds()
{
  char seqa[401], seqb[561];
  scoring_t scoring;
  const unsigned int rates[] = {0, 10, 50, 150, 400};
  size_t it, i, j, checked = 0, built = 0;
  for(it = 0; it < 300; it++)
  {
    const char *alpha = it % 4 == 3 ? "acgtn" : "acgt";
    const int match = 1 + (int)test_rand(3);
    scoring_init(&scoring, match, -(int)test_rand(6), -(int)test_rand(8),
                 -(int)test_rand(4), false, false, false, false, false,
                 it % 5 == 4);
    test_rand_seq(seqa, 20 + test_rand(280), alpha);
    if(it % 9 == 8) test_rand_seq(seqb, 20 + test_rand(280), alpha);
    else test_mutate(seqb, seqa, rates[it % 5], alpha);
    const size_t m = strlen(seqa), n = strlen(seqb);
    if(it % 2) for(j = 0; j < n; j += 3) seqb[j] = (char)toupper(seqb[j]);
    sub_table_t t;
    ASSERT(sub_table_build(&t, &scoring, seqa, m, seqb, n));
    score_t smax = INT_MIN;
    for(i = 0; i < t.n * t.n; i++) smax = MAX2(smax, t.scores[i]);
    const size_t k = 2 + test_rand(5);
    seeds_t *sd = seeds_build(&t, &scoring, seqa, m, seqb, n, smax, k);
    if(sd == NULL) {
      if(it < 40)
        printf("seeds: none for match %d mismatch %d open %d extend %d "
               "case %d alphabet %s k %zu letters %zu\n", match,
               scoring.mismatch, scoring.gap_open, scoring.gap_extend,
               (int)(it % 5 == 4), alpha, k, t.n);
      sub_table_free(&t);
      continue;
    }
    built++;
    // F[step][i][j]: the best score of aligning a[i..] with b[j..] after a
    // step of that type into (i, j): 0 a pair, 1 a gap taking a letter of
    // a, 2 one taking a letter of b
    const int64_t NEG = INT64_MIN / 4;
    const int64_t o = scoring.gap_open, e = scoring.gap_extend;
    const size_t W = n + 1;
    int64_t *F = malloc(3 * (m + 1) * W * sizeof(int64_t));
    ASSERT(F != NULL);
    #define FS(s, i, j) F[((s) * (m + 1) + (i)) * W + (j)]
    for(i = m + 1; i-- > 0; ) for(j = n + 1; j-- > 0; ) {
      int64_t pair = NEG, ga = NEG, gb = NEG, s;
      if(i == m && j == n) { FS(0,i,j) = FS(1,i,j) = FS(2,i,j) = 0; continue; }
      if(i < m && j < n) {
        const int sc = t.scores[t.code[(unsigned char)seqb[j]] * t.n +
                                t.code[(unsigned char)seqa[i]]];
        pair = sc + FS(0,i+1,j+1);
      }
      if(i < m) ga = FS(1,i+1,j) + e;
      if(j < n) gb = FS(2,i,j+1) + e;
      for(s = 0; s < 3; s++)
        FS(s,i,j) = MAX3(pair, ga + (s == 1 ? 0 : o), gb + (s == 2 ? 0 : o));
    }
    for(i = 0; i <= m; i++) for(j = 0; j <= n; j++) {
      const int64_t exact = MAX3(FS(0,i,j), FS(1,i,j), FS(2,i,j));
      const int64_t bound = seeds_bound(sd, i, j);
      if(bound < exact)
        printf("seeds: cell %zu,%zu bound %lld under %lld (%zu x %zu)\n",
               i, j, (long long)bound, (long long)exact, m, n);
      ASSERT(bound >= exact);
      checked++;
    }
    #undef FS
    free(F);
    seeds_free(sd);
    sub_table_free(&t);
  }
  ASSERT(built > 150);
  ASSERT(checked > 100000);
}

void test_simd_fill(void)
{
  SUITE_START("Striped SIMD fill");
  RUN(simd_test_diff());
  RUN(simd_test_seeds());
  RUN(simd_test_compact_band(true));
  RUN(simd_test_compact_band(false));
  RUN(simd_test_compact());
  RUN(simd_test_tiny_blocks());
  RUN(simd_test_fill_blocks());
  RUN(simd_test_matches_scalar());
  RUN(simd_test_odd_penalties());
  RUN(simd_test_directional_pairs());
  SUITE_END();
}
