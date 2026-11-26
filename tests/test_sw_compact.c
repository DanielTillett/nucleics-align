/*
 tests/test_sw_compact.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 Compact Smith-Waterman: rectangles, strips, narrow lanes and fast hits give the full search's hits
 */

#include "test_common.h"


/* Large compact Smith-Waterman alignments (forced here for small ones) find
 * the best hit from scans and a fill of just its rectangle: the same hits
 * as the full search, first and later ones, local hits in random flanks */
void simd_test_compact_rect()
{
  size_t min_cells = alignment_block_min_cells;
  size_t min_cols = alignment_block_min_cols;
  size_t band_cells = alignment_band_block_min_cells;
  size_t band_cols = alignment_band_block_min_cols;
  int max_threads = alignment_max_threads;
  size_t rect_cells = smith_waterman_rect_min_cells;
  sw_aligner_t *full = smith_waterman_new(), *rect = smith_waterman_new();
  alignment_t *r_full = alignment_create(256), *r_rect = alignment_create(256);
  char core[701], hit[1401], seqa[2001], seqb[3001];
  scoring_t scoring;
  const unsigned int rates[] = {0, 10, 60, 200};
  size_t i, h, rects = 0, hits = 0;

  smith_waterman_set_compact(rect, true);
  smith_waterman_rect_min_cells = 0;
  alignment_block_min_cells = 1;
  alignment_block_min_cols = 8;
  alignment_band_block_min_cells = 1;
  alignment_band_block_min_cols = 8;

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
      int extend = (int)test_rand(6) - 3;
      scoring_init(&scoring, 1 + (int)test_rand(4), -(int)test_rand(8),
                   (int)test_rand(9) - 4 - extend, extend, false, false,
                   false, false, false, false);
    }
    else scoring_init(&scoring, 1 + (int)test_rand(4), -(int)test_rand(8),
                      -(int)test_rand(8), -(int)test_rand(4), false, false,
                      false, false, false, false);
    // A related piece somewhere in each of two longer sequences
    test_rand_seq(core, 1 + test_rand(700), alpha);
    test_mutate(hit, core, rates[i % 4], alpha);
    size_t fa = test_rand(600), fb = test_rand(600), ta = test_rand(600);
    test_rand_seq(seqa, fa, alpha);
    strcat(seqa, core);
    test_rand_seq(seqa + strlen(seqa), ta, alpha);
    test_rand_seq(seqb, fb, alpha);
    strcat(seqb, hit);
    test_rand_seq(seqb + strlen(seqb), test_rand(600), alpha);

    score_t min_score = 1 + (score_t)test_rand(3) * (score_t)test_rand(20);
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
  // fill). Hits that start in the top half of the rows so far are filled
  // whole (see alignment_sw_scan_short()): about 220 of the 600 use one
#ifdef ALN_STRIPED
  ASSERT(rects > 150);
#else
  ASSERT(rects == 0);
#endif
  ASSERT(hits > 500);

  smith_waterman_rect_min_cells = rect_cells;
  alignment_block_min_cells = min_cells;
  alignment_block_min_cols = min_cols;
  alignment_band_block_min_cells = band_cells;
  alignment_band_block_min_cols = band_cols;
  alignment_max_threads = max_threads;
  smith_waterman_free(full);
  smith_waterman_free(rect);
  alignment_free(r_full);
  alignment_free(r_rect);
}

/* Strip alignments (forced here, with short strips): the best hit from
 * rows kept by the scan and strips filled again, over the rectangle or the
 * whole matrix, is the same as the full search's, and later hits too */
void simd_test_strips()
{
  size_t min_cells = alignment_block_min_cells;
  size_t min_cols = alignment_block_min_cols;
  size_t band_cells = alignment_band_block_min_cells;
  size_t band_cols = alignment_band_block_min_cols;
  int max_threads = alignment_max_threads;
  size_t rect_cells = smith_waterman_rect_min_cells;
  size_t strip_cells = smith_waterman_strip_min_cells;
  sw_aligner_t *full = smith_waterman_new(), *strips = smith_waterman_new();
  alignment_t *r_full = alignment_create(256), *r_strips = alignment_create(256);
  char core[701], hit[1401], seqa[2001], seqb[3001];
  scoring_t scoring;
  const unsigned int rates[] = {0, 10, 60, 200};
  size_t i, h, hits = 0, used = smith_waterman_strip_alignments;

  smith_waterman_set_compact(strips, true);
  smith_waterman_strip_min_cells = 0;
  alignment_block_min_cells = 1;
  alignment_block_min_cols = 8;
  alignment_band_block_min_cells = 1;
  alignment_band_block_min_cols = 8;

  for(i = 0; i < 400; i++)
  {
    alignment_max_threads = 1 + (int)(i % 3) * 2;
    // Over the rectangle, or the whole matrix
    smith_waterman_rect_min_cells = i % 2 ? 0 : SIZE_MAX;
    smith_waterman_strip_rows = 1 + test_rand(20);
    const char *alpha = i % 2 ? "acgt" : "ac";
    if(i % 7 == 6) {
      scoring_system_BLOSUM62(&scoring);
      alpha = "ARNDCQEGHILKMFPSTWYV";
    }
    else scoring_init(&scoring, 1 + (int)test_rand(4), -(int)test_rand(8),
                      -(int)test_rand(8), -(int)test_rand(4), false, false,
                      false, false, false, false);
    test_rand_seq(core, 1 + test_rand(700), alpha);
    test_mutate(hit, core, rates[i % 4], alpha);
    test_rand_seq(seqa, test_rand(600), alpha);
    strcat(seqa, core);
    test_rand_seq(seqa + strlen(seqa), test_rand(600), alpha);
    test_rand_seq(seqb, test_rand(600), alpha);
    strcat(seqb, hit);
    test_rand_seq(seqb + strlen(seqb), test_rand(600), alpha);

    score_t min_score = 1 + (score_t)test_rand(3) * (score_t)test_rand(20);
    smith_waterman_set_min_score(full, min_score);
    smith_waterman_set_min_score(strips, min_score);
    smith_waterman_align(seqa, seqb, &scoring, full);
    smith_waterman_align(seqa, seqb, &scoring, strips);
    size_t max_hits = i % 3 == 0 ? 4 : 1;
    for(h = 0; h < max_hits; h++) {
      int got_full = smith_waterman_fetch(full, r_full);
      int got_strips = smith_waterman_fetch(strips, r_strips);
      ASSERT(got_full == got_strips);
      if(!got_full || !got_strips) break;
      ASSERT(alignments_equal(r_full, r_strips));
      hits++;
    }
  }
#ifdef ALN_STRIPED
  ASSERT(smith_waterman_strip_alignments - used > 300); // strips were used
#endif
  ASSERT(hits > 350);

  smith_waterman_strip_rows = 0;
  smith_waterman_strip_min_cells = strip_cells;
  smith_waterman_rect_min_cells = rect_cells;
  alignment_block_min_cells = min_cells;
  alignment_block_min_cols = min_cols;
  alignment_band_block_min_cells = band_cells;
  alignment_band_block_min_cols = band_cols;
  alignment_max_threads = max_threads;
  smith_waterman_free(full);
  smith_waterman_free(strips);
  alignment_free(r_full);
  alignment_free(r_strips);
}

/* The rectangle's first scan gives up on a hit spanning the sequences, and
 * otherwise finds what the full scan finds */
void simd_test_scan_short()
{
  static char a[3001], b[6001];
  scoring_t scoring;
  score_t best, best_s;
  size_t x, y, xs, ys, i;
  scoring_init(&scoring, 1, -2, -4, -1, false, false, false, false, false,
               false);
  for(i = 0; i < 20; i++) {
    // A shared sequence (a hit from start to end), then unrelated ones
    test_rand_seq(a, 1000 + test_rand(2000), "acgt");
    if(i % 2) test_rand_seq(b, 1000 + test_rand(2000), "acgt");
    else test_mutate(b, a, 20, "acgt");
    const size_t la = strlen(a), lb = strlen(b);
    bool ok = alignment_sw_scan(a, la, b, lb, &scoring, 0, NULL, &best, &x, &y);
    bool ok_s = alignment_sw_scan_short(a, la, b, lb, &scoring, NULL,
                                        &best_s, &xs, &ys);
#ifdef ALN_STRIPED
    ASSERT(ok);
    ASSERT(ok_s == (i % 2 == 1)); // gives up only on the shared one
#endif
    if(ok && ok_s) ASSERT(best == best_s && x == xs && y == ys);
  }
}

/* Smith-Waterman fills and scans with 8- or 16-bit lanes (bits) give the
 * same results as 32-bit ones and the score matrices, including when
 * scores grow past the lanes' range part-way (the fill starts again with
 * wider lanes: 8-bit fills with 16-bit lanes, or with 32-bit ones in some
 * pairs), on one thread and on column blocks */
void simd_test_narrow(int bits)
{
  size_t min_cells = alignment_block_min_cells;
  size_t min_cols = alignment_block_min_cols;
  size_t band_cells = alignment_band_block_min_cells;
  size_t band_cols = alignment_band_block_min_cols;
  int max_threads = alignment_max_threads;
  sw_aligner_t *full = smith_waterman_new(), *compact = smith_waterman_new();
  alignment_t *r_full = alignment_create(256), *r_compact = alignment_create(256);
  char seqa[1201], seqb[2201], core[1001];
  scoring_t scoring;
  const unsigned int rates[] = {0, 10, 60, 300};
  const int top = bits == 8 ? INT8_MAX : INT16_MAX;
  size_t i, past = 0, within = 0;
  size_t fills0, retries0, fills1, retries1;

  smith_waterman_set_compact(compact, true);
  alignment_block_min_cells = 1;
  alignment_block_min_cols = 8;
  alignment_band_block_min_cells = 1;
  alignment_band_block_min_cols = 8;
  if(bits == 8) alignment_byte_counts(&fills0, &retries0);
  else alignment_narrow_counts(&fills0, &retries0);

  for(i = 0; i < 400; i++)
  {
    alignment_max_threads = 1 + (int)(i % 4);
    // The lanes under test; 8-bit fills that pass their range go on with
    // 32-bit lanes in some pairs
    alignment_byte_enabled = bits == 8;
    alignment_narrow_enabled = bits == 16 || i % 3 != 0;
    // Large scores in some, so the best passes the range less the highest
    // score (with 8 bits, scores small enough for the lanes in most)
    int match = bits == 8 ? (i % 5 ? 1 + (int)test_rand(4) : 1 + (int)test_rand(60))
                          : (i % 2 ? 1 + (int)test_rand(400) : 1 + (int)test_rand(5));
    // With 8 bits, penalties past the lanes' range in some (16-bit lanes
    // then, or 32-bit ones)
    int pen = bits == 8 && i % 7 != 0 ? 40 : 300;
    scoring_init(&scoring, match, -(int)test_rand((unsigned int)match * 2 + 2),
                 -(int)test_rand((unsigned int)pen),
                 -(int)test_rand((unsigned int)pen / 3), false, false,
                 false, false, false, false);
    test_rand_seq(core, test_rand(1000), i % 3 ? "acgt" : "ac");
    test_rand_seq(seqa, test_rand(100), "acgt");
    strcat(seqa, core);
    test_mutate(seqb, core, rates[i % 4], "acgt");
    test_rand_seq(seqb + strlen(seqb), test_rand(100), "acgt");

    smith_waterman_set_min_score(full, 1);
    smith_waterman_set_min_score(compact, 1);
    smith_waterman_align(seqa, seqb, &scoring, full);
    smith_waterman_align(seqa, seqb, &scoring, compact);
    int got_full = smith_waterman_fetch(full, r_full);
    int got_compact = smith_waterman_fetch(compact, r_compact);
    ASSERT(got_full == got_compact);
    if(got_full && got_compact) {
      ASSERT(alignments_equal(r_full, r_compact));
      if(r_full->score > top - match) past++;
      else within++;
    }

    // Scans: the best cell, then the highest cells with that score, with
    // the narrow lanes and with 32-bit ones only
    score_t bestn, best32, tn, t32;
    size_t xn, yn, x32, y32, txn, tyn, tx32, ty32;
    size_t la = strlen(seqa), lb = strlen(seqb);
    bool narrow = alignment_narrow_enabled;
    bool okn = alignment_sw_scan(seqa, la, seqb, lb, &scoring, 0, NULL,
                                 &bestn, &xn, &yn);
    alignment_narrow_enabled = alignment_byte_enabled = false;
    bool ok32 = alignment_sw_scan(seqa, la, seqb, lb, &scoring, 0, NULL,
                                  &best32, &x32, &y32);
#ifdef ALN_STRIPED
    ASSERT(okn == (la > 0 && lb > 0));
#endif
    ASSERT(okn == ok32);
    if(okn && ok32) {
      ASSERT(bestn == best32 && xn == x32 && yn == y32);
      if(best32 > 0) {
        ok32 = alignment_sw_scan(seqa, la, seqb, lb, &scoring, best32, NULL,
                                 &t32, &tx32, &ty32);
        alignment_byte_enabled = bits == 8;
        alignment_narrow_enabled = narrow;
        okn = alignment_sw_scan(seqa, la, seqb, lb, &scoring, best32, NULL,
                                &tn, &txn, &tyn);
        ASSERT(okn == ok32);
        ASSERT(!okn || (tn == t32 && txn == tx32 && tyn == ty32));
      }
    }
  }
  alignment_narrow_enabled = alignment_byte_enabled = true;

  // Edges of the 8-bit range. Identical sequences whose cells are all
  // multiples of 32 step from 96 (just over the limit of 95) straight to
  // 128, past 127, so the range must be checked at the right row; other
  // match scores near the limit. Then penalties beyond -128 (8-bit lanes
  // would read -300 as -44), which need wider lanes: on short sequences, so
  // the best stays in the 8-bit range and no retry would hide a wrong score
  // (b: a with its 21st base deleted, or changed)
  static const int edges[][5] = {{32, -128, -32, -32, 0},
                                 {43, -50, -50, -25, 0},
                                 {61, -50, -50, -25, 0},
                                 {3, -2, -299, -1, 1}, {3, -2, -2, -300, 1},
                                 {3, -300, -50, -50, 2}};
  // One thread: on narrow column blocks the best can jump where a hit
  // enters a block, which starts the fill again with 16-bit lanes and
  // would hide a wrong 8-bit score
  alignment_max_threads = 1;
  for(i = 0; bits == 8 && i < sizeof(edges) / sizeof(edges[0]); i++) {
    scoring_init(&scoring, edges[i][0], edges[i][1], edges[i][2], edges[i][3],
                 false, false, false, false, false, false);
    test_rand_seq(seqa, edges[i][4] ? 40 : 400, "acgt");
    strcpy(seqb, seqa);
    if(edges[i][4] == 1) memmove(seqb + 20, seqb + 21, strlen(seqb + 20));
    if(edges[i][4] == 2) seqb[20] = seqa[20] == 'a' ? 'c' : 'a';
    smith_waterman_align(seqa, seqb, &scoring, full);
    smith_waterman_align(seqa, seqb, &scoring, compact);
    ASSERT(smith_waterman_fetch(full, r_full) &&
           smith_waterman_fetch(compact, r_compact) &&
           alignments_equal(r_full, r_compact));
    score_t bestn, best32;
    size_t xn, yn, x32, y32;
    size_t la = strlen(seqa), lb = strlen(seqb);
    bool okn = alignment_sw_scan(seqa, la, seqb, lb, &scoring, 0, NULL,
                                 &bestn, &xn, &yn);
    alignment_byte_enabled = alignment_narrow_enabled = false;
    bool ok32 = alignment_sw_scan(seqa, la, seqb, lb, &scoring, 0, NULL,
                                  &best32, &x32, &y32);
    alignment_byte_enabled = alignment_narrow_enabled = true;
    ASSERT(okn == ok32 &&
           (!okn || (bestn == best32 && xn == x32 && yn == y32)));
  }

  // Both sides of the range were tested
  ASSERT(past > 50);
  ASSERT(within > 50);
  // The narrow fills ran, and every best hit past the range was filled
  // again with wider lanes (the scans may stop earlier)
  if(bits == 8) alignment_byte_counts(&fills1, &retries1);
  else alignment_narrow_counts(&fills1, &retries1);
#ifdef ALN_STRIPED
  ASSERT(fills1 - fills0 > 600); // a fill and a scan for most pairs
  ASSERT(retries1 - retries0 >= past);
#else
  ASSERT(fills1 == fills0 && retries1 == retries0);
#endif

  alignment_block_min_cells = min_cells;
  alignment_block_min_cols = min_cols;
  alignment_band_block_min_cells = band_cells;
  alignment_band_block_min_cols = band_cols;
  alignment_max_threads = max_threads;
  smith_waterman_free(full);
  smith_waterman_free(compact);
  alignment_free(r_full);
  alignment_free(r_compact);
}

/* A gap scoring above 0 can open a local hit, which the rectangle's
 * reverse scan would miss: the compact aligner gives the same hit (the
 * checked aligner refuses such scoring; the unchecked one is the way in) */
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

void test_sw_fastlocal(void)
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
  ASSERT(smith_waterman_fetch(sw, r1) == 1);
  ASSERT(r1->score == 6);

  // 2. Codex multi-hit counterexample: ACAA vs ACCCAA (+2/-1/0/-1)
  scoring_t sc2;
  scoring_init(&sc2, 2, -1, 0, -1, false, false, false, false, false, true);
  smith_waterman_align("ACAA", "ACCCAA", &sc2, sw);
  ASSERT(smith_waterman_fetch(sw, r1) == 1);
  ASSERT(r1->score == 6);
  int f2 = smith_waterman_fetch(sw, r2);
  ASSERT(f2 == 0 || (f2 == 1 && r2->score <= r1->score));

  // 3. 3000 bp pair: exact optimal score match vs default Gotoh + 2nd hit
  scoring_t sc3;
  scoring_init(&sc3, 1, -2, -4, -1, false, false, false, false, false, true);
  char *a = malloc(3100), *b = malloc(6001);
  ASSERT(a && b);
  test_rand_seq(a, 3000, "ACGT");
  test_mutate(b, a, 50, "ACGT");

  smith_waterman_fast_rect = false;
  smith_waterman_align(a, b, &sc3, sw);
  ASSERT(smith_waterman_fetch(sw, rdef) == 1);

  smith_waterman_fast_rect = true;
  smith_waterman_align(a, b, &sc3, sw);
  ASSERT(smith_waterman_fetch(sw, r1) == 1);
  ASSERT(r1->score == rdef->score);
  f2 = smith_waterman_fetch(sw, r2);
  ASSERT(f2 == 0 || (f2 == 1 && r2->score <= r1->score));

  smith_waterman_fast_rect = false;
  smith_waterman_rect_min_cells = saved_min;
  alignment_free(r1);
  alignment_free(r2);
  alignment_free(rdef);
  smith_waterman_free(sw);
  free(a);
  free(b);
}

void test_sw_compact(void)
{
  SUITE_START("Compact Smith-Waterman");
  RUN(simd_test_compact_rect());
  RUN(simd_test_scan_short());
  RUN(simd_test_strips());
  RUN(simd_test_rect_gap_score());
  RUN(simd_test_narrow(16));
  RUN(simd_test_narrow(8));
  RUN(test_sw_fastlocal());
  SUITE_END();
}
