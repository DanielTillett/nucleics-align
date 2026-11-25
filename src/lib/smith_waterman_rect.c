/*
 smith_waterman_rect.c
 O(len_b) rectangle scan, strip checkpoint reconstruction, and fast-local path.
 */

#include "smith_waterman_internal.h"

// The best hit of the region a[x0 .. x0+rx) against b[y0 .. y0+ry) (of the
// whole sequences) without its traceback bytes (see
// smith_waterman_strip_min_cells): one scan finds the best cell, keeping
// every strip_rows-th row; sw_follow_strips() walks the hit. want: the best
// score a rectangle's region must have, at its last cell (0: any). Sets
// *used where it applies (else nothing is changed); ALN_ERR_NOMEM if out
// of memory
aln_status_t sw_strips_fill(sw_aligner_t *sw, size_t x0, size_t y0,
                               size_t rx, size_t ry, const scoring_t *scoring,
                               score_t want, bool *used)
{
  aligner_t *aligner = &sw->aligner;
  *used = false;
  if(rx == 0 || ry == 0) return ALN_OK; // no cells: no hit, nothing to do
  // Rows per strip near sqrt(12 x rows): the strip's bytes (k bytes per
  // column) and the rows kept (12 bytes per column each) then take about
  // the same memory
  size_t k = smith_waterman_strip_rows;
  if(k == 0) {
    for(k = 16; k * k < 12 * ry; k++) {}
  }
  k = MIN2(k, ry);
  const size_t kept = ry / k + 1;
  if(kept > SIZE_MAX / 3 / sizeof(score_t) / (rx + 1))
    return ALN_OK;
  score_t *ckpt = aln_malloc(kept * 3 * (rx + 1) * sizeof(score_t));
  if(ckpt == NULL) return ALN_ERR_NOMEM;
  score_t best;
  size_t bx, by;
  if(!alignment_sw_rows(sw->seq_a + x0, rx, sw->seq_b + y0, ry, scoring,
                        &aligner->opts, k, ckpt, &best, &bx, &by) ||
     (want > 0 && (best != want || bx != rx || by != ry))) {
    free(ckpt);
    return ALN_OK;
  }
  free(sw->ckpt);
  sw->ckpt = ckpt;
  sw->strip_rows = k;
  sw_describe(sw, scoring, sw->seq_a + x0, sw->seq_b + y0, rx, ry, best,
              by * (rx + 1) + bx, x0, y0);
  sw->strips = true;
  sw->full_best = (y0 + by) * (sw->len_a + 1) + x0 + bx;
  __atomic_add_fetch(&smith_waterman_strip_alignments, 1, __ATOMIC_RELAXED);
  *used = true;
  return ALN_OK;
}

// Follow the best hit of a strip alignment back from its last cell, as
// _follow_hit() follows traceback bytes: from the bottom, each strip of
// rows gets its bytes again (up to the hit's column) when the walk enters
// it. Returns 1, or -1 if out of memory
int sw_follow_strips(sw_aligner_t *sw, alignment_t *result)
{
  const aligner_t *aligner = &sw->aligner;
  const scoring_t *scoring = aligner->scoring;
  const sub_table_t *t = &sw->table;
  const size_t w = aligner->score_width, k = sw->strip_rows;
  const char *a = aligner->seq_a, *b = aligner->seq_b;
  const size_t end_x = aligner->best_index % w, end_y = aligner->best_index / w;
  size_t x = end_x, y = end_y, n = 0, s0 = 0, xs = 0;
  enum Matrix m = MATCH;
  score_t score = aligner->best_score;

  uint8_t *ops = aln_malloc(end_x + end_y + 1); // the matrix of each step
  uint8_t *trace = aln_malloc(k * end_x);
  if(ops == NULL || trace == NULL) {
    free(ops);
    free(trace);
    return -1;
  }

  while(score > 0) {
    if(n == 0 || y < s0) {
      // The strip holding row y, from the row kept above it
      const size_t kept = (y - 1) / k;
      s0 = kept * k + 1;
      xs = x;
      if(!alignment_sw_trace_rows(a, w - 1, b, aligner->score_height - 1,
                                  scoring, t, &aligner->opts,
                                  aligner->best_score,
                                  sw->ckpt + kept * 3 * w, s0, y, xs, trace)) {
        free(ops);
        free(trace);
        return -1;
      }
    }
    // As _sw_reverse_move() steps on traceback bytes
    const int next = (trace[(y - s0) * xs + x - 1] >> (2 * (int)m)) & 3;
    ops[n++] = (uint8_t)m;
    sw_trace_step(t, scoring, a, b, next, &m, &x, &y, &score,
                  "_follow_strips");
  }
  free(trace);

  if(!alignment_reserve(result, n)) {
    free(ops);
    return -1;
  }
  // From the first cell: the steps in reverse
  size_t i, xa = x, yb = y;
  for(i = 0; i < n; i++) {
    const enum Matrix op = (enum Matrix)ops[n-1-i];
    if(op != GAP_A) xa++;
    if(op != GAP_B) yb++;
    result->result_a[i] = op == GAP_A ? '-' : a[xa-1];
    result->result_b[i] = op == GAP_B ? '-' : b[yb-1];
  }
  free(ops);
  result->result_a[n] = result->result_b[n] = '\0';
  result->length = n;
  result->score = aligner->best_score;
  result->pos_a = x;
  result->pos_b = y;
  result->len_a = end_x - x;
  result->len_b = end_y - y;
  return 1;
}

// Verify that a global alignment of a[x0..x0+rx) vs b[y0..y0+ry) is a valid
// optimal local hit scoring exactly best with strictly positive local prefixes
static bool _sw_verify_local_hit(const sub_table_t *t, const scoring_t *scoring,
                                 const alignment_t *res, size_t rx, size_t ry,
                                 score_t best)
{
  if(res->length == 0 || res->score != best) return false;
  if(res->result_a[0] == '-' || res->result_b[0] == '-') return false;
  if(res->result_a[res->length - 1] == '-' ||
     res->result_b[res->length - 1] == '-') return false;
  const long long open = (long long)scoring->gap_open + scoring->gap_extend;
  const long long ext = (long long)scoring->gap_extend;
  long long sc = 0;
  size_t ia = 0, ib = 0, k;
  int st = 0; // 0: match, 1: gap_a, 2: gap_b
  for(k = 0; k < res->length; k++) {
    char ca = res->result_a[k], cb = res->result_b[k];
    if(ca == '-' && cb == '-') return false;
    if(ca == '-') {
      sc += (st == 1 ? ext : open);
      st = 1;
      ib++;
    } else if(cb == '-') {
      sc += (st == 2 ? ext : open);
      st = 2;
      ia++;
    } else {
      sc += t->scores[t->code[(unsigned char)cb] * t->n +
                      t->code[(unsigned char)ca]];
      st = 0;
      ia++;
      ib++;
    }
    if(sc <= 0) return false;
  }
  return ia == rx && ib == ry && sc == (long long)best;
}

bool sw_align_fast_rect(sw_aligner_t *sw, const char *a, const char *b,
                                size_t len_a, size_t len_b,
                                const scoring_t *scoring)
{
  aligner_t *aligner = &sw->aligner;
  score_t best = 0, best2 = 0;
  size_t bx = 0, by = 0, rx = 0, ry = 0, i;

  if(scoring->gap_extend > 0 || scoring->gap_open + scoring->gap_extend > 0 ||
     scoring->no_start_gap_penalty || scoring->no_end_gap_penalty ||
     scoring->no_gaps_in_a || scoring->no_gaps_in_b ||
     scoring->no_gaps_in_longest || scoring->no_mismatches)
    return false;

  if(!alignment_sw_scan_hint(a, len_a, b, len_b, scoring, 0, sw->min_score,
                             &aligner->opts, &best, &bx, &by))
    return false;

  if(best < sw->min_score) {
    // No hit: nothing to fill, but describe the sequences as a fill would
    // (best_index 0: not an earlier alignment's cell)
    sw_describe(sw, scoring, a, b, len_a, len_b, best, 0, 0, 0);
    return true;
  }

  // Reverse scan with target=0 locates a genuine jointly-located optimal start
  char *ra = aln_malloc(bx + by);
  if(ra == NULL) return false;
  char *rb = ra + bx;
  for(i = 0; i < bx; i++) ra[i] = a[bx - 1 - i];
  for(i = 0; i < by; i++) rb[i] = b[by - 1 - i];
  bool ok = alignment_sw_scan_hint(ra, bx, rb, by, scoring, 0, best,
                                   &aligner->opts, &best2, &rx, &ry);
  free(ra);
  if(!ok || best2 != best || rx == 0 || ry == 0) return false;

  const size_t x0 = bx - rx, y0 = by - ry;
  const size_t wfa_budget = aligner->opts.wfa_budget != 0
                            ? aligner->opts.wfa_budget
                            : alignment_wfa_default_budget(rx, ry);

  aln_status_t st = alignment_wfa_align(&sw->table, scoring, a + x0, rx,
                                        b + y0, ry, wfa_budget, &sw->fast_res);
  if(st != ALN_OK) {
    st = hirschberg_align_checked(a + x0, b + y0, rx, ry, scoring,
                                  &aligner->opts, &sw->fast_res, NULL);
  }
  if(st != ALN_OK ||
     !_sw_verify_local_hit(&sw->table, scoring, &sw->fast_res, rx, ry, best))
    return false;

  sw->fast_res.pos_a = x0;
  sw->fast_res.pos_b = y0;
  sw->fast_res.len_a = rx;
  sw->fast_res.len_b = ry;
  sw->fast_hit = true;
  sw->full_best = by * (len_a + 1) + bx;
  sw_describe(sw, scoring, a, b, len_a, len_b, best, sw->full_best, x0, y0);
  return true;
}

// Best hit of a large compact alignment in little memory (as SSW finds it,
// but with the same hit as the full search). A scan finds the best score
// and its first cell (bx, by). A local path ending there and scoring best
// is a best path of the reversed prefixes, so a scan of a[0..bx) and
// b[0..by) reversed finds where every best path ending at (bx, by) starts:
// the rectangle from the lowest start column and row to (bx, by) holds them
// all. A fill of just that rectangle gives exact values along every best
// path, and no cell outside it could tie with one on them, so the
// traceback takes the same steps. The rectangle is as big as the hit: this
// saves memory and time when hits are short against the sequences.
// That needs every hit to start with a pair: the traceback stops at a score
// of 0, so a hit can open with a gap only if a gap can score above 0.
// Returns false where the scans do not apply, or if out of memory (the
// whole fill then reports it)
bool sw_align_rect(sw_aligner_t *sw, const char *a, const char *b,
                           size_t len_a, size_t len_b,
                           const scoring_t *scoring)
{
  aligner_t *aligner = &sw->aligner;
  score_t best, best2;
  size_t bx, by, rx, ry, i;

  // Only where no gap scores above 0, and where the full matrices would fit
  // (aligner_align() reports sequences too long for them)
  if(scoring->gap_extend > 0 || scoring->gap_open + scoring->gap_extend > 0 ||
     len_b+1 > SIZE_MAX / (len_a+1) / sizeof(score_t))
    return false;

  // Gives up on a hit that spans the sequences: the whole fill is quicker
  if(!alignment_sw_scan_short(a, len_a, b, len_b, scoring, &aligner->opts,
                              &best, &bx, &by))
    return false;

  if(best < sw->min_score) {
    // No hit: nothing to fill, but describe the sequences as a fill would
    // (best_index 0: not an earlier alignment's cell)
    sw_describe(sw, scoring, a, b, len_a, len_b, best, 0, 0, 0);
    return true;
  }

  // A hit scoring best has at least best / smax aligned pairs. If that is
  // most of the prefix up to (bx, by), along its longer side, fill the
  // whole prefix (exact too: nothing after (bx, by) changes it) rather than
  // scan for the rectangle. (Along the shorter side alone, a 100-letter
  // hit at (100, 100000) would have the whole 100 x 100000 prefix filled)
  const score_t smax = sw->table.smax;
  if(smax > 0 && (double)best >= 0.5 * smax * (double)MAX2(bx, by)) {
    rx = bx;
    ry = by;
  }
  else
  {
    // Reversed prefixes up to the end of the best hit
    char *ra = aln_malloc(bx + by);
    if(ra == NULL) return false;
    char *rb = ra + bx;
    for(i = 0; i < bx; i++) ra[i] = a[bx-1-i];
    for(i = 0; i < by; i++) rb[i] = b[by-1-i];
    bool ok = alignment_sw_scan(ra, bx, rb, by, scoring, best, &aligner->opts,
                                &best2, &rx, &ry);
    free(ra);
    if(!ok || rx == 0 || ry == 0) return false;
  }

  // Reversed cell (rx, ry) is the start pair (bx-rx+1, by-ry+1)
  const size_t x0 = bx - rx, y0 = by - ry;
  if((double)rx * (double)ry >= (double)aligner->opts.sw_strip_min_cells) {
    bool used;
    if(sw_strips_fill(sw, x0, y0, rx, ry, scoring, best, &used) != ALN_OK)
      return false; // the whole fill reports running out of memory
    if(used) return true;
  }
  if(aligner_align_checked(aligner, a + x0, b + y0, rx, ry, scoring,
                           1) != ALN_OK ||
     !aligner->traced || aligner->best_score != best ||
     aligner->best_index != ry * (rx + 1) + rx) return false;

  sw->rect = true;
  sw->off_a = x0;
  sw->off_b = y0;
  sw->full_best = by * (len_a + 1) + bx;
  return true;
}

