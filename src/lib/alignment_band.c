/*
 alignment_band.c
 Adaptive strip-pruned banded Needleman-Wunsch fill and trace strip caching.
 */

#define _DEFAULT_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdatomic.h>
#include <unistd.h>

#include "alignment_fill_internal.h"
#include "alignment_seeds.h"
#include "alignment_macros.h"
#include "alignment_alloc.h"
#include "alignment_testing.h"

// A path through diagonal d has at least I = |d| + |n-m-d| gap characters
// and at most (m+n-I)/2 aligned pairs, each scoring at most smax; opening
// gaps only lowers the score (gap_open <= 0)
size_t alignment_band_width(score_t smax, score_t gap_extend,
                            size_t m, size_t n, score_t lower)
{
  const double s = smax, h = gap_extend;
  const double gaps = m > n ? (double)(m - n) : (double)(n - m);
  const double step = 2 * (s - 2*h); // each diagonal further out
  const double room = s * ((double)m + n - gaps) + 2*h*gaps - 2.0*lower;
  if(step <= 0) return SIZE_MAX;
  if(room < 0) return 0;
  return room / step >= 1e15 ? SIZE_MAX : (size_t)(room / step);
}


#ifdef ALN_STRIPED
// Rows of a band's strip starting at row y0, and the columns c0+1 .. c1 its
// rows reach on diagonals x - y from dlo to dhi
static void alignment_strip(size_t y0, size_t strip, size_t len_i,
                            size_t len_j, ptrdiff_t dlo, ptrdiff_t dhi,
                            size_t *y1, size_t *c0, size_t *c1)
{
  *y1 = MIN2(len_j, y0 + strip - 1);
  ptrdiff_t lo = (ptrdiff_t)y0 + dlo, hi = (ptrdiff_t)*y1 + dhi;
  *c0 = lo <= 1 ? 0 : (size_t)lo - 1;
  *c1 = hi >= (ptrdiff_t)len_i ? len_i : (size_t)hi;
}

// A guess at the best global alignment score, meant to be no lower than it,
// for deciding whether a band is worth trying (the alignment never depends
// on it): words of a few letters of a are looked for near where they would
// lie in b, and the share found (each found with chance about p^K, for a
// share p of letters that match) gives p; the rest count as mismatches,
// with no gaps but one for the difference in length. Returns false for
// short sequences, where a first, narrow band costs little
static bool alignment_score_guess(const sub_table_t *t,
                                  const scoring_t *scoring,
                                  const char *a, size_t len_a,
                                  const char *b, size_t len_b, double *guess)
{
  enum { K = 12, SAMPLES = 128, WINDOW = 256 };
  if(len_a < 64 * K || len_b < 64 * K) return false;
  size_t s, i, j, k, found = 0;
  // Letters match as the scoring says (so soft-masked lower case against
  // upper case counts under case-insensitive scoring), through the table
  const unsigned char *code = t->code;
  const bool *is_match = t->is_match;
  const size_t n = t->n;
  for(s = 0; s < SAMPLES; s++) {
    i = (len_a - K) * s / (SAMPLES - 1);
    const size_t at = (size_t)((double)i * (double)(len_b - K) /
                               (double)(len_a - K));
    const size_t lo = at > WINDOW ? at - WINDOW : 0;
    const size_t hi = MIN2(len_b - K, at + WINDOW);
    for(j = lo; j <= hi; j++) {
      for(k = 0; k < K; k++)
        if(!is_match[code[(unsigned char)b[j+k]] * n +
                     code[(unsigned char)a[i+k]]]) break;
      if(k == K) { found++; break; }
    }
  }
  // p^K is the share found, taken high so the guess stays high (a band
  // wrongly skipped costs more than one wrongly tried: 50 kb at 3%, 0.16
  // -> 0.19 s): two of the largest standard errors 128 samples can have
  // (0.5 / sqrt(128) = 0.044), and one sample more
  const double f = MIN2(1.0, (double)found / SAMPLES + 0.1);
  double plo = 0, phi = 1, p, pk;
  for(s = 0; s < 40; s++) {
    p = (plo + phi) / 2;
    for(pk = 1, i = 0; i < K; i++) pk *= p;
    if(pk < f) plo = p; else phi = p;
  }
  p = (plo + phi) / 2;
  // Mean scores of matching and of mismatching pairs
  double match = 0, mismatch = 0;
  size_t nmatch = 0, nmismatch = 0;
  for(i = 0; i < t->n * t->n; i++) {
    if(t->is_match[i]) { match += t->scores[i]; nmatch++; }
    else { mismatch += t->scores[i]; nmismatch++; }
  }
  if(nmatch == 0) return false;
  match /= (double)nmatch;
  mismatch = nmismatch ? mismatch / (double)nmismatch : match;
  const size_t L = MIN2(len_a, len_b);
  const size_t D = len_a > len_b ? len_a - len_b : len_b - len_a;
  *guess = (double)L * (p * match + (1 - p) * mismatch) +
           (D ? (double)scoring->gap_open + (double)D * scoring->gap_extend
              : 0);
  return true;
}

// A band kept in strips (see alignment_band_keep_max_cells): what filling a
// strip again needs. Strip s holds rows 1 + s*strip .. and columns c0[s] ..
// c1 (the band's, see alignment_fill_banded()); its top row (row y0-1),
// columns c0[s] .. c1, is kept at rows + off[s], 3 x (c1 - c0[s] + 1) as
// MATCH, GAP_A, GAP_B
typedef struct {
  size_t strip, nstrips; // strip: the most rows in one
  size_t *y0s; // strip s holds rows y0s[s] .. y0s[s+1] - 1
  size_t *c0;
  score_t min;
  size_t *off;
  score_t *rows;
  score_t *col0; // column 0, 3 x score_height, for strips that reach it
  score_t *top; // 3 x score_width: a strip's top row by column
  score_t *lcol; // 3 x (strip + 1): its left column
  size_t loaded, loaded_y, loaded_x; // the strip in trace, to row, column
  sub_table_t table; // the pair's scores, built on the first reload
  bool has_table;
} trace_strips_t;

void trace_strips_free(void *ts_ptr)
{
  trace_strips_t *ts = ts_ptr;
  if(ts == NULL) return;
  free(ts->off);
  free(ts->y0s);
  free(ts->c0);
  free(ts->rows);
  free(ts->col0);
  free(ts->top);
  free(ts->lcol);
  if(ts->has_table) sub_table_free(&ts->table);
  free(ts);
}

// Fill strip s of a band again, rows y0 .. y (y in the strip) and its
// columns up to x, into trace from its kept top row
static bool trace_strips_fill(aligner_t *aligner, trace_strips_t *ts,
                              size_t s, size_t y, size_t x)
{
  const size_t w = aligner->score_width, h = aligner->score_height;
  const size_t len_i = w - 1, len_j = h - 1;
  const size_t y0 = ts->y0s[s], span = ts->off[s+1] - ts->off[s];
  const size_t c0 = ts->c0[s], c1 = MIN2(c0 + span / 3 - 1, x);
  const size_t width = c1 - c0;
  size_t r, m;
  // Its top row, by column (the kept columns run to the strip's own c1)
  for(m = 0; m < 3; m++)
    memcpy(ts->top + m*w + c0, ts->rows + ts->off[s] + m * (span / 3),
           (width + 1) * sizeof(score_t));
  const score_t *top[3] = {ts->top, ts->top + w, ts->top + 2*w};
  for(r = y0; r <= y; r++)
    aligner->trace_rows[r] = (ptrdiff_t)((r - y0) * width) -
                             (ptrdiff_t)(c0 + 1);
  // The scores table once, not per strip (building it scans both whole
  // sequences, and a 1 Mb traceback reloads hundreds of strips)
  if(!ts->has_table) {
    if(!sub_table_build(&ts->table, aligner->scoring, aligner->seq_a, len_i,
                        aligner->seq_b, len_j)) return false;
    ts->has_table = true;
  }
  alignment_block_t tmpl = {
    .aligner = aligner, .table = &ts->table, .min = ts->min,
    .c0 = c0, .width = width, .y0 = y0, .y1 = y,
    .r0 = {top[MATCH], top[GAP_A], top[GAP_B]}, .lstride = 1,
    .band = true,
    .trace = aligner->trace, .trace_rows = aligner->trace_rows};
  for(m = 0; m < 3; m++) {
    if(c0 == 0) tmpl.lc[m] = ts->col0 + m*h;
    else {
      // As the fill set it: the cell above the strip, then min
      score_t *lcm = ts->lcol + m * (ts->strip + 1);
      lcm[0] = top[m][c0];
      for(r = 1; r <= y - y0 + 1; r++) lcm[r] = ts->min;
      tmpl.lc[m] = lcm;
      tmpl.ly0 = y0 - 1;
    }
  }
  bool filled = alignment_fill_blocks(&tmpl);
  if(filled) {
    ts->loaded = s;
    ts->loaded_y = y;
    ts->loaded_x = c1;
  }
  else ts->loaded = ts->nstrips; // none
  return filled;
}

bool alignment_trace_load(aligner_t *aligner, size_t x, size_t y)
{
  trace_strips_t *ts = aligner->trace_strips;
  if(ts == NULL || y == 0) return true;
  // The strip holding row y: the last starting at or above it
  size_t lo = 0, hi = ts->nstrips - 1;
  while(lo < hi) {
    const size_t mid = (lo + hi + 1) / 2;
    if(ts->y0s[mid] <= y) lo = mid; else hi = mid - 1;
  }
  const size_t s = lo;
  if(s == ts->loaded && y <= ts->loaded_y && x <= ts->loaded_x) return true;
  return trace_strips_fill(aligner, ts, s, y, x);
}

// The diagonals x - y the strips below row y need, from row y's cells
// (columns c0+1 .. c1, plus real column 0 when c0 == 0) and a lower bound lb
// on the best score. A cell scoring v can end no better than its bound,
// v + smax min(i, j) + e |i - j| (i, j the letters left of each sequence,
// e the gap extension). A step towards diagonal len_i - len_j does not
// raise the bound and one away displaces a matched pair and adds two gap
// extensions, lowering it by at least smax - 2e, so from a cell whose bound
// is lb + k, cells with bounds of lb or more lie only on the diagonals from
// it to len_i - len_j and k / (smax - 2e) beyond. A cell with a bound under
// lb is on no best path, and could not tie with a cell on one, so leaving it
// out changes no cell on a best path (the computational volume A*PA2 narrows
// its band to: Groot Koerkamp 2024, WABI). Narrows *lo .. *hi to those
// diagonals; no change if the gap extension is 0, or no cell's bound reaches lb
static void alignment_band_narrow(score_t *const row[3],
                                  score_t *const col0[3], size_t y,
                                  size_t c0, size_t c1,
                                  size_t len_i, size_t len_j, score_t smax,
                                  int gap_extend, score_t lb,
                                  ptrdiff_t *lo, ptrdiff_t *hi)
{
  if(gap_extend >= 0) return;
  const long long e2 = (long long)smax - 2LL * gap_extend;
  const ptrdiff_t D = (ptrdiff_t)len_i - (ptrdiff_t)len_j;
  ptrdiff_t nlo = PTRDIFF_MAX, nhi = PTRDIFF_MIN;
  size_t x;
  if(c0 == 0) {
    const long long v = MAX3(col0[MATCH][y], col0[GAP_A][y], col0[GAP_B][y]);
    const size_t ri = len_i, rj = len_j - y;
    const long long bound = v + (long long)smax * (long long)MIN2(ri, rj) +
                            (long long)gap_extend *
                            (long long)(ri > rj ? ri - rj : rj - ri);
    if(bound >= lb) {
      const ptrdiff_t d = -(ptrdiff_t)y;
      const ptrdiff_t k = (ptrdiff_t)MIN2((bound - lb) / e2,
                                          (long long)(len_i + len_j));
      nlo = MIN2(nlo, MIN2(d, D) - k);
      nhi = MAX2(nhi, MAX2(d, D) + k);
    }
  }
  for(x = c0 + 1; x <= c1; x++) {
    const long long v = MAX3(row[MATCH][x], row[GAP_A][x], row[GAP_B][x]);
    const size_t ri = len_i - x, rj = len_j - y;
    const long long bound = v + (long long)smax * (long long)MIN2(ri, rj) +
                            (long long)gap_extend *
                            (long long)(ri > rj ? ri - rj : rj - ri);
    if(bound < lb) continue;
    const ptrdiff_t d = (ptrdiff_t)x - (ptrdiff_t)y;
    const ptrdiff_t k = (ptrdiff_t)MIN2((bound - lb) / e2,
                                        (long long)(len_i + len_j));
    nlo = MIN2(nlo, MIN2(d, D) - k);
    nhi = MAX2(nhi, MAX2(d, D) + k);
  }
  if(nlo > nhi) return;
  *lo = MAX2(*lo, nlo);
  *hi = MIN2(*hi, nhi);
}

// Compact Needleman-Wunsch fill of only the cells near the main diagonals.
// A path through a cell far from them needs many gap characters, which caps
// its score (alignment_band_width()). The rows run in strips, each over
// just the columns its rows' band reaches, with min for the column on its
// left and for top-row cells the strip above did not fill: cells next to
// those come out too low, never too high, and cells whose best path stays
// in the band are exact. If the best score found (a real alignment's)
// needs no wider band, every path scoring as much is inside it, so no cell
// outside could tie with one on a best path, and the traceback takes the
// same steps as with every cell filled. Otherwise the band is widened and
// filled again. Returns false, with nothing done, once a band would not be
// faster than the whole fill on its threads. A band of more than
// band_keep_max_cells cells keeps one strip's bytes at a time, and each
// strip's top row (aligner->trace_strips, see alignment_trace_load()), in
// strips of at most about sqrt(12 x rows) rows
fill_result_t alignment_fill_banded(aligner_t *aligner,
                                           const sub_table_t *table,
                                           score_t min, score_t *const row0[3],
                                           score_t *const col0[3])
{
  const scoring_t *scoring = aligner->scoring;
  const size_t w = aligner->score_width, h = aligner->score_height;
  const size_t len_i = w-1, len_j = h-1;
  const ptrdiff_t D = (ptrdiff_t)len_i - (ptrdiff_t)len_j;
  // The whole fill's threads, or one for a narrow band
  const double full_threads =
    (double)alignment_fill_threads(&aligner->opts, len_i, len_j, true);
  size_t i, m, y0, y1, c0, c1;

  if(!aligner->opts.band || scoring->gap_open > 0 ||
     scoring->gap_extend > 0) return FILL_SKIPPED;

  score_t smax = 0;
  for(i = 0; i < table->n * table->n; i++) smax = MAX2(smax, table->scores[i]);

  // A band must cost under this much of the whole fill, whose difference-
  // encoded lanes (where they apply) take about 0.6 of the time per cell
  // of the band's absolute ones: banding at 20 kb and 3-20% divergence
  // found the band too narrow and filled the whole matrix anyway, 0.033 ->
  // 0.044 s
  score_t floor, ceiling;
  const double pays = aligner->opts.diff &&
                      alignment_diff_bits(table, scoring, &floor, &ceiling)
                      ? 0.35 : 0.6;

  // Quick reject for small single-threaded inputs where even the narrowest
  // allowed band with the minimum strip height cannot beat a full compact fill
  if(full_threads <= 1.0) {
    const double min_cols = (double)(D < 0 ? -D : D) +
                            2.0 * (double)aligner->opts.band_min_width + 1.0;
    const size_t min_strip = aligner->opts.band_strip_rows
                             ? aligner->opts.band_strip_rows : 64;
    if(min_cols + (double)min_strip >= pays * (double)len_i)
      return FILL_SKIPPED;
  }

  // The best score itself where the wavefront algorithm finds it within
  // its budget (close pairs): then the band it needs is filled once (50 kb
  // at 0.1%: 0.040 -> 0.018 s). Otherwise skip the band if even a high
  // guess at the best score needs one too wide to pay (the cost rule
  // below): the first, narrow band would only find that out (20 kb at
  // 3-20% divergence: 0.044 -> 0.033 s)
  double guess;
  score_t exact;
  const size_t budget = aligner->opts.band_wfa_budget;
  const bool have_exact = budget > 0 && len_i <= SIZE_MAX / budget &&
    alignment_wfa_score(table, scoring, aligner->seq_a, len_i,
                        aligner->seq_b, len_j, budget * len_i, &exact);
  if(have_exact) guess = exact;
  if((have_exact ||
      alignment_score_guess(table, scoring, aligner->seq_a, len_i,
                            aligner->seq_b, len_j, &guess)) &&
     guess > (double)INT_MIN && guess < (double)INT_MAX) {
    const size_t need = alignment_band_width(smax, scoring->gap_extend,
                                             len_i, len_j, (score_t)guess);
    if(need >= len_i) return FILL_SKIPPED;
    // The least any band at least that wide could cost by the rule below. A
    // wider band can get more threads and so cost less than a narrower one
    // (50 kb at 3%), so: on one thread, its columns per row; on several,
    // at least twice the columns per block, and its columns shared by as
    // many threads as the whole fill has
    const double gcols = (double)(D < 0 ? -D : D) + 2.0*need + 1;
    const size_t gstrip = aligner->opts.band_strip_rows
                          ? aligner->opts.band_strip_rows
                          : MAX2(64, (size_t)(gcols / 4));
    const double one = gcols + (double)gstrip;
    const double several =
      MAX2(2.0 * (double)aligner->opts.band_block_min_cols,
           one / MAX2(1.0, full_threads / 2));
    if(MIN2(one, several) >=
       pays * (double)len_i / MAX2(1.0, full_threads / 2)) return FILL_SKIPPED;
  }

  size_t bw = aligner->opts.band_min_width;
  if(have_exact) {
    bw = MAX2(bw, alignment_band_width(smax, scoring->gap_extend,
                                       len_i, len_j, exact));
    atomic_fetch_add(&aln_stat_band_wfa_sized, 1);
  }
  // A lower bound on the best score, once known, narrows strips
  bool have_lb = have_exact;
  score_t lb = have_exact ? exact : 0;
  score_t *top = NULL, *lcol = NULL;
  fill_result_t done = FILL_SKIPPED;
  int widened = 0; // times this band was widened to get more threads
  // Seeds narrow strips further (alignment_seeds.c), built once a lower
  // bound is known
  seeds_t *seeds = NULL;
  bool seeds_tried = false;

  while(done == FILL_SKIPPED)
  {
    // Columns per row, plus the strip height the window grows by
    if(bw >= len_i) break;
    const double cols = (double)(D < 0 ? -D : D) + 2.0*bw + 1;
    size_t strip = aligner->opts.band_strip_rows
                   ? aligner->opts.band_strip_rows
                   : MAX2(64, (size_t)(cols / 4));
    double keep_max = (double)aligner->opts.band_keep_max_cells;
    if(keep_max == 0) {
      static atomic_size_t cached_keep_max = 0;
      size_t km = atomic_load_explicit(&cached_keep_max, memory_order_relaxed);
      if(km == 0) {
        long pages = sysconf(_SC_PHYS_PAGES), page = sysconf(_SC_PAGESIZE);
        km = pages > 0 && page > 0 && (size_t)pages <= SIZE_MAX / (size_t)page
             ? ((size_t)pages * (size_t)page) / 4
             : ((size_t)1 << 30);
        if(km == 0) km = (size_t)1 << 30;
        atomic_store_explicit(&cached_keep_max, km, memory_order_relaxed);
      }
      keep_max = (double)km;
    }
    const bool keep = (double)len_j * (cols + (double)strip) > keep_max;
    if(keep && !aligner->opts.band_strip_rows) {
      size_t most;
      for(most = 64; most * most < 12 * len_j; most++) {}
      strip = MIN2(strip, most);
    }
    // A strip's own rows add as many columns to its band as its diagonals
    // span, so once narrowing leaves few diagonals, strips get shorter
    // (about half the diagonals, at least min_rows): the most strips there
    // can be sizes the kept strips' tables
    const size_t min_rows = MIN2((size_t)32, strip);
    const size_t max_strips = len_j / min_rows + 2;
    // A fill on t threads runs about max(1, t/2) times faster (16 threads:
    // 6-8x); a narrow band runs on one
    const double band_threads =
      (double)alignment_fill_threads(&aligner->opts,
                                     (size_t)cols + strip, strip, false);
    if((cols + strip) / MAX2(1.0, band_threads / 2) >=
       pays * (double)len_i / MAX2(1.0, full_threads / 2)) {
      if(widened < 8 && band_threads < full_threads && bw < len_i) {
        bw += MAX2(bw / 4, (size_t)1);
        widened++;
        continue;
      }
      break;
    }
    widened = 0;

    atomic_fetch_add(&aln_stat_band_tries, 1);
    const ptrdiff_t dlo = MIN2(0, D) - (ptrdiff_t)bw;
    const ptrdiff_t dhi = MAX2(0, D) + (ptrdiff_t)bw;

    // Each row keeps the bytes of its strip's columns
    if(aligner->trace_rows_capacity < h) {
      free(aligner->trace_rows);
      aligner->trace_rows = aln_malloc(h * sizeof(ptrdiff_t));
      aligner->trace_rows_capacity = aligner->trace_rows != NULL ? h : 0;
    }
    if(top == NULL) top = aln_malloc(3 * w * sizeof(score_t));
    score_t *lc = aln_realloc(lcol, 3 * (strip + 1) * sizeof(score_t));
    if(lc != NULL) lcol = lc;
    if(aligner->trace_rows == NULL || top == NULL || lc == NULL) {
      done = FILL_NOMEM;
      break;
    }
    // Room for every strip's bytes, or kept strips: the biggest strip's
    // bytes, and every strip's top row (with column 0, which the first
    // strips use). Strips narrowed by a lower bound (below) need less
    size_t cells = 0, kept = 0;
    for(y0 = 1; y0 <= len_j; y0 = y1 + 1) {
      alignment_strip(y0, strip, len_i, len_j, dlo, dhi, &y1, &c0, &c1);
      if(keep) {
        cells = MAX2(cells, (y1 - y0 + 1) * (c1 - c0));
        kept += 3 * (c1 - c0 + 1);
      }
      else cells += (y1 - y0 + 1) * (c1 - c0);
    }
    trace_strips_free(aligner->trace_strips);
    aligner->trace_strips = NULL;
    trace_strips_t *ts = NULL;
    if(keep) {
      ts = aln_calloc(1, sizeof(trace_strips_t));
      if(ts != NULL) {
        ts->strip = strip;
        ts->nstrips = ts->loaded = max_strips;
        ts->c0 = aln_malloc(max_strips * sizeof(size_t));
        ts->y0s = aln_malloc((max_strips + 1) * sizeof(size_t));
        ts->min = min;
        ts->off = aln_malloc((max_strips + 1) * sizeof(size_t));
        ts->rows = aln_malloc(MAX2(kept, (size_t)3) * sizeof(score_t));
        ts->col0 = aln_malloc(3 * h * sizeof(score_t));
        ts->top = aln_malloc(3 * w * sizeof(score_t));
        ts->lcol = aln_malloc(3 * (strip + 1) * sizeof(score_t));
      }
      if(ts == NULL || ts->off == NULL || ts->c0 == NULL || ts->y0s == NULL ||
         ts->rows == NULL || ts->col0 == NULL || ts->top == NULL ||
         ts->lcol == NULL) {
        trace_strips_free(ts);
        done = FILL_NOMEM;
        break;
      }
      for(m = 0; m < 3; m++)
        memcpy(ts->col0 + m*h, col0[m], h * sizeof(score_t));
    }
    if(!alignment_trace_capacity(aligner, cells)) {
      trace_strips_free(ts);
      done = FILL_NOMEM;
      break;
    }

    // Row 0, then each strip's last row for the next
    score_t *topm[3] = {top, top + w, top + 2*w}, *lcm[3];
    for(m = 0; m < 3; m++) {
      memcpy(topm[m], row0[m], w * sizeof(score_t));
      lcm[m] = lcol + m * (strip + 1);
    }

    size_t jr_prev = len_i; // the top row is right up to here
    size_t s = 0, used = 0, c0_prev = 0, rows_cap = MAX2(kept, (size_t)3);
    ptrdiff_t sdlo = dlo, sdhi = dhi; // this strip's diagonals
    for(y0 = 1; y0 <= len_j; y0 = y1 + 1, s++)
    {
      // Not left of the strip above's left column: a path only moves right,
      // so one reaching a cell there passes left of the strip above too
      sdlo = MAX2(sdlo, (ptrdiff_t)c0_prev + 1 - (ptrdiff_t)y0);
      const size_t span_d = sdhi >= sdlo ? (size_t)(sdhi - sdlo) + 1 : 1;
      const size_t rows = s == 0 ? strip
                          : MIN2(strip, MAX2(min_rows, span_d / 2));
      alignment_strip(y0, rows, len_i, len_j, sdlo, sdhi, &y1, &c0, &c1);
      c0_prev = c0;
      // Columns the strip above did not reach
      for(i = jr_prev + 1; i <= c1; i++)
        topm[MATCH][i] = topm[GAP_A][i] = topm[GAP_B][i] = min;
      jr_prev = c1;
      if(keep) {
        // The strip's top row, and its bytes from the start of trace
        const size_t span = c1 - c0 + 1;
        if(s == 0) ts->off[0] = 0;
        ts->c0[s] = c0;
        ts->y0s[s] = y0;
        if(ts->off[s] + 3 * span > rows_cap) {
          const size_t cap = MAX2(2 * rows_cap, ts->off[s] + 3 * span);
          score_t *grown = aln_realloc(ts->rows, cap * sizeof(score_t));
          if(grown == NULL) break; // out of memory: y0 <= len_j below
          ts->rows = grown;
          rows_cap = cap;
        }
        for(m = 0; m < 3; m++)
          memcpy(ts->rows + ts->off[s] + m * span, topm[m] + c0,
                 span * sizeof(score_t));
        ts->off[s+1] = ts->off[s] + 3 * span;
        for(i = y0; i <= y1; i++)
          aligner->trace_rows[i] = (ptrdiff_t)((i - y0) * (c1 - c0)) -
                                   (ptrdiff_t)(c0 + 1);
      }
      else for(i = y0; i <= y1; i++, used += c1 - c0)
        aligner->trace_rows[i] = (ptrdiff_t)used - (ptrdiff_t)(c0 + 1);

      alignment_block_t tmpl = {
        .aligner = aligner, .table = table, .min = min,
        .c0 = c0, .width = c1 - c0, .y0 = y0, .y1 = y1,
        .r0 = {topm[MATCH], topm[GAP_A], topm[GAP_B]}, .lstride = 1,
        .band = true,
        .trace = aligner->trace, .trace_rows = aligner->trace_rows,
        .lastrow = {topm[MATCH], topm[GAP_A], topm[GAP_B]}};
      for(m = 0; m < 3; m++) {
        if(c0 == 0) tmpl.lc[m] = col0[m]; // the real column 0
        else {
          // Column c0 is outside the band, apart from the cell above the
          // strip (diagonally before its first cell)
          lcm[m][0] = topm[m][c0];
          for(i = 1; i <= y1 - y0 + 1; i++) lcm[m][i] = min;
          tmpl.lc[m] = lcm[m];
          tmpl.ly0 = y0 - 1;
        }
      }
      if(!alignment_fill_blocks(&tmpl)) break;
      // Column c0 of the last row is outside the band, for a strip below
      // with the same left column
      if(c0 > 0)
        topm[MATCH][c0] = topm[GAP_A][c0] = topm[GAP_B][c0] = min;
      sdlo = dlo;
      sdhi = dhi;
      if(have_lb) {
        alignment_band_narrow(topm, col0, y1, c0, c1, len_i, len_j, smax,
                              scoring->gap_extend, lb, &sdlo, &sdhi);
        if(aligner->opts.band_seeds && !seeds_tried) {
          seeds_tried = true;
          seeds = seeds_build(table, scoring, aligner->seq_a, len_i,
                              aligner->seq_b, len_j, smax,
                              aligner->opts.seed_length);
        }
        if(seeds != NULL && y1 < len_j) {
          // Over the rows the next strip would have at these diagonals
          // (it may get fewer: the same cells, in part)
          const size_t next_d = sdhi >= sdlo ? (size_t)(sdhi - sdlo) + 1 : 1;
          const size_t next = MIN2(strip, MAX2(min_rows, next_d / 2));
          seeds_narrow(seeds, topm, col0, y1, c0, c1,
                       MIN2(len_j, y1 + next), smax, scoring->gap_extend,
                       lb, &sdlo, &sdhi);
        }
      }
    }
    if(y0 <= len_j) {
      trace_strips_free(ts);
      done = FILL_NOMEM;
      break;
    }

    score_t best = MAX3(aligner->final_scores[MATCH],
                        aligner->final_scores[GAP_A],
                        aligner->final_scores[GAP_B]);
    size_t need = alignment_band_width(smax, scoring->gap_extend,
                                       len_i, len_j, best);
    if(need <= bw) {
      done = FILL_DONE;
      if(ts != NULL) {
        // The last strip's bytes are in trace, for the start of the walk
        ts->nstrips = s;
        ts->y0s[s] = len_j + 1;
        ts->loaded = s - 1;
        ts->loaded_y = len_j;
        ts->loaded_x = len_i;
        aligner->trace_strips = ts;
      }
    }
    else {
      trace_strips_free(ts);
      bw = MAX2(need, 2*bw);
      have_lb = true; // a real alignment's score
      lb = best;
    }
  }

  seeds_free(seeds);
  free(top);
  free(lcol);
  aligner->trace_banded = done == FILL_DONE;
  return done;
}
#else
void trace_strips_free(void *ts_ptr) { (void)ts_ptr; }
bool alignment_trace_load(aligner_t *aligner, size_t x, size_t y)
{
  (void)aligner; (void)x; (void)y;
  return false;
}
#endif
