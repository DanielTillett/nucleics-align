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

// Rows per strip of a band of cols columns per row: the setting, or a
// quarter of the columns and at least 64
static size_t band_strip_rows(const alignment_options_t *o, double cols)
{
  return o->band_strip_rows ? o->band_strip_rows
                            : MAX2(64, (size_t)(cols / 4));
}

// What a band costs against the whole fill: the pair's sizes, the whole
// fill's threads and the fraction of the whole fill's time per cell a band
// must stay under
typedef struct {
  size_t len_i, len_j;
  ptrdiff_t D;         // len_i - len_j, the last cell's diagonal
  double full_threads; // the whole fill's threads
  double pays;         // a band must cost under this much of the whole fill
  score_t smax;
} band_cost_t;

// A fill on t threads runs about max(1, t/2) times faster (16 threads:
// 6-8x); a narrow band runs on one
static double threads_speedup(double threads)
{
  return MAX2(1.0, threads / 2);
}

static void band_cost_init(const aligner_t *aligner, const sub_table_t *table,
                           band_cost_t *c)
{
  c->len_i = aligner->score_width - 1;
  c->len_j = aligner->score_height - 1;
  c->D = (ptrdiff_t)c->len_i - (ptrdiff_t)c->len_j;
  // The whole fill's threads, or one for a narrow band
  c->full_threads =
    (double)alignment_fill_threads(&aligner->opts, c->len_i, c->len_j, true);
  // A band must cost under this much of the whole fill, whose difference-
  // encoded lanes (where they apply) take about 0.6 of the time per cell
  // of the band's absolute ones: banding at 20 kb and 3-20% divergence
  // found the band too narrow and filled the whole matrix anyway, 0.033 ->
  // 0.044 s
  score_t floor, ceiling;
  c->pays = aligner->opts.diff &&
            alignment_diff_bits(table, aligner->scoring, &floor, &ceiling)
            ? 0.35 : 0.6;
  c->smax = table->smax;
}

// Columns per row of a band of half-width bw
static double band_cols(const band_cost_t *c, size_t bw)
{
  return (double)(c->D < 0 ? -c->D : c->D) + 2.0 * (double)bw + 1;
}

// The first band's half-width (*bw): from the best score itself where the
// wavefront algorithm finds it within its budget (close pairs; *have_exact
// and *exact), so the band it needs is filled once (50 kb at 0.1%: 0.040
// -> 0.018 s), else the narrowest allowed. False where no band could pay:
// even the narrowest with the shortest strips against a full compact fill
// on one thread, or one wide enough for a high guess at the best score
// (the first, narrow band would only find that out: 20 kb at 3-20%
// divergence, 0.044 -> 0.033 s)
static bool band_first_width(const aligner_t *aligner, const sub_table_t *table,
                             const band_cost_t *c, size_t *bw,
                             bool *have_exact, score_t *exact)
{
  const scoring_t *scoring = aligner->scoring;
  const size_t len_i = c->len_i, len_j = c->len_j;

  if(c->full_threads <= 1.0) {
    const double min_cols = band_cols(c, aligner->opts.band_min_width);
    const size_t min_strip = band_strip_rows(&aligner->opts, 0);
    if(min_cols + (double)min_strip >= c->pays * (double)len_i) return false;
  }

  double guess;
  const size_t budget = aligner->opts.band_wfa_budget;
  *have_exact = budget > 0 && len_i <= SIZE_MAX / budget &&
    alignment_wfa_score(table, scoring, aligner->seq_a, len_i,
                        aligner->seq_b, len_j, budget * len_i, exact);
  if(*have_exact) guess = *exact;
  if((*have_exact ||
      alignment_score_guess(table, scoring, aligner->seq_a, len_i,
                            aligner->seq_b, len_j, &guess)) &&
     guess > (double)INT_MIN && guess < (double)INT_MAX) {
    const size_t need = alignment_band_width(c->smax, scoring->gap_extend,
                                             len_i, len_j, (score_t)guess);
    if(need >= len_i) return false;
    // The least any band at least that wide could cost by band_plan()'s
    // rule. A wider band can get more threads and so cost less than a
    // narrower one (50 kb at 3%), so: on one thread, its columns per row;
    // on several, at least twice the columns per block, and its columns
    // shared by as many threads as the whole fill has
    const double gcols = band_cols(c, need);
    const size_t gstrip = band_strip_rows(&aligner->opts, gcols);
    const double one = gcols + (double)gstrip;
    const double several =
      MAX2(2.0 * (double)aligner->opts.band_block_min_cols,
           one / threads_speedup(c->full_threads));
    if(MIN2(one, several) >=
       c->pays * (double)len_i / threads_speedup(c->full_threads)) return false;
  }

  *bw = aligner->opts.band_min_width;
  if(*have_exact) {
    *bw = MAX2(*bw, alignment_band_width(c->smax, scoring->gap_extend,
                                         len_i, len_j, *exact));
    atomic_fetch_add(&aln_stat_band_wfa_sized, 1);
  }
  return true;
}

// The strips of one band: their height, whether only one strip's bytes
// are kept at a time (see alignment_band_keep_max_cells), the shortest a
// narrowed strip gets, and the most strips there can be
typedef struct {
  size_t strip, min_rows, max_strips;
  bool keep;
} band_plan_t;

// Cells a band may keep every byte of: the setting, or a quarter of the
// physical memory (found once)
static double band_keep_max(const aligner_t *aligner)
{
  if(aligner->opts.band_keep_max_cells != 0)
    return (double)aligner->opts.band_keep_max_cells;
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
  return (double)km;
}

typedef enum { PLAN_FILL, PLAN_WIDEN, PLAN_STOP } band_verdict_t;

// Plan the strips of a band of half-width *bw: PLAN_FILL with the plan;
// PLAN_WIDEN with *bw raised, where the band would not pay but a wider one
// could get more threads and so cost less (*widened counts those in a row,
// up to 8); PLAN_STOP once no band pays
static band_verdict_t band_plan(const aligner_t *aligner, const band_cost_t *c,
                                size_t *bw, int *widened, band_plan_t *p)
{
  if(*bw >= c->len_i) return PLAN_STOP;
  // Columns per row, plus the strip height the window grows by
  const double cols = band_cols(c, *bw);
  size_t strip = band_strip_rows(&aligner->opts, cols);
  const bool keep = (double)c->len_j * (cols + (double)strip) > band_keep_max(aligner);
  if(keep && !aligner->opts.band_strip_rows) {
    // Kept strips: about sqrt(12 x rows) rows, so a strip's bytes and the
    // rows kept take about the same memory
    size_t most;
    for(most = 64; most * most < 12 * c->len_j; most++) {}
    strip = MIN2(strip, most);
  }
  // A strip's own rows add as many columns to its band as its diagonals
  // span, so once narrowing leaves few diagonals, strips get shorter
  // (about half the diagonals, at least min_rows): the most strips there
  // can be sizes the kept strips' tables
  const size_t min_rows = MIN2((size_t)32, strip);
  const double band_threads =
    (double)alignment_fill_threads(&aligner->opts,
                                   (size_t)cols + strip, strip, false);
  if((cols + strip) / threads_speedup(band_threads) >=
     c->pays * (double)c->len_i / threads_speedup(c->full_threads)) {
    if(*widened < 8 && band_threads < c->full_threads && *bw < c->len_i) {
      *bw += MAX2(*bw / 4, (size_t)1);
      (*widened)++;
      return PLAN_WIDEN;
    }
    return PLAN_STOP;
  }
  *widened = 0;
  p->strip = strip;
  p->keep = keep;
  p->min_rows = min_rows;
  p->max_strips = c->len_j / min_rows + 2;
  return PLAN_FILL;
}

// The buffers of one band fill, beside the aligner's traceback bytes and
// row offsets: the top row of the strip being filled (top, the three
// matrices by column), the column left of a strip (lcol) and, with keep,
// the table of kept strips (ts, with room for rows_cap scores of their
// top rows)
typedef struct {
  score_t *top, *lcol;
  trace_strips_t *ts;
  size_t rows_cap;
} band_bufs_t;

// Room for the band: every strip's bytes, or kept strips: the biggest
// strip's bytes, and every strip's top row (with column 0, which the first
// strips use). Strips narrowed by a lower bound need less. False if out of
// memory (ts freed; top and lcol stay for the caller to free)
static bool band_alloc(aligner_t *aligner, const band_cost_t *c,
                       const band_plan_t *p, score_t min,
                       score_t *const col0[3], ptrdiff_t dlo, ptrdiff_t dhi,
                       band_bufs_t *bufs)
{
  const size_t w = c->len_i + 1, h = c->len_j + 1;
  size_t y0, y1, c0, c1, m;

  // Each row keeps the bytes of its strip's columns
  if(aligner->trace_rows_capacity < h) {
    free(aligner->trace_rows);
    aligner->trace_rows = aln_malloc(h * sizeof(ptrdiff_t));
    aligner->trace_rows_capacity = aligner->trace_rows != NULL ? h : 0;
  }
  if(bufs->top == NULL) bufs->top = aln_malloc(3 * w * sizeof(score_t));
  score_t *lc = aln_realloc(bufs->lcol, 3 * (p->strip + 1) * sizeof(score_t));
  if(lc != NULL) bufs->lcol = lc;
  if(aligner->trace_rows == NULL || bufs->top == NULL || lc == NULL)
    return false;

  size_t cells = 0, kept = 0;
  for(y0 = 1; y0 <= c->len_j; y0 = y1 + 1) {
    alignment_strip(y0, p->strip, c->len_i, c->len_j, dlo, dhi, &y1, &c0, &c1);
    if(p->keep) {
      cells = MAX2(cells, (y1 - y0 + 1) * (c1 - c0));
      kept += 3 * (c1 - c0 + 1);
    }
    else cells += (y1 - y0 + 1) * (c1 - c0);
  }
  trace_strips_free(aligner->trace_strips);
  aligner->trace_strips = NULL;
  bufs->ts = NULL;
  bufs->rows_cap = MAX2(kept, (size_t)3);
  if(p->keep) {
    trace_strips_t *ts = aln_calloc(1, sizeof(trace_strips_t));
    if(ts != NULL) {
      ts->strip = p->strip;
      ts->nstrips = ts->loaded = p->max_strips;
      ts->c0 = aln_malloc(p->max_strips * sizeof(size_t));
      ts->y0s = aln_malloc((p->max_strips + 1) * sizeof(size_t));
      ts->min = min;
      ts->off = aln_malloc((p->max_strips + 1) * sizeof(size_t));
      ts->rows = aln_malloc(bufs->rows_cap * sizeof(score_t));
      ts->col0 = aln_malloc(3 * h * sizeof(score_t));
      ts->top = aln_malloc(3 * w * sizeof(score_t));
      ts->lcol = aln_malloc(3 * (p->strip + 1) * sizeof(score_t));
    }
    if(ts == NULL || ts->off == NULL || ts->c0 == NULL || ts->y0s == NULL ||
       ts->rows == NULL || ts->col0 == NULL || ts->top == NULL ||
       ts->lcol == NULL) {
      trace_strips_free(ts);
      return false;
    }
    for(m = 0; m < 3; m++)
      memcpy(ts->col0 + m*h, col0[m], h * sizeof(score_t));
    bufs->ts = ts;
  }
  if(!alignment_trace_capacity(aligner, cells)) {
    trace_strips_free(bufs->ts);
    bufs->ts = NULL;
    return false;
  }
  return true;
}

// One band fill: what every strip needs, and what the fill learns
typedef struct {
  aligner_t *aligner;
  const sub_table_t *table;
  const band_cost_t *cost;
  score_t min;
  score_t *const *row0, *const *col0; // row 0 and column 0 by matrix
  ptrdiff_t dlo, dhi; // the band's diagonals
  band_plan_t plan;
  band_bufs_t bufs;
  // A lower bound on the best score, once known, narrows strips
  bool have_lb;
  score_t lb;
  // Seeds narrow strips further (alignment_seeds.c), built once a lower
  // bound is known
  seeds_t *seeds;
  bool seeds_tried;
} band_fill_t;

// Fill the band's strips from the top: each over the columns its rows'
// diagonals reach, from the strip above's last row (min where that strip
// did not reach), narrowed by the lower bound and the seeds. Each strip's
// bytes go to the aligner's trace (with keep, over the last strip's; the
// table ts keeps each strip's top row). *strips gets how many there were.
// False if out of memory part way
static bool band_fill_strips(band_fill_t *f, size_t *strips)
{
  aligner_t *aligner = f->aligner;
  const scoring_t *scoring = aligner->scoring;
  const sub_table_t *table = f->table;
  const size_t len_i = f->cost->len_i, len_j = f->cost->len_j;
  const size_t w = len_i + 1;
  const score_t min = f->min, smax = f->cost->smax;
  const size_t strip = f->plan.strip, min_rows = f->plan.min_rows;
  const bool keep = f->plan.keep;
  trace_strips_t *ts = f->bufs.ts;
  size_t i, m, y0, y1, c0, c1;

  // Row 0, then each strip's last row for the next
  score_t *topm[3] = {f->bufs.top, f->bufs.top + w, f->bufs.top + 2*w}, *lcm[3];
  for(m = 0; m < 3; m++) {
    memcpy(topm[m], f->row0[m], w * sizeof(score_t));
    lcm[m] = f->bufs.lcol + m * (strip + 1);
  }

  size_t jr_prev = len_i; // the top row is right up to here
  size_t s = 0, used = 0, c0_prev = 0, rows_cap = f->bufs.rows_cap;
  ptrdiff_t sdlo = f->dlo, sdhi = f->dhi; // this strip's diagonals
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
      if(c0 == 0) tmpl.lc[m] = f->col0[m]; // the real column 0
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
    sdlo = f->dlo;
    sdhi = f->dhi;
    if(f->have_lb) {
      alignment_band_narrow(topm, f->col0, y1, c0, c1, len_i, len_j, smax,
                            scoring->gap_extend, f->lb, &sdlo, &sdhi);
      if(aligner->opts.band_seeds && !f->seeds_tried) {
        f->seeds_tried = true;
        f->seeds = seeds_build(table, scoring, aligner->seq_a, len_i,
                               aligner->seq_b, len_j, smax,
                               aligner->opts.seed_length);
      }
      if(f->seeds != NULL && y1 < len_j) {
        // Over the rows the next strip would have at these diagonals
        // (it may get fewer: the same cells, in part)
        const size_t next_d = sdhi >= sdlo ? (size_t)(sdhi - sdlo) + 1 : 1;
        const size_t next = MIN2(strip, MAX2(min_rows, next_d / 2));
        seeds_narrow(f->seeds, topm, f->col0, y1, c0, c1,
                     MIN2(len_j, y1 + next), smax, scoring->gap_extend,
                     f->lb, &sdlo, &sdhi);
      }
    }
  }
  *strips = s;
  return y0 > len_j;
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
// strips of at most about sqrt(12 x rows) rows.
// In turn: band_cost_init() and band_first_width() decide whether any band
// could pay and how wide the first is; then for each width band_plan() its
// strips, band_alloc() their room and band_fill_strips() the fill, and the
// best score found says whether the band held (done) or how much wider the
// next must be
fill_result_t alignment_fill_banded(aligner_t *aligner,
                                           const sub_table_t *table,
                                           score_t min, score_t *const row0[3],
                                           score_t *const col0[3])
{
  const scoring_t *scoring = aligner->scoring;
  if(!aligner->opts.band || scoring->gap_open > 0 ||
     scoring->gap_extend > 0) return FILL_SKIPPED;

  band_cost_t cost;
  band_cost_init(aligner, table, &cost);
  size_t bw;
  bool have_exact;
  score_t exact;
  if(!band_first_width(aligner, table, &cost, &bw, &have_exact, &exact))
    return FILL_SKIPPED;

  band_fill_t f = {aligner, table, &cost, min, row0, col0, 0, 0,
                   {0, 0, 0, false}, {NULL, NULL, NULL, 0},
                   have_exact, have_exact ? exact : 0, NULL, false};
  fill_result_t done = FILL_SKIPPED;
  int widened = 0; // times this band was widened to get more threads

  while(done == FILL_SKIPPED)
  {
    const band_verdict_t verdict = band_plan(aligner, &cost, &bw, &widened,
                                             &f.plan);
    if(verdict == PLAN_WIDEN) continue;
    if(verdict == PLAN_STOP) break;

    atomic_fetch_add(&aln_stat_band_tries, 1);
    f.dlo = MIN2(0, cost.D) - (ptrdiff_t)bw;
    f.dhi = MAX2(0, cost.D) + (ptrdiff_t)bw;
    if(!band_alloc(aligner, &cost, &f.plan, min, col0, f.dlo, f.dhi,
                   &f.bufs)) {
      done = FILL_NOMEM;
      break;
    }
    size_t strips;
    if(!band_fill_strips(&f, &strips)) {
      trace_strips_free(f.bufs.ts);
      f.bufs.ts = NULL;
      done = FILL_NOMEM;
      break;
    }

    score_t best = MAX3(aligner->final_scores[MATCH],
                        aligner->final_scores[GAP_A],
                        aligner->final_scores[GAP_B]);
    size_t need = alignment_band_width(cost.smax, scoring->gap_extend,
                                       cost.len_i, cost.len_j, best);
    if(need <= bw) {
      done = FILL_DONE;
      if(f.bufs.ts != NULL) {
        // The last strip's bytes are in trace, for the start of the walk
        trace_strips_t *ts = f.bufs.ts;
        ts->nstrips = strips;
        ts->y0s[strips] = cost.len_j + 1;
        ts->loaded = strips - 1;
        ts->loaded_y = cost.len_j;
        ts->loaded_x = cost.len_i;
        aligner->trace_strips = ts;
        f.bufs.ts = NULL;
      }
    }
    else {
      trace_strips_free(f.bufs.ts);
      f.bufs.ts = NULL;
      bw = MAX2(need, 2*bw);
      f.have_lb = true; // a real alignment's score
      f.lb = best;
    }
  }

  seeds_free(f.seeds);
  free(f.bufs.top);
  free(f.bufs.lcol);
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
