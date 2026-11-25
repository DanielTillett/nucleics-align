/*
 hirschberg.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdatomic.h>
#include <unistd.h>

#include "hirschberg_internal.h"

static void hb_emit_del(hb_state_t *st, const char *a, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++, st->pos++) {
    st->out_a[st->pos] = a[i];
    st->out_b[st->pos] = '-';
    st->ops[st->pos] = GAP_B;
  }
}

static void hb_emit_ins(hb_state_t *st, const char *b, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++, st->pos++) {
    st->out_a[st->pos] = '-';
    st->out_b[st->pos] = b[i];
    st->ops[st->pos] = GAP_A;
  }
}

// Column-block threading of a pass: blocks of columns run on separate
// threads, each one row behind the block to its left. Tests lower these
size_t hirschberg_block_min_cells = 4000000; // rows x columns
size_t hirschberg_block_min_cols = 1024; // columns per block
int hirschberg_max_threads = 0; // see hirschberg.h; 0: half the CPUs


// Banded passes (see hb_diff()). Tests lower the last two
bool hirschberg_band_enabled = true;
size_t hirschberg_band_strip_rows = 0; // rows per strip; 0: automatic
size_t hirschberg_band_min_width = 64; // first band to try with no bound
size_t hirschberg_banded_passes = 0;

// Rows reused down the recursion (see hb_diff())
size_t hirschberg_reuse_depth = 2;
size_t hirschberg_reuse_min_rows = 32;
size_t hirschberg_reused_rows = 0;


static void hb_chain_free(hb_chain_t *c)
{
  while(c->len > 0) free(c->h[--c->len]);
}

// Drop the last row of c
static void hb_chain_pop(hb_chain_t *c)
{
  free(c->h[--c->len]);
}

// Rows for a pass of `rows` rows over `width` columns to save: row r for
// the nearest sub-problem below it, then the rows of the sub-problems below
// that one, as long as they have enough rows. A forward pass serves the
// first half of each (r/2 rows next), a reverse pass the second (r - r/2).
// An allocation that fails only saves fewer rows
static void hb_chain_plan(const hb_state_t *st, hb_chain_t *c, size_t rows,
                          size_t r, bool reverse, size_t width)
{
  const size_t depth = MIN2(st->opts.hb_reuse_depth, (size_t)HB_CHAIN_MAX);
  const size_t min_rows = MAX2(st->opts.hb_reuse_min_rows, (size_t)2);
  size_t k;
  c->len = 0;
  // r falls at each step (r - r/2 < r for r >= 2) and stays below rows
  for(; c->len < depth && r >= min_rows && r < rows;
      r = reverse ? r - r/2 : r/2)
  {
    score_t *t = aln_malloc(2 * (width+1) * sizeof(score_t));
    if(t == NULL) break;
    c->h[c->len] = t;
    c->d[c->len] = t + width+1;
    c->rows[c->len++] = r;
  }
  // into rising order
  for(k = 0; k < c->len / 2; k++) {
    size_t o = c->len - 1 - k, r2 = c->rows[k];
    score_t *h = c->h[k], *d = c->d[k];
    c->h[k] = c->h[o]; c->d[k] = c->d[o]; c->rows[k] = c->rows[o];
    c->h[o] = h; c->d[o] = d; c->rows[o] = r2;
  }
}

// Keep columns 0 .. width (in pass order) of the rows of c, for a
// sub-problem that narrow: what the rows held down the recursion add up to
// stays within a few rows of the whole problem. If out of memory, drops them
static void hb_chain_trim(hb_chain_t *c, size_t width)
{
  size_t k;
  for(k = 0; k < c->len; k++) {
    score_t *t = aln_malloc(2 * (width+1) * sizeof(score_t));
    if(t == NULL) { hb_chain_free(c); return; }
    memcpy(t, c->h[k], (width+1) * sizeof(score_t));
    memcpy(t + width+1, c->d[k], (width+1) * sizeof(score_t));
    free(c->h[k]);
    c->h[k] = t;
    c->d[k] = t + width+1;
  }
}

// Copy the last row of c into a pass's rows h, d over columns 0 .. n (the
// reverse rows the other way round), and drop it
static void hb_chain_take(hb_chain_t *c, size_t n, bool reverse,
                          score_t *h, score_t *d)
{
  const size_t k = c->len - 1;
  size_t p;
  for(p = 0; p <= n; p++) {
    size_t j = reverse ? n - p : p;
    h[j] = c->h[k][p];
    d[j] = c->d[k][p];
  }
  hb_chain_pop(c);
  __atomic_add_fetch(&hirschberg_reused_rows, 1, __ATOMIC_RELAXED);
}

size_t hirschberg_striped_mismatches = 0;
size_t hirschberg_part16_max_cols = 0, hirschberg_parts16 = 0;
bool hirschberg_myers_enabled = true;
size_t hirschberg_myers_scores = 0;

// The score by Myers' bit-vector algorithm (alignment_myers_score()) where
// the scoring is unit cost and tests have not turned it off
static bool hb_myers(const hb_state_t *st, const char *a, size_t len_a,
                     const char *b, size_t len_b, score_t *score)
{
  if(!hirschberg_myers_enabled ||
     !alignment_myers_score(&st->table, st->scoring, a, len_a, b, len_b,
                            score)) return false;
  __atomic_add_fetch(&hirschberg_myers_scores, 1, __ATOMIC_RELAXED);
  return true;
}

#ifdef HB_STRIPED
// As hb_myers(), with the banded phase only (declines on divergent pairs,
// which the striped passes then take on several threads)
static bool hb_myers_banded(const hb_state_t *st, const char *a, size_t len_a,
                            const char *b, size_t len_b, score_t *score)
{
  if(!hirschberg_myers_enabled ||
     !alignment_myers_score_banded(&st->table, st->scoring, a, len_a, b, len_b,
                                   score)) return false;
  __atomic_add_fetch(&hirschberg_myers_scores, 1, __ATOMIC_RELAXED);
  return true;
}
#endif

bool hirschberg_simd_available(const alignment_options_t *options)
{
#ifdef HB_STRIPED
  const bool striped = options ? options->striped : alignment_striped_enabled;
  return striped && HB_FAST_VMAX();
#else
  (void)options;
  return false;
#endif
}


score_t hb_midpoint(const hb_state_t *st, size_t n, size_t *midj,
                           int *type)
{
  const score_t g = st->scoring->gap_open;
  const score_t *cc = st->cc, *dd = st->dd, *rr = st->rr, *ss = st->ss;
  score_t midc = cc[0] + rr[0], c;
  size_t j;
  *midj = 0;
  *type = 1;

  for(j = 0; j <= n; j++)
  {
    c = cc[j] + rr[j];
    if(c > midc) { midc = c; *midj = j; *type = 1; }
    c = dd[j] + ss[j] - g; // gap open counted on both sides
    if(c > midc) { midc = c; *midj = j; *type = 2; }
  }
  return midc;
}

// Forward and reverse passes of hb_diff(), into cc, dd (row imid from the
// top) and rr, ss (from the bottom), unless out of memory (st->nomem).
// lower: see hb_diff(). fwd / rev: rows saved by the passes above (see
// hb_diff()), the last one this sub-problem's own row imid from the top /
// bottom if not empty. A pass over all columns takes that row instead of
// making it, or else saves rows for the sub-problems below. Afterwards fwd
// / rev hold the rows for those below
static void hb_passes(hb_state_t *st, const char *a, const char *b,
                      size_t m, size_t n, size_t imid,
                      score_t tb, score_t te, score_t lower,
                      hb_chain_t *fwd, hb_chain_t *rev)
{
  #ifdef HB_STRIPED
  const bool simd = hb_simd_ok(st, m, n);
  if(simd && hb_passes_banded(st, a, b, m, n, imid, tb, te, lower)) {
    // Banded rows are not exact, so not saved; the ones below still are
    if(fwd->len > 0) hb_chain_pop(fwd);
    if(rev->len > 0) hb_chain_pop(rev);
    return;
  }
  #endif
  const bool do_f = fwd->len == 0, do_r = rev->len == 0;
  const size_t rrows = m - imid; // rows of the reverse pass
  if(!do_f) hb_chain_take(fwd, n, false, st->cc, st->dd);
  if(!do_r) hb_chain_take(rev, n, true, st->rr, st->ss);
  // Unit-cost Myers split rows on one thread (faster than unbanded striped
  // SIMD on divergent unit-cost pairs; with gap_open == 0, type 2 in
  // hb_midpoint() requires strict '>' and cannot win over type 1)
  #ifdef HB_STRIPED
  const bool one_thread = st->pass_threads <= 1;
  #else
  const bool one_thread = true;
  #endif
  if(hirschberg_myers_enabled && tb == 0 && te == 0 && one_thread &&
     (!do_f || alignment_myers_row(&st->table, st->scoring, a, imid, b, n,
                                   false, st->cc)) &&
     (!do_r || alignment_myers_row(&st->table, st->scoring, a + imid, rrows,
                                   b, n, true, st->rr))) {
    if(do_f) memcpy(st->dd, st->cc, (n + 1) * sizeof(score_t));
    if(do_r) memcpy(st->ss, st->rr, (n + 1) * sizeof(score_t));
    return;
  }
  if(do_f) hb_chain_plan(st, fwd, imid, imid / 2, false, n);
  if(do_r) hb_chain_plan(st, rev, rrows, rrows - rrows / 2, true, n);
  #ifdef HB_STRIPED
  if(simd) {
    hb_passes_striped(st, a, b, m, n, imid, tb, te, 0, -1, do_f, do_r,
                      fwd, rev);
    if(hirschberg_check_striped && !st->nomem)
      hb_check_striped(st, a, b, m, n, imid, tb, te, false);
    return;
  }
  #endif
  (void)lower;
  hb_passes_scalar(st, a, b, m, n, imid, tb, te, do_f, do_r, fwd, rev);
}

// Align a[0..m) with b[0..n), appending to the output; stops, with the
// output incomplete, if out of memory (st->nomem).
// tb / te: score of opening a deletion that touches the start / end of this
// sub-problem; 0 if it continues a deletion already paid for by the caller.
// lower: a score some alignment of this sub-problem reaches (from the
// caller's split), or HB_NO_BOUND. It narrows the band of the SIMD passes.
// fwd / rev: rows saved for this sub-problem and those below it (taken over
// and freed here). A forward pass's row i depends only on tb, a[0..i) and
// the columns up to its own, so row i of this sub-problem's forward pass,
// cut to the first n_L + 1 columns, is row i of the forward pass of the
// first half (same corner, same tb, a and b starting at the same place).
// Likewise a reverse pass's row serves the second half, which shares the
// bottom-right corner and te. A pass over all columns saves the rows its
// half and the halves of that half will need, if they split at the same
// rows (a deletion across the middle row shifts them by one; such rows go
// unused). Banded rows are lower than exact off the best paths, so only
// unbanded passes save rows
static void hb_diff(hb_state_t *st, const char *a, const char *b,
                    size_t m, size_t n, score_t tb, score_t te,
                    score_t lower, hb_chain_t *fwd, hb_chain_t *rev)
{
  const scoring_t *scoring = st->scoring;
  const score_t g = scoring->gap_open, h = scoring->gap_extend;
  size_t j;

  if(n == 0 || m <= 1) { hb_chain_free(fwd); hb_chain_free(rev); }
  if(n == 0) { hb_emit_del(st, a, m); return; }
  if(m == 0) { hb_emit_ins(st, b, n); return; }

  if(m == 1)
  {
    // Either delete a[0] and insert all of b, or align a[0] with some b[j-1]
    score_t best = MAX2(tb, te) + h + hb_gap(scoring, n), c;
    size_t best_j = 0;

    for(j = 1; j <= n; j++)
    {
      c = hb_gap(scoring, j-1) + hb_sub(scoring, a[0], b[j-1]) +
          hb_gap(scoring, n-j);
      if(c > best) { best = c; best_j = j; }
    }

    if(best_j == 0)
    {
      // Put the deletion next to whichever neighbouring deletion it joins
      if(tb >= te) { hb_emit_del(st, a, 1); hb_emit_ins(st, b, n); }
      else         { hb_emit_ins(st, b, n); hb_emit_del(st, a, 1); }
    }
    else
    {
      hb_emit_ins(st, b, best_j-1);
      st->out_a[st->pos] = a[0];
      st->out_b[st->pos] = b[best_j-1];
      st->ops[st->pos] = MATCH;
      st->pos++;
      hb_emit_ins(st, b+best_j, n-best_j);
    }
    return;
  }

  const score_t *cc = st->cc, *dd = st->dd, *rr = st->rr, *ss = st->ss;
  size_t imid = m / 2, midj;
  int type;

  hb_chain_t none = {.len = 0};

  // Saved rows other than this sub-problem's own are of no use
  if(fwd->len > 0 && fwd->rows[fwd->len-1] != imid) hb_chain_free(fwd);
  if(rev->len > 0 && rev->rows[rev->len-1] != m - imid) hb_chain_free(rev);

  hb_passes(st, a, b, m, n, imid, tb, te, lower, fwd, rev);
  if(st->nomem) { hb_chain_free(fwd); hb_chain_free(rev); return; }
  hb_midpoint(st, n, &midj, &type);
  hb_chain_trim(fwd, midj);
  hb_chain_trim(rev, n - midj);

  // The best score of each half is known now: a bound for its own band
  if(type == 1)
  {
    score_t left = cc[midj], right = rr[midj];
    hb_diff(st, a, b, imid, midj, tb, g, left, fwd, &none);
    if(st->nomem) { hb_chain_free(rev); return; }
    hb_diff(st, a+imid, b+midj, m-imid, n-midj, g, te, right, &none, rev);
  }
  else
  {
    // Without a[imid-1] (opening or extending the deletion) and with the
    // deletion's end free
    score_t left = dd[midj] - g - h, right = ss[midj] - g - h;
    hb_diff(st, a, b, imid-1, midj, tb, 0, left, fwd, &none);
    if(st->nomem) { hb_chain_free(rev); return; }
    hb_emit_del(st, a+imid-1, 2);
    hb_diff(st, a+imid+1, b+midj, m-imid-1, n-midj, 0, te, right, &none, rev);
  }
}

score_t hirschberg_score_alignment(const char *result_a, const char *result_b,
                                   const scoring_t *scoring)
{
  const score_t g = scoring->gap_open, h = scoring->gap_extend;
  enum Matrix prev = MATCH;
  score_t score = 0;
  size_t i;

  for(i = 0; result_a[i] != '\0' && result_b[i] != '\0'; i++)
  {
    if(result_a[i] == '-')
    {
      score += (prev == GAP_A ? h : g + h);
      prev = GAP_A;
    }
    else if(result_b[i] == '-')
    {
      score += (prev == GAP_B ? h : g + h);
      prev = GAP_B;
    }
    else
    {
      score += hb_sub(scoring, result_a[i], result_b[i]);
      prev = MATCH;
    }
  }

  return score;
}

// The score of the alignment written, from the steps taken rather than the
// '-' characters, so sequences containing '-' or '\0' score as
// needleman_wunsch_align() scores them
static score_t hb_score_ops(const hb_state_t *st, const scoring_t *scoring)
{
  const score_t g = scoring->gap_open, h = scoring->gap_extend;
  int prev = MATCH;
  score_t score = 0;
  size_t i;

  for(i = 0; i < st->pos; i++)
  {
    if(st->ops[i] == MATCH) score += hb_sub(scoring, st->out_a[i], st->out_b[i]);
    else score += (st->ops[i] == prev ? h : g + h);
    prev = st->ops[i];
  }

  return score;
}

bool hirschberg_align(const char *a, const char *b,
                      const scoring_t *scoring, alignment_t *result)
{
  return hirschberg_align2(a, b, strlen(a), strlen(b), scoring, result);
}

// Options that need the full matrices. A positive gap_open (opening a gap
// earns points) breaks the Myers-Miller split: its scores then differ from
// Needleman-Wunsch, so it is refused too
static bool hb_unsupported(const scoring_t *scoring)
{
  return scoring->gap_open > 0 ||
         scoring->no_start_gap_penalty || scoring->no_end_gap_penalty ||
         scoring->no_gaps_in_a || scoring->no_gaps_in_b ||
         scoring->no_gaps_in_longest || scoring->no_mismatches;
}

static void hb_state_free(hb_state_t *st);

// Returns false if out of memory, with nothing to free
static bool hb_state_init(hb_state_t *st, const char *a, const char *b,
                          size_t len_a, size_t len_b, const scoring_t *scoring,
                          const alignment_options_t *opts)
{
  memset(st, 0, sizeof(*st));
  st->scoring = scoring;
  alignment_options_resolve(opts, &st->opts);

  score_t *rows = aln_malloc(4 * (len_b+1) * sizeof(score_t));
  st->cc = rows;
  st->b = b;
  st->row_b = aln_malloc((len_b+1) * sizeof(size_t));
  if(rows == NULL || st->row_b == NULL ||
     !sub_table_build(&st->table, scoring, a, len_a, b, len_b)) {
    hb_state_free(st);
    return false;
  }
  size_t j;
  for(j = 0; j < len_b; j++)
    st->row_b[j] = st->table.code[(unsigned char)b[j]] * st->table.n;

  // Every value the passes compute, and the sum of two at the split, must
  // fit in an int. Each is a best score, so at least that of deleting and
  // inserting everything (two gaps), less one more gap open and a gap
  // extension or substitution; and at most the best substitution score per
  // step. Penalties large enough to break this are refused
  #define HB_DABS(x) ((x) < 0 ? -(double)(x) : (double)(x))
  {
    double go = HB_DABS(scoring->gap_open), ge = HB_DABS(scoring->gap_extend);
    double sub = 0, len = (double)len_a + (double)len_b + 2;
    size_t c;
    for(c = 0; c < st->table.n * st->table.n; c++)
      sub = MAX2(sub, HB_DABS(st->table.scores[c]));
    double lo = 3*go + len*ge + sub, hi = len*sub;
    st->in_range = 2*lo + go < (double)INT_MAX && 2*hi < (double)INT_MAX;
  }

#ifdef HB_STRIPED
  // Three striped rows plus a profile row per character code
  size_t seg = (len_b + ALN_VL - 1) / ALN_VL + 1;
  // HB_NEG must stay below every real score, with room to add penalties.
  // No path through the matrices (including padding columns) has more than
  // len_a + len_b + ALN_VL steps, each scoring at most |gap_open| +
  // |gap_extend| or the largest substitution score in magnitude. Beyond a
  // conservative bound use the scalar passes
  {
    double gap = HB_DABS(scoring->gap_open) + HB_DABS(scoring->gap_extend);
    st->smax = st->table.smax;
    double sub = MAX2((double)st->smax, -(double)st->table.smin);
    double steps = (double)len_a + (double)len_b + ALN_VL + 2;
    st->simd_ok = steps * (gap + sub) < (double)(1 << 27);

    // 16-bit blocks (hb_block_pass16()): a cell of a W-column block lies
    // within W (smax + e) + p of the reference (e and p the gap extension
    // and open penalties, negated), the cells above it within one more row
    // (p + e + sub), and a sum taken before the next max adds at most p + e
    // + sub again
    double per_col = (double)st->smax + HB_DABS(scoring->gap_extend);
    double fixed = 4 * gap + 2 * sub;
    st->w16 = 0;
    if(st->opts.narrow && fixed < HB16_RANGE) {
      double w = per_col > 0 ? (HB16_RANGE - fixed) / per_col : 1e18;
      st->w16 = w >= 1e15 ? (size_t)1e15 : (size_t)w / HB_VL16 * HB_VL16;
    }
    // Narrower parts for tests (never below the smallest used)
    if(hirschberg_part16_max_cols > 0 && st->w16 >= 4 * HB_VL16)
      st->w16 = MAX2(MIN2(st->w16, hirschberg_part16_max_cols / HB_VL16 *
                                   HB_VL16), 4 * HB_VL16);
  }

  // Threads per pass: the forward and reverse passes already run at once
  st->pass_threads = st->opts.hb_max_threads;
  if(st->pass_threads <= 0) {
    static atomic_int cached_cpus = 0;
    int cpus = atomic_load_explicit(&cached_cpus, memory_order_relaxed);
    if(cpus <= 0) {
      long n = sysconf(_SC_NPROCESSORS_ONLN);
      cpus = n > 0 ? (int)n : 1;
      atomic_store_explicit(&cached_cpus, cpus, memory_order_relaxed);
    }
    st->pass_threads = cpus > 1 ? cpus / 2 : 1;
  }
  if(st->pass_threads > ALN_MAX_THREADS) st->pass_threads = ALN_MAX_THREADS;

  st->vbuf = NULL;
  st->vbuf2 = NULL;
  st->vbuf_bytes = (3 + st->table.n) * seg * sizeof(vscore_t);
#endif

  st->dd = rows + (len_b+1);
  st->rr = rows + 2*(len_b+1);
  st->ss = rows + 3*(len_b+1);
  return true;
}

static void hb_state_free(hb_state_t *st)
{
  free(st->cc); // all four rows
  free(st->row_b);
#ifdef HB_STRIPED
  free(st->vbuf);
  free(st->vbuf2);
#endif
  sub_table_free(&st->table);
}

bool hirschberg_align2(const char *a, const char *b,
                       size_t len_a, size_t len_b,
                       const scoring_t *scoring, alignment_t *result)
{
  return hirschberg_align_opts(a, b, len_a, len_b, scoring, NULL, result);
}

// The unchecked aligners' result: true for ALN_OK, false for
// ALN_ERR_UNSUPPORTED (the scoring, or penalties too large), and the
// error printed and the program ended for anything else
static bool hb_unchecked(aln_status_t status, const aln_error_t *e)
{
  if(status == ALN_OK) return true;
  if(status == ALN_ERR_UNSUPPORTED) return false;
  alignment_fprint_error(stderr, e);
  exit(EXIT_FAILURE);
}

// hirschberg_align_opts(), returning ALN_ERR_UNSUPPORTED where it returns
// false, or ALN_ERR_NOMEM with no alignment in result
static aln_status_t hb_align(const char *a, const char *b,
                             size_t len_a, size_t len_b,
                             const scoring_t *scoring,
                             const alignment_options_t *opts,
                             alignment_t *result)
{
  if(hb_unsupported(scoring)) return ALN_ERR_UNSUPPORTED;

  hb_state_t st;
  if(!hb_state_init(&st, a, b, len_a, len_b, scoring, opts))
    return ALN_ERR_NOMEM;
  if(!st.in_range) { hb_state_free(&st); return ALN_ERR_UNSUPPORTED; }

  st.ops = aln_malloc(len_a + len_b + 1);
  if(st.ops == NULL || !alignment_reserve(result, len_a + len_b)) {
    free(st.ops);
    hb_state_free(&st);
    return ALN_ERR_NOMEM;
  }
  st.out_a = result->result_a;
  st.out_b = result->result_b;

  hb_chain_t fwd = {.len = 0}, rev = {.len = 0};
  hb_diff(&st, a, b, len_a, len_b, scoring->gap_open, scoring->gap_open,
          HB_NO_BOUND, &fwd, &rev);
  hb_state_free(&st);
  if(st.nomem) {
    free(st.ops);
    result->result_a[0] = result->result_b[0] = '\0';
    result->length = 0;
    return ALN_ERR_NOMEM;
  }

  result->result_a[st.pos] = '\0';
  result->result_b[st.pos] = '\0';
  result->length = st.pos;
  result->score = hb_score_ops(&st, scoring);
  free(st.ops);
  result->pos_a = result->pos_b = 0;
  result->len_a = len_a;
  result->len_b = len_b;
  return ALN_OK;
}

bool hirschberg_align_opts(const char *a, const char *b,
                           size_t len_a, size_t len_b,
                           const scoring_t *scoring,
                           const alignment_options_t *opts,
                           alignment_t *result)
{
  aln_error_t e = {ALN_OK, 0, 0};
  return hb_unchecked(hirschberg_align_checked(a, b, len_a, len_b, scoring,
                                               opts, result, &e), &e);
}

bool hirschberg_score(const char *a, const char *b,
                      const scoring_t *scoring, score_t *score)
{
  return hirschberg_score2(a, b, strlen(a), strlen(b), scoring, score);
}

static aln_status_t hb_score(const char *a, const char *b,
                             size_t len_a, size_t len_b,
                             const scoring_t *scoring,
                             const alignment_options_t *opts, score_t *score);

// The checks the aligner makes before it would stop the program, in its
// order: options it refuses (and sizes it could not hold), then the
// scoring, then the range of scores
static aln_status_t hb_checked(const char *a, const char *b,
                               size_t len_a, size_t len_b,
                               const scoring_t *scoring, aln_error_t *err)
{
  aln_error_t e = {ALN_OK, 0, 0};
  if(hb_unsupported(scoring)) e.status = ALN_ERR_UNSUPPORTED;
  // Rows of len_b + 1 scores, and the alignment's len_a + len_b steps
  else if(len_b >= SIZE_MAX / (4 * sizeof(score_t)) - 1 ||
          len_a >= SIZE_MAX / 2 || len_b >= SIZE_MAX / 2)
    e.status = ALN_ERR_TOO_LONG;
  else alignment_scoring_complete(scoring, a, len_a, b, len_b, &e);
  if(err != NULL) *err = e;
  return e.status;
}

aln_status_t hirschberg_align_checked(const char *a, const char *b,
                                      size_t len_a, size_t len_b,
                                      const scoring_t *scoring,
                                      const alignment_options_t *opts,
                                      alignment_t *result, aln_error_t *err)
{
  aln_status_t status = hb_checked(a, b, len_a, len_b, scoring, err);
  if(status == ALN_OK) {
    // ALN_ERR_UNSUPPORTED: penalties too large
    status = hb_align(a, b, len_a, len_b, scoring, opts, result);
    if(err != NULL) err->status = status;
  }
  return status;
}

aln_status_t hirschberg_score_checked(const char *a, const char *b,
                                      size_t len_a, size_t len_b,
                                      const scoring_t *scoring,
                                      const alignment_options_t *opts,
                                      score_t *score, aln_error_t *err)
{
  aln_status_t status = hb_checked(a, b, len_a, len_b, scoring, err);
  if(status == ALN_OK) {
    status = hb_score(a, b, len_a, len_b, scoring, opts, score);
    if(err != NULL) err->status = status;
  }
  return status;
}

bool hirschberg_score2(const char *a, const char *b,
                       size_t len_a, size_t len_b,
                       const scoring_t *scoring, score_t *score)
{
  return hirschberg_score_opts(a, b, len_a, len_b, scoring, NULL, score);
}

bool hirschberg_score_opts(const char *a, const char *b,
                           size_t len_a, size_t len_b,
                           const scoring_t *scoring,
                           const alignment_options_t *opts, score_t *score)
{
  aln_error_t e = {ALN_OK, 0, 0};
  return hb_unchecked(hirschberg_score_checked(a, b, len_a, len_b, scoring,
                                               opts, score, &e), &e);
}

// hirschberg_score_opts(), returning ALN_ERR_UNSUPPORTED where it returns
// false, or ALN_ERR_NOMEM with score unset
static aln_status_t hb_score(const char *a, const char *b,
                             size_t len_a, size_t len_b,
                             const scoring_t *scoring,
                             const alignment_options_t *opts, score_t *score)
{
  if(hb_unsupported(scoring)) return ALN_ERR_UNSUPPORTED;

  hb_state_t st;
  if(!hb_state_init(&st, a, b, len_a, len_b, scoring, opts))
    return ALN_ERR_NOMEM;
  if(!st.in_range) { hb_state_free(&st); return ALN_ERR_UNSUPPORTED; }
  const score_t g = scoring->gap_open;

  size_t midj;
  int type;

  // Unit-cost scoring: banded Myers beats speculative striped bands at any
  // thread count; on divergent pairs with multiple threads, let multi-threaded
  // striped SIMD handle the full pass instead of single-threaded full Myers
  bool scored;
  #ifdef HB_STRIPED
  const bool simd = hb_simd_ok(&st, len_a, len_b);
  if(st.pass_threads <= 1) scored = hb_myers(&st, a, len_a, b, len_b, score);
  else scored = hb_myers_banded(&st, a, len_a, b, len_b, score);
  #else
  scored = hb_myers(&st, a, len_a, b, len_b, score);
  #endif

  if(scored)
  {
  }
  #ifdef HB_STRIPED
  // Banded forward and reverse passes of hb_diff(), which run at the same
  // time: every alignment crosses the middle row, so the best crossing is
  // the best score
  else if(simd &&
          hb_passes_banded(&st, a, b, len_a, len_b, len_a / 2, g, g, HB_NO_BOUND))
  {
    if(!st.nomem) *score = hb_midpoint(&st, len_b, &midj, &type);
  }
  // Otherwise one forward pass over every row of a: cc[len_b] is the best
  // score, ending in any state
  else if(simd && hb_ensure_vbuf(&st, false))
  {
    hb_pass_striped(&st, a, 1, len_a, b, 1, 0, len_b, g, st.cc, st.dd, 1,
                    st.vbuf, NULL, st.pass_threads);
    if(!st.nomem) *score = st.cc[len_b];
  }
  #endif
  else if(!st.nomem)
  {
    // with imid = len_a the reverse pass has no rows
    hb_passes_scalar(&st, a, b, len_a, len_b, len_a, g, g, true, false,
                     NULL, NULL);
    *score = st.cc[len_b];
  }
  (void)midj; (void)type;

  hb_state_free(&st);
  return st.nomem ? ALN_ERR_NOMEM : ALN_OK;
}
