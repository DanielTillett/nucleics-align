/*
 hirschberg.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// Linear space global alignment with affine gaps, following the 'diff'
// procedure of Myers & Miller (1988) Comput Appl Biosci 4:11-17, written to
// maximise a score rather than minimise a cost.
// A gap of length k scores gap_open + k*gap_extend, as in needleman_wunsch.c.
// Deleting a character of a means a gap in b, inserting one of b a gap in a.

#include <stdlib.h>
#include <stddef.h> // ptrdiff_t
#include <stdio.h>
#include <string.h>
#include <limits.h> // INT_MAX
#include <pthread.h>
#include <sched.h> // sched_yield
#include <stdatomic.h>
#include <unistd.h> // sysconf

#include "alignment_macros.h"
#include "alignment_simd.h"
#include "hirschberg.h"
#include "alignment_testing.h"
#include "alignment_alloc.h"

// SIMD passes only where the CPU has a native vector max: without one
// (plain SSE2, the x86-64 baseline) the scalar loops are faster here.
// With runtime dispatch (ALN_DISPATCH) that is checked when running
#if defined(ALN_STRIPED) && defined(ALN_NATIVE_VMAX)
#  define HB_STRIPED 1
#  define HB_FAST_VMAX() 1
#elif defined(ALN_STRIPED) && defined(ALN_DISPATCH)
#  define HB_STRIPED 1
#  define HB_FAST_VMAX() __builtin_cpu_supports("sse4.1")
#endif

typedef struct
{
  const scoring_t *scoring;
  sub_table_t table;
  const char *b; // whole of seq b, to index row_b
  size_t *row_b; // table row offset for each character of b
  score_t *cc, *dd, *rr, *ss; // forward and reverse rows, len_b+1 each
  char *out_a, *out_b;
  size_t pos; // next position to write in out_a, out_b
  uint8_t *ops; // step at each position (MATCH, GAP_A, GAP_B), for the score
  bool in_range; // no score or sum of two can overflow, see hb_state_init()
  bool nomem; // an allocation failed: the rows are not usable (hb_nomem())
  alignment_options_t opts; // settings in effect
#ifdef HB_STRIPED
  void *vbuf; // striped rows and score profile, see hb_pass_striped()
  void *vbuf2; // the same for a pass on a second thread
  size_t vbuf_bytes; // lazy allocation size for vbuf / vbuf2
  bool simd_ok; // all scores far enough from HB_NEG and overflow
  int pass_threads; // threads for the passes of one step (see hb_passes_striped())
  score_t smax; // best substitution score, at least 0, see hb_band()
  size_t w16; // most columns (lanes x segments) of a 16-bit block; 0: none
#endif
} hb_state_t;

static inline score_t hb_gap(const scoring_t *scoring, size_t k)
{
  return k == 0 ? 0 : scoring->gap_open + (score_t)k * scoring->gap_extend;
}

static inline score_t hb_sub(const scoring_t *scoring, char a, char b)
{
  int score;
  bool is_match;
  scoring_lookup(scoring, a, b, &score, &is_match);
  return score;
}

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

// No score known for a sub-problem yet, see hb_diff()
#define HB_NO_BOUND SCORE_MIN

// An allocation failed, perhaps on one of two threads running passes at
// once; read after they are joined
static inline void hb_nomem(hb_state_t *st)
{
  __atomic_store_n(&st->nomem, true, __ATOMIC_RELAXED);
}

// Banded passes (see hb_diff()). Tests lower the last two
bool hirschberg_band_enabled = true;
size_t hirschberg_band_strip_rows = 0; // rows per strip; 0: automatic
size_t hirschberg_band_min_width = 64; // first band to try with no bound
size_t hirschberg_banded_passes = 0;

// Rows reused down the recursion (see hb_diff())
size_t hirschberg_reuse_depth = 2;
size_t hirschberg_reuse_min_rows = 32;
size_t hirschberg_reused_rows = 0;

// Rows a pass saved for sub-problems below it (see hb_diff()): for k <
// len, row rows[k] of a forward pass from the sub-problem's top-left corner,
// or of a reverse pass from its bottom-right corner, as h[k][p] (best score)
// and d[k][p] (ending in a deletion) for p = 0 .. the width in pass order:
// column p of a forward pass, column n - p of a reverse one. rows[] rise
// with k; the last is the row the nearest sub-problem needs
#define HB_CHAIN_MAX 8
typedef struct {
  score_t *h[HB_CHAIN_MAX], *d[HB_CHAIN_MAX];
  size_t rows[HB_CHAIN_MAX], len;
} hb_chain_t;

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

#ifdef HB_STRIPED

// Below any real score, with room to add penalties without overflow
#define HB_NEG (SCORE_MIN / 4)

// Sub-problems with at least this many cells run their two passes on two
// threads; below it, starting a thread costs more than it saves
#define HB_THREAD_CELLS 1e6


// What a block passes to the block on its right, for each row r = 0..rows:
// H[r] best score in its last column, F[r] the insertion gap score entering
// the next column. done counts the rows published so far (64-byte aligned
// to avoid false sharing between adjacent blocks)
typedef struct {
  _Alignas(64) score_t *H, *F;
  atomic_size_t done;
  atomic_uint waiters;
  pthread_mutex_t mu;
  pthread_cond_t cv;
} hb_boundary_t;

#define HB_PUB_ROWS 16

static inline void hb_publish_row(hb_boundary_t *right, size_t val)
{
  atomic_store_explicit(&right->done, val, memory_order_release);
  atomic_thread_fence(memory_order_seq_cst);
  if(atomic_load_explicit(&right->waiters, memory_order_relaxed) > 0) {
    pthread_mutex_lock(&right->mu);
    pthread_cond_broadcast(&right->cv);
    pthread_mutex_unlock(&right->mu);
  }
}

// Wait until the block on the left has published row r
static inline void hb_wait_row(hb_boundary_t *left, size_t r, size_t *seen)
{
  if(*seen > r) return;
  unsigned spins = 0;
  size_t cur;
  while((cur = atomic_load_explicit(&left->done, memory_order_acquire)) <= r) {
    if(spins < 512) {
#if defined(__x86_64__) || defined(__i386__)
      __builtin_ia32_pause();
#endif
      spins++;
      continue;
    }
    pthread_mutex_lock(&left->mu);
    atomic_fetch_add_explicit(&left->waiters, 1, memory_order_seq_cst);
    atomic_thread_fence(memory_order_seq_cst);
    while((cur = atomic_load_explicit(&left->done, memory_order_acquire)) <= r) {
      pthread_cond_wait(&left->cv, &left->mu);
    }
    atomic_fetch_sub_explicit(&left->waiters, 1, memory_order_relaxed);
    pthread_mutex_unlock(&left->mu);
    break;
  }
  *seen = cur;
}

// Where a strip of a banded pass starts (see hb_pass_banded()): its top row
// and the column on its left. A whole pass (NULL) starts at row 0
typedef struct {
  const score_t *H, *D; // top row, best score and ending in a deletion,
                        // indexed by column of the pass
  score_t col0; // running deletion score in the left column, for column 0
  bool neg_left; // left column outside the band: HB_NEG in every row
  bool real; // every cell of the top row is a real score (no HB_NEG)
} hb_top_t;

// 16-bit lanes, twice as many per vector, for hb_block_pass16()
#define HB_VL16 (2 * ALN_VL)
typedef int16_t v16_t __attribute__((vector_size(HB_VL16 * sizeof(int16_t))));

// Stands for -infinity in 16-bit lanes: below every real score there (see
// hb_state_init(): those stay within HB16_RANGE of 0), with room for a
// penalty to be added
#define HB16_NEG (-16000)
#define HB16_RANGE 12000
#define HB16_MAX_PARTS 1024

ALN_INLINE v16_t v16splat(int16_t x) { return (v16_t){0} + x; }

ALN_INLINE v16_t v16max(v16_t a, v16_t b)
{
  v16_t gt = a > b;
  return (a & gt) | (b & ~gt);
}

ALN_INLINE v16_t v16shift(v16_t v, int16_t x)
{
#if HB_VL16 == 32
  return __builtin_shufflevector(v, v16splat(x), 32, 0, 1, 2, 3, 4, 5, 6, 7,
                                 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30);
#elif HB_VL16 == 16
  return __builtin_shufflevector(v, v16splat(x), 16, 0, 1, 2, 3, 4, 5, 6, 7,
                                 8, 9, 10, 11, 12, 13, 14);
#else
  return __builtin_shufflevector(v, v16splat(x), 8, 0, 1, 2, 3, 4, 5, 6);
#endif
}

ALN_INLINE bool v16any(v16_t v)
{
  int16_t r = 0;
  int l;
  for(l = 0; l < HB_VL16; l++) r |= v[l];
  return r != 0;
}

// One part of a 16-bit block (see hb_block_pass16()): columns c0 ..
// c0+width-1 of the pass, its rows and profile, and its reference, the
// score of the column on its left in the row last filled
typedef struct {
  size_t c0, width, seg;
  v16_t *H, *Hn, *D, *prof;
  score_t ref;
} hb16_part_t;

// The next row of part s, given its left column's score in that row (col0)
// and the insertion gap entering it (fstart), from the row of the sequence
// whose character has code `code`. Gives the same for the part on the right:
// the score of the last column (*hlast) and the gap leaving it (*fout)
ALN_INLINE void hb16_row(hb16_part_t *s, size_t code, score_t col0,
                         score_t fstart, score_t g, score_t h,
                         score_t *hlast, score_t *fout)
{
  const size_t seg = s->seg;
  const v16_t vg = v16splat((int16_t)g), vh = v16splat((int16_t)h);
  const v16_t vneg = v16splat(HB16_NEG);
  const v16_t *pr = s->prof + code * seg;
  v16_t *H = s->H, *Hn = s->Hn, *D = s->D, vdiag, vf, vd, vc;
  size_t k, l;

  // The new reference less the old: off the row above
  const v16_t vdelta = v16splat((int16_t)(col0 - s->ref)), vhd = vh - vdelta;
  const long long frel = MAX2((long long)fstart - col0, (long long)HB16_NEG);

  vdiag = v16shift(H[seg-1], 0) - vdelta; // the old reference less itself
  vf = v16shift(vneg, (int16_t)frel);
  for(k = 0; k < seg; k++)
  {
    vd = v16max(D[k], H[k] + vg) + vhd;
    vc = v16max(v16max(vd, vf), vdiag + pr[k]);
    vdiag = H[k] - vdelta;
    D[k] = vd;
    Hn[k] = vc;
    vf = v16max(vf, vc + vg) + vh;
  }

  long long fo = vf[HB_VL16-1];

  // Lazy F with the carries found first (see hb_block_pass())
  vf = v16shift(vf, HB16_NEG);
  for(l = 2; l < HB_VL16; l++) {
    long long carry = (long long)vf[l-1] + (long long)seg * h;
    if(carry > vf[l]) vf[l] = (int16_t)carry;
  }
  fo = MAX2(fo, (long long)vf[HB_VL16-1] + (long long)seg * h);
  k = 0;
  while(v16any(vf > Hn[k] + vg)) {
    Hn[k] = v16max(Hn[k], vf);
    if(++k == seg) break;
    vf = v16max(vf + vh, vneg);
  }

  *hlast = Hn[seg-1][HB_VL16-1] + col0;
  *fout = fo <= HB16_NEG ? HB_NEG : (score_t)(fo + col0);
  s->H = Hn;
  s->Hn = H;
  s->ref = col0;
}

// hb_block_pass() with 16-bit lanes: the same values, for a block whose top
// row (if any) is all real scores. The block is cut into parts of at most
// st->w16 columns (a multiple of HB_VL16 but the last), filled one row at a
// time, left to right; each part hands the next its last column's score
// and the insertion gap leaving it, as blocks do. A part keeps each row
// less one 32-bit score, the row's cell in the column on its left, and
// every cell of the row lies within HB16_RANGE of it: a cell is at least
// that cell plus a gap across the part, and at most the score of the best
// path from the left column down to it, which crosses the part's columns
// with at most one gap open, a substitution score per column and a gap
// extension per row it moves down more than across, while the left column
// loses one per row (see hb_state_init()). Moving to the next row's
// reference takes the difference off the cells above, in the gap extension
// for the deletion and once more for the diagonal
ALN_TARGET_CLONES
static void hb_block_pass16(hb_state_t *st, const char *a, ptrdiff_t astep,
                            size_t rows, const char *b, ptrdiff_t bstep,
                            size_t c0, size_t width, score_t t0,
                            score_t *oc, score_t *od, ptrdiff_t ostep,
                            void *buf, hb_boundary_t *left,
                            hb_boundary_t *right, const hb_top_t *top)
{
  const scoring_t *scoring = st->scoring;
  const sub_table_t *t = &st->table;
  const score_t g = scoring->gap_open, h = scoring->gap_extend;
  const size_t nparts = (width + st->w16 - 1) / st->w16;
  const size_t pw = ((width + nparts - 1) / nparts + HB_VL16 - 1) /
                    HB_VL16 * HB_VL16;
  hb16_part_t parts[nparts];
  v16_t *next = buf;
  size_t r, s, k, l, p, c;
  if(nparts > 1) __atomic_add_fetch(&hirschberg_parts16, 1, __ATOMIC_RELAXED);

  // Each part's rows and profile, then its top row: row 0, or the strip's
  // top row, less the reference (padding lanes copy the last column, so
  // they stay in range too)
  score_t col0 = t0, fstart, hl, fo;
  for(s = 0; s < nparts; s++) {
    hb16_part_t *q = &parts[s];
    q->c0 = c0 + s * pw;
    q->width = s + 1 < nparts ? pw : width - s * pw;
    q->seg = (q->width + HB_VL16 - 1) / HB_VL16;
    const size_t seg = q->seg;
    q->H = next; q->Hn = q->H + seg; q->D = q->Hn + seg; q->prof = q->D + seg;
    next = q->prof + t->n * seg;

    memset(q->prof, 0, t->n * seg * sizeof(v16_t));
    for(p = 0; p < q->width; p++) {
      size_t row = t->code[(unsigned char)b[(ptrdiff_t)(q->c0+p) * bstep]] *
                   t->n;
      l = p / seg; k = p - l*seg;
      for(c = 0; c < t->n; c++)
        q->prof[c*seg + k][l] = (int16_t)t->scores[row + c];
    }

    if(s == 0 && left != NULL) {
      size_t seen0 = 0;
      hb_wait_row(left, 0, &seen0);
      q->ref = left->H[0];
    }
    else if(top != NULL) q->ref = top->H[q->c0];
    else q->ref = q->c0 == 0 ? 0 : g + (score_t)q->c0 * h;
    if(s == 0 && left == NULL && top != NULL) col0 = top->col0;

    for(k = 0; k < seg; k++)
      for(l = 0; l < HB_VL16; l++) {
        p = q->c0 + MIN2(l*seg + k, q->width - 1) + 1;
        score_t hv = top == NULL ? g + (score_t)p * h : top->H[p];
        score_t dv = top == NULL ? hv + g : top->D[p];
        q->H[k][l] = (int16_t)(hv - q->ref);
        q->D[k][l] = (int16_t)(dv - q->ref);
      }
  }

  if(right != NULL) {
    right->H[0] = top == NULL ? g + (score_t)(c0 + width) * h
                              : top->H[c0+width];
    right->F[0] = HB_NEG; // not used
    hb_publish_row(right, 1);
  }

  size_t left_seen = 0;
  for(r = 0; r < rows; r++)
  {
    const size_t code = t->code[(unsigned char)a[(ptrdiff_t)r * astep]];

    if(left == NULL) {
      col0 += h;
      fstart = col0 + g + h;
    } else {
      hb_wait_row(left, r+1, &left_seen);
      col0 = left->H[r+1];
      fstart = left->F[r+1];
    }
    hl = col0; fo = fstart;
    for(s = 0; s < nparts; s++)
      hb16_row(&parts[s], code, hl, fo, g, h, &hl, &fo);

    if(right != NULL) {
      right->H[r+1] = hl;
      right->F[r+1] = fo;
      if(((r + 1) & (HB_PUB_ROWS - 1)) == 0 || r + 1 == rows)
        hb_publish_row(right, r+2);
    }
  }

  if(left == NULL) oc[(ptrdiff_t)c0 * ostep] = od[(ptrdiff_t)c0 * ostep] = col0;
  for(s = 0; s < nparts; s++) {
    const hb16_part_t *q = &parts[s];
    for(l = 0; l < HB_VL16; l++)
      for(k = 0, p = l*q->seg; k < q->seg && p < q->width; k++, p++) {
        oc[(ptrdiff_t)(q->c0+p+1) * ostep] = q->H[k][l] + q->ref;
        od[(ptrdiff_t)(q->c0+p+1) * ostep] = q->D[k][l] + q->ref;
      }
  }
}

// Striped SIMD (Farrar 2007) version of one pass of hb_diff() over the
// columns c0 .. c0+width-1 of the pass: gives exactly the same values as the
// scalar loops. Rows run over `rows` characters of a starting at a[0] and
// moving by astep (+1 forward, -1 reverse); columns run over characters of b
// the same way (bstep), b pointing at column 0 of the pass. t0 is tb or te.
// The last row is written to oc[p*ostep] (best score) and od[p*ostep]
// (ending in a deletion) for column p in pass order; the first block also
// writes column p = 0.
// The insertion (horizontal) gap feeds back into the best score of the same
// row, so Farrar's lazy F loop carries it across segment boundaries (with
// the carries found first, see the loop). That is
// exact only for gap_open <= 0 and gap_extend <= 0, which the caller checks.
// left / right are the boundaries shared with neighbouring blocks (NULL at
// the ends). A block with a right neighbour must have width a multiple of
// ALN_VL, so its last column is the last cell of its last lane.
// top, if not NULL, gives the top row and left column of a strip instead
ALN_TARGET_CLONES
static void hb_block_pass(hb_state_t *st, const char *a, ptrdiff_t astep,
                          size_t rows, const char *b, ptrdiff_t bstep,
                          size_t c0, size_t width, score_t t0,
                          score_t *oc, score_t *od, ptrdiff_t ostep,
                          void *buf, hb_boundary_t *left, hb_boundary_t *right,
                          const hb_top_t *top)
{
  // Twice the lanes where the scores fit (see hb_block_pass16())
  // (parts of at least 4 vectors' columns keep its rows and profiles within
  // the space this block's 32-bit ones take; at most HB16_MAX_PARTS parts,
  // as their list is on the stack)
  if((top == NULL || top->real) && st->w16 >= 4 * HB_VL16 &&
     (width + st->w16 - 1) / st->w16 <= HB16_MAX_PARTS &&
     (right == NULL || width % HB_VL16 == 0)) {
    hb_block_pass16(st, a, astep, rows, b, bstep, c0, width, t0, oc, od,
                    ostep, buf, left, right, top);
    return;
  }

  const scoring_t *scoring = st->scoring;
  const sub_table_t *t = &st->table;
  const score_t g = scoring->gap_open, h = scoring->gap_extend;
  const vscore_t vg = vsplat(g), vh = vsplat(h), vneg = vsplat(HB_NEG);
  const size_t seg = (width + ALN_VL - 1) / ALN_VL;
  vscore_t *H = buf, *Hn = H + seg, *D = Hn + seg, *prof = D + seg;
  vscore_t *tmp, vdiag, vf, vd, vc;
  size_t r, k, l, p, c;

  // prof[c*seg + k] lane l: score of a character with code c against the
  // b character in column c0+l*seg+k (0 for padding past the block)
  memset(prof, 0, t->n * seg * sizeof(vscore_t));
  for(p = 0; p < width; p++) {
    size_t row = t->code[(unsigned char)b[(ptrdiff_t)(c0+p) * bstep]] * t->n;
    l = p / seg; k = p - l*seg;
    for(c = 0; c < t->n; c++) prof[c*seg + k][l] = t->scores[row + c];
  }

  // Row 0: only insertions; or the strip's top row (padding: HB_NEG)
  for(k = 0; k < seg; k++)
    for(l = 0; l < ALN_VL; l++) {
      p = l*seg + k;
      if(top == NULL) {
        H[k][l] = g + (score_t)(c0 + p + 1) * h;
        D[k][l] = H[k][l] + g;
      } else {
        H[k][l] = p < width ? top->H[c0+p+1] : HB_NEG;
        D[k][l] = p < width ? top->D[c0+p+1] : HB_NEG;
      }
    }
  if(right != NULL) {
    right->H[0] = top == NULL ? g + (score_t)(c0 + width) * h
                              : top->H[c0+width];
    right->F[0] = HB_NEG; // not used
    hb_publish_row(right, 1);
  }

  // The column left of the block: column 0 (only deletions) for the first
  // block, otherwise the last column of the block on the left
  score_t col0 = t0, col0_prev = 0, fstart;
  const bool neg_left = top != NULL && top->neg_left;
  // (with a left block, its last column comes from the boundary instead:
  // the left block may be writing top->H[c0] as the strip's output)
  if(top != NULL && left == NULL) { col0 = top->col0; col0_prev = top->H[c0]; }

  size_t left_seen = 0;
  for(r = 0; r < rows; r++)
  {
    const vscore_t *pr = prof +
      (size_t)t->code[(unsigned char)a[(ptrdiff_t)r * astep]] * seg;

    if(neg_left && left == NULL) {
      col0 = fstart = HB_NEG;
    } else if(left == NULL) {
      col0 += h;
      fstart = col0 + g + h;
    } else {
      hb_wait_row(left, r+1, &left_seen);
      col0_prev = left->H[r];
      col0 = left->H[r+1];
      fstart = left->F[r+1];
    }

    vdiag = vshift(H[seg-1], col0_prev);
    vf = vshift(vneg, fstart); // lanes after 0 are fixed below
    for(k = 0; k < seg; k++)
    {
      vd = vmax(D[k], H[k] + vg) + vh;
      vc = vmax(vmax(vd, vf), vdiag + pr[k]);
      vdiag = H[k];
      D[k] = vd;
      Hn[k] = vc;
      vf = vmax(vf, vc + vg) + vh;
    }

    // Insertion gap leaving the last lane: into the next block
    score_t fout = vf[ALN_VL-1];

    // Lazy F with the carries found first: what lane l passes to lane l+1
    // is the better of its own and what it was passed, less seg gap
    // extensions (a scan over the lanes, as in alignment_diff.inc). Each
    // lane then takes its carry along its cells while it could still raise
    // a best score, once, rather than the carries going round the lanes
    // until they stop (up to ALN_VL-1 times). A carry stopped early can
    // never win later: the gap leaving that cell already beats it
    vf = vshift(vf, HB_NEG);
    for(l = 2; l < ALN_VL; l++) {
      long long carry = (long long)vf[l-1] + (long long)seg * h;
      if(carry > vf[l]) vf[l] = (score_t)carry;
    }
    fout = (score_t)MAX2((long long)fout,
                         (long long)vf[ALN_VL-1] + (long long)seg * h);
    k = 0;
    while(vany(vf > Hn[k] + vg)) {
      Hn[k] = vmax(Hn[k], vf);
      if(++k == seg) break;
      vf = vmax(vf + vh, vneg);
    }

    if(right != NULL) {
      right->H[r+1] = Hn[seg-1][ALN_VL-1];
      right->F[r+1] = fout;
      if(((r + 1) & (HB_PUB_ROWS - 1)) == 0 || r + 1 == rows)
        hb_publish_row(right, r+2);
    }

    tmp = H; H = Hn; Hn = tmp;
    col0_prev = col0;
  }

  if(left == NULL) oc[(ptrdiff_t)c0 * ostep] = od[(ptrdiff_t)c0 * ostep] = col0;
  for(l = 0; l < ALN_VL; l++)
    for(k = 0, p = l*seg; k < seg && p < width; k++, p++) {
      oc[(ptrdiff_t)(c0+p+1) * ostep] = H[k][l];
      od[(ptrdiff_t)(c0+p+1) * ostep] = D[k][l];
    }
}

// One column block of a pass, run on its own thread
typedef struct {
  hb_state_t *st; const char *a; ptrdiff_t astep; size_t rows;
  const char *b; ptrdiff_t bstep; size_t c0, width; score_t t0;
  score_t *oc, *od; ptrdiff_t ostep; void *buf;
  hb_boundary_t *left, *right; const hb_top_t *top;
} hb_block_args_t;

static void* hb_block_thread(void *ptr)
{
  hb_block_args_t *p = ptr;
  hb_block_pass(p->st, p->a, p->astep, p->rows, p->b, p->bstep, p->c0,
                p->width, p->t0, p->oc, p->od, p->ostep, p->buf,
                p->left, p->right, p->top);
  return NULL;
}

// One pass of hb_diff() over the n columns c0+1 .. c0+n: as one block, or
// split into column blocks on several threads when it is big enough (see
// hb_block_pass()). A whole pass has c0 = 0 and top = NULL
static void hb_pass_striped(hb_state_t *st, const char *a, ptrdiff_t astep,
                            size_t rows, const char *b, ptrdiff_t bstep,
                            size_t c0, size_t n, score_t t0,
                            score_t *oc, score_t *od, ptrdiff_t ostep,
                            void *buf, const hb_top_t *top, int threads)
{
  // At least one vector of columns per block, so no block is empty
  const size_t min_cols = MAX2(st->opts.hb_block_min_cols, (size_t)ALN_VL);
  size_t nblocks = 1, i;
  if(threads > 1 && (double)rows * n >= st->opts.hb_block_min_cells &&
     n >= 2 * min_cols)
  {
    nblocks = MIN2((size_t)threads, n / min_cols);
  }

  if(nblocks <= 1) {
    hb_block_pass(st, a, astep, rows, b, bstep, c0, n, t0, oc, od, ostep, buf,
                  NULL, NULL, top);
    return;
  }

  // Blocks of equal width (a multiple of ALN_VL, or of HB_VL16 so they can
  // take 16-bit lanes); the last takes the rest
  const size_t unit = n / nblocks >= HB_VL16 ? HB_VL16 : ALN_VL;
  size_t width = (n / nblocks) / unit * unit;
  size_t seg = (n - width * (nblocks-1) + ALN_VL - 1) / ALN_VL; // widest
  size_t bufsize = (3 + st->table.n) * seg * sizeof(vscore_t);
  hb_block_args_t args[nblocks];
  hb_boundary_t bounds[nblocks-1];
  pthread_t thread_ids[nblocks];
  bool started[nblocks];
  void *bufs = aln_aligned_alloc(sizeof(vscore_t), bufsize * nblocks);
  score_t *edges = aln_malloc(2 * (nblocks-1) * (rows+1) * sizeof(score_t));
  if(bufs == NULL || edges == NULL) {
    free(bufs);
    free(edges);
    hb_nomem(st);
    return;
  }

  for(i = 0; i + 1 < nblocks; i++) {
    bounds[i].H = edges + 2*i*(rows+1);
    bounds[i].F = bounds[i].H + (rows+1);
    atomic_init(&bounds[i].done, 0);
    atomic_init(&bounds[i].waiters, 0);
    if(pthread_mutex_init(&bounds[i].mu, NULL) != 0) {
      for(size_t j = 0; j < i; j++) {
        pthread_cond_destroy(&bounds[j].cv);
        pthread_mutex_destroy(&bounds[j].mu);
      }
      free(bufs);
      free(edges);
      hb_nomem(st);
      return;
    }
    if(pthread_cond_init(&bounds[i].cv, NULL) != 0) {
      pthread_mutex_destroy(&bounds[i].mu);
      for(size_t j = 0; j < i; j++) {
        pthread_cond_destroy(&bounds[j].cv);
        pthread_mutex_destroy(&bounds[j].mu);
      }
      free(bufs);
      free(edges);
      hb_nomem(st);
      return;
    }
  }

  for(i = 0; i < nblocks; i++) {
    hb_block_args_t blk = {st, a, astep, rows, b, bstep, c0 + i * width,
                           i + 1 < nblocks ? width : n - i * width, t0,
                           oc, od, ostep, (char*)bufs + i * bufsize,
                           i > 0 ? &bounds[i-1] : NULL,
                           i + 1 < nblocks ? &bounds[i] : NULL, top};
    args[i] = blk;
  }

  // Block 0 runs here; the rest on their own threads. A block whose thread
  // cannot start runs here too, after the blocks to its left
  for(i = 1; i < nblocks; i++)
    started[i] = pthread_create(&thread_ids[i], NULL, hb_block_thread,
                                &args[i]) == 0;
  hb_block_thread(&args[0]);
  for(i = 1; i < nblocks; i++) {
    if(started[i]) pthread_join(thread_ids[i], NULL);
    else hb_block_thread(&args[i]);
  }
  for(i = 0; i + 1 < nblocks; i++) {
    pthread_cond_destroy(&bounds[i].cv);
    pthread_mutex_destroy(&bounds[i].mu);
  }

  free(bufs);
  free(edges);
}

// Banded version of hb_pass_striped() over all n columns: only cells on
// diagonals dlo <= j - i <= dhi (row i, column j of the pass) need to be
// right. The rows are cut into strips; each strip runs over just the columns
// its rows' diagonals reach, with HB_NEG for the column on its left and for
// top-row cells the strip above did not compute. Cells next to those come
// out too low, never too high, and cells whose best path stays inside the
// band are exact. Columns of the last row outside the strip get HB_NEG
static void hb_pass_banded(hb_state_t *st, const char *a, ptrdiff_t astep,
                           size_t rows, const char *b, ptrdiff_t bstep,
                           size_t n, score_t t0,
                           score_t *oc, score_t *od, ptrdiff_t ostep,
                           void *buf, ptrdiff_t dlo, ptrdiff_t dhi,
                           int threads)
{
  const score_t g = st->scoring->gap_open, h = st->scoring->gap_extend;
  size_t strip = st->opts.hb_band_strip_rows;
  if(strip == 0) strip = MAX2(256, (size_t)(dhi - dlo + 1) / 4);

  // Last row of the strip above, by column; the first strip starts at row 0
  score_t *top_h = aln_malloc(2 * (n+1) * sizeof(score_t));
  if(top_h == NULL) {
    hb_nomem(st);
    return;
  }
  score_t *top_d = top_h + n+1;
  size_t i0, j, hs, c0 = 0, jr = n, jr_prev = n;
  top_h[0] = top_d[0] = 0;
  for(j = 1; j <= n; j++) {
    top_h[j] = g + (score_t)j * h;
    top_d[j] = top_h[j] + g;
  }

  for(i0 = 0; i0 < rows; i0 += hs)
  {
    hs = MIN2(strip, rows - i0);
    // Columns the band reaches in rows i0+1 .. i0+hs
    ptrdiff_t lo = (ptrdiff_t)i0 + 1 + dlo, hi = (ptrdiff_t)(i0 + hs) + dhi;
    c0 = lo <= 1 ? 0 : (size_t)lo - 1;
    jr = hi >= (ptrdiff_t)n ? n : (size_t)hi;
    // Columns the strip above did not reach
    for(j = jr_prev + 1; j <= jr; j++) top_h[j] = top_d[j] = HB_NEG;
    jr_prev = jr;

    hb_top_t top = {top_h, top_d, i0 == 0 ? t0 : top_h[0], c0 > 0, false};
    hb_pass_striped(st, a + (ptrdiff_t)i0 * astep, astep, hs, b, bstep,
                    c0, jr - c0, t0, top_h, top_d, 1, buf, &top, threads);
  }

  for(j = 0; j <= n; j++) {
    bool in = j >= c0 && j <= jr;
    oc[(ptrdiff_t)j * ostep] = in ? top_h[j] : HB_NEG;
    od[(ptrdiff_t)j * ostep] = in ? top_d[j] : HB_NEG;
  }
  free(top_h);
}

// hb_pass_striped() over all n columns that also writes the rows of save
// (NULL: none) on the way: it runs in pieces, each starting from the row the
// one before wrote, as a strip of hb_pass_banded() starts from the strip
// above
static void hb_pass_saving(hb_state_t *st, const char *a, ptrdiff_t astep,
                           size_t rows, const char *b, ptrdiff_t bstep,
                           size_t n, score_t t0,
                           score_t *oc, score_t *od, ptrdiff_t ostep,
                           void *buf, int threads, const hb_chain_t *save)
{
  hb_top_t top;
  const hb_top_t *tp = NULL;
  size_t done = 0, k;
  for(k = 0; save != NULL && k < save->len; k++) {
    hb_pass_striped(st, a + (ptrdiff_t)done * astep, astep,
                    save->rows[k] - done, b, bstep, 0, n, t0,
                    save->h[k], save->d[k], 1, buf, tp, threads);
    if(__atomic_load_n(&st->nomem, __ATOMIC_RELAXED)) return;
    top = (hb_top_t){save->h[k], save->d[k], save->h[k][0], false, true};
    tp = &top;
    done = save->rows[k];
  }
  hb_pass_striped(st, a + (ptrdiff_t)done * astep, astep, rows - done,
                  b, bstep, 0, n, t0, oc, od, ostep, buf, tp, threads);
}

// Arguments for running hb_pass_saving() or hb_pass_banded() on another
// thread: banded when dlo <= dhi
typedef struct {
  hb_state_t *st; const char *a; ptrdiff_t astep; size_t rows;
  const char *b; ptrdiff_t bstep; size_t n; score_t t0;
  score_t *oc, *od; ptrdiff_t ostep; void *buf; ptrdiff_t dlo, dhi;
  int threads; const hb_chain_t *save;
} hb_pass_args_t;

static void* hb_pass_thread(void *ptr)
{
  hb_pass_args_t *p = ptr;
  if(p->dlo <= p->dhi)
    hb_pass_banded(p->st, p->a, p->astep, p->rows, p->b, p->bstep, p->n,
                   p->t0, p->oc, p->od, p->ostep, p->buf, p->dlo, p->dhi,
                   p->threads);
  else
    hb_pass_saving(p->st, p->a, p->astep, p->rows, p->b, p->bstep, p->n,
                   p->t0, p->oc, p->od, p->ostep, p->buf, p->threads,
                   p->save);
  return NULL;
}

// Half-width of the band an m x n sub-problem needs for a path scoring at
// least `lower` (see alignment_band_width())
static size_t hb_band(const hb_state_t *st, size_t m, size_t n,
                      score_t lower)
{
  return alignment_band_width(st->smax, st->scoring->gap_extend, m, n, lower);
}

// Is a band of half-width w over n columns narrow enough to be worth it?
static bool hb_band_pays(const hb_state_t *st, size_t m, size_t n, size_t w)
{
  if(!st->opts.hb_band || w >= n) return false;
  double cols = (double)(m > n ? m - n : n - m) + 2.0*w + 1;
  double strip = st->opts.hb_band_strip_rows
                 ? (double)st->opts.hb_band_strip_rows : MAX2(256, cols / 4);
  return cols + strip < 0.6 * (double)n;
}

#endif /* HB_STRIPED */

// Forward and reverse passes of hb_diff() without SIMD, or just one of them
// (do_f, do_r). Each also writes the rows of its save chain (NULL: none)
static void hb_passes_scalar(hb_state_t *st, const char *a, const char *b,
                             size_t m, size_t n, size_t imid,
                             score_t tb, score_t te, bool do_f, bool do_r,
                             const hb_chain_t *fsave, const hb_chain_t *rsave)
{
  const score_t g = st->scoring->gap_open, h = st->scoring->gap_extend;
  score_t *cc = st->cc, *dd = st->dd, *rr = st->rr, *ss = st->ss;
  score_t s, c, e, t;
  size_t i, j, k;
  // Substitution scores: column of the table for a[i], row for each b[j]
  const int *col;
  const size_t *row_b = st->row_b + (b - st->b);

  // Forward pass: best scores for a[0..imid) against b[0..j),
  // cc ending in any state, dd ending in a deletion
  if(do_f)
  {
    cc[0] = 0;
    t = g;
    for(j = 1; j <= n; j++) { t += h; cc[j] = t; dd[j] = t + g; }

    t = tb;
    for(i = 1, k = 0; i <= imid; i++)
    {
      col = st->table.scores + st->table.code[(unsigned char)a[i-1]];
      s = cc[0];
      t += h;
      c = t;
      cc[0] = c;
      e = t + g;
      for(j = 1; j <= n; j++)
      {
        e = MAX2(e, c + g) + h;
        dd[j] = MAX2(dd[j], cc[j] + g) + h;
        c = MAX3(dd[j], e, s + col[row_b[j-1]]);
        s = cc[j];
        cc[j] = c;
      }
      if(fsave != NULL && k < fsave->len && i == fsave->rows[k]) {
        memcpy(fsave->h[k], cc, (n+1) * sizeof(score_t));
        memcpy(fsave->d[k], dd, (n+1) * sizeof(score_t));
        fsave->d[k][0] = cc[0];
        k++;
      }
    }
    dd[0] = cc[0];
  }

  // Reverse pass: best scores for a[imid..m) against b[j..n),
  // rr starting in any state, ss starting with a deletion
  if(do_r)
  {
    rr[n] = 0;
    t = g;
    for(j = n; j-- > 0; ) { t += h; rr[j] = t; ss[j] = t + g; }

    t = te;
    for(i = m, k = 0; i-- > imid; )
    {
      col = st->table.scores + st->table.code[(unsigned char)a[i]];
      s = rr[n];
      t += h;
      c = t;
      rr[n] = c;
      e = t + g;
      for(j = n; j-- > 0; )
      {
        e = MAX2(e, c + g) + h;
        ss[j] = MAX2(ss[j], rr[j] + g) + h;
        c = MAX3(ss[j], e, s + col[row_b[j]]);
        s = rr[j];
        rr[j] = c;
      }
      // in pass order: column n first
      if(rsave != NULL && k < rsave->len && m - i == rsave->rows[k]) {
        for(j = 0; j <= n; j++) {
          rsave->h[k][n-j] = rr[j];
          rsave->d[k][n-j] = ss[j];
        }
        rsave->d[k][0] = rr[n];
        k++;
      }
    }
    ss[n] = rr[n];
  }
}

bool hirschberg_check_striped = false;
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

bool hirschberg_simd_available(void)
{
#ifdef HB_STRIPED
  return alignment_striped_enabled && HB_FAST_VMAX();
#else
  return false;
#endif
}

#ifdef HB_STRIPED
static score_t hb_midpoint(const hb_state_t *st, size_t n, size_t *midj,
                           int *type);

// Testing aid: redo the passes without SIMD and count any difference. A
// banded row may be lower than the scalar one but never higher, and must
// give the same split
static void hb_check_striped(hb_state_t *st, const char *a, const char *b,
                             size_t m, size_t n, size_t imid,
                             score_t tb, score_t te, bool banded)
{
  size_t bytes = (n+1) * sizeof(score_t), j;
  score_t *copy = malloc(4 * bytes);
  if(copy == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  size_t midj, midj2;
  int type, type2;
  score_t midc = hb_midpoint(st, n, &midj, &type);
  memcpy(copy, st->cc, bytes);
  memcpy(copy + (n+1), st->dd, bytes);
  memcpy(copy + 2*(n+1), st->rr, bytes);
  memcpy(copy + 3*(n+1), st->ss, bytes);
  hb_passes_scalar(st, a, b, m, n, imid, tb, te, true, true, NULL, NULL);
  bool differs;
  if(banded) {
    differs = hb_midpoint(st, n, &midj2, &type2) != midc ||
              midj2 != midj || type2 != type;
    for(j = 0; j <= n; j++)
      differs |= copy[j] > st->cc[j] || copy[n+1+j] > st->dd[j] ||
                 copy[2*(n+1)+j] > st->rr[j] || copy[3*(n+1)+j] > st->ss[j];
  } else {
    differs = memcmp(copy, st->cc, bytes) ||
              memcmp(copy + (n+1), st->dd, bytes) ||
              memcmp(copy + 2*(n+1), st->rr, bytes) ||
              memcmp(copy + 3*(n+1), st->ss, bytes);
  }
  // Alignments on several threads may count at once
  if(differs)
    __atomic_add_fetch(&hirschberg_striped_mismatches, 1, __ATOMIC_RELAXED);
  free(copy);
}
#endif

// Where does the best path cross row imid, given the forward and reverse
// rows? Returns its score; *type 2 means a deletion spans the row: a[imid-1]
// and a[imid] are both deleted. The first best column wins
static score_t hb_midpoint(const hb_state_t *st, size_t n, size_t *midj,
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

#ifdef HB_STRIPED
// Forward and reverse passes of hb_diff() with SIMD, over diagonals
// dlo .. dhi only when dlo <= dhi (see hb_pass_banded()), or just one of
// them (do_f, do_r). Passes over all columns also write the rows of their
// save chain (NULL: none). The passes read shared, unchanging data and write
// different rows, so big ones run at the same time
static bool hb_ensure_vbuf(hb_state_t *st, bool need_second)
{
  if(st->vbuf == NULL) {
    st->vbuf = aln_aligned_alloc(sizeof(vscore_t), st->vbuf_bytes);
    if(st->vbuf == NULL) { hb_nomem(st); return false; }
  }
  if(need_second && st->vbuf2 == NULL) {
    st->vbuf2 = aln_aligned_alloc(sizeof(vscore_t), st->vbuf_bytes);
    if(st->vbuf2 == NULL) { hb_nomem(st); return false; }
  }
  return true;
}

static void hb_passes_striped(hb_state_t *st, const char *a, const char *b,
                              size_t m, size_t n, size_t imid,
                              score_t tb, score_t te,
                              ptrdiff_t dlo, ptrdiff_t dhi,
                              bool do_f, bool do_r, const hb_chain_t *fsave,
                              const hb_chain_t *rsave)
{
  if(!hb_ensure_vbuf(st, do_r)) return;
  // In the reverse pass, diagonal j - i becomes (n - m) - (j - i): the band
  // around 0 .. n-m maps onto itself
  hb_pass_args_t fwd = {st, a, 1, imid, b, 1, n, tb,
                        st->cc, st->dd, 1, st->vbuf, dlo, dhi,
                        st->pass_threads, fsave};
  hb_pass_args_t rev = {st, a+m-1, -1, m-imid, b+n-1, -1, n, te,
                        st->rr+n, st->ss+n, -1, st->vbuf2, dlo, dhi,
                        st->pass_threads, rsave};
  pthread_t thread;
  // Passes run at the same time share the threads: with the full count
  // each, twice as many threads as cores slowed 50 kb unrelated pairs
  // (0.31 -> 0.39 s) and tripled the CPU time. With one thread, they run
  // one after the other
  if(do_f && do_r && st->pass_threads >= 2 &&
     (double)m * (double)n >= HB_THREAD_CELLS)
  {
    fwd.threads = rev.threads = MAX2(1, st->pass_threads / 2);
    if(pthread_create(&thread, NULL, hb_pass_thread, &rev) == 0) {
      hb_pass_thread(&fwd);
      pthread_join(thread, NULL);
      return;
    }
    fwd.threads = rev.threads = st->pass_threads; // one after the other
  }
  if(do_f) hb_pass_thread(&fwd);
  if(do_r) hb_pass_thread(&rev);
}
#endif

#ifdef HB_STRIPED
// Can the SIMD passes be used on an m x n sub-problem? Building the score
// profile costs O(n x alphabet), so only on those with enough rows to pay
// for it
static bool hb_simd_ok(const hb_state_t *st, size_t m, size_t n)
{
  return st->opts.striped && st->simd_ok && HB_FAST_VMAX() &&
         st->scoring->gap_open <= 0 && st->scoring->gap_extend <= 0 &&
         n >= 2*ALN_VL && m >= 32 && m >= 2*st->table.n;
}

// Banded forward and reverse passes of hb_diff() (see hb_passes()). Cells
// too far from the main diagonals to be on a path scoring `lower` are
// skipped. Banded rows are never too high and are exact along every best
// path inside the band, so if the best score found (a real alignment's)
// needs no wider band, the band held every best path and the split is the
// same as without a band. Otherwise widen it and try again; with no bound
// to start from, try a narrow band first: it costs little, and its score
// shows how wide the band must be. Returns false, with nothing
// usable in the rows, once a band would be too wide to pay
static bool hb_passes_banded(hb_state_t *st, const char *a, const char *b,
                             size_t m, size_t n, size_t imid,
                             score_t tb, score_t te, score_t lower)
{
  size_t w = lower == HB_NO_BOUND
             ? st->opts.hb_band_min_width
             : hb_band(st, m, n, lower);
  const ptrdiff_t D = (ptrdiff_t)n - (ptrdiff_t)m;
  size_t midj;
  int type;

  while(hb_band_pays(st, m, n, w))
  {
    hb_passes_striped(st, a, b, m, n, imid, tb, te,
                      MIN2(0, D) - (ptrdiff_t)w, MAX2(0, D) + (ptrdiff_t)w,
                      true, true, NULL, NULL);
    if(st->nomem) return true; // hb_diff() stops
    if(hirschberg_check_striped)
      __atomic_add_fetch(&hirschberg_banded_passes, 1, __ATOMIC_RELAXED);
    size_t need = hb_band(st, m, n, hb_midpoint(st, n, &midj, &type));
    if(need <= w) {
      if(hirschberg_check_striped)
        hb_check_striped(st, a, b, m, n, imid, tb, te, true);
      return true;
    }
    w = MAX2(need, 2*w);
  }
  return false;
}
#endif

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
  if(hirschberg_myers_enabled && tb == 0 && te == 0 &&
#ifdef HB_STRIPED
     st->pass_threads <= 1 &&
#endif
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
  if(opts != NULL) st->opts = *opts;
  else alignment_options_default(&st->opts);

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
    double sub = 0;
    size_t c;
    st->smax = 0;
    for(c = 0; c < st->table.n * st->table.n; c++) {
      sub = MAX2(sub, HB_DABS(st->table.scores[c]));
      st->smax = MAX2(st->smax, st->table.scores[c]);
    }
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

// Exits if out of memory (ALN_ERR_NOMEM)
static bool hb_exit_nomem(aln_status_t status)
{
  if(status == ALN_ERR_NOMEM) {
    aln_error_t e = {status, 0, 0};
    alignment_fprint_error(stderr, &e);
    exit(EXIT_FAILURE);
  }
  return status == ALN_OK;
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
  return hb_exit_nomem(hb_align(a, b, len_a, len_b, scoring, opts, result));
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
  return hb_exit_nomem(hb_score(a, b, len_a, len_b, scoring, opts, score));
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
  if((
#ifdef HB_STRIPED
      st.pass_threads <= 1
        ? hb_myers(&st, a, len_a, b, len_b, score)
        : (hirschberg_myers_enabled &&
           alignment_myers_score_banded(&st.table, scoring, a, len_a, b, len_b,
                                        score))
#else
      hb_myers(&st, a, len_a, b, len_b, score)
#endif
     ))
  {
  }
  #ifdef HB_STRIPED
  // Banded forward and reverse passes of hb_diff(), which run at the same
  // time: every alignment crosses the middle row, so the best crossing is
  // the best score
  else if(hb_simd_ok(&st, len_a, len_b) &&
          hb_passes_banded(&st, a, b, len_a, len_b, len_a / 2, g, g, HB_NO_BOUND))
  {
    if(!st.nomem) *score = hb_midpoint(&st, len_b, &midj, &type);
  }
  // Otherwise one forward pass over every row of a: cc[len_b] is the best
  // score, ending in any state
  else if(hb_simd_ok(&st, len_a, len_b) && hb_ensure_vbuf(&st, false))
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
