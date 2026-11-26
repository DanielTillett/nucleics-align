/*
 hirschberg_pass.c
 Forward/reverse SIMD, banded, and scalar row passes for Hirschberg alignment.
 */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <pthread.h>
#include <stdatomic.h>

#include "hirschberg_internal.h"

#ifdef HB_STRIPED

// What a block passes to the block on its right, for each row r = 0..rows:
// H[r] best score in its last column, F[r] the insertion gap score entering
// the next column. done counts the rows published so far (64-byte aligned
// to avoid false sharing between adjacent blocks)
typedef struct hb_boundary {
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
typedef struct hb_top {
  const score_t *H, *D; // top row, best score and ending in a deletion,
                        // indexed by column of the pass
  score_t col0; // running deletion score in the left column, for column 0
  bool neg_left; // left column outside the band: HB_NEG in every row
  bool real; // every cell of the top row is a real score (no HB_NEG)
} hb_top_t;

// 16-bit lanes, twice as many per vector, for hb_block_pass16()
typedef int16_t v16_t __attribute__((vector_size(HB_VL16 * sizeof(int16_t))));

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
static void hb_block_pass16(const hb_block_args_t *blk)
{
  hb_state_t *const st = blk->st;
  const char *const a = blk->a, *const b = blk->b;
  const ptrdiff_t astep = blk->astep, bstep = blk->bstep, ostep = blk->ostep;
  const size_t rows = blk->rows, c0 = blk->c0, width = blk->width;
  const score_t t0 = blk->t0;
  score_t *const oc = blk->oc, *const od = blk->od;
  void *const buf = blk->buf;
  hb_boundary_t *const left = blk->left, *const right = blk->right;
  const hb_top_t *const top = blk->top;
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
static void hb_block_pass(const hb_block_args_t *blk)
{
  hb_state_t *const st = blk->st;
  const char *const a = blk->a, *const b = blk->b;
  const ptrdiff_t astep = blk->astep, bstep = blk->bstep, ostep = blk->ostep;
  const size_t rows = blk->rows, c0 = blk->c0, width = blk->width;
  const score_t t0 = blk->t0;
  score_t *const oc = blk->oc, *const od = blk->od;
  void *const buf = blk->buf;
  hb_boundary_t *const left = blk->left, *const right = blk->right;
  const hb_top_t *const top = blk->top;
  // Twice the lanes where the scores fit (see hb_block_pass16())
  // (parts of at least 4 vectors' columns keep its rows and profiles within
  // the space this block's 32-bit ones take; at most HB16_MAX_PARTS parts,
  // as their list is on the stack)
  if((top == NULL || top->real) && st->w16 >= 4 * HB_VL16 &&
     (width + st->w16 - 1) / st->w16 <= HB16_MAX_PARTS &&
     (right == NULL || width % HB_VL16 == 0)) {
    hb_block_pass16(blk);
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
static void* hb_block_thread(void *ptr)
{
  hb_block_pass(ptr);
  return NULL;
}

// One pass of hb_diff() over the n columns c0+1 .. c0+n: as one block, or
// split into column blocks on several threads when it is big enough (see
// hb_block_pass()). A whole pass has c0 = 0 and top = NULL
void hb_pass_striped(const hb_block_args_t *pass, int threads)
{
  hb_state_t *const st = pass->st;
  const size_t rows = pass->rows, c0 = pass->c0, n = pass->width;
  // At least one vector of columns per block, so no block is empty
  const size_t min_cols = MAX2(st->opts.hb_block_min_cols, (size_t)ALN_VL);
  size_t nblocks = 1, i;
  if(threads > 1 && (double)rows * n >= st->opts.hb_block_min_cells &&
     n >= 2 * min_cols)
  {
    nblocks = MIN2((size_t)threads, n / min_cols);
  }

  if(nblocks <= 1) {
    hb_block_pass(pass);
    return;
  }

  // Blocks of equal width (a multiple of ALN_VL, or of HB_VL16 so they can
  // take 16-bit lanes); the last takes the rest
  const size_t unit = n / nblocks >= HB_VL16 ? HB_VL16 : ALN_VL;
  size_t width = (n / nblocks) / unit * unit;
  size_t seg = (n - width * (nblocks-1) + ALN_VL - 1) / ALN_VL; // widest
  size_t bufsize = (3 + st->table.n) * seg * sizeof(vscore_t);
  // Per-block state on the heap: this runs on a pass thread, whose stack
  // (512 KB on macOS) thousands of blocks would not fit
  hb_block_args_t *args = aln_malloc(nblocks * sizeof(*args));
  hb_boundary_t *bounds = aln_malloc((nblocks-1) * sizeof(*bounds));
  pthread_t *thread_ids = aln_malloc(nblocks * sizeof(*thread_ids));
  bool *started = aln_malloc(nblocks * sizeof(*started));
  void *bufs = aln_aligned_alloc(sizeof(vscore_t), bufsize * nblocks);
  score_t *edges = aln_malloc(2 * (nblocks-1) * (rows+1) * sizeof(score_t));
  if(args == NULL || bounds == NULL || thread_ids == NULL || started == NULL ||
     bufs == NULL || edges == NULL) {
    free(args); free(bounds); free(thread_ids); free(started);
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
    const bool mu_ok = pthread_mutex_init(&bounds[i].mu, NULL) == 0;
    if(!mu_ok || pthread_cond_init(&bounds[i].cv, NULL) != 0) {
      if(mu_ok) pthread_mutex_destroy(&bounds[i].mu);
      for(size_t j = 0; j < i; j++) {
        pthread_cond_destroy(&bounds[j].cv);
        pthread_mutex_destroy(&bounds[j].mu);
      }
      free(args); free(bounds); free(thread_ids); free(started);
      free(bufs);
      free(edges);
      hb_nomem(st);
      return;
    }
  }

  for(i = 0; i < nblocks; i++) {
    args[i] = *pass;
    args[i].c0 = c0 + i * width;
    args[i].width = i + 1 < nblocks ? width : n - i * width;
    args[i].buf = (char*)bufs + i * bufsize;
    args[i].left = i > 0 ? &bounds[i-1] : NULL;
    args[i].right = i + 1 < nblocks ? &bounds[i] : NULL;
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

  free(args); free(bounds); free(thread_ids); free(started);
  free(bufs);
  free(edges);
}

// One whole pass of hb_diff() over the n columns 1 .. n (as
// hb_block_args_t, without the block fields), on up to threads threads:
// banded over the diagonals dlo .. dhi when dlo <= dhi
// (hb_pass_banded()), or saving the rows of save on the way
// (hb_pass_saving())
typedef struct {
  hb_state_t *st; const char *a; ptrdiff_t astep; size_t rows;
  const char *b; ptrdiff_t bstep; size_t n; score_t t0;
  score_t *oc, *od; ptrdiff_t ostep; void *buf; ptrdiff_t dlo, dhi;
  int threads; const hb_chain_t *save;
} hb_pass_args_t;

// The block a pass runs as: the pass's fields, from row `from` of a, over
// rows rows and the columns c0+1 .. c0+width, writing oc / od by column
// with stride ostep (the pass's own, or 1 for a strip's top row)
static hb_block_args_t hb_pass_block(const hb_pass_args_t *p, size_t from,
                                     size_t rows, size_t c0, size_t width,
                                     score_t *oc, score_t *od,
                                     ptrdiff_t ostep, const hb_top_t *top)
{
  hb_block_args_t blk = {p->st, p->a + (ptrdiff_t)from * p->astep, p->astep,
                         rows, p->b, p->bstep, c0, width, p->t0, oc, od,
                         ostep, p->buf, NULL, NULL, top};
  return blk;
}

// Banded version of hb_pass_striped() over all n columns: only cells on
// diagonals dlo <= j - i <= dhi (row i, column j of the pass) need to be
// right. The rows are cut into strips; each strip runs over just the columns
// its rows' diagonals reach, with HB_NEG for the column on its left and for
// top-row cells the strip above did not compute. Cells next to those come
// out too low, never too high, and cells whose best path stays inside the
// band are exact. Columns of the last row outside the strip get HB_NEG
static void hb_pass_banded(const hb_pass_args_t *p)
{
  hb_state_t *const st = p->st;
  const size_t rows = p->rows, n = p->n;
  const score_t t0 = p->t0;
  score_t *const oc = p->oc, *const od = p->od;
  const ptrdiff_t ostep = p->ostep, dlo = p->dlo, dhi = p->dhi;
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
    const hb_block_args_t strip_pass =
      hb_pass_block(p, i0, hs, c0, jr - c0, top_h, top_d, 1, &top);
    hb_pass_striped(&strip_pass, p->threads);
    // Out of memory in a strip: the ones below would build on stale rows
    if(__atomic_load_n(&st->nomem, __ATOMIC_RELAXED)) {
      free(top_h);
      return;
    }
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
static void hb_pass_saving(const hb_pass_args_t *p)
{
  const hb_chain_t *save = p->save;
  hb_top_t top;
  const hb_top_t *tp = NULL;
  size_t done = 0, k;
  for(k = 0; save != NULL && k < save->len; k++) {
    const hb_block_args_t piece =
      hb_pass_block(p, done, save->rows[k] - done, 0, p->n,
                    save->h[k], save->d[k], 1, tp);
    hb_pass_striped(&piece, p->threads);
    if(__atomic_load_n(&p->st->nomem, __ATOMIC_RELAXED)) return;
    top = (hb_top_t){save->h[k], save->d[k], save->h[k][0], false, true};
    tp = &top;
    done = save->rows[k];
  }
  const hb_block_args_t rest =
    hb_pass_block(p, done, p->rows - done, 0, p->n, p->oc, p->od, p->ostep,
                  tp);
  hb_pass_striped(&rest, p->threads);
}

static void* hb_pass_thread(void *ptr)
{
  const hb_pass_args_t *p = ptr;
  if(p->dlo <= p->dhi) hb_pass_banded(p);
  else hb_pass_saving(p);
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
void hb_passes_scalar(hb_state_t *st, const char *a, const char *b,
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

#ifdef HB_STRIPED
// Testing aid: redo the passes without SIMD and count any difference. A
// banded row may be lower than the scalar one but never higher, and must
// give the same split
void hb_check_striped(hb_state_t *st, const char *a, const char *b,
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

bool hb_ensure_vbuf(hb_state_t *st, bool need_second)
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

void hb_passes_striped(hb_state_t *st, const char *a, const char *b,
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

// Can the SIMD passes be used on an m x n sub-problem? Building the score
// profile costs O(n x alphabet), so only on those with enough rows to pay
// for it
bool hb_simd_ok(const hb_state_t *st, size_t m, size_t n)
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
bool hb_passes_banded(hb_state_t *st, const char *a, const char *b,
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
