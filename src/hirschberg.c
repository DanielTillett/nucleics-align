/*
 hirschberg.c
 url: https://github.com/DanielTillett/seq-align
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
#include <pthread.h>
#include <sched.h> // sched_yield
#include <stdatomic.h>
#include <unistd.h> // sysconf

#include "alignment_macros.h"
#include "alignment_simd.h"
#include "hirschberg.h"

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
#ifdef HB_STRIPED
  void *vbuf; // striped rows and score profile, see hb_pass_striped()
  void *vbuf2; // the same for a pass on a second thread
  bool simd_ok; // all scores far enough from HB_NEG and overflow
  int pass_threads; // threads for one pass, see hb_pass_striped()
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
  }
}

static void hb_emit_ins(hb_state_t *st, const char *b, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++, st->pos++) {
    st->out_a[st->pos] = '-';
    st->out_b[st->pos] = b[i];
  }
}

// Column-block threading of a pass: blocks of columns run on separate
// threads, each one row behind the block to its left. Tests lower these
size_t hirschberg_block_min_cells = 4000000; // rows x columns
size_t hirschberg_block_min_cols = 1024; // columns per block
int hirschberg_max_threads = 0; // threads per pass; 0: half the CPUs

#ifdef HB_STRIPED

// Below any real score, with room to add penalties without overflow
#define HB_NEG (SCORE_MIN / 4)

// Sub-problems with at least this many cells run their two passes on two
// threads; below it, starting a thread costs more than it saves
#define HB_THREAD_CELLS 1e6


// What a block passes to the block on its right, for each row r = 0..rows:
// H[r] best score in its last column, F[r] the insertion gap score entering
// the next column. done counts the rows published so far
typedef struct {
  score_t *H, *F;
  atomic_size_t done;
} hb_boundary_t;

// Wait until the block on the left has published row r
static inline void hb_wait_row(hb_boundary_t *left, size_t r)
{
  unsigned spins = 0;
  while(atomic_load_explicit(&left->done, memory_order_acquire) <= r)
    if(++spins > 256) { spins = 0; sched_yield(); }
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
// row, so Farrar's lazy F loop carries it across segment boundaries. That is
// exact only for gap_open <= 0 and gap_extend <= 0, which the caller checks.
// left / right are the boundaries shared with neighbouring blocks (NULL at
// the ends). A block with a right neighbour must have width a multiple of
// ALN_VL, so its last column is the last cell of its last lane.
ALN_TARGET_CLONES
static void hb_block_pass(hb_state_t *st, const char *a, ptrdiff_t astep,
                          size_t rows, const char *b, ptrdiff_t bstep,
                          size_t c0, size_t width, score_t t0,
                          score_t *oc, score_t *od, ptrdiff_t ostep,
                          void *buf, hb_boundary_t *left, hb_boundary_t *right)
{
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

  // Row 0: only insertions
  for(k = 0; k < seg; k++)
    for(l = 0; l < ALN_VL; l++) {
      H[k][l] = g + (score_t)(c0 + l*seg + k + 1) * h;
      D[k][l] = H[k][l] + g;
    }
  if(right != NULL) {
    right->H[0] = g + (score_t)(c0 + width) * h;
    right->F[0] = HB_NEG; // not used
    atomic_store_explicit(&right->done, 1, memory_order_release);
  }

  // The column left of the block: column 0 (only deletions) for the first
  // block, otherwise the last column of the block on the left
  score_t col0 = t0, col0_prev = 0, fstart;

  for(r = 0; r < rows; r++)
  {
    const vscore_t *pr = prof +
      (size_t)t->code[(unsigned char)a[(ptrdiff_t)r * astep]] * seg;

    if(left == NULL) {
      col0 += h;
      fstart = col0 + g + h;
    } else {
      hb_wait_row(left, r+1);
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

    // Lazy F: carry the insertion gap across segment boundaries while it
    // could still raise a best score
    vf = vshift(vf, HB_NEG);
    k = 0;
    while(vany(vf > Hn[k] + vg)) {
      Hn[k] = vmax(Hn[k], vf);
      vf = vmax(vf + vh, vneg);
      if(++k == seg) {
        fout = MAX2(fout, vf[ALN_VL-1]);
        k = 0;
        vf = vshift(vf, HB_NEG);
      }
    }

    if(right != NULL) {
      right->H[r+1] = Hn[seg-1][ALN_VL-1];
      right->F[r+1] = fout;
      atomic_store_explicit(&right->done, r+2, memory_order_release);
    }

    tmp = H; H = Hn; Hn = tmp;
    col0_prev = col0;
  }

  if(left == NULL) oc[0] = od[0] = col0;
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
  hb_boundary_t *left, *right;
} hb_block_args_t;

static void* hb_block_thread(void *ptr)
{
  hb_block_args_t *p = ptr;
  hb_block_pass(p->st, p->a, p->astep, p->rows, p->b, p->bstep, p->c0,
                p->width, p->t0, p->oc, p->od, p->ostep, p->buf,
                p->left, p->right);
  return NULL;
}

// One pass of hb_diff() over n columns: as one block, or split into column
// blocks on several threads when it is big enough (see hb_block_pass())
static void hb_pass_striped(hb_state_t *st, const char *a, ptrdiff_t astep,
                            size_t rows, const char *b, ptrdiff_t bstep,
                            size_t n, score_t t0,
                            score_t *oc, score_t *od, ptrdiff_t ostep,
                            void *buf)
{
  size_t nblocks = 1, i;
  if(st->pass_threads > 1 && (double)rows * n >= hirschberg_block_min_cells &&
     n >= 2 * hirschberg_block_min_cols)
  {
    nblocks = MIN2((size_t)st->pass_threads, n / hirschberg_block_min_cols);
  }

  if(nblocks <= 1) {
    hb_block_pass(st, a, astep, rows, b, bstep, 0, n, t0, oc, od, ostep, buf,
                  NULL, NULL);
    return;
  }

  // Blocks of equal width (a multiple of ALN_VL); the last takes the rest
  size_t width = (n / nblocks) / ALN_VL * ALN_VL;
  size_t seg = (n - width * (nblocks-1) + ALN_VL - 1) / ALN_VL; // widest
  size_t bufsize = (3 + st->table.n) * seg * sizeof(vscore_t);
  hb_block_args_t args[nblocks];
  hb_boundary_t bounds[nblocks-1];
  pthread_t threads[nblocks];
  bool started[nblocks];
  void *bufs = aligned_alloc(sizeof(vscore_t), bufsize * nblocks);
  score_t *edges = malloc(2 * (nblocks-1) * (rows+1) * sizeof(score_t));
  if(bufs == NULL || edges == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  for(i = 0; i + 1 < nblocks; i++) {
    bounds[i].H = edges + 2*i*(rows+1);
    bounds[i].F = bounds[i].H + (rows+1);
    atomic_init(&bounds[i].done, 0);
  }

  for(i = 0; i < nblocks; i++) {
    hb_block_args_t blk = {st, a, astep, rows, b, bstep, i * width,
                           i + 1 < nblocks ? width : n - i * width, t0,
                           oc, od, ostep, (char*)bufs + i * bufsize,
                           i > 0 ? &bounds[i-1] : NULL,
                           i + 1 < nblocks ? &bounds[i] : NULL};
    args[i] = blk;
  }

  // Block 0 runs here; the rest on their own threads. A block whose thread
  // cannot start runs here too, after the blocks to its left
  for(i = 1; i < nblocks; i++)
    started[i] = pthread_create(&threads[i], NULL, hb_block_thread,
                                &args[i]) == 0;
  hb_block_thread(&args[0]);
  for(i = 1; i < nblocks; i++) {
    if(started[i]) pthread_join(threads[i], NULL);
    else hb_block_thread(&args[i]);
  }

  free(bufs);
  free(edges);
}

// Arguments for running hb_pass_striped() on another thread
typedef struct {
  hb_state_t *st; const char *a; ptrdiff_t astep; size_t rows;
  const char *b; ptrdiff_t bstep; size_t n; score_t t0;
  score_t *oc, *od; ptrdiff_t ostep; void *buf;
} hb_pass_args_t;

static void* hb_pass_thread(void *ptr)
{
  hb_pass_args_t *p = ptr;
  hb_pass_striped(p->st, p->a, p->astep, p->rows, p->b, p->bstep, p->n,
                  p->t0, p->oc, p->od, p->ostep, p->buf);
  return NULL;
}

#endif /* HB_STRIPED */

// Forward and reverse passes of hb_diff() without SIMD
static void hb_passes_scalar(hb_state_t *st, const char *a, const char *b,
                             size_t m, size_t n, size_t imid,
                             score_t tb, score_t te)
{
  const score_t g = st->scoring->gap_open, h = st->scoring->gap_extend;
  score_t *cc = st->cc, *dd = st->dd, *rr = st->rr, *ss = st->ss;
  score_t s, c, e, t;
  size_t i, j;

  // Forward pass: best scores for a[0..imid) against b[0..j),
  // cc ending in any state, dd ending in a deletion
  cc[0] = 0;
  t = g;
  for(j = 1; j <= n; j++) { t += h; cc[j] = t; dd[j] = t + g; }

  t = tb;
  // Substitution scores: column of the table for a[i], row for each b[j]
  const int *col;
  const size_t *row_b = st->row_b + (b - st->b);

  for(i = 1; i <= imid; i++)
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
  }
  dd[0] = cc[0];

  // Reverse pass: best scores for a[imid..m) against b[j..n),
  // rr starting in any state, ss starting with a deletion
  rr[n] = 0;
  t = g;
  for(j = n; j-- > 0; ) { t += h; rr[j] = t; ss[j] = t + g; }

  t = te;
  for(i = m; i-- > imid; )
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
  }
  ss[n] = rr[n];
}

bool hirschberg_check_striped = false;
size_t hirschberg_striped_mismatches = 0;

#ifdef HB_STRIPED
// Testing aid: redo the passes without SIMD and count any difference
static void hb_check_striped(hb_state_t *st, const char *a, const char *b,
                             size_t m, size_t n, size_t imid,
                             score_t tb, score_t te)
{
  size_t bytes = (n+1) * sizeof(score_t);
  score_t *copy = malloc(4 * bytes);
  if(copy == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  memcpy(copy, st->cc, bytes);
  memcpy(copy + (n+1), st->dd, bytes);
  memcpy(copy + 2*(n+1), st->rr, bytes);
  memcpy(copy + 3*(n+1), st->ss, bytes);
  hb_passes_scalar(st, a, b, m, n, imid, tb, te);
  if(memcmp(copy, st->cc, bytes) || memcmp(copy + (n+1), st->dd, bytes) ||
     memcmp(copy + 2*(n+1), st->rr, bytes) ||
     memcmp(copy + 3*(n+1), st->ss, bytes))
  {
    hirschberg_striped_mismatches++;
  }
  free(copy);
}
#endif

// Align a[0..m) with b[0..n), appending to the output.
// tb / te: score of opening a deletion that touches the start / end of this
// sub-problem; 0 if it continues a deletion already paid for by the caller.
static void hb_diff(hb_state_t *st, const char *a, const char *b,
                    size_t m, size_t n, score_t tb, score_t te)
{
  const scoring_t *scoring = st->scoring;
  const score_t g = scoring->gap_open, h = scoring->gap_extend;
  size_t j;

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
      st->pos++;
      hb_emit_ins(st, b+best_j, n-best_j);
    }
    return;
  }

  score_t *cc = st->cc, *dd = st->dd, *rr = st->rr, *ss = st->ss;
  score_t c;
  size_t imid = m / 2;

  #ifdef HB_STRIPED
  // Building the score profile costs O(n x alphabet), so only use SIMD on
  // sub-problems with enough rows to pay for it
  if(alignment_striped_enabled && st->simd_ok && HB_FAST_VMAX() &&
     g <= 0 && h <= 0 &&
     n >= 2*ALN_VL && m >= 32 && m >= 2*st->table.n)
  {
    // The passes read shared, unchanging data and write different rows, so
    // big ones run at the same time
    hb_pass_args_t rev = {st, a+m-1, -1, m-imid, b+n-1, -1, n, te,
                          rr+n, ss+n, -1, st->vbuf2};
    pthread_t thread;
    if((double)m * (double)n >= HB_THREAD_CELLS &&
       pthread_create(&thread, NULL, hb_pass_thread, &rev) == 0)
    {
      hb_pass_striped(st, a, 1, imid, b, 1, n, tb, cc, dd, 1, st->vbuf);
      pthread_join(thread, NULL);
    }
    else
    {
      hb_pass_striped(st, a, 1, imid, b, 1, n, tb, cc, dd, 1, st->vbuf);
      hb_pass_thread(&rev);
    }
    if(hirschberg_check_striped) hb_check_striped(st, a, b, m, n, imid, tb, te);
  }
  else
  #endif
  {
    hb_passes_scalar(st, a, b, m, n, imid, tb, te);
  }

  // Where does the optimal path cross row imid? Type 2 means a deletion
  // spans the row: a[imid-1] and a[imid] are both deleted
  score_t midc = cc[0] + rr[0];
  size_t midj = 0;
  int type = 1;

  for(j = 0; j <= n; j++)
  {
    c = cc[j] + rr[j];
    if(c > midc) { midc = c; midj = j; type = 1; }
    c = dd[j] + ss[j] - g; // gap open counted on both sides
    if(c > midc) { midc = c; midj = j; type = 2; }
  }

  if(type == 1)
  {
    hb_diff(st, a, b, imid, midj, tb, g);
    hb_diff(st, a+imid, b+midj, m-imid, n-midj, g, te);
  }
  else
  {
    hb_diff(st, a, b, imid-1, midj, tb, 0);
    hb_emit_del(st, a+imid-1, 2);
    hb_diff(st, a+imid+1, b+midj, m-imid-1, n-midj, 0, te);
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

static void hb_state_init(hb_state_t *st, const char *a, const char *b,
                          size_t len_a, size_t len_b, const scoring_t *scoring)
{
  st->scoring = scoring;
  st->out_a = st->out_b = NULL;
  st->pos = 0;

  score_t *rows = malloc(4 * (len_b+1) * sizeof(score_t));
  if(rows == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  sub_table_build(&st->table, scoring, a, len_a, b, len_b);
  st->b = b;
  st->row_b = malloc((len_b+1) * sizeof(size_t));
  if(st->row_b == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  size_t j;
  for(j = 0; j < len_b; j++)
    st->row_b[j] = st->table.code[(unsigned char)b[j]] * st->table.n;

#ifdef HB_STRIPED
  // Three striped rows plus a profile row per character code
  size_t seg = (len_b + ALN_VL - 1) / ALN_VL + 1;
  // HB_NEG must stay below every real score, with room to add penalties.
  // No path through the matrices (including padding columns) has more than
  // len_a + len_b + ALN_VL steps, each scoring at most |gap_open| +
  // |gap_extend| or the largest substitution score in magnitude. Beyond a
  // conservative bound use the scalar passes
  {
    #define HB_DABS(x) ((x) < 0 ? -(double)(x) : (double)(x))
    double gap = HB_DABS(scoring->gap_open) + HB_DABS(scoring->gap_extend);
    double sub = 0;
    size_t c;
    for(c = 0; c < st->table.n * st->table.n; c++)
      sub = MAX2(sub, HB_DABS(st->table.scores[c]));
    double steps = (double)len_a + (double)len_b + ALN_VL + 2;
    st->simd_ok = steps * (gap + sub) < (double)(1 << 27);
  }

  // Threads per pass: the forward and reverse passes already run at once
  st->pass_threads = hirschberg_max_threads;
  if(st->pass_threads <= 0) {
    long cpus = sysconf(_SC_NPROCESSORS_ONLN);
    st->pass_threads = cpus > 1 ? (int)(cpus / 2) : 1;
  }

  st->vbuf = aligned_alloc(sizeof(vscore_t),
                           (3 + st->table.n) * seg * sizeof(vscore_t));
  st->vbuf2 = aligned_alloc(sizeof(vscore_t),
                            (3 + st->table.n) * seg * sizeof(vscore_t));
  if(st->vbuf == NULL || st->vbuf2 == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
#endif

  st->cc = rows;
  st->dd = rows + (len_b+1);
  st->rr = rows + 2*(len_b+1);
  st->ss = rows + 3*(len_b+1);
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
  if(hb_unsupported(scoring)) return false;

  alignment_ensure_capacity(result, len_a + len_b);

  hb_state_t st;
  hb_state_init(&st, a, b, len_a, len_b, scoring);
  st.out_a = result->result_a;
  st.out_b = result->result_b;

  hb_diff(&st, a, b, len_a, len_b, scoring->gap_open, scoring->gap_open);
  hb_state_free(&st);

  result->result_a[st.pos] = '\0';
  result->result_b[st.pos] = '\0';
  result->length = st.pos;
  result->score = hirschberg_score_alignment(result->result_a,
                                             result->result_b, scoring);
  return true;
}

bool hirschberg_score(const char *a, const char *b,
                      const scoring_t *scoring, score_t *score)
{
  return hirschberg_score2(a, b, strlen(a), strlen(b), scoring, score);
}

bool hirschberg_score2(const char *a, const char *b,
                       size_t len_a, size_t len_b,
                       const scoring_t *scoring, score_t *score)
{
  if(hb_unsupported(scoring)) return false;

  hb_state_t st;
  hb_state_init(&st, a, b, len_a, len_b, scoring);

  // One forward pass over every row of a: cc[len_b] is then the best score
  // of aligning all of a with all of b, ending in any state
  const score_t g = scoring->gap_open, h = scoring->gap_extend;
  #ifdef HB_STRIPED
  if(alignment_striped_enabled && st.simd_ok && HB_FAST_VMAX() &&
     g <= 0 && h <= 0 && len_b >= 2*ALN_VL)
  {
    hb_pass_striped(&st, a, 1, len_a, b, 1, len_b, g, st.cc, st.dd, 1,
                    st.vbuf);
  }
  else
  #endif
  {
    // with imid = len_a the reverse pass has no rows
    hb_passes_scalar(&st, a, b, len_a, len_b, len_a, g, g);
  }
  (void)h;

  *score = st.cc[len_b];
  hb_state_free(&st);
  return true;
}
