/*
 alignment.c
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Feb 2015
 */

// Turn on debugging output by defining SEQ_ALIGN_VERBOSE
//#define SEQ_ALIGN_VERBOSE

#define _DEFAULT_SOURCE // madvise

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower
#include <assert.h>
#include <limits.h> // INT_MAX
#include <pthread.h>
#include <sched.h> // sched_yield
#include <stdatomic.h>
#include <unistd.h> // sysconf
#ifdef __linux__
#  include <sys/mman.h> // madvise
#endif

#include "alignment.h"
#include "alignment_macros.h"
#include "alignment_simd.h"

const char align_col_mismatch[] = "\033[92m"; // Mismatch (GREEN)
const char align_col_indel[] = "\033[91m"; // Insertion / deletion (RED)
// Pink used by SmithWaterman local alignment for printing surrounding bases
const char align_col_context[] = "\033[95m";
const char align_col_stop[] = "\033[0m";

// Which sequences may not contain gaps once no_gaps_in_longest is resolved
static inline void alignment_no_gaps(const scoring_t *scoring,
                                     size_t len_a, size_t len_b,
                                     bool *no_gaps_in_a, bool *no_gaps_in_b)
{
  *no_gaps_in_a = scoring->no_gaps_in_a ||
                  (scoring->no_gaps_in_longest && len_a >= len_b);
  *no_gaps_in_b = scoring->no_gaps_in_b ||
                  (scoring->no_gaps_in_longest && len_a < len_b);
}

void sub_table_build(sub_table_t *t, const scoring_t *scoring,
                            const char *a, size_t len_a,
                            const char *b, size_t len_b)
{
  bool in_a[256] = {false}, in_b[256] = {false};
  unsigned char bytes[256];
  size_t i, j;

  for(i = 0; i < len_a; i++) in_a[(unsigned char)a[i]] = true;
  for(i = 0; i < len_b; i++) in_b[(unsigned char)b[i]] = true;

  for(i = 0, t->n = 0; i < 256; i++) {
    t->code[i] = 0;
    if(in_a[i] || in_b[i]) {
      t->code[i] = (unsigned char)t->n;
      bytes[t->n++] = (unsigned char)i;
    }
  }

  size_t cells = t->n * t->n + 1;
  t->scores = malloc(cells * sizeof(int));
  t->is_match = malloc(cells * sizeof(bool));
  if(t->scores == NULL || t->is_match == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  // Only look up pairs that can be aligned (a character of a against one of
  // b): scoring_lookup() exits on pairs a substitution file does not define
  for(j = 0; j < t->n; j++) {
    for(i = 0; i < t->n; i++) {
      if(in_a[bytes[i]] && in_b[bytes[j]]) {
        scoring_lookup(scoring, (char)bytes[i], (char)bytes[j],
                       &t->scores[j*t->n+i], &t->is_match[j*t->n+i]);
      } else {
        t->scores[j*t->n+i] = 0; // never read
        t->is_match[j*t->n+i] = false;
      }
    }
  }
}

void sub_table_free(sub_table_t *t)
{
  free(t->scores);
  free(t->is_match);
}

//
// Striped SIMD fill (Farrar 2007, Bioinformatics 23:156-161)
//
// Each row of the matrices runs along seq_a. A row is split into ALN_VL
// segments of seg cells; lane l of vector k holds cell l*seg+k. match and
// gap_a only depend on the previous row, so they are computed a vector at a
// time. gap_b depends on the cell to its left: it is first computed within
// each segment, then corrected by Farrar's 'lazy F' loop, which carries
// values across segment boundaries until they stop improving any cell.
// Every value is identical to the scalar recurrence, including the clamp at
// min, so traceback and Smith-Waterman hit order are unchanged.
// Vector helpers are in alignment_simd.h.

bool alignment_striped_enabled = true;

// Column-block threading of the striped fill: blocks of columns run on
// separate threads, each one row behind the block to its left. Tests lower
// these
size_t alignment_block_min_cells = 4000000; // rows x columns
size_t alignment_block_min_cols = 1024; // columns per block
int alignment_max_threads = 0; // 0: half the CPUs

// Banded compact fills (see alignment_fill_banded())
bool alignment_band_enabled = true;
size_t alignment_band_strip_rows = 0; // 0: automatic
size_t alignment_band_min_width = 64;

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

// Wait until the block on the left has finished row y
static inline void alignment_wait_row(atomic_size_t *left_done, size_t y)
{
  unsigned spins = 0;
  while(atomic_load_explicit(left_done, memory_order_acquire) < y)
    if(++spins > 256) { spins = 0; sched_yield(); }
}

// One column block of the striped fill: the columns of seq_a c0 ..
// c0+width-1 (matrix columns c0+1 .. c0+width). r0[m] is row 0 of matrix m
// (MATCH, GAP_A, GAP_B) by matrix column; lc[m] is matrix column c0 by row,
// lstride apart: the initialised first column for the first block and the
// last column of the block on the left otherwise. With a left block, row y
// is only started once left_done reaches y; done is set to y once row y is
// written.
// Without trace the three matrices are written. With trace each cell gets
// its traceback byte instead, rc[m] gets the block's last column and, with
// track_best, best_* the first cell (lowest x, then y) with the highest
// match score.
// Rows y0 .. y1 are filled (1 .. the last row for a whole fill): r0 is then
// row y0-1, and lc[m][(y - ly0) * lstride] is the left column in row y. A
// banded fill (see alignment_fill_banded()) stores the traceback byte of
// cell (x, y) at trace[trace_rows[y] + x], and writes the strip's last row
// to lastrow[m] by matrix column
typedef struct {
  aligner_t *aligner; const sub_table_t *table; score_t min;
  size_t c0, width; atomic_size_t *left_done, *done;
  const score_t *r0[3], *lc[3]; size_t lstride;
  uint8_t *trace; score_t *rc[3];
  bool track_best; score_t best_score; size_t best_x, best_y;
  size_t y0, y1, ly0;
  const ptrdiff_t *trace_rows;
  score_t *lastrow[3];
  // Smith-Waterman scans (alignment_sw_scan()): no traceback bytes, just
  // the best cell, or with target > 0 the highest column and row of cells
  // with that match score
  bool score_only;
  score_t target;
  size_t tgt_x, tgt_y;
} alignment_block_t;

// Lanes of a where mask is set, of b elsewhere
ALN_INLINE vscore_t vsel(vscore_t mask, vscore_t a, vscore_t b)
{
  return (a & mask) | (b & ~mask);
}

// Matrix the traceback steps to from a cell scoring cur, given each
// predecessor's score plus the penalty for stepping from it: the first of
// GAP_A, GAP_B and MATCH that gives cur, as alignment_reverse_move() checks
ALN_INLINE vscore_t vtrace(vscore_t cur, vscore_t from_a, vscore_t from_b,
                              vscore_t from_m)
{
  vscore_t c = vsplat(ALN_TRACE_NONE);
  c = vsel(from_m == cur, vsplat(MATCH), c);
  c = vsel(from_b == cur, vsplat(GAP_B), c);
  return vsel(from_a == cur, vsplat(GAP_A), c);
}

ALN_TARGET_CLONES
static void alignment_fill_striped(alignment_block_t *p)
{
  aligner_t *aligner = p->aligner;
  const sub_table_t *table = p->table;
  const scoring_t *scoring = aligner->scoring;
  const score_t min = p->min;
  const size_t c0 = p->c0, ls = p->lstride;
  const bool full = p->trace == NULL && !p->score_only;
  score_t *match_scores = full ? aligner->match_scores + c0 : NULL;
  score_t *gap_a_scores = full ? aligner->gap_a_scores + c0 : NULL;
  score_t *gap_b_scores = full ? aligner->gap_b_scores + c0 : NULL;
  const score_t *lcM = p->lc[MATCH], *lcA = p->lc[GAP_A], *lcB = p->lc[GAP_B];
  const char *seq_a = aligner->seq_a + c0;
  const size_t w = aligner->score_width;
  const size_t len_i = p->width, ly0 = p->ly0;
  const size_t seg = (len_i + ALN_VL - 1) / ALN_VL;
  const size_t n = table->n;
  size_t c, k, l, y, i;

  const vscore_t vmin = vsplat(min);
  const vscore_t vopen = vsplat(scoring->gap_open + scoring->gap_extend);
  const vscore_t vext = vsplat(scoring->gap_extend);

  // profile[c*seg + k] lane l: score of seq_a[l*seg+k] against character c
  // of seq_b (0 past the end of seq_a; those cells are never read back)
  size_t nvec = (n + 7) * seg;
  vscore_t *buf = aligned_alloc(sizeof(vscore_t), nvec * sizeof(vscore_t));
  if(buf == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  vscore_t *profile = buf;
  vscore_t *pM = buf + n*seg, *pA = pM + seg, *pB = pA + seg; // previous row
  vscore_t *cM = pB + seg, *cA = cM + seg, *cB = cA + seg; // current row
  vscore_t *tb = cB + seg; // traceback bytes of the current row
  vscore_t *tmp;

  for(c = 0; c < n; c++) {
    for(k = 0; k < seg; k++) {
      for(l = 0; l < ALN_VL; l++) {
        i = l*seg + k;
        profile[c*seg+k][l] = i < len_i
          ? table->scores[c*n + table->code[(unsigned char)seq_a[i]]]
          : 0;
      }
    }
  }

  // Row y0-1 (row 0 for a whole fill) in striped order
  for(k = 0; k < seg; k++) {
    for(l = 0; l < ALN_VL; l++) {
      i = l*seg + k;
      pM[k][l] = i < len_i ? p->r0[MATCH][c0+1+i] : min;
      pA[k][l] = i < len_i ? p->r0[GAP_A][c0+1+i] : min;
      pB[k][l] = i < len_i ? p->r0[GAP_B][c0+1+i] : min;
    }
  }
  if(!full)
    for(c = 0; c < 3; c++) p->rc[c][p->y0-1] = p->r0[c][c0+len_i];

  // Column of each lane's first cell, to leave out cells past the end of
  // seq_a when looking for the best cell
  vscore_t vfirst;
  for(l = 0; l < ALN_VL; l++) vfirst[l] = (score_t)(l*seg);
  const vscore_t vlen = vsplat((score_t)len_i);

  for(y = p->y0; y <= p->y1; y++)
  {
    const size_t row = y*w;
    const vscore_t *prof = profile +
      (size_t)table->code[(unsigned char)aligner->seq_b[y-1]] * seg;

    // Column c0 (initialised, or the left block's) feeds lane 0 at k = 0
    if(p->left_done != NULL) alignment_wait_row(p->left_done, y);
    const size_t lp = (y-1-ly0)*ls, lc = (y-ly0)*ls;
    const score_t lm_prev = lcM[lp], la_prev = lcA[lp], lb_prev = lcB[lp];
    const score_t lm = lcM[lc], la = lcA[lc], lb = lcB[lc];
    score_t h0_prev = MAX3(lm_prev, la_prev, lb_prev);
    score_t x0 = MAX2(lm, la);
    score_t b0 = lb; // min in the real column 0

    // match from the diagonal, gap_a from the cell above
    vscore_t vdiag = vshift(vmax(vmax(pM[seg-1], pA[seg-1]), pB[seg-1]),
                            h0_prev);
    for(k = 0; k < seg; k++) {
      cM[k] = vmax(vdiag + prof[k], vmin);
      cA[k] = vmax(vmax(vmax(pM[k], pB[k]) + vopen, pA[k] + vext), vmin);
      vdiag = vmax(vmax(pM[k], pA[k]), pB[k]);
    }

    // gap_b within each segment; lanes other than 0 start without the
    // gap_b value carried from the previous segment
    vscore_t vx = vshift(vmax(cM[seg-1], cA[seg-1]), x0);
    vscore_t vb = vshift(vmin, b0);
    for(k = 0; k < seg; k++) {
      cB[k] = vmax(vmax(vx + vopen, vb + vext), vmin);
      vx = vmax(cM[k], cA[k]);
      vb = cB[k];
    }

    // Lazy F: carry gap_b across segment boundaries while it improves a cell
    vscore_t vf = vmax(vshift(cB[seg-1], min) + vext, vmin);
    k = 0;
    while(vany(vf > cB[k])) {
      cB[k] = vmax(cB[k], vf);
      vf = vmax(vf + vext, vmin);
      if(++k == seg) { k = 0; vf = vmax(vshift(vf, min), vmin); }
    }

    if(full)
    {
      // Back to row-major for traceback
      for(l = 0; l < ALN_VL; l++) {
        for(k = 0, i = l*seg; k < seg && i < len_i; k++, i++) {
          match_scores[row+1+i] = cM[k][l];
          gap_a_scores[row+1+i] = cA[k][l];
          gap_b_scores[row+1+i] = cB[k][l];
        }
      }
    }
    else
    {
      if(!p->score_only)
      {
        // Traceback bytes, from each cell's predecessors: diagonal (match),
        // above (gap_a) and left (gap_b)
        vscore_t dM = vshift(pM[seg-1], lm_prev);
        vscore_t dA = vshift(pA[seg-1], la_prev);
        vscore_t dB = vshift(pB[seg-1], lb_prev);
        vscore_t eM = vshift(cM[seg-1], lm), eA = vshift(cA[seg-1], la);
        vscore_t eB = vshift(cB[seg-1], lb);
        for(k = 0; k < seg; k++) {
          tb[k] = vtrace(cM[k], dA + prof[k], dB + prof[k], dM + prof[k]) |
                  vtrace(cA[k], pA[k] + vext, pB[k] + vopen, pM[k] + vopen) << 2 |
                  vtrace(cB[k], eA + vopen, eB + vext, eM + vopen) << 4;
          dM = pM[k]; dA = pA[k]; dB = pB[k];
          eM = cM[k]; eA = cA[k]; eB = cB[k];
        }
        uint8_t *trow = p->trace + (p->trace_rows != NULL
                                    ? (size_t)(p->trace_rows[y] +
                                               (ptrdiff_t)(c0 + 1))
                                    : row + c0 + 1);
        for(l = 0; l < ALN_VL; l++)
          for(k = 0, i = l*seg; k < seg && i < len_i; k++, i++)
            trow[i] = (uint8_t)tb[k][l];
      }

      i = len_i - 1;
      p->rc[MATCH][y] = cM[i % seg][i / seg];
      p->rc[GAP_A][y] = cA[i % seg][i / seg];
      p->rc[GAP_B][y] = cB[i % seg][i / seg];

      if(p->track_best || p->target > 0)
      {
        vscore_t vm = vmin;
        for(k = 0; k < seg; k++)
          vm = vmax(vm, vsel(vfirst + (score_t)k < vlen, cM[k], vmin));
        score_t rmax = vm[0];
        for(l = 1; l < ALN_VL; l++) rmax = MAX2(rmax, vm[l]);
        if(p->track_best && rmax > 0 && rmax >= p->best_score) {
          // Lowest column in this row with that score: the first k in each
          // lane (a vector pass, as the best rises in most rows of similar
          // sequences), then the first lane with one
          const vscore_t vr = vsplat(rmax), vnone = vsplat((score_t)seg);
          vscore_t vk = vnone;
          for(k = 0; k < seg; k++)
            vk = vsel((cM[k] == vr) & (vk == vnone) &
                      (vfirst + (score_t)k < vlen), vsplat((score_t)k), vk);
          for(l = 0, i = len_i; l < ALN_VL && i == len_i; l++)
            if(vk[l] < (score_t)seg) i = l*seg + (size_t)vk[l];
          if(rmax > p->best_score || c0+1+i < p->best_x) {
            p->best_score = rmax; p->best_x = c0+1+i; p->best_y = y;
          }
        }
        if(p->target > 0 && rmax >= p->target) {
          // Highest column in this row with the target score (none score
          // more, see alignment_sw_scan())
          for(l = ALN_VL, i = len_i; l-- > 0 && i == len_i; )
            for(k = seg; k-- > 0; )
              if(l*seg+k < len_i && cM[k][l] == p->target) {
                i = l*seg+k;
                break;
              }
          if(i < len_i) {
            p->tgt_x = MAX2(p->tgt_x, c0+1+i);
            p->tgt_y = y;
          }
          p->best_score = MAX2(p->best_score, rmax);
        }
      }
    }
    if(p->done != NULL) atomic_store_explicit(p->done, y, memory_order_release);

    tmp = pM; pM = cM; cM = tmp;
    tmp = pA; pA = cA; cA = tmp;
    tmp = pB; pB = cB; cB = tmp;
  }

  // The last row, for the strip below
  if(p->lastrow[MATCH] != NULL) {
    for(l = 0; l < ALN_VL; l++) {
      for(k = 0, i = l*seg; k < seg && i < len_i; k++, i++) {
        p->lastrow[MATCH][c0+1+i] = pM[k][l];
        p->lastrow[GAP_A][c0+1+i] = pA[k][l];
        p->lastrow[GAP_B][c0+1+i] = pB[k][l];
      }
    }
  }

  free(buf);
}

static void* alignment_block_thread(void *ptr)
{
  alignment_fill_striped(ptr);
  return NULL;
}

// Threads (column blocks) for a striped fill of width x rows cells
static size_t alignment_fill_threads(size_t width, size_t rows)
{
  int threads = alignment_max_threads;
  if(threads <= 0) {
    long cpus = sysconf(_SC_NPROCESSORS_ONLN);
    threads = cpus > 1 ? (int)(cpus / 2) : 1;
  }
  // At least one vector of columns per block, so no block is empty
  const size_t min_cols = MAX2(alignment_block_min_cols, (size_t)ALN_VL);
  if(threads > 1 && (double)width * rows >= alignment_block_min_cells &&
     width >= 2 * min_cols)
    return MIN2((size_t)threads, width / min_cols);
  return 1;
}

// Striped fill of the columns tmpl->c0+1 .. tmpl->c0+tmpl->width over rows
// tmpl->y0 .. tmpl->y1 (the whole matrix, or a strip of a banded fill): as
// one block, or split into column blocks on several threads when it is big
// enough. tmpl gives the settings for the first block; the others read
// their left column from the block before. With trace, the last cell's
// scores go to final_scores (once its row and column are filled) and the
// best cell (track_best) to best_score and best_index
static void alignment_fill_blocks(alignment_block_t *tmpl)
{
  aligner_t *aligner = tmpl->aligner;
  const size_t w = aligner->score_width, h = aligner->score_height;
  const size_t len_i = tmpl->width, len_j = tmpl->y1 - tmpl->y0 + 1;
  const size_t nblocks = alignment_fill_threads(len_i, len_j);
  size_t i, m;

  // Blocks of equal width (a multiple of ALN_VL); the last takes the rest
  size_t width = (len_i / nblocks) / ALN_VL * ALN_VL;
  alignment_block_t args[nblocks];
  atomic_size_t done[nblocks];
  pthread_t thread_ids[nblocks];
  bool started[nblocks];

  // With trace, each block's last column goes here for the next block
  score_t *cols = NULL;
  if(tmpl->trace != NULL || tmpl->score_only) {
    if(h > SIZE_MAX / sizeof(score_t) / 3 / nblocks) cols = NULL;
    else cols = malloc(nblocks * 3 * h * sizeof(score_t));
    if(cols == NULL) {
      fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  for(i = 0; i < nblocks; i++) {
    atomic_init(&done[i], 0);
    alignment_block_t *b = &args[i];
    *b = *tmpl;
    b->c0 = tmpl->c0 + i * width;
    b->width = i + 1 < nblocks ? width : len_i - i * width;
    b->left_done = i > 0 ? &done[i-1] : NULL;
    b->done = i + 1 < nblocks ? &done[i] : NULL;
    b->best_score = 0;
    b->best_x = b->best_y = 0;
    for(m = 0; m < 3; m++) {
      if(cols != NULL) {
        b->rc[m] = cols + (i*3 + m) * h; // by row
        if(i > 0) { b->lc[m] = args[i-1].rc[m]; b->lstride = 1; b->ly0 = 0; }
      }
      else b->lc[m] = tmpl->lc[m] + (b->c0 - tmpl->c0); // the matrices
    }
  }

  // Block 0 runs here; the rest on their own threads. A block whose thread
  // cannot start runs here too, after the blocks to its left
  for(i = 1; i < nblocks; i++)
    started[i] = pthread_create(&thread_ids[i], NULL, alignment_block_thread,
                                &args[i]) == 0;
  alignment_block_thread(&args[0]);
  for(i = 1; i < nblocks; i++) {
    if(started[i]) pthread_join(thread_ids[i], NULL);
    else alignment_block_thread(&args[i]);
  }

  if(cols != NULL) {
    if(tmpl->y1 == h-1 && tmpl->c0 + tmpl->width == w-1)
      for(m = 0; m < 3; m++)
        aligner->final_scores[m] = args[nblocks-1].rc[m][h-1];
    // Blocks are in column order, so on a tie the earlier block's cell wins
    aligner->best_score = 0;
    aligner->best_index = 0;
    for(i = 0; i < nblocks; i++) {
      if(args[i].best_score > aligner->best_score) {
        aligner->best_score = args[i].best_score;
        aligner->best_index = args[i].best_y * w + args[i].best_x;
      }
    }
    // Target scans: the highest column and row, and the best score seen
    tmpl->tgt_x = tmpl->tgt_y = 0;
    tmpl->best_score = 0;
    for(i = 0; i < nblocks; i++) {
      tmpl->tgt_x = MAX2(tmpl->tgt_x, args[i].tgt_x);
      tmpl->tgt_y = MAX2(tmpl->tgt_y, args[i].tgt_y);
      tmpl->best_score = MAX2(tmpl->best_score, args[i].best_score);
    }
    free(cols);
  }
}

// Smith-Waterman scans of a against b in O(len_b) memory, for the best hit
// of smith_waterman.c. With target 0: the best match score and its first
// cell (lowest x, then y; *x, *y count from 1), as a full fill would rank
// it. With target > 0, for a scan in which no cell scores more than target:
// the highest x and the highest y of cells scoring target (0 if none).
// Returns false, with nothing found, where the striped fill does not apply
// or (target > 0) a cell scores more than target
bool alignment_sw_scan(const char *a, size_t len_a, const char *b, size_t len_b,
                       const scoring_t *scoring, score_t target,
                       score_t *best, size_t *x, size_t *y)
{
  bool no_gaps_in_a, no_gaps_in_b;
  alignment_no_gaps(scoring, len_a, len_b, &no_gaps_in_a, &no_gaps_in_b);
  if(!alignment_striped_enabled || len_a == 0 || len_b == 0 ||
     len_a >= INT_MAX || len_b >= INT_MAX ||
     len_a + len_b > SIZE_MAX / (3 * sizeof(score_t)) ||
     scoring->no_start_gap_penalty || scoring->no_end_gap_penalty ||
     scoring->no_mismatches || no_gaps_in_a || no_gaps_in_b) return false;

  aligner_t tmp;
  aligner_init(&tmp);
  tmp.scoring = scoring;
  tmp.seq_a = a;
  tmp.seq_b = b;
  tmp.score_width = len_a + 1;
  tmp.score_height = len_b + 1;

  // Row 0 and column 0 are all 0 in Smith-Waterman
  score_t *zeros = calloc(len_a + len_b + 2, sizeof(score_t));
  if(zeros == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  sub_table_t table;
  sub_table_build(&table, scoring, a, len_a, b, len_b);
  alignment_block_t tmpl = {
    .aligner = &tmp, .table = &table, .min = 0,
    .r0 = {zeros, zeros, zeros}, .lc = {zeros, zeros, zeros}, .lstride = 0,
    .width = len_a, .y0 = 1, .y1 = len_b,
    .score_only = true, .track_best = target == 0, .target = target};
  alignment_fill_blocks(&tmpl);
  sub_table_free(&table);
  free(zeros);

  if(target == 0) {
    *best = tmp.best_score;
    *x = tmp.best_index % tmp.score_width;
    *y = tmp.best_index / tmp.score_width;
    return true;
  }
  if(tmpl.best_score > target) return false;
  *best = tmpl.best_score;
  *x = tmpl.tgt_x;
  *y = tmpl.tgt_y;
  return true;
}

#else

bool alignment_sw_scan(const char *a, size_t len_a, const char *b, size_t len_b,
                       const scoring_t *scoring, score_t target,
                       score_t *best, size_t *x, size_t *y)
{
  (void)a; (void)len_a; (void)b; (void)len_b; (void)scoring; (void)target;
  (void)best; (void)x; (void)y;
  return false;
}

#endif /* ALN_STRIPED */

// Lowest score a cell may take (cells are clamped to it)
static score_t alignment_min_score(const scoring_t *scoring, char is_sw)
{
  int gap_open_penalty = scoring->gap_extend + scoring->gap_open;
  int gap_extend_penalty = scoring->gap_extend;

  // Penalties may have been changed since scoring_init() (e.g. by the command
  // line parser), so take the most negative penalty from the current values
  int worst_penalty = MIN3(scoring->min_penalty, gap_open_penalty, gap_extend_penalty);
  worst_penalty = MIN3(worst_penalty, scoring->match, scoring->mismatch);

  return is_sw ? 0 : SCORE_MIN + abs(worst_penalty);
}

// Fill in traceback matrix
static void alignment_fill_matrices(aligner_t *aligner, char is_sw)
{
  score_t *match_scores = aligner->match_scores;
  score_t *gap_a_scores = aligner->gap_a_scores;
  score_t *gap_b_scores = aligner->gap_b_scores;
  const scoring_t *scoring = aligner->scoring;
  size_t score_width = aligner->score_width;
  size_t score_height = aligner->score_height;
  size_t i, j;

  int gap_open_penalty = scoring->gap_extend + scoring->gap_open;
  int gap_extend_penalty = scoring->gap_extend;

  const score_t min = alignment_min_score(scoring, is_sw);

  size_t seq_i, seq_j, len_i = score_width-1, len_j = score_height-1;

  bool no_gaps_in_a, no_gaps_in_b;
  alignment_no_gaps(scoring, len_i, len_j, &no_gaps_in_a, &no_gaps_in_b);
  size_t index, index_left, index_up, index_upleft;

  sub_table_t table;
  sub_table_build(&table, scoring, aligner->seq_a, len_i,
                  aligner->seq_b, len_j);

  // [0][0]
  match_scores[0] = 0;
  gap_a_scores[0] = 0;
  gap_b_scores[0] = 0;

  if(is_sw)
  {
    for(i = 1; i < score_width; i++)
      match_scores[i] = gap_a_scores[i] = gap_b_scores[i] = 0;
    for(j = 1, index = score_width; j < score_height; j++, index += score_width)
      match_scores[index] = gap_a_scores[index] = gap_b_scores[index] = min;
  }
  else
  {
    // work along first row -> [i][0]
    for(i = 1; i < score_width; i++)
    {
      match_scores[i] = min;

      // Think carefully about which way round these are
      gap_a_scores[i] = min;
      gap_b_scores[i] = scoring->no_start_gap_penalty ? 0
                        : scoring->gap_open + (int)i * scoring->gap_extend;
    }

    // work down first column -> [0][j]
    for(j = 1, index = score_width; j < score_height; j++, index += score_width)
    {
      match_scores[index] = min;

      // Think carefully about which way round these are
      gap_a_scores[index] = scoring->no_start_gap_penalty ? 0
                            : scoring->gap_open + (int)j * scoring->gap_extend;
      gap_b_scores[index] = min;
    }
  }

  #ifdef ALN_STRIPED
  // The striped fill handles the uniform recurrence only; options that
  // change it for particular cells use the scalar loop below
  if(alignment_striped_enabled && len_i > 0 &&
     !scoring->no_end_gap_penalty && !scoring->no_mismatches &&
     !no_gaps_in_a && !no_gaps_in_b)
  {
    alignment_block_t tmpl = {
      .aligner = aligner, .table = &table, .min = min,
      .r0 = {match_scores, gap_a_scores, gap_b_scores},
      .lc = {match_scores, gap_a_scores, gap_b_scores},
      .lstride = score_width, .width = len_i, .y0 = 1, .y1 = len_j};
    alignment_fill_blocks(&tmpl);
    sub_table_free(&table);
    return;
  }
  #endif

  // start at position [1][1]
  index_upleft = 0;
  index_up = 1;
  index_left = score_width;
  index = score_width+1;

  for(seq_j = 0; seq_j < len_j; seq_j++)
  {
    size_t row = table.code[(unsigned char)aligner->seq_b[seq_j]] * table.n;
    const int *sub_row = table.scores + row;
    const bool *match_row = table.is_match + row;

    for(seq_i = 0; seq_i < len_i; seq_i++)
    {
      // Update match_scores[i][j] with position [i-1][j-1]
      // substitution penalty
      unsigned char code_a = table.code[(unsigned char)aligner->seq_a[seq_i]];
      bool is_match = match_row[code_a];
      int substitution_penalty = sub_row[code_a];

      if(scoring->no_mismatches && !is_match)
      {
        match_scores[index] = min;
      }
      else
      {
        // substitution
        // 1) continue alignment
        // 2) close gap in seq_a
        // 3) close gap in seq_b
        // max(x+s, y+s, z+s) == max(x,y,z)+s, with one add instead of three
        match_scores[index]
          = MAX2(MAX3(match_scores[index_upleft],
                      gap_a_scores[index_upleft],
                      gap_b_scores[index_upleft]) + substitution_penalty,
                 min);
      }

      // Long arithmetic since some INTs are set to min and penalty is -ve
      // (adding as ints would cause an integer overflow)

      // Update gap_a_scores[i][j] from position [i][j-1]
      if(seq_i == len_i-1 && scoring->no_end_gap_penalty)
      {
        gap_a_scores[index] = MAX3(match_scores[index_up],
                                   gap_a_scores[index_up],
                                   gap_b_scores[index_up]);
      }
      else if(!no_gaps_in_a || seq_i == len_i-1)
      {
        gap_a_scores[index]
          = MAX3(MAX2(match_scores[index_up],
                      gap_b_scores[index_up]) + gap_open_penalty,
                 gap_a_scores[index_up] + gap_extend_penalty,
                 min);
      }
      else
        gap_a_scores[index] = min;

      // Update gap_b_scores[i][j] from position [i-1][j]
      if(seq_j == len_j-1 && scoring->no_end_gap_penalty)
      {
        gap_b_scores[index] = MAX3(match_scores[index_left],
                                   gap_a_scores[index_left],
                                   gap_b_scores[index_left]);
      }
      else if(!no_gaps_in_b || seq_j == len_j-1)
      {
        gap_b_scores[index]
          = MAX3(MAX2(match_scores[index_left],
                      gap_a_scores[index_left]) + gap_open_penalty,
                 gap_b_scores[index_left] + gap_extend_penalty,
                 min);
      }
      else
        gap_b_scores[index] = min;

      index++;
      index_left++;
      index_up++;
      index_upleft++;
    }

    index++;
    index_left++;
    index_up++;
    index_upleft++;
  }

  sub_table_free(&table);
}

// Allocate memory for one score matrix. On Linux, large matrices are 2 MB
// aligned and marked for transparent huge pages: with 4 KB pages a 5 kb x
// 5 kb alignment takes ~72,000 page faults and many TLB misses, with huge
// pages ~320 (whole program 0.14 -> 0.08 s; repeated alignments ~30% faster)
static void* alignment_alloc_matrix(size_t mem)
{
#ifdef __linux__
  const size_t huge = (size_t)2 << 20;
  if(mem >= huge) {
    mem = (mem + huge - 1) / huge * huge; // aligned_alloc needs a multiple
    void *ptr = aligned_alloc(huge, mem);
    if(ptr != NULL) madvise(ptr, mem, MADV_HUGEPAGE); // only a hint
    return ptr;
  }
#endif
  return malloc(mem);
}

#ifdef ALN_STRIPED
// Room for the traceback bytes of `cells` cells
static void alignment_trace_capacity(aligner_t *aligner, size_t cells)
{
  if(aligner->trace_capacity < cells) {
    free(aligner->trace);
    aligner->trace_capacity = ROUNDUP2POW(cells);
    aligner->trace = alignment_alloc_matrix(aligner->trace_capacity);
    if(aligner->trace == NULL) {
      fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

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
// faster than the whole fill on its threads
static bool alignment_fill_banded(aligner_t *aligner, const sub_table_t *table,
                                  score_t min, score_t *const row0[3],
                                  score_t *const col0[3])
{
  const scoring_t *scoring = aligner->scoring;
  const size_t w = aligner->score_width, h = aligner->score_height;
  const size_t len_i = w-1, len_j = h-1;
  const ptrdiff_t D = (ptrdiff_t)len_i - (ptrdiff_t)len_j;
  // The whole fill's threads, or one for a narrow band
  const double full_threads = (double)alignment_fill_threads(len_i, len_j);
  size_t i, m, y0, y1, c0, c1;

  if(!alignment_band_enabled || scoring->gap_open > 0 ||
     scoring->gap_extend > 0) return false;

  score_t smax = 0;
  for(i = 0; i < table->n * table->n; i++) smax = MAX2(smax, table->scores[i]);

  size_t bw = alignment_band_min_width;
  score_t *top = NULL, *lcol = NULL;
  bool done = false;

  while(!done)
  {
    // Columns per row, plus the strip height the window grows by
    if(bw >= len_i) break;
    const double cols = (double)(D < 0 ? -D : D) + 2.0*bw + 1;
    const size_t strip = alignment_band_strip_rows ? alignment_band_strip_rows
                         : MAX2(256, (size_t)(cols / 4));
    // A fill on t threads runs about max(1, t/2) times faster (16 threads:
    // 6-8x); a narrow band runs on one
    const double band_threads =
      (double)alignment_fill_threads((size_t)cols + strip, strip);
    if((cols + strip) / MAX2(1.0, band_threads / 2) >=
       0.6 * (double)len_i / MAX2(1.0, full_threads / 2)) break;

    const ptrdiff_t dlo = MIN2(0, D) - (ptrdiff_t)bw;
    const ptrdiff_t dhi = MAX2(0, D) + (ptrdiff_t)bw;

    // Each row keeps the bytes of its strip's columns
    if(aligner->trace_rows_capacity < h) {
      free(aligner->trace_rows);
      aligner->trace_rows_capacity = h;
      aligner->trace_rows = malloc(h * sizeof(ptrdiff_t));
    }
    if(top == NULL) {
      top = malloc(3 * w * sizeof(score_t));
      lcol = malloc(3 * (strip + 1) * sizeof(score_t));
    }
    else lcol = realloc(lcol, 3 * (strip + 1) * sizeof(score_t));
    if(aligner->trace_rows == NULL || top == NULL || lcol == NULL) {
      fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    size_t cells = 0;
    for(y0 = 1; y0 <= len_j; y0 = y1 + 1) {
      alignment_strip(y0, strip, len_i, len_j, dlo, dhi, &y1, &c0, &c1);
      for(i = y0; i <= y1; i++, cells += c1 - c0)
        aligner->trace_rows[i] = (ptrdiff_t)cells - (ptrdiff_t)(c0 + 1);
    }
    alignment_trace_capacity(aligner, cells);

    // Row 0, then each strip's last row for the next
    score_t *topm[3] = {top, top + w, top + 2*w}, *lcm[3];
    for(m = 0; m < 3; m++) {
      memcpy(topm[m], row0[m], w * sizeof(score_t));
      lcm[m] = lcol + m * (strip + 1);
    }

    size_t jr_prev = len_i; // the top row is right up to here
    for(y0 = 1; y0 <= len_j; y0 = y1 + 1)
    {
      alignment_strip(y0, strip, len_i, len_j, dlo, dhi, &y1, &c0, &c1);
      // Columns the strip above did not reach
      for(i = jr_prev + 1; i <= c1; i++)
        topm[MATCH][i] = topm[GAP_A][i] = topm[GAP_B][i] = min;
      jr_prev = c1;

      alignment_block_t tmpl = {
        .aligner = aligner, .table = table, .min = min,
        .c0 = c0, .width = c1 - c0, .y0 = y0, .y1 = y1,
        .r0 = {topm[MATCH], topm[GAP_A], topm[GAP_B]}, .lstride = 1,
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
      alignment_fill_blocks(&tmpl);
    }

    score_t best = MAX3(aligner->final_scores[MATCH],
                        aligner->final_scores[GAP_A],
                        aligner->final_scores[GAP_B]);
    size_t need = alignment_band_width(smax, scoring->gap_extend,
                                       len_i, len_j, best);
    if(need <= bw) done = true;
    else bw = MAX2(need, 2*bw);
  }

  free(top);
  free(lcol);
  aligner->trace_banded = done;
  return done;
}
#endif

// Compact mode: run the striped fill keeping one traceback byte per cell
// instead of the score matrices. Returns false (nothing done) where the
// striped fill does not apply, with free start gaps (the traceback changes
// penalties along the first row and column), or where the first row and
// column fall below min: working out the traceback bytes adds a penalty to
// each predecessor, which could then overflow (the fill adds it only to the
// best one)
static bool alignment_fill_trace(aligner_t *aligner, char is_sw)
{
#ifdef ALN_STRIPED
  const scoring_t *scoring = aligner->scoring;
  const size_t w = aligner->score_width, h = aligner->score_height;
  const size_t len_i = w-1, len_j = h-1;
  size_t i, m;

  bool no_gaps_in_a, no_gaps_in_b;
  alignment_no_gaps(scoring, len_i, len_j, &no_gaps_in_a, &no_gaps_in_b);
  if(!alignment_striped_enabled || len_i == 0 || len_j == 0 ||
     len_i >= INT_MAX || len_j >= INT_MAX ||
     w + h > SIZE_MAX / (3 * sizeof(score_t)) ||
     scoring->no_start_gap_penalty ||
     scoring->no_end_gap_penalty || scoring->no_mismatches ||
     no_gaps_in_a || no_gaps_in_b) return false;

  const score_t min = alignment_min_score(scoring, is_sw);
  if(!is_sw) {
    // Gap scores along the first row and column run from these two
    long long edge1 = (long long)scoring->gap_open + scoring->gap_extend;
    long long edgen = (long long)scoring->gap_open +
                      (long long)MAX2(len_i, len_j) * scoring->gap_extend;
    if(MIN2(edge1, edgen) < min || MAX2(edge1, edgen) > INT_MAX) return false;
  }

  // Row 0 and column 0, as alignment_fill_matrices() sets them
  score_t *edges = malloc(3 * (w + h) * sizeof(score_t));
  if(edges == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  score_t *row0[3], *col0[3];
  for(m = 0; m < 3; m++) {
    row0[m] = edges + m * w;
    col0[m] = edges + 3 * w + m * h;
    row0[m][0] = col0[m][0] = 0;
  }
  for(i = 1; i < w; i++) {
    row0[MATCH][i] = row0[GAP_A][i] = is_sw ? 0 : min;
    row0[GAP_B][i] = is_sw ? 0 : scoring->gap_open + (int)i * scoring->gap_extend;
  }
  for(i = 1; i < h; i++) {
    col0[MATCH][i] = col0[GAP_B][i] = min;
    col0[GAP_A][i] = is_sw ? min : scoring->gap_open + (int)i * scoring->gap_extend;
  }

  sub_table_t table;
  sub_table_build(&table, scoring, aligner->seq_a, len_i, aligner->seq_b, len_j);
  aligner->trace_banded = false;
  if(is_sw || !alignment_fill_banded(aligner, &table, min, row0, col0))
  {
    alignment_trace_capacity(aligner, w * h);
    alignment_block_t tmpl = {
      .aligner = aligner, .table = &table, .min = min,
      .r0 = {row0[MATCH], row0[GAP_A], row0[GAP_B]},
      .lc = {col0[MATCH], col0[GAP_A], col0[GAP_B]}, .lstride = 1,
      .trace = aligner->trace, .track_best = is_sw,
      .width = len_i, .y0 = 1, .y1 = len_j};
    alignment_fill_blocks(&tmpl);
  }
  sub_table_free(&table);
  free(edges);
  return true;
#else
  (void)aligner; (void)is_sw;
  return false;
#endif
}

void aligner_align(aligner_t *aligner,
                   const char *seq_a, const char *seq_b,
                   size_t len_a, size_t len_b,
                   const scoring_t *scoring, char is_sw)
{
  aligner->scoring = scoring;
  aligner->seq_a = seq_a;
  aligner->seq_b = seq_b;
  aligner->score_width = len_a+1;
  aligner->score_height = len_b+1;

  if(aligner->score_height > SIZE_MAX / aligner->score_width / sizeof(score_t)) {
    fprintf(stderr, "%s:%i: Sequences too long to align\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  aligner->traced = aligner->compact && alignment_fill_trace(aligner, is_sw);
  if(aligner->traced) return;

  size_t new_capacity = aligner->score_width * aligner->score_height;

  if(aligner->capacity < new_capacity)
  {
    aligner->capacity = ROUNDUP2POW(new_capacity);
    size_t mem = sizeof(score_t) * aligner->capacity;
    // The fill overwrites every cell it uses, so old contents are not kept
    free(aligner->match_scores);
    free(aligner->gap_a_scores);
    free(aligner->gap_b_scores);
    aligner->match_scores = alignment_alloc_matrix(mem);
    aligner->gap_a_scores = alignment_alloc_matrix(mem);
    aligner->gap_b_scores = alignment_alloc_matrix(mem);
    if(!aligner->match_scores || !aligner->gap_a_scores || !aligner->gap_b_scores) {
      fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  alignment_fill_matrices(aligner, is_sw);
}

void aligner_destroy(aligner_t *aligner)
{
  free(aligner->trace);
  free(aligner->trace_rows);
  if(aligner->capacity > 0) {
    free(aligner->match_scores);
    free(aligner->gap_a_scores);
    free(aligner->gap_b_scores);
  }
}


alignment_t* alignment_create(size_t capacity)
{
  capacity = ROUNDUP2POW(capacity);
  alignment_t *result = malloc(sizeof(alignment_t));
  if(result == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  result->result_a = malloc(sizeof(char)*capacity);
  result->result_b = malloc(sizeof(char)*capacity);
  if(result->result_a == NULL || result->result_b == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  result->capacity = capacity;
  result->length = 0;
  result->result_a[0] = result->result_b[0] = '\0';
  result->pos_a = result->pos_b = result->len_a = result->len_b = 0;
  result->score = 0;
  return result;
}

void alignment_ensure_capacity(alignment_t* result, size_t strlength)
{
  size_t capacity = strlength+1;
  if(result->capacity < capacity)
  {
    capacity = ROUNDUP2POW(capacity);
    result->result_a = realloc(result->result_a, sizeof(char)*capacity);
    result->result_b = realloc(result->result_b, sizeof(char)*capacity);
    result->capacity = capacity;
    if(result->result_a == NULL || result->result_b == NULL) {
      fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

void alignment_free(alignment_t* result)
{
  free(result->result_a);
  free(result->result_b);
  free(result);
}


// Backtrack through scoring matrices
void alignment_reverse_move(enum Matrix *curr_matrix, score_t *curr_score,
                            size_t *score_x, size_t *score_y,
                            size_t *arr_index, const aligner_t *aligner)
{
  size_t seq_x = (*score_x)-1, seq_y = (*score_y)-1;
  size_t len_i = aligner->score_width-1, len_j = aligner->score_height-1;

  bool is_match;
  int match_penalty;
  const scoring_t *scoring = aligner->scoring;

  bool no_gaps_in_a, no_gaps_in_b;
  alignment_no_gaps(scoring, len_i, len_j, &no_gaps_in_a, &no_gaps_in_b);

  scoring_lookup(scoring, aligner->seq_a[seq_x], aligner->seq_b[seq_y],
                 &match_penalty, &is_match);

  int gap_a_open_penalty, gap_b_open_penalty;
  int gap_a_extend_penalty, gap_b_extend_penalty;

  gap_a_open_penalty = gap_b_open_penalty = scoring->gap_extend + scoring->gap_open;
  gap_a_extend_penalty = gap_b_extend_penalty = scoring->gap_extend;

  // Free gaps at the ends
  if(scoring->no_end_gap_penalty) {
    if(*score_x == len_i) gap_a_open_penalty = gap_a_extend_penalty = 0;
    if(*score_y == len_j) gap_b_open_penalty = gap_b_extend_penalty = 0;
  }
  if(scoring->no_start_gap_penalty) {
    if(*score_x == 0) gap_a_open_penalty = gap_a_extend_penalty = 0;
    if(*score_y == 0) gap_b_open_penalty = gap_b_extend_penalty = 0;
  }

  int prev_match_penalty, prev_gap_a_penalty, prev_gap_b_penalty;

  switch(*curr_matrix)
  {
    case MATCH:
      prev_match_penalty = match_penalty;
      prev_gap_a_penalty = match_penalty;
      prev_gap_b_penalty = match_penalty;
      (*score_x)--;
      (*score_y)--;
      (*arr_index) -= aligner->score_width + 1;
      break;

    case GAP_A:
      prev_match_penalty = gap_a_open_penalty;
      prev_gap_a_penalty = gap_a_extend_penalty;
      prev_gap_b_penalty = gap_a_open_penalty;
      (*score_y)--;
      (*arr_index) -= aligner->score_width;
      break;

    case GAP_B:
      prev_match_penalty = gap_b_open_penalty;
      prev_gap_a_penalty = gap_b_open_penalty;
      prev_gap_b_penalty = gap_b_extend_penalty;
      (*score_x)--;
      (*arr_index)--;
      break;

    default:
      fprintf(stderr, "Program error: invalid matrix in get_reverse_move()\n");
      fprintf(stderr, "Please submit a bug report to: https://github.com/DanielTillett/nucleics-align/issues\n");
      exit(EXIT_FAILURE);
  }

  // *arr_index = ARR_2D_INDEX(aligner->score_width, *score_x, *score_y);

  if((!no_gaps_in_a || *score_x == 0 || *score_x == len_i) &&
     aligner->gap_a_scores[*arr_index] + prev_gap_a_penalty == *curr_score)
  {
    *curr_matrix = GAP_A;
    *curr_score = aligner->gap_a_scores[*arr_index];
  }
  else if((!no_gaps_in_b || *score_y == 0 || *score_y == len_j) &&
          aligner->gap_b_scores[*arr_index] + prev_gap_b_penalty == *curr_score)
  {
    *curr_matrix = GAP_B;
    *curr_score = aligner->gap_b_scores[*arr_index];
  }
  else if(aligner->match_scores[*arr_index] + prev_match_penalty == *curr_score)
  {
    *curr_matrix = MATCH;
    *curr_score = aligner->match_scores[*arr_index];
  }
  else
  {
    alignment_print_matrices(aligner);

    fprintf(stderr, "[%s:%zu,%zu]: %i [ismatch: %i] '%c' '%c'\n",
            MATRIX_NAME(*curr_matrix), *score_x, *score_y, *curr_score,
            is_match, aligner->seq_a[seq_x], aligner->seq_b[seq_y]);
    fprintf(stderr, " Penalties match: %i gap_open: %i gap_extend: %i\n",
            prev_match_penalty, prev_gap_a_penalty, prev_gap_b_penalty);
    fprintf(stderr, " Expected MATCH: %i GAP_A: %i GAP_B: %i\n",
            aligner->match_scores[*arr_index],
            aligner->gap_a_scores[*arr_index],
            aligner->gap_b_scores[*arr_index]);

    fprintf(stderr,
"Program error: traceback fail (get_reverse_move)\n"
"This may be due to an integer overflow if your sequences are long or scores\n"
"are large. If this is the case using smaller scores or shorter sequences may\n"
"work around this problem.  \n"
"  If you think this is a bug, please report it to: https://github.com/DanielTillett/nucleics-align/issues\n");
    exit(EXIT_FAILURE);
  }
}


void alignment_print_matrices(const aligner_t *aligner)
{
  const score_t* match_scores = aligner->match_scores;
  const score_t* gap_a_scores = aligner->gap_a_scores;
  const score_t* gap_b_scores = aligner->gap_b_scores;

  size_t i, j;

  printf("seq_a: %.*s\nseq_b: %.*s\n",
         (int)aligner->score_width-1, aligner->seq_a,
         (int)aligner->score_height-1, aligner->seq_b);

  printf("match_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    printf("%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      printf("\t%3i", (int)ARR_LOOKUP(match_scores, aligner->score_width, i, j));
    }
    putc('\n', stdout);
  }
  printf("gap_a_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    printf("%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      printf("\t%3i", (int)ARR_LOOKUP(gap_a_scores, aligner->score_width, i, j));
    }
    putc('\n', stdout);
  }
  printf("gap_b_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    printf("%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      printf("\t%3i", (int)ARR_LOOKUP(gap_b_scores, aligner->score_width, i, j));
    }
    putc('\n', stdout);
  }

  printf("match: %i mismatch: %i gapopen: %i gapexend: %i\n",
         aligner->scoring->match, aligner->scoring->mismatch,
         aligner->scoring->gap_open, aligner->scoring->gap_extend);
  printf("\n");
}

void alignment_colour_fprint_against(FILE *out, const char *alignment_a,
                                     const char *alignment_b,
                                     char case_sensitive)
{
  size_t i;
  char red = 0, green = 0;

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_b[i] == '-')
    {
      if(!red)
      {
        fputs(align_col_indel, out);
        red = 1;
      }
    }
    else if(red)
    {
      red = 0;
      fputs(align_col_stop, out);
    }

    if(((case_sensitive && alignment_a[i] != alignment_b[i]) ||
        (!case_sensitive && tolower((unsigned char)alignment_a[i]) !=
                            tolower((unsigned char)alignment_b[i]))) &&
       alignment_a[i] != '-' && alignment_b[i] != '-')
    {
      if(!green)
      {
        fputs(align_col_mismatch, out);
        green = 1;
      }
    }
    else if(green)
    {
      green = 0;
      fputs(align_col_stop, out);
    }

    putc(alignment_a[i], out);
  }

  if(green || red)
  {
    // Stop all colours
    fputs(align_col_stop, out);
  }
}

// Order of alignment_a / alignment_b is not important
void alignment_fprint_spacer(FILE *out, const char* alignment_a,
                             const char* alignment_b,
                             const scoring_t* scoring)
{
  size_t i;

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '-' || alignment_b[i] == '-')
    {
      putc(' ', out);
    }
    else if(alignment_a[i] == alignment_b[i] ||
            (!scoring->case_sensitive &&
             tolower((unsigned char)alignment_a[i]) ==
             tolower((unsigned char)alignment_b[i])))
    {
      putc('|', out);
    }
    else
    {
      putc('*', out);
    }
  }
}

void alignment_colour_print_against(const char *alignment_a,
                                    const char *alignment_b,
                                    char case_sensitive)
{
  alignment_colour_fprint_against(stdout, alignment_a, alignment_b,
                                  case_sensitive);
}

void alignment_print_spacer(const char* alignment_a, const char* alignment_b,
                            const scoring_t* scoring)
{
  alignment_fprint_spacer(stdout, alignment_a, alignment_b, scoring);
}
