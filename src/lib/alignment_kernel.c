/*
 alignment_kernel.c
 Striped SIMD fill kernels, column-block threading, and Smith-Waterman scans.
 */

#define _DEFAULT_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <pthread.h>
#include <sched.h>
#include <stdatomic.h>
#include <unistd.h>

#include "alignment_fill_internal.h"
#include "alignment_macros.h"
#include "alignment_alloc.h"
#include "alignment_testing.h"

#ifdef ALN_STRIPED

// Publish progress every ALN_PUB_ROWS rows (plus the final row) to cut
// cross-core atomic traffic while keeping the pipeline full
#define ALN_PUB_ROWS 16

static inline void aln_publish_row(aln_done_t *d, size_t y)
{
  atomic_store_explicit(&d->v, y, memory_order_release);
  atomic_thread_fence(memory_order_seq_cst);
  if(atomic_load_explicit(&d->waiters, memory_order_relaxed) > 0) {
    pthread_mutex_lock(&d->mu);
    pthread_cond_broadcast(&d->cv);
    pthread_mutex_unlock(&d->mu);
  }
}

static inline void aln_wake_right(aln_done_t *d)
{
  if(d == NULL) return;
  // Publish SIZE_MAX so a late waiter registering after an abort never parks
  atomic_store_explicit(&d->v, SIZE_MAX, memory_order_release);
  atomic_thread_fence(memory_order_seq_cst);
  if(atomic_load_explicit(&d->waiters, memory_order_relaxed) > 0) {
    pthread_mutex_lock(&d->mu);
    pthread_cond_broadcast(&d->cv);
    pthread_mutex_unlock(&d->mu);
  }
}

// Wait until the block on the left has finished row y. Caches the last
// observed row in *seen so the common case needs no atomic load. Uses a
// bounded spin followed by a condition-variable wait to cut wasted CPU time.
// Returns false if a block stopped meanwhile (see alignment_block_t.stop)
static inline bool alignment_wait_row(aln_done_t *left_done, size_t y,
                                      atomic_bool *stop, size_t *seen)
{
  if(*seen >= y) return true;
  unsigned spins = 0;
  size_t cur;
  while((cur = atomic_load_explicit(&left_done->v, memory_order_acquire)) < y) {
    if(stop != NULL &&
       atomic_load_explicit(stop, memory_order_relaxed)) return false;
    if(spins < 512) {
#if defined(__x86_64__) || defined(__i386__)
      __builtin_ia32_pause();
#endif
      spins++;
      continue;
    }
    pthread_mutex_lock(&left_done->mu);
    atomic_fetch_add_explicit(&left_done->waiters, 1, memory_order_seq_cst);
    atomic_thread_fence(memory_order_seq_cst);
    while((cur = atomic_load_explicit(&left_done->v, memory_order_acquire)) < y &&
          (stop == NULL || !atomic_load_explicit(stop, memory_order_acquire))) {
      pthread_cond_wait(&left_done->cv, &left_done->mu);
    }
    atomic_fetch_sub_explicit(&left_done->waiters, 1, memory_order_relaxed);
    pthread_mutex_unlock(&left_done->mu);
    break;
  }
  if(stop != NULL && atomic_load_explicit(stop, memory_order_acquire))
    return false;
  if(cur < y) return false;
  *seen = cur;
  return true;
}

#define ALN_K_SUFFIX 32
#define ALN_K_S score_t
#define ALN_K_VL ALN_VL
#define ALN_K_NARROW 0
#define ALN_K_BYTE 0
#include "alignment_kernel.inc"

#define ALN_K_SUFFIX 16
#define ALN_K_S int16_t
#define ALN_K_VL (2 * ALN_VL)
#define ALN_K_NARROW 1
#define ALN_K_BYTE 0
#include "alignment_kernel.inc"

#define ALN_K_SUFFIX 8
#define ALN_K_S int8_t
#define ALN_K_VL (4 * ALN_VL)
#define ALN_K_NARROW 1
#define ALN_K_BYTE 1
#include "alignment_kernel.inc"

#define ALN_D_SUFFIX 8
#define ALN_D_S int8_t
#define ALN_D_VL (4 * ALN_VL)
#include "alignment_diff.inc"

#define ALN_D_SUFFIX 16
#define ALN_D_S int16_t
#define ALN_D_VL (2 * ALN_VL)
#include "alignment_diff.inc"

// Lanes for a difference-encoded Needleman-Wunsch fill with this table:
// 8 or 16 bits, or 0 where neither applies, with the floor that stands for
// -infinity in *floor and the most a cell lies above the one above it in
// *ceiling. With p the gap open penalty (gap_open + gap_extend)
// and e the gap extension penalty, both negated (gaps must not score above
// 0), and s+ and s- the highest and the most negative substitution scores
// (at least 0): a cell and the cell above or beside it differ by at most
// s+ + p, so every real difference lies within K = s- + 2 s+ + 5 p + e of
// 0. The floor must stay below every real value after one step adds to it
// (at most s+ + p), and no sum may leave the lane type: the most taken off
// in one go is e + s+ + p
int alignment_diff_bits(const sub_table_t *table,
                               const scoring_t *scoring, score_t *floor,
                               score_t *ceiling)
{
  if(scoring->gap_open > 0 || scoring->gap_extend > 0) return 0;
  const long long p = -((long long)scoring->gap_open + scoring->gap_extend);
  const long long e = -(long long)scoring->gap_extend;
  long long smax = 0, smin = 0;
  size_t i;
  for(i = 0; i < table->n * table->n; i++) {
    smax = MAX2(smax, (long long)table->scores[i]);
    smin = MAX2(smin, -(long long)table->scores[i]);
  }
  const long long K = smin + 2*smax + 5*p + e, step = smax + p;
  *ceiling = (score_t)step; // no cell is more than this above the one above
  if(K + step < 100 && e + step <= 28) { *floor = -100; return 8; }
  if(K + step < 16000 && e + step <= 16000) { *floor = -16000; return 16; }
  return 0;
}

// Whether a Smith-Waterman fill of len_i columns with this table may use
// lanes of `bits` bits (8 or 16), and the highest row best (*limit) that
// still keeps the next row in range. Cells are then never below 0 (min) nor
// above the best match score so far: no gap scores above 0, so gap cells
// take at most the score of a match cell before them, and a match cell at
// most the best of the row above plus the highest substitution score. Sums
// of a cell (at least 0) and a penalty stay in range as every penalty is
// at least the lane type's lowest value (-16000 for 16 bits, so a sum of
// two stays in range too). 16-bit lanes also hold segment numbers
bool alignment_narrow_limit(const sub_table_t *table,
                                   const scoring_t *scoring, score_t min,
                                   size_t len_i, int bits, score_t *limit)
{
  const long long open = (long long)scoring->gap_open + scoring->gap_extend;
  const long long ext = scoring->gap_extend;
  const long long lo = bits == 8 ? INT8_MIN : -16000;
  const score_t hi = bits == 8 ? INT8_MAX : INT16_MAX;
  size_t i;
  if(min != 0 || open > 0 || ext > 0 || open < lo || ext < lo ||
     (bits == 16 && len_i / (2 * ALN_VL) + 1 >= INT16_MAX)) return false;
  score_t smax = 0;
  for(i = 0; i < table->n * table->n; i++) {
    if(table->scores[i] < lo || table->scores[i] > -lo) return false;
    smax = MAX2(smax, table->scores[i]);
  }
  *limit = hi - smax;
  return *limit > 0;
}

// The narrowest lanes a Smith-Waterman fill of len_i columns may start
// with (tmpl->bits, and narrow_limit), allowed by opts: 8 bits, then 16,
// else 32 (bits 0). A target scan needs the target in range
void alignment_lanes(alignment_block_t *tmpl,
                            const alignment_options_t *opts,
                            const scoring_t *scoring, size_t len_i)
{
  const sub_table_t *table = tmpl->table;
  tmpl->bits = 0;
  if(opts->byte &&
     alignment_narrow_limit(table, scoring, tmpl->min, len_i, 8,
                            &tmpl->narrow_limit) &&
     tmpl->target <= tmpl->narrow_limit) tmpl->bits = 8;
  else if(opts->narrow &&
          alignment_narrow_limit(table, scoring, tmpl->min, len_i, 16,
                                 &tmpl->narrow_limit) &&
          tmpl->target <= tmpl->narrow_limit) tmpl->bits = 16;
}

static void* alignment_block_thread(void *ptr)
{
  alignment_block_t *p = ptr;
  if(p->diff == 8) alignment_fill_diff_8(p);
  else if(p->diff == 16) alignment_fill_diff_16(p);
  else if(p->bits == 8) alignment_fill_striped_8(p);
  else if(p->bits == 16) alignment_fill_striped_16(p);
  else alignment_fill_striped_32(p);
  // Always wake any downstream sleeper on normal or early (OOM/overflow/stop) exit
  aln_wake_right(p->done);
  return NULL;
}

static int alignment_online_cpus(void)
{
  static atomic_int cached = 0;
  int c = atomic_load_explicit(&cached, memory_order_relaxed);
  if(c <= 0) {
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    c = n > 0 ? (int)n : 1;
    atomic_store_explicit(&cached, c, memory_order_relaxed);
  }
  return c;
}

// Threads (column blocks) for a striped fill of width x rows cells
// A fill across the whole width (whole: every fill but a band's strips)
// takes blocks of at least twice block_min_cols columns: thinner blocks
// finish a row too soon after the block on their left, so threads mostly
// wait (20 kb: 16 threads no faster than 10, with 0.79 -> 0.44 s of CPU).
// A band's strips take their own, smaller minimums (band_block_min_cols,
// band_block_min_cells): narrowed strips have few columns and few rows,
// and would otherwise run on one thread (1 Mb at 1%: 8.5 -> 2.4 s)
size_t alignment_fill_threads(const alignment_options_t *o,
                                    size_t width, size_t rows, bool whole)
{
  int threads = o->max_threads;
  if(threads <= 0) {
    int cpus = alignment_online_cpus();
    threads = cpus > 1 ? cpus / 2 : 1;
  }
  if(threads > ALN_MAX_THREADS) threads = ALN_MAX_THREADS;
  // At least one vector of columns per block, so no block is empty
  const size_t min_cols = MAX2(whole ? 2 * o->block_min_cols
                                     : o->band_block_min_cols, (size_t)ALN_VL);
  const size_t min_cells = whole ? o->block_min_cells
                                 : o->band_block_min_cells;
  if(threads > 1 && (double)width * rows >= (double)min_cells &&
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
// Returns false if out of memory, with nothing filled
bool alignment_fill_blocks(alignment_block_t *tmpl)
{
  aligner_t *aligner = tmpl->aligner;
  const size_t w = aligner->score_width, h = aligner->score_height;
  const size_t len_i = tmpl->width, len_j = tmpl->y1 - tmpl->y0 + 1;
  const size_t nblocks = alignment_fill_threads(&aligner->opts, len_i, len_j,
                                                !tmpl->band);
  size_t i, m;

  // Blocks of equal width (a multiple of ALN_VL); the last takes the rest
  size_t width = (len_i / nblocks) / ALN_VL * ALN_VL;

  // Per-block state on the heap: thousands of blocks would not fit a small
  // thread stack (macOS gives threads 512 KB)
  alignment_block_t *args = aln_malloc(nblocks * sizeof(*args));
  aln_done_t *done = aln_malloc(nblocks * sizeof(*done));
  pthread_t *thread_ids = aln_malloc(nblocks * sizeof(*thread_ids));
  bool *started = aln_malloc(nblocks * sizeof(*started));
  if(args == NULL || done == NULL || thread_ids == NULL || started == NULL) {
    free(args); free(done); free(thread_ids); free(started);
    return false;
  }

  atomic_bool stop, nomem, gave_up;
  atomic_init(&stop, false);
  atomic_init(&nomem, false);
  atomic_init(&gave_up, false);
  if(tmpl->bits == 8) atomic_fetch_add(&aln_stat_byte_fills, 1);
  if(tmpl->bits == 16) atomic_fetch_add(&aln_stat_narrow_fills, 1);
  if(tmpl->diff != 0) atomic_fetch_add(&aln_stat_diff_fills, 1);

  // With trace, each block's last column goes here for the next block:
  // the rows of this fill (a band's strip, or the whole matrix) and the
  // row above them, not the whole matrix height (a 1 Mb pair's strips
  // would each map 12 MB per block for a few hundred rows)
  score_t *cols = NULL;
  const size_t col_rows = len_j + 1; // rows y0-1 .. y1
  if(tmpl->trace != NULL || tmpl->score_only) {
    if(col_rows > SIZE_MAX / sizeof(score_t) / 3 / nblocks) cols = NULL;
    else cols = aln_malloc(nblocks * 3 * col_rows * sizeof(score_t));
    if(cols == NULL) {
      free(args); free(done); free(thread_ids); free(started);
      return false;
    }
  }

  for(i = 0; i < nblocks; i++) {
    atomic_init(&done[i].v, 0);
    atomic_init(&done[i].waiters, 0);
    if(pthread_mutex_init(&done[i].mu, NULL) != 0) {
      for(size_t j = 0; j < i; j++) {
        pthread_cond_destroy(&done[j].cv);
        pthread_mutex_destroy(&done[j].mu);
      }
      free(cols);
      free(args); free(done); free(thread_ids); free(started);
      return false;
    }
    if(pthread_cond_init(&done[i].cv, NULL) != 0) {
      pthread_mutex_destroy(&done[i].mu);
      for(size_t j = 0; j < i; j++) {
        pthread_cond_destroy(&done[j].cv);
        pthread_mutex_destroy(&done[j].mu);
      }
      free(cols);
      free(args); free(done); free(thread_ids); free(started);
      return false;
    }
    alignment_block_t *b = &args[i];
    *b = *tmpl;
    b->c0 = tmpl->c0 + i * width;
    b->width = i + 1 < nblocks ? width : len_i - i * width;
    b->trace_boff = tmpl->diff != 0
                    ? i * (aligner->trace_seg0 * aligner->trace_vl)
                    : 0;
    b->left_done = i > 0 ? &done[i-1] : NULL;
    b->done = i + 1 < nblocks ? &done[i] : NULL;
    b->best_score = 0;
    b->best_x = b->best_y = 0;
    b->stop = &stop;
    b->nomem = &nomem;
    b->gave_up = &gave_up;
    b->rc_y0 = tmpl->y0 - 1;
    for(m = 0; m < 3; m++) {
      if(cols != NULL) {
        b->rc[m] = cols + (i*3 + m) * col_rows; // by row, from row y0-1
        if(i > 0) {
          b->lc[m] = args[i-1].rc[m];
          b->lstride = 1;
          b->ly0 = tmpl->y0 - 1;
        }
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
  for(i = 0; i < nblocks; i++) {
    pthread_cond_destroy(&done[i].cv);
    pthread_mutex_destroy(&done[i].mu);
  }
  free(done);
  free(thread_ids);
  free(started);

  if(atomic_load(&nomem)) {
    free(cols);
    free(args);
    return false;
  }
  tmpl->given_up = atomic_load(&gave_up);
  if(tmpl->given_up) {
    free(cols);
    free(args);
    return true;
  }

  // Scores grew past 8- or 16-bit lanes: all again with wider ones (16 bits
  // after 8 where they are allowed, else 32)
  if(tmpl->bits != 0 && atomic_load(&stop)) {
    free(cols);
    free(args);
    if(tmpl->bits == 8) {
      atomic_fetch_add(&aln_stat_byte_retries, 1);
      const scoring_t *scoring = aligner->scoring;
      tmpl->bits = aligner->opts.narrow &&
                   alignment_narrow_limit(tmpl->table, scoring, tmpl->min,
                                          len_i, 16, &tmpl->narrow_limit) &&
                   tmpl->target <= tmpl->narrow_limit ? 16 : 0;
    }
    else {
      atomic_fetch_add(&aln_stat_narrow_retries, 1);
      tmpl->bits = 0;
    }
    return alignment_fill_blocks(tmpl);
  }

  if(cols != NULL) {
    if(tmpl->y1 == h-1 && tmpl->c0 + tmpl->width == w-1)
      for(m = 0; m < 3; m++)
        aligner->final_scores[m] = args[nblocks-1].rc[m][h-1 - (tmpl->y0-1)];
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
  free(args);
  return true;
}

// Smith-Waterman scans of a against b in O(len_b) memory, for the best hit
// of smith_waterman.c. With target 0: the best match score and its first
// cell (lowest x, then y; *x, *y count from 1), as a full fill would rank
// it. With target > 0, for a scan in which no cell scores more than target:
// the highest x and the highest y of cells scoring target (0 if none).
// Returns false, with nothing found, where the striped fill does not apply
// or (target > 0) a cell scores more than target, or if out of memory.
// opts may be NULL (the defaults). short_hit: see alignment_sw_scan_short();
// ckpt, every: see alignment_sw_rows()
static bool alignment_sw_scan_any(const char *a, size_t len_a,
                                  const char *b, size_t len_b,
                                  const scoring_t *scoring, score_t target,
                                  score_t min_expected,
                                  const alignment_options_t *opts,
                                  bool short_hit, score_t *ckpt, size_t every,
                                  score_t *best, size_t *x, size_t *y)
{
  aligner_t tmp;
  aligner_init(&tmp);
  if(opts != NULL) tmp.opts = *opts;
  else alignment_options_default(&tmp.opts);

  bool no_gaps_in_a, no_gaps_in_b;
  alignment_no_gaps(scoring, len_a, len_b, &no_gaps_in_a, &no_gaps_in_b);
  if(!tmp.opts.striped || len_a == 0 || len_b == 0 ||
     len_a >= INT_MAX || len_b >= INT_MAX ||
     len_a + len_b > SIZE_MAX / (3 * sizeof(score_t)) ||
     scoring->no_start_gap_penalty || scoring->no_end_gap_penalty ||
     scoring->no_mismatches || no_gaps_in_a || no_gaps_in_b) return false;
  if(alignment_scores_may_overflow(scoring, a, len_a, b, len_b, true))
    return false; // aligner_align() reports it

  tmp.scoring = scoring;
  tmp.seq_a = a;
  tmp.seq_b = b;
  tmp.score_width = len_a + 1;
  tmp.score_height = len_b + 1;

  // Row 0 and column 0 are all 0 in Smith-Waterman: row 0 is read across
  // the columns, column 0 through one cell (stride 0)
  score_t *zeros = aln_calloc(len_a + 1, sizeof(score_t));
  sub_table_t table;
  if(zeros == NULL || !sub_table_build(&table, scoring, a, len_a, b, len_b)) {
    free(zeros);
    return false;
  }
  alignment_block_t tmpl = {
    .aligner = &tmp, .table = &table, .min = 0,
    .r0 = {zeros, zeros, zeros}, .lc = {zeros, zeros, zeros}, .lstride = 0,
    .width = len_a, .y0 = 1, .y1 = len_b,
    .score_only = true, .track_best = target == 0, .target = target,
    .ckpt = ckpt, .ckpt_every = every};
  alignment_lanes(&tmpl, &tmp.opts, scoring, len_a);
  if(min_expected > 0 && tmpl.bits != 0 && min_expected > tmpl.narrow_limit) {
    if(tmpl.bits == 8 && tmp.opts.narrow &&
       alignment_narrow_limit(&table, scoring, tmpl.min, len_a, 16,
                              &tmpl.narrow_limit) &&
       MAX2(tmpl.target, min_expected) <= tmpl.narrow_limit)
      tmpl.bits = 16;
    else
      tmpl.bits = 0;
  }
  // Row 0 is all 0 in Smith-Waterman (the fill writes the others; column 0
  // is never read back, see alignment_sw_trace_rows())
  if(ckpt != NULL) memset(ckpt, 0, 3 * (len_a + 1) * sizeof(score_t));
  if(short_hit) {
    // A hit scoring at least half the best substitution score per row of
    // b it has crossed (see _sw_align_rect() in smith_waterman.c)
    score_t smax = 0;
    size_t i;
    for(i = 0; i < table.n * table.n; i++) smax = MAX2(smax, table.scores[i]);
    tmpl.giveup_rate = 0.5 * MAX2(smax, 1);
  }
  bool filled = alignment_fill_blocks(&tmpl);
  sub_table_free(&table);
  free(zeros);
  if(!filled || tmpl.given_up) return false;

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

bool alignment_sw_scan(const char *a, size_t len_a, const char *b, size_t len_b,
                       const scoring_t *scoring, score_t target,
                       const alignment_options_t *opts,
                       score_t *best, size_t *x, size_t *y)
{
  return alignment_sw_scan_any(a, len_a, b, len_b, scoring, target, 0, opts,
                               false, NULL, 0, best, x, y);
}

bool alignment_sw_scan_hint(const char *a, size_t len_a,
                            const char *b, size_t len_b,
                            const scoring_t *scoring, score_t target,
                            score_t min_expected,
                            const alignment_options_t *opts,
                            score_t *best, size_t *x, size_t *y)
{
  return alignment_sw_scan_any(a, len_a, b, len_b, scoring, target,
                               min_expected, opts, false, NULL, 0, best, x, y);
}

bool alignment_sw_scan_short(const char *a, size_t len_a,
                             const char *b, size_t len_b,
                             const scoring_t *scoring,
                             const alignment_options_t *opts,
                             score_t *best, size_t *x, size_t *y)
{
  return alignment_sw_scan_any(a, len_a, b, len_b, scoring, 0, 0, opts, true,
                               NULL, 0, best, x, y);
}

bool alignment_sw_rows(const char *a, size_t len_a,
                       const char *b, size_t len_b,
                       const scoring_t *scoring,
                       const alignment_options_t *opts,
                       size_t every, score_t *ckpt,
                       score_t *best, size_t *x, size_t *y)
{
  if(every == 0) return false;
  return alignment_sw_scan_any(a, len_a, b, len_b, scoring, 0, 0, opts, false,
                               ckpt, every, best, x, y);
}

bool alignment_sw_trace_rows(const char *a, size_t len_a,
                             const char *b, size_t len_b,
                             const scoring_t *scoring,
                             const sub_table_t *table,
                             const alignment_options_t *opts,
                             score_t max_score,
                             const score_t *row, size_t y0, size_t y1,
                             size_t xmax, uint8_t *trace)
{
  aligner_t tmp;
  aligner_init(&tmp);
  if(opts != NULL) tmp.opts = *opts;
  else alignment_options_default(&tmp.opts);
  if(!tmp.opts.striped || xmax == 0 || xmax > len_a || y0 == 0 ||
     y0 > y1 || y1 > len_b) return false;

  tmp.scoring = scoring;
  tmp.seq_a = a;
  tmp.seq_b = b;
  tmp.score_width = len_a + 1;
  tmp.score_height = y1 + 1; // the fill's columns hold rows up to y1

  // Column 0 is all 0 in Smith-Waterman. Cell (x, y) of the strip at
  // trace[(y - y0) * xmax + x - 1]
  const size_t w = len_a + 1;
  score_t zero = 0;
  ptrdiff_t *trace_rows = aln_malloc((y1 + 1) * sizeof(ptrdiff_t));
  if(trace_rows == NULL) return false;
  sub_table_t local_table;
  const sub_table_t *use_table = table;
  if(use_table == NULL) {
    if(!sub_table_build(&local_table, scoring, a, len_a, b, len_b)) {
      free(trace_rows);
      return false;
    }
    use_table = &local_table;
  }
  size_t r;
  for(r = y0; r <= y1; r++)
    trace_rows[r] = (ptrdiff_t)((r - y0) * xmax) - 1;
  alignment_block_t tmpl = {
    .aligner = &tmp, .table = use_table, .min = 0,
    .r0 = {row, row + w, row + 2*w}, .lc = {&zero, &zero, &zero},
    .lstride = 0, .width = xmax, .y0 = y0, .y1 = y1,
    .trace = trace, .trace_rows = trace_rows};
  // Use 8- or 16-bit lanes when the known global max_score fits safely,
  // each where the options allow it (as alignment_lanes() chooses them)
  score_t lim = 0;
  if(max_score > 0) {
    if(tmp.opts.byte &&
       alignment_narrow_limit(use_table, scoring, 0, xmax, 8, &lim) &&
       max_score <= lim) {
      tmpl.bits = 8;
      tmpl.narrow_limit = lim;
    } else if(tmp.opts.narrow &&
              alignment_narrow_limit(use_table, scoring, 0, xmax, 16, &lim) &&
              max_score <= lim) {
      tmpl.bits = 16;
      tmpl.narrow_limit = lim;
    }
  }
  bool filled = alignment_fill_blocks(&tmpl);
  if(table == NULL) sub_table_free(&local_table);
  free(trace_rows);
  return filled;
}

#else

bool alignment_sw_scan(const char *a, size_t len_a, const char *b, size_t len_b,
                       const scoring_t *scoring, score_t target,
                       const alignment_options_t *opts,
                       score_t *best, size_t *x, size_t *y)
{
  (void)a; (void)len_a; (void)b; (void)len_b; (void)scoring; (void)target;
  (void)opts;
  (void)best; (void)x; (void)y;
  return false;
}

bool alignment_sw_scan_hint(const char *a, size_t len_a,
                            const char *b, size_t len_b,
                            const scoring_t *scoring, score_t target,
                            score_t max_hint, const alignment_options_t *opts,
                            score_t *best, size_t *x, size_t *y)
{
  (void)max_hint;
  return alignment_sw_scan(a, len_a, b, len_b, scoring, target, opts, best, x, y);
}

bool alignment_sw_scan_short(const char *a, size_t len_a,
                             const char *b, size_t len_b,
                             const scoring_t *scoring,
                             const alignment_options_t *opts,
                             score_t *best, size_t *x, size_t *y)
{
  return alignment_sw_scan(a, len_a, b, len_b, scoring, 0, opts, best, x, y);
}

bool alignment_sw_rows(const char *a, size_t len_a,
                       const char *b, size_t len_b,
                       const scoring_t *scoring,
                       const alignment_options_t *opts,
                       size_t every, score_t *ckpt,
                       score_t *best, size_t *x, size_t *y)
{
  (void)a; (void)len_a; (void)b; (void)len_b; (void)scoring; (void)opts;
  (void)every; (void)ckpt; (void)best; (void)x; (void)y;
  return false;
}

bool alignment_sw_trace_rows(const char *a, size_t len_a,
                             const char *b, size_t len_b,
                             const scoring_t *scoring,
                             const sub_table_t *table,
                             const alignment_options_t *opts,
                             score_t max_score,
                             const score_t *row, size_t y0, size_t y1,
                             size_t xmax, uint8_t *trace)
{
  (void)a; (void)len_a; (void)b; (void)len_b; (void)scoring; (void)table;
  (void)opts; (void)max_score;
  (void)row; (void)y0; (void)y1; (void)xmax; (void)trace;
  return false;
}

#endif /* ALN_STRIPED */
