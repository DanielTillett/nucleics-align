/*
 alignment.c
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Feb 2015
 */

// Turn on debugging output by defining NUCLEICS_ALIGN_VERBOSE
//#define NUCLEICS_ALIGN_VERBOSE

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
#include "hirschberg.h" // defaults for alignment_options_t
#include "smith_waterman.h"
#include "alignment_testing.h"
#include "alignment_alloc.h"

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

bool sub_table_build(sub_table_t *t, const scoring_t *scoring,
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
  t->scores = aln_malloc(cells * sizeof(int));
  t->is_match = aln_malloc(cells * sizeof(bool));
  if(t->scores == NULL || t->is_match == NULL) {
    sub_table_free(t);
    return false;
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
  return true;
}

void sub_table_free(sub_table_t *t)
{
  free(t->scores);
  free(t->is_match);
  t->scores = NULL;
  t->is_match = NULL;
}

//
// Striped SIMD fill (Farrar 2007, Bioinformatics 23:156-161)
//
// Each row of the matrices runs along seq_a. A row is split into ALN_VL
// segments of seg cells; lane l of vector k holds cell l*seg+k. match and
// gap_a only depend on the previous row, so they are computed a vector at a
// time. gap_b depends on the cell to its left: it is first computed within
// each segment, then corrected by Farrar's 'lazy F' loop, which carries
// values across segment boundaries while they improve a cell (what each
// lane passes on is found first, with a scan over the lanes).
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
size_t alignment_band_keep_max_cells = 0; // 0: a quarter of the memory
size_t alignment_band_wfa_budget = 16;

// 8- and 16-bit lanes for Smith-Waterman fills (see alignment_narrow_limit())
bool alignment_narrow_enabled = true, alignment_byte_enabled = true;
static atomic_size_t alignment_narrow_fills, alignment_narrow_retries;
static atomic_size_t alignment_byte_fills, alignment_byte_retries;

// Difference-encoded Needleman-Wunsch fills (see alignment_diff_bits())
bool alignment_diff_enabled = true;
bool alignment_linear_fallback = false;
static atomic_size_t alignment_diff_fills;
// Band fills (each width tried), and bands sized by alignment_wfa_score()
static atomic_size_t alignment_band_tries, alignment_band_wfa_sized;

void alignment_band_counts(size_t *tries, size_t *wfa_sized)
{
  *tries = atomic_load(&alignment_band_tries);
  *wfa_sized = atomic_load(&alignment_band_wfa_sized);
}

size_t alignment_diff_count(void)
{
  return atomic_load(&alignment_diff_fills);
}

void alignment_narrow_counts(size_t *fills, size_t *retries)
{
  *fills = atomic_load(&alignment_narrow_fills);
  *retries = atomic_load(&alignment_narrow_retries);
}

void alignment_byte_counts(size_t *fills, size_t *retries)
{
  *fills = atomic_load(&alignment_byte_fills);
  *retries = atomic_load(&alignment_byte_retries);
}

// Tests: the allocation to fail (see alignment_alloc.h)
atomic_size_t alignment_alloc_fail_in;

void alignment_alloc_fail(size_t n)
{
  atomic_store(&alignment_alloc_fail_in, n);
}

size_t alignment_alloc_fail_pending(void)
{
  return atomic_load(&alignment_alloc_fail_in);
}

void alignment_options_default(alignment_options_t *opts)
{
  opts->max_threads = alignment_max_threads;
  opts->striped = alignment_striped_enabled;
  opts->band = alignment_band_enabled;
  opts->narrow = alignment_narrow_enabled;
  opts->byte = alignment_byte_enabled;
  opts->diff = alignment_diff_enabled;
  opts->linear_fallback = alignment_linear_fallback;
  opts->block_min_cells = alignment_block_min_cells;
  opts->block_min_cols = alignment_block_min_cols;
  opts->band_strip_rows = alignment_band_strip_rows;
  opts->band_min_width = alignment_band_min_width;
  opts->band_keep_max_cells = alignment_band_keep_max_cells;
  opts->band_wfa_budget = alignment_band_wfa_budget;
  opts->sw_rect_min_cells = smith_waterman_rect_min_cells;
  opts->sw_strip_min_cells = smith_waterman_strip_min_cells;
  opts->hb_max_threads = hirschberg_max_threads;
  opts->hb_band = hirschberg_band_enabled;
  opts->hb_block_min_cells = hirschberg_block_min_cells;
  opts->hb_block_min_cols = hirschberg_block_min_cols;
  opts->hb_band_strip_rows = hirschberg_band_strip_rows;
  opts->hb_band_min_width = hirschberg_band_min_width;
  opts->hb_reuse_depth = hirschberg_reuse_depth;
  opts->hb_reuse_min_rows = hirschberg_reuse_min_rows;
}

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

// Wait until the block on the left has finished row y. Returns false if a
// block stopped meanwhile (see alignment_block_t.stop; stop may be NULL)
static inline bool alignment_wait_row(atomic_size_t *left_done, size_t y,
                                      atomic_bool *stop)
{
  unsigned spins = 0;
  while(atomic_load_explicit(left_done, memory_order_acquire) < y) {
    if(stop != NULL &&
       atomic_load_explicit(stop, memory_order_relaxed)) return false;
    if(++spins > 256) { spins = 0; sched_yield(); }
  }
  return true;
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
  // A scan for the best cell (track_best) gives up, setting *stop and
  // *gave_up, once the best so far reaches giveup_rate per row filled,
  // from a sixteenth of the rows on (0: never; see alignment_sw_scan_short())
  double giveup_rate;
  atomic_bool *gave_up;
  // A strip of a band, which takes thinner blocks than other fills (see
  // alignment_fill_threads())
  bool band;
  bool given_up; // the fill gave up (set by alignment_fill_blocks())
  // Scans: row k*ckpt_every of the three matrices goes to ckpt (see
  // alignment_sw_rows()); NULL for none
  score_t *ckpt;
  size_t ckpt_every;
  // Fill with 8- or 16-bit lanes (4 or 2 times as many; bits 8 or 16, 0
  // for 32-bit lanes) where alignment_narrow_limit() allows: blocks stop
  // once a row's best match score passes narrow_limit, setting *stop, and
  // the fill is done again with wider lanes
  int bits;
  score_t narrow_limit;
  // Needleman-Wunsch traceback bytes from difference-encoded lanes of diff
  // bits (8 or 16; 0: absolute scores), with floor for -infinity and
  // ceiling above any real difference (see alignment_diff.inc)
  int diff;
  score_t floor, ceiling;
  // Shared by the blocks of a fill: *stop is set when one stops early (the
  // others then stop too), and also *nomem if it ran out of memory
  atomic_bool *stop, *nomem;
} alignment_block_t;

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
static int alignment_diff_bits(const sub_table_t *table,
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
static bool alignment_narrow_limit(const sub_table_t *table,
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
static void alignment_lanes(alignment_block_t *tmpl,
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
  return NULL;
}

// Threads (column blocks) for a striped fill of width x rows cells
// A fill across the whole width (whole: every fill but a band's strips)
// takes blocks of at least twice block_min_cols columns: thinner blocks
// finish a row too soon after the block on their left, so threads mostly
// wait (20 kb: 16 threads no faster than 10, with 0.79 -> 0.44 s of CPU).
// A band's strips keep the minimum, or narrow bands would lose their
// threads and stop paying (50 kb at 10%: 0.15 -> 0.29 s)
static size_t alignment_fill_threads(const alignment_options_t *o,
                                    size_t width, size_t rows, bool whole)
{
  int threads = o->max_threads;
  if(threads <= 0) {
    long cpus = sysconf(_SC_NPROCESSORS_ONLN);
    threads = cpus > 1 ? (int)(cpus / 2) : 1;
  }
  // At least one vector of columns per block, so no block is empty
  const size_t min_cols = MAX2(whole ? 2 * o->block_min_cols
                                     : o->block_min_cols, (size_t)ALN_VL);
  if(threads > 1 && (double)width * rows >= o->block_min_cells &&
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
static bool alignment_fill_blocks(alignment_block_t *tmpl)
{
  aligner_t *aligner = tmpl->aligner;
  const size_t w = aligner->score_width, h = aligner->score_height;
  const size_t len_i = tmpl->width, len_j = tmpl->y1 - tmpl->y0 + 1;
  const size_t nblocks = alignment_fill_threads(&aligner->opts, len_i, len_j,
                                                !tmpl->band);
  size_t i, m;

  // Blocks of equal width (a multiple of ALN_VL); the last takes the rest
  size_t width = (len_i / nblocks) / ALN_VL * ALN_VL;
  alignment_block_t args[nblocks];
  atomic_size_t done[nblocks];
  pthread_t thread_ids[nblocks];
  bool started[nblocks];

  atomic_bool stop, nomem, gave_up;
  atomic_init(&stop, false);
  atomic_init(&nomem, false);
  atomic_init(&gave_up, false);
  if(tmpl->bits == 8) atomic_fetch_add(&alignment_byte_fills, 1);
  if(tmpl->bits == 16) atomic_fetch_add(&alignment_narrow_fills, 1);
  if(tmpl->diff != 0) atomic_fetch_add(&alignment_diff_fills, 1);

  // With trace, each block's last column goes here for the next block
  score_t *cols = NULL;
  if(tmpl->trace != NULL || tmpl->score_only) {
    if(h > SIZE_MAX / sizeof(score_t) / 3 / nblocks) cols = NULL;
    else cols = aln_malloc(nblocks * 3 * h * sizeof(score_t));
    if(cols == NULL) return false;
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
    b->stop = &stop;
    b->nomem = &nomem;
    b->gave_up = &gave_up;
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

  if(atomic_load(&nomem)) {
    free(cols);
    return false;
  }
  tmpl->given_up = atomic_load(&gave_up);
  if(tmpl->given_up) {
    free(cols);
    return true;
  }

  // Scores grew past 8- or 16-bit lanes: all again with wider ones (16 bits
  // after 8 where they are allowed, else 32)
  if(tmpl->bits != 0 && atomic_load(&stop)) {
    free(cols);
    if(tmpl->bits == 8) {
      atomic_fetch_add(&alignment_byte_retries, 1);
      const scoring_t *scoring = aligner->scoring;
      tmpl->bits = aligner->opts.narrow &&
                   alignment_narrow_limit(tmpl->table, scoring, tmpl->min,
                                          len_i, 16, &tmpl->narrow_limit) &&
                   tmpl->target <= tmpl->narrow_limit ? 16 : 0;
    }
    else {
      atomic_fetch_add(&alignment_narrow_retries, 1);
      tmpl->bits = 0;
    }
    return alignment_fill_blocks(tmpl);
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

  // Row 0 and column 0 are all 0 in Smith-Waterman
  score_t *zeros = aln_calloc(len_a + len_b + 2, sizeof(score_t));
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
  return alignment_sw_scan_any(a, len_a, b, len_b, scoring, target, opts,
                               false, NULL, 0, best, x, y);
}

bool alignment_sw_scan_short(const char *a, size_t len_a,
                             const char *b, size_t len_b,
                             const scoring_t *scoring,
                             const alignment_options_t *opts,
                             score_t *best, size_t *x, size_t *y)
{
  return alignment_sw_scan_any(a, len_a, b, len_b, scoring, 0, opts, true,
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
  return alignment_sw_scan_any(a, len_a, b, len_b, scoring, 0, opts, false,
                               ckpt, every, best, x, y);
}

bool alignment_sw_trace_rows(const char *a, size_t len_a,
                             const char *b, size_t len_b,
                             const scoring_t *scoring,
                             const alignment_options_t *opts,
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
  sub_table_t table;
  if(trace_rows == NULL ||
     !sub_table_build(&table, scoring, a, len_a, b, len_b)) {
    free(trace_rows);
    return false;
  }
  size_t r;
  for(r = y0; r <= y1; r++)
    trace_rows[r] = (ptrdiff_t)((r - y0) * xmax) - 1;
  // 32-bit lanes: row y0-1 may hold any score
  alignment_block_t tmpl = {
    .aligner = &tmp, .table = &table, .min = 0,
    .r0 = {row, row + w, row + 2*w}, .lc = {&zero, &zero, &zero},
    .lstride = 0, .width = xmax, .y0 = y0, .y1 = y1,
    .trace = trace, .trace_rows = trace_rows};
  bool filled = alignment_fill_blocks(&tmpl);
  sub_table_free(&table);
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
                             const alignment_options_t *opts,
                             const score_t *row, size_t y0, size_t y1,
                             size_t xmax, uint8_t *trace)
{
  (void)a; (void)len_a; (void)b; (void)len_b; (void)scoring; (void)opts;
  (void)row; (void)y0; (void)y1; (void)xmax; (void)trace;
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

// Fill in traceback matrix. Returns false if out of memory
static bool alignment_fill_matrices(aligner_t *aligner, char is_sw)
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
  if(!sub_table_build(&table, scoring, aligner->seq_a, len_i,
                      aligner->seq_b, len_j)) return false;

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
    // The first row's gaps start the alignment, and end it too where seq_b
    // is empty (the first row is then the last); the first column's the
    // same with seq_a
    const bool free_row = scoring->no_start_gap_penalty ||
                          (scoring->no_end_gap_penalty && len_j == 0);
    const bool free_col = scoring->no_start_gap_penalty ||
                          (scoring->no_end_gap_penalty && len_i == 0);

    // work along first row -> [i][0]
    for(i = 1; i < score_width; i++)
    {
      match_scores[i] = min;

      // Think carefully about which way round these are
      gap_a_scores[i] = min;
      gap_b_scores[i] = free_row ? 0
                        : scoring->gap_open + (int)i * scoring->gap_extend;
    }

    // work down first column -> [0][j]
    for(j = 1, index = score_width; j < score_height; j++, index += score_width)
    {
      match_scores[index] = min;

      // Think carefully about which way round these are
      gap_a_scores[index] = free_col ? 0
                            : scoring->gap_open + (int)j * scoring->gap_extend;
      gap_b_scores[index] = min;
    }
  }

  #ifdef ALN_STRIPED
  // The striped fill handles the uniform recurrence only; options that
  // change it for particular cells use the scalar loop below
  if(aligner->opts.striped && len_i > 0 &&
     !scoring->no_end_gap_penalty && !scoring->no_mismatches &&
     !no_gaps_in_a && !no_gaps_in_b)
  {
    alignment_block_t tmpl = {
      .aligner = aligner, .table = &table, .min = min,
      .r0 = {match_scores, gap_a_scores, gap_b_scores},
      .lc = {match_scores, gap_a_scores, gap_b_scores},
      .lstride = score_width, .width = len_i, .y0 = 1, .y1 = len_j};
    bool filled = alignment_fill_blocks(&tmpl);
    sub_table_free(&table);
    return filled;
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
  return true;
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
    void *ptr = aln_aligned_alloc(huge, mem);
    if(ptr != NULL) madvise(ptr, mem, MADV_HUGEPAGE); // only a hint
    return ptr;
  }
#endif
  return aln_malloc(mem);
}

#ifdef ALN_STRIPED
// Room for the traceback bytes of `cells` cells. Returns false if out of
// memory
static bool alignment_trace_capacity(aligner_t *aligner, size_t cells)
{
  if(aligner->trace_capacity < cells) {
    free(aligner->trace);
    aligner->trace_capacity = ROUNDUP2POW(cells);
    aligner->trace = alignment_alloc_matrix(aligner->trace_capacity);
    if(aligner->trace == NULL) {
      aligner->trace_capacity = 0;
      return false;
    }
  }
  return true;
}

// What a compact fill did
typedef enum { FILL_SKIPPED, FILL_DONE, FILL_NOMEM } fill_result_t;

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
  size_t s, i, j, found = 0;
  for(s = 0; s < SAMPLES; s++) {
    i = (len_a - K) * s / (SAMPLES - 1);
    const size_t at = (size_t)((double)i * (double)(len_b - K) /
                               (double)(len_a - K));
    const size_t lo = at > WINDOW ? at - WINDOW : 0;
    const size_t hi = MIN2(len_b - K, at + WINDOW);
    for(j = lo; j <= hi; j++)
      if(memcmp(a + i, b + j, K) == 0) { found++; break; }
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
  size_t strip, nstrips;
  size_t *c0;
  score_t min;
  size_t *off;
  score_t *rows;
  score_t *col0; // column 0, 3 x score_height, for strips that reach it
  score_t *top; // 3 x score_width: a strip's top row by column
  score_t *lcol; // 3 x (strip + 1): its left column
  size_t loaded, loaded_y, loaded_x; // the strip in trace, to row, column
} trace_strips_t;

static void trace_strips_free(trace_strips_t *ts)
{
  if(ts == NULL) return;
  free(ts->off);
  free(ts->c0);
  free(ts->rows);
  free(ts->col0);
  free(ts->top);
  free(ts->lcol);
  free(ts);
}

// Fill strip s of a band again, rows y0 .. y (y in the strip) and its
// columns up to x, into trace from its kept top row
static bool trace_strips_fill(aligner_t *aligner, trace_strips_t *ts,
                              size_t s, size_t y, size_t x)
{
  const size_t w = aligner->score_width, h = aligner->score_height;
  const size_t len_i = w - 1, len_j = h - 1;
  const size_t y0 = 1 + s * ts->strip, span = ts->off[s+1] - ts->off[s];
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
  sub_table_t table;
  if(!sub_table_build(&table, aligner->scoring, aligner->seq_a, len_i,
                      aligner->seq_b, len_j)) return false;
  alignment_block_t tmpl = {
    .aligner = aligner, .table = &table, .min = ts->min,
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
  sub_table_free(&table);
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
  const size_t s = (y - 1) / ts->strip;
  if(s == ts->loaded && y <= ts->loaded_y && x <= ts->loaded_x) return true;
  return trace_strips_fill(aligner, ts, s, y, x);
}

// The diagonals x - y the strips below row y need, from row y's cells
// (columns c0+1 .. c1; the rest are outside the band) and a lower bound lb
// on the best score. A cell scoring v can end no better than its bound,
// v + smax min(i, j) + e |i - j| (i, j the letters left of each sequence,
// e the gap extension). A step towards diagonal len_i - len_j does not
// raise the bound and one away lowers it by at least -2e, so from a cell
// whose bound is lb + k, cells with bounds of lb or more lie only on the
// diagonals from it to len_i - len_j and k / -2e beyond. A cell with a
// bound under lb is on no best path, and could not tie with a cell on one,
// so leaving it out changes no cell on a best path (the computational
// volume A*PA2 narrows its band to: Groot Koerkamp 2024, WABI). Narrows
// *lo .. *hi to those diagonals; no change if the gap extension is 0, or
// no cell's bound reaches lb
static void alignment_band_narrow(score_t *const row[3], size_t y,
                                  size_t c0, size_t c1,
                                  size_t len_i, size_t len_j, score_t smax,
                                  int gap_extend, score_t lb,
                                  ptrdiff_t *lo, ptrdiff_t *hi)
{
  if(gap_extend >= 0) return;
  const long long e2 = -2LL * gap_extend;
  const ptrdiff_t D = (ptrdiff_t)len_i - (ptrdiff_t)len_j;
  ptrdiff_t nlo = PTRDIFF_MAX, nhi = PTRDIFF_MIN;
  size_t x;
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
static fill_result_t alignment_fill_banded(aligner_t *aligner,
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
    const double several = MAX2(2.0 * (double)aligner->opts.block_min_cols,
                                one / MAX2(1.0, full_threads / 2));
    if(MIN2(one, several) >=
       pays * (double)len_i / MAX2(1.0, full_threads / 2)) return FILL_SKIPPED;
  }

  size_t bw = aligner->opts.band_min_width;
  if(have_exact) {
    bw = MAX2(bw, alignment_band_width(smax, scoring->gap_extend,
                                       len_i, len_j, exact));
    atomic_fetch_add(&alignment_band_wfa_sized, 1);
  }
  // A lower bound on the best score, once known, narrows strips
  bool have_lb = have_exact;
  score_t lb = have_exact ? exact : 0;
  score_t *top = NULL, *lcol = NULL;
  fill_result_t done = FILL_SKIPPED;
  int widened = 0; // times this band was widened to get more threads

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
      // A quarter of the physical memory (or 1 GB if unknown): the strips
      // cost a second fill of the band as the traceback reaches them
      long pages = sysconf(_SC_PHYS_PAGES), page = sysconf(_SC_PAGESIZE);
      keep_max = pages > 0 && page > 0 ? (double)pages * (double)page / 4
                                       : (double)((size_t)1 << 30);
    }
    const bool keep = (double)len_j * (cols + (double)strip) > keep_max;
    if(keep && !aligner->opts.band_strip_rows) {
      // At most about sqrt(12 x rows) rows, so a wide band's one strip of
      // bytes stays small; a narrow band keeps its usual strips (taller
      // ones would reach further across: 300 kb at 0.1%, 0.68 -> 2.0 s),
      // and their top rows take only about 60 bytes per row
      size_t most;
      for(most = 64; most * most < 12 * len_j; most++) {}
      strip = MIN2(strip, most);
    }
    // A fill on t threads runs about max(1, t/2) times faster (16 threads:
    // 6-8x); a narrow band runs on one
    const double band_threads =
      (double)alignment_fill_threads(&aligner->opts,
                                     (size_t)cols + strip, strip, false);
    if((cols + strip) / MAX2(1.0, band_threads / 2) >=
       pays * (double)len_i / MAX2(1.0, full_threads / 2)) {
      // A band too small for its strips to be split across all the threads
      // can cost more than a wider one that is (50 kb at 3%: 1,677
      // diagonals each side ran on one thread, and the whole matrix was
      // filled instead). So, while it has fewer threads than the whole
      // fill, try it a quarter wider, a few times: a wider band holds every
      // path the narrower one would, so its proof, and the alignment, are
      // unchanged
      if(widened < 8 && band_threads < full_threads && bw < len_i) {
        bw += MAX2(bw / 4, (size_t)1);
        widened++;
        continue;
      }
      break;
    }
    widened = 0;

    atomic_fetch_add(&alignment_band_tries, 1);
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
    size_t cells = 0, kept = 0, nstrips = 0;
    for(y0 = 1; y0 <= len_j; y0 = y1 + 1, nstrips++) {
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
        ts->nstrips = ts->loaded = nstrips;
        ts->c0 = aln_malloc(nstrips * sizeof(size_t));
        ts->min = min;
        ts->off = aln_malloc((nstrips + 1) * sizeof(size_t));
        ts->rows = aln_malloc(kept * sizeof(score_t));
        ts->col0 = aln_malloc(3 * h * sizeof(score_t));
        ts->top = aln_malloc(3 * w * sizeof(score_t));
        ts->lcol = aln_malloc(3 * (strip + 1) * sizeof(score_t));
      }
      if(ts == NULL || ts->off == NULL || ts->c0 == NULL ||
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
    size_t s = 0, used = 0, c0_prev = 0;
    ptrdiff_t sdlo = dlo, sdhi = dhi; // this strip's diagonals
    for(y0 = 1; y0 <= len_j; y0 = y1 + 1, s++)
    {
      // Not left of the strip above's left column: a path only moves right,
      // so one reaching a cell there passes left of the strip above too
      sdlo = MAX2(sdlo, (ptrdiff_t)c0_prev + 1 - (ptrdiff_t)y0);
      alignment_strip(y0, strip, len_i, len_j, sdlo, sdhi, &y1, &c0, &c1);
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
      if(have_lb)
        alignment_band_narrow(topm, y1, c0, c1, len_i, len_j, smax,
                              scoring->gap_extend, lb, &sdlo, &sdhi);
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
        ts->loaded = nstrips - 1;
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

  free(top);
  free(lcol);
  aligner->trace_banded = done == FILL_DONE;
  return done;
}
#endif

// Compact mode: run the striped fill keeping one traceback byte per cell
// instead of the score matrices. Returns false (nothing done) where the
// striped fill does not apply, with free start gaps (the traceback changes
// penalties along the first row and column), or where the first row and
// column fall below min: working out the traceback bytes adds a penalty to
// each predecessor, which could then overflow (the fill adds it only to the
// best one). FILL_NOMEM: out of memory, with no trace
static fill_result_t alignment_fill_trace(aligner_t *aligner, char is_sw)
{
#ifdef ALN_STRIPED
  const scoring_t *scoring = aligner->scoring;
  const size_t w = aligner->score_width, h = aligner->score_height;
  const size_t len_i = w-1, len_j = h-1;
  size_t i, m;

  bool no_gaps_in_a, no_gaps_in_b;
  alignment_no_gaps(scoring, len_i, len_j, &no_gaps_in_a, &no_gaps_in_b);
  if(!aligner->opts.striped || len_i == 0 || len_j == 0 ||
     len_i >= INT_MAX || len_j >= INT_MAX ||
     w + h > SIZE_MAX / (3 * sizeof(score_t)) ||
     scoring->no_start_gap_penalty ||
     scoring->no_end_gap_penalty || scoring->no_mismatches ||
     no_gaps_in_a || no_gaps_in_b) return FILL_SKIPPED;

  const score_t min = alignment_min_score(scoring, is_sw);
  if(!is_sw) {
    // Gap scores along the first row and column run from these two
    long long edge1 = (long long)scoring->gap_open + scoring->gap_extend;
    long long edgen = (long long)scoring->gap_open +
                      (long long)MAX2(len_i, len_j) * scoring->gap_extend;
    if(MIN2(edge1, edgen) < min || MAX2(edge1, edgen) > INT_MAX)
      return FILL_SKIPPED;
  }

  // Row 0 and column 0, as alignment_fill_matrices() sets them
  score_t *edges = aln_malloc(3 * (w + h) * sizeof(score_t));
  if(edges == NULL) return FILL_NOMEM;
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
  if(!sub_table_build(&table, scoring, aligner->seq_a, len_i,
                      aligner->seq_b, len_j)) {
    free(edges);
    return FILL_NOMEM;
  }
  aligner->trace_banded = false;
  trace_strips_free(aligner->trace_strips);
  aligner->trace_strips = NULL;
  fill_result_t done = FILL_SKIPPED;
  if(!is_sw) done = alignment_fill_banded(aligner, &table, min, row0, col0);
  if(done == FILL_SKIPPED && !alignment_trace_capacity(aligner, w * h))
    done = FILL_NOMEM;
  else if(done == FILL_SKIPPED)
  {
    alignment_block_t tmpl = {
      .aligner = aligner, .table = &table, .min = min,
      .r0 = {row0[MATCH], row0[GAP_A], row0[GAP_B]},
      .lc = {col0[MATCH], col0[GAP_A], col0[GAP_B]}, .lstride = 1,
      .trace = aligner->trace, .track_best = is_sw,
      .width = len_i, .y0 = 1, .y1 = len_j};
    if(is_sw) alignment_lanes(&tmpl, &aligner->opts, scoring, len_i);
    else if(aligner->opts.diff)
      tmpl.diff = alignment_diff_bits(&table, scoring, &tmpl.floor,
                                      &tmpl.ceiling);
    done = alignment_fill_blocks(&tmpl) ? FILL_DONE : FILL_NOMEM;
  }
  sub_table_free(&table);
  free(edges);
  return done;
#else
  (void)aligner; (void)is_sw;
  return FILL_SKIPPED;
#endif
}

// The score matrices of len_a x len_b would not fit in memory addresses
static bool alignment_too_long(size_t len_a, size_t len_b)
{
  return len_a == SIZE_MAX || len_b == SIZE_MAX ||
         len_b + 1 > SIZE_MAX / (len_a + 1) / sizeof(score_t);
}

bool alignment_scores_may_overflow(const scoring_t *scoring,
                                   const char *a, size_t len_a,
                                   const char *b, size_t len_b, bool is_sw)
{
  // Upwards: a path takes at most len_a + len_b steps, each adding at most
  // top, the highest score of a gap or of a pair in the sequences
  bool in_a[256] = {false}, in_b[256] = {false};
  size_t i, j;
  for(i = 0; i < len_a; i++) in_a[(unsigned char)a[i]] = true;
  for(i = 0; i < len_b; i++) in_b[(unsigned char)b[i]] = true;
  double top = MAX2((double)scoring->gap_open + scoring->gap_extend,
                    (double)scoring->gap_extend);
  double low = 0; // the lowest score of a pair in the sequences
  for(j = 0; j < 256; j++) {
    if(!in_b[j]) continue;
    for(i = 0; i < 256; i++) {
      if(!in_a[i] || !scoring_has_pair(scoring, (char)i, (char)j)) continue;
      int score;
      bool is_match;
      scoring_lookup(scoring, (char)i, (char)j, &score, &is_match);
      top = MAX2(top, (double)score);
      low = MIN2(low, (double)score);
    }
  }
  const double half = (double)INT_MAX / 2;
  if(((double)len_a + (double)len_b + 1) * top > half) return true;

  // Downwards: cells are clamped at the min score (see
  // alignment_min_score()). That leaves Needleman-Wunsch's first row and
  // column, gaps of every length up to the longer sequence, which must be
  // ints, and its score, at least that of a pair at each position of the
  // shorter sequence and one gap
  if(!is_sw) {
    const double open = scoring->gap_open, ext = scoring->gap_extend;
    const double longest = (double)MAX2(len_a, len_b);
    const double shortest = (double)MIN2(len_a, len_b);
    const double rest = longest - shortest;
    const double edges[] = {longest * ext, open + ext, open + longest * ext,
                            shortest * low + (rest > 0 ? open + rest * ext : 0)};
    for(i = 0; i < sizeof(edges) / sizeof(edges[0]); i++)
      if(edges[i] < (double)INT_MIN || edges[i] > (double)INT_MAX) return true;
  }
  return false;
}

bool alignment_scoring_complete(const scoring_t *scoring,
                                const char *a, size_t len_a,
                                const char *b, size_t len_b,
                                aln_error_t *err)
{
  bool in_a[256] = {false}, in_b[256] = {false};
  size_t i, j;
  for(i = 0; i < len_a; i++) in_a[(unsigned char)a[i]] = true;
  for(i = 0; i < len_b; i++) in_b[(unsigned char)b[i]] = true;

  // In the order sub_table_build() looks pairs up, so the pair reported is
  // the one scoring_lookup() would stop at
  for(j = 0; j < 256; j++) {
    if(!in_b[j]) continue;
    for(i = 0; i < 256; i++) {
      if(in_a[i] && !scoring_has_pair(scoring, (char)i, (char)j)) {
        if(err != NULL) {
          err->status = ALN_ERR_SCORING;
          err->a = scoring->case_sensitive ? (char)i : (char)tolower((int)i);
          err->b = scoring->case_sensitive ? (char)j : (char)tolower((int)j);
        }
        return false;
      }
    }
  }
  return true;
}

aln_status_t alignment_check(const scoring_t *scoring,
                             const char *a, size_t len_a,
                             const char *b, size_t len_b,
                             bool scoring_first, aln_error_t *err)
{
  aln_error_t e = {ALN_OK, 0, 0};
  bool ok = !scoring_first ||
            alignment_scoring_complete(scoring, a, len_a, b, len_b, &e);
  if(ok && alignment_too_long(len_a, len_b)) {
    e.status = ALN_ERR_TOO_LONG;
    ok = false;
  }
  if(ok && alignment_scores_may_overflow(scoring, a, len_a, b, len_b,
                                         scoring_first)) {
    e.status = ALN_ERR_OVERFLOW;
    ok = false;
  }
  if(ok && !scoring_first)
    alignment_scoring_complete(scoring, a, len_a, b, len_b, &e);
  if(err != NULL) *err = e;
  return e.status;
}

void alignment_fprint_error(FILE *out, const aln_error_t *err)
{
  switch(err->status) {
    case ALN_OK: break;
    case ALN_ERR_SCORING:
      fprintf(out, "Error: Unknown character pair (%c,%c) and "
                   "match/mismatch have not been set\n", err->a, err->b);
      break;
    case ALN_ERR_TOO_LONG:
      fprintf(out, "Error: Sequences too long to align\n");
      break;
    case ALN_ERR_OVERFLOW:
      fprintf(out, "Error: Scores could overflow: use smaller scores or "
                   "shorter sequences\n");
      break;
    case ALN_ERR_UNSUPPORTED:
      fprintf(out, "Error: linear space does not support these scoring "
                   "options (or penalties this large)\n");
      break;
    case ALN_ERR_NOMEM:
      fprintf(out, "Error: Out of memory\n");
      break;
  }
}

void aligner_align(aligner_t *aligner,
                   const char *seq_a, const char *seq_b,
                   size_t len_a, size_t len_b,
                   const scoring_t *scoring, char is_sw)
{
  aln_error_t e = {ALN_OK, 0, 0};
  e.status = aligner_align_checked(aligner, seq_a, seq_b, len_a, len_b,
                                   scoring, is_sw);
  if(e.status != ALN_OK) {
    alignment_fprint_error(stderr, &e);
    exit(EXIT_FAILURE);
  }
}

aln_status_t aligner_align_checked(aligner_t *aligner,
                                   const char *seq_a, const char *seq_b,
                                   size_t len_a, size_t len_b,
                                   const scoring_t *scoring, char is_sw)
{
  aligner->scoring = scoring;
  if(aligner->options != NULL) aligner->opts = *aligner->options;
  else alignment_options_default(&aligner->opts);
  aligner->seq_a = seq_a;
  aligner->seq_b = seq_b;
  aligner->score_width = len_a+1;
  aligner->score_height = len_b+1;

  // The same checks as the *_checked() aligners
  aligner->traced = false;
  if(alignment_too_long(len_a, len_b)) return ALN_ERR_TOO_LONG;
  if(alignment_scores_may_overflow(scoring, seq_a, len_a, seq_b, len_b,
                                   is_sw)) return ALN_ERR_OVERFLOW;

  if(aligner->compact) {
    fill_result_t done = alignment_fill_trace(aligner, is_sw);
    if(done == FILL_NOMEM) return ALN_ERR_NOMEM;
    aligner->traced = done == FILL_DONE;
    if(aligner->traced) return ALN_OK;
  }

  size_t new_capacity = aligner->score_width * aligner->score_height;

  if(aligner->capacity < new_capacity)
  {
    aligner->capacity = ROUNDUP2POW(new_capacity);
    // Rounding up can take the bytes past SIZE_MAX: no rounding then
    if(aligner->capacity > SIZE_MAX / sizeof(score_t))
      aligner->capacity = new_capacity;
    size_t mem = sizeof(score_t) * aligner->capacity;
    // The fill overwrites every cell it uses, so old contents are not kept
    free(aligner->match_scores);
    free(aligner->gap_a_scores);
    free(aligner->gap_b_scores);
    aligner->match_scores = alignment_alloc_matrix(mem);
    aligner->gap_a_scores = alignment_alloc_matrix(mem);
    aligner->gap_b_scores = alignment_alloc_matrix(mem);
    if(!aligner->match_scores || !aligner->gap_a_scores || !aligner->gap_b_scores) {
      free(aligner->match_scores);
      free(aligner->gap_a_scores);
      free(aligner->gap_b_scores);
      aligner->match_scores = NULL;
      aligner->gap_a_scores = aligner->gap_b_scores = NULL;
      aligner->capacity = 0;
      return ALN_ERR_NOMEM;
    }
  }

  return alignment_fill_matrices(aligner, is_sw) ? ALN_OK : ALN_ERR_NOMEM;
}

void aligner_destroy(aligner_t *aligner)
{
  trace_strips_free(aligner->trace_strips);
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
  if(!alignment_reserve(result, strlength)) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

bool alignment_reserve(alignment_t* result, size_t strlength)
{
  size_t capacity = strlength+1;
  if(result->capacity < capacity)
  {
    capacity = ROUNDUP2POW(capacity);
    // Each string keeps its contents if the other cannot grow
    char *a = aln_realloc(result->result_a, sizeof(char)*capacity);
    if(a == NULL) return false;
    result->result_a = a;
    char *b = aln_realloc(result->result_b, sizeof(char)*capacity);
    if(b == NULL) return false;
    result->result_b = b;
    result->capacity = capacity;
  }
  return true;
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
