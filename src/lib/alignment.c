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
#include <stdatomic.h>
#include <unistd.h> // sysconf

#include "alignment.h"
#include "alignment_macros.h"
#include "alignment_simd.h"
#include "hirschberg.h" // defaults for alignment_options_t
#include "smith_waterman.h"
#include "alignment_testing.h"
#include "alignment_alloc.h"

#include "alignment_fill_internal.h"

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
// A band's strips: narrowed strips have few columns, and their blocks run
// only for the strip's rows, so they split finer than a whole fill
size_t alignment_band_block_min_cells = 262144;
size_t alignment_band_block_min_cols = 256;

// Banded compact fills (see alignment_fill_banded())
bool alignment_band_enabled = true;
size_t alignment_band_strip_rows = 0; // 0: automatic
size_t alignment_band_min_width = 64;
bool alignment_seeds_enabled = true;
size_t alignment_seed_length = 0; // 0: automatic
size_t alignment_band_keep_max_cells = 0; // 0: a quarter of the memory
size_t alignment_band_wfa_budget = 16;

// 8- and 16-bit lanes for Smith-Waterman fills (see alignment_narrow_limit())
bool alignment_narrow_enabled = true, alignment_byte_enabled = true;
atomic_size_t aln_stat_narrow_fills, aln_stat_narrow_retries;
atomic_size_t aln_stat_byte_fills, aln_stat_byte_retries;

// Difference-encoded Needleman-Wunsch fills (see alignment_diff_bits())
bool alignment_diff_enabled = true;
bool alignment_linear_fallback = false;
size_t alignment_linear_above = 0;
bool alignment_wfa_alignment = false;
size_t alignment_wfa_budget = 0;

size_t alignment_wfa_default_budget(size_t len_a, size_t len_b)
{
  const size_t per_letter = 256, max_cells = (size_t)1 << 28;
  const size_t total = len_a <= SIZE_MAX - len_b ? len_a + len_b : SIZE_MAX;
  return total <= max_cells / per_letter ? per_letter * total : max_cells;
}
atomic_size_t aln_stat_diff_fills;
// Band fills (each width tried), and bands sized by alignment_wfa_score()
atomic_size_t aln_stat_band_tries, aln_stat_band_wfa_sized;

void alignment_band_counts(size_t *tries, size_t *wfa_sized)
{
  *tries = atomic_load(&aln_stat_band_tries);
  *wfa_sized = atomic_load(&aln_stat_band_wfa_sized);
}

size_t alignment_diff_count(void)
{
  return atomic_load(&aln_stat_diff_fills);
}

void alignment_narrow_counts(size_t *fills, size_t *retries)
{
  *fills = atomic_load(&aln_stat_narrow_fills);
  *retries = atomic_load(&aln_stat_narrow_retries);
}

void alignment_byte_counts(size_t *fills, size_t *retries)
{
  *fills = atomic_load(&aln_stat_byte_fills);
  *retries = atomic_load(&aln_stat_byte_retries);
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
  opts->linear_above = alignment_linear_above;
  opts->wfa_align = alignment_wfa_alignment;
  opts->wfa_budget = alignment_wfa_budget;
  opts->sw_fast_rect = smith_waterman_fast_rect;
  opts->block_min_cells = alignment_block_min_cells;
  opts->block_min_cols = alignment_block_min_cols;
  opts->band_block_min_cells = alignment_band_block_min_cells;
  opts->band_block_min_cols = alignment_band_block_min_cols;
  opts->band_strip_rows = alignment_band_strip_rows;
  opts->band_min_width = alignment_band_min_width;
  opts->band_seeds = alignment_seeds_enabled;
  opts->seed_length = alignment_seed_length;
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
  // column (unless free), gaps of every length up to the longer sequence,
  // and its score, all of which must stay within [min, INT_MAX]
  if(!is_sw) {
    const double floor = (double)alignment_min_score(scoring, false);
    const double open = scoring->gap_open, ext = scoring->gap_extend;
    const bool free_row = scoring->no_start_gap_penalty ||
                          (scoring->no_end_gap_penalty && len_b == 0);
    const bool free_col = scoring->no_start_gap_penalty ||
                          (scoring->no_end_gap_penalty && len_a == 0);
    const double edge_len = MAX2(free_row ? 0.0 : (double)len_a,
                                 free_col ? 0.0 : (double)len_b);
    const double longest = (double)MAX2(len_a, len_b);
    const double shortest = (double)MIN2(len_a, len_b);
    const double rest = longest - shortest;
    const bool free_tail = scoring->no_start_gap_penalty ||
                           scoring->no_end_gap_penalty;
    double score_lb;
    if(scoring->no_mismatches) {
      if(scoring->no_start_gap_penalty && scoring->no_end_gap_penalty)
        score_lb = 0;
      else if(free_tail)
        score_lb = shortest > 0 ? open + shortest * ext : 0;
      else
        score_lb = (len_a > 0 ? open + (double)len_a * ext : 0) +
                   (len_b > 0 ? open + (double)len_b * ext : 0);
    }
    else {
      score_lb = shortest * low +
                 (rest > 0 && !free_tail ? open + rest * ext : 0);
    }
    const double edges[] = {
      edge_len > 0 ? edge_len * ext : 0,
      edge_len > 0 ? open + ext : 0,
      edge_len > 0 ? open + edge_len * ext : 0,
      score_lb
    };
    // Strictly above the floor: a cell at the floor could not be told from
    // an unreachable one (clamped there), and the traceback would fail
    for(i = 0; i < sizeof(edges) / sizeof(edges[0]); i++)
      if(edges[i] <= floor || edges[i] > (double)INT_MAX) return true;
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

bool alignment_sw_gaps_unsupported(const scoring_t *scoring)
{
  // A gap scoring above 0 would grow out of the empty cells (clamped at
  // 0), reporting a gap where there is no alignment
  return scoring->gap_extend > 0 || scoring->gap_open + scoring->gap_extend > 0;
}

aln_status_t alignment_check(const scoring_t *scoring,
                             const char *a, size_t len_a,
                             const char *b, size_t len_b,
                             bool is_sw, aln_error_t *err)
{
  // Smith-Waterman checks the scoring pairs first, Needleman-Wunsch last
  aln_error_t e = {ALN_OK, 0, 0};
  bool ok = !is_sw ||
            alignment_scoring_complete(scoring, a, len_a, b, len_b, &e);
  if(ok && is_sw && alignment_sw_gaps_unsupported(scoring)) {
    e.status = ALN_ERR_UNSUPPORTED;
    ok = false;
  }
  if(ok && alignment_too_long(len_a, len_b)) {
    e.status = ALN_ERR_TOO_LONG;
    ok = false;
  }
  if(ok && alignment_scores_may_overflow(scoring, a, len_a, b, len_b, is_sw)) {
    e.status = ALN_ERR_OVERFLOW;
    ok = false;
  }
  if(ok && !is_sw)
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
      fprintf(out, "Error: the aligner does not support these scoring "
                   "options or penalties (linear space refuses some; "
                   "Smith-Waterman gaps must not score above 0)\n");
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
    aligner->capacity = alignment_grow_capacity(new_capacity);
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


alignment_t* alignment_try_create(size_t capacity)
{
  capacity = capacity < 16 ? 16 : ROUNDUP2POW(capacity);
  alignment_t *result = aln_calloc(1, sizeof(alignment_t));
  if(result == NULL) return NULL;
  result->result_a = aln_malloc(sizeof(char) * capacity);
  result->result_b = aln_malloc(sizeof(char) * capacity);
  if(result->result_a == NULL || result->result_b == NULL) {
    free(result->result_a);
    free(result->result_b);
    free(result);
    return NULL;
  }
  result->capacity = capacity;
  result->result_a[0] = result->result_b[0] = '\0';
  return result;
}

alignment_t* alignment_create(size_t capacity)
{
  alignment_t *result = alignment_try_create(capacity);
  if(result == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
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
  if(strlength == SIZE_MAX) return false; // its '\0' would not fit
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
  if(result == NULL) return;
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
     (long long)aligner->gap_a_scores[*arr_index] + prev_gap_a_penalty == *curr_score)
  {
    *curr_matrix = GAP_A;
    *curr_score = aligner->gap_a_scores[*arr_index];
  }
  else if((!no_gaps_in_b || *score_y == 0 || *score_y == len_j) &&
          (long long)aligner->gap_b_scores[*arr_index] + prev_gap_b_penalty == *curr_score)
  {
    *curr_matrix = GAP_B;
    *curr_score = aligner->gap_b_scores[*arr_index];
  }
  else if((long long)aligner->match_scores[*arr_index] + prev_match_penalty == *curr_score)
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

