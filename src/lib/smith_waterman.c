/*
 smith_waterman.c
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */

#include "smith_waterman_internal.h"

size_t smith_waterman_rect_min_cells = (size_t)1 << 24;
size_t smith_waterman_strip_min_cells = (size_t)1 << 28;
bool smith_waterman_fast_rect = false;
size_t smith_waterman_strip_rows = 0;
size_t smith_waterman_strip_alignments = 0;

// One step back through the matrices. Gives the same result as
// alignment_reverse_move() (same checks, same order), but uses the
// substitution table and fixed gap penalties, so it only applies when no
// option changes the penalties or forbids gaps for particular cells
static inline void _sw_reverse_move(const sw_aligner_t *sw,
                                    enum Matrix *curr_matrix,
                                    score_t *curr_score,
                                    size_t *score_x, size_t *score_y,
                                    size_t *arr_index)
{
  const aligner_t *aligner = &sw->aligner;

  if(aligner->traced)
  {
    // The traceback byte says which matrix the step goes to; the score is
    // the current one less the penalty of that step
    const scoring_t *scoring = aligner->scoring;
    const score_t open = scoring->gap_open + scoring->gap_extend;
    const score_t ext = scoring->gap_extend;
    int next = alignment_trace_next(aligner, *arr_index, *curr_matrix);
    score_t step;
    if(next == ALN_TRACE_NONE) {
      fprintf(stderr, "Program error: traceback fail (_sw_reverse_move)\n");
      exit(EXIT_FAILURE);
    }
    if(*curr_matrix == MATCH) {
      const sub_table_t *t = &sw->table;
      unsigned char ca = (unsigned char)aligner->seq_a[*score_x-1];
      unsigned char cb = (unsigned char)aligner->seq_b[*score_y-1];
      step = t->scores[t->code[cb]*t->n + t->code[ca]];
      *arr_index -= aligner->score_width + 1;
    }
    else if(*curr_matrix == GAP_A) {
      step = next == GAP_A ? ext : open;
      *arr_index -= aligner->score_width;
    }
    else {
      step = next == GAP_B ? ext : open;
      *arr_index -= 1;
    }
    if(*curr_matrix != GAP_A) (*score_x)--;
    if(*curr_matrix != GAP_B) (*score_y)--;
    *curr_matrix = (enum Matrix)next;
    *curr_score -= step;
    return;
  }

  if(!sw->fast_walk || (*curr_matrix != MATCH && *curr_matrix != GAP_A &&
                        *curr_matrix != GAP_B))
  {
    alignment_reverse_move(curr_matrix, curr_score, score_x, score_y,
                           arr_index, aligner);
    return;
  }

  const scoring_t *scoring = aligner->scoring;
  const score_t open = scoring->gap_open + scoring->gap_extend;
  const score_t ext = scoring->gap_extend;
  score_t from_a, from_b, from_m; // penalty for arriving from each matrix
  size_t prev;

  if(*curr_matrix == MATCH)
  {
    const sub_table_t *t = &sw->table;
    unsigned char ca = (unsigned char)aligner->seq_a[*score_x-1];
    unsigned char cb = (unsigned char)aligner->seq_b[*score_y-1];
    from_a = from_b = from_m = t->scores[t->code[cb]*t->n + t->code[ca]];
    prev = *arr_index - aligner->score_width - 1;
  }
  else if(*curr_matrix == GAP_A)
  {
    from_a = ext; from_b = open; from_m = open;
    prev = *arr_index - aligner->score_width;
  }
  else
  {
    from_a = open; from_b = ext; from_m = open;
    prev = *arr_index - 1;
  }

  enum Matrix next;
  score_t next_score;

  if(aligner->gap_a_scores[prev] + from_a == *curr_score) {
    next = GAP_A; next_score = aligner->gap_a_scores[prev];
  }
  else if(aligner->gap_b_scores[prev] + from_b == *curr_score) {
    next = GAP_B; next_score = aligner->gap_b_scores[prev];
  }
  else if(aligner->match_scores[prev] + from_m == *curr_score) {
    next = MATCH; next_score = aligner->match_scores[prev];
  }
  else {
    // No valid predecessor: the general version reports it and exits
    alignment_reverse_move(curr_matrix, curr_score, score_x, score_y,
                           arr_index, aligner);
    return;
  }

  if(*curr_matrix != GAP_A) (*score_x)--;
  if(*curr_matrix != GAP_B) (*score_y)--;
  *arr_index = prev;
  *curr_matrix = next;
  *curr_score = next_score;
}

// Sort indices by their matrix values
// Struct used to pass data to sort_match_indices
typedef struct
{
  score_t *match_scores;
  size_t score_width;
} MatrixSort;

// Compare two candidate cell indices: descending score, then ascending x, then y
static inline int cmp_match_index(size_t ia, size_t ib, const MatrixSort *ctx)
{
  const score_t *match_scores = ctx->match_scores;
  size_t score_width = ctx->score_width;

  long diff = (long)match_scores[ib] - match_scores[ia];
  if(diff == 0) {
    size_t xa = ia % score_width, xb = ib % score_width;
    if(xa != xb) return (xa > xb) - (xa < xb);
    return (ia > ib) - (ia < ib);
  }
  return diff > 0 ? 1 : -1;
}

// In-place O(n log n) heapsort for candidate hit indices (portable C11, no sort_r)
static void sw_sort_match_indices(size_t *arr, size_t n, const MatrixSort *ctx)
{
  if(n < 2) return;
  for(size_t start = n / 2; start > 0; start--) {
    size_t root = start - 1;
    while(root * 2 + 1 < n) {
      size_t child = root * 2 + 1;
      if(child + 1 < n && cmp_match_index(arr[child], arr[child + 1], ctx) < 0)
        child++;
      if(cmp_match_index(arr[root], arr[child], ctx) < 0) {
        size_t tmp = arr[root]; arr[root] = arr[child]; arr[child] = tmp;
        root = child;
      } else break;
    }
  }
  for(size_t end = n - 1; end > 0; end--) {
    size_t tmp = arr[end]; arr[end] = arr[0]; arr[0] = tmp;
    size_t root = 0;
    while(root * 2 + 1 < end) {
      size_t child = root * 2 + 1;
      if(child + 1 < end && cmp_match_index(arr[child], arr[child + 1], ctx) < 0)
        child++;
      if(cmp_match_index(arr[root], arr[child], ctx) < 0) {
        size_t sw = arr[root]; arr[root] = arr[child]; arr[child] = sw;
        root = child;
      } else break;
    }
  }
}

static bool _try_init_history(sw_history_t *hist)
{
  hist->hits_capacity = 256;
  size_t mem = hist->hits_capacity * sizeof(*(hist->sorted_match_indices));
  hist->sorted_match_indices = aln_malloc(mem);
  if(hist->sorted_match_indices == NULL) return false;
  if(!bitset_alloc(&hist->match_scores_mask, 256)) {
    free(hist->sorted_match_indices);
    hist->sorted_match_indices = NULL;
    return false;
  }
  return true;
}

sw_aligner_t* smith_waterman_try_new(void)
{
  sw_aligner_t *sw = aln_calloc(1, sizeof(sw_aligner_t));
  if(sw == NULL) return NULL;
  if(!_try_init_history(&sw->history)) {
    free(sw);
    return NULL;
  }
  sw->min_score = 1;
  return sw;
}

// Room for one bit per cell in the used-cell mask. Returns false if out of
// memory
static bool _ensure_mask_capacity(sw_history_t *hist, size_t arr_size)
{
  return arr_size <= hist->match_scores_mask.l ||
         bitset_set_length(&hist->match_scores_mask, ROUNDUP2POW(arr_size));
}

// Room for arr_size cells in the sorted index (only needed past the first
// hit). Returns false if out of memory
static bool _ensure_index_capacity(sw_history_t *hist, size_t arr_size)
{
  if(arr_size > hist->hits_capacity) {
    size_t capacity = ROUNDUP2POW(arr_size);
    size_t mem = capacity * sizeof(*(hist->sorted_match_indices));
    size_t *indices = NULL;
    if(capacity <= SIZE_MAX / sizeof(*(hist->sorted_match_indices)))
      indices = aln_realloc(hist->sorted_match_indices, mem);
    if(indices == NULL) return false;
    hist->sorted_match_indices = indices;
    hist->hits_capacity = capacity;
  }
  return true;
}

// After an error: no hits to fetch until aligned again
static void _no_hits(sw_history_t *hist)
{
  hist->num_of_hits = hist->next_hit = 0;
  hist->indexed = hist->first_done = true;
}

static aln_status_t _sw_align(const char *a, const char *b,
                              size_t len_a, size_t len_b,
                              const scoring_t *scoring, sw_aligner_t *sw);
static int _sw_fetch(sw_aligner_t *sw, alignment_t *result,
                     aln_status_t *status);

sw_aligner_t* smith_waterman_new(void)
{
  sw_aligner_t *sw = smith_waterman_try_new();
  if(sw == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  return sw;
}

void smith_waterman_free(sw_aligner_t *sw)
{
  free(sw->ckpt);
  free(sw->fast_res.result_a);
  free(sw->fast_res.result_b);
  aligner_destroy(&(sw->aligner));
  sub_table_free(&sw->table);
  bitset_dealloc(&sw->history.match_scores_mask);
  free(sw->history.sorted_match_indices);
  free(sw);
}

void smith_waterman_set_min_score(sw_aligner_t *sw, score_t min_score)
{
  sw->min_score = MAX2(min_score, 1);
}

void smith_waterman_set_compact(sw_aligner_t *sw, bool compact)
{
  sw->compact = compact;
}

void smith_waterman_set_options(sw_aligner_t *sw,
                                const alignment_options_t *opts)
{
  sw->aligner.options = opts;
}

aligner_t* smith_waterman_get_aligner(sw_aligner_t *sw)
{
  return &sw->aligner;
}

void smith_waterman_align(const char *a, const char *b,
                          const scoring_t *scoring, sw_aligner_t *sw)
{
  smith_waterman_align2(a, b, strlen(a), strlen(b), scoring, sw);
}

// The best hit of the region a[x0 .. x0+rx) against b[y0 .. y0+ry) (of the
// whole sequences) without its traceback bytes (see
// smith_waterman_strip_min_cells): one scan finds the best cell, keeping
// every strip_rows-th row; sw_follow_strips() walks the hit. want: the best
// score a rectangle's region must have, at its last cell (0: any). Sets
// *used where it applies (else nothing is changed); ALN_ERR_NOMEM if out
// of memory

aln_status_t smith_waterman_align_checked(const char *seq_a, const char *seq_b,
                                          size_t len_a, size_t len_b,
                                          const scoring_t *scoring,
                                          sw_aligner_t *sw, aln_error_t *err)
{
  // Whatever happens, the hits of the pair before are gone (they refer to
  // sequences the caller may have freed)
  _no_hits(&sw->history);
  aln_status_t status = alignment_check(scoring, seq_a, len_a, seq_b, len_b,
                                        true, err);
  if(status == ALN_OK) {
    status = _sw_align(seq_a, seq_b, len_a, len_b, scoring, sw);
    if(err != NULL) err->status = status;
  }
  return status;
}

void smith_waterman_align2(const char *a, const char *b,
                           size_t len_a, size_t len_b,
                           const scoring_t *scoring, sw_aligner_t *sw)
{
  aln_error_t e = {ALN_OK, 0, 0};
  e.status = _sw_align(a, b, len_a, len_b, scoring, sw);
  if(e.status != ALN_OK) {
    alignment_fprint_error(stderr, &e);
    exit(EXIT_FAILURE);
  }
}

// Returns ALN_ERR_TOO_LONG or ALN_ERR_OVERFLOW (see aligner_align()), or
// ALN_ERR_NOMEM; after an error there are no hits to fetch
static aln_status_t _sw_align(const char *a, const char *b,
                              size_t len_a, size_t len_b,
                              const scoring_t *scoring, sw_aligner_t *sw)
{
  aligner_t *aligner = &sw->aligner;
  sw_history_t *hist = &sw->history;
  sw->seq_a = a;
  sw->seq_b = b;
  sw->len_a = len_a;
  sw->len_b = len_b;
  sw->rect = sw->strips = sw->fast_hit = false;
  free(sw->ckpt);
  sw->ckpt = NULL;
  _no_hits(hist);
  sub_table_free(&sw->table);
  if(!sub_table_build(&sw->table, scoring, a, len_a, b, len_b))
    return ALN_ERR_NOMEM;
  aligner->compact = sw->compact;
  if(aligner->options != NULL) aligner->opts = *aligner->options;
  else alignment_options_default(&aligner->opts);
  if(sw->compact && aligner->opts.sw_fast_rect &&
     (double)len_a * (double)len_b >= (double)aligner->opts.sw_rect_min_cells &&
     sw_align_fast_rect(sw, a, b, len_a, len_b, scoring))
  {
    hist->indexed = hist->first_done = false;
    return ALN_OK;
  }
  if(!sw->compact ||
     (double)len_a * (double)len_b < (double)aligner->opts.sw_rect_min_cells ||
     !sw_align_rect(sw, a, b, len_a, len_b, scoring))
  {
    sw->rect = false;
    bool used = false;
    if(sw->compact && (double)len_a * (double)len_b >=
                      (double)aligner->opts.sw_strip_min_cells) {
      aln_status_t status = sw_strips_fill(sw, 0, 0, len_a, len_b, scoring, 0,
                                       &used);
      if(status != ALN_OK) return status;
    }
    if(!used) {
      aln_status_t status = aligner_align_checked(aligner, a, b, len_a, len_b,
                                                  scoring, 1);
      if(status != ALN_OK) return status;
    }
  }

  sw->fast_walk = !scoring->no_start_gap_penalty &&
                  !scoring->no_end_gap_penalty &&
                  !scoring->no_gaps_in_a && !scoring->no_gaps_in_b &&
                  !scoring->no_gaps_in_longest;

  // No cell is walked when the scan found no hit, or (strips) before a
  // second hit is fetched
  size_t arr_size = sw->rect && (aligner->best_score < sw->min_score ||
                                 sw->strips) ? 0 :
                    aligner->score_width * aligner->score_height;
  if(!_ensure_mask_capacity(hist, arr_size)) return ALN_ERR_NOMEM;

  // Clear the visited-cell mask for every cell (bits past them are never
  // read). Size is in bytes: clearing only one byte per 32-bit word left
  // stale bits that silently dropped hits on reuse
  memset(hist->match_scores_mask.b, 0,
         (arr_size + 31) / 32 * sizeof(hist->match_scores_mask.b[0]));

  // Set number of hits
  hist->num_of_hits = hist->next_hit = 0;
  hist->indexed = hist->first_done = false;
  return ALN_OK;
}

// The first hit in index order: the highest-scoring cell, then lowest
// position in seq_a, then seq_b. Its walk always succeeds, as no cell is used
// yet. Returns false if no cell scores at least min_score
static bool _best_cell(const sw_aligner_t *sw, size_t *best)
{
  const aligner_t *aligner = &sw->aligner;

  // A compact fill found it already
  if(aligner->traced) {
    *best = aligner->best_index;
    return aligner->best_score > 0 && aligner->best_score >= sw->min_score;
  }

  const score_t *scores = aligner->match_scores;
  const size_t w = aligner->score_width, h = aligner->score_height;
  size_t x, y, row, bx = w;
  score_t max_score = 0;

  for(row = 0; row < w * h; row++)
    if(scores[row] > max_score) max_score = scores[row];
  if(max_score < sw->min_score) return false;

  // Lowest x, then y (rows are scanned in order, so the first in a column)
  for(y = 0, row = 0; y < h; y++, row += w) {
    for(x = 0; x < bx; x++) {
      if(scores[row+x] == max_score) { bx = x; *best = row + x; break; }
    }
  }
  return true;
}

// Order every cell scoring at least min_score for smith_waterman_fetch().
// Returns false if out of memory
static bool _index_hits(sw_aligner_t *sw)
{
  aligner_t *aligner = &sw->aligner;
  sw_history_t *hist = &sw->history;
  const size_t arr_size = aligner->score_width * aligner->score_height;

  // Hits may start from cells scoring at least lo. Cells scoring less are
  // only reached after every higher one, so leaving them out changes no
  // hit that scores at least lo
  const score_t *scores = aligner->match_scores;
  const score_t lo = sw->min_score;
  size_t pos, x, y, w = aligner->score_width, h = aligner->score_height;
  size_t total = 0, count;
  score_t s, max_score = 0;
  for(pos = 0; pos < arr_size; pos++) {
    if(scores[pos] > max_score) max_score = scores[pos];
    total += (scores[pos] >= lo);
  }
  // Room for those cells only: one entry per cell of the matrix would be
  // 34 GB for a 50 kb pair
  if(!_ensure_index_capacity(hist, total)) return false;

  // Order cells by score (highest first), then position in seq_a, then seq_b.
  // Scores are bounded integers, so a counting sort replaces the comparison
  // sort; it also makes the order of tied hits the same on every platform.
  // With few candidates (e.g. a high min_score) a comparison sort of just
  // those is cheaper than the counting sort's column-major pass
  size_t *starts = NULL;
  if(total * 16 >= arr_size && (size_t)max_score <= arr_size)
    starts = aln_calloc((size_t)max_score + 2, sizeof(size_t));

  if(starts == NULL)
  {
    for(pos = 0; pos < arr_size; pos++) {
      if(scores[pos] >= lo)
        hist->sorted_match_indices[hist->num_of_hits++] = pos;
    }
    MatrixSort tmp_struct = {aligner->match_scores, aligner->score_width};
    sw_sort_match_indices(hist->sorted_match_indices, hist->num_of_hits,
                          &tmp_struct);
    return true;
  }

  // starts[s] = number of candidate cells scoring more than s
  for(pos = 0; pos < arr_size; pos++)
    if(scores[pos] >= lo) starts[scores[pos]]++;
  total = 0;
  for(s = max_score; s >= lo; s--) { count = starts[s]; starts[s] = total; total += count; }
  hist->num_of_hits = total;

  // Column-major pass so ties are ordered by x, then y
  for(x = 0; x < w; x++)
    for(y = 0, pos = x; y < h; y++, pos += w)
      if(scores[pos] >= lo)
        hist->sorted_match_indices[starts[scores[pos]]++] = pos;

  free(starts);
  return true;
}

// Return 1 if alignment was found, 0 otherwise, -1 if out of memory (its
// cells are then used, and the hit lost)
static int _follow_hit(sw_aligner_t* sw, size_t arr_index,
                       alignment_t* result)
{
  const aligner_t *aligner = &(sw->aligner);
  const sw_history_t *hist = &(sw->history);

  // Follow path through matrix
  size_t score_x = (size_t)ARR_2D_X(arr_index, aligner->score_width);
  size_t score_y = (size_t)ARR_2D_Y(arr_index, aligner->score_width);

  // Local alignments always (start and) end with a match
  enum Matrix curr_matrix = MATCH;
  score_t curr_score = aligner->traced ? aligner->best_score // the only hit
                                       : aligner->match_scores[arr_index];

  // Store end arr_index and (x,y) coords for later
  size_t end_arr_index = arr_index;
  size_t end_score_x = score_x;
  size_t end_score_y = score_y;
  score_t end_score = curr_score;

  size_t length;

  for(length = 0; ; length++)
  {
    // Stop at a cell already used by an earlier hit, else mark it as used
    uint32_t *word = &hist->match_scores_mask.b[arr_index >> 5];
    uint32_t bit = 1U << (arr_index & 31);
    if(*word & bit) return 0;
    *word |= bit;

    if(curr_score == 0) break;


    // The current matrix fixes which cell comes next, so check it before
    // working out the move: most rejected walks stop after a step or two
    size_t next = arr_index - (curr_matrix == GAP_B ? 1 :
                               curr_matrix == GAP_A ? aligner->score_width :
                                                      aligner->score_width + 1);
    if(bitset32_get(hist->match_scores_mask.b, next)) return 0;

    // Find out where to go next
    _sw_reverse_move(sw, &curr_matrix, &curr_score,
                     &score_x, &score_y, &arr_index);
  }

  // We got a result!
  // Allocate memory for the result
  if(!alignment_reserve(result, length)) return -1;
  result->length = length;

  // Jump back to the end of the alignment
  arr_index = end_arr_index;
  score_x = end_score_x;
  score_y = end_score_y;
  curr_matrix = MATCH;
  curr_score = end_score;

  // Now follow backwards again to create alignment!
  size_t i;

  for(i = length-1; curr_score > 0; i--)
  {
    switch(curr_matrix)
    {
      case MATCH:
        result->result_a[i] = aligner->seq_a[score_x-1];
        result->result_b[i] = aligner->seq_b[score_y-1];
        break;

      case GAP_A:
        result->result_a[i] = '-';
        result->result_b[i] = aligner->seq_b[score_y-1];
        break;

      case GAP_B:
        result->result_a[i] = aligner->seq_a[score_x-1];
        result->result_b[i] = '-';
        break;

      default:
      fprintf(stderr, "Program error: invalid matrix in _follow_hit()\n");
      fprintf(stderr, "Please submit a bug report to: https://github.com/DanielTillett/nucleics-align/issues\n");
      exit(EXIT_FAILURE);
    }

    _sw_reverse_move(sw, &curr_matrix, &curr_score,
                     &score_x, &score_y, &arr_index);
  }

  result->result_a[length] = '\0';
  result->result_b[length] = '\0';

  result->score = end_score;

  result->pos_a = score_x;
  result->pos_b = score_y;

  result->len_a = end_score_x - score_x;
  result->len_b = end_score_y - score_y;

  return 1;
}

int smith_waterman_fetch(sw_aligner_t *sw, alignment_t *result)
{
  aln_error_t e = {ALN_OK, 0, 0};
  int found = smith_waterman_fetch_checked(sw, result, &e);
  if(e.status != ALN_OK) {
    alignment_fprint_error(stderr, &e);
    exit(EXIT_FAILURE);
  }
  return found;
}

int smith_waterman_fetch_checked(sw_aligner_t *sw, alignment_t *result,
                                 aln_error_t *err)
{
  aln_status_t status = ALN_OK;
  int found = _sw_fetch(sw, result, &status);
  if(status != ALN_OK) _no_hits(&sw->history);
  if(err != NULL) {
    err->status = status;
    err->a = err->b = 0;
  }
  return found;
}

static int _sw_fetch(sw_aligner_t *sw, alignment_t *result,
                     aln_status_t *status)
{
  sw_history_t *hist = &(sw->history);
  int found;

  if(!hist->first_done) {
    hist->first_done = true;
    if(sw->fast_hit) {
      if(!alignment_reserve(result, sw->fast_res.length)) {
        *status = ALN_ERR_NOMEM;
        return 0;
      }
      memcpy(result->result_a, sw->fast_res.result_a, sw->fast_res.length + 1);
      memcpy(result->result_b, sw->fast_res.result_b, sw->fast_res.length + 1);
      result->length = sw->fast_res.length;
      result->pos_a = sw->fast_res.pos_a;
      result->pos_b = sw->fast_res.pos_b;
      result->len_a = sw->fast_res.len_a;
      result->len_b = sw->fast_res.len_b;
      result->score = sw->fast_res.score;
      return 1;
    }
    size_t best = 0;
    if(!_best_cell(sw, &best)) { hist->indexed = true; return 0; }
    // Cannot be refused: no cell is used yet
    if((sw->strips ? sw_follow_strips(sw, result)
                   : _follow_hit(sw, best, result)) < 0) {
      *status = ALN_ERR_NOMEM;
      return 0;
    }
    result->pos_a += sw->rect ? sw->off_a : 0;
    result->pos_b += sw->rect ? sw->off_b : 0;
    return 1;
  }

  // The best cell is first in the index and now used, so it is skipped
  if(!hist->indexed) {
    hist->indexed = true;
    aligner_t *aligner = &sw->aligner;
    if(sw->rect) {
      // Fill the whole matrices, and walk the first hit again in them to
      // mark the cells it used
      sw->rect = sw->strips = false;
      free(sw->ckpt);
      sw->ckpt = NULL;
      aligner->compact = false;
      *status = aligner_align_checked(aligner, sw->seq_a, sw->seq_b,
                                      sw->len_a, sw->len_b,
                                      aligner->scoring, 1);
      aligner->compact = true;
      if(*status != ALN_OK) return 0;
      size_t arr_size = aligner->score_width * aligner->score_height;
      if(!_ensure_mask_capacity(hist, arr_size)) {
        *status = ALN_ERR_NOMEM;
        return 0;
      }
      memset(hist->match_scores_mask.b, 0,
             (arr_size + 31) / 32 * sizeof(hist->match_scores_mask.b[0]));
      if(sw->fast_hit) {
        // Mark the exact cells of the fast_res hit already returned to the caller
        size_t w = aligner->score_width;
        size_t x = sw->fast_res.pos_a, y = sw->fast_res.pos_b;
        size_t idx = y * w + x;
        hist->match_scores_mask.b[idx >> 5] |= 1U << (idx & 31);
        for(size_t k = 0; k < sw->fast_res.length; k++) {
          if(sw->fast_res.result_a[k] != '-') x++;
          if(sw->fast_res.result_b[k] != '-') y++;
          idx = y * w + x;
          hist->match_scores_mask.b[idx >> 5] |= 1U << (idx & 31);
        }
        sw->fast_hit = false;
      } else {
        alignment_t first;
        memset(&first, 0, sizeof(first));
        found = _follow_hit(sw, sw->full_best, &first);
        free(first.result_a);
        free(first.result_b);
        if(found < 0) {
          *status = ALN_ERR_NOMEM;
          return 0;
        }
      }
    }
    else if(aligner->traced) {
      // Later hits need the score matrices: fill them (the used cells stay
      // marked)
      aligner->compact = false;
      *status = aligner_align_checked(aligner, aligner->seq_a, aligner->seq_b,
                                      aligner->score_width-1,
                                      aligner->score_height-1,
                                      aligner->scoring, 1);
      aligner->compact = true;
      if(*status != ALN_OK) return 0;
    }
    if(!_index_hits(sw)) {
      *status = ALN_ERR_NOMEM;
      return 0;
    }
  }

  while(hist->next_hit < hist->num_of_hits)
  {
    size_t arr_index = hist->sorted_match_indices[hist->next_hit++];

    if(!bitset32_get(hist->match_scores_mask.b, arr_index)) {
      found = _follow_hit(sw, arr_index, result);
      if(found < 0) {
        *status = ALN_ERR_NOMEM;
        return 0;
      }
      if(found) return 1;
    }
  }

  return 0;
}
