/*
 smith_waterman.c
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */

// Turn on debugging output by defining NUCLEICS_ALIGN_VERBOSE
//#define NUCLEICS_ALIGN_VERBOSE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sort_r/sort_r.h"

#include "smith_waterman.h"
#include "hirschberg.h"
#include "alignment_macros.h"
#include "alignment_alloc.h"
#include "alignment_testing.h"

typedef struct {
  uint32_t *b; size_t l, s; // l is bits, s in uint32_t
} BitSet;

static inline BitSet* bitset_alloc(BitSet *bs, size_t l) {
  bs->s = l/32 + (l%32 != 0); // not (l+31)/32, which wraps for huge l
  bs->l = l;
  bs->b = calloc(bs->s, sizeof(bs->b[0]));
  return bs->b ? bs : NULL;
}

static inline void bitset_dealloc(BitSet *bs) {
  free(bs->b);
  memset(bs, 0, sizeof(*bs));
}

static inline BitSet* bitset_set_length(BitSet *bs, size_t l) {
  size_t ss = l/32 + (l%32 != 0);
  if(ss > bs->s) {
    uint32_t *b = aln_realloc(bs->b, ss*sizeof(bs->b[0]));
    if(b == NULL) return NULL;
    memset(b+bs->s, 0, (ss-bs->s)*sizeof(b[0])); // zero new memory
    bs->b = b;
    bs->s = ss;
  }
  bs->l = l;
  return bs;
}

// For iterating through local alignments. The sorted index of hits is only
// built if a second hit is fetched: the first is just the best cell
typedef struct
{
  BitSet match_scores_mask;
  size_t *sorted_match_indices, hits_capacity, num_of_hits, next_hit;
  bool indexed, first_done;
} sw_history_t;

// Store alignment here
struct sw_aligner_t
{
  aligner_t aligner;
  sw_history_t history;
  // For following hits back through the matrices without scoring_lookup()
  sub_table_t table;
  bool fast_walk;
  score_t min_score; // lowest score a hit may start from, at least 1
  bool compact; // fill traceback bytes only, until a second hit is fetched
  // Large compact alignments (see _sw_align_rect()): the aligner holds only
  // the rectangle of the best hits, which starts at off_a, off_b
  bool rect;
  size_t off_a, off_b, full_best; // full_best: best cell in the whole matrix
  const char *seq_a, *seq_b; // the whole sequences
  size_t len_a, len_b;
  // Strip alignments (see _sw_strips()): rect is set too, and the aligner
  // describes the region but holds no traceback bytes; ckpt holds every
  // strip_rows-th row of its matrices
  bool strips;
  score_t *ckpt;
  size_t strip_rows;
  // Opt-in fast best-hit reconstruction (smith_waterman_fast_rect)
  bool fast_hit;
  alignment_t fast_res;
};

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

// Function passed to sort_r
static int sort_match_indices(const void *aa, const void *bb, void *arg)
{
  const size_t *a = aa;
  const size_t *b = bb;
  const MatrixSort *tmp = arg;

  // Recover variables from the struct
  const score_t *match_scores = tmp->match_scores;
  size_t score_width = tmp->score_width;

  long diff = (long)match_scores[*b] - match_scores[*a];

  // Then by position on seq_a, then seq_b (same order as the counting sort)
  if(diff == 0) {
    size_t xa = *a % score_width, xb = *b % score_width;
    if(xa != xb) return (xa > xb) - (xa < xb);
    return (*a > *b) - (*a < *b);
  }
  else return diff > 0 ? 1 : -1;
}

static void _init_history(sw_history_t *hist)
{
  hist->hits_capacity = 256;
  size_t mem = hist->hits_capacity * sizeof(*(hist->sorted_match_indices));
  hist->sorted_match_indices = malloc(mem);
  if(!bitset_alloc(&hist->match_scores_mask, 256) || !hist->sorted_match_indices) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

// Room for one bit per cell in the used-cell mask. Returns false if out of
// memory
static bool _ensure_mask_capacity(sw_history_t *hist, size_t arr_size)
{
  return arr_size <= hist->match_scores_mask.l ||
         bitset_set_length(&hist->match_scores_mask, ROUNDUP2POW(arr_size));
}

// Room for every cell in the sorted index (only needed past the first hit).
// Returns false if out of memory
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
  sw_aligner_t *sw = calloc(1, sizeof(sw_aligner_t));
  if(sw == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  _init_history(&(sw->history));
  sw->min_score = 1;
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
// every strip_rows-th row; _follow_strips() walks the hit. want: the best
// score a rectangle's region must have, at its last cell (0: any). Sets
// *used where it applies (else nothing is changed); ALN_ERR_NOMEM if out
// of memory
static aln_status_t _sw_strips(sw_aligner_t *sw, size_t x0, size_t y0,
                               size_t rx, size_t ry, const scoring_t *scoring,
                               score_t want, bool *used)
{
  aligner_t *aligner = &sw->aligner;
  *used = false;
  // Rows per strip near sqrt(12 x rows): the strip's bytes (k bytes per
  // column) and the rows kept (12 bytes per column each) then take about
  // the same memory
  size_t k = smith_waterman_strip_rows;
  if(k == 0) {
    for(k = 16; k * k < 12 * ry; k++) {}
  }
  k = MIN2(k, ry);
  const size_t kept = ry / k + 1;
  if(rx == 0 || ry == 0 || kept > SIZE_MAX / 3 / sizeof(score_t) / (rx + 1))
    return ALN_OK;
  score_t *ckpt = aln_malloc(kept * 3 * (rx + 1) * sizeof(score_t));
  if(ckpt == NULL) return ALN_ERR_NOMEM;
  score_t best;
  size_t bx, by;
  if(!alignment_sw_rows(sw->seq_a + x0, rx, sw->seq_b + y0, ry, scoring,
                        &aligner->opts, k, ckpt, &best, &bx, &by) ||
     (want > 0 && (best != want || bx != rx || by != ry))) {
    free(ckpt);
    return ALN_OK;
  }
  free(sw->ckpt);
  sw->ckpt = ckpt;
  sw->strip_rows = k;
  aligner->scoring = scoring;
  aligner->seq_a = sw->seq_a + x0;
  aligner->seq_b = sw->seq_b + y0;
  aligner->score_width = rx + 1;
  aligner->score_height = ry + 1;
  aligner->traced = true;
  aligner->trace_banded = false;
  aligner->best_score = best;
  aligner->best_index = by * (rx + 1) + bx;
  sw->strips = sw->rect = true;
  sw->off_a = x0;
  sw->off_b = y0;
  sw->full_best = (y0 + by) * (sw->len_a + 1) + x0 + bx;
  __atomic_add_fetch(&smith_waterman_strip_alignments, 1, __ATOMIC_RELAXED);
  *used = true;
  return ALN_OK;
}

// Follow the best hit of a strip alignment back from its last cell, as
// _follow_hit() follows traceback bytes: from the bottom, each strip of
// rows gets its bytes again (up to the hit's column) when the walk enters
// it. Returns 1, or -1 if out of memory
static int _follow_strips(sw_aligner_t *sw, alignment_t *result)
{
  const aligner_t *aligner = &sw->aligner;
  const scoring_t *scoring = aligner->scoring;
  const sub_table_t *t = &sw->table;
  const score_t open = scoring->gap_open + scoring->gap_extend;
  const score_t ext = scoring->gap_extend;
  const size_t w = aligner->score_width, k = sw->strip_rows;
  const char *a = aligner->seq_a, *b = aligner->seq_b;
  const size_t end_x = aligner->best_index % w, end_y = aligner->best_index / w;
  size_t x = end_x, y = end_y, n = 0, s0 = 0, xs = 0;
  enum Matrix m = MATCH;
  score_t score = aligner->best_score;

  uint8_t *ops = aln_malloc(end_x + end_y + 1); // the matrix of each step
  uint8_t *trace = aln_malloc(k * end_x);
  if(ops == NULL || trace == NULL) {
    free(ops);
    free(trace);
    return -1;
  }

  while(score > 0) {
    if(n == 0 || y < s0) {
      // The strip holding row y, from the row kept above it
      const size_t kept = (y - 1) / k;
      s0 = kept * k + 1;
      xs = x;
      if(!alignment_sw_trace_rows(a, w - 1, b, aligner->score_height - 1,
                                  scoring, t, &aligner->opts,
                                  aligner->best_score,
                                  sw->ckpt + kept * 3 * w, s0, y, xs, trace)) {
        free(ops);
        free(trace);
        return -1;
      }
    }
    // As _sw_reverse_move() steps on traceback bytes
    int next = (trace[(y - s0) * xs + x - 1] >> (2 * (int)m)) & 3;
    score_t step;
    if(next == ALN_TRACE_NONE) {
      fprintf(stderr, "Program error: traceback fail (_follow_strips)\n");
      exit(EXIT_FAILURE);
    }
    if(m == MATCH)
      step = t->scores[t->code[(unsigned char)b[y-1]] * t->n +
                       t->code[(unsigned char)a[x-1]]];
    else if(m == GAP_A) step = next == GAP_A ? ext : open;
    else step = next == GAP_B ? ext : open;
    ops[n++] = (uint8_t)m;
    if(m != GAP_A) x--;
    if(m != GAP_B) y--;
    m = (enum Matrix)next;
    score -= step;
  }
  free(trace);

  if(!alignment_reserve(result, n)) {
    free(ops);
    return -1;
  }
  // From the first cell: the steps in reverse
  size_t i, xa = x, yb = y;
  for(i = 0; i < n; i++) {
    const enum Matrix op = (enum Matrix)ops[n-1-i];
    if(op != GAP_A) xa++;
    if(op != GAP_B) yb++;
    result->result_a[i] = op == GAP_A ? '-' : a[xa-1];
    result->result_b[i] = op == GAP_B ? '-' : b[yb-1];
  }
  free(ops);
  result->result_a[n] = result->result_b[n] = '\0';
  result->length = n;
  result->score = aligner->best_score;
  result->pos_a = x;
  result->pos_b = y;
  result->len_a = end_x - x;
  result->len_b = end_y - y;
  return 1;
}

// Best hit of a large compact alignment in little memory (as SSW finds it,
// but with the same hit as the full search). A scan finds the best score
// and its first cell (bx, by). A local path ending there and scoring best
// is a best path of the reversed prefixes, so a scan of a[0..bx) and
// b[0..by) reversed finds where every best path ending at (bx, by) starts:
// the rectangle from the lowest start column and row to (bx, by) holds them
// all. A fill of just that rectangle gives exact values along every best
// path, and no cell outside it could tie with one on them, so the
// traceback takes the same steps. The rectangle is as big as the hit: this
// saves memory and time when hits are short against the sequences.
// That needs every hit to start with a pair: the traceback stops at a score
// of 0, so a hit can open with a gap only if a gap can score above 0.
// Returns false where the scans do not apply, or if out of memory (the
// whole fill then reports it)
// Verify that a global alignment of a[x0..x0+rx) vs b[y0..y0+ry) is a valid
// optimal local hit scoring exactly best with strictly positive local prefixes
static bool _sw_verify_local_hit(const sub_table_t *t, const scoring_t *scoring,
                                 const alignment_t *res, size_t rx, size_t ry,
                                 score_t best)
{
  if(res->length == 0 || res->score != best) return false;
  if(res->result_a[0] == '-' || res->result_b[0] == '-') return false;
  if(res->result_a[res->length - 1] == '-' ||
     res->result_b[res->length - 1] == '-') return false;
  const long long open = (long long)scoring->gap_open + scoring->gap_extend;
  const long long ext = (long long)scoring->gap_extend;
  long long sc = 0;
  size_t ia = 0, ib = 0, k;
  int st = 0; // 0: match, 1: gap_a, 2: gap_b
  for(k = 0; k < res->length; k++) {
    char ca = res->result_a[k], cb = res->result_b[k];
    if(ca == '-' && cb == '-') return false;
    if(ca == '-') {
      sc += (st == 1 ? ext : open);
      st = 1;
      ib++;
    } else if(cb == '-') {
      sc += (st == 2 ? ext : open);
      st = 2;
      ia++;
    } else {
      sc += t->scores[t->code[(unsigned char)cb] * t->n +
                      t->code[(unsigned char)ca]];
      st = 0;
      ia++;
      ib++;
    }
    if(sc <= 0) return false;
  }
  return ia == rx && ib == ry && sc == (long long)best;
}

static bool _sw_align_fast_rect(sw_aligner_t *sw, const char *a, const char *b,
                                size_t len_a, size_t len_b,
                                const scoring_t *scoring)
{
  aligner_t *aligner = &sw->aligner;
  score_t best = 0, best2 = 0;
  size_t bx = 0, by = 0, rx = 0, ry = 0, i;

  if(scoring->gap_extend > 0 || scoring->gap_open + scoring->gap_extend > 0 ||
     scoring->no_start_gap_penalty || scoring->no_end_gap_penalty ||
     scoring->no_gaps_in_a || scoring->no_gaps_in_b ||
     scoring->no_gaps_in_longest || scoring->no_mismatches)
    return false;

  if(!alignment_sw_scan(a, len_a, b, len_b, scoring, 0, &aligner->opts,
                        &best, &bx, &by))
    return false;

  if(best < sw->min_score) {
    aligner->scoring = scoring;
    aligner->seq_a = a;
    aligner->seq_b = b;
    aligner->score_width = len_a + 1;
    aligner->score_height = len_b + 1;
    aligner->traced = true;
    aligner->best_score = best;
    sw->rect = true;
    sw->off_a = sw->off_b = 0;
    return true;
  }

  // Reverse scan with target=0 locates a genuine jointly-located optimal start
  char *ra = aln_malloc(bx + by);
  if(ra == NULL) return false;
  char *rb = ra + bx;
  for(i = 0; i < bx; i++) ra[i] = a[bx - 1 - i];
  for(i = 0; i < by; i++) rb[i] = b[by - 1 - i];
  bool ok = alignment_sw_scan(ra, bx, rb, by, scoring, 0, &aligner->opts,
                              &best2, &rx, &ry);
  free(ra);
  if(!ok || best2 != best || rx == 0 || ry == 0) return false;

  const size_t x0 = bx - rx, y0 = by - ry;
  const size_t total = rx <= SIZE_MAX - ry ? rx + ry : SIZE_MAX;
  size_t wfa_budget = alignment_wfa_budget;
  if(wfa_budget == 0) {
    const size_t max_c = (size_t)1 << 28;
    wfa_budget = total <= max_c / 256 ? 256 * total : max_c;
  }

  aln_status_t st = alignment_wfa_align(&sw->table, scoring, a + x0, rx,
                                        b + y0, ry, wfa_budget, &sw->fast_res);
  if(st != ALN_OK) {
    st = hirschberg_align_checked(a + x0, b + y0, rx, ry, scoring,
                                  &aligner->opts, &sw->fast_res, NULL);
  }
  if(st != ALN_OK ||
     !_sw_verify_local_hit(&sw->table, scoring, &sw->fast_res, rx, ry, best))
    return false;

  sw->fast_res.pos_a = x0;
  sw->fast_res.pos_b = y0;
  sw->fast_res.len_a = rx;
  sw->fast_res.len_b = ry;
  sw->fast_hit = true;
  sw->rect = true;
  sw->off_a = x0;
  sw->off_b = y0;
  sw->full_best = by * (len_a + 1) + bx;
  aligner->scoring = scoring;
  aligner->seq_a = a;
  aligner->seq_b = b;
  aligner->score_width = len_a + 1;
  aligner->score_height = len_b + 1;
  aligner->traced = true;
  aligner->best_score = best;
  return true;
}

static bool _sw_align_rect(sw_aligner_t *sw, const char *a, const char *b,
                           size_t len_a, size_t len_b,
                           const scoring_t *scoring)
{
  aligner_t *aligner = &sw->aligner;
  score_t best, best2;
  size_t bx, by, rx, ry, i;

  // Only where no gap scores above 0, and where the full matrices would fit
  // (aligner_align() reports sequences too long for them)
  if(scoring->gap_extend > 0 || scoring->gap_open + scoring->gap_extend > 0 ||
     len_b+1 > SIZE_MAX / (len_a+1) / sizeof(score_t))
    return false;

  // Gives up on a hit that spans the sequences: the whole fill is quicker
  if(!alignment_sw_scan_short(a, len_a, b, len_b, scoring, &aligner->opts,
                              &best, &bx, &by))
    return false;

  if(best < sw->min_score) {
    // No hit: nothing to fill, but describe the sequences as a fill would
    aligner->scoring = scoring;
    aligner->seq_a = a;
    aligner->seq_b = b;
    aligner->score_width = len_a+1;
    aligner->score_height = len_b+1;
    aligner->traced = true;
    aligner->best_score = best;
    sw->rect = true;
    sw->off_a = sw->off_b = 0;
    return true;
  }

  // A hit scoring best has at least best / smax aligned pairs. If that is
  // most of the prefix up to (bx, by), fill the whole prefix (exact too:
  // nothing after (bx, by) changes it) rather than scan for the rectangle
  const sub_table_t *t = &sw->table;
  score_t smax = 0;
  for(i = 0; i < t->n * t->n; i++) smax = MAX2(smax, t->scores[i]);
  if(smax > 0 && (double)best >= 0.5 * smax * (double)MIN2(bx, by)) {
    rx = bx;
    ry = by;
  }
  else
  {
    // Reversed prefixes up to the end of the best hit
    char *ra = aln_malloc(bx + by);
    if(ra == NULL) return false;
    char *rb = ra + bx;
    for(i = 0; i < bx; i++) ra[i] = a[bx-1-i];
    for(i = 0; i < by; i++) rb[i] = b[by-1-i];
    bool ok = alignment_sw_scan(ra, bx, rb, by, scoring, best, &aligner->opts,
                                &best2, &rx, &ry);
    free(ra);
    if(!ok || rx == 0 || ry == 0) return false;
  }

  // Reversed cell (rx, ry) is the start pair (bx-rx+1, by-ry+1)
  const size_t x0 = bx - rx, y0 = by - ry;
  if((double)rx * (double)ry >= (double)aligner->opts.sw_strip_min_cells) {
    bool used;
    if(_sw_strips(sw, x0, y0, rx, ry, scoring, best, &used) != ALN_OK)
      return false; // the whole fill reports running out of memory
    if(used) return true;
  }
  if(aligner_align_checked(aligner, a + x0, b + y0, rx, ry, scoring,
                           1) != ALN_OK ||
     !aligner->traced || aligner->best_score != best ||
     aligner->best_index != ry * (rx + 1) + rx) return false;

  sw->rect = true;
  sw->off_a = x0;
  sw->off_b = y0;
  sw->full_best = by * (len_a + 1) + bx;
  return true;
}

aln_status_t smith_waterman_align_checked(const char *seq_a, const char *seq_b,
                                          size_t len_a, size_t len_b,
                                          const scoring_t *scoring,
                                          sw_aligner_t *sw, aln_error_t *err)
{
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
  if(sw->compact && smith_waterman_fast_rect &&
     (double)len_a * (double)len_b >= (double)aligner->opts.sw_rect_min_cells &&
     _sw_align_fast_rect(sw, a, b, len_a, len_b, scoring))
  {
    hist->indexed = hist->first_done = false;
    return ALN_OK;
  }
  if(!sw->compact ||
     (double)len_a * (double)len_b < (double)aligner->opts.sw_rect_min_cells ||
     !_sw_align_rect(sw, a, b, len_a, len_b, scoring))
  {
    sw->rect = false;
    bool used = false;
    if(sw->compact && (double)len_a * (double)len_b >=
                      (double)aligner->opts.sw_strip_min_cells) {
      aln_status_t status = _sw_strips(sw, 0, 0, len_a, len_b, scoring, 0,
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
  if(!_ensure_index_capacity(hist, arr_size)) return false;

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
    sort_r(hist->sorted_match_indices, hist->num_of_hits,
           sizeof(size_t), sort_match_indices, &tmp_struct);
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
    if((sw->strips ? _follow_strips(sw, result)
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
