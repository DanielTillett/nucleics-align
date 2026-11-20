/*
 needleman_wunsch.c
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "needleman_wunsch.h"
#include "hirschberg.h"

nw_aligner_t* needleman_wunsch_new(void)
{
  nw_aligner_t *nw = calloc(1, sizeof(nw_aligner_t));
  if(nw == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  return nw;
}

void needleman_wunsch_free(nw_aligner_t *nw)
{
  aligner_destroy(nw);
  free(nw);
}

void needleman_wunsch_align(const char *a, const char *b,
                            const scoring_t *scoring,
                            nw_aligner_t *nw, alignment_t *result)
{
  needleman_wunsch_align2(a, b, strlen(a), strlen(b), scoring, nw, result);
}

// Gaps before the first aligned pair, then shift the alignment (written
// backwards from next_char) to the start of the buffers
static aln_status_t _nw_align(const char *a, const char *b,
                              size_t len_a, size_t len_b,
                              const scoring_t *scoring,
                              nw_aligner_t *nw, alignment_t *result);

// Too big for the matrices (status): align in linear space instead, if the
// options ask (alignment_linear_fallback) and it supports the scoring.
// Returns true if that worked
static bool _nw_linear(const char *a, const char *b, size_t len_a,
                       size_t len_b, const scoring_t *scoring,
                       const nw_aligner_t *nw, alignment_t *result,
                       aln_status_t status)
{
  alignment_options_t defaults;
  const alignment_options_t *opts = nw->options;
  if(opts == NULL) {
    alignment_options_default(&defaults);
    opts = &defaults;
  }
  return (status == ALN_ERR_TOO_LONG || status == ALN_ERR_NOMEM) &&
         opts->linear_fallback &&
         hirschberg_align_checked(a, b, len_a, len_b, scoring, opts, result,
                                  NULL) == ALN_OK;
}

static void _nw_finish(const nw_aligner_t *nw, alignment_t *result,
                       size_t score_x, size_t score_y, size_t next_char,
                       size_t longest_alignment)
{
  char *alignment_a = result->result_a, *alignment_b = result->result_b;

  // Gap in A
  while(score_y > 0)
  {
    alignment_a[next_char] = '-';
    alignment_b[next_char] = nw->seq_b[score_y-1];
    next_char--;
    score_y--;
  }

  // Gap in B
  while(score_x > 0)
  {
    alignment_a[next_char] = nw->seq_a[score_x-1];
    alignment_b[next_char] = '-';
    next_char--;
    score_x--;
  }

  // Shift alignment strings back into 0th position in char arrays
  size_t first_char = next_char+1;
  size_t alignment_len = longest_alignment - first_char;

  // Use memmove
  memmove(alignment_a, alignment_a+first_char, alignment_len);
  memmove(alignment_b, alignment_b+first_char, alignment_len);

  alignment_a[alignment_len] = '\0';
  alignment_b[alignment_len] = '\0';

  result->length = alignment_len;
  // A global alignment covers both sequences (a result reused from
  // Smith-Waterman would otherwise keep its hit's positions)
  result->pos_a = result->pos_b = 0;
  result->len_a = nw->score_width - 1;
  result->len_b = nw->score_height - 1;
}

// The traceback below, from the traceback bytes of a compact fill. Returns
// false at a step with no valid predecessor, for the matrices to report
// Returns 1, or 0 if a step has no predecessor (the caller then fills the
// matrices, to report it), or -1 if out of memory (filling a band's strip
// again, see alignment_trace_load())
static int _nw_follow_trace(nw_aligner_t *nw, alignment_t *result)
{
  const size_t w = nw->score_width;
  size_t longest_alignment = w-1 + nw->score_height-1;
  size_t next_char = longest_alignment-1;

  enum Matrix curr_matrix = MATCH;
  score_t curr_score = nw->final_scores[MATCH];
  if(nw->final_scores[GAP_B] >= curr_score) {
    curr_matrix = GAP_B;
    curr_score = nw->final_scores[GAP_B];
  }
  if(nw->final_scores[GAP_A] >= curr_score) {
    curr_matrix = GAP_A;
    curr_score = nw->final_scores[GAP_A];
  }
  result->score = curr_score;

  size_t score_x = w-1, score_y = nw->score_height-1;
  for(; score_x > 0 && score_y > 0; next_char--)
  {
    result->result_a[next_char] = curr_matrix == GAP_A ? '-' : nw->seq_a[score_x-1];
    result->result_b[next_char] = curr_matrix == GAP_B ? '-' : nw->seq_b[score_y-1];

    if(!alignment_trace_load(nw, score_x, score_y)) return -1;
    int next = alignment_trace_next(nw,
                                    alignment_trace_index(nw, score_x, score_y),
                                    curr_matrix);
    if(next == ALN_TRACE_NONE) return 0;
    if(curr_matrix != GAP_A) score_x--;
    if(curr_matrix != GAP_B) score_y--;
    curr_matrix = (enum Matrix)next;
  }

  _nw_finish(nw, result, score_x, score_y, next_char, longest_alignment);
  return 1;
}

aln_status_t needleman_wunsch_align_checked(const char *a, const char *b,
                                            size_t len_a, size_t len_b,
                                            const scoring_t *scoring,
                                            nw_aligner_t *nw,
                                            alignment_t *result,
                                            aln_error_t *err)
{
  aln_status_t status = alignment_check(scoring, a, len_a, b, len_b, false,
                                        err);
  if(status == ALN_OK) {
    status = _nw_align(a, b, len_a, len_b, scoring, nw, result);
    if(err != NULL) err->status = status;
  }
  if(_nw_linear(a, b, len_a, len_b, scoring, nw, result, status)) {
    status = ALN_OK;
    if(err != NULL) err->status = status;
  }
  return status;
}

void needleman_wunsch_align2(const char *a, const char *b,
                             size_t len_a, size_t len_b,
                             const scoring_t *scoring,
                             nw_aligner_t *nw, alignment_t *result)
{
  aln_error_t e = {ALN_OK, 0, 0};
  e.status = _nw_align(a, b, len_a, len_b, scoring, nw, result);
  if(_nw_linear(a, b, len_a, len_b, scoring, nw, result, e.status))
    e.status = ALN_OK;
  if(e.status != ALN_OK) {
    alignment_fprint_error(stderr, &e);
    exit(EXIT_FAILURE);
  }
}

// Returns ALN_ERR_TOO_LONG or ALN_ERR_OVERFLOW (see aligner_align()), or
// ALN_ERR_NOMEM with no alignment in result
static aln_status_t _nw_align(const char *a, const char *b,
                              size_t len_a, size_t len_b,
                              const scoring_t *scoring,
                              nw_aligner_t *nw, alignment_t *result)
{
  aln_status_t status;
  #ifdef SEQ_ALIGN_VERBOSE
    // The diagnostics below need the score matrices
    bool compact = nw->compact;
    nw->compact = false;
    status = aligner_align_checked(nw, a, b, len_a, len_b, scoring, 0);
    nw->compact = compact;
  #else
    status = aligner_align_checked(nw, a, b, len_a, len_b, scoring, 0);
  #endif
  if(status != ALN_OK) return status;

  // note: longest_alignment = strlen(seq_a) + strlen(seq_b)
  size_t longest_alignment = nw->score_width-1 + nw->score_height-1;
  if(!alignment_reserve(result, longest_alignment)) return ALN_ERR_NOMEM;

  if(nw->traced) {
    const int followed = _nw_follow_trace(nw, result);
    if(followed > 0) return ALN_OK;
    if(followed < 0) return ALN_ERR_NOMEM;
    // Fill the matrices, so the traceback below reports the failed step
    nw->compact = false;
    status = aligner_align_checked(nw, a, b, len_a, len_b, scoring, 0);
    nw->compact = true;
    if(status != ALN_OK) return status;
  }

  // work backwards re-tracing optimal alignment, then shift sequences into place

  // Position of next alignment character in buffer (working backwards)
  size_t next_char = longest_alignment-1;

  size_t arr_size = nw->score_width * nw->score_height;

  // Get max score (and therefore current matrix)
  enum Matrix curr_matrix = MATCH;
  score_t curr_score = nw->match_scores[arr_size-1];

  if(nw->gap_b_scores[arr_size-1] >= curr_score)
  {
    curr_matrix = GAP_B;
    curr_score = nw->gap_b_scores[arr_size-1];
  }

  if(nw->gap_a_scores[arr_size-1] >= curr_score)
  {
    curr_matrix = GAP_A;
    curr_score = nw->gap_a_scores[arr_size-1];
  }

  #ifdef SEQ_ALIGN_VERBOSE
    alignment_print_matrices(nw);
  #endif

  result->score = curr_score;
  char *alignment_a = result->result_a, *alignment_b = result->result_b;

  // coords in score matrices
  size_t score_x = nw->score_width-1, score_y = nw->score_height-1;
  size_t arr_index = arr_size - 1;

  for(; score_x > 0 && score_y > 0; next_char--)
  {
    #ifdef SEQ_ALIGN_VERBOSE
    printf("matrix: %s (%zu,%zu) score: %i\n",
           MATRIX_NAME(curr_matrix), score_x-1, score_y-1, curr_score);
    #endif

    switch(curr_matrix)
    {
      case MATCH:
        alignment_a[next_char] = nw->seq_a[score_x-1];
        alignment_b[next_char] = nw->seq_b[score_y-1];
        break;

      case GAP_A:
        alignment_a[next_char] = '-';
        alignment_b[next_char] = nw->seq_b[score_y-1];
        break;

      case GAP_B:
        alignment_a[next_char] = nw->seq_a[score_x-1];
        alignment_b[next_char] = '-';
        break;

      default:
        fprintf(stderr, "Program error: invalid matrix number\n");
        fprintf(stderr, "Please submit a bug report to: https://github.com/DanielTillett/nucleics-align/issues\n");
        exit(EXIT_FAILURE);
    }

    if(score_x > 0 && score_y > 0)
    {
      alignment_reverse_move(&curr_matrix, &curr_score,
                             &score_x, &score_y, &arr_index, nw);
    }
  }

  _nw_finish(nw, result, score_x, score_y, next_char, longest_alignment);
  return ALN_OK;
}
