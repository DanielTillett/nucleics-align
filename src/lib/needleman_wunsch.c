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

#include "alignment_macros.h"
#include "alignment_alloc.h"
#include "needleman_wunsch.h"
#include "hirschberg.h"

nw_aligner_t* needleman_wunsch_try_new(void)
{
  return aln_calloc(1, sizeof(nw_aligner_t));
}

nw_aligner_t* needleman_wunsch_new(void)
{
  nw_aligner_t *nw = needleman_wunsch_try_new();
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

// The alignment with the full matrices (compact: traceback bytes)
static aln_status_t _nw_align(const char *a, const char *b,
                              size_t len_a, size_t len_b,
                              const scoring_t *scoring,
                              nw_aligner_t *nw, alignment_t *result);

// The aligner's settings: its own, or the defaults (into *defaults)
static const alignment_options_t *_nw_opts(const nw_aligner_t *nw,
                                           alignment_options_t *defaults)
{
  if(nw->options != NULL) return nw->options;
  alignment_options_default(defaults);
  return defaults;
}

// Too big for the matrices (status): align in linear space instead, if the
// options ask (alignment_linear_fallback) and it supports the scoring.
// Returns true if that worked
static bool _nw_linear(const char *a, const char *b, size_t len_a,
                       size_t len_b, const scoring_t *scoring,
                       const nw_aligner_t *nw, alignment_t *result,
                       aln_status_t status)
{
  alignment_options_t defaults;
  const alignment_options_t *opts = _nw_opts(nw, &defaults);
  return (status == ALN_ERR_TOO_LONG || status == ALN_ERR_NOMEM) &&
         opts->linear_fallback &&
         hirschberg_align_checked(a, b, len_a, len_b, scoring, opts, result,
                                  NULL) == ALN_OK;
}

// By the wavefront algorithm, if the options ask (alignment_wfa_alignment)
// and it applies: the scoring suits it and the pair is similar enough.
// ALN_OK if that worked, ALN_ERR_NOMEM if it ran out of memory (the
// matrices need more, so there is no point trying them), otherwise
// ALN_ERR_UNSUPPORTED: the other aligners align it
static aln_status_t _nw_wfa(const char *a, const char *b, size_t len_a,
                            size_t len_b, const scoring_t *scoring,
                            const nw_aligner_t *nw, alignment_t *result,
                            bool checked)
{
  alignment_options_t defaults;
  const alignment_options_t *opts = _nw_opts(nw, &defaults);
  if(!opts->wfa_align) return ALN_ERR_UNSUPPORTED;
  // (a pair the scoring lacks is left to the matrices, which report it;
  // the checked entry point has already scanned for one)
  if(!checked && !alignment_scoring_complete(scoring, a, len_a, b, len_b, NULL))
    return ALN_ERR_UNSUPPORTED;
  sub_table_t table;
  if(!sub_table_build(&table, scoring, a, len_a, b, len_b))
    return ALN_ERR_NOMEM;
  const size_t budget = opts->wfa_budget != 0
                        ? opts->wfa_budget
                        : alignment_wfa_default_budget(len_a, len_b);
  aln_status_t status = alignment_wfa_align(&table, scoring, a, len_a, b,
                                            len_b, budget, result);
  sub_table_free(&table);
  return status == ALN_OK || status == ALN_ERR_NOMEM ? status
                                                     : ALN_ERR_UNSUPPORTED;
}

// A long pair (the longer sequence at least alignment_linear_above
// letters): align in linear space, if the options ask and it supports the
// scoring. ALN_OK if that worked, ALN_ERR_NOMEM if it ran out of memory,
// otherwise ALN_ERR_UNSUPPORTED: the matrices align it
static aln_status_t _nw_long(const char *a, const char *b, size_t len_a,
                             size_t len_b, const scoring_t *scoring,
                             const nw_aligner_t *nw, alignment_t *result)
{
  alignment_options_t defaults;
  const alignment_options_t *opts = _nw_opts(nw, &defaults);
  if(opts->linear_above == 0 || MAX2(len_a, len_b) < opts->linear_above)
    return ALN_ERR_UNSUPPORTED;
  aln_status_t status = hirschberg_align_checked(a, b, len_a, len_b, scoring,
                                                 opts, result, NULL);
  return status == ALN_OK || status == ALN_ERR_NOMEM ? status
                                                     : ALN_ERR_UNSUPPORTED;
}

// The fast paths that go before the matrices: ALN_OK or ALN_ERR_NOMEM
// settle it, ALN_ERR_UNSUPPORTED leaves it to the matrices. checked: the
// input has passed alignment_check()
static aln_status_t _nw_fast_paths(const char *a, const char *b,
                                   size_t len_a, size_t len_b,
                                   const scoring_t *scoring,
                                   const nw_aligner_t *nw,
                                   alignment_t *result, bool checked)
{
  aln_status_t status = _nw_wfa(a, b, len_a, len_b, scoring, nw, result,
                                checked);
  if(status == ALN_ERR_UNSUPPORTED)
    status = _nw_long(a, b, len_a, len_b, scoring, nw, result);
  return status;
}

// Gaps before the first aligned pair, then shift the alignment (written
// backwards from next_char) to the start of the buffers. (With two empty
// sequences longest_alignment is 0 and next_char has wrapped to SIZE_MAX:
// first_char = next_char + 1 wraps back to 0, and the alignment is empty)
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

// The traceback from the traceback bytes of a compact fill. Returns 1, or
// 0 if a step has no predecessor (the caller then fills the matrices, to
// report it), or -1 if out of memory (filling a band's strip again, see
// alignment_trace_load())
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
  if(status != ALN_OK) return status;
  status = _nw_fast_paths(a, b, len_a, len_b, scoring, nw, result, true);
  if(status == ALN_ERR_UNSUPPORTED) {
    status = _nw_align(a, b, len_a, len_b, scoring, nw, result);
    if(_nw_linear(a, b, len_a, len_b, scoring, nw, result, status))
      status = ALN_OK;
  }
  if(err != NULL) err->status = status;
  return status;
}

void needleman_wunsch_align2(const char *a, const char *b,
                             size_t len_a, size_t len_b,
                             const scoring_t *scoring,
                             nw_aligner_t *nw, alignment_t *result)
{
  // As needleman_wunsch_align_checked() without alignment_check(): a pair
  // the scoring lacks ends the program from the table build instead
  aln_error_t e = {ALN_OK, 0, 0};
  e.status = _nw_fast_paths(a, b, len_a, len_b, scoring, nw, result, false);
  if(e.status == ALN_ERR_UNSUPPORTED) {
    e.status = _nw_align(a, b, len_a, len_b, scoring, nw, result);
    if(_nw_linear(a, b, len_a, len_b, scoring, nw, result, e.status))
      e.status = ALN_OK;
  }
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
  #ifdef NUCLEICS_ALIGN_VERBOSE
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

  #ifdef NUCLEICS_ALIGN_VERBOSE
    alignment_print_matrices(nw);
  #endif

  result->score = curr_score;
  char *alignment_a = result->result_a, *alignment_b = result->result_b;

  // coords in score matrices
  size_t score_x = nw->score_width-1, score_y = nw->score_height-1;
  size_t arr_index = arr_size - 1;

  for(; score_x > 0 && score_y > 0; next_char--)
  {
    #ifdef NUCLEICS_ALIGN_VERBOSE
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

    alignment_reverse_move(&curr_matrix, &curr_score,
                           &score_x, &score_y, &arr_index, nw);
  }

  _nw_finish(nw, result, score_x, score_y, next_char, longest_alignment);
  return ALN_OK;
}
