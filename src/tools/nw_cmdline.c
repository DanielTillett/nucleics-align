/*
 tools/nw_cmdline.h
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _DEFAULT_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower

// Alignment scoring and loading
#include "alignment_cmdline.h"

#include "needleman_wunsch.h"
#include "hirschberg.h"

// For this run
static cmdline_t *cmd;
static scoring_t scoring;

static void nw_set_default_scoring(void)
{
  scoring_system_default(&scoring);
}

static void* nw_aligner_new(const cmdline_t *c)
{
  nw_aligner_t *nw = needleman_wunsch_new();
  // Traceback bytes instead of score matrices, unless they are printed
  nw->compact = !c->print_matrices;
  return nw;
}

static void nw_aligner_free(void *nw)
{
  needleman_wunsch_free(nw);
}

static void nw_set_options(void *nw, const alignment_options_t *opts)
{
  ((nw_aligner_t*)nw)->options = opts;
}

// Write an error as the library would before exiting, or --linear's own
// message. Returns 0 for none, else EXIT_FAILURE
static int nw_error(const aln_error_t *e, FILE *err)
{
  if(e->status == ALN_OK) return 0;
  if(e->status == ALN_ERR_UNSUPPORTED)
    fprintf(err, "Error: --linear does not support these scoring options "
                 "(or penalties this large)\n");
  else alignment_fprint_error(err, e);
  fflush(err);
  return EXIT_FAILURE;
}

// Returns 0, or EXIT_FAILURE with the error written to err
static int nw_align(tool_worker_t *wk, const char *seq_a, const char *seq_b,
                    FILE *err)
{
  nw_aligner_t *nw = wk->aligner;
  const size_t len_a = strlen(seq_a), len_b = strlen(seq_b);
  aln_error_t e;
  if(!cmd->linear_space)
    needleman_wunsch_align_checked(seq_a, seq_b, len_a, len_b, &scoring,
                                   nw, wk->result, &e);
  else
    hirschberg_align_checked(seq_a, seq_b, len_a, len_b, &scoring,
                             nw->options, wk->result, &e);
  return nw_error(&e, err);
}

static int align_zam(tool_worker_t *wk, const char *seq_a, const char *seq_b,
                     FILE *out, FILE *err)
{
  alignment_t *result = wk->result;
  int status = nw_align(wk, seq_a, seq_b, err);
  if(status != 0) return status;

  // Swap '-' for '_'
  int i;
  for(i = 0; result->result_a[i] != '\0'; i++)
  {
    if(result->result_a[i] == '-') result->result_a[i] = '_';
    if(result->result_b[i] == '-') result->result_b[i] = '_';
  }

  int num_of_mismatches = 0;
  int num_of_indels = 0;

  // Print branch 1 and spacer
  fprintf(out, "Br1:%s\n    ", result->result_a);

  for(i = 0; result->result_a[i] != '\0'; i++)
  {
    if(result->result_a[i] == '_' || result->result_b[i] == '_')
    {
      putc(' ', out);
      num_of_indels++;
    }
    else if((scoring.case_sensitive && result->result_a[i] != result->result_b[i]) ||
            tolower((unsigned char)result->result_a[i]) !=
            tolower((unsigned char)result->result_b[i]))
    {
      putc('*', out);
      num_of_mismatches++;
    }
    else
    {
      putc('|', out);
    }
  }

  // Print branch 2 and mismatch indel numbers
  fprintf(out, "\nBr2:%s\n%i %i\n\n", result->result_b,
          num_of_mismatches, num_of_indels);
  return 0;
}

// Align a pair and print it to out (errors to err). Returns 0, or a status
// to exit with. Pairs are not numbered
static int align(tool_worker_t *wk, const char *seq_a, const char *seq_b,
                 const char *seq_a_name, const char *seq_b_name,
                 size_t number, FILE *out, FILE *err)
{
  nw_aligner_t *nw = wk->aligner;
  alignment_t *result = wk->result;
  (void)number;
  int status;

  if(cmd->score_only)
  {
    // One linear-memory pass; options it cannot handle use the full matrices
    score_t score;
    const size_t len_a = strlen(seq_a), len_b = strlen(seq_b);
    aln_error_t e;
    if(hirschberg_score_checked(seq_a, seq_b, len_a, len_b, &scoring,
                                nw->options, &score, &e)
       == ALN_ERR_UNSUPPORTED) {
      needleman_wunsch_align_checked(seq_a, seq_b, len_a, len_b, &scoring,
                                     nw, result, &e);
      score = result->score;
    }
    if(e.status != ALN_OK) return nw_error(&e, err);
    if(cmd->print_fasta && seq_a_name != NULL) fprintf(out, "%s\n", seq_a_name);
    if(cmd->print_fasta && seq_b_name != NULL) fprintf(out, "%s\n", seq_b_name);
    fprintf(out, "score: %i\n\n", score);
    fflush(out);
    return 0;
  }

  if(cmd->zam_stle_output)
  {
    status = align_zam(wk, seq_a, seq_b, out, err);
    fflush(out);
    return status;
  }

  if((status = nw_align(wk, seq_a, seq_b, err)) != 0) return status;

  if(cmd->print_matrices)
  {
    // Only one pair at a time, so out is stdout
    alignment_print_matrices(nw);
  }

  if(cmd->print_fasta && seq_a_name != NULL)
  {
    fputs(seq_a_name, out);
    putc('\n', out);
  }

  if(cmd->print_fasta && cmd->print_pretty && seq_b_name != NULL)
  {
    fputs(seq_b_name, out);
    putc('\n', out);
  }

  if(cmd->print_colour)
  {
    // Print alignment line 1
    alignment_colour_fprint_against(out, result->result_a, result->result_b,
                                    scoring.case_sensitive);
  }
  else
  {
    fputs(result->result_a, out);
  }
  putc('\n', out);

  if(cmd->print_pretty)
  {
    // Print spacer
    alignment_fprint_spacer(out, result->result_a, result->result_b, &scoring);
    putc('\n', out);
  }
  else if(cmd->print_fasta && seq_b_name != NULL)
  {
    fputs(seq_b_name, out);
    putc('\n', out);
  }

  if(cmd->print_colour)
  {
    // Print alignment line 2
    alignment_colour_fprint_against(out, result->result_b, result->result_a,
                                    scoring.case_sensitive);
  }
  else
  {
    fputs(result->result_b, out);
  }
  putc('\n', out);

  if(cmd->print_scores) {
    fprintf(out, "score: %i\n", result->score);
  }

  putc('\n', out);
  fflush(out);
  return 0;
}

int main(int argc, char* argv[])
{
  #ifdef SEQ_ALIGN_VERBOSE
  printf("VERBOSE: on\n");
  #endif

  nw_set_default_scoring();
  cmd = cmdline_new(argc, argv, &scoring, SEQ_ALIGN_NW_CMD);

  static const cmd_tool_t tool = {nw_aligner_new, nw_aligner_free,
                                  nw_set_options, align, NULL, NULL};
  return cmdline_run_tool(cmd, &tool);
}
