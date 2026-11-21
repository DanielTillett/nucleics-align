/*
 tools/sw_cmdline.c
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

// my utility functions
#include "seq_file/seq_file.h"

// Alignment scoring and loading
#include "alignment_scoring_load.h"
#include "alignment_cmdline.h"
#include "alignment_macros.h"

#include "smith_waterman.h"
#include "hirschberg.h"

static cmdline_t *cmd;
static scoring_t scoring;

static size_t alignment_index = 0;
static bool wait_on_keystroke = 0;

static void sw_set_default_scoring(void)
{
  scoring_system_default(&scoring);

  // Change slightly
  scoring.match = 2;
  scoring.mismatch = -2;
  scoring.gap_open = -2;
  scoring.gap_extend = -1;
}

static void* sw_aligner_new(const cmdline_t *c)
{
  sw_aligner_t *sw = smith_waterman_new();
  // Just the best hit: traceback bytes instead of score matrices
  smith_waterman_set_compact(sw, c->max_hits_per_alignment_set &&
                                 c->max_hits_per_alignment == 1 &&
                                 !c->print_matrices);
  return sw;
}

static void sw_aligner_free(void *sw)
{
  smith_waterman_free(sw);
}

static void sw_set_options(void *sw, const alignment_options_t *opts)
{
  smith_waterman_set_options(sw, opts);
}

// Interactive input: wait for a key after each hit
static void sw_on_stdin(void)
{
  wait_on_keystroke = 1;
}

// Print one line of an alignment
static void print_alignment_part(FILE *out, const char* seq1, const char* seq2,
                                 size_t pos, size_t len,
                                 const char* context_str,
                                 size_t spaces_left, size_t spaces_right,
                                 size_t context_left, size_t context_right)
{
  size_t i;
  fprintf(out, "  ");

  for(i = 0; i < spaces_left; i++) fprintf(out, " ");

  if(context_left > 0)
  {
    if(cmd->print_colour) fputs(align_col_context, out);
    fprintf(out, "%.*s", (int)context_left, context_str+pos-context_left);
    if(cmd->print_colour) fputs(align_col_stop, out);
  }

  if(cmd->print_colour)
    alignment_colour_fprint_against(out, seq1, seq2, scoring.case_sensitive);
  else
    fputs(seq1, out);

  if(context_right > 0)
  {
    if(cmd->print_colour) fputs(align_col_context, out);
    fprintf(out, "%.*s", (int)context_right, context_str+pos+len);
    if(cmd->print_colour) fputs(align_col_stop, out);
  }

  for(i = 0; i < spaces_right; i++) putc(' ', out);

  fprintf(out, "  [pos: %zu; len: %zu]\n", pos, len);
}

static char get_next_hit(void)
{
  if(!wait_on_keystroke)
    return 1;

  int r = 0;

  char response = 0;
  char next_hit = 0;

  while(!response)
  {
    printf("next [h]it or [a]lignment: ");
    fflush(stdout);

    while((r = getc(stdin)) != -1 && r != '\n' && r != '\r')
    {
      if(r == 'h' || r == 'H')
      {
        next_hit = 1;
        response = 1;
      }
      else if(r == 'a' || r == 'A')
      {
        next_hit = 0;
        response = 1;
      }
    }

    if(r == -1)
    {
      // We're done here
      putc('\n', stdout);
      exit(EXIT_SUCCESS);
    }
  }

  return next_hit;
}

// Align two sequences against each other to find local alignments between
// them, and print them to out (errors to err) as alignment number index
static int align(tool_worker_t *wk, const char *seq_a, const char *seq_b,
                  const char *seq_a_name, const char *seq_b_name,
                  size_t index, FILE *out, FILE *err)
{
  sw_aligner_t *sw = wk->aligner;
  alignment_t *result = wk->result;

  if((seq_a_name != NULL || seq_b_name != NULL) && wait_on_keystroke)
  {
    fprintf(err, "Error: Interactive input takes seq only "
                 "(no FASTA/FASTQ) '%s:%s'\n", seq_a_name, seq_b_name);
    fflush(err);
    exit(EXIT_FAILURE);
  }

  // Check both arguments have length > 0
  if(seq_a[0] == '\0' || seq_b[0] == '\0')
  {
    fprintf(err, "Error: Sequences must have length > 0\n");
    fflush(err);

    if(cmd->print_fasta && seq_a_name != NULL && seq_b_name != NULL)
    {
      fprintf(err, "%s\n%s\n", seq_a_name, seq_b_name);
    }

    fflush(err);

    return 0;
  }

  size_t len_a = strlen(seq_a), len_b = strlen(seq_b);

  // If min_score hasn't been set, set a limit based on the lengths of seqs
  // or zero if we're running interactively
  score_t min_score = cmd->min_score_set ? cmd->min_score
                      : wait_on_keystroke ? 0
                      : scoring.match * MAX2(0.2 * MIN2(len_a, len_b), 2);

  #ifdef NUCLEICS_ALIGN_VERBOSE
  fprintf(out, "min_score: %i\n", min_score);
  #endif

  // Hits below min_score are never printed, so do not index them
  smith_waterman_set_min_score(sw, min_score);
  aln_error_t e;
  if(smith_waterman_align_checked(seq_a, seq_b, len_a, len_b, &scoring, sw,
                                  &e) != ALN_OK) {
    alignment_fprint_error(err, &e);
    fflush(err);
    return EXIT_FAILURE;
  }

  aligner_t *aligner = smith_waterman_get_aligner(sw);

  fprintf(out, "== Alignment %zu lengths (%zu, %zu):\n", index, len_a, len_b);

  if(cmd->print_matrices)
  {
    // Only one pair at a time, so out is stdout
    alignment_print_matrices(aligner);
  }

  // seqA
  if(cmd->print_fasta && seq_a_name != NULL)
  {
    fputs(seq_a_name, out);
    putc('\n', out);
  }

  if(cmd->print_seq)
  {
    fputs(seq_a, out);
    putc('\n', out);
  }

  // seqB
  if(cmd->print_fasta && seq_b_name != NULL)
  {
    fputs(seq_b_name, out);
    putc('\n', out);
  }

  if(cmd->print_seq)
  {
    fputs(seq_b, out);
    putc('\n', out);
  }

  putc('\n', out);

  fflush(out);

  size_t hit_index = 0;

  // For print context
  size_t context_left = 0, context_right = 0;
  size_t left_spaces_a = 0, left_spaces_b = 0;
  size_t right_spaces_a = 0, right_spaces_b = 0;


  // Check the hit limit before fetching: finding a hit that would not be
  // printed can mean following many rejected paths
  while(get_next_hit() &&
        (!cmd->max_hits_per_alignment_set ||
         hit_index < cmd->max_hits_per_alignment) &&
        smith_waterman_fetch_checked(sw, result, &e) &&
        result->score >= min_score)
  {
    fprintf(out, "hit %zu.%zu score: %i\n", index, hit_index++, result->score);

    if(cmd->print_context)
    {
      // Calculate number of characters of context to print either side
      context_left = MAX2(result->pos_a, result->pos_b);
      context_left = MIN2(context_left, cmd->print_context);

      size_t rem_a = len_a - (result->pos_a + result->len_a);
      size_t rem_b = len_b - (result->pos_b + result->len_b);

      context_right = MAX2(rem_a, rem_b);
      context_right = MIN2(context_right, cmd->print_context);

      left_spaces_a = (context_left > result->pos_a)
                      ? context_left - result->pos_a : 0;

      left_spaces_b = (context_left > result->pos_b)
                      ? context_left - result->pos_b : 0;

      right_spaces_a = (context_right > rem_a) ? context_right - rem_a : 0;
      right_spaces_b = (context_right > rem_b) ? context_right - rem_b : 0;
    }

    #ifdef NUCLEICS_ALIGN_VERBOSE
    fprintf(out, "context left = %zu; right = %zu spacing: [%zu,%zu] [%zu,%zu]\n",
            context_left, context_right,
            left_spaces_a, right_spaces_a,
            left_spaces_b, right_spaces_b);
    #endif

    // seq a
    print_alignment_part(out, result->result_a, result->result_b,
                         result->pos_a, result->len_a,
                         seq_a,
                         left_spaces_a, right_spaces_a,
                         context_left-left_spaces_a,
                         context_right-right_spaces_a);

    if(cmd->print_pretty)
    {
      fputs("  ", out);

      size_t max_left_spaces = MAX2(left_spaces_a, left_spaces_b);
      size_t max_right_spaces = MAX2(right_spaces_a, right_spaces_b);
      size_t spacer;

      // Print spaces for lefthand spacing
      for(spacer = 0; spacer < max_left_spaces; spacer++)
      {
        putc(' ', out);
      }

      // Print dots for lefthand context sequence
      for(spacer = 0; spacer < context_left-max_left_spaces; spacer++)
      {
        putc('.', out);
      }

      alignment_fprint_spacer(out, result->result_a, result->result_b,
                              &scoring);

      // Print dots for righthand context sequence
      for(spacer = 0; spacer < context_right-max_right_spaces; spacer++)
      {
        putc('.', out);
      }

      // Print spaces for righthand spacing
      for(spacer = 0; spacer < max_right_spaces; spacer++)
      {
        putc(' ', out);
      }

      putc('\n', out);
    }

    // seq b
    print_alignment_part(out, result->result_b, result->result_a,
                         result->pos_b, result->len_b,
                         seq_b,
                         left_spaces_b, right_spaces_b,
                         context_left-left_spaces_b,
                         context_right-right_spaces_b);

    fprintf(out, "\n");

    // Flush output here
    fflush(out);
  }

  // Out of memory fetching a hit: the hits before it are printed
  if(e.status != ALN_OK) {
    alignment_fprint_error(err, &e);
    fflush(err);
    return EXIT_FAILURE;
  }

  fputs("==\n", out);
  fflush(out);
  return 0;
}

// The number of the next alignment: pairs that align() rejects as empty
// (it reads them as strings) get none
static size_t number_pair(const char *seq_a, const char *seq_b)
{
  return seq_a[0] == '\0' || seq_b[0] == '\0' ? 0 : alignment_index++;
}

int main(int argc, char* argv[])
{
  #ifdef NUCLEICS_ALIGN_VERBOSE
  printf("VERBOSE: on\n");
  #endif

  sw_set_default_scoring();
  cmd = cmdline_new(argc, argv, &scoring, NUCLEICS_ALIGN_SW_CMD);

  static const cmd_tool_t tool = {sw_aligner_new, sw_aligner_free,
                                  sw_set_options, align, number_pair,
                                  sw_on_stdin};
  return cmdline_run_tool(cmd, &tool);
}
