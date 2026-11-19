/*
 tools/nw_cmdline.h
 url: https://github.com/DanielTillett/seq-align
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
cmdline_t *cmd;
scoring_t scoring;

// An aligner and somewhere to put its result: one for each thread
typedef struct
{
  nw_aligner_t *nw;
  alignment_t *result;
} nw_worker_t;

static void nw_set_default_scoring(void)
{
  scoring_system_default(&scoring);
}

static void* nw_worker_new(void)
{
  nw_worker_t *wk = malloc(sizeof(nw_worker_t));
  if(wk == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  wk->nw = needleman_wunsch_new();
  // Traceback bytes instead of score matrices, unless they are printed
  wk->nw->compact = !cmd->print_matrices;
  wk->result = alignment_create(256);
  return wk;
}

static void nw_worker_free(void *ptr)
{
  nw_worker_t *wk = ptr;
  needleman_wunsch_free(wk->nw);
  alignment_free(wk->result);
  free(wk);
}

// Returns 0, or EXIT_FAILURE if --linear cannot align with these options
static int nw_align(nw_worker_t *wk, const char *seq_a, const char *seq_b,
                    FILE *err)
{
  if(!cmd->linear_space)
    needleman_wunsch_align(seq_a, seq_b, &scoring, wk->nw, wk->result);
  else if(!hirschberg_align(seq_a, seq_b, &scoring, wk->result))
  {
    fprintf(err, "Error: --linear does not support these scoring options "
                 "(or penalties this large)\n");
    return EXIT_FAILURE;
  }
  return 0;
}

static int align_zam(nw_worker_t *wk, const char *seq_a, const char *seq_b,
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
// to exit with
static int align(nw_worker_t *wk, const char *seq_a, const char *seq_b,
                 const char *seq_a_name, const char *seq_b_name,
                 FILE *out, FILE *err)
{
  alignment_t *result = wk->result;
  int status;

  if(cmd->score_only)
  {
    // One linear-memory pass; options it cannot handle use the full matrices
    score_t score;
    if(!hirschberg_score(seq_a, seq_b, &scoring, &score)) {
      needleman_wunsch_align(seq_a, seq_b, &scoring, wk->nw, result);
      score = result->score;
    }
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
    alignment_print_matrices(wk->nw);
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

static int align_pair_run(void *wk, read_t *read1, read_t *read2,
                          size_t number, FILE *out, FILE *err)
{
  (void)number;
  return align(wk, read1->seq.b, read2->seq.b,
               (read1->name.end == 0 ? NULL : read1->name.b),
               (read2->name.end == 0 ? NULL : read2->name.b), out, err);
}

// One pair at a time, on this thread
static nw_worker_t *main_worker;

static void align_pair_from_file(read_t *read1, read_t *read2)
{
  int status = align_pair_run(main_worker, read1, read2, 0, stdout, stderr);
  if(status != 0) exit(status);
}

int main(int argc, char* argv[])
{
  #ifdef SEQ_ALIGN_VERBOSE
  printf("VERBOSE: on\n");
  #endif

  nw_set_default_scoring();
  cmd = cmdline_new(argc, argv, &scoring, SEQ_ALIGN_NW_CMD);

  // Threads within an alignment, and for pairs from files at once
  unsigned int threads = cmdline_threads(cmd);
  alignment_max_threads = hirschberg_max_threads = (int)threads;

  // Align!
  main_worker = nw_worker_new();

  if(cmd->seq1 != NULL)
  {
    // Align seq1 and seq2 pair passed on the command line
    int status = align(main_worker, cmd->seq1, cmd->seq2, NULL, NULL,
                       stdout, stderr);
    if(status != 0) exit(status);
  }

  // Align from files
  size_t i, num_of_file_pairs = cmdline_get_num_of_file_pairs(cmd);
  for(i = 0; i < num_of_file_pairs; i++)
  {
    const char *file1 = cmdline_get_file1(cmd, i);
    const char *file2 = cmdline_get_file2(cmd, i);
    if(file1 != NULL && *file1 == '\0' && file2 == NULL) {
      file1 = "-";
    }
    // Interactive input, and printed matrices, one pair at a time
    if(threads > 1 && !cmd->interactive && !cmd->print_matrices)
      align_from_file_threaded(file1, file2, !cmd->interactive, threads,
                               &scoring,
                               NULL, align_pair_run, nw_worker_new,
                               nw_worker_free);
    else
      align_from_file(file1, file2, &align_pair_from_file, !cmd->interactive);
  }

  // Free memory for storing alignment results
  nw_worker_free(main_worker);

  cmdline_free(cmd);

  return EXIT_SUCCESS;
}
