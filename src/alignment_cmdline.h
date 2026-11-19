/*
 alignment_cmdline.h
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */

#ifndef ALIGNMENT_CMDLINE_HEADER_SEEN
#define ALIGNMENT_CMDLINE_HEADER_SEEN

// request decent POSIX version
#ifndef _XOPEN_SOURCE
  #define _XOPEN_SOURCE 700
#endif
#ifndef _DEFAULT_SOURCE
  #define _DEFAULT_SOURCE
#endif

#include <stdarg.h> // required for va_list
#include <stdbool.h>
#include "seq_file/seq_file.h"
#include "alignment.h"

enum SeqAlignCmdType {SEQ_ALIGN_SW_CMD, SEQ_ALIGN_NW_CMD, SEQ_ALIGN_LCS_CMD};

typedef struct
{
  // file inputs
  size_t file_list_length, file_list_capacity;
  char **file_paths1, **file_paths2;

  // All values initially 0
  bool case_sensitive;
  int match, mismatch, gap_open, gap_extend;

  // SW specific
  score_t min_score;
  unsigned int print_context, max_hits_per_alignment;
  bool min_score_set, max_hits_per_alignment_set;
  bool print_seq;

  // NW specific
  bool freestartgap_set, freeendgap_set;
  bool print_matrices, print_scores;
  bool zam_stle_output;
  bool linear_space; // Hirschberg / Myers-Miller
  bool score_only; // print the score without the alignment

  // Turns off zlib for stdin
  bool interactive;

  // --threads (0: not set, half the CPUs)
  unsigned int threads;

  // General output
  bool print_fasta, print_pretty, print_colour;

  // Experimental
  bool no_gaps_in1, no_gaps_in2;
  bool no_mismatches;

  // Pair of sequences to align
  const char *seq1, *seq2;
} cmdline_t;

char parse_entire_int(char *str, int *result);
char parse_entire_uint(char *str, unsigned int *result);

cmdline_t* cmdline_new(int argc, char **argv, scoring_t *scoring,
                       enum SeqAlignCmdType cmd_type);
void cmdline_free(cmdline_t* cmd);

void cmdline_add_files(cmdline_t* cmd, char* p1, char* p2);
size_t cmdline_get_num_of_file_pairs(cmdline_t* cmd);
char* cmdline_get_file1(cmdline_t* cmd, size_t i);
char* cmdline_get_file2(cmdline_t* cmd, size_t i);

void align_from_file(const char *path1, const char *path2,
                     void (align)(read_t *r1, read_t *r2),
                     bool use_zlib);

// Threads for the command line tools: --threads, or half the CPUs
unsigned int cmdline_threads(const cmdline_t *cmd);

// Like align_from_file(), but aligns several pairs at once on up to
// `threads` threads. For each pair, in input order on this thread, number()
// gives a number for it (may be NULL: the pair's position); then run()
// aligns it on a worker thread with that worker's state (from new_worker())
// and the settings to use (opts, with its thread counts set: one beside
// other pairs, `threads` alone; opts may be NULL for the defaults),
// writing to out and err instead of stdout and stderr, and returns 0 or a
// status to exit with. No global setting is changed. Each pair's output is
// written as soon as the pairs before it are, so it is the same as one at a
// time. A pair of more than opts->block_min_cells cells runs on its own
// with all the threads. run() should return errors in the input (see the
// *_checked() aligners) rather than exit, so the pairs before are written
void align_from_file_threaded(const char *path1, const char *path2,
                              bool use_zlib, unsigned int threads,
                              size_t (*number)(read_t *r1, read_t *r2),
                              const alignment_options_t *opts,
                              int (*run)(void *worker, read_t *r1, read_t *r2,
                                         size_t number,
                                         const alignment_options_t *opts,
                                         FILE *out, FILE *err),
                              void* (*new_worker)(void),
                              void (*free_worker)(void *worker));

#endif
