/*
 alignment_scoring_load.c
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013

 Scoring files for the command line tools: the (perhaps gzipped) file is
 read whole, parsed by the library (alignment_scoring_parse.c), and an
 error ends the program naming the file and line
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "alignment_scoring_load.h"

static void out_of_memory(void) __attribute__((noreturn));
static void out_of_memory(void)
{
  fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);
}

// The whole of file into memory (len bytes), or exit on a read error (a
// damaged gzip file, say)
static char *read_whole_file(gzFile file, size_t *len)
{
  size_t cap = 1 << 16, n = 0;
  char *buf = malloc(cap);
  if(buf == NULL) out_of_memory();
  for(;;) {
    if(n == cap) {
      if(cap > SIZE_MAX / 2) out_of_memory();
      cap *= 2;
      char *nb = realloc(buf, cap);
      if(nb == NULL) out_of_memory();
      buf = nb;
    }
    const size_t room = cap - n;
    int got = gzread(file, buf + n, room > INT_MAX ? INT_MAX : (unsigned)room);
    if(got <= 0) {
      // The end of the file, or a damaged one (which may read as the end
      // of the data rather than an error): not the same thing
      int errnum = 0;
      const char *msg = gzerror(file, &errnum);
      if(got < 0 || (errnum != Z_OK && errnum != Z_STREAM_END)) {
        fprintf(stderr, "Error: reading scoring file: %s\n",
                msg != NULL && *msg != '\0' ? msg : "read error");
        exit(EXIT_FAILURE);
      }
      break;
    }
    n += (size_t)got;
  }
  *len = n;
  return buf;
}

static void loading_error(scoring_parse_t status, const char *file_path,
                          int line_num, bool is_matrix) __attribute__((noreturn));
static void loading_error(scoring_parse_t status, const char *file_path,
                          int line_num, bool is_matrix)
{
  if(status == SCORING_PARSE_NOMEM) out_of_memory();
  fprintf(stderr, "Error: substitution %s : %s\n",
          is_matrix ? "matrix" : "pairs", scoring_parse_strerror(status));
  if(file_path != NULL) fprintf(stderr, "File: %s\n", file_path);
  if(line_num != -1) fprintf(stderr, "Line: %i\n", line_num);
  exit(EXIT_FAILURE);
}

// Line numbers reported in errors are 1-based
void align_scoring_load_matrix(gzFile file, const char* file_path,
                               scoring_t* scoring)
{
  size_t len;
  char *text = read_whole_file(file, &len);
  int line;
  scoring_parse_t status = scoring_parse_matrix(scoring, text, len, &line);
  free(text);
  if(status != SCORING_PARSE_OK) loading_error(status, file_path, line, true);
}

void align_scoring_load_pairwise(gzFile file, const char* file_path,
                                 scoring_t* scoring)
{
  size_t len;
  char *text = read_whole_file(file, &len);
  int line;
  scoring_parse_t status = scoring_parse_pairs(scoring, text, len, &line);
  free(text);
  if(status != SCORING_PARSE_OK) loading_error(status, file_path, line, false);
}
