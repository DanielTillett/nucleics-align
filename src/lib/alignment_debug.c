/*
 alignment_debug.c
 Matrix debug printing and coloured alignment output helpers.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "alignment.h"
#include "alignment_macros.h"

const char align_col_mismatch[] = "\033[92m"; // Mismatch (GREEN)
const char align_col_indel[] = "\033[91m"; // Insertion / deletion (RED)
// Pink used by SmithWaterman local alignment for printing surrounding bases
const char align_col_context[] = "\033[95m";
const char align_col_stop[] = "\033[0m";


void alignment_fprint_matrices(FILE *fp, const aligner_t *aligner)
{
  if(fp == NULL || aligner == NULL) return;
  const score_t* match_scores = aligner->match_scores;
  const score_t* gap_a_scores = aligner->gap_a_scores;
  const score_t* gap_b_scores = aligner->gap_b_scores;

  size_t i, j;

  fprintf(fp, "seq_a: %.*s\nseq_b: %.*s\n",
          (int)aligner->score_width-1, aligner->seq_a,
          (int)aligner->score_height-1, aligner->seq_b);

  fprintf(fp, "match_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    fprintf(fp, "%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      fprintf(fp, "\t%3i", (int)ARR_LOOKUP(match_scores, aligner->score_width, i, j));
    }
    putc('\n', fp);
  }
  fprintf(fp, "gap_a_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    fprintf(fp, "%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      fprintf(fp, "\t%3i", (int)ARR_LOOKUP(gap_a_scores, aligner->score_width, i, j));
    }
    putc('\n', fp);
  }
  fprintf(fp, "gap_b_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    fprintf(fp, "%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      fprintf(fp, "\t%3i", (int)ARR_LOOKUP(gap_b_scores, aligner->score_width, i, j));
    }
    putc('\n', fp);
  }

  fprintf(fp, "match: %i mismatch: %i gapopen: %i gapexend: %i\n",
          aligner->scoring->match, aligner->scoring->mismatch,
          aligner->scoring->gap_open, aligner->scoring->gap_extend);
  fprintf(fp, "\n");
}

void alignment_print_matrices(const aligner_t *aligner)
{
  alignment_fprint_matrices(stdout, aligner);
}

void alignment_colour_fprint_against(FILE *out, const char *alignment_a,
                                     const char *alignment_b,
                                     char case_sensitive)
{
  size_t i;
  char red = 0, green = 0;

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_b[i] == '-')
    {
      if(!red)
      {
        fputs(align_col_indel, out);
        red = 1;
      }
    }
    else if(red)
    {
      red = 0;
      fputs(align_col_stop, out);
    }

    if(((case_sensitive && alignment_a[i] != alignment_b[i]) ||
        (!case_sensitive && tolower((unsigned char)alignment_a[i]) !=
                            tolower((unsigned char)alignment_b[i]))) &&
       alignment_a[i] != '-' && alignment_b[i] != '-')
    {
      if(!green)
      {
        fputs(align_col_mismatch, out);
        green = 1;
      }
    }
    else if(green)
    {
      green = 0;
      fputs(align_col_stop, out);
    }

    putc(alignment_a[i], out);
  }

  if(green || red)
  {
    // Stop all colours
    fputs(align_col_stop, out);
  }
}

// Order of alignment_a / alignment_b is not important
void alignment_fprint_spacer(FILE *out, const char* alignment_a,
                             const char* alignment_b,
                             const scoring_t* scoring)
{
  size_t i;

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '-' || alignment_b[i] == '-')
    {
      putc(' ', out);
    }
    else if(alignment_a[i] == alignment_b[i] ||
            (!scoring->case_sensitive &&
             tolower((unsigned char)alignment_a[i]) ==
             tolower((unsigned char)alignment_b[i])))
    {
      putc('|', out);
    }
    else
    {
      putc('*', out);
    }
  }
}

void alignment_colour_print_against(const char *alignment_a,
                                    const char *alignment_b,
                                    char case_sensitive)
{
  alignment_colour_fprint_against(stdout, alignment_a, alignment_b,
                                  case_sensitive);
}

void alignment_print_spacer(const char* alignment_a, const char* alignment_b,
                            const scoring_t* scoring)
{
  alignment_fprint_spacer(stdout, alignment_a, alignment_b, scoring);
}
