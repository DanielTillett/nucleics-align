/*
 alignment_scoring_load.c
 url: https://github.com/DanielTillett/seq-align
 original original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _DEFAULT_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower isspace

#include "string_buffer/string_buffer.h"

#include "alignment_cmdline.h"
#include "alignment_scoring_load.h"

#define lc(c,case_sensitive) ((case_sensitive) ? (c) : (char)tolower((unsigned char)(c)))

static void _loading_error(const char* err_msg, const char* file_path,
                           int line_num, char is_matrix)
  __attribute__((noreturn));

static void _loading_error(const char* err_msg, const char* file_path,
                           int line_num, char is_matrix)
{
  if(is_matrix) fprintf(stderr, "Error: substitution matrix : %s\n", err_msg);
  else fprintf(stderr, "Error: substitution pairs : %s\n", err_msg);

  if(file_path != NULL) fprintf(stderr, "File: %s\n", file_path);
  if(line_num != -1) fprintf(stderr, "Line: %i\n", line_num);

  exit(EXIT_FAILURE);
}

// Line numbers reported in errors are 1-based
void align_scoring_load_matrix(gzFile file, const char* file_path,
                               scoring_t* scoring, char case_sensitive)
{
  StrBuf* sbuf = strbuf_new(500);
  size_t read_length;
  int line_num = 0;

  // Read first line (column headings)
  while((read_length = strbuf_reset_gzreadline(sbuf, file)) > 0)
  {
    line_num++;
    strbuf_chomp(sbuf);

    if(sbuf->end > 0 && sbuf->b[0] != '#' && // line is not empty, not comment
       !string_is_all_whitespace(sbuf->b)) // and not whitespace
    {
      // Read first line

      if(sbuf->end < 2)
      {
        _loading_error("Too few column headings", file_path, line_num, 1);
      }

      break;
    }
  }

  if(read_length == 0)
  {
    // Reached the end of the file without finding the column headings
    _loading_error("Empty file", file_path, -1, 1);
  }

  // If the separator character is whitespace,
  // the set of whitespace characters is used
  char sep = sbuf->b[0];

  if((sep >= '0' && sep <= '9') || sep == '-')
  {
    _loading_error("Numbers (0-9) and dashes (-) do not make good separators",
                   file_path, line_num, 1);
  }

  char* characters = malloc(sbuf->end);
  int num_of_chars = 0;

  if(characters == NULL)
  {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  if(isspace((unsigned char)sep))
  {
    char* next = sbuf->b;

    while((next = string_next_nonwhitespace(next+1)) != NULL)
    {
      characters[num_of_chars++] = lc(*next, case_sensitive);
    }

    // Now read lines below
    while((read_length = strbuf_reset_gzreadline(sbuf, file)) > 0)
    {
      line_num++;
      strbuf_chomp(sbuf);

      char* from_char_pos = string_next_nonwhitespace(sbuf->b);

      if(from_char_pos == NULL || sbuf->b[0] == '#')
      {
        // skip this line
        continue;
      }

      char from_char = lc(*from_char_pos, case_sensitive);
      char to_char;

      char* score_txt = from_char_pos+1;
      int score;

      int i;
      for(i = 0; i < num_of_chars; i++)
      {
        to_char = characters[i];

        if(*score_txt == '\0')
        {
          _loading_error("Too few columns on row", file_path, line_num, 1);
        }

        if(!isspace((unsigned char)*score_txt))
        {
          _loading_error("Expected whitespace between elements - found character",
                         file_path, line_num, 1);
        }

        score_txt = string_next_nonwhitespace(score_txt+1);

        if(score_txt == NULL)
        {
          _loading_error("Too few columns on row", file_path, line_num, 1);
        }

        char* strtol_last_char_ptr = score_txt;
        score = (int)strtol(strtol_last_char_ptr, &strtol_last_char_ptr, 10);

        // If pointer to end of number string hasn't moved -> error
        if(strtol_last_char_ptr == score_txt)
        {
          _loading_error("Missing number value on line", file_path, line_num, 1);
        }

        scoring_add_mutation(scoring, from_char, to_char, score);

        score_txt = strtol_last_char_ptr;
      }

      if(*score_txt != '\0' && !string_is_all_whitespace(score_txt))
      {
        _loading_error("Too many columns on row", file_path, line_num, 1);
      }
    }
  }
  else
  {
    size_t i;

    // Header is <sep><char><sep><char>...
    for(i = 0; i + 1 < sbuf->end; i += 2)
    {
      if(sbuf->b[i] != sep)
      {
        _loading_error("Separator missing from line", file_path, line_num, 1);
      }

      characters[num_of_chars++] = lc(sbuf->b[i+1], case_sensitive);
    }

    if(i < sbuf->end)
    {
      // e.g. ",A,C," - a separator with no heading after it
      _loading_error("Empty column heading", file_path, line_num, 1);
    }

    int score;

    // Read rows: <char><sep><score><sep><score>...
    while((read_length = strbuf_reset_gzreadline(sbuf, file)) > 0)
    {
      line_num++;
      strbuf_chomp(sbuf);

      if(sbuf->b[0] == '#' || string_is_all_whitespace(sbuf->b))
      {
        // skip this line
        continue;
      }

      char from_char = lc(sbuf->b[0], case_sensitive);

      // Skip the row's character
      char* str_pos = sbuf->b+1;

      int to_char_index = 0;
      char to_char;

      while(*str_pos != '\0')
      {
        if(to_char_index >= num_of_chars)
        {
          _loading_error("Too many columns on row", file_path, line_num, 1);
        }

        to_char = characters[to_char_index++];

        if(*str_pos != sep)
        {
          _loading_error("Separator missing from line", file_path, line_num, 1);
        }

        // Move past separator
        str_pos++;

        char* after_num_str = str_pos;
        score = (int)strtol(str_pos, &after_num_str, 10);

        // If pointer to end of number string hasn't moved -> error
        if(str_pos == after_num_str)
        {
          _loading_error("Missing number value on line", file_path, line_num, 1);
        }

        scoring_add_mutation(scoring, from_char, to_char, score);

        str_pos = after_num_str;
      }

      if(to_char_index < num_of_chars)
      {
        _loading_error("Too few columns on row", file_path, line_num, 1);
      }
    }
  }

  free(characters);
  strbuf_free(sbuf);
}


void align_scoring_load_pairwise(gzFile file, const char* file_path,
                                 scoring_t* scoring, char case_sensitive)
{
  StrBuf* sbuf = strbuf_new(200);
  size_t read_length;
  int line_num = 0;

  char a, b;
  int score;

  int num_pairs_added = 0;

  while((read_length = strbuf_reset_gzreadline(sbuf, file)) > 0)
  {
    line_num++;
    strbuf_chomp(sbuf);

    if(sbuf->end > 0 && sbuf->b[0] != '#' && // line is not empty, not comment
       !string_is_all_whitespace(sbuf->b)) // and not whitespace
    {
      if(sbuf->end < 5)
      {
        _loading_error("Line too short", file_path, line_num, 0);
      }

      if(isspace((unsigned char)sbuf->b[1]))
      {
        // split by whitespace
        a = sbuf->b[0];

        size_t char2_pos;

        for(char2_pos = 1;
            sbuf->b[char2_pos] != '\0' && isspace((unsigned char)sbuf->b[char2_pos]);
            char2_pos++);

        if(char2_pos+2 >= sbuf->end ||
           !isspace((unsigned char)sbuf->b[char2_pos+1]))
        {
          _loading_error("Line too short", file_path, line_num, 0);
        }

        b = sbuf->b[char2_pos];

        if(!parse_entire_int(sbuf->b+char2_pos+2, &score))
        {
          _loading_error("Invalid number", file_path, line_num, 0);
        }
      }
      else
      {
        if(sbuf->b[1] != sbuf->b[3])
        {
          _loading_error("Inconsistent separators used", file_path, line_num, 0);
        }

        a = sbuf->b[0];
        b = sbuf->b[2];

        if(!parse_entire_int(sbuf->b + 4, &score))
        {
          _loading_error("Invalid number", file_path, line_num, 0);
        }
      }

      scoring_add_mutation(scoring, lc(a, case_sensitive), lc(b, case_sensitive),
                           score);
      num_pairs_added++;
    }
  }

  strbuf_free(sbuf);

  if(num_pairs_added == 0)
  {
    _loading_error("No pairs added from file (file empty?)",
                   file_path, -1, 0);
  }
}
