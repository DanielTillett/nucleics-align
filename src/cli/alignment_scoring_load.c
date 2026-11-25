/*
 alignment_scoring_load.c
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
#include <ctype.h> // tolower isspace

#include "alignment_cmdline.h"
#include "alignment_scoring_load.h"

#define lc(c,case_sensitive) ((case_sensitive) ? (c) : (char)tolower((unsigned char)(c)))

// Minimal growable line buffer for reading gzFile lines without external deps
typedef struct {
  char *b;
  size_t end;
  size_t cap;
} line_buf_t;

static void line_buf_init(line_buf_t *lb, size_t init_cap)
{
  lb->cap = init_cap > 0 ? init_cap : 256;
  lb->b = malloc(lb->cap);
  if(lb->b == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  lb->b[0] = '\0';
  lb->end = 0;
}

static void line_buf_free(line_buf_t *lb)
{
  free(lb->b);
  lb->b = NULL;
  lb->end = lb->cap = 0;
}

static size_t line_buf_gzreadline(line_buf_t *lb, gzFile file)
{
  lb->end = 0;
  lb->b[0] = '\0';
  for(;;) {
    if(lb->cap - lb->end < 128) {
      size_t ncap = lb->cap * 2;
      char *nb = realloc(lb->b, ncap);
      if(nb == NULL) {
        fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }
      lb->b = nb;
      lb->cap = ncap;
    }
    int avail = (lb->cap - lb->end > (size_t)INT_MAX) ? INT_MAX : (int)(lb->cap - lb->end);
    if(gzgets(file, lb->b + lb->end, avail) == NULL) {
      // The end of the file, or a damaged one: not the same thing
      int errnum = 0;
      const char *msg = gzerror(file, &errnum);
      if(errnum != Z_OK && errnum != Z_STREAM_END) {
        fprintf(stderr, "Error: reading scoring file: %s\n",
                msg != NULL && *msg != '\0' ? msg : "read error");
        exit(EXIT_FAILURE);
      }
      return lb->end;
    }
    size_t got = strlen(lb->b + lb->end);
    lb->end += got;
    if(got == 0 || lb->b[lb->end - 1] == '\n')
      return lb->end;
  }
}

static void line_buf_chomp(line_buf_t *lb)
{
  while(lb->end > 0 && (lb->b[lb->end - 1] == '\n' || lb->b[lb->end - 1] == '\r'))
    lb->b[--lb->end] = '\0';
}

static bool str_is_all_whitespace(const char *s)
{
  while(*s != '\0') {
    if(!isspace((unsigned char)*s)) return false;
    s++;
  }
  return true;
}

static char *str_next_nonwhitespace(char *s)
{
  while(*s != '\0' && isspace((unsigned char)*s)) s++;
  return *s == '\0' ? NULL : s;
}

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
  line_buf_t sbuf_val;
  line_buf_t *sbuf = &sbuf_val;
  line_buf_init(sbuf, 500);
  size_t read_length;
  int line_num = 0;

  // Read first line (column headings)
  while((read_length = line_buf_gzreadline(sbuf, file)) > 0)
  {
    line_num++;
    line_buf_chomp(sbuf);

    if(sbuf->end > 0 && sbuf->b[0] != '#' && // line is not empty, not comment
       !str_is_all_whitespace(sbuf->b)) // and not whitespace
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

    while((next = str_next_nonwhitespace(next+1)) != NULL)
    {
      characters[num_of_chars++] = lc(*next, case_sensitive);
    }

    // Now read lines below
    while((read_length = line_buf_gzreadline(sbuf, file)) > 0)
    {
      line_num++;
      line_buf_chomp(sbuf);

      char* from_char_pos = str_next_nonwhitespace(sbuf->b);

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

        score_txt = str_next_nonwhitespace(score_txt+1);

        if(score_txt == NULL)
        {
          _loading_error("Too few columns on row", file_path, line_num, 1);
        }

        char* strtol_last_char_ptr = score_txt;
        long num = strtol(strtol_last_char_ptr, &strtol_last_char_ptr, 10);

        // If pointer to end of number string hasn't moved -> error
        if(strtol_last_char_ptr == score_txt)
        {
          _loading_error("Missing number value on line", file_path, line_num, 1);
        }
        if(num < -SCORE_ABS_MAX || num > SCORE_ABS_MAX)
          _loading_error("Score out of range", file_path, line_num, 1);
        score = (int)num;

        scoring_add_mutation(scoring, from_char, to_char, score);

        score_txt = strtol_last_char_ptr;
      }

      if(*score_txt != '\0' && !str_is_all_whitespace(score_txt))
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
    while((read_length = line_buf_gzreadline(sbuf, file)) > 0)
    {
      line_num++;
      line_buf_chomp(sbuf);

      if(sbuf->b[0] == '#' || str_is_all_whitespace(sbuf->b))
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
        long num = strtol(str_pos, &after_num_str, 10);

        // If pointer to end of number string hasn't moved -> error
        if(str_pos == after_num_str)
        {
          _loading_error("Missing number value on line", file_path, line_num, 1);
        }
        if(num < -SCORE_ABS_MAX || num > SCORE_ABS_MAX)
          _loading_error("Score out of range", file_path, line_num, 1);
        score = (int)num;

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
  line_buf_free(sbuf);
}


void align_scoring_load_pairwise(gzFile file, const char* file_path,
                                 scoring_t* scoring, char case_sensitive)
{
  line_buf_t sbuf_val;
  line_buf_t *sbuf = &sbuf_val;
  line_buf_init(sbuf, 200);
  size_t read_length;
  int line_num = 0;

  char a, b;
  int score;

  int num_pairs_added = 0;

  while((read_length = line_buf_gzreadline(sbuf, file)) > 0)
  {
    line_num++;
    line_buf_chomp(sbuf);

    if(sbuf->end > 0 && sbuf->b[0] != '#' && // line is not empty, not comment
       !str_is_all_whitespace(sbuf->b)) // and not whitespace
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

      if(score < -SCORE_ABS_MAX || score > SCORE_ABS_MAX)
        _loading_error("Score out of range", file_path, line_num, 0);

      scoring_add_mutation(scoring, lc(a, case_sensitive), lc(b, case_sensitive),
                           score);
      num_pairs_added++;
    }
  }

  line_buf_free(sbuf);

  if(num_pairs_added == 0)
  {
    _loading_error("No pairs added from file (file empty?)",
                   file_path, -1, 0);
  }
}
