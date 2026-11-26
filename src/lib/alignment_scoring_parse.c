/*
 alignment_scoring_parse.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026

 The two scoring file formats in scoring/, read from memory: a substitution
 matrix (a header line of the column characters, separated by whitespace or
 by one separator character that also starts the line, then one row per
 character), and pairs (one "a b score" or "a,b,score" per line). Lines
 starting with # and blank lines are skipped. Both the command line tools
 (src/cli/alignment_scoring_load.c) and the public API (na_api.c) read the
 formats through here, so one grammar serves both
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>

#include "alignment_scoring.h"
#include "alignment_alloc.h"

const char *scoring_parse_strerror(scoring_parse_t status)
{
  switch(status) {
    case SCORING_PARSE_OK:            return "No error";
    case SCORING_PARSE_NOMEM:         return "Out of memory";
    case SCORING_PARSE_EMPTY:         return "Empty file";
    case SCORING_PARSE_FEW_HEADINGS:  return "Too few column headings";
    case SCORING_PARSE_SEPARATOR:
      return "Numbers (0-9) and dashes (-) do not make good separators";
    case SCORING_PARSE_SEP_MISSING:   return "Separator missing from line";
    case SCORING_PARSE_EMPTY_HEADING: return "Empty column heading";
    case SCORING_PARSE_FEW_COLUMNS:   return "Too few columns on row";
    case SCORING_PARSE_EXPECTED_WS:
      return "Expected whitespace between elements - found character";
    case SCORING_PARSE_NUMBER:        return "Missing number value on line";
    case SCORING_PARSE_RANGE:         return "Score out of range";
    case SCORING_PARSE_MANY_COLUMNS:  return "Too many columns on row";
    case SCORING_PARSE_NO_ROWS:       return "No rows in matrix";
    case SCORING_PARSE_SHORT_LINE:    return "Line too short";
    case SCORING_PARSE_SEPARATORS:    return "Inconsistent separators used";
    case SCORING_PARSE_INVALID_NUMBER: return "Invalid number";
    case SCORING_PARSE_NO_PAIRS:
      return "No pairs added from file (file empty?)";
  }
  return "Unknown error";
}

// The next line of text (CR LF or LF ended; a last line needs neither):
// its start and length without the line ending, or NULL at the end
static const char *next_line(const char *p, const char *end,
                             const char **line_start, size_t *line_len)
{
  if(p >= end) return NULL;
  *line_start = p;
  const char *nl = memchr(p, '\n', (size_t)(end - p));
  const char *lend = nl ? nl : end;
  size_t len = (size_t)(lend - p);
  if(len > 0 && p[len - 1] == '\r') len--;
  *line_len = len;
  return nl ? nl + 1 : end;
}

static bool span_all_ws(const char *s, size_t len)
{
  for(size_t i = 0; i < len; i++)
    if(!isspace((unsigned char)s[i])) return false;
  return true;
}

static char fold_ch(char c, bool cs)
{
  return cs ? c : (char)tolower((unsigned char)c);
}

// A score token: whitespace around it is skipped
typedef enum { NUM_OK, NUM_NOT_A_NUMBER, NUM_TOO_BIG } num_t;

static num_t parse_num(const char *s, size_t len, long *out)
{
  if(memchr(s, '\0', len) != NULL) return NUM_NOT_A_NUMBER;
  while(len > 0 && isspace((unsigned char)*s)) { s++; len--; }
  while(len > 0 && isspace((unsigned char)s[len - 1])) len--;
  if(len == 0 || len >= 64) return len == 0 ? NUM_NOT_A_NUMBER : NUM_TOO_BIG;
  char tmp[64];
  memcpy(tmp, s, len);
  tmp[len] = '\0';
  char *ep = tmp;
  errno = 0;
  long v = strtol(tmp, &ep, 10);
  if(ep != tmp + len) return NUM_NOT_A_NUMBER;
  if(errno == ERANGE) return NUM_TOO_BIG;
  *out = v;
  return NUM_OK;
}

// A matrix score: anything that is not a number is NUMBER, a number too
// big for a score (or for a long) is RANGE
static scoring_parse_t matrix_score(const char *s, size_t len, int *score)
{
  long v;
  switch(parse_num(s, len, &v)) {
    case NUM_NOT_A_NUMBER: return SCORING_PARSE_NUMBER;
    case NUM_TOO_BIG: return SCORING_PARSE_RANGE;
    default: break;
  }
  if(v < -SCORE_ABS_MAX || v > SCORE_ABS_MAX) return SCORING_PARSE_RANGE;
  *score = (int)v;
  return SCORING_PARSE_OK;
}

// A pairs score: a number an int cannot hold is not a valid one
// (INVALID_NUMBER, as the tools have always said); one it can hold but a
// score cannot is RANGE
static scoring_parse_t pairs_score(const char *s, size_t len, int *score)
{
  long v;
  if(parse_num(s, len, &v) != NUM_OK || v < INT_MIN || v > INT_MAX)
    return SCORING_PARSE_INVALID_NUMBER;
  if(v < -SCORE_ABS_MAX || v > SCORE_ABS_MAX) return SCORING_PARSE_RANGE;
  *score = (int)v;
  return SCORING_PARSE_OK;
}

// Both parsers read into a copy of the scoring and swap it in only at the
// end, so nothing changes on an error. Every failure sets the status and
// the line and jumps to the cleanup
#define PARSE_FAIL(st, ln) \
  do { status = (st); fail_line = (ln); goto fail; } while(0)

scoring_parse_t scoring_parse_matrix(scoring_t *scoring, const char *text,
                                     size_t len, int *line)
{
  scoring_parse_t status = SCORING_PARSE_OK, st;
  int fail_line = -1;
  scoring_t *tmp = NULL;
  char *chars = NULL;

  if((tmp = aln_malloc(sizeof(*tmp))) == NULL) PARSE_FAIL(SCORING_PARSE_NOMEM, -1);
  *tmp = *scoring;
  const bool cs = tmp->case_sensitive;
  const char *p = text, *end = text + len;
  const char *ls = NULL;
  size_t llen = 0;
  int line_num = 0, rows_added = 0, score;

  // The header line: the first that is not blank or a comment
  bool found = false;
  while((p = next_line(p, end, &ls, &llen)) != NULL) {
    line_num++;
    if(llen > 0 && ls[0] != '#' && !span_all_ws(ls, llen)) {
      if(llen < 2) PARSE_FAIL(SCORING_PARSE_FEW_HEADINGS, line_num);
      found = true;
      break;
    }
  }
  if(!found) PARSE_FAIL(SCORING_PARSE_EMPTY, -1);

  // Its first character is the separator: whitespace (then any run of
  // whitespace separates), or one character
  char sep = ls[0];
  if((sep >= '0' && sep <= '9') || sep == '-')
    PARSE_FAIL(SCORING_PARSE_SEPARATOR, line_num);

  if((chars = aln_malloc(llen)) == NULL) PARSE_FAIL(SCORING_PARSE_NOMEM, -1);
  size_t nchars = 0;

  if(isspace((unsigned char)sep)) {
    for(size_t i = 1; i < llen; i++)
      if(!isspace((unsigned char)ls[i])) chars[nchars++] = fold_ch(ls[i], cs);
    // Rows: the row's character, then a score per column, whitespace between
    while((p = next_line(p, end, &ls, &llen)) != NULL) {
      line_num++;
      if(llen == 0 || ls[0] == '#' || span_all_ws(ls, llen)) continue;
      size_t idx = 0;
      while(idx < llen && isspace((unsigned char)ls[idx])) idx++;
      char from_c = fold_ch(ls[idx++], cs);
      for(size_t c = 0; c < nchars; c++) {
        if(idx >= llen) PARSE_FAIL(SCORING_PARSE_FEW_COLUMNS, line_num);
        if(!isspace((unsigned char)ls[idx]))
          PARSE_FAIL(SCORING_PARSE_EXPECTED_WS, line_num);
        while(idx < llen && isspace((unsigned char)ls[idx])) idx++;
        if(idx >= llen) PARSE_FAIL(SCORING_PARSE_FEW_COLUMNS, line_num);
        size_t start = idx;
        while(idx < llen && !isspace((unsigned char)ls[idx])) idx++;
        if((st = matrix_score(ls + start, idx - start, &score)) != SCORING_PARSE_OK)
          PARSE_FAIL(st, line_num);
        scoring_add_mutation(tmp, from_c, chars[c], score);
      }
      if(!span_all_ws(ls + idx, llen - idx))
        PARSE_FAIL(SCORING_PARSE_MANY_COLUMNS, line_num);
      rows_added++;
    }
  } else {
    // Header: <sep><char><sep><char>...
    size_t i;
    for(i = 0; i + 1 < llen; i += 2) {
      if(ls[i] != sep) PARSE_FAIL(SCORING_PARSE_SEP_MISSING, line_num);
      chars[nchars++] = fold_ch(ls[i + 1], cs);
    }
    if(i < llen) PARSE_FAIL(SCORING_PARSE_EMPTY_HEADING, line_num);
    // Rows: <char><sep><score><sep><score>...
    while((p = next_line(p, end, &ls, &llen)) != NULL) {
      line_num++;
      if(llen == 0 || ls[0] == '#' || span_all_ws(ls, llen)) continue;
      char from_c = fold_ch(ls[0], cs);
      size_t idx = 1;
      for(size_t c = 0; c < nchars; c++) {
        if(idx >= llen) PARSE_FAIL(SCORING_PARSE_FEW_COLUMNS, line_num);
        if(ls[idx] != sep) PARSE_FAIL(SCORING_PARSE_SEP_MISSING, line_num);
        idx++;
        size_t start = idx;
        while(idx < llen && ls[idx] != sep) idx++;
        if((st = matrix_score(ls + start, idx - start, &score)) != SCORING_PARSE_OK)
          PARSE_FAIL(st, line_num);
        scoring_add_mutation(tmp, from_c, chars[c], score);
      }
      if(idx < llen) PARSE_FAIL(SCORING_PARSE_MANY_COLUMNS, line_num);
      rows_added++;
    }
  }

  if(rows_added == 0) PARSE_FAIL(SCORING_PARSE_NO_ROWS, -1);
  free(chars);
  *scoring = *tmp;
  free(tmp);
  if(line != NULL) *line = -1;
  return SCORING_PARSE_OK;

fail:
  free(chars);
  free(tmp);
  if(line != NULL) *line = fail_line;
  return status;
}

scoring_parse_t scoring_parse_pairs(scoring_t *scoring, const char *text,
                                    size_t len, int *line)
{
  scoring_parse_t status = SCORING_PARSE_OK, st;
  int fail_line = -1;
  scoring_t *tmp = NULL;

  if((tmp = aln_malloc(sizeof(*tmp))) == NULL) PARSE_FAIL(SCORING_PARSE_NOMEM, -1);
  *tmp = *scoring;
  const bool cs = tmp->case_sensitive;
  const char *p = text, *end = text + len;
  const char *ls = NULL;
  size_t llen = 0;
  int line_num = 0, added = 0;

  while((p = next_line(p, end, &ls, &llen)) != NULL) {
    line_num++;
    if(llen == 0 || ls[0] == '#' || span_all_ws(ls, llen)) continue;
    if(llen < 5) PARSE_FAIL(SCORING_PARSE_SHORT_LINE, line_num);
    char a, b;
    int score = 0;
    if(isspace((unsigned char)ls[1])) {
      // a <whitespace> b <whitespace> score
      a = ls[0];
      size_t pos = 1;
      while(pos < llen && isspace((unsigned char)ls[pos])) pos++;
      if(pos + 2 >= llen || !isspace((unsigned char)ls[pos + 1]))
        PARSE_FAIL(SCORING_PARSE_SHORT_LINE, line_num);
      b = ls[pos];
      if((st = pairs_score(ls + pos + 2, llen - (pos + 2), &score)) != SCORING_PARSE_OK)
        PARSE_FAIL(st, line_num);
    } else {
      // a<sep>b<sep>score
      if(ls[1] != ls[3]) PARSE_FAIL(SCORING_PARSE_SEPARATORS, line_num);
      a = ls[0];
      b = ls[2];
      if((st = pairs_score(ls + 4, llen - 4, &score)) != SCORING_PARSE_OK)
        PARSE_FAIL(st, line_num);
    }
    scoring_add_mutation(tmp, fold_ch(a, cs), fold_ch(b, cs), score);
    added++;
  }

  if(added == 0) PARSE_FAIL(SCORING_PARSE_NO_PAIRS, -1);
  *scoring = *tmp;
  free(tmp);
  if(line != NULL) *line = -1;
  return SCORING_PARSE_OK;

fail:
  free(tmp);
  if(line != NULL) *line = fail_line;
  return status;
}
