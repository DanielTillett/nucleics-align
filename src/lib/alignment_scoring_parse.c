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
 formats through here, so one grammar serves both: the tools' own, as
 their line reader has read it since 2013. Matrix scores are read with
 strtol(), so a score may carry a sign and leading spaces, and what
 follows it is checked as a separator or the end of the row; a pairs
 score must be the rest of its line
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

// The next line, as the tools' line reader gave it: a C string in *buf
// (grown as needed) holding the text up to the line feed, with the
// carriage returns and line feeds at its end removed and anything from a
// NUL byte on dropped. Returns where the line after starts, or NULL at the
// end (or if out of memory: *nomem)
static const char *next_line(const char *p, const char *end, char **buf,
                             size_t *cap, bool *nomem)
{
  if(p >= end) return NULL;
  const char *nl = memchr(p, '\n', (size_t)(end - p));
  const size_t len = (size_t)((nl ? nl : end) - p);
  if(len + 1 > *cap) {
    char *grown = realloc(*buf, len + 1);
    if(grown == NULL) { *nomem = true; return NULL; }
    *buf = grown;
    *cap = len + 1;
  }
  memcpy(*buf, p, len);
  (*buf)[len] = '\0';
  size_t n = strlen(*buf);
  while(n > 0 && ((*buf)[n-1] == '\r' || (*buf)[n-1] == '\n')) (*buf)[--n] = '\0';
  return nl ? nl + 1 : end;
}

static bool all_ws(const char *s)
{
  for(; *s != '\0'; s++) if(!isspace((unsigned char)*s)) return false;
  return true;
}

// The first character at or after s that is not whitespace, or NULL
static const char *next_nonws(const char *s)
{
  while(*s != '\0' && isspace((unsigned char)*s)) s++;
  return *s == '\0' ? NULL : s;
}

static char fold_ch(char c, bool cs)
{
  return cs ? c : (char)tolower((unsigned char)c);
}

// A matrix score: the number strtol() reads at s (*after: where it ends);
// NUMBER where there is none, RANGE past +-SCORE_ABS_MAX
static scoring_parse_t prefix_score(const char *s, const char **after,
                                    int *score)
{
  char *e = (char *)s;
  long v = strtol(s, &e, 10);
  *after = e;
  if(e == s) return SCORING_PARSE_NUMBER;
  if(v < -SCORE_ABS_MAX || v > SCORE_ABS_MAX) return SCORING_PARSE_RANGE;
  *score = (int)v;
  return SCORING_PARSE_OK;
}

// A pairs score: the whole of s must be one int (INVALID_NUMBER
// otherwise), within +-SCORE_ABS_MAX (RANGE)
static scoring_parse_t entire_score(const char *s, int *score)
{
  if(*s == '\0') return SCORING_PARSE_INVALID_NUMBER;
  char *e = (char *)s;
  errno = 0;
  long v = strtol(s, &e, 10);
  if(errno == ERANGE || v > INT_MAX || v < INT_MIN || *e != '\0')
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
#define NEXT_LINE() next_line(p, end, &line, &cap, &nomem)

scoring_parse_t scoring_parse_matrix(scoring_t *scoring, const char *text,
                                     size_t len, int *line_out)
{
  scoring_parse_t status = SCORING_PARSE_OK, st;
  int fail_line = -1;
  scoring_t *tmp = NULL;
  char *chars = NULL, *line = NULL;
  size_t cap = 0;
  bool nomem = false;

  if((tmp = aln_malloc(sizeof(*tmp))) == NULL) PARSE_FAIL(SCORING_PARSE_NOMEM, -1);
  *tmp = *scoring;
  const bool cs = tmp->case_sensitive;
  const char *p = text, *end = text + len;
  int line_num = 0, rows_added = 0, score;
  size_t nchars = 0, i;

  // The header line: the first that is not empty, a comment or blank
  bool found = false;
  while((p = NEXT_LINE()) != NULL) {
    line_num++;
    if(line[0] != '\0' && line[0] != '#' && !all_ws(line)) {
      if(strlen(line) < 2) PARSE_FAIL(SCORING_PARSE_FEW_HEADINGS, line_num);
      found = true;
      break;
    }
  }
  if(nomem) PARSE_FAIL(SCORING_PARSE_NOMEM, -1);
  if(!found) PARSE_FAIL(SCORING_PARSE_EMPTY, -1);

  // Its first character is the separator: whitespace (then any run of
  // whitespace separates), or one character
  const char sep = line[0];
  if((sep >= '0' && sep <= '9') || sep == '-')
    PARSE_FAIL(SCORING_PARSE_SEPARATOR, line_num);
  const size_t hlen = strlen(line);
  if((chars = aln_malloc(hlen)) == NULL) PARSE_FAIL(SCORING_PARSE_NOMEM, -1);

  if(isspace((unsigned char)sep)) {
    const char *next = line;
    while((next = next_nonws(next + 1)) != NULL) chars[nchars++] = fold_ch(*next, cs);
    // Rows: the row's character, then a score per column, whitespace between
    while((p = NEXT_LINE()) != NULL) {
      line_num++;
      const char *from = next_nonws(line);
      if(from == NULL || line[0] == '#') continue;
      const char from_c = fold_ch(*from, cs);
      const char *s = from + 1;
      for(i = 0; i < nchars; i++) {
        if(*s == '\0') PARSE_FAIL(SCORING_PARSE_FEW_COLUMNS, line_num);
        if(!isspace((unsigned char)*s)) PARSE_FAIL(SCORING_PARSE_EXPECTED_WS, line_num);
        if((s = next_nonws(s + 1)) == NULL) PARSE_FAIL(SCORING_PARSE_FEW_COLUMNS, line_num);
        if((st = prefix_score(s, &s, &score)) != SCORING_PARSE_OK) PARSE_FAIL(st, line_num);
        scoring_add_mutation(tmp, from_c, chars[i], score);
      }
      if(!all_ws(s)) PARSE_FAIL(SCORING_PARSE_MANY_COLUMNS, line_num);
      rows_added++;
    }
  } else {
    // Header: <sep><char><sep><char>...
    for(i = 0; i + 1 < hlen; i += 2) {
      if(line[i] != sep) PARSE_FAIL(SCORING_PARSE_SEP_MISSING, line_num);
      chars[nchars++] = fold_ch(line[i + 1], cs);
    }
    // e.g. ",A,C," - a separator with no heading after it
    if(i < hlen) PARSE_FAIL(SCORING_PARSE_EMPTY_HEADING, line_num);
    // Rows: <char><sep><score><sep><score>...
    while((p = NEXT_LINE()) != NULL) {
      line_num++;
      if(line[0] == '#' || all_ws(line)) continue;
      const char from_c = fold_ch(line[0], cs);
      const char *s = line + 1;
      size_t col = 0;
      while(*s != '\0') {
        if(col >= nchars) PARSE_FAIL(SCORING_PARSE_MANY_COLUMNS, line_num);
        if(*s != sep) PARSE_FAIL(SCORING_PARSE_SEP_MISSING, line_num);
        s++;
        if((st = prefix_score(s, &s, &score)) != SCORING_PARSE_OK) PARSE_FAIL(st, line_num);
        scoring_add_mutation(tmp, from_c, chars[col++], score);
      }
      if(col < nchars) PARSE_FAIL(SCORING_PARSE_FEW_COLUMNS, line_num);
      rows_added++;
    }
  }
  if(nomem) PARSE_FAIL(SCORING_PARSE_NOMEM, -1);

  if(rows_added == 0) PARSE_FAIL(SCORING_PARSE_NO_ROWS, -1);
  free(line);
  free(chars);
  *scoring = *tmp;
  free(tmp);
  if(line_out != NULL) *line_out = -1;
  return SCORING_PARSE_OK;

fail:
  free(line);
  free(chars);
  free(tmp);
  if(line_out != NULL) *line_out = fail_line;
  return status;
}

scoring_parse_t scoring_parse_pairs(scoring_t *scoring, const char *text,
                                    size_t len, int *line_out)
{
  scoring_parse_t status = SCORING_PARSE_OK, st;
  int fail_line = -1;
  scoring_t *tmp = NULL;
  char *line = NULL;
  size_t cap = 0;
  bool nomem = false;

  if((tmp = aln_malloc(sizeof(*tmp))) == NULL) PARSE_FAIL(SCORING_PARSE_NOMEM, -1);
  *tmp = *scoring;
  const bool cs = tmp->case_sensitive;
  const char *p = text, *end = text + len;
  int line_num = 0, added = 0;

  while((p = NEXT_LINE()) != NULL) {
    line_num++;
    if(line[0] == '\0' || line[0] == '#' || all_ws(line)) continue;
    const size_t n = strlen(line);
    if(n < 5) PARSE_FAIL(SCORING_PARSE_SHORT_LINE, line_num);
    char a, b;
    int score = 0;
    if(isspace((unsigned char)line[1])) {
      // a <whitespace> b <whitespace> score
      a = line[0];
      size_t pos = 1;
      while(line[pos] != '\0' && isspace((unsigned char)line[pos])) pos++;
      if(pos + 2 >= n || !isspace((unsigned char)line[pos + 1]))
        PARSE_FAIL(SCORING_PARSE_SHORT_LINE, line_num);
      b = line[pos];
      if((st = entire_score(line + pos + 2, &score)) != SCORING_PARSE_OK)
        PARSE_FAIL(st, line_num);
    } else {
      // a<sep>b<sep>score
      if(line[1] != line[3]) PARSE_FAIL(SCORING_PARSE_SEPARATORS, line_num);
      a = line[0];
      b = line[2];
      if((st = entire_score(line + 4, &score)) != SCORING_PARSE_OK)
        PARSE_FAIL(st, line_num);
    }
    scoring_add_mutation(tmp, fold_ch(a, cs), fold_ch(b, cs), score);
    added++;
  }
  if(nomem) PARSE_FAIL(SCORING_PARSE_NOMEM, -1);

  if(added == 0) PARSE_FAIL(SCORING_PARSE_NO_PAIRS, -1);
  free(line);
  *scoring = *tmp;
  free(tmp);
  if(line_out != NULL) *line_out = -1;
  return SCORING_PARSE_OK;

fail:
  free(line);
  free(tmp);
  if(line_out != NULL) *line_out = fail_line;
  return status;
}
