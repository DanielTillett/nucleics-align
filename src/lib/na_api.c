/*
 src/lib/na_api.c
 Implementation of the public Nucleics Align v2.0 C11 API (include/nucleics_align.h).
 Zero exit() calls, zero zlib dependencies, zero mutable global state mutations.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include "nucleics_align.h"
#include "alignment.h"
#include "needleman_wunsch.h"
#include "smith_waterman.h"
#include "hirschberg.h"
#include "alignment_alloc.h"

struct na_scoring {
  scoring_t s;
};

struct na_alignment {
  alignment_t *a;
};

struct na_aligner {
  nw_aligner_t *nw;
  sw_aligner_t *sw;
  alignment_options_t opts;
  na_method_t method;
  score_t min_local_score;
  bool local_active;
  bool has_full_matrices;
  bool last_was_local;
  na_error_t last_err;
};

static void na_set_err(na_error_t *err, na_status_t st, char a, char b,
                       int line, const char *msg)
{
  if(err == NULL) return;
  err->status = st;
  err->a = a;
  err->b = b;
  err->line = line;
  if(msg != NULL) {
    size_t n = strlen(msg);
    if(n >= sizeof(err->msg)) n = sizeof(err->msg) - 1;
    memcpy(err->msg, msg, n);
    err->msg[n] = '\0';
  } else {
    err->msg[0] = '\0';
  }
}

static na_status_t map_status(aln_status_t st)
{
  switch(st) {
    case ALN_OK:              return NA_OK;
    case ALN_ERR_NOMEM:       return NA_ERR_NOMEM;
    case ALN_ERR_TOO_LONG:    return NA_ERR_TOO_LONG;
    case ALN_ERR_OVERFLOW:    return NA_ERR_OVERFLOW;
    case ALN_ERR_SCORING:     return NA_ERR_SCORING_PAIR;
    case ALN_ERR_UNSUPPORTED: return NA_ERR_UNSUPPORTED;
    default:                  return NA_ERR_INTERNAL;
  }
}

static na_status_t sync_err(na_aligner_t *al, aln_status_t st, const aln_error_t *ae)
{
  na_status_t ns = map_status(st);
  if(al != NULL) {
    if(ae != NULL)
      na_set_err(&al->last_err, ns, ae->a, ae->b, -1, NULL);
    else
      na_set_err(&al->last_err, ns, 0, 0, -1, NULL);
  }
  return ns;
}

const char *na_version(void)
{
  return NA_VERSION_STR;
}

const char *na_strerror(na_status_t status)
{
  switch(status) {
    case NA_OK:               return "ok";
    case NA_DONE:             return "done";
    case NA_ERR_INVALID_ARG:  return "invalid argument";
    case NA_ERR_SCORING_PAIR: return "no substitution score for character pair";
    case NA_ERR_TOO_LONG:     return "sequences too long";
    case NA_ERR_OVERFLOW:     return "alignment score may overflow";
    case NA_ERR_UNSUPPORTED:  return "unsupported scoring options for method";
    case NA_ERR_NOMEM:        return "out of memory";
    case NA_ERR_PARSE:        return "parse error";
    case NA_ERR_INTERNAL:     return "internal error";
    default:                  return "unknown error";
  }
}

size_t na_error_format(const na_error_t *err, char *buf, size_t buf_size)
{
  if(err == NULL || buf == NULL || buf_size == 0) return 0;
  int n = 0;
  if(err->status == NA_ERR_SCORING_PAIR) {
    n = snprintf(buf, buf_size,
                 "Error: no score for aligning '%c' against '%c'",
                 err->a, err->b);
  } else if(err->status == NA_ERR_PARSE) {
    if(err->line > 0)
      n = snprintf(buf, buf_size, "Error (line %d): %s", err->line,
                   err->msg[0] ? err->msg : "parse error");
    else
      n = snprintf(buf, buf_size, "Error: %s",
                   err->msg[0] ? err->msg : "parse error");
  } else if(err->msg[0] != '\0') {
    n = snprintf(buf, buf_size, "Error: %s", err->msg);
  } else {
    n = snprintf(buf, buf_size, "Error: %s", na_strerror(err->status));
  }
  if(n < 0) { buf[0] = '\0'; return 0; }
  return (size_t)n >= buf_size ? buf_size - 1 : (size_t)n;
}

/* =========================================================================
 * Scoring
 * ========================================================================= */

static bool valid_score(int v)
{
  return v >= -SCORE_ABS_MAX && v <= SCORE_ABS_MAX;
}

static int na_strcasecmp(const char *a, const char *b)
{
  while(*a && *b) {
    int ca = tolower((unsigned char)*a++);
    int cb = tolower((unsigned char)*b++);
    if(ca != cb) return ca - cb;
  }
  return tolower((unsigned char)*a) - tolower((unsigned char)*b);
}

na_scoring_t *na_scoring_new(int match, int mismatch, int gap_open, int gap_extend)
{
  if(!valid_score(match) || !valid_score(mismatch) ||
     !valid_score(gap_open) || !valid_score(gap_extend))
    return NULL;
  na_scoring_t *sc = aln_calloc(1, sizeof(*sc));
  if(sc == NULL) return NULL;
  scoring_init(&sc->s, match, mismatch, gap_open, gap_extend,
               false, false, false, false, false, false);
  return sc;
}

na_scoring_t *na_scoring_new_preset(const char *name)
{
  if(name == NULL) return NULL;
  na_scoring_t *sc = aln_calloc(1, sizeof(*sc));
  if(sc == NULL) return NULL;
  if(na_strcasecmp(name, "DNA_DEFAULT") == 0 || na_strcasecmp(name, "DEFAULT") == 0) {
    scoring_system_default(&sc->s);
    return sc;
  }
  if(!scoring_system_by_name(&sc->s, name)) {
    free(sc);
    return NULL;
  }
  return sc;
}

void na_scoring_free(na_scoring_t *sc)
{
  free(sc);
}

const char *na_scoring_preset_names(void)
{
  return scoring_system_names();
}

na_status_t na_scoring_set_flags(na_scoring_t *sc, unsigned flags)
{
  if(sc == NULL) return NA_ERR_INVALID_ARG;
  const unsigned valid_mask = NA_FREE_START_GAP | NA_FREE_END_GAP |
                              NA_NO_GAPS_A | NA_NO_GAPS_B |
                              NA_NO_GAPS_LONGEST | NA_NO_MISMATCHES |
                              NA_CASE_SENSITIVE;
  if((flags & ~valid_mask) != 0) return NA_ERR_INVALID_ARG;

  sc->s.no_start_gap_penalty = (flags & NA_FREE_START_GAP) != 0;
  sc->s.no_end_gap_penalty   = (flags & NA_FREE_END_GAP) != 0;
  sc->s.no_gaps_in_a         = (flags & NA_NO_GAPS_A) != 0;
  sc->s.no_gaps_in_b         = (flags & NA_NO_GAPS_B) != 0;
  sc->s.no_gaps_in_longest   = (flags & NA_NO_GAPS_LONGEST) != 0;
  sc->s.no_mismatches        = (flags & NA_NO_MISMATCHES) != 0;
  bool cs = (flags & NA_CASE_SENSITIVE) != 0;
  if(!cs && sc->s.case_sensitive)
    return NA_ERR_INVALID_ARG;
  if(cs && !sc->s.case_sensitive)
    scoring_make_case_sensitive(&sc->s);
  return NA_OK;
}

na_status_t na_scoring_set_match_mismatch(na_scoring_t *sc, bool enable,
                                          int match, int mismatch)
{
  if(sc == NULL || !valid_score(match) || !valid_score(mismatch))
    return NA_ERR_INVALID_ARG;
  sc->s.use_match_mismatch = enable;
  sc->s.match = match;
  sc->s.mismatch = mismatch;
  if(match < sc->s.min_penalty) sc->s.min_penalty = match;
  if(mismatch < sc->s.min_penalty) sc->s.min_penalty = mismatch;
  if(match > sc->s.max_penalty) sc->s.max_penalty = match;
  if(mismatch > sc->s.max_penalty) sc->s.max_penalty = mismatch;
  return NA_OK;
}

na_status_t na_scoring_set_gaps(na_scoring_t *sc, int gap_open, int gap_extend)
{
  if(sc == NULL || !valid_score(gap_open) || !valid_score(gap_extend))
    return NA_ERR_INVALID_ARG;
  sc->s.gap_open = gap_open;
  sc->s.gap_extend = gap_extend;
  int go = gap_open + gap_extend;
  if(go < sc->s.min_penalty) sc->s.min_penalty = go;
  if(gap_extend < sc->s.min_penalty) sc->s.min_penalty = gap_extend;
  if(go > sc->s.max_penalty) sc->s.max_penalty = go;
  if(gap_extend > sc->s.max_penalty) sc->s.max_penalty = gap_extend;
  return NA_OK;
}

na_status_t na_scoring_set_pair(na_scoring_t *sc, char a, char b, int score)
{
  if(sc == NULL || !valid_score(score)) return NA_ERR_INVALID_ARG;
  scoring_add_mutation(&sc->s, a, b, score);
  return NA_OK;
}

na_status_t na_scoring_set_wildcard(na_scoring_t *sc, char c, int score)
{
  if(sc == NULL || !valid_score(score)) return NA_ERR_INVALID_ARG;
  scoring_add_wildcard(&sc->s, c, score);
  return NA_OK;
}

na_status_t na_scoring_lookup(const na_scoring_t *sc, char a, char b,
                              int *score_out, bool *is_match_out)
{
  if(sc == NULL || score_out == NULL || is_match_out == NULL)
    return NA_ERR_INVALID_ARG;
  if(!scoring_has_pair(&sc->s, a, b))
    return NA_ERR_SCORING_PAIR;
  scoring_lookup(&sc->s, a, b, score_out, is_match_out);
  return NA_OK;
}

/* --- In-memory scoring parsers (transactional, zero exit(), zero zlib) --- */

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

static bool parse_int_span(const char *s, size_t len, int *out)
{
  if(memchr(s, '\0', len) != NULL) return false;
  while(len > 0 && isspace((unsigned char)*s)) { s++; len--; }
  while(len > 0 && isspace((unsigned char)s[len - 1])) len--;
  if(len == 0 || len >= 32) return false;
  char tmp[32];
  memcpy(tmp, s, len);
  tmp[len] = '\0';
  char *ep = tmp;
  long v = strtol(tmp, &ep, 10);
  if(ep != tmp + len) return false;
  if(v < -SCORE_ABS_MAX || v > SCORE_ABS_MAX) return false;
  *out = (int)v;
  return true;
}

na_status_t na_scoring_parse_matrix(na_scoring_t *sc, const char *text,
                                    size_t len, na_error_t *err)
{
  if(sc == NULL || text == NULL) {
    na_set_err(err, NA_ERR_INVALID_ARG, 0, 0, -1, "invalid argument");
    return NA_ERR_INVALID_ARG;
  }
  scoring_t *tmp = aln_malloc(sizeof(*tmp));
  if(tmp == NULL) {
    na_set_err(err, NA_ERR_NOMEM, 0, 0, -1, "out of memory");
    return NA_ERR_NOMEM;
  }
  *tmp = sc->s;
  const bool cs = tmp->case_sensitive;
  const char *p = text, *end = text + len;
  const char *ls = NULL;
  size_t llen = 0;
  int line_num = 0, rows_added = 0;

  // Find header line
  bool found = false;
  while((p = next_line(p, end, &ls, &llen)) != NULL) {
    line_num++;
    if(llen > 0 && ls[0] != '#' && !span_all_ws(ls, llen)) {
      if(llen < 2) {
        free(tmp);
        na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "too few column headings");
        return NA_ERR_PARSE;
      }
      found = true;
      break;
    }
  }
  if(!found) {
    free(tmp);
    na_set_err(err, NA_ERR_PARSE, 0, 0, -1, "empty matrix");
    return NA_ERR_PARSE;
  }

  char sep = ls[0];
  if((sep >= '0' && sep <= '9') || sep == '-') {
    free(tmp);
    na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "invalid separator");
    return NA_ERR_PARSE;
  }

  char *chars = aln_malloc(llen);
  if(chars == NULL) {
    free(tmp);
    na_set_err(err, NA_ERR_NOMEM, 0, 0, -1, "out of memory");
    return NA_ERR_NOMEM;
  }
  size_t nchars = 0;

  if(isspace((unsigned char)sep)) {
    for(size_t i = 1; i < llen; i++) {
      if(!isspace((unsigned char)ls[i])) {
        chars[nchars++] = fold_ch(ls[i], cs);
      }
    }
    while((p = next_line(p, end, &ls, &llen)) != NULL) {
      line_num++;
      if(llen == 0 || ls[0] == '#' || span_all_ws(ls, llen)) continue;
      size_t idx = 0;
      while(idx < llen && isspace((unsigned char)ls[idx])) idx++;
      char from_c = fold_ch(ls[idx++], cs);
      for(size_t c = 0; c < nchars; c++) {
        if(idx >= llen || !isspace((unsigned char)ls[idx])) {
          free(chars); free(tmp);
          na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "too few columns or missing whitespace");
          return NA_ERR_PARSE;
        }
        while(idx < llen && isspace((unsigned char)ls[idx])) idx++;
        size_t start = idx;
        while(idx < llen && !isspace((unsigned char)ls[idx])) idx++;
        int v = 0;
        if(idx == start || !parse_int_span(ls + start, idx - start, &v)) {
          free(chars); free(tmp);
          na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "invalid or out-of-range score");
          return NA_ERR_PARSE;
        }
        scoring_add_mutation(tmp, from_c, chars[c], v);
      }
      if(!span_all_ws(ls + idx, llen - idx)) {
        free(chars); free(tmp);
        na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "too many columns on row");
        return NA_ERR_PARSE;
      }
      rows_added++;
    }
  } else {
    size_t i;
    for(i = 0; i + 1 < llen; i += 2) {
      if(ls[i] != sep) {
        free(chars); free(tmp);
        na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "separator missing from header");
        return NA_ERR_PARSE;
      }
      chars[nchars++] = fold_ch(ls[i + 1], cs);
    }
    if(i < llen) {
      free(chars); free(tmp);
      na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "empty column heading");
      return NA_ERR_PARSE;
    }
    while((p = next_line(p, end, &ls, &llen)) != NULL) {
      line_num++;
      if(llen == 0 || ls[0] == '#' || span_all_ws(ls, llen)) continue;
      char from_c = fold_ch(ls[0], cs);
      size_t idx = 1;
      for(size_t c = 0; c < nchars; c++) {
        if(idx >= llen || ls[idx] != sep) {
          free(chars); free(tmp);
          na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "separator missing or too few columns");
          return NA_ERR_PARSE;
        }
        idx++;
        size_t start = idx;
        while(idx < llen && ls[idx] != sep) idx++;
        int v = 0;
        if(!parse_int_span(ls + start, idx - start, &v)) {
          free(chars); free(tmp);
          na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "invalid or out-of-range score");
          return NA_ERR_PARSE;
        }
        scoring_add_mutation(tmp, from_c, chars[c], v);
      }
      if(idx < llen) {
        free(chars); free(tmp);
        na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "too many columns on row");
        return NA_ERR_PARSE;
      }
      rows_added++;
    }
  }

  free(chars);
  if(rows_added == 0) {
    free(tmp);
    na_set_err(err, NA_ERR_PARSE, 0, 0, -1, "no rows in matrix");
    return NA_ERR_PARSE;
  }
  sc->s = *tmp;
  free(tmp);
  na_set_err(err, NA_OK, 0, 0, -1, NULL);
  return NA_OK;
}

na_status_t na_scoring_parse_pairs(na_scoring_t *sc, const char *text,
                                   size_t len, na_error_t *err)
{
  if(sc == NULL || text == NULL) {
    na_set_err(err, NA_ERR_INVALID_ARG, 0, 0, -1, "invalid argument");
    return NA_ERR_INVALID_ARG;
  }
  scoring_t *tmp = aln_malloc(sizeof(*tmp));
  if(tmp == NULL) {
    na_set_err(err, NA_ERR_NOMEM, 0, 0, -1, "out of memory");
    return NA_ERR_NOMEM;
  }
  *tmp = sc->s;
  const bool cs = tmp->case_sensitive;
  const char *p = text, *end = text + len;
  const char *ls = NULL;
  size_t llen = 0;
  int line_num = 0, added = 0;

  while((p = next_line(p, end, &ls, &llen)) != NULL) {
    line_num++;
    if(llen == 0 || ls[0] == '#' || span_all_ws(ls, llen)) continue;
    if(llen < 5) {
      free(tmp);
      na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "line too short");
      return NA_ERR_PARSE;
    }
    char a, b;
    int score = 0;
    if(isspace((unsigned char)ls[1])) {
      a = ls[0];
      size_t pos = 1;
      while(pos < llen && isspace((unsigned char)ls[pos])) pos++;
      if(pos + 2 >= llen || !isspace((unsigned char)ls[pos + 1])) {
        free(tmp);
        na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "line too short");
        return NA_ERR_PARSE;
      }
      b = ls[pos];
      if(!parse_int_span(ls + pos + 2, llen - (pos + 2), &score)) {
        free(tmp);
        na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "invalid or out-of-range score");
        return NA_ERR_PARSE;
      }
    } else {
      if(ls[1] != ls[3]) {
        free(tmp);
        na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "inconsistent separators");
        return NA_ERR_PARSE;
      }
      a = ls[0];
      b = ls[2];
      if(!parse_int_span(ls + 4, llen - 4, &score)) {
        free(tmp);
        na_set_err(err, NA_ERR_PARSE, 0, 0, line_num, "invalid or out-of-range score");
        return NA_ERR_PARSE;
      }
    }
    scoring_add_mutation(tmp, fold_ch(a, cs), fold_ch(b, cs), score);
    added++;
  }

  if(added == 0) {
    free(tmp);
    na_set_err(err, NA_ERR_PARSE, 0, 0, -1, "no pairs added");
    return NA_ERR_PARSE;
  }
  sc->s = *tmp;
  free(tmp);
  na_set_err(err, NA_OK, 0, 0, -1, NULL);
  return NA_OK;
}

/* =========================================================================
 * Alignment Result
 * ========================================================================= */

na_alignment_t *na_alignment_new(void)
{
  na_alignment_t *res = aln_calloc(1, sizeof(*res));
  if(res == NULL) return NULL;
  res->a = alignment_try_create(256);
  if(res->a == NULL) { free(res); return NULL; }
  return res;
}

void na_alignment_free(na_alignment_t *res)
{
  if(res == NULL) return;
  alignment_free(res->a);
  free(res);
}

const char *na_alignment_a(const na_alignment_t *res)
{
  return (res && res->a) ? res->a->result_a : "";
}

const char *na_alignment_b(const na_alignment_t *res)
{
  return (res && res->a) ? res->a->result_b : "";
}

size_t na_alignment_length(const na_alignment_t *res)
{
  return (res && res->a) ? res->a->length : 0;
}

na_score_t na_alignment_score(const na_alignment_t *res)
{
  return (res && res->a) ? (na_score_t)res->a->score : 0;
}

void na_alignment_span(const na_alignment_t *res,
                       size_t *pos_a, size_t *len_a,
                       size_t *pos_b, size_t *len_b)
{
  if(res == NULL || res->a == NULL) return;
  if(pos_a) *pos_a = res->a->pos_a;
  if(len_a) *len_a = res->a->len_a;
  if(pos_b) *pos_b = res->a->pos_b;
  if(len_b) *len_b = res->a->len_b;
}

/* =========================================================================
 * Aligner Context
 * ========================================================================= */

na_aligner_t *na_aligner_new(void)
{
  na_aligner_t *al = aln_calloc(1, sizeof(*al));
  if(al == NULL) return NULL;
  al->nw = needleman_wunsch_try_new();
  al->sw = smith_waterman_try_new();
  if(al->nw == NULL || al->sw == NULL) {
    if(al->nw) needleman_wunsch_free(al->nw);
    if(al->sw) smith_waterman_free(al->sw);
    free(al);
    return NULL;
  }
  alignment_options_default(&al->opts);
  al->method = NA_METHOD_AUTO;
  al->min_local_score = 1;
  al->nw->compact = true;
  al->nw->options = &al->opts;
  smith_waterman_set_compact(al->sw, true);
  smith_waterman_set_options(al->sw, &al->opts);
  return al;
}

void na_aligner_free(na_aligner_t *al)
{
  if(al == NULL) return;
  needleman_wunsch_free(al->nw);
  smith_waterman_free(al->sw);
  free(al);
}

void na_aligner_set_threads(na_aligner_t *al, unsigned threads)
{
  if(al == NULL) return;
  int t = threads > (unsigned)INT_MAX ? INT_MAX : (int)threads;
  al->opts.max_threads = t;
  al->opts.hb_max_threads = t;
}

void na_aligner_set_method(na_aligner_t *al, na_method_t method)
{
  if(al == NULL) return;
  al->method = method;
  al->opts.wfa_align = (method == NA_METHOD_WFA);
  al->nw->compact = (method != NA_METHOD_FULL);
}

void na_aligner_set_linear_above(na_aligner_t *al, size_t max_len)
{
  if(al == NULL) return;
  al->opts.linear_above = max_len;
}

void na_aligner_set_linear_fallback(na_aligner_t *al, bool enable)
{
  if(al == NULL) return;
  al->opts.linear_fallback = enable;
}

void na_aligner_set_wfa_budget(na_aligner_t *al, size_t budget)
{
  if(al == NULL) return;
  al->opts.wfa_budget = budget;
}

void na_aligner_set_min_local_score(na_aligner_t *al, na_score_t min_score)
{
  if(al == NULL) return;
  al->min_local_score = min_score < 1 ? 1 : (score_t)min_score;
  smith_waterman_set_min_score(al->sw, al->min_local_score);
}

void na_aligner_set_fast_local(na_aligner_t *al, bool enable)
{
  if(al == NULL) return;
  al->opts.sw_fast_rect = enable;
}

const na_error_t *na_aligner_last_error(const na_aligner_t *al)
{
  return al ? &al->last_err : NULL;
}

/* =========================================================================
 * Global & Local Alignment Dispatch
 * ========================================================================= */

static void clear_alignment(alignment_t *a)
{
  if(a == NULL) return;
  a->length = 0;
  a->score = 0;
  a->pos_a = a->pos_b = a->len_a = a->len_b = 0;
  if(a->result_a) a->result_a[0] = '\0';
  if(a->result_b) a->result_b[0] = '\0';
}

na_status_t na_global_align(na_aligner_t *al,
                            const char *a, size_t len_a,
                            const char *b, size_t len_b,
                            const na_scoring_t *sc,
                            na_alignment_t *out)
{
  if(al != NULL) { al->local_active = false; al->has_full_matrices = false; }
  if(out != NULL) clear_alignment(out->a);
  if(al == NULL || sc == NULL || out == NULL || out->a == NULL ||
     (a == NULL && len_a > 0) || (b == NULL && len_b > 0)) {
    if(al) na_set_err(&al->last_err, NA_ERR_INVALID_ARG, 0, 0, -1, "invalid argument");
    return NA_ERR_INVALID_ARG;
  }
  if(sc->s.no_mismatches &&
     (sc->s.no_gaps_in_a || sc->s.no_gaps_in_b || sc->s.no_gaps_in_longest)) {
    na_set_err(&al->last_err, NA_ERR_UNSUPPORTED, 0, 0, -1,
               "cannot combine NA_NO_MISMATCHES with NA_NO_GAPS_*");
    return NA_ERR_UNSUPPORTED;
  }
  const char *sa = a ? a : "";
  const char *sb = b ? b : "";
  al->last_was_local = false;
  aln_error_t ae = {0};

  if(al->method == NA_METHOD_LINEAR) {
    aln_status_t st = hirschberg_align_checked(sa, sb, len_a, len_b,
                                               &sc->s, &al->opts, out->a, &ae);
    return sync_err(al, st, &ae);
  }

  alignment_options_t eff = al->opts;
  if(al->method == NA_METHOD_FULL) {
    eff.linear_above = 0;
    eff.linear_fallback = false;
    eff.wfa_align = false;
  }
  al->nw->compact = (al->method != NA_METHOD_FULL);
  al->nw->options = &eff;
  aln_status_t st = needleman_wunsch_align_checked(sa, sb, len_a, len_b,
                                                   &sc->s, al->nw, out->a, &ae);
  al->nw->options = &al->opts;
  if(st == ALN_OK && al->method == NA_METHOD_FULL &&
     al->nw->match_scores != NULL && !al->nw->traced)
    al->has_full_matrices = true;
  return sync_err(al, st, &ae);
}

na_status_t na_global_score(na_aligner_t *al,
                            const char *a, size_t len_a,
                            const char *b, size_t len_b,
                            const na_scoring_t *sc,
                            na_score_t *out_score)
{
  if(al != NULL) { al->local_active = false; al->has_full_matrices = false; }
  if(al == NULL || sc == NULL || out_score == NULL ||
     (a == NULL && len_a > 0) || (b == NULL && len_b > 0)) {
    if(al) na_set_err(&al->last_err, NA_ERR_INVALID_ARG, 0, 0, -1, "invalid argument");
    return NA_ERR_INVALID_ARG;
  }
  if(sc->s.no_mismatches &&
     (sc->s.no_gaps_in_a || sc->s.no_gaps_in_b || sc->s.no_gaps_in_longest)) {
    na_set_err(&al->last_err, NA_ERR_UNSUPPORTED, 0, 0, -1,
               "cannot combine NA_NO_MISMATCHES with NA_NO_GAPS_*");
    return NA_ERR_UNSUPPORTED;
  }
  const char *sa = a ? a : "";
  const char *sb = b ? b : "";
  al->last_was_local = false;
  aln_error_t ae = {0};
  score_t s = 0;

  aln_status_t st = hirschberg_score_checked(sa, sb, len_a, len_b,
                                             &sc->s, &al->opts, &s, &ae);
  if(st == ALN_OK) {
    *out_score = (na_score_t)s;
    return sync_err(al, ALN_OK, &ae);
  }
  if(st == ALN_ERR_UNSUPPORTED) {
    alignment_t *tmp = alignment_try_create(64);
    if(tmp == NULL) return sync_err(al, ALN_ERR_NOMEM, NULL);
    al->nw->compact = true;
    al->nw->options = &al->opts;
    st = needleman_wunsch_align_checked(sa, sb, len_a, len_b,
                                        &sc->s, al->nw, tmp, &ae);
    if(st == ALN_OK) *out_score = (na_score_t)tmp->score;
    alignment_free(tmp);
    return sync_err(al, st, &ae);
  }
  return sync_err(al, st, &ae);
}

na_status_t na_local_align(na_aligner_t *al,
                           const char *a, size_t len_a,
                           const char *b, size_t len_b,
                           const na_scoring_t *sc)
{
  if(al != NULL) { al->local_active = false; al->has_full_matrices = false; }
  if(al == NULL || sc == NULL ||
     (a == NULL && len_a > 0) || (b == NULL && len_b > 0)) {
    if(al) na_set_err(&al->last_err, NA_ERR_INVALID_ARG, 0, 0, -1, "invalid argument");
    return NA_ERR_INVALID_ARG;
  }
  const char *sa = a ? a : "";
  const char *sb = b ? b : "";
  al->last_was_local = true;
  smith_waterman_set_compact(al->sw, al->method != NA_METHOD_FULL);
  smith_waterman_set_min_score(al->sw, al->min_local_score);
  smith_waterman_set_options(al->sw, &al->opts);

  aln_error_t ae = {0};
  aln_status_t st = smith_waterman_align_checked(sa, sb, len_a, len_b,
                                                 &sc->s, al->sw, &ae);
  if(st == ALN_OK) {
    al->local_active = true;
    const aligner_t *raw = smith_waterman_get_aligner(al->sw);
    if(al->method == NA_METHOD_FULL && raw != NULL &&
       raw->match_scores != NULL && !raw->traced)
      al->has_full_matrices = true;
  }
  return sync_err(al, st, &ae);
}

na_status_t na_local_next(na_aligner_t *al, na_alignment_t *out)
{
  if(out != NULL) clear_alignment(out->a);
  if(al == NULL || out == NULL || out->a == NULL) {
    if(al) na_set_err(&al->last_err, NA_ERR_INVALID_ARG, 0, 0, -1, "invalid argument");
    return NA_ERR_INVALID_ARG;
  }
  if(!al->local_active) {
    na_set_err(&al->last_err, NA_DONE, 0, 0, -1, NULL);
    return NA_DONE;
  }
  aln_error_t ae = {0};
  int got = smith_waterman_fetch_checked(al->sw, out->a, &ae);
  if(got > 0) return sync_err(al, ALN_OK, &ae);
  al->local_active = false;
  if(ae.status != ALN_OK) return sync_err(al, ae.status, &ae);
  na_set_err(&al->last_err, NA_DONE, 0, 0, -1, NULL);
  return NA_DONE;
}

na_status_t na_score_alignment(const char *res_a, const char *res_b, size_t len,
                               const na_scoring_t *sc, na_score_t *out_score)
{
  if(res_a == NULL || res_b == NULL || sc == NULL || out_score == NULL)
    return NA_ERR_INVALID_ARG;
  for(size_t i = 0; i < len; i++) {
    if(res_a[i] != '-' && res_b[i] != '-') {
      if(!scoring_has_pair(&sc->s, res_a[i], res_b[i]))
        return NA_ERR_SCORING_PAIR;
    }
  }
  // Letters of each sequence, to tell the DP's free edges: a gap in a (a
  // column consuming b) is in the first column while no letter of a has
  // been consumed, and in the last column once all have; a gap in b is in
  // the first row while no letter of b has been consumed, and in the last
  // row once all have. Those gaps are free under the start/end flags
  size_t letters_a = 0, letters_b = 0;
  for(size_t i = 0; i < len; i++) {
    letters_a += res_a[i] != '-';
    letters_b += res_b[i] != '-';
  }
  const bool free_start = sc->s.no_start_gap_penalty;
  const bool free_end = sc->s.no_end_gap_penalty;

  long long total = 0;
  size_t x = 0, y = 0; // letters of a and b consumed so far
  bool in_a = false, in_b = false;
  for(size_t i = 0; i < len; i++) {
    char ca = res_a[i], cb = res_b[i];
    if(ca == '-' && cb == '-') return NA_ERR_INVALID_ARG;
    if(ca == '-') {
      const bool free_gap = (free_start && x == 0) ||
                            (free_end && x == letters_a);
      if(!free_gap)
        total += in_a ? sc->s.gap_extend : (sc->s.gap_open + sc->s.gap_extend);
      in_a = true; in_b = false;
      y++;
    } else if(cb == '-') {
      const bool free_gap = (free_start && y == 0) ||
                            (free_end && y == letters_b);
      if(!free_gap)
        total += in_b ? sc->s.gap_extend : (sc->s.gap_open + sc->s.gap_extend);
      in_b = true; in_a = false;
      x++;
    } else {
      int s = 0; bool m = false;
      scoring_lookup(&sc->s, ca, cb, &s, &m);
      if(sc->s.no_mismatches && !m) return NA_ERR_SCORING_PAIR;
      total += s;
      in_a = false; in_b = false;
      x++; y++;
    }
  }
  if(total < INT32_MIN || total > INT32_MAX) return NA_ERR_OVERFLOW;
  *out_score = (na_score_t)total;
  return NA_OK;
}

void na_debug_fprint_matrices(const na_aligner_t *al, FILE *fp)
{
  if(al == NULL || fp == NULL || !al->has_full_matrices) return;
  const aligner_t *raw = al->last_was_local ? smith_waterman_get_aligner(al->sw) : al->nw;
  if(raw == NULL || raw->match_scores == NULL || raw->traced) return;
  alignment_fprint_matrices(fp, raw);
}
