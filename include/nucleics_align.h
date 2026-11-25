/*
 include/nucleics_align.h
 Nucleics Align v2.0 — Public C11 / C++ API
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 */

#ifndef NUCLEICS_ALIGN_V2_HEADER_SEEN
#define NUCLEICS_ALIGN_V2_HEADER_SEEN

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

#define NA_VERSION_MAJOR 2
#define NA_VERSION_MINOR 0
#define NA_VERSION_PATCH 0
/* "2.0.0", from the three above */
#define NA_VERSION_STR_(x) #x
#define NA_VERSION_STR__(a, b, c) NA_VERSION_STR_(a.b.c)
#define NA_VERSION_STR NA_VERSION_STR__(NA_VERSION_MAJOR, NA_VERSION_MINOR, \
                                        NA_VERSION_PATCH)
/* One number for #if checks: 2.0.0 is 0x020000 */
#define NA_VERSION_NUMBER ((NA_VERSION_MAJOR << 16) | (NA_VERSION_MINOR << 8) \
                           | NA_VERSION_PATCH)

#ifdef __cplusplus
extern "C" {
#endif

typedef int32_t na_score_t;

/* Status codes returned by all fallible API calls. */
typedef enum {
  NA_OK = 0,
  NA_DONE,               /* Local hit iterator exhausted / no qualifying hit */
  NA_ERR_INVALID_ARG,    /* NULL pointer, invalid flag, or out-of-range value */
  NA_ERR_SCORING_PAIR,   /* Sequence contains a character pair with no score */
  NA_ERR_TOO_LONG,       /* Sequences too long for address space */
  NA_ERR_OVERFLOW,       /* Score could overflow int32_t */
  NA_ERR_UNSUPPORTED,    /* Requested method does not support scoring options */
  NA_ERR_NOMEM,          /* Out of memory */
  NA_ERR_PARSE,          /* Syntax or range error while parsing scoring text */
  NA_ERR_INTERNAL        /* Internal engine error */
} na_status_t;

/* Optional diagnostic detail populated on error. */
typedef struct {
  na_status_t status;
  char a;                /* Offending seq_a char (NA_ERR_SCORING_PAIR) */
  char b;                /* Offending seq_b char (NA_ERR_SCORING_PAIR) */
  int line;              /* 1-based line number (NA_ERR_PARSE), -1 if none */
  char msg[96];          /* Short human-readable error message */
} na_error_t;

const char *na_version(void);
const char *na_strerror(na_status_t status);
size_t      na_error_format(const na_error_t *err, char *buf, size_t buf_size);

/* =========================================================================
 * Scoring configuration (opaque; safe to share read-only across threads
 * once set up)
 * Gap penalty formula: a gap of length k >= 1 scores gap_open + k * gap_extend.
 * Scores must lie within -INT_MAX/4 .. INT_MAX/4 (NA_ERR_INVALID_ARG).
 * ========================================================================= */
typedef struct na_scoring na_scoring_t;

enum {
  NA_FREE_START_GAP  = 1u << 0,  /* a gap along the first row/column is free */
  NA_FREE_END_GAP    = 1u << 1,  /* a gap along the last row/column is free */
  NA_NO_GAPS_A       = 1u << 2,  /* no gap in a (global alignment only) */
  NA_NO_GAPS_B       = 1u << 3,  /* no gap in b */
  NA_NO_GAPS_LONGEST = 1u << 4,  /* no gap in the longer sequence */
  NA_NO_MISMATCHES   = 1u << 5,  /* no mismatched pair (not with NA_NO_GAPS_*) */
  NA_CASE_SENSITIVE  = 1u << 6   /* 'a' and 'A' differ; cannot be unset again */
};

/* Match/mismatch and gap scores; or a preset: one of
 * na_scoring_preset_names() ("PAM30", "PAM70", "BLOSUM80", "BLOSUM62",
 * "DNA_HYBRIDIZATION"), or "DEFAULT" / "DNA_DEFAULT" for the tools' default
 * DNA scoring. NULL for an unknown name or out of memory. */
na_scoring_t *na_scoring_new(int match, int mismatch, int gap_open, int gap_extend);
na_scoring_t *na_scoring_new_preset(const char *name);
void          na_scoring_free(na_scoring_t *sc);
const char   *na_scoring_preset_names(void);

/* Replaces every flag with the set given. NA_CASE_SENSITIVE is one-way: a
 * call without it after it was set is refused (NA_ERR_INVALID_ARG) and
 * changes nothing. Pairs and wildcards set before it are copied to every
 * mix of cases; those set after are exact. */
na_status_t na_scoring_set_flags(na_scoring_t *sc, unsigned flags);
na_status_t na_scoring_set_match_mismatch(na_scoring_t *sc, bool enable,
                                          int match, int mismatch);
na_status_t na_scoring_set_gaps(na_scoring_t *sc, int gap_open, int gap_extend);
/* A score for one pair, or for a character against anything (two
 * wildcards score the lower of their scores) */
na_status_t na_scoring_set_pair(na_scoring_t *sc, char a, char b, int score);
na_status_t na_scoring_set_wildcard(na_scoring_t *sc, char c, int score);
/* The file formats in scoring/, from memory; nothing changes on an error
 * (err->line says where it is) */
na_status_t na_scoring_parse_matrix(na_scoring_t *sc, const char *text,
                                    size_t len, na_error_t *err);
na_status_t na_scoring_parse_pairs(na_scoring_t *sc, const char *text,
                                   size_t len, na_error_t *err);
/* NA_ERR_SCORING_PAIR for a pair the scoring does not cover */
na_status_t na_scoring_lookup(const na_scoring_t *sc, char a, char b,
                              int *score_out, bool *is_match_out);

/* =========================================================================
 * Alignment result buffer (opaque; reusable across calls)
 * ========================================================================= */
typedef struct na_alignment na_alignment_t;

na_alignment_t *na_alignment_new(void);
void            na_alignment_free(na_alignment_t *res);

/* The aligned strings ('-' for a gap) belong to res: valid until the next
 * call that writes to res, or frees it. Copy them to keep them. With res
 * NULL: "", 0, no span. The span is the 0-based start and the letters
 * consumed in each sequence (0, len_a, 0, len_b for a global alignment).
 * Do not pass a result's own strings as input to a call writing to it. */
const char *na_alignment_a(const na_alignment_t *res);
const char *na_alignment_b(const na_alignment_t *res);
size_t      na_alignment_length(const na_alignment_t *res);
na_score_t  na_alignment_score(const na_alignment_t *res);
void        na_alignment_span(const na_alignment_t *res,
                              size_t *pos_a, size_t *len_a,
                              size_t *pos_b, size_t *len_b);

/* =========================================================================
 * Aligner context (opaque; one per thread, reusable across calls)
 * ========================================================================= */
typedef struct na_aligner na_aligner_t;

/* Every method gives the same optimal score. AUTO and FULL give the same
 * alignment (the documented tie-breaking); LINEAR and WFA, and AUTO once
 * linear_above or linear_fallback sends a pair to linear memory, may give
 * another alignment with that score. na_global_score() ignores the method. */
typedef enum {
  NA_METHOD_AUTO = 0,    /* Exact DP, one traceback byte per cell, banded
                            and SIMD where it applies */
  NA_METHOD_FULL,        /* The three full score matrices (12 bytes per
                            cell); needed for na_debug_fprint_matrices() */
  NA_METHOD_LINEAR,      /* Linear memory (Hirschberg / Myers-Miller);
                            refuses the free-gap, no-gap and no-mismatch
                            flags and a positive gap open: NA_ERR_UNSUPPORTED */
  NA_METHOD_WFA          /* Bidirectional wavefront alignment within the
                            wfa_budget, else as AUTO */
} na_method_t;

na_aligner_t *na_aligner_new(void);
void          na_aligner_free(na_aligner_t *al);

/* Settings apply to the calls after them. threads: within one alignment
 * (0: half the CPUs; at most 4096). linear_above: under AUTO, a pair whose
 * longer sequence has at least this many letters is aligned in linear
 * memory (0: never). linear_fallback: under AUTO, linear memory when the
 * matrices would not fit or ran out of memory. wfa_budget: diagonal
 * updates NA_METHOD_WFA may spend (0: min(2^28, 256 * (len_a + len_b))).
 * min_local_score: the lowest hit score na_local_next() fetches (default
 * 1; values below 1 are clamped to 1; the command line tools default to
 * match * max(0.2 * shorter length, 2) instead). fast_local: see
 * na_local_align(). */
void na_aligner_set_threads(na_aligner_t *al, unsigned threads);
void na_aligner_set_method(na_aligner_t *al, na_method_t method);
void na_aligner_set_linear_above(na_aligner_t *al, size_t max_len);
void na_aligner_set_linear_fallback(na_aligner_t *al, bool enable);
void na_aligner_set_wfa_budget(na_aligner_t *al, size_t budget);
void na_aligner_set_min_local_score(na_aligner_t *al, na_score_t min_score);
void na_aligner_set_fast_local(na_aligner_t *al, bool enable);
/* The status and detail of the aligner's last call; valid until its next */
const na_error_t *na_aligner_last_error(const na_aligner_t *al);

/* Global alignment (Needleman-Wunsch / Hirschberg / BiWFA / Myers).
 * Sequences are byte strings of the lengths given (empty, and NUL bytes,
 * are fine). On an error out is cleared (length 0, score 0). Both calls
 * end any local iteration in progress. */
na_status_t na_global_align(na_aligner_t *al,
                            const char *a, size_t len_a,
                            const char *b, size_t len_b,
                            const na_scoring_t *sc,
                            na_alignment_t *out);

/* The score alone, in linear memory (Myers' bit-vector algorithm for
 * unit costs); flags linear memory refuses are scored with the compact
 * matrices instead */
na_status_t na_global_score(na_aligner_t *al,
                            const char *a, size_t len_a,
                            const char *b, size_t len_b,
                            const na_scoring_t *sc,
                            na_score_t *out_score);

/* Local alignment (Smith-Waterman): na_local_next() then yields the hits
 * best first (NA_OK with a hit, NA_DONE after the last), no two sharing a
 * cell, down to min_local_score. a, b and sc must stay valid and unchanged
 * until the next na_local_align(), na_global_align(), na_global_score() or
 * na_aligner_free(), each of which ends the iteration. The first hit of a
 * long pair needs little memory; the second fills the score matrices (12
 * bytes per cell). fast_local reconstructs the first hit with wavefront or
 * linear-memory alignment, verified against the exact score: the same
 * score, perhaps another equal-scoring alignment, and as the cells it
 * marks as used differ, later hits may differ too. */
na_status_t na_local_align(na_aligner_t *al,
                           const char *a, size_t len_a,
                           const char *b, size_t len_b,
                           const na_scoring_t *sc);

na_status_t na_local_next(na_aligner_t *al, na_alignment_t *out);

/* Utilities */
/* Score a pair of aligned strings as the aligner would, free start/end
 * gap flags included. A column of two gaps is NA_ERR_INVALID_ARG; a pair
 * without a score, or a mismatch under NA_NO_MISMATCHES, NA_ERR_SCORING_PAIR */
na_status_t na_score_alignment(const char *res_a, const char *res_b, size_t len,
                               const na_scoring_t *sc, na_score_t *out_score);

/* The three score matrices of the last NA_METHOD_FULL alignment (nothing
 * after any other method). Reads the sequences and scoring given to that
 * call, which must still be valid. */
void na_debug_fprint_matrices(const na_aligner_t *al, FILE *fp);

#ifdef __cplusplus
}
#endif

#endif /* NUCLEICS_ALIGN_V2_HEADER_SEEN */
