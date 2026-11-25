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
#define NA_VERSION_STR   "2.0.0"

/* Legacy version macros for CLI compatibility */
#define NUCLEICS_ALIGN_VERSION_STR NA_VERSION_STR
#define NUCLEICS_ALIGN_VERSION     0x200

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
 * Scoring configuration (opaque; safe to share read-only across threads)
 * Gap penalty formula: a gap of length k >= 1 scores gap_open + k * gap_extend.
 * ========================================================================= */
typedef struct na_scoring na_scoring_t;

enum {
  NA_FREE_START_GAP  = 1u << 0,
  NA_FREE_END_GAP    = 1u << 1,
  NA_NO_GAPS_A       = 1u << 2,
  NA_NO_GAPS_B       = 1u << 3,
  NA_NO_GAPS_LONGEST = 1u << 4,
  NA_NO_MISMATCHES   = 1u << 5,
  NA_CASE_SENSITIVE  = 1u << 6
};

na_scoring_t *na_scoring_new(int match, int mismatch, int gap_open, int gap_extend);
na_scoring_t *na_scoring_new_preset(const char *name);
void          na_scoring_free(na_scoring_t *sc);
const char   *na_scoring_preset_names(void);

na_status_t na_scoring_set_flags(na_scoring_t *sc, unsigned flags);
na_status_t na_scoring_set_match_mismatch(na_scoring_t *sc, bool enable,
                                          int match, int mismatch);
na_status_t na_scoring_set_gaps(na_scoring_t *sc, int gap_open, int gap_extend);
na_status_t na_scoring_set_pair(na_scoring_t *sc, char a, char b, int score);
na_status_t na_scoring_set_wildcard(na_scoring_t *sc, char c, int score);
na_status_t na_scoring_parse_matrix(na_scoring_t *sc, const char *text,
                                    size_t len, na_error_t *err);
na_status_t na_scoring_parse_pairs(na_scoring_t *sc, const char *text,
                                   size_t len, na_error_t *err);
na_status_t na_scoring_lookup(const na_scoring_t *sc, char a, char b,
                              int *score_out, bool *is_match_out);

/* =========================================================================
 * Alignment result buffer (opaque; reusable across calls)
 * ========================================================================= */
typedef struct na_alignment na_alignment_t;

na_alignment_t *na_alignment_new(void);
void            na_alignment_free(na_alignment_t *res);

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

typedef enum {
  NA_METHOD_AUTO = 0,    /* Exact compact/banded NW (default tie-breaking) */
  NA_METHOD_FULL,        /* Full 3-matrix DP (required for debug matrix dump) */
  NA_METHOD_LINEAR,      /* O(n)-memory Hirschberg / Myers-Miller */
  NA_METHOD_WFA          /* Bidirectional Wavefront Alignment (falls back if needed) */
} na_method_t;

na_aligner_t *na_aligner_new(void);
void          na_aligner_free(na_aligner_t *al);

void na_aligner_set_threads(na_aligner_t *al, unsigned threads);
void na_aligner_set_method(na_aligner_t *al, na_method_t method);
void na_aligner_set_linear_above(na_aligner_t *al, size_t max_len);
void na_aligner_set_linear_fallback(na_aligner_t *al, bool enable);
void na_aligner_set_wfa_budget(na_aligner_t *al, size_t budget);
void na_aligner_set_min_local_score(na_aligner_t *al, na_score_t min_score);
void na_aligner_set_fast_local(na_aligner_t *al, bool enable);
const na_error_t *na_aligner_last_error(const na_aligner_t *al);

/* Global alignment (Needleman-Wunsch / Hirschberg / BiWFA / Myers) */
na_status_t na_global_align(na_aligner_t *al,
                            const char *a, size_t len_a,
                            const char *b, size_t len_b,
                            const na_scoring_t *sc,
                            na_alignment_t *out);

na_status_t na_global_score(na_aligner_t *al,
                            const char *a, size_t len_a,
                            const char *b, size_t len_b,
                            const na_scoring_t *sc,
                            na_score_t *out_score);

/* Local alignment (Smith-Waterman): hits are yielded best-first.
 * Sequence buffers (a, b) and scoring (sc) must remain valid until the next
 * na_local_align() / na_global_align() call or na_aligner_free(). */
na_status_t na_local_align(na_aligner_t *al,
                           const char *a, size_t len_a,
                           const char *b, size_t len_b,
                           const na_scoring_t *sc);

na_status_t na_local_next(na_aligner_t *al, na_alignment_t *out);

/* Utilities */
na_status_t na_score_alignment(const char *res_a, const char *res_b, size_t len,
                               const na_scoring_t *sc, na_score_t *out_score);

void na_debug_fprint_matrices(const na_aligner_t *al, FILE *fp);

#ifdef __cplusplus
}
#endif

#endif /* NUCLEICS_ALIGN_V2_HEADER_SEEN */
