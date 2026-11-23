/*
 alignment.h
 original author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 May 2013
 */

#ifndef ALIGNMENT_HEADER_SEEN
#define ALIGNMENT_HEADER_SEEN

#include <string.h> // memset
#include <stddef.h> // ptrdiff_t
#include <stdio.h> // FILE
#include <stdint.h>
#include <stdbool.h>
#include "alignment_scoring.h"

#ifndef ROUNDUP2POW
  #define ROUNDUP2POW(x) _rndup2pow64(x)
  static inline size_t _rndup2pow64(unsigned long long x) {
    // long long >=64 bits guaranteed in C99
    if(x <= 1) return 1; // without this, 0 would round up to 0
    // No power of two that big fits: leave it (an allocation then fails)
    if(x > ((unsigned long long)SIZE_MAX >> 1) + 1)
      return x < SIZE_MAX ? (size_t)x : SIZE_MAX;
    --x; x|=x>>1; x|=x>>2; x|=x>>4; x|=x>>8; x|=x>>16; x|=x>>32; ++x;
    return x;
  }
#endif

// Settings for an alignment: threads, and switches and sizes for tuning and
// tests. None changes a result. alignment_options_default() gives the
// process-wide defaults (the alignment_*, hirschberg_* and
// smith_waterman_*_min_cells variables), which apply to any aligner
// without options of its own; set options to use others, e.g. different
// ones on different threads
typedef struct
{
  int max_threads; // full-matrix fills (0: half the CPUs)
  bool striped, band, narrow, byte, diff;
  // Needleman-Wunsch too big for the matrices (too long, or out of memory):
  // align in linear space instead (see alignment_linear_fallback); and
  // pairs whose longer sequence has at least linear_above letters (0:
  // none), see alignment_linear_above. The settings that can change a
  // result
  bool linear_fallback;
  size_t linear_above;
  // Needleman-Wunsch by the wavefront algorithm where it applies (see
  // alignment_wfa_alignment): the same score, perhaps another alignment
  bool wfa_align;
  size_t block_min_cells, block_min_cols, band_strip_rows, band_min_width;
  size_t band_keep_max_cells; // see alignment_band_keep_max_cells
  size_t band_wfa_budget; // see alignment_band_wfa_budget
  size_t sw_rect_min_cells; // see smith_waterman_rect_min_cells
  size_t sw_strip_min_cells; // see smith_waterman_strip_min_cells
  // Linear space (hirschberg.h)
  int hb_max_threads; // 0: half the CPUs
  bool hb_band;
  size_t hb_block_min_cells, hb_block_min_cols;
  size_t hb_band_strip_rows, hb_band_min_width;
  size_t hb_reuse_depth, hb_reuse_min_rows; // see hirschberg_reuse_depth
} alignment_options_t;

typedef struct
{
  const scoring_t* scoring;
  const char *seq_a, *seq_b;
  size_t score_width, score_height; // width=len(seq_a)+1, height=len(seq_b)+1
  // Settings for this aligner's alignments (NULL: the defaults). opts holds
  // those in effect during an alignment
  const alignment_options_t *options;
  alignment_options_t opts;
  score_t *match_scores, *gap_a_scores, *gap_b_scores;
  size_t capacity;
  // Set compact before aligning to keep one traceback byte per cell instead
  // of the three score matrices, where the options allow (12x less memory).
  // traced then says the bytes were made and the matrices were not filled
  bool compact, traced;
  uint8_t *trace; // per cell: for matrix m, bits 2m..2m+1 = matrix stepped to
  size_t trace_capacity;
  score_t final_scores[3]; // traced: MATCH, GAP_A, GAP_B at the last cell
  size_t best_index; // traced Smith-Waterman: first cell with best_score
  score_t best_score;
  // Banded trace (Needleman-Wunsch): only cells near the main diagonals
  // were filled, and cell (x, y)'s byte is at trace[trace_rows[y] + x]
  bool trace_banded;
  ptrdiff_t *trace_rows;
  size_t trace_rows_capacity;
  // A band too big to keep (see alignment_band_keep_max_cells): trace and
  // trace_rows hold one strip of rows at a time, which alignment_trace_load()
  // fills again as needed (private; NULL otherwise)
  void *trace_strips;
} aligner_t;

// Store alignment result here
typedef struct
{
  char *result_a, *result_b;
  size_t capacity, length;
  size_t pos_a, pos_b; // position of first base (0-based)
  size_t len_a, len_b; // number of bases in alignment
  score_t score;
} alignment_t;

// Substitution scores for every pair of characters present in two
// sequences, so DP inner loops do an array load instead of scoring_lookup()
typedef struct
{
  unsigned char code[256]; // byte -> row/column in scores
  size_t n; // number of distinct bytes
  int *scores; // scores[code_b*n + code_a]
  bool *is_match;
} sub_table_t;

// Matrix names
enum Matrix { MATCH,GAP_A,GAP_B };

#define MATRIX_NAME(x) ((x) == MATCH ? "MATCH" : ((x) == GAP_A ? "GAP_A" : "GAP_B"))

// Ties: every traceback (alignment_reverse_move(), the traceback bytes of
// every fill, and the walks that read them) steps to the first predecessor,
// in the order GAP_A, GAP_B, MATCH, whose score plus the step's gives the
// current score; Needleman-Wunsch starts in the last cell's best matrix in
// the same order. Any new traceback must keep this to give the same
// alignments (see the README's notes on equal-scoring alignments)

// Traceback byte value for a step with no valid predecessor
#define ALN_TRACE_NONE 3

// Matrix the step back from matrix m at cell index goes to (traced only)
static inline int alignment_trace_next(const aligner_t *aligner, size_t index,
                                       enum Matrix m)
{
  return (aligner->trace[index] >> (2 * (int)m)) & 3;
}

// Index of the traceback byte of cell (x, y)
static inline size_t alignment_trace_index(const aligner_t *aligner,
                                           size_t x, size_t y)
{
  return aligner->trace_banded
         ? (size_t)(aligner->trace_rows[y] + (ptrdiff_t)x)
         : y * aligner->score_width + x;
}

#ifdef __cplusplus
extern "C" {
#endif

// Printing colour codes
extern const char align_col_mismatch[], align_col_indel[], align_col_context[],
                  align_col_stop[];

#define aligner_init(a) (memset(a, 0, sizeof(aligner_t)))
void aligner_align(aligner_t *aligner,
                   const char *seq_a, const char *seq_b,
                   size_t len_a, size_t len_b,
                   const scoring_t *scoring, char is_sw);
void aligner_destroy(aligner_t *aligner);

// Errors that would otherwise end the program. The *_checked() aligners
// (needleman_wunsch.h, smith_waterman.h, hirschberg.h) report them instead:
// errors in the input before any work, running out of memory as it happens
typedef enum {
  ALN_OK = 0,
  ALN_ERR_SCORING,     // a character pair has no score: err->a, err->b
  ALN_ERR_TOO_LONG,    // the score matrices would not fit in memory addresses
  ALN_ERR_OVERFLOW,    // scores could pass INT_MAX (see SCORE_ABS_MAX)
  ALN_ERR_UNSUPPORTED, // linear space refuses these options or penalties
  ALN_ERR_NOMEM        // out of memory; nothing to fetch until aligned again
} aln_status_t;

typedef struct {
  aln_status_t status;
  char a, b; // ALN_ERR_SCORING: the pair, as scoring_lookup() reports it
} aln_error_t;

// Whether every pair of a character of a with one of b has a score. If not,
// err (may be NULL) gets the first pair scoring_lookup() would stop at
bool alignment_scoring_complete(const scoring_t *scoring,
                                const char *a, size_t len_a,
                                const char *b, size_t len_b,
                                aln_error_t *err);

// Whether a score of a full-matrix alignment of a and b (is_sw: local) could
// pass the range of int, given the scores of the pairs in the sequences
bool alignment_scores_may_overflow(const scoring_t *scoring,
                                   const char *a, size_t len_a,
                                   const char *b, size_t len_b, bool is_sw);

// The first error a full-matrix alignment of a and b would stop at: the
// scoring first (as Smith-Waterman checks it) or last (Needleman-Wunsch).
// err may be NULL
aln_status_t alignment_check(const scoring_t *scoring,
                             const char *a, size_t len_a,
                             const char *b, size_t len_b,
                             bool scoring_first, aln_error_t *err);

// The message for an error, as the aligners print it before exiting
void alignment_fprint_error(FILE *out, const aln_error_t *err);

// As aligner_align(), but returns ALN_ERR_TOO_LONG, ALN_ERR_OVERFLOW or
// ALN_ERR_NOMEM (with no matrices or trace to read) rather than exiting
aln_status_t aligner_align_checked(aligner_t *aligner,
                                   const char *seq_a, const char *seq_b,
                                   size_t len_a, size_t len_b,
                                   const scoring_t *scoring, char is_sw);

// Returns false if out of memory, with nothing to free
bool sub_table_build(sub_table_t *t, const scoring_t *scoring,
                     const char *a, size_t len_a,
                     const char *b, size_t len_b);
void sub_table_free(sub_table_t *t);

// The best global alignment's score by the wavefront algorithm (see
// alignment_wfa.c), for a table built for these sequences. Returns false,
// setting nothing, unless both are non-empty and the scoring has one match
// and one mismatch score, gap penalties and no other rules; or past budget
// diagonal updates; or if out of memory
bool alignment_wfa_score(const sub_table_t *table, const scoring_t *scoring,
                         const char *seq_a, size_t len_a,
                         const char *seq_b, size_t len_b,
                         size_t budget, score_t *score);

// The best global alignment itself by the wavefront algorithm, keeping
// every wavefront (at most budget cells), for a table built for these
// sequences: the same score as needleman_wunsch_align(), but where several
// alignments share it, perhaps another of them. ALN_ERR_UNSUPPORTED if the
// scoring is not one match and one mismatch score with gap penalties, or
// past the budget; ALN_ERR_NOMEM if out of memory
aln_status_t alignment_wfa_align(const sub_table_t *table,
                                 const scoring_t *scoring,
                                 const char *seq_a, size_t len_a,
                                 const char *seq_b, size_t len_b,
                                 size_t budget, alignment_t *result);

// The best global alignment's score under unit costs (match 0, mismatch and
// gap letter -1, gap open 0: minus the edit distance) by Myers' bit-vector
// algorithm (see alignment_myers.c), for a table built for these
// sequences. Returns false, setting nothing, unless the scoring is unit
// cost for every pair in the table and seq_a is not empty; or if out of
// memory
bool alignment_myers_score(const sub_table_t *table, const scoring_t *scoring,
                           const char *seq_a, size_t len_a,
                           const char *seq_b, size_t len_b, score_t *score);

// Full bottom-row global scores under unit costs for Hirschberg --linear:
// rev=false fills cc[0..len_b] for prefixes of seq_b; rev=true fills
// cc[len_b..0] for suffixes of seq_b against reversed seq_a
bool alignment_myers_row(const sub_table_t *table, const scoring_t *scoring,
                         const char *seq_a, size_t len_a,
                         const char *seq_b, size_t len_b,
                         bool rev, score_t *cc);

// Use the striped SIMD fill where it applies (default true). It gives the
// same matrices as the scalar fill; turning it off is for tests/benchmarks
extern bool alignment_striped_enabled;

// Large striped fills split the columns into blocks, one thread each: at
// least alignment_block_min_cells cells and alignment_block_min_cols columns
// per block, on up to alignment_max_threads threads (0: half the CPUs)
extern size_t alignment_block_min_cells, alignment_block_min_cols;
extern int alignment_max_threads;

// Compact Needleman-Wunsch fills skip cells that cannot be on a best path,
// with the same result (see alignment_band_width()). alignment_band_enabled
// (default true) turns that off; alignment_band_strip_rows (0: automatic)
// and alignment_band_min_width (first band tried, default 64) are for
// tuning and tests
extern bool alignment_band_enabled;
extern size_t alignment_band_strip_rows, alignment_band_min_width;
// Bands of more than this many cells (0, the default: a quarter of the
// physical memory, one traceback byte per cell) keep only one strip of
// rows' bytes at a time, and each strip's top row: the traceback fills a
// strip again when it reaches it (see alignment_trace_load()). Same
// alignment, in at most about sqrt(12 x rows) rows' worth of bytes (1 Mb at
// 1% divergence: 35 GB -> 350 MB), but the band is filled about twice. For
// tuning and tests
extern size_t alignment_band_keep_max_cells;
// Before a band, the best score is worked out by the wavefront algorithm
// (alignment_wfa_score()) where the scoring allows, in up to this many
// diagonal updates per letter of a (default 16; 0: none).
// If it finishes, the score gives the band that passes on its first fill,
// or shows that no band pays; if not, bands start narrow and widen as
// before. Same alignment either way. For tuning and tests
extern size_t alignment_band_wfa_budget;

// Smith-Waterman fills and scans use 8-bit lanes (4 times as many) or
// 16-bit lanes (twice as many) where the scores allow, and start again
// with wider lanes if they grow too big (same result): 8, then 16, then 32
// bits. The linear-space passes over every column (hirschberg.h) use
// 16-bit lanes where the scores allow, which they check before starting.
// alignment_byte_enabled and alignment_narrow_enabled (default true) turn
// off the 8- and 16-bit lanes, for tests
extern bool alignment_narrow_enabled, alignment_byte_enabled;

// Needleman-Wunsch traceback-byte fills that are not banded keep each cell
// as its difference from the cell above, in 8- or 16-bit lanes where the
// scores allow (same bytes, so the same alignment). alignment_diff_enabled
// (default true) turns that off, for tests
extern bool alignment_diff_enabled;

// Needleman-Wunsch alignments too big for the score matrices or traceback
// bytes (ALN_ERR_TOO_LONG, or ALN_ERR_NOMEM) are made in linear space
// instead (hirschberg_align_opts()) where it supports the scoring: the same
// optimal score, but it may choose another of several equal-scoring
// alignments. Default false
extern bool alignment_linear_fallback;
// Needleman-Wunsch alignments whose longer sequence has at least this many
// letters are made in linear space (hirschberg_align_opts()) where it
// supports the scoring, and with the matrices otherwise: faster and in far
// less memory for long pairs, with the same score, but it may choose
// another of several equal-scoring alignments. Default 0 (never)
extern size_t alignment_linear_above;
// Needleman-Wunsch alignments are made by the wavefront algorithm
// (alignment_wfa_align()) where the scoring is one match and one mismatch
// score with gap penalties and the pair is similar enough to finish within
// 64 cells per letter (at most 2^25 cells, about 400 MB), and with the
// matrices otherwise: much faster for near-identical pairs, with the same
// score, but it may choose another of several equal-scoring alignments.
// Default false
extern bool alignment_wfa_alignment;
// Diagonal-update cap for --wfa / alignment_wfa_alignment (0: automatic,
// min(2^28, 256 * (len_a + len_b))). Separate from alignment_band_wfa_budget
extern size_t alignment_wfa_budget;

// The default settings, from the variables above (and their hirschberg_*
// and smith_waterman_* counterparts) as they are now
void alignment_options_default(alignment_options_t *opts);

// Half-width w of the band a global alignment of m x n needs for a path
// scoring at least `lower`: such a path stays on diagonals from
// MIN2(0, n-m) - w to MAX2(0, n-m) + w. smax is the best substitution score
// (at least 0); gap penalties must be 0 or less. SIZE_MAX: no useful band
size_t alignment_band_width(score_t smax, score_t gap_extend,
                            size_t m, size_t n, score_t lower);

// Smith-Waterman scans in O(len_b) memory (see alignment.c), used by
// smith_waterman.c for the best hit of large alignments
bool alignment_sw_scan(const char *a, size_t len_a, const char *b, size_t len_b,
                       const scoring_t *scoring, score_t target,
                       const alignment_options_t *opts,
                       score_t *best, size_t *x, size_t *y);
// The scan for the best cell (target 0), but giving up (returning false)
// once the best so far reaches half the best substitution score per row
// filled, from a sixteenth of the rows on: a hit that rises that fast from
// near the start spans most of the sequences, so finding where it starts
// would not pay
bool alignment_sw_scan_short(const char *a, size_t len_a,
                             const char *b, size_t len_b,
                             const scoring_t *scoring,
                             const alignment_options_t *opts,
                             score_t *best, size_t *x, size_t *y);
// The scan for the best cell (target 0), also keeping every every-th row
// of the three matrices: row k*every (k = 0, 1, ...) at ckpt + k*3*(len_a +
// 1) as MATCH, GAP_A then GAP_B scores by column (0 .. len_a). ckpt needs
// room for len_b/every + 1 rows. For alignment_sw_trace_rows()
bool alignment_sw_rows(const char *a, size_t len_a,
                       const char *b, size_t len_b,
                       const scoring_t *scoring,
                       const alignment_options_t *opts,
                       size_t every, score_t *ckpt,
                       score_t *best, size_t *x, size_t *y);
// Traceback bytes of rows y0 .. y1 and columns 1 .. xmax of a compact
// Smith-Waterman fill, from row y0-1 (row: as alignment_sw_rows() keeps
// it): cell (x, y) at trace[(y - y0) * xmax + x - 1], the same byte as a
// fill of the whole matrix gives it (cells right of xmax cannot change
// those up to it). Returns false if out of memory or where the striped
// fill does not apply
bool alignment_sw_trace_rows(const char *a, size_t len_a,
                             const char *b, size_t len_b,
                             const scoring_t *scoring,
                             const sub_table_t *table,
                             const alignment_options_t *opts,
                             score_t max_score,
                             const score_t *row, size_t y0, size_t y1,
                             size_t xmax, uint8_t *trace);

// Before reading the traceback byte of cell (x, y) of a band kept in strips
// (aligner->trace_strips set), fill again the strip holding row y, up to
// column x, from the row kept above it: the same bytes as the whole band
// had. Does nothing otherwise, or if that strip is loaded. Returns false if
// out of memory
bool alignment_trace_load(aligner_t *aligner, size_t x, size_t y);

// Constructors/Destructors for alignment. alignment_create() and
// alignment_ensure_capacity() exit if out of memory; alignment_reserve()
// returns false, leaving the result as it was
alignment_t* alignment_create(size_t capacity);
void alignment_ensure_capacity(alignment_t* result, size_t strlength);
bool alignment_reserve(alignment_t* result, size_t strlength);
void alignment_free(alignment_t* result);

void alignment_reverse_move(enum Matrix *curr_matrix, score_t *curr_score,
                            size_t *score_x, size_t *score_y,
                            size_t *arr_index, const aligner_t *aligner);

// Printing
void alignment_print_matrices(const aligner_t *aligner);

void alignment_colour_print_against(const char *alignment_a,
                                    const char *alignment_b,
                                    char case_sensitive);

void alignment_print_spacer(const char* alignment_a, const char* alignment_b,
                            const scoring_t* scoring);

// The same, to out instead of stdout
void alignment_colour_fprint_against(FILE *out, const char *alignment_a,
                                     const char *alignment_b,
                                     char case_sensitive);
void alignment_fprint_spacer(FILE *out, const char* alignment_a,
                             const char* alignment_b,
                             const scoring_t* scoring);

#ifdef __cplusplus
}
#endif

#endif
