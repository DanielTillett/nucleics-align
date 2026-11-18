/*
 alignment.h
 original author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/DanielTillett/seq-align
 license: BSD 3-Clause, see LICENSE
 May 2013
 */

#ifndef ALIGNMENT_HEADER_SEEN
#define ALIGNMENT_HEADER_SEEN

#include <string.h> // memset
#include "alignment_scoring.h"

#ifndef ROUNDUP2POW
  #define ROUNDUP2POW(x) _rndup2pow64(x)
  static inline size_t _rndup2pow64(unsigned long long x) {
    // long long >=64 bits guaranteed in C99
    if(x <= 1) return 1; // without this, 0 would round up to 0
    --x; x|=x>>1; x|=x>>2; x|=x>>4; x|=x>>8; x|=x>>16; x|=x>>32; ++x;
    return x;
  }
#endif

typedef struct
{
  const scoring_t* scoring;
  const char *seq_a, *seq_b;
  size_t score_width, score_height; // width=len(seq_a)+1, height=len(seq_b)+1
  score_t *match_scores, *gap_a_scores, *gap_b_scores;
  size_t capacity;
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

void sub_table_build(sub_table_t *t, const scoring_t *scoring,
                     const char *a, size_t len_a,
                     const char *b, size_t len_b);
void sub_table_free(sub_table_t *t);

// Use the striped SIMD fill where it applies (default true). It gives the
// same matrices as the scalar fill; turning it off is for tests/benchmarks
extern bool alignment_striped_enabled;

// Large striped fills split the columns into blocks, one thread each: at
// least alignment_block_min_cells cells and alignment_block_min_cols columns
// per block, on up to alignment_max_threads threads (0: half the CPUs)
extern size_t alignment_block_min_cells, alignment_block_min_cols;
extern int alignment_max_threads;

// Constructors/Destructors for alignment
alignment_t* alignment_create(size_t capacity);
void alignment_ensure_capacity(alignment_t* result, size_t strlength);
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

#ifdef __cplusplus
}
#endif

#endif
