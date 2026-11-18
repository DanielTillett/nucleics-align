/*
 alignment_scoring.h
 url: https://github.com/DanielTillett/seq-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */

#ifndef ALIGNMENT_SCORING_HEADER_SEEN
#define ALIGNMENT_SCORING_HEADER_SEEN

#include <inttypes.h>
#include <stdbool.h>
#include <limits.h> // INT_MIN

typedef int score_t;
#define SCORE_MIN INT_MIN

typedef struct
{
  int gap_open, gap_extend;

  // Needleman Wunsch only
  // Turn these on to turn off penalties for gaps at the start/end of alignment
  bool no_start_gap_penalty, no_end_gap_penalty;

  // Turn at most one of these on at a time to prevent gaps/mismatches
  bool no_gaps_in_a, no_gaps_in_b, no_mismatches;

  // Apply no_gaps_in_a or no_gaps_in_b to whichever sequence is longer
  // (the first sequence if they are the same length)
  bool no_gaps_in_longest;

  // If swap_score not set, should we use match/mismatch values?
  bool use_match_mismatch;
  int match, mismatch;

  bool case_sensitive;

  // Array of characters that match to everything with the same penalty (i.e. 'N's)
  uint32_t wildcards[256/32], swap_set[256][256/32];
  score_t wildscores[256], swap_scores[256][256];
  int min_penalty, max_penalty; // min, max {match/mismatch,gapopen etc.}
} scoring_t;

#ifndef bitset32_get
  #define bitset32_get(arr,idx)   (((arr)[(idx)>>5] >> ((idx)&31)) & 0x1U)
  #define bitset32_set(arr,idx)   ((arr)[(idx)>>5] |=   (1U<<((idx)&31)))
  #define bitset32_clear(arr,idx) ((arr)[(idx)>>5] &=  ~(1U<<((idx)&31)))
#endif

// Characters are cast to unsigned char so bytes >= 0x80 index tables safely
#define get_wildcard_bit(scoring,c) bitset32_get((scoring)->wildcards,(unsigned char)(c))
#define set_wildcard_bit(scoring,c) bitset32_set((scoring)->wildcards,(unsigned char)(c))

#define get_swap_bit(scoring,a,b)   bitset32_get((scoring)->swap_set[(unsigned char)(a)],(unsigned char)(b))
#define set_swap_bit(scoring,a,b)   bitset32_set((scoring)->swap_set[(unsigned char)(a)],(unsigned char)(b))

#ifdef __cplusplus
extern "C" {
#endif

void scoring_init(scoring_t* scoring, int match, int mismatch,
                  int gap_open, int gap_extend,
                  bool no_start_gap_penalty, bool no_end_gap_penalty,
                  bool no_gaps_in_a, bool no_gaps_in_b,
                  bool no_mismatches, bool case_sensitive);

#define scoring_is_wildcard(scoring,c) (get_wildcard_bit(scoring,c))

void scoring_add_wildcard(scoring_t* scoring, char c, int s);

void scoring_add_mutation(scoring_t* scoring, char a, char b, int score);

void scoring_print(const scoring_t* scoring);

void scoring_lookup(const scoring_t* scoring, char a, char b,
                    int *score, bool *is_match);

// Some scoring systems
void scoring_system_PAM30(scoring_t *scoring);
void scoring_system_PAM70(scoring_t *scoring);
void scoring_system_BLOSUM80(scoring_t *scoring);
void scoring_system_BLOSUM62(scoring_t *scoring);
void scoring_system_DNA_hybridization(scoring_t *scoring);
void scoring_system_default(scoring_t *scoring); // DNA/RNA default

#ifdef __cplusplus
}
#endif

#endif
