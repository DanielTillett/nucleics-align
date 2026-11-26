/*
 alignment_scoring.h
 url: https://github.com/DanielTillett/nucleics-align
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

// Scores (match, mismatch, gap open and extend, substitution and wildcard
// scores) should lie within -SCORE_ABS_MAX .. SCORE_ABS_MAX, so a few added
// together cannot overflow. The command line tools and the scoring file
// loaders refuse others
#define SCORE_ABS_MAX (INT_MAX / 4)

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

// Recompute min_penalty and max_penalty from the current scores (after
// changing match, mismatch or the gap penalties in place; scoring_add_*()
// keep them up to date themselves)
void scoring_update_penalty_range(scoring_t* scoring);

#define scoring_is_wildcard(scoring,c) (get_wildcard_bit(scoring,c))

void scoring_add_wildcard(scoring_t* scoring, char c, int s);

void scoring_add_mutation(scoring_t* scoring, char a, char b, int score);

// Compare characters case sensitively from now on. Scores already set (in
// lower case, as scoring_add_*() store them while case insensitive) are
// copied to every mix of cases, so 'A' against 'R' or 'r' scores as 'a'
// against 'r' did. 'A' against 'a' gets no score from them (a mismatch
// where match/mismatch scoring is used)
void scoring_make_case_sensitive(scoring_t* scoring);


// The score of a against b, and whether they count as a match. Exits on a
// pair the scoring has no score for; scoring_try_lookup() returns false
// instead, setting nothing
void scoring_lookup(const scoring_t* scoring, char a, char b,
                    int *score, bool *is_match);
bool scoring_try_lookup(const scoring_t* scoring, char a, char b,
                        int *score, bool *is_match);

// Does scoring_lookup() find a score for (a, b)? (It exits if not)
bool scoring_has_pair(const scoring_t* scoring, char a, char b);

// Some scoring systems
void scoring_system_PAM30(scoring_t *scoring);
void scoring_system_PAM70(scoring_t *scoring);
void scoring_system_BLOSUM80(scoring_t *scoring);
void scoring_system_BLOSUM62(scoring_t *scoring);
void scoring_system_DNA_hybridization(scoring_t *scoring);
void scoring_system_default(scoring_t *scoring); // DNA/RNA default

// Set one of the systems above by name, ignoring case (PAM30, PAM70,
// BLOSUM80, BLOSUM62, DNA_HYBRIDIZATION); false if there is none
bool scoring_system_by_name(scoring_t *scoring, const char *name);
// Their names, separated by '|'
const char* scoring_system_names(void);

// The scoring file formats in scoring/, read from memory (see
// alignment_scoring_parse.c): a substitution matrix, or pairs. Nothing
// changes on an error; *line (may be NULL) then gets the 1-based line it
// is on, or -1 for one of the file as a whole
typedef enum {
  SCORING_PARSE_OK = 0,
  SCORING_PARSE_NOMEM,
  SCORING_PARSE_EMPTY,          // matrix: no header line
  SCORING_PARSE_FEW_HEADINGS,   // matrix: a header of one character
  SCORING_PARSE_SEPARATOR,      // matrix: a digit or - as the separator
  SCORING_PARSE_SEP_MISSING,    // matrix: a column without its separator
  SCORING_PARSE_EMPTY_HEADING,  // matrix: a separator with nothing after it
  SCORING_PARSE_FEW_COLUMNS,    // matrix: a row ends before its last column
  SCORING_PARSE_EXPECTED_WS,    // matrix: a character where whitespace
                                // should separate scores
  SCORING_PARSE_NUMBER,         // matrix: a score that is not a number
  SCORING_PARSE_RANGE,          // a score beyond +-SCORE_ABS_MAX
  SCORING_PARSE_MANY_COLUMNS,   // matrix: more scores than headings
  SCORING_PARSE_NO_ROWS,        // matrix: headings only
  SCORING_PARSE_SHORT_LINE,     // pairs: too short for a pair and a score
  SCORING_PARSE_SEPARATORS,     // pairs: the two separators differ
  SCORING_PARSE_INVALID_NUMBER, // pairs: a score that is not an int
  SCORING_PARSE_NO_PAIRS        // pairs: nothing added
} scoring_parse_t;

scoring_parse_t scoring_parse_matrix(scoring_t *scoring, const char *text,
                                     size_t len, int *line);
scoring_parse_t scoring_parse_pairs(scoring_t *scoring, const char *text,
                                    size_t len, int *line);
// A short message for each status (the command line tools' wording)
const char *scoring_parse_strerror(scoring_parse_t status);

#ifdef __cplusplus
}
#endif

#endif
