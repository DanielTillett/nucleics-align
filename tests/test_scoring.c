/*
 tests/test_scoring.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 Scoring tables: penalty ranges, case, the built-in systems
 */

#include "test_common.h"


// min_penalty / max_penalty follow the scores as they are changed: a harsh
// mismatch replaced by a mild one no longer sets the DP floor
void scoring_test_penalty_range()
{
  scoring_t scoring;
  scoring_init(&scoring, 1, -500000000, -4, -1, false, false, false, false,
               false, false);
  ASSERT(scoring.min_penalty == -500000000);
  scoring.mismatch = -2;
  scoring_update_penalty_range(&scoring);
  ASSERT(scoring.min_penalty == -5 && scoring.max_penalty == 1);
  // Pairs and wildcards count, and only while match/mismatch are in use
  scoring_add_mutation(&scoring, 'a', 'c', -9);
  scoring_add_wildcard(&scoring, 'n', 7);
  scoring_update_penalty_range(&scoring);
  ASSERT(scoring.min_penalty == -9 && scoring.max_penalty == 7);
  scoring.use_match_mismatch = false;
  scoring.match = 100;
  scoring_update_penalty_range(&scoring);
  ASSERT(scoring.min_penalty == -9 && scoring.max_penalty == 7);
}

/* scoring_make_case_sensitive() keeps scores already set, for upper-case
 * characters too, and tells a character from its other case */
void scoring_test_case_sensitive()
{
  scoring_t scoring;
  int s1, s2;
  bool m1, m2;

  scoring_system_BLOSUM62(&scoring);
  scoring_lookup(&scoring, 'w', 'r', &s1, &m1);
  scoring_make_case_sensitive(&scoring);
  scoring_lookup(&scoring, 'W', 'R', &s2, &m2);
  ASSERT(s1 == s2 && m1 == m2);
  scoring_lookup(&scoring, 'w', 'r', &s2, &m2);
  ASSERT(s1 == s2 && m1 == m2);
  scoring_lookup(&scoring, 'W', 'r', &s2, &m2);
  ASSERT(s1 == s2 && m1 == m2);
  scoring_lookup(&scoring, 'w', 'R', &s2, &m2);
  ASSERT(s1 == s2 && m1 == m2);
  scoring_lookup(&scoring, 'W', 'w', &s2, &m2); // match/mismatch
  ASSERT(s2 == scoring.mismatch && !m2);

  scoring_system_default(&scoring);
  scoring_add_wildcard(&scoring, 'n', 0);
  scoring_make_case_sensitive(&scoring);
  scoring_lookup(&scoring, 'A', 'a', &s1, &m1);
  ASSERT(s1 == scoring.mismatch && !m1);
  scoring_lookup(&scoring, 'N', 'C', &s1, &m1);
  ASSERT(s1 == 0 && m1);
}

/* The two matrix forms read through the shared parser (as the tools'
 * --substitution_matrix reads them), case-folded: whitespace separated as
 * in the files in scoring/, and a single-character separator */
void scoring_test_load_matrix()
{
  scoring_t scoring;
  int score, line;
  bool is_match;
  const char *ws = "# comment\n   A  C\nA  5 -3\nC -3  7\n";
  const char *comma = ",A,C\nA,5,-3\nC,-3,7\n";

  scoring_init(&scoring, 0, 0, -4, -1, false, false, false, false, false, false);
  ASSERT(scoring_parse_matrix(&scoring, ws, strlen(ws), &line) == SCORING_PARSE_OK);
  scoring_lookup(&scoring, 'a', 'c', &score, &is_match);
  ASSERT(score == -3);
  scoring_lookup(&scoring, 'C', 'C', &score, &is_match);
  ASSERT(score == 7);

  scoring_init(&scoring, 0, 0, -4, -1, false, false, false, false, false, false);
  ASSERT(scoring_parse_matrix(&scoring, comma, strlen(comma), &line) == SCORING_PARSE_OK);
  scoring_lookup(&scoring, 'a', 'a', &score, &is_match);
  ASSERT(score == 5);
  scoring_lookup(&scoring, 'c', 'a', &score, &is_match);
  ASSERT(score == -3);
  scoring_lookup(&scoring, 'c', 'c', &score, &is_match);
  ASSERT(score == 7);
}

/* Every status of the scoring file parsers, with the line it names, and
 * nothing changed by a failed parse (the tools' grammar: a matrix score is
 * a number prefix checked against what follows it, a pairs score must be
 * the rest of its line) */
void scoring_test_parse_statuses()
{
  static const struct {
    bool matrix;
    const char *text;
    scoring_parse_t status;
    int line;
  } cases[] = {
    {true,  " A C\nA 1 -1\nC -1 1\n",       SCORING_PARSE_OK, -1},
    {true,  ",A,C\r\nA,1,-1\r\nC,-1,+1\r\n", SCORING_PARSE_OK, -1},
    {true,  "# c\n\n  \n",                   SCORING_PARSE_EMPTY, -1},
    {true,  "A\n",                           SCORING_PARSE_FEW_HEADINGS, 1},
    {true,  "1A1C\nA111\n",                  SCORING_PARSE_SEPARATOR, 1},
    {true,  ",A;C\nA,1,1\n",                 SCORING_PARSE_SEP_MISSING, 1},
    {true,  ",A,C\nA;1,1\n",                 SCORING_PARSE_SEP_MISSING, 2},
    {true,  ",A,C,\nA,1,1\n",                SCORING_PARSE_EMPTY_HEADING, 1},
    {true,  " A C\nA 1\n",                   SCORING_PARSE_FEW_COLUMNS, 2},
    {true,  ",A,C\nA,1\n",                   SCORING_PARSE_FEW_COLUMNS, 2},
    {true,  " A C\nA 1x 2\n",                SCORING_PARSE_EXPECTED_WS, 2},
    {true,  " A\nA x\n",                     SCORING_PARSE_NUMBER, 2},
    {true,  ",A\nA,\n",                      SCORING_PARSE_NUMBER, 2},
    {true,  " A\nA 99999999999\n",           SCORING_PARSE_RANGE, 2},
    {true,  " A\n# c\nA 1 2\n",              SCORING_PARSE_MANY_COLUMNS, 3},
    {true,  ",A\nA,1,2\n",                   SCORING_PARSE_MANY_COLUMNS, 2},
    {true,  " A C\n",                        SCORING_PARSE_NO_ROWS, -1},
    {false, "A C 3\nC,A,-3\n",               SCORING_PARSE_OK, -1},
    {false, "A C 3\nAC 3\n",                 SCORING_PARSE_SHORT_LINE, 2},
    {false, "A,C;3\n",                       SCORING_PARSE_SEPARATORS, 1},
    {false, "A C 3 \n",                      SCORING_PARSE_INVALID_NUMBER, 1},
    {false, "A C 3000000000\n",              SCORING_PARSE_INVALID_NUMBER, 1},
    {false, "A C 2000000000\n",              SCORING_PARSE_RANGE, 1},
    {false, "# c\n\n",                       SCORING_PARSE_NO_PAIRS, -1},
  };
  size_t i;
  for(i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    scoring_t scoring, before;
    scoring_init(&scoring, 1, -1, -4, -1, false, false, false, false, false,
                 false);
    memcpy(&before, &scoring, sizeof(scoring));
    int line = 0;
    const char *t = cases[i].text;
    scoring_parse_t st = cases[i].matrix
      ? scoring_parse_matrix(&scoring, t, strlen(t), &line)
      : scoring_parse_pairs(&scoring, t, strlen(t), &line);
    ASSERT(st == cases[i].status);
    ASSERT(line == cases[i].line);
    if(st != SCORING_PARSE_OK)
      ASSERT(memcmp(&scoring, &before, sizeof(scoring)) == 0);
    ASSERT(strcmp(scoring_parse_strerror(st), "Unknown error") != 0);
  }

  // Out of memory: the copy of the scoring is the first allocation
  scoring_t scoring, before;
  scoring_init(&scoring, 1, -1, -4, -1, false, false, false, false, false,
               false);
  memcpy(&before, &scoring, sizeof(scoring));
  int line = 0;
  alignment_alloc_fail(1);
  ASSERT(scoring_parse_pairs(&scoring, "A C 3\n", 6, &line) ==
         SCORING_PARSE_NOMEM);
  alignment_alloc_fail(0);
  ASSERT(memcmp(&scoring, &before, sizeof(scoring)) == 0);
}

/* Every built-in scoring system is found by name, in any case; an unknown
 * name is refused; the one found is the table itself */
void scoring_test_system_by_name()
{
  scoring_t scoring, direct;
  char names[128], lower[64];
  size_t n = 0, i;
  strcpy(names, scoring_system_names());
  for(char *name = strtok(names, "|"); name != NULL; name = strtok(NULL, "|")) {
    ASSERT(scoring_system_by_name(&scoring, name));
    for(i = 0; name[i] != '\0'; i++) lower[i] = (char)tolower((unsigned char)name[i]);
    lower[i] = '\0';
    ASSERT(scoring_system_by_name(&scoring, lower));
    n++;
  }
  ASSERT(n == 5);
  ASSERT(!scoring_system_by_name(&scoring, "PAM31"));
  ASSERT(!scoring_system_by_name(&scoring, ""));
  int s1, s2;
  bool m1, m2;
  scoring_system_BLOSUM62(&direct);
  ASSERT(scoring_system_by_name(&scoring, "blosum62"));
  scoring_lookup(&scoring, 'W', 'W', &s1, &m1);
  scoring_lookup(&direct, 'W', 'W', &s2, &m2);
  ASSERT(s1 == s2 && s1 == 11 && m1 && m2);
}

void test_scoring(void)
{
  SUITE_START("Scoring");
  RUN(scoring_test_load_matrix());
  RUN(scoring_test_parse_statuses());
  RUN(scoring_test_penalty_range());
  RUN(scoring_test_case_sensitive());
  RUN(scoring_test_system_by_name());
  SUITE_END();
}
