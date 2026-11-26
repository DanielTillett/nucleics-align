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
  RUN(scoring_test_penalty_range());
  RUN(scoring_test_case_sensitive());
  RUN(scoring_test_system_by_name());
  SUITE_END();
}
