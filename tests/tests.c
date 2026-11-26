/*
 tests/tests.c
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Feb 2015

 The unit tests' main: one suite per tests/test_*.c file (see
 tests/test_common.h for the harness). Build and run with make test
 */

#include <time.h>
#include "test_common.h"

int main(int argc, char **argv)
{
  if(argc != 1)
  {
    printf("  Unused args '%s..'\n", argv[1]);
    printf("Usage: ./nucleics_align_tests\n");
    exit(EXIT_FAILURE);
  }

  // The seed every test's random data comes from (with its name):
  // NUCLEICS_ALIGN_TEST_SEED repeats a run
  const char *seed_env = getenv("NUCLEICS_ALIGN_TEST_SEED");
  unsigned int seed = seed_env != NULL
                      ? (unsigned int)strtoul(seed_env, NULL, 10)
                      : (unsigned int)time(NULL) + (unsigned int)getpid();
  test_set_seed(seed);
  printf("  Random seed: %u (set NUCLEICS_ALIGN_TEST_SEED to repeat)\n", seed);

  printf("  Test Nucleics Align C library:\n\n");

  test_nw();
  test_hirschberg();
  test_simd_fill();
  test_simd_hirschberg();
  test_sw_compact();
  test_sw();
  test_buffers();
  test_scoring();
  test_options();
  test_errors();
  test_wfa();
  test_regressions();

  printf("\n");
  printf(" %i / %i suites failed\n", suites_failed, suites_run);
  printf(" %i / %i suites empty\n", suites_empty, suites_run);
  printf(" %i / %i tests failed\n", tests_failed, tests_run);

  printf("\n THE END.\n");

  return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
