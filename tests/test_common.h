/*
 tests/test_common.h
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026

 The unit tests' harness: suites, the ASSERT that counts and continues,
 one random generator reseeded for every test from the run's seed, random
 sequences and mutations, alignment comparisons, and a snapshot of the
 library's tuning globals taken before each test and restored after it
 */

#ifndef TEST_COMMON_HEADER_SEEN
#define TEST_COMMON_HEADER_SEEN

// request decent POSIX version (open_memstream, getpid)
#define _XOPEN_SOURCE 700

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <limits.h>
#include <ctype.h>
#include <unistd.h>
#include <pthread.h>

#include "alignment.h"
#include "alignment_macros.h"
#include "alignment_simd.h" // ALN_STRIPED
#include "alignment_testing.h" // counters, test hooks and the OOM switch
#include "needleman_wunsch.h"
#include "smith_waterman.h"
#include "hirschberg.h"

//
// Suites and assertions
//

// Counts (read by main for the summary)
extern int tests_run, tests_failed;
extern int suites_run, suites_failed, suites_empty;

// A failed assertion names the test it is in and the run's seed, and the
// test goes on (a crash would hide the other failures)
#define QUOTE(str) #str
#define ASSERT(x) do { test_counted(); if(!(x)) \
  test_failed(__FILE__, __LINE__, QUOTE(x)); } while(0)

// Checks on how a fill was made that hold only where the striped fill
// exists (the scalar-only build, -DALN_NO_SIMD, still checks the results)
#ifdef ALN_STRIPED
#  define ASSERT_STRIPED(x) ASSERT(x)
#else
#  define ASSERT_STRIPED(x) ((void)0)
#endif

void test_counted(void);
void test_failed(const char *file, int line, const char *expr);

// A suite: SUITE_START(name); RUN(test()); ...; SUITE_END();
void SUITE_START(const char *name);
void SUITE_END(void);

// One test of a suite: the generator is reseeded from the run's seed and
// the test's name (so a test's random data does not depend on the tests
// before it), and the tuning globals are saved before and restored after
// (a test that left one changed is reported)
#define RUN(call) do { test_begin(#call); call; test_end(#call); } while(0)
void test_begin(const char *name);
void test_end(const char *name);

// The run's seed (NUCLEICS_ALIGN_TEST_SEED, or the time and process id)
void test_set_seed(unsigned int seed);
unsigned int test_seed(void);

//
// Random data (deterministic from the seed)
//

// 0 .. n-1
unsigned int test_rand(unsigned int n);
// len letters of alphabet, NUL-terminated (seq holds len + 1 bytes)
void test_rand_seq(char *seq, size_t len, const char *alphabet);
// A copy of in with about per_mille substitutions, a quarter as many
// deletions and insertions, and now and then a deleted run of 20. out
// holds at least 3/2 the length of in, plus one
void test_mutate(char *out, const char *in, unsigned int per_mille,
                 const char *alphabet);

//
// Comparisons
//

// Removing the gaps from an alignment string gives back the sequence
bool test_ungapped_equals(const char *aln, const char *seq);
// The same score, strings and span
bool alignments_equal(const alignment_t *x, const alignment_t *y);

// The striped SIMD fill must give exactly the scalar matrices
void simd_check(const char *a, const char *b, const scoring_t *scoring,
                char is_sw, aligner_t *scalar, aligner_t *striped);
// Linear-space alignment must give the same alignment with or without SIMD
// (a positive gap open is refused by both)
void simd_hb_check(const char *a, const char *b, const scoring_t *scoring,
                   alignment_t *scalar, alignment_t *striped);

//
// Suites (one per test_*.c file)
//

void test_nw(void);
void test_hirschberg(void);
void test_simd_fill(void);
void test_simd_hirschberg(void);
void test_sw_compact(void);
void test_sw(void);
void test_buffers(void);
void test_scoring(void);
void test_options(void);
void test_errors(void);
void test_wfa(void);
void test_regressions(void);

#endif /* TEST_COMMON_HEADER_SEEN */
