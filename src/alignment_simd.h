/*
 alignment_simd.h
 url: https://github.com/DanielTillett/seq-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// Small vector helpers for the striped (Farrar 2007) dynamic programming in
// alignment.c and hirschberg.c. Uses GCC/clang vector extensions rather than
// intrinsics so it builds on x86 and ARM; with -mavx2 (or -march=native) the
// compilers emit full-width instructions such as vpmaxsd.
// ALN_STRIPED is defined only if the compiler supports them.

#ifndef ALIGNMENT_SIMD_HEADER_SEEN
#define ALIGNMENT_SIMD_HEADER_SEEN

#include <stdbool.h>
#include "alignment_scoring.h" // score_t

#if defined(__has_builtin)
#  if __has_builtin(__builtin_shufflevector)
#    define ALN_STRIPED 1
#  endif
#endif

#ifdef ALN_STRIPED

// 8 lanes with AVX2, otherwise 4 (native width for SSE2 and NEON)
#ifdef __AVX2__
#  define ALN_VL 8
#else
#  define ALN_VL 4
#endif
// Whether the target has a vector max instruction for 32-bit integers.
// Without one (plain SSE2, the x86-64 baseline) vmax() takes four
// instructions, which still pays off in the Needleman-Wunsch/Smith-Waterman
// fill but not in hirschberg.c
#if defined(__SSE4_1__) || defined(__AVX2__) || defined(__ARM_NEON) ||     defined(__aarch64__)
#  define ALN_NATIVE_VMAX 1
#endif

typedef score_t vscore_t __attribute__((vector_size(ALN_VL * sizeof(score_t))));

static inline vscore_t vsplat(score_t x) { return (vscore_t){0} + x; }

static inline vscore_t vmax(vscore_t a, vscore_t b)
{
  vscore_t gt = a > b; // all bits set in lanes where a > b
  return (a & gt) | (b & ~gt);
}

// Move every lane up by one; lane 0 gets x
static inline vscore_t vshift(vscore_t v, score_t x)
{
#if ALN_VL == 8
  return __builtin_shufflevector(v, vsplat(x), 8, 0, 1, 2, 3, 4, 5, 6);
#else
  return __builtin_shufflevector(v, vsplat(x), 4, 0, 1, 2);
#endif
}

static inline bool vany(vscore_t v)
{
  score_t r = 0;
  int l;
  for(l = 0; l < ALN_VL; l++) r |= v[l];
  return r != 0;
}

#endif /* ALN_STRIPED */

#endif /* ALIGNMENT_SIMD_HEADER_SEEN */
