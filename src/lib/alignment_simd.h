/*
 alignment_simd.h
 url: https://github.com/DanielTillett/nucleics-align
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

// Whether the target has a vector max instruction for 32-bit integers.
// Without one (plain SSE2, the x86-64 baseline) vmax() takes four
// instructions, which still pays off in the Needleman-Wunsch/Smith-Waterman
// fill but not in hirschberg.c
#if defined(__SSE4_1__) || defined(__AVX2__) || defined(__ARM_NEON) ||     defined(__aarch64__)
#  define ALN_NATIVE_VMAX 1
#endif

// x86-64 Linux builds without -mavx2: compile the SIMD kernels for AVX2,
// SSE4.1 and the baseline, and let the loader pick one for the CPU
// (GCC/clang function multi-versioning, which needs ifunc support).
// Not under ThreadSanitizer: the ifunc resolvers run before it starts and
// crash the program
#if defined(__has_feature)
#  if __has_feature(thread_sanitizer)
#    define ALN_TSAN 1
#  endif
#endif
#if defined(__SANITIZE_THREAD__)
#  define ALN_TSAN 1
#endif
#if defined(__x86_64__) && defined(__linux__) && !defined(__AVX2__) && \
    !defined(ALN_TSAN) && defined(__has_attribute)
#  if __has_attribute(target_clones)
#    define ALN_DISPATCH 1
#    define ALN_TARGET_CLONES \
       __attribute__((target_clones("avx2", "sse4.1", "default")))
#  endif
#endif
#ifndef ALN_TARGET_CLONES
#  define ALN_TARGET_CLONES
#endif

// 8 lanes with AVX2, otherwise 4 (native width for SSE2 and NEON). GCC
// dispatch builds use 8 in every version: the SSE4.1 and baseline versions
// run each vector as two halves, which still beats 4 lanes (fewer segments,
// so the lazy F loop wraps less; 20 kb unrelated score pass, one thread,
// SSE4.1: 0.21 -> 0.11 s). clang rejects 8-lane vectors in the versions
// without AVX, so its dispatch builds keep 4
// AVX-512 targets (e.g. NATIVE=1 on Zen 4/5 or Ice Lake and later) use 16:
// a 20 kb unrelated score pass on one thread 0.10 -> 0.055 s. Not in the
// dispatch builds, where every version shares one lane count and 16 lanes
// run 10x slower on AVX2
#ifndef ALN_VL
#  if defined(__AVX512F__)
#    define ALN_VL 16
#  elif defined(__AVX2__) || (defined(ALN_DISPATCH) && !defined(__clang__))
#    define ALN_VL 8
#  else
#    define ALN_VL 4
#  endif
#endif

typedef score_t vscore_t __attribute__((vector_size(ALN_VL * sizeof(score_t))));

// Helpers that take or return vectors are always inlined, even without
// optimisation: an 8-lane vector passed between the AVX2 version of a
// function and a helper compiled for the baseline would be passed in
// different registers on each side (a crash in DEBUG=1 builds)
#define ALN_INLINE static inline __attribute__((always_inline))

ALN_INLINE vscore_t vsplat(score_t x) { return (vscore_t){0} + x; }

ALN_INLINE vscore_t vmax(vscore_t a, vscore_t b)
{
  vscore_t gt = a > b; // all bits set in lanes where a > b
  return (a & gt) | (b & ~gt);
}

// Move every lane up by one; lane 0 gets x
ALN_INLINE vscore_t vshift(vscore_t v, score_t x)
{
#if ALN_VL == 16
  return __builtin_shufflevector(v, vsplat(x), 16, 0, 1, 2, 3, 4, 5, 6, 7, 8,
                                 9, 10, 11, 12, 13, 14);
#elif ALN_VL == 8
  return __builtin_shufflevector(v, vsplat(x), 8, 0, 1, 2, 3, 4, 5, 6);
#else
  return __builtin_shufflevector(v, vsplat(x), 4, 0, 1, 2);
#endif
}

ALN_INLINE bool vany(vscore_t v)
{
  score_t r = 0;
  int l;
  for(l = 0; l < ALN_VL; l++) r |= v[l];
  return r != 0;
}

#endif /* ALN_STRIPED */

#endif /* ALIGNMENT_SIMD_HEADER_SEEN */
