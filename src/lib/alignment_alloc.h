/*
 alignment_alloc.h
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// Private: the allocations an alignment makes, which the *_checked()
// aligners recover from. Tests make one fail on purpose (see
// alignment_alloc_fail() in alignment_testing.h)

#ifndef ALIGNMENT_ALLOC_HEADER_SEEN
#define ALIGNMENT_ALLOC_HEADER_SEEN

#include <stdlib.h>
#include <stdbool.h>
#include <stdatomic.h>

extern atomic_size_t alignment_alloc_fail_in;

// Whether this allocation is the one a test asked to fail
static inline bool aln_alloc_fails(void)
{
  size_t n = atomic_load_explicit(&alignment_alloc_fail_in,
                                  memory_order_relaxed);
  while(n > 0 && !atomic_compare_exchange_weak(&alignment_alloc_fail_in,
                                               &n, n - 1)) {}
  return n == 1;
}

static inline void* aln_malloc(size_t size)
{
  return aln_alloc_fails() ? NULL : malloc(size);
}

static inline void* aln_calloc(size_t n, size_t size)
{
  return aln_alloc_fails() ? NULL : calloc(n, size);
}

static inline void* aln_realloc(void *ptr, size_t size)
{
  return aln_alloc_fails() ? NULL : realloc(ptr, size);
}

static inline void* aln_aligned_alloc(size_t align, size_t size)
{
  return aln_alloc_fails() ? NULL : aligned_alloc(align, size);
}

#endif /* ALIGNMENT_ALLOC_HEADER_SEEN */
