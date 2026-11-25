/*
 hirschberg_internal.h
 Internal cross-module seam for hirschberg.c and hirschberg_pass.c.
 */

#ifndef HIRSCHBERG_INTERNAL_HEADER_SEEN
#define HIRSCHBERG_INTERNAL_HEADER_SEEN

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdatomic.h>

#include "alignment_macros.h"
#include "alignment_simd.h"
#include "hirschberg.h"
#include "alignment_testing.h"
#include "alignment_alloc.h"

#if defined(ALN_STRIPED) && defined(ALN_NATIVE_VMAX)
#  define HB_STRIPED 1
#  define HB_FAST_VMAX() 1
#elif defined(ALN_STRIPED) && defined(ALN_DISPATCH)
#  define HB_STRIPED 1
#  define HB_FAST_VMAX() __builtin_cpu_supports("sse4.1")
#endif

#define HB_NO_BOUND SCORE_MIN
#define HB_CHAIN_MAX 8

#ifdef HB_STRIPED
#  define HB_NEG (SCORE_MIN / 4)
#  define HB_THREAD_CELLS 1e6
#  define HB_VL16 (2 * ALN_VL)
#  define HB16_NEG (-16000)
#  define HB16_RANGE 12000

typedef struct hb_top hb_top_t;
#endif

typedef struct {
  score_t *h[HB_CHAIN_MAX], *d[HB_CHAIN_MAX];
  size_t rows[HB_CHAIN_MAX], len;
} hb_chain_t;

typedef struct {
  const scoring_t *scoring;
  sub_table_t table;
  const char *b;
  size_t *row_b;
  score_t *cc, *dd, *rr, *ss;
  char *out_a, *out_b;
  size_t pos;
  uint8_t *ops;
  bool in_range;
  bool nomem;
  alignment_options_t opts;
#ifdef HB_STRIPED
  void *vbuf;
  void *vbuf2;
  size_t vbuf_bytes;
  bool simd_ok;
  int pass_threads;
  score_t smax;
  size_t w16;
#endif
} hb_state_t;

static inline score_t hb_gap(const scoring_t *scoring, size_t k)
{
  return k == 0 ? 0 : scoring->gap_open + (score_t)k * scoring->gap_extend;
}

static inline score_t hb_sub(const scoring_t *scoring, char a, char b)
{
  int score;
  bool is_match;
  scoring_lookup(scoring, a, b, &score, &is_match);
  return score;
}

static inline void hb_nomem(hb_state_t *st)
{
  __atomic_store_n(&st->nomem, true, __ATOMIC_RELAXED);
}

// Defined in hirschberg.c
score_t hb_midpoint(const hb_state_t *st, size_t n, size_t *midj, int *type);

// Defined in hirschberg_pass.c
void hb_passes_scalar(hb_state_t *st, const char *a, const char *b,
                      size_t m, size_t n, size_t imid,
                      score_t tb, score_t te,
                      bool do_f, bool do_r,
                      const hb_chain_t *fsave, const hb_chain_t *rsave);
#ifdef HB_STRIPED
bool hb_ensure_vbuf(hb_state_t *st, bool need_second);
bool hb_simd_ok(const hb_state_t *st, size_t m, size_t n);
void hb_pass_striped(hb_state_t *st, const char *a, ptrdiff_t astep,
                     size_t rows, const char *b, ptrdiff_t bstep,
                     size_t c0, size_t n, score_t t0,
                     score_t *oc, score_t *od, ptrdiff_t ostep,
                     void *buf, const hb_top_t *top, int threads);
void hb_check_striped(hb_state_t *st, const char *a, const char *b,
                      size_t m, size_t n, size_t imid,
                      score_t tb, score_t te, bool banded);
void hb_passes_striped(hb_state_t *st, const char *a, const char *b,
                       size_t m, size_t n, size_t imid,
                       score_t tb, score_t te, ptrdiff_t dlo, ptrdiff_t dhi,
                       bool do_f, bool do_r,
                       const hb_chain_t *fsave, const hb_chain_t *rsave);
bool hb_passes_banded(hb_state_t *st, const char *a, const char *b,
                      size_t m, size_t n, size_t imid,
                      score_t tb, score_t te, score_t lower);
#endif

#endif /* HIRSCHBERG_INTERNAL_HEADER_SEEN */
