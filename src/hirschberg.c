/*
 hirschberg.c
 url: https://github.com/DanielTillett/seq-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// Linear space global alignment with affine gaps, following the 'diff'
// procedure of Myers & Miller (1988) Comput Appl Biosci 4:11-17, written to
// maximise a score rather than minimise a cost.
// A gap of length k scores gap_open + k*gap_extend, as in needleman_wunsch.c.
// Deleting a character of a means a gap in b, inserting one of b a gap in a.

#include <stdlib.h>
#include <stddef.h> // ptrdiff_t
#include <stdio.h>
#include <string.h>

#include "alignment_macros.h"
#include "alignment_simd.h"
#include "hirschberg.h"

// SIMD passes only where the target has a native vector max: without one
// (plain SSE2, the x86-64 baseline) the scalar loops are faster here
#if defined(ALN_STRIPED) && defined(ALN_NATIVE_VMAX)
#  define HB_STRIPED 1
#endif

typedef struct
{
  const scoring_t *scoring;
  sub_table_t table;
  const char *b; // whole of seq b, to index row_b
  size_t *row_b; // table row offset for each character of b
  score_t *cc, *dd, *rr, *ss; // forward and reverse rows, len_b+1 each
  char *out_a, *out_b;
  size_t pos; // next position to write in out_a, out_b
#ifdef HB_STRIPED
  void *vbuf; // striped rows and score profile, see hb_pass_striped()
  bool simd_ok; // all scores far enough from HB_NEG and overflow
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

static void hb_emit_del(hb_state_t *st, const char *a, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++, st->pos++) {
    st->out_a[st->pos] = a[i];
    st->out_b[st->pos] = '-';
  }
}

static void hb_emit_ins(hb_state_t *st, const char *b, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++, st->pos++) {
    st->out_a[st->pos] = '-';
    st->out_b[st->pos] = b[i];
  }
}

#ifdef HB_STRIPED

// Below any real score, with room to add penalties without overflow
#define HB_NEG (SCORE_MIN / 4)

// Striped SIMD (Farrar 2007) version of one pass of hb_diff(): gives exactly
// the same values as the scalar loops. Rows run over `rows` characters of a
// starting at a[0] and moving by astep (+1 forward, -1 reverse); columns run
// over n characters of b the same way (bstep). t0 is tb or te. The last row
// is written to oc[p*ostep] (best score) and od[p*ostep] (ending in a
// deletion) for column p = 0..n in pass order.
// The insertion (horizontal) gap feeds back into the best score of the same
// row, so Farrar's lazy F loop carries it across segment boundaries. That is
// exact only for gap_open <= 0 and gap_extend <= 0, which the caller checks.
static void hb_pass_striped(hb_state_t *st, const char *a, ptrdiff_t astep,
                            size_t rows, const char *b, ptrdiff_t bstep,
                            size_t n, score_t t0,
                            score_t *oc, score_t *od, ptrdiff_t ostep)
{
  const scoring_t *scoring = st->scoring;
  const sub_table_t *t = &st->table;
  const score_t g = scoring->gap_open, h = scoring->gap_extend;
  const vscore_t vg = vsplat(g), vh = vsplat(h), vneg = vsplat(HB_NEG);
  const size_t seg = (n + ALN_VL - 1) / ALN_VL;
  vscore_t *H = st->vbuf, *Hn = H + seg, *D = Hn + seg, *prof = D + seg;
  vscore_t *tmp, vdiag, vf, vd, vc;
  size_t r, k, l, p, c;

  // prof[c*seg + k] lane l: score of a character with code c against the
  // b character in column l*seg+k (0 for padding past the end of b)
  memset(prof, 0, t->n * seg * sizeof(vscore_t));
  for(p = 0; p < n; p++) {
    size_t row = t->code[(unsigned char)b[(ptrdiff_t)p * bstep]] * t->n;
    l = p / seg; k = p - l*seg;
    for(c = 0; c < t->n; c++) prof[c*seg + k][l] = t->scores[row + c];
  }

  // Row 0: only insertions
  for(k = 0; k < seg; k++)
    for(l = 0; l < ALN_VL; l++) {
      H[k][l] = g + (score_t)(l*seg + k + 1) * h;
      D[k][l] = H[k][l] + g;
    }

  score_t col0 = t0, col0_prev = 0; // column 0: only deletions

  for(r = 0; r < rows; r++)
  {
    const vscore_t *pr = prof +
      (size_t)t->code[(unsigned char)a[(ptrdiff_t)r * astep]] * seg;
    col0 += h;

    vdiag = vshift(H[seg-1], col0_prev);
    vf = vshift(vneg, col0 + g + h); // lanes after 0 are fixed below
    for(k = 0; k < seg; k++)
    {
      vd = vmax(D[k], H[k] + vg) + vh;
      vc = vmax(vmax(vd, vf), vdiag + pr[k]);
      vdiag = H[k];
      D[k] = vd;
      Hn[k] = vc;
      vf = vmax(vf, vc + vg) + vh;
    }

    // Lazy F: carry the insertion gap across segment boundaries while it
    // could still raise a best score
    vf = vshift(vf, HB_NEG);
    k = 0;
    while(vany(vf > Hn[k] + vg)) {
      Hn[k] = vmax(Hn[k], vf);
      vf = vmax(vf + vh, vneg);
      if(++k == seg) { k = 0; vf = vshift(vf, HB_NEG); }
    }

    tmp = H; H = Hn; Hn = tmp;
    col0_prev = col0;
  }

  oc[0] = od[0] = col0;
  for(l = 0; l < ALN_VL; l++)
    for(k = 0, p = l*seg; k < seg && p < n; k++, p++) {
      oc[(ptrdiff_t)(p+1) * ostep] = H[k][l];
      od[(ptrdiff_t)(p+1) * ostep] = D[k][l];
    }
}

#endif /* HB_STRIPED */

// Forward and reverse passes of hb_diff() without SIMD
static void hb_passes_scalar(hb_state_t *st, const char *a, const char *b,
                             size_t m, size_t n, size_t imid,
                             score_t tb, score_t te)
{
  const score_t g = st->scoring->gap_open, h = st->scoring->gap_extend;
  score_t *cc = st->cc, *dd = st->dd, *rr = st->rr, *ss = st->ss;
  score_t s, c, e, t;
  size_t i, j;

  // Forward pass: best scores for a[0..imid) against b[0..j),
  // cc ending in any state, dd ending in a deletion
  cc[0] = 0;
  t = g;
  for(j = 1; j <= n; j++) { t += h; cc[j] = t; dd[j] = t + g; }

  t = tb;
  // Substitution scores: column of the table for a[i], row for each b[j]
  const int *col;
  const size_t *row_b = st->row_b + (b - st->b);

  for(i = 1; i <= imid; i++)
  {
    col = st->table.scores + st->table.code[(unsigned char)a[i-1]];
    s = cc[0];
    t += h;
    c = t;
    cc[0] = c;
    e = t + g;
    for(j = 1; j <= n; j++)
    {
      e = MAX2(e, c + g) + h;
      dd[j] = MAX2(dd[j], cc[j] + g) + h;
      c = MAX3(dd[j], e, s + col[row_b[j-1]]);
      s = cc[j];
      cc[j] = c;
    }
  }
  dd[0] = cc[0];

  // Reverse pass: best scores for a[imid..m) against b[j..n),
  // rr starting in any state, ss starting with a deletion
  rr[n] = 0;
  t = g;
  for(j = n; j-- > 0; ) { t += h; rr[j] = t; ss[j] = t + g; }

  t = te;
  for(i = m; i-- > imid; )
  {
    col = st->table.scores + st->table.code[(unsigned char)a[i]];
    s = rr[n];
    t += h;
    c = t;
    rr[n] = c;
    e = t + g;
    for(j = n; j-- > 0; )
    {
      e = MAX2(e, c + g) + h;
      ss[j] = MAX2(ss[j], rr[j] + g) + h;
      c = MAX3(ss[j], e, s + col[row_b[j]]);
      s = rr[j];
      rr[j] = c;
    }
  }
  ss[n] = rr[n];
}

bool hirschberg_check_striped = false;
size_t hirschberg_striped_mismatches = 0;

#ifdef HB_STRIPED
// Testing aid: redo the passes without SIMD and count any difference
static void hb_check_striped(hb_state_t *st, const char *a, const char *b,
                             size_t m, size_t n, size_t imid,
                             score_t tb, score_t te)
{
  size_t bytes = (n+1) * sizeof(score_t);
  score_t *copy = malloc(4 * bytes);
  if(copy == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  memcpy(copy, st->cc, bytes);
  memcpy(copy + (n+1), st->dd, bytes);
  memcpy(copy + 2*(n+1), st->rr, bytes);
  memcpy(copy + 3*(n+1), st->ss, bytes);
  hb_passes_scalar(st, a, b, m, n, imid, tb, te);
  if(memcmp(copy, st->cc, bytes) || memcmp(copy + (n+1), st->dd, bytes) ||
     memcmp(copy + 2*(n+1), st->rr, bytes) ||
     memcmp(copy + 3*(n+1), st->ss, bytes))
  {
    hirschberg_striped_mismatches++;
  }
  free(copy);
}
#endif

// Align a[0..m) with b[0..n), appending to the output.
// tb / te: score of opening a deletion that touches the start / end of this
// sub-problem; 0 if it continues a deletion already paid for by the caller.
static void hb_diff(hb_state_t *st, const char *a, const char *b,
                    size_t m, size_t n, score_t tb, score_t te)
{
  const scoring_t *scoring = st->scoring;
  const score_t g = scoring->gap_open, h = scoring->gap_extend;
  size_t j;

  if(n == 0) { hb_emit_del(st, a, m); return; }
  if(m == 0) { hb_emit_ins(st, b, n); return; }

  if(m == 1)
  {
    // Either delete a[0] and insert all of b, or align a[0] with some b[j-1]
    score_t best = MAX2(tb, te) + h + hb_gap(scoring, n), c;
    size_t best_j = 0;

    for(j = 1; j <= n; j++)
    {
      c = hb_gap(scoring, j-1) + hb_sub(scoring, a[0], b[j-1]) +
          hb_gap(scoring, n-j);
      if(c > best) { best = c; best_j = j; }
    }

    if(best_j == 0)
    {
      // Put the deletion next to whichever neighbouring deletion it joins
      if(tb >= te) { hb_emit_del(st, a, 1); hb_emit_ins(st, b, n); }
      else         { hb_emit_ins(st, b, n); hb_emit_del(st, a, 1); }
    }
    else
    {
      hb_emit_ins(st, b, best_j-1);
      st->out_a[st->pos] = a[0];
      st->out_b[st->pos] = b[best_j-1];
      st->pos++;
      hb_emit_ins(st, b+best_j, n-best_j);
    }
    return;
  }

  score_t *cc = st->cc, *dd = st->dd, *rr = st->rr, *ss = st->ss;
  score_t c;
  size_t imid = m / 2;

  #ifdef HB_STRIPED
  // Building the score profile costs O(n x alphabet), so only use SIMD on
  // sub-problems with enough rows to pay for it
  if(alignment_striped_enabled && st->simd_ok && g <= 0 && h <= 0 &&
     n >= 2*ALN_VL && m >= 32 && m >= 2*st->table.n)
  {
    hb_pass_striped(st, a, 1, imid, b, 1, n, tb, cc, dd, 1);
    hb_pass_striped(st, a+m-1, -1, m-imid, b+n-1, -1, n, te, rr+n, ss+n, -1);
    if(hirschberg_check_striped) hb_check_striped(st, a, b, m, n, imid, tb, te);
  }
  else
  #endif
  {
    hb_passes_scalar(st, a, b, m, n, imid, tb, te);
  }

  // Where does the optimal path cross row imid? Type 2 means a deletion
  // spans the row: a[imid-1] and a[imid] are both deleted
  score_t midc = cc[0] + rr[0];
  size_t midj = 0;
  int type = 1;

  for(j = 0; j <= n; j++)
  {
    c = cc[j] + rr[j];
    if(c > midc) { midc = c; midj = j; type = 1; }
    c = dd[j] + ss[j] - g; // gap open counted on both sides
    if(c > midc) { midc = c; midj = j; type = 2; }
  }

  if(type == 1)
  {
    hb_diff(st, a, b, imid, midj, tb, g);
    hb_diff(st, a+imid, b+midj, m-imid, n-midj, g, te);
  }
  else
  {
    hb_diff(st, a, b, imid-1, midj, tb, 0);
    hb_emit_del(st, a+imid-1, 2);
    hb_diff(st, a+imid+1, b+midj, m-imid-1, n-midj, 0, te);
  }
}

score_t hirschberg_score_alignment(const char *result_a, const char *result_b,
                                   const scoring_t *scoring)
{
  const score_t g = scoring->gap_open, h = scoring->gap_extend;
  enum Matrix prev = MATCH;
  score_t score = 0;
  size_t i;

  for(i = 0; result_a[i] != '\0' && result_b[i] != '\0'; i++)
  {
    if(result_a[i] == '-')
    {
      score += (prev == GAP_A ? h : g + h);
      prev = GAP_A;
    }
    else if(result_b[i] == '-')
    {
      score += (prev == GAP_B ? h : g + h);
      prev = GAP_B;
    }
    else
    {
      score += hb_sub(scoring, result_a[i], result_b[i]);
      prev = MATCH;
    }
  }

  return score;
}

bool hirschberg_align(const char *a, const char *b,
                      const scoring_t *scoring, alignment_t *result)
{
  return hirschberg_align2(a, b, strlen(a), strlen(b), scoring, result);
}

bool hirschberg_align2(const char *a, const char *b,
                       size_t len_a, size_t len_b,
                       const scoring_t *scoring, alignment_t *result)
{
  if(scoring->no_start_gap_penalty || scoring->no_end_gap_penalty ||
     scoring->no_gaps_in_a || scoring->no_gaps_in_b ||
     scoring->no_gaps_in_longest || scoring->no_mismatches)
  {
    return false;
  }

  alignment_ensure_capacity(result, len_a + len_b);

  hb_state_t st;
  st.scoring = scoring;
  st.out_a = result->result_a;
  st.out_b = result->result_b;
  st.pos = 0;

  score_t *rows = malloc(4 * (len_b+1) * sizeof(score_t));
  if(rows == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  sub_table_build(&st.table, scoring, a, len_a, b, len_b);
  st.b = b;
  st.row_b = malloc((len_b+1) * sizeof(size_t));
  if(st.row_b == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  size_t j;
  for(j = 0; j < len_b; j++)
    st.row_b[j] = st.table.code[(unsigned char)b[j]] * st.table.n;

#ifdef HB_STRIPED
  // Three striped rows plus a profile row per character code
  size_t seg = (len_b + ALN_VL - 1) / ALN_VL + 1;
  // HB_NEG must stay below every real score, with room to add penalties.
  // No path through the matrices (including padding columns) has more than
  // len_a + len_b + ALN_VL steps, each scoring at most |gap_open| +
  // |gap_extend| or the largest substitution score in magnitude. Beyond a
  // conservative bound use the scalar passes
  {
    #define HB_DABS(x) ((x) < 0 ? -(double)(x) : (double)(x))
    double gap = HB_DABS(scoring->gap_open) + HB_DABS(scoring->gap_extend);
    double sub = 0;
    size_t c;
    for(c = 0; c < st.table.n * st.table.n; c++)
      sub = MAX2(sub, HB_DABS(st.table.scores[c]));
    double steps = (double)len_a + (double)len_b + ALN_VL + 2;
    st.simd_ok = steps * (gap + sub) < (double)(1 << 27);
  }

  st.vbuf = aligned_alloc(sizeof(vscore_t),
                          (3 + st.table.n) * seg * sizeof(vscore_t));
  if(st.vbuf == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
#endif

  st.cc = rows;
  st.dd = rows + (len_b+1);
  st.rr = rows + 2*(len_b+1);
  st.ss = rows + 3*(len_b+1);

  hb_diff(&st, a, b, len_a, len_b, scoring->gap_open, scoring->gap_open);
  free(rows);
  free(st.row_b);
#ifdef HB_STRIPED
  free(st.vbuf);
#endif
  sub_table_free(&st.table);

  result->result_a[st.pos] = '\0';
  result->result_b[st.pos] = '\0';
  result->length = st.pos;
  result->score = hirschberg_score_alignment(result->result_a,
                                             result->result_b, scoring);
  return true;
}
