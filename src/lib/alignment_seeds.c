/*
 alignment_seeds.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// A bound on the score any alignment of the rest of a and b can add from a
// cell, from seeds: the words of k letters of a at every kth position. A
// seed with no exact copy in b forces an error inside it; between seeds
// matched exactly, in order, a change of diagonal costs gap letters. This
// is the gap-chaining seed heuristic of A*PA (Groot Koerkamp & Ivanov
// 2024, Bioinformatics 40:btae032), with affine scores converted to
// penalties as in alignment_wfa.c (Eizenga & Paten 2022): a mismatch
// 2(l - mm), a gap open -2o, a gap letter l - 2e, for match score l, the
// best mismatch score mm, gap open o and extension e; twice the score of
// an alignment is l x (its letters) less its penalty. Each term below is a
// lower bound on the penalty of every alignment, so the score bound never
// excludes a best path (alignment_band.c narrows strips with it).
//
// A seed costs at least c when broken, c the least penalty of one error
// seen from the seeds: a mismatch, a gap in b inside the seed, or a run of
// r gapped letters of a, which breaks at most ceil((r-1)/k) + 1 seeds.
// Between consecutive matched seeds the penalty is at least the gap cost
// of the diagonal change, and at least c per seed skipped: max of the two.
// Chains are followed over the next W seeds; seeds further on count c
// each. The gap the path is in at the cell may go on without opening, so
// one gap open is given back. The bound is also never above the plain one
// (all letters left matched but the gap the lengths force)

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "alignment_seeds.h"
#include "alignment_alloc.h"
#include "alignment_macros.h"

#define SEEDS_WINDOW 64 // seeds ahead a chain is followed over
#define SEEDS_MAX_COPIES 32 // seeds with more copies in b are not charged
#define SEEDS_SAMPLE 8 // columns between bounds evaluated in a row
#define SEEDS_NONE (INT64_MAX / 4)

struct seeds_t {
  size_t k, S, M;
  size_t *sa;     // seed s starts at a[sa[s]]
  size_t *mstart; // its copies in b: mj[mstart[s] .. mstart[s+1]-1]
  size_t *mj;
  int64_t *G;     // least chain penalty from a copy to the end
  int64_t *Gs;    // least c x seed + G over copies of seeds from s on
  int64_t c, open, gamma, l, lip; // lip: most h changes per column
  int gap_extend;
  ptrdiff_t D;    // the end's diagonal, len_a - len_b
  size_t len_a, len_b;
  // Penalties sampled along the last bottom row, every SEEDS_SAMPLE
  // columns from column from, for the next strip's top row
  int64_t *last; size_t last_from, last_to, last_n, last_y;
};

typedef struct { uint64_t key; uint32_t pos; } kmer_t;

static int kmer_cmp(const void *pa, const void *pb)
{
  const kmer_t *a = pa, *b = pb;
  return a->key < b->key ? -1 : a->key > b->key ? 1 :
         a->pos < b->pos ? -1 : a->pos > b->pos;
}

// Penalty of moving from one diagonal to another
static inline int64_t seeds_gap(const seeds_t *sd, ptrdiff_t delta)
{
  if(delta == 0) return 0;
  return sd->open + sd->gamma * (int64_t)(delta < 0 ? -delta : delta);
}

seeds_t *seeds_build(const sub_table_t *t, const scoring_t *sc,
                     const char *a, size_t len_a,
                     const char *b, size_t len_b, score_t smax, size_t k)
{
  const size_t n = t->n;
  if(n == 0 || n > 64 || sc->gap_open > 0 || sc->gap_extend > 0) return NULL;
  // Letters must fall into classes: a letter of a and one of b match
  // exactly when they are in one class, and then score smax. The table
  // holds only pairs of a letter of a and one of b, so a's letters are
  // grouped by the set of b's letters they match (two sets that share a
  // letter must be the same), and the class of a letter is the lowest code
  // in its set: for one of b, of the set holding it. Letters of a matching
  // none of b get a class of their own
  bool in_a[64] = {false}, in_b[64] = {false};
  size_t i, x, y;
  for(i = 0; i < len_a; i++) in_a[t->code[(unsigned char)a[i]]] = true;
  for(i = 0; i < len_b; i++) in_b[t->code[(unsigned char)b[i]]] = true;
  uint64_t set[64] = {0};
  int cls[64], mm = INT_MIN;
  for(x = 0; x < n; x++) cls[x] = -1;
  for(x = 0; x < n; x++) {
    if(!in_a[x]) continue;
    for(y = 0; y < n; y++) {
      if(!in_b[y]) continue;
      if(t->is_match[y*n + x]) {
        if(t->scores[y*n + x] != smax) return NULL;
        set[x] |= (uint64_t)1 << y;
      }
      else mm = MAX2(mm, t->scores[y*n + x]);
    }
  }
  if(mm == INT_MIN || mm >= smax) return NULL;
  for(x = 0; x < n; x++) {
    if(!in_a[x]) continue;
    size_t x2;
    for(x2 = 0; x2 < x; x2++)
      if(in_a[x2] && (set[x] & set[x2]) && set[x] != set[x2]) return NULL;
    int id = (int)n + (int)x;
    if(set[x]) { id = 0; while(!((set[x] >> id) & 1)) id++; }
    cls[x] = id;
    for(y = 0; y < n; y++) {
      if(!((set[x] >> y) & 1)) continue;
      if(cls[y] >= 0 && cls[y] != id) return NULL; // y in a as well
      cls[y] = id;
    }
  }
  for(y = 0; y < n; y++) if(cls[y] < 0) cls[y] = (int)y;
  int bits = 1;
  while(((size_t)1 << bits) < 2 * n) bits++;
  if(k == 0) {
    // Enough letters that a seed has few copies in b by chance: about one
    // for random letters
    k = 8;
    while(k < 15 && k * bits < 60 && ((size_t)1 << (bits * k)) < len_b) k++;
  }
  if(k < 2 || k * bits > 62 || len_a < k || len_b < k) return NULL;

  seeds_t *sd = aln_calloc(1, sizeof(*sd));
  if(sd == NULL) return NULL;
  sd->k = k;
  sd->l = smax;
  sd->D = (ptrdiff_t)len_a - (ptrdiff_t)len_b;
  sd->len_a = len_a;
  sd->len_b = len_b;
  sd->gap_extend = sc->gap_extend;
  sd->open = -2LL * sc->gap_open;
  sd->gamma = (int64_t)smax - 2LL * sc->gap_extend;
  // A run of r gapped letters of a breaks at most ceil((r-1)/k) + 1 seeds:
  // one letter one seed, two letters up to two, and the longest runs k
  // letters a seed, so the least cost per seed is at r = 1, r = 2 or in
  // the limit, whichever is least
  int64_t c = MIN2(2LL * (smax - mm), sd->open + sd->gamma);
  c = MIN2(c, (sd->open + 2 * sd->gamma) / 2);
  c = MIN2(c, (int64_t)k * sd->gamma);
  sd->c = c;
  sd->lip = MAX2(sd->open + sd->gamma, c);
  if(c <= 0) { free(sd); return NULL; }

  // b's words, sorted
  const size_t nb = len_b - k + 1;
  kmer_t *bk = aln_malloc(nb * sizeof(kmer_t));
  if(bk == NULL) { free(sd); return NULL; }
  size_t j; uint64_t key = 0;
  const uint64_t mask = (((uint64_t)1 << (bits * k - 1)) << 1) - 1;
  for(j = 0; j < len_b; j++) {
    key = ((key << bits) | (uint64_t)cls[t->code[(unsigned char)b[j]]]) & mask;
    if(j + 1 >= k) { bk[j+1-k].key = key; bk[j+1-k].pos = (uint32_t)(j+1-k); }
  }
  qsort(bk, nb, sizeof(kmer_t), kmer_cmp);

  // a's seeds and their copies
  const size_t smaxn = len_a / k;
  size_t mcap = smaxn * 2 + 16, M = 0, S = 0, s;
  sd->sa = aln_malloc(smaxn * sizeof(size_t));
  sd->mstart = aln_malloc((smaxn + 1) * sizeof(size_t));
  sd->mj = aln_malloc(mcap * sizeof(size_t));
  if(sd->sa == NULL || sd->mstart == NULL || sd->mj == NULL) {
    free(bk); seeds_free(sd); return NULL;
  }
  for(s = 0; s < smaxn; s++) {
    const size_t i = s * k;
    key = 0;
    for(x = 0; x < k; x++)
      key = (key << bits) | (uint64_t)cls[t->code[(unsigned char)a[i+x]]];
    size_t lo = 0, hi = nb;
    while(lo < hi) {
      const size_t mid = (lo + hi) / 2;
      if(bk[mid].key < key) lo = mid + 1; else hi = mid;
    }
    size_t cnt = 0;
    while(lo + cnt < nb && bk[lo+cnt].key == key) cnt++;
    if(cnt > SEEDS_MAX_COPIES) continue;
    if(M + cnt > mcap) {
      mcap = (M + cnt) * 2;
      size_t *mj = aln_realloc(sd->mj, mcap * sizeof(size_t));
      if(mj == NULL) { free(bk); seeds_free(sd); return NULL; }
      sd->mj = mj;
    }
    sd->sa[S] = i;
    sd->mstart[S] = M;
    for(x = 0; x < cnt; x++) sd->mj[M++] = bk[lo+x].pos;
    S++;
  }
  free(bk);
  sd->mstart[S] = M;
  sd->S = S;
  sd->M = M;

  // Chains from each copy, by decreasing seed
  sd->G = aln_malloc((M + 1) * sizeof(int64_t));
  sd->Gs = aln_malloc((S + 2) * sizeof(int64_t));
  if(sd->G == NULL || sd->Gs == NULL) { seeds_free(sd); return NULL; }
  sd->Gs[S] = sd->Gs[S+1] = SEEDS_NONE;
  for(s = S; s-- > 0; ) {
    int64_t gs = sd->Gs[s+1];
    size_t m;
    for(m = sd->mstart[s]; m < sd->mstart[s+1]; m++) {
      const ptrdiff_t dm = (ptrdiff_t)sd->sa[s] - (ptrdiff_t)sd->mj[m];
      int64_t g = MAX2(seeds_gap(sd, sd->D - dm), c * (int64_t)(S - s - 1));
      size_t s2;
      for(s2 = s + 1; s2 < S && s2 <= s + SEEDS_WINDOW; s2++) {
        const int64_t skip = c * (int64_t)(s2 - s - 1);
        size_t m2;
        for(m2 = sd->mstart[s2]; m2 < sd->mstart[s2+1]; m2++) {
          if(sd->mj[m2] < sd->mj[m] + k) continue; // before this copy ends
          const ptrdiff_t dm2 = (ptrdiff_t)sd->sa[s2] - (ptrdiff_t)sd->mj[m2];
          g = MIN2(g, MAX2(seeds_gap(sd, dm2 - dm), skip) + sd->G[m2]);
        }
      }
      if(s + SEEDS_WINDOW + 1 < S)
        g = MIN2(g, sd->Gs[s + SEEDS_WINDOW + 1] - c * (int64_t)(s + 1));
      sd->G[m] = g;
      gs = MIN2(gs, c * (int64_t)s + g);
    }
    sd->Gs[s] = gs;
  }
  return sd;
}

void seeds_free(seeds_t *sd)
{
  if(sd == NULL) return;
  free(sd->sa);
  free(sd->mstart);
  free(sd->mj);
  free(sd->G);
  free(sd->Gs);
  free(sd->last);
  free(sd);
}

// Least penalty of a chain from cell (x, y): from the first seed at or
// after x, copies at or after row y
static int64_t seeds_h(const seeds_t *sd, size_t x, size_t y)
{
  const ptrdiff_t dw = (ptrdiff_t)x - (ptrdiff_t)y;
  size_t lo = 0, hi = sd->S;
  while(lo < hi) {
    const size_t mid = (lo + hi) / 2;
    if(sd->sa[mid] < x) lo = mid + 1; else hi = mid;
  }
  const size_t s0 = lo;
  int64_t h = seeds_gap(sd, sd->D - dw);
  if(s0 >= sd->S) return h;
  h = MAX2(h, sd->c * (int64_t)(sd->S - s0));
  size_t s;
  for(s = s0; s < sd->S && s < s0 + SEEDS_WINDOW; s++) {
    const int64_t skip = sd->c * (int64_t)(s - s0);
    size_t m;
    for(m = sd->mstart[s]; m < sd->mstart[s+1]; m++) {
      if(sd->mj[m] < y) continue;
      const ptrdiff_t dm = (ptrdiff_t)sd->sa[s] - (ptrdiff_t)sd->mj[m];
      h = MIN2(h, MAX2(seeds_gap(sd, dm - dw), skip) + sd->G[m]);
    }
  }
  if(s0 + SEEDS_WINDOW < sd->S)
    h = MIN2(h, sd->Gs[s0 + SEEDS_WINDOW] - sd->c * (int64_t)s0);
  return h;
}

// The score bound from a penalty bound h at (x, y)
static int64_t seeds_score(const seeds_t *sd, size_t x, size_t y, int64_t h)
{
  const int64_t ri = (int64_t)(sd->len_a - x), rj = (int64_t)(sd->len_b - y);
  const int64_t plain = sd->l * MIN2(ri, rj) +
                        (int64_t)sd->gap_extend * (ri > rj ? ri - rj : rj - ri);
  h -= sd->open;
  const int64_t num = sd->l * (ri + rj) - h;
  const int64_t ours = num >= 0 ? num / 2 : -((-num + 1) / 2);
  return MIN2(ours, plain);
}

int64_t seeds_bound(const seeds_t *sd, size_t x, size_t y)
{
  return seeds_score(sd, x, y, seeds_h(sd, x, y));
}

// Penalty bounds along row y for columns from .. to, exact every
// SEEDS_SAMPLE columns (and at to) and, between, from the nearer samples:
// h moves by at most lip per column. Returns the samples (n of them, from
// column from) for keeping
static int64_t *seeds_row(const seeds_t *sd, size_t y, size_t from,
                          size_t to, size_t *n_out, int64_t *h)
{
  const size_t n = (to - from) / SEEDS_SAMPLE + 2;
  int64_t *samples = aln_malloc(n * sizeof(int64_t));
  if(samples == NULL) return NULL;
  size_t i, x;
  for(i = 0; i < n; i++) {
    x = MIN2(to, from + i * SEEDS_SAMPLE);
    samples[i] = seeds_h(sd, x, y);
  }
  for(x = from; x <= to; x++) {
    const size_t i = (x - from) / SEEDS_SAMPLE, xi = from + i * SEEDS_SAMPLE;
    int64_t v = samples[i] - sd->lip * (int64_t)(x - xi);
    if(i + 1 < n) {
      const size_t xn = MIN2(to, xi + SEEDS_SAMPLE);
      v = MAX2(v, samples[i+1] - sd->lip * (int64_t)(xn - x));
    }
    h[x - from] = v;
  }
  *n_out = n;
  return samples;
}

// Penalty bounds along the row the last call sampled, from its samples
static bool seeds_row_kept(const seeds_t *sd, size_t y, size_t from,
                           size_t to, int64_t *h)
{
  if(sd->last == NULL || sd->last_y != y || from < sd->last_from ||
     to > sd->last_to) return false;
  size_t x;
  for(x = from; x <= to; x++) {
    const size_t i = (x - sd->last_from) / SEEDS_SAMPLE;
    const size_t xi = sd->last_from + i * SEEDS_SAMPLE;
    int64_t v = sd->last[i] - sd->lip * (int64_t)(x - xi);
    if(i + 1 < sd->last_n) {
      const size_t xn = MIN2(sd->last_to, xi + SEEDS_SAMPLE);
      v = MAX2(v, sd->last[i+1] - sd->lip * (int64_t)(xn - x));
    }
    h[x - from] = v;
  }
  return true;
}

void seeds_narrow(seeds_t *sd, score_t *const row[3],
                  score_t *const col0[3], size_t y1, size_t c0, size_t c1,
                  size_t y2, score_t smax, int gap_extend, score_t lb,
                  ptrdiff_t *lo, ptrdiff_t *hi)
{
  if(gap_extend >= 0 || y2 <= y1 || *lo > *hi) return;
  const int64_t e = gap_extend, e2 = (int64_t)smax - 2 * e;
  const int64_t NEG = INT64_MIN / 4;
  const size_t dy = y2 - y1;
  // Row y1's cells on a best path, and how far below its diagonals a path
  // from one can stray (two gap letters per diagonal, a matched pair less)
  const size_t xlo = c0 == 0 ? 0 : c0 + 1, xhi = c1;
  const size_t nx = xhi - xlo + 1;
  int64_t *v = aln_malloc(nx * sizeof(int64_t));
  int64_t *h = aln_malloc(nx * sizeof(int64_t));
  if(v == NULL || h == NULL) { free(v); free(h); return; }
  if(!seeds_row_kept(sd, y1, xlo, xhi, h)) {
    size_t n;
    int64_t *samples = seeds_row(sd, y1, xlo, xhi, &n, h);
    if(samples == NULL) { free(v); free(h); return; }
    free(samples);
  }
  ptrdiff_t dmin = PTRDIFF_MAX, dmax = PTRDIFF_MIN;
  int64_t tmax = 0;
  size_t x;
  for(x = xlo; x <= xhi; x++) {
    const int64_t val = x == 0
      ? MAX3(col0[MATCH][y1], col0[GAP_A][y1], col0[GAP_B][y1])
      : MAX3(row[MATCH][x], row[GAP_A][x], row[GAP_B][x]);
    const int64_t bound = val + seeds_score(sd, x, y1, h[x - xlo]);
    if(bound >= lb) {
      v[x - xlo] = val;
      const ptrdiff_t d = (ptrdiff_t)x - (ptrdiff_t)y1;
      dmin = MIN2(dmin, d);
      dmax = MAX2(dmax, d);
      tmax = MAX2(tmax, (bound - lb) / e2);
    }
    else v[x - xlo] = NEG;
  }
  free(h);
  if(dmin > dmax) { free(v); return; }

  // Row y2's cells of the strip's diagonals: each scores at most the best
  // of row y1's cells plus dy rows down and dx columns right at most
  // smax min(dx, dy) + e |dx - dy|: for dx >= dy a running best of
  // v(x) - e x, under dy a sliding window of v(x) - (smax - e) x
  const ptrdiff_t ylo = (ptrdiff_t)y2 + *lo, yhi = (ptrdiff_t)y2 + *hi;
  const size_t x2lo = ylo < 0 ? 0 : (size_t)ylo;
  const size_t x2hi = yhi >= (ptrdiff_t)sd->len_a ? sd->len_a : (size_t)yhi;
  if(x2lo > x2hi) { free(v); return; }
  const size_t n2 = x2hi - x2lo + 1;
  int64_t *ub = aln_malloc(n2 * sizeof(int64_t));
  int64_t *h2 = aln_malloc(n2 * sizeof(int64_t));
  size_t *dq = aln_malloc(nx * sizeof(size_t));
  if(ub == NULL || h2 == NULL || dq == NULL) {
    free(v); free(ub); free(h2); free(dq);
    return;
  }
  const int64_t slope = (int64_t)smax - e;
  int64_t far = NEG;
  size_t xa = xlo, xb = xlo, qh = 0, qt = 0, i2;
  for(i2 = 0; i2 < n2; i2++) {
    const size_t xp = x2lo + i2;
    while(xa <= xhi && xa + dy <= xp) {
      if(v[xa - xlo] != NEG) far = MAX2(far, v[xa - xlo] - e * (int64_t)xa);
      xa++;
    }
    while(xb <= xhi && xb <= xp) {
      if(v[xb - xlo] != NEG) {
        const int64_t key = v[xb - xlo] - slope * (int64_t)xb;
        while(qt > qh &&
              v[dq[qt-1] - xlo] - slope * (int64_t)dq[qt-1] <= key) qt--;
        dq[qt++] = xb;
      }
      xb++;
    }
    while(qt > qh && dq[qh] + dy <= xp) qh++;
    int64_t best = NEG;
    if(far != NEG)
      best = far + e * ((int64_t)xp - (int64_t)dy) + (int64_t)smax * (int64_t)dy;
    if(qt > qh) {
      const size_t xq = dq[qh];
      const int64_t near = v[xq - xlo] + slope * ((int64_t)xp - (int64_t)xq) +
                           e * (int64_t)dy;
      best = MAX2(best, near);
    }
    ub[i2] = best;
  }
  free(dq);
  free(v);
  size_t n;
  int64_t *samples = seeds_row(sd, y2, x2lo, x2hi, &n, h2);
  if(samples == NULL) { free(ub); free(h2); return; }
  for(i2 = 0; i2 < n2; i2++) {
    if(ub[i2] == NEG) continue;
    const size_t xp = x2lo + i2;
    const int64_t bound = ub[i2] + seeds_score(sd, xp, y2, h2[i2]);
    if(bound >= lb) {
      const ptrdiff_t d = (ptrdiff_t)xp - (ptrdiff_t)y2;
      dmin = MIN2(dmin, d);
      dmax = MAX2(dmax, d);
      tmax = MAX2(tmax, (bound - lb) / e2);
    }
  }
  free(ub);
  free(h2);
  *lo = MAX2(*lo, dmin - (ptrdiff_t)tmax);
  *hi = MIN2(*hi, dmax + (ptrdiff_t)tmax);
  free(sd->last);
  sd->last = samples;
  sd->last_from = x2lo;
  sd->last_to = x2hi;
  sd->last_n = n;
  sd->last_y = y2;
}
