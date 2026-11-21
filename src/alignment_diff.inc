/*
 alignment_diff.inc
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// Private to alignment.c: the striped fill of Needleman-Wunsch traceback
// bytes with difference-encoded lanes (after Suzuki & Kasahara 2018, BMC
// Bioinformatics 19:45, as ksw2 uses them), included once for each lane
// width. Scores grow with the sequences, but the difference between a cell
// and the cell above it does not, so 8- or 16-bit lanes hold every cell of
// any length. Row y keeps, for each column j:
//   m, a, b   M, A and B of (y, j) less H(y-1, j), the cell above
//   oM, oA, oB  M, A and B of (y, j) less H(y, j): 0 for the best, else below
//   v         H(y, j) - H(y, j-1)
// with H the best of the three. Every traceback byte comes from comparing
// values with the same reference, so the bytes are exactly those of the
// 32-bit fill (alignment_kernel.inc) and the alignment is the same.
// Gap scores must not be above 0. The cells that stand for -infinity in
// the absolute fill (M and A of row 0, M and B of column 0) take p->floor:
// it is below every real difference by more than one step can raise it, so
// it never ties with or beats a real value (see alignment_diff_bits()).
// One block (alignment_block_t) fills rows y0 .. y1 of its columns, with
// trace set (not banded: no trace_rows or lastrow) and neither track_best
// nor score_only. It reads its left column and writes its last (rc) as
// absolute scores, so blocks of either fill can sit side by side.
// Before including, define:
//   ALN_D_SUFFIX  8 or 16: the alignment_kernel.inc instance whose vector
//                 type and helpers to use (int8_t or int16_t lanes)
//   ALN_D_S       its lane type
//   ALN_D_VL      its lanes per vector

#define ALN_D_CAT2(a, b) a##_##b
#define ALN_D_CAT(a, b) ALN_D_CAT2(a, b)
#define DN(name) ALN_D_CAT(name, ALN_D_SUFFIX)
#define DS ALN_D_S
#define DV DN(kv)

ALN_TARGET_CLONES
static void DN(alignment_fill_diff)(alignment_block_t *p)
{
  aligner_t *aligner = p->aligner;
  const sub_table_t *table = p->table;
  const scoring_t *scoring = aligner->scoring;
  const size_t c0 = p->c0, ls = p->lstride, ly0 = p->ly0;
  const score_t *lcM = p->lc[MATCH], *lcA = p->lc[GAP_A], *lcB = p->lc[GAP_B];
  const char *seq_a = aligner->seq_a + c0;
  const size_t w = aligner->score_width, len_i = p->width;
  const size_t seg = (len_i + ALN_D_VL - 1) / ALN_D_VL;
  const size_t n = table->n;
  const score_t lo = p->floor, up = p->ceiling;
  const score_t open = scoring->gap_open + scoring->gap_extend;
  size_t c, k, l, y, i;

  const DV vF = DN(ksplat)((DS)lo);
  const DV vopen = DN(ksplat)((DS)open);
  const DV vext = DN(ksplat)((DS)scoring->gap_extend);
  const DV vzero = DN(ksplat)(0);
  const DV codeA = DN(ksplat)(GAP_A), codeB = DN(ksplat)(GAP_B);
  const DV codeM = DN(ksplat)(MATCH);

  // profile[c*seg + k] lane l: score of seq_a[l*seg+k] against character c
  // of seq_b (0 past the end of seq_a); then the rows' vectors
  DV *buf = aln_aligned_alloc(sizeof(DV), (n + 12) * seg * sizeof(DV));
  if(buf == NULL) {
    atomic_store_explicit(p->nomem, true, memory_order_relaxed);
    atomic_store_explicit(p->stop, true, memory_order_relaxed);
    return;
  }
  DV *profile = buf;
  DV *pOM = buf + n*seg, *pOA = pOM + seg, *pOB = pOA + seg, *pV = pOB + seg;
  DV *cOM = pV + seg, *cOA = cOM + seg, *cOB = cOA + seg, *cV = cOB + seg;
  DV *cm = cV + seg, *ca = cm + seg, *cb = ca + seg, *tb = cb + seg;
  DV *tmp;

  for(c = 0; c < n; c++) {
    for(k = 0; k < seg; k++) {
      for(l = 0; l < ALN_D_VL; l++) {
        i = l*seg + k;
        profile[c*seg+k][l] = i < len_i
          ? (DS)table->scores[c*n + table->code[(unsigned char)seq_a[i]]]
          : 0;
      }
    }
  }

  // Row y0-1 from its absolute scores (r0); columns past the end of seq_a
  // start level (all 0), which keeps them in range and changes no real
  // cell (none lies to their right)
  // The first column's v is set each row, once the left column is ready
  // (it may be the left block's, still being written)
  score_t hl = 0; // H(y0-1) of the column on the left
  score_t hfirst = 0, hlast = 0; // H(y-1) of the block's first, last column
  // H(y-1) of each lane's last column (past the end of seq_a: of the last
  // real one, as those columns start level), for the carries between lanes
  score_t hlane[ALN_D_VL];
  for(i = 0; i < len_i; i++) {
    const size_t x = c0 + 1 + i;
    const score_t hm = p->r0[MATCH][x], ha = p->r0[GAP_A][x];
    const score_t hb = p->r0[GAP_B][x], h = MAX3(hm, ha, hb);
    k = i % seg;
    l = i / seg;
    pOM[k][l] = (DS)MAX2(hm - h, lo);
    pOA[k][l] = (DS)MAX2(ha - h, lo);
    pOB[k][l] = (DS)MAX2(hb - h, lo);
    pV[k][l] = (DS)(i == 0 ? 0 : h - hl);
    if(i == 0) hfirst = h;
    if(k == seg - 1) hlane[l] = h;
    hlast = h;
    hl = h;
  }
  for(l = (len_i - 1) / seg; l < ALN_D_VL; l++)
    if(l * seg + seg > len_i) hlane[l] = hlast;
  for(; i < seg * ALN_D_VL; i++) {
    k = i % seg;
    l = i / seg;
    pOM[k][l] = pOA[k][l] = pOB[k][l] = pV[k][l] = 0;
  }
  for(c = 0; c < 3; c++) p->rc[c][p->y0-1] = p->r0[c][c0+len_i];

  const size_t kl = (len_i - 1) % seg, ll = (len_i - 1) / seg; // last cell

  for(y = p->y0; y <= p->y1; y++)
  {
    const DV *prof = profile +
      (size_t)table->code[(unsigned char)aligner->seq_b[y-1]] * seg;

    // Column c0 (initialised, or the left block's) feeds lane 0 at k = 0
    if(p->left_done != NULL && !alignment_wait_row(p->left_done, y, p->stop))
      break;
    const size_t lp = (y-1-ly0)*ls, lc = (y-ly0)*ls;
    const score_t lm_prev = lcM[lp], la_prev = lcA[lp], lb_prev = lcB[lp];
    const score_t lm = lcM[lc], la = lcA[lc], lb = lcB[lc];
    const score_t hc0p = MAX3(lm_prev, la_prev, lb_prev); // H(y-1, c0)
    // Column c0 in row y less H(y-1, c0), and its offsets in row y-1
    const score_t em0 = MAX2(lm - hc0p, lo), ea0 = MAX2(la - hc0p, lo);
    const score_t eb0 = MAX2(lb - hc0p, lo);
    const score_t h0 = MAX3(lm, la, lb) - hc0p;
    const score_t da0 = MAX2(la_prev - hc0p, lo);
    const score_t db0 = MAX2(lb_prev - hc0p, lo);
    pV[0][0] = (DS)(hfirst - hc0p); // H(y-1, c0+1) - H(y-1, c0)

    // Match from the diagonal and gap_a from above, both less the cell
    // above; gap_b from the left within each lane, less the cell above:
    // the left cell's own reference is one column back, so v comes off.
    // Lanes other than 0 start without what the lane before carries over
    // (-infinity), which the lazy F loop adds
    DV vx = DN(kshift)(vF, (DS)MAX2(em0, ea0));
    DV vb = DN(kshift)(vF, (DS)eb0);
    for(k = 0; k < seg; k++) {
      cm[k] = DN(kmax)(prof[k] - pV[k], vF);
      ca[k] = DN(kmax)(DN(kmax)(pOA[k] + vext,
                                DN(kmax)(pOM[k], pOB[k]) + vopen), vF);
      const DV t = DN(kmax)(DN(kmax)(vx + vopen, vb + vext), vF);
      cb[k] = DN(kmax)(t - pV[k], vF);
      vx = DN(kmax)(cm[k], ca[k]);
      vb = cb[k];
    }

    // Lazy F with the carries found first. What lane l passes to lane l+1
    // is the better of its own (from its last column) and what it was
    // passed, less seg gap extensions (one per column it crosses): a scan
    // over the lanes, in absolute scores. Each lane then takes its carry,
    // as its first column's gap_b from the left, along its columns while
    // it improves a cell. A carry at the floor (none, or so low it lost)
    // never wins: the lane's own gap_b starts above the floor and each
    // column takes the same off both (lane 0's comes from column c0, a
    // real cell, in the pass above)
    const long long none = LLONG_MIN / 4;
    long long cin = none;
    DV vf = vF;
    for(l = 0; l + 1 < ALN_D_VL; l++) {
      long long out = (long long)vx[l] + open;
      if(vb[l] != (DS)lo) out = MAX2(out, (long long)vb[l] + scoring->gap_extend);
      out += hlane[l];
      cin = cin == none ? out
                        : MAX2(out, cin + (long long)seg * scoring->gap_extend);
      // Relative to H(y-1) of lane l+1's first column (any real carry is
      // at most a cell above it, so within the lanes)
      long long rel = cin - hlane[l] - pV[0][l+1];
      vf[l+1] = (DS)MIN2(MAX2(rel, (long long)lo), (long long)up);
    }
    k = 0;
    while(DN(kany)(vf > cb[k])) {
      cb[k] = DN(kmax)(cb[k], vf);
      if(++k == seg) break;
      vf = DN(kmax)(DN(kmax)(vf + vext, vF) - pV[k], vF);
    }

    // Each cell's best, offsets, v, and traceback bytes from its
    // predecessors: diagonal (match), above (gap_a) and left (gap_b)
    const DV hend = DN(kmax)(DN(kmax)(cm[seg-1], ca[seg-1]), cb[seg-1]);
    DV hleft = DN(kshift)(hend, (DS)h0);
    DV dA = DN(kshift)(pOA[seg-1], (DS)da0);
    DV dB = DN(kshift)(pOB[seg-1], (DS)db0);
    DV eM = DN(kshift)(cm[seg-1], (DS)em0);
    DV eA = DN(kshift)(ca[seg-1], (DS)ea0);
    DV eB = DN(kshift)(cb[seg-1], (DS)eb0);
    for(k = 0; k < seg; k++) {
      const DV h = DN(kmax)(DN(kmax)(cm[k], ca[k]), cb[k]);
      // The diagonal cell's best state, gap_a before gap_b before match
      DV fM = DN(ksel)(dB == vzero, codeB, codeM);
      fM = DN(ksel)(dA == vzero, codeA, fM);
      // gap_b against values one column back: add v back
      const DV bc = cb[k] + pV[k];
      tb[k] = fM |
              DN(ktrace)(ca[k], pOA[k] + vext, pOB[k] + vopen,
                         pOM[k] + vopen) << 2 |
              DN(ktrace)(bc, eA + vopen, eB + vext, eM + vopen) << 4;
      cV[k] = pV[k] + h - hleft;
      cOM[k] = DN(kmax)(cm[k] - h, vF);
      cOA[k] = DN(kmax)(ca[k] - h, vF);
      cOB[k] = DN(kmax)(cb[k] - h, vF);
      hleft = h;
      dA = pOA[k]; dB = pOB[k];
      eM = cm[k]; eA = ca[k]; eB = cb[k];
    }

    uint8_t *trow = p->trace + y*w + c0 + 1;
    for(l = 0; l < ALN_D_VL; l++)
      for(k = 0, i = l*seg; k < seg && i < len_i; k++, i++)
        trow[i] = (uint8_t)tb[k][l];

    // The last column's scores for the block on the right (and the final
    // cell), and the first and last columns' H for the next row
    p->rc[MATCH][y] = hlast + cm[kl][ll];
    p->rc[GAP_A][y] = hlast + ca[kl][ll];
    p->rc[GAP_B][y] = hlast + cb[kl][ll];
    hlast += MAX3((score_t)cm[kl][ll], (score_t)ca[kl][ll],
                  (score_t)cb[kl][ll]);
    hfirst += MAX3((score_t)cm[0][0], (score_t)ca[0][0], (score_t)cb[0][0]);
    for(l = 0; l < ALN_D_VL; l++) hlane[l] += hend[l];

    if(p->done != NULL) atomic_store_explicit(p->done, y, memory_order_release);

    tmp = pOM; pOM = cOM; cOM = tmp;
    tmp = pOA; pOA = cOA; cOA = tmp;
    tmp = pOB; pOB = cOB; cOB = tmp;
    tmp = pV; pV = cV; cV = tmp;
  }

  free(buf);
}

#undef DV
#undef DS
#undef DN
#undef ALN_D_CAT
#undef ALN_D_CAT2
#undef ALN_D_SUFFIX
#undef ALN_D_S
#undef ALN_D_VL
