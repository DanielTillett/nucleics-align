/*
 alignment_kernel.inc
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// Private to alignment.c: the striped fill, included once for each lane
// width. One block (alignment_block_t) fills rows y0 .. y1 of its columns
// and, by which fields are set, gives one of these outputs:
//   full matrices  trace == NULL, !score_only: the three score matrices
//   traceback      trace != NULL: a byte per cell (trace_rows: banded, by
//                  strip), the last column (rc) and, with track_best, the
//                  first best match cell; lastrow for the strip below
//   score scan     score_only: the last column and the best cell
//                  (track_best) or the cells with score target (target > 0)
// Before including, define:
//   ALN_K_SUFFIX  appended to the names defined here (e.g. 32)
//   ALN_K_S       lane type (score_t, or int16_t)
//   ALN_K_VL      lanes per vector (a plain number or a sum of them, for #if)
//   ALN_K_NARROW  1 if lanes are narrower than score_t: rows stop, setting
//                 *p->overflow, once a row's best match score passes
//                 p->narrow_limit (see alignment_narrow_limit())
// alignment_block_t, ALN_INLINE and the vector helpers must already be
// defined. The macros are undefined at the end, ready for the next width

#define ALN_K_CAT2(a, b) a##_##b
#define ALN_K_CAT(a, b) ALN_K_CAT2(a, b)
#define KN(name) ALN_K_CAT(name, ALN_K_SUFFIX)
#define KS ALN_K_S
#define KV KN(kv)

typedef KS KV __attribute__((vector_size(ALN_K_VL * sizeof(KS))));

ALN_INLINE KV KN(ksplat)(KS x) { return (KV){0} + x; }

ALN_INLINE KV KN(kmax)(KV a, KV b)
{
  KV gt = a > b; // all bits set in lanes where a > b
  return (a & gt) | (b & ~gt);
}

// Lanes of a where mask is set, of b elsewhere
ALN_INLINE KV KN(ksel)(KV mask, KV a, KV b)
{
  return (a & mask) | (b & ~mask);
}

// Move every lane up by one; lane 0 gets x
ALN_INLINE KV KN(kshift)(KV v, KS x)
{
#if ALN_K_VL == 32
  return __builtin_shufflevector(v, KN(ksplat)(x), 32, 0, 1, 2, 3, 4, 5, 6,
                                 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                                 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                 30);
#elif ALN_K_VL == 16
  return __builtin_shufflevector(v, KN(ksplat)(x), 16, 0, 1, 2, 3, 4, 5, 6,
                                 7, 8, 9, 10, 11, 12, 13, 14);
#elif ALN_K_VL == 8
  return __builtin_shufflevector(v, KN(ksplat)(x), 8, 0, 1, 2, 3, 4, 5, 6);
#else
  return __builtin_shufflevector(v, KN(ksplat)(x), 4, 0, 1, 2);
#endif
}

ALN_INLINE bool KN(kany)(KV v)
{
  KS r = 0;
  int l;
  for(l = 0; l < ALN_K_VL; l++) r |= v[l];
  return r != 0;
}

// Matrix the traceback steps to from a cell scoring cur, given each
// predecessor's score plus the penalty for stepping from it: the first of
// GAP_A, GAP_B and MATCH that gives cur, as alignment_reverse_move() checks
ALN_INLINE KV KN(ktrace)(KV cur, KV from_a, KV from_b, KV from_m)
{
  KV c = KN(ksplat)(ALN_TRACE_NONE);
  c = KN(ksel)(from_m == cur, KN(ksplat)(MATCH), c);
  c = KN(ksel)(from_b == cur, KN(ksplat)(GAP_B), c);
  return KN(ksel)(from_a == cur, KN(ksplat)(GAP_A), c);
}

ALN_TARGET_CLONES
static void KN(alignment_fill_striped)(alignment_block_t *p)
{
  aligner_t *aligner = p->aligner;
  const sub_table_t *table = p->table;
  const scoring_t *scoring = aligner->scoring;
  const score_t min = p->min;
  const size_t c0 = p->c0, ls = p->lstride;
  const bool full = p->trace == NULL && !p->score_only;
  score_t *match_scores = full ? aligner->match_scores + c0 : NULL;
  score_t *gap_a_scores = full ? aligner->gap_a_scores + c0 : NULL;
  score_t *gap_b_scores = full ? aligner->gap_b_scores + c0 : NULL;
  const score_t *lcM = p->lc[MATCH], *lcA = p->lc[GAP_A], *lcB = p->lc[GAP_B];
  const char *seq_a = aligner->seq_a + c0;
  const size_t w = aligner->score_width;
  const size_t len_i = p->width, ly0 = p->ly0;
  const size_t seg = (len_i + ALN_K_VL - 1) / ALN_K_VL;
  const size_t n = table->n;
  size_t c, k, l, y, i;

  const KV vmin = KN(ksplat)((KS)min);
  const KV vopen = KN(ksplat)((KS)(scoring->gap_open + scoring->gap_extend));
  const KV vext = KN(ksplat)((KS)scoring->gap_extend);

  // profile[c*seg + k] lane l: score of seq_a[l*seg+k] against character c
  // of seq_b (0 past the end of seq_a; those cells are never read back)
  size_t nvec = (n + 7) * seg;
  KV *buf = aligned_alloc(sizeof(KV), nvec * sizeof(KV));
  if(buf == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  KV *profile = buf;
  KV *pM = buf + n*seg, *pA = pM + seg, *pB = pA + seg; // previous row
  KV *cM = pB + seg, *cA = cM + seg, *cB = cA + seg; // current row
  KV *tb = cB + seg; // traceback bytes of the current row
  KV *tmp;

  for(c = 0; c < n; c++) {
    for(k = 0; k < seg; k++) {
      for(l = 0; l < ALN_K_VL; l++) {
        i = l*seg + k;
        profile[c*seg+k][l] = i < len_i
          ? (KS)table->scores[c*n + table->code[(unsigned char)seq_a[i]]]
          : 0;
      }
    }
  }

  // Row y0-1 (row 0 for a whole fill) in striped order
  for(k = 0; k < seg; k++) {
    for(l = 0; l < ALN_K_VL; l++) {
      i = l*seg + k;
      pM[k][l] = (KS)(i < len_i ? p->r0[MATCH][c0+1+i] : min);
      pA[k][l] = (KS)(i < len_i ? p->r0[GAP_A][c0+1+i] : min);
      pB[k][l] = (KS)(i < len_i ? p->r0[GAP_B][c0+1+i] : min);
    }
  }
  if(!full)
    for(c = 0; c < 3; c++) p->rc[c][p->y0-1] = p->r0[c][c0+len_i];

  // Segments of each lane that are real cells (lanes past the end of seq_a
  // are left out when looking for the best cell)
  KV vlim;
  for(l = 0; l < ALN_K_VL; l++)
    vlim[l] = (KS)(l*seg >= len_i ? 0 : MIN2(seg, len_i - l*seg));

  for(y = p->y0; y <= p->y1; y++)
  {
    const size_t row = y*w;
    const KV *prof = profile +
      (size_t)table->code[(unsigned char)aligner->seq_b[y-1]] * seg;

    // Column c0 (initialised, or the left block's) feeds lane 0 at k = 0
    if(p->left_done != NULL && !alignment_wait_row(p->left_done, y, p->overflow))
      break;
    const size_t lp = (y-1-ly0)*ls, lc = (y-ly0)*ls;
    const score_t lm_prev = lcM[lp], la_prev = lcA[lp], lb_prev = lcB[lp];
    const score_t lm = lcM[lc], la = lcA[lc], lb = lcB[lc];
    KS h0_prev = (KS)MAX3(lm_prev, la_prev, lb_prev);
    KS x0 = (KS)MAX2(lm, la);
    KS b0 = (KS)lb; // min in the real column 0

    // One pass: match from the diagonal, gap_a from the cell above, gap_b
    // from the cell on the left within each lane (lanes other than 0 start
    // without what the lane before carries over, which the lazy F loop
    // adds) and, where wanted, each lane's best match cell (its first, as
    // only a higher score replaces it)
    const bool track = p->track_best || p->target > 0;
    KV vdiag = KN(kshift)(KN(kmax)(KN(kmax)(pM[seg-1], pA[seg-1]), pB[seg-1]),
                          h0_prev);
    KV vx = KN(kshift)(vmin, x0);
    KV vb = KN(kshift)(vmin, b0);
    KV vbest = vmin, vbk = KN(ksplat)(0);
    for(k = 0; k < seg; k++) {
      cM[k] = KN(kmax)(vdiag + prof[k], vmin);
      cA[k] = KN(kmax)(KN(kmax)(KN(kmax)(pM[k], pB[k]) + vopen, pA[k] + vext),
                       vmin);
      vdiag = KN(kmax)(KN(kmax)(pM[k], pA[k]), pB[k]);
      cB[k] = KN(kmax)(KN(kmax)(vx + vopen, vb + vext), vmin);
      vx = KN(kmax)(cM[k], cA[k]);
      vb = cB[k];
      if(track) {
        const KV vk = KN(ksplat)((KS)k);
        const KV gt = (cM[k] > vbest) & (vk < vlim);
        vbest = KN(ksel)(gt, cM[k], vbest);
        vbk = KN(ksel)(gt, vk, vbk);
      }
    }

    // Lazy F: carry gap_b across segment boundaries while it improves a cell,
    // starting with what each lane's last cell passes to the next lane
    KV vf = KN(kmax)(KN(kmax)(KN(kshift)(vx, (KS)min) + vopen,
                              KN(kshift)(vb, (KS)min) + vext), vmin);
    k = 0;
    while(KN(kany)(vf > cB[k])) {
      cB[k] = KN(kmax)(cB[k], vf);
      vf = KN(kmax)(vf + vext, vmin);
      if(++k == seg) { k = 0; vf = KN(kmax)(KN(kshift)(vf, (KS)min), vmin); }
    }

    if(full)
    {
      // Back to row-major for traceback
      for(l = 0; l < ALN_K_VL; l++) {
        for(k = 0, i = l*seg; k < seg && i < len_i; k++, i++) {
          match_scores[row+1+i] = cM[k][l];
          gap_a_scores[row+1+i] = cA[k][l];
          gap_b_scores[row+1+i] = cB[k][l];
        }
      }
    }
    else
    {
      if(!p->score_only)
      {
        // Traceback bytes, from each cell's predecessors: diagonal (match),
        // above (gap_a) and left (gap_b)
        KV dM = KN(kshift)(pM[seg-1], (KS)lm_prev);
        KV dA = KN(kshift)(pA[seg-1], (KS)la_prev);
        KV dB = KN(kshift)(pB[seg-1], (KS)lb_prev);
        KV eM = KN(kshift)(cM[seg-1], (KS)lm);
        KV eA = KN(kshift)(cA[seg-1], (KS)la);
        KV eB = KN(kshift)(cB[seg-1], (KS)lb);
        for(k = 0; k < seg; k++) {
          tb[k] = KN(ktrace)(cM[k], dA + prof[k], dB + prof[k], dM + prof[k]) |
                  KN(ktrace)(cA[k], pA[k] + vext, pB[k] + vopen,
                             pM[k] + vopen) << 2 |
                  KN(ktrace)(cB[k], eA + vopen, eB + vext, eM + vopen) << 4;
          dM = pM[k]; dA = pA[k]; dB = pB[k];
          eM = cM[k]; eA = cA[k]; eB = cB[k];
        }
        uint8_t *trow = p->trace + (p->trace_rows != NULL
                                    ? (size_t)(p->trace_rows[y] +
                                               (ptrdiff_t)(c0 + 1))
                                    : row + c0 + 1);
        for(l = 0; l < ALN_K_VL; l++)
          for(k = 0, i = l*seg; k < seg && i < len_i; k++, i++)
            trow[i] = (uint8_t)tb[k][l];
      }

      i = len_i - 1;
      p->rc[MATCH][y] = cM[i % seg][i / seg];
      p->rc[GAP_A][y] = cA[i % seg][i / seg];
      p->rc[GAP_B][y] = cB[i % seg][i / seg];

      if(track)
      {
        score_t rmax = vbest[0];
        for(l = 1; l < ALN_K_VL; l++) rmax = MAX2(rmax, (score_t)vbest[l]);
#if ALN_K_NARROW
        // The next row could pass the lane type's range: stop, to be filled
        // again with wider lanes. Stop early too where the best so far,
        // rising as it has, would pass it by the last row (it rises about
        // steadily along a hit spanning the sequences), checking from a
        // sixteenth of the rows on: an early wrong guess only costs time
        const size_t rows = p->y1 - p->y0 + 1, done_rows = y - p->y0 + 1;
        if(rmax > p->narrow_limit ||
           (p->track_best && done_rows * 16 >= rows &&
            (double)MAX2(rmax, p->best_score) * (double)rows >
            (double)p->narrow_limit * (double)done_rows)) {
          atomic_store_explicit(p->overflow, true, memory_order_relaxed);
          break;
        }
#endif
        if(p->track_best && rmax > 0 && rmax >= p->best_score) {
          // Lowest column in this row with that score: the first lane whose
          // best has it, at that lane's first cell with it
          for(l = 0; (score_t)vbest[l] != rmax; l++) {}
          i = l*seg + (size_t)vbk[l];
          if(rmax > p->best_score || c0+1+i < p->best_x) {
            p->best_score = rmax; p->best_x = c0+1+i; p->best_y = y;
          }
        }
        if(p->target > 0 && rmax >= p->target) {
          // Highest column in this row with the target score (none score
          // more, see alignment_sw_scan())
          for(l = ALN_K_VL, i = len_i; l-- > 0 && i == len_i; )
            for(k = seg; k-- > 0; )
              if(l*seg+k < len_i && cM[k][l] == p->target) {
                i = l*seg+k;
                break;
              }
          if(i < len_i) {
            p->tgt_x = MAX2(p->tgt_x, c0+1+i);
            p->tgt_y = y;
          }
          p->best_score = MAX2(p->best_score, rmax);
        }
      }
    }
    if(p->done != NULL) atomic_store_explicit(p->done, y, memory_order_release);

    tmp = pM; pM = cM; cM = tmp;
    tmp = pA; pA = cA; cA = tmp;
    tmp = pB; pB = cB; cB = tmp;
  }

  // The last row, for the strip below
  if(p->lastrow[MATCH] != NULL) {
    for(l = 0; l < ALN_K_VL; l++) {
      for(k = 0, i = l*seg; k < seg && i < len_i; k++, i++) {
        p->lastrow[MATCH][c0+1+i] = pM[k][l];
        p->lastrow[GAP_A][c0+1+i] = pA[k][l];
        p->lastrow[GAP_B][c0+1+i] = pB[k][l];
      }
    }
  }

  free(buf);
}

#undef KV
#undef KS
#undef KN
#undef ALN_K_CAT
#undef ALN_K_CAT2
#undef ALN_K_SUFFIX
#undef ALN_K_S
#undef ALN_K_VL
#undef ALN_K_NARROW
