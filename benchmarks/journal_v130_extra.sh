#!/usr/bin/env bash
# Supplementary journal-grade benchmarks for Nucleics Align v1.3.0:
#   A. In-process C protein benchmark (Nucleics vs parasail scan & striped, zero startup overhead)
#   B. parasail AVX2 DNA global traceback (20 kb & 50 kb, 5 seeds)
#   C. Banded ksw2 (at Nucleics certified band width, 20 kb & 50 kb 1%, 5 seeds)
#   D. parasail 32-bit exact local + Nucleics --fastlocal (20 kb 5 seeds + 197 kb Mpox)
#   E. Previously skipped 50 kb cells (WFA2 BiWFA rand, Nucleics --wfa 10% & rand)
#   F. Thread scaling curve (1, 2, 4, 8, 16 threads on 20 kb 10%, 50 kb 10%, 50 kb rand)
#   G. Long-indel & unequal-length DNA workload (20 kb vs 18.5 kb, 50-200 bp indels, 3 seeds)
#
# Usage: bash benchmarks/journal_v130_extra.sh [work_dir] [threads]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:-/tmp/nucleics_compare}"
THREADS="${2:-8}"
mkdir -p "$WORK/src" "$WORK/bin" "$WORK/long_indel"

NW="$ROOT/bin/needleman_wunsch"
SW="$ROOT/bin/smith_waterman"
WFA2="$WORK/bin/run_wfa2"
KSW2="$WORK/bin/run_ksw2"

echo "# host: $(uname -a)"
echo "# cpu:  $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')"
echo "# governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)"
echo "# smt:  $(cat /sys/devices/system/cpu/smt/active 2>/dev/null || echo unknown)"
echo "# gcc:  $(gcc --version | head -n1)"
echo "# nucleics-align: $(git -C "$ROOT" rev-parse --short HEAD)"
echo "# nw_sha256: $(sha256sum "$NW" | awk '{print $1}')"
echo "# sw_sha256: $(sha256sum "$SW" | awk '{print $1}')"

# Ensure parasail config.h exists
cat > "$WORK/src/config.h" <<'EOF'
#define HAVE_POSIX_MEMALIGN 1
#define HAVE_AVX2 1
#define HAVE_SSE41 1
#define HAVE_SSE2 1
EOF

# --- Build A: In-process protein benchmark (Nucleics Align C API vs parasail scan & striped) ---
cat > "$WORK/src/inproc_prot.c" <<'EOF'
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "needleman_wunsch.h"
#include "hirschberg.h"
#include "parasail.h"
#include "parasail/matrices/blosum62.h"

static char *read_fa(const char *path, size_t *out_len) {
  FILE *f = fopen(path, "r");
  if (!f) { perror(path); exit(1); }
  size_t cap = 4096, len = 0;
  char *buf = malloc(cap);
  char line[8192];
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '>') continue;
    for (char *p = line; *p; p++) {
      if (*p <= ' ') continue;
      if (len + 2 > cap) { cap *= 2; buf = realloc(buf, cap); }
      buf[len++] = *p;
    }
  }
  fclose(f);
  buf[len] = '\0';
  *out_len = len;
  return buf;
}

static double now_ms(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

static int cmp_dbl(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

static void summarize(const char *tag, double *vals, int n, int score) {
  qsort(vals, (size_t)n, sizeof(double), cmp_dbl);
  double med = vals[n / 2];
  printf("  %-22s med=%7.3f ms [%6.3f-%6.3f]  score=%d\n",
         tag, med, vals[0], vals[n - 1], score);
}

int main(int argc, char **argv) {
  if (argc < 4) return 1;
  const char *name = argv[1];
  size_t la = 0, lb = 0;
  char *a = read_fa(argv[2], &la);
  char *b = read_fa(argv[3], &lb);
  int reps = (argc >= 5) ? atoi(argv[4]) : 15;
  double *t = malloc((size_t)reps * sizeof(double));

  printf("-- %s (%zu / %zu aa, %d in-process reps) --\n", name, la, lb, reps);

  /* 1. parasail scan (1t) */
  int sc_ps = 0;
  for (int r = -1; r < reps; r++) {
    double t0 = now_ms();
    parasail_result_t *res = parasail_nw_trace_scan_avx2_256_16(
        a, (int)la, b, (int)lb, 11, 1, &parasail_blosum62);
    parasail_cigar_t *cig = parasail_result_get_cigar(res, a, (int)la, b, (int)lb, &parasail_blosum62);
    char *cstr = parasail_cigar_decode(cig);
    sc_ps = res->score;
    free(cstr);
    parasail_cigar_free(cig);
    parasail_result_free(res);
    double dt = now_ms() - t0;
    if (r >= 0) t[r] = dt;
  }
  summarize("parasail_scan_1t", t, reps, sc_ps);

  /* 2. parasail striped (1t) */
  int sc_pst = 0;
  for (int r = -1; r < reps; r++) {
    double t0 = now_ms();
    parasail_result_t *res = parasail_nw_trace_striped_avx2_256_16(
        a, (int)la, b, (int)lb, 11, 1, &parasail_blosum62);
    parasail_cigar_t *cig = parasail_result_get_cigar(res, a, (int)la, b, (int)lb, &parasail_blosum62);
    char *cstr = parasail_cigar_decode(cig);
    sc_pst = res->score;
    free(cstr);
    parasail_cigar_free(cig);
    parasail_result_free(res);
    double dt = now_ms() - t0;
    if (r >= 0) t[r] = dt;
  }
  if (sc_pst != sc_ps) { fprintf(stderr, "SCORE MISMATCH parasail_striped\n"); return 2; }
  summarize("parasail_striped_1t", t, reps, sc_pst);

  scoring_t scoring;
  scoring_system_BLOSUM62(&scoring);
  alignment_options_t o1, o8;
  alignment_options_default(&o1); o1.max_threads = 1; o1.hb_max_threads = 1;
  alignment_options_default(&o8); o8.max_threads = 8; o8.hb_max_threads = 8;
  aln_error_t err;

  /* 3. Nucleics full (1t) */
  int sc_na = 0;
  for (int r = -1; r < reps; r++) {
    nw_aligner_t *nw = needleman_wunsch_new();
    nw->compact = true;
    nw->options = &o1;
    alignment_t *aln = alignment_create(la + lb + 16);
    double t0 = now_ms();
    if (needleman_wunsch_align_checked(a, b, la, lb, &scoring, nw, aln, &err) != ALN_OK) return 3;
    sc_na = aln->score;
    double dt = now_ms() - t0;
    alignment_free(aln);
    needleman_wunsch_free(nw);
    if (r >= 0) t[r] = dt;
  }
  if (sc_na != sc_ps) { fprintf(stderr, "SCORE MISMATCH na_full_1t\n"); return 2; }
  summarize("na_full_1t", t, reps, sc_na);

  /* 4. Nucleics --linear (1t) */
  int sc_lin1 = 0;
  for (int r = -1; r < reps; r++) {
    alignment_t *aln = alignment_create(la + lb + 16);
    double t0 = now_ms();
    if (hirschberg_align_checked(a, b, la, lb, &scoring, &o1, aln, &err) != ALN_OK) return 3;
    sc_lin1 = aln->score;
    double dt = now_ms() - t0;
    alignment_free(aln);
    if (r >= 0) t[r] = dt;
  }
  if (sc_lin1 != sc_ps) { fprintf(stderr, "SCORE MISMATCH na_linear_1t\n"); return 2; }
  summarize("na_linear_1t", t, reps, sc_lin1);

  /* 5. Nucleics --linear (8t) */
  int sc_lin8 = 0;
  for (int r = -1; r < reps; r++) {
    alignment_t *aln = alignment_create(la + lb + 16);
    double t0 = now_ms();
    if (hirschberg_align_checked(a, b, la, lb, &scoring, &o8, aln, &err) != ALN_OK) return 3;
    sc_lin8 = aln->score;
    double dt = now_ms() - t0;
    alignment_free(aln);
    if (r >= 0) t[r] = dt;
  }
  if (sc_lin8 != sc_ps) { fprintf(stderr, "SCORE MISMATCH na_linear_8t\n"); return 2; }
  summarize("na_linear_8t", t, reps, sc_lin8);

  free(t); free(a); free(b);
  return 0;
}
EOF

gcc -O3 -mavx2 -std=c11 -D_POSIX_C_SOURCE=200809L \
  -I "$WORK/src" -I "$WORK/src/parasail" -I "$WORK/src/parasail/src" \
  -I "$ROOT/src" -I "$ROOT/libs" \
  "$WORK/src/inproc_prot.c" \
  "$WORK/src/parasail/src/nw_trace_scan_avx2_256_16.c" \
  "$WORK/src/parasail/src/nw_trace_striped_avx2_256_16.c" \
  "$WORK/src/parasail/src/memory.c" \
  "$WORK/src/parasail/src/memory_avx2.c" \
  "$WORK/src/parasail/src/cigar.c" \
  "$WORK/src/parasail/src/matrix_lookup.c" \
  "$ROOT/src/libalign.a" "$ROOT/libs/string_buffer/libstrbuf.a" \
  -lpthread -lz -o "$WORK/bin/inproc_prot"

# --- Build B: parasail AVX2 32-bit global DNA (scan) & 32-bit local DNA (striped) ---
cat > "$WORK/src/run_parasail_dna.c" <<'EOF'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parasail.h"

static char *read_fa_pair(FILE *f, int *out_len) {
  size_t cap = 4096, len = 0; char *buf = malloc(cap); char line[8192]; int started = 0;
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '>') { if (started) { fseek(f, -(long)strlen(line), SEEK_CUR); break; } started = 1; continue; }
    started = 1;
    for (char *p = line; *p; p++) {
      if (*p <= ' ') continue;
      if (len + 2 > cap) { cap *= 2; buf = realloc(buf, cap); }
      buf[len++] = *p;
    }
  }
  buf[len] = '\0'; *out_len = (int)len; return buf;
}

int main(int argc, char **argv) {
  if (argc < 3) return 1;
  const char *mode = argv[1];
  FILE *f = fopen(argv[2], "r");
  if (!f) { perror(argv[2]); return 1; }
  int la = 0, lb = 0;
  char *a = read_fa_pair(f, &la), *b = read_fa_pair(f, &lb);
  fclose(f);
  parasail_matrix_t *mat;
  parasail_result_t *res;
  if (strcmp(mode, "nw") == 0) {
    /* +1 / -2 / -4 / -1 -> parasail open=5, extend=1 */
    mat = parasail_matrix_create("ACGTN", 1, -2);
    res = parasail_nw_trace_scan_avx2_256_32(a, la, b, lb, 5, 1, mat);
  } else {
    /* local +2 / -2 / -2 / -1 -> parasail open=3, extend=1 */
    mat = parasail_matrix_create("ACGTN", 2, -2);
    res = parasail_sw_trace_striped_avx2_256_32(a, la, b, lb, 3, 1, mat);
  }
  if (!res || parasail_result_is_saturated(res)) {
    fprintf(stderr, "parasail saturated or failed\n"); return 2;
  }
  parasail_cigar_t *cig = parasail_result_get_cigar(res, a, la, b, lb, mat);
  char *cstr = parasail_cigar_decode(cig);
  printf("score: %d\ncigar_len: %zu\n", res->score, strlen(cstr));
  free(cstr);
  parasail_cigar_free(cig);
  parasail_result_free(res);
  parasail_matrix_free(mat);
  free(a); free(b);
  return 0;
}
EOF

gcc -O3 -mavx2 -I "$WORK/src" \
  -I "$WORK/src/parasail" -I "$WORK/src/parasail/src" \
  "$WORK/src/run_parasail_dna.c" \
  "$WORK/src/parasail/src/nw_trace_scan_avx2_256_32.c" \
  "$WORK/src/parasail/src/sw_trace_striped_avx2_256_32.c" \
  "$WORK/src/parasail/src/memory.c" \
  "$WORK/src/parasail/src/memory_avx2.c" \
  "$WORK/src/parasail/src/cigar.c" \
  "$WORK/src/parasail/src/matrix_lookup.c" \
  -o "$WORK/bin/run_parasail_dna"

# --- Build C: Banded ksw2 wrapper (accepts explicit band width w) ---
cat > "$WORK/src/run_ksw2_banded.c" <<'EOF'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "ksw2.h"
static uint8_t *read_fa_enc(FILE *f, int *out_len) {
  size_t cap = 4096, len = 0; uint8_t *buf = malloc(cap); char line[8192]; int started = 0;
  static const uint8_t tbl[256] = {['A']=0,['a']=0,['C']=1,['c']=1,['G']=2,['g']=2,['T']=3,['t']=3};
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '>') { if (started) { fseek(f, -(long)strlen(line), SEEK_CUR); break; } started = 1; continue; }
    started = 1;
    for (char *p = line; *p; p++) { if (*p <= ' ') continue; if (len+1>cap){cap*=2;buf=realloc(buf,cap);} buf[len++]=tbl[(unsigned char)*p]; }
  }
  *out_len = (int)len; return buf;
}
int main(int argc, char **argv) {
  if (argc < 3) return 1;
  FILE *f = fopen(argv[1], "r");
  int w = atoi(argv[2]);
  int la=0, lb=0; uint8_t *a = read_fa_enc(f, &la), *b = read_fa_enc(f, &lb); fclose(f);
  int8_t mat[25]; for (int i=0;i<5;i++) for (int j=0;j<5;j++) mat[i*5+j] = (i==j && i<4) ? 1 : -2;
  ksw_extz_t ez; memset(&ez, 0, sizeof(ez));
  ksw_extz2_sse(NULL, la, a, lb, b, 5, mat, 4, 1, w, -1, 0, 0, &ez);
  printf("score: %d\nn_cigar: %d\n", ez.score, ez.n_cigar);
  free(ez.cigar); free(a); free(b);
  return 0;
}
EOF

gcc -O3 -msse4.1 -I "$WORK/src/ksw2" \
  "$WORK/src/run_ksw2_banded.c" "$WORK/src/ksw2/ksw2_extz2_sse.c" \
  -o "$WORK/bin/run_ksw2_banded"

# --- Run all benchmark sections via Python driver ---
python3 - "$ROOT" "$WORK" "$THREADS" <<'PY'
import os, random, re, statistics, subprocess, sys, time

ROOT, WORK, THREADS = sys.argv[1], sys.argv[2], int(sys.argv[3])
NW = os.path.join(ROOT, "bin/needleman_wunsch")
SW = os.path.join(ROOT, "bin/smith_waterman")
WFA2 = os.path.join(WORK, "bin/run_wfa2")
KSW2 = os.path.join(WORK, "bin/run_ksw2")
KSW2_B = os.path.join(WORK, "bin/run_ksw2_banded")
PARA_DNA = os.path.join(WORK, "bin/run_parasail_dna")
INPROC = os.path.join(WORK, "bin/inproc_prot")

def run_one(cmd):
    t0 = time.perf_counter()
    p = subprocess.run(["/usr/bin/time", "-v"] + cmd, capture_output=True, text=True, check=True)
    dt = (time.perf_counter() - t0) * 1000.0
    rss = int(re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", p.stderr).group(1)) / 1024.0
    usr = float(re.search(r"User time \(seconds\):\s*([\d.]+)", p.stderr).group(1))
    sys_t = float(re.search(r"System time \(seconds\):\s*([\d.]+)", p.stderr).group(1))
    m_sc = re.findall(r"(?:^|\s)score:\s*(-?\d+)", p.stdout)
    sc = int(m_sc[-1]) if m_sc else None
    return dt, (usr + sys_t) * 1000.0, rss, sc

def run_seeds(tag, cmd_fn, seeds, expected_scores=None):
    ws, cs, rs, scs = [], [], [], []
    for i, s in enumerate(seeds):
        w, c, r, sc = run_one(cmd_fn(s))
        ws.append(w); cs.append(c); rs.append(r); scs.append(sc)
        if expected_scores is not None:
            assert sc == expected_scores[i], f"{tag} seed {s}: {sc} != {expected_scores[i]}"
    print(f"  {tag:20s} wall={statistics.median(ws):8.2f} ms [{min(ws):7.2f}-{max(ws):7.2f}]  cpu={statistics.median(cs):7.1f} ms  rss={statistics.median(rs):7.1f} MiB  scores={scs}")
    return scs

print("\n=== A. IN-PROCESS PROTEIN BENCHMARK (BLOSUM62 -10/-1, ZERO STARTUP OVERHEAD) ===")
for name, f1, f2 in [
    ("Human/Mouse TP53", "h_tp53.fa", "m_tp53.fa"),
    ("Human POLR2A/Yeast Rpo21", "h_polr2a.fa", "y_rpo21.fa"),
    ("Human/Mouse BRCA1", "h_brca1.fa", "m_brca1.fa"),
    ("Human/Mouse MTOR", "h_mtor.fa", "m_mtor.fa"),
    ("Human/Mouse ATM", "h_atm.fa", "m_atm.fa"),
]:
    out = subprocess.check_output([INPROC, name, os.path.join(WORK, "proteins", f1), os.path.join(WORK, "proteins", f2), "15"], text=True)
    print(out, end="")

print("\n=== B. PARASAIL AVX2 32-BIT GLOBAL DNA (20 kb & 50 kb, 5 SEEDS, +1/-2/-4/-1) ===")
for length in [20000, 50000]:
    for div in ["1", "10", "rand"]:
        print(f"-- Global {length//1000} kb {div} --")
        # Get exact scores from Nucleics --scoreonly
        exp = []
        for s in range(1, 6):
            fa = os.path.join(WORK, "pairs", f"{length}_{div}_s{s}.fa")
            _, _, _, sc = run_one([NW, "--scoreonly", "--threads", "8", "--file", fa])
            exp.append(sc)
        run_seeds("parasail_scan32_1t", lambda s: [PARA_DNA, "nw", os.path.join(WORK, "pairs", f"{length}_{div}_s{s}.fa")], range(1, 6), exp)

print("\n=== C. BANDED KSW2 AT CERTIFIED BAND WIDTH (20 kb & 50 kb 1%, 5 SEEDS) ===")
for length in [20000, 50000]:
    print(f"-- Banded ksw2 {length//1000} kb 1% --")
    exp = []
    bands = []
    for s in range(1, 6):
        fa = os.path.join(WORK, "pairs", f"{length}_1_s{s}.fa")
        _, _, _, sc = run_one([NW, "--scoreonly", "--threads", "8", "--file", fa])
        exp.append(sc)
        # Certified band width: b = floor((s_max*min(m,n) + g_ext*|D| - L) / (s_max - 2*g_ext))
        # Read lengths
        seqs = []
        cur = []
        for line in open(fa):
            if line.startswith(">"):
                if cur: seqs.append("".join(cur)); cur = []
            else:
                cur.append(line.strip())
        if cur: seqs.append("".join(cur))
        m, n = len(seqs[0]), len(seqs[1])
        b = (1 * min(m, n) - 1 * abs(m - n) - sc) // 3
        bands.append(b)
    print(f"  certified_bands={bands}")
    run_seeds("ksw2_banded_1t", lambda s: [KSW2_B, os.path.join(WORK, "pairs", f"{length}_1_s{s}.fa"), str(bands[s-1])], range(1, 6), exp)

print("\n=== D. PARASAIL 32-BIT LOCAL & NUCLEICS --fastlocal (20 kb, 5 SEEDS, +2/-2/-2/-1) ===")
for div in ["1", "10", "rand"]:
    print(f"-- Local 20 kb {div} --")
    exp = []
    for s in range(1, 6):
        fa = os.path.join(WORK, "pairs", f"20000_{div}_s{s}.fa")
        _, _, _, sc = run_one([SW, "--maxhits", "1", "--minscore", "1", "--threads", "8", "--file", fa])
        exp.append(sc)
    run_seeds("parasail_sw32_1t", lambda s: [PARA_DNA, "sw", os.path.join(WORK, "pairs", f"20000_{div}_s{s}.fa")], range(1, 6), exp)
    run_seeds("sw_fastlocal_1t", lambda s: [SW, "--maxhits", "1", "--minscore", "1", "--fastlocal", "--threads", "1", "--file", os.path.join(WORK, "pairs", f"20000_{div}_s{s}.fa")], range(1, 6), exp)
    run_seeds("sw_fastlocal_8t", lambda s: [SW, "--maxhits", "1", "--minscore", "1", "--fastlocal", "--threads", "8", "--file", os.path.join(WORK, "pairs", f"20000_{div}_s{s}.fa")], range(1, 6), exp)

# Also check --fastlocal on real Mpox 197 kb if cached
mpox_a = "/tmp/seq_real_bench/mpox_ref.fa"
mpox_b = "/tmp/seq_real_bench/mpox_2022.fa"
if os.path.exists(mpox_a) and os.path.exists(mpox_b):
    print("-- Local Mpox 197 kb --")
    w1, c1, r1, sc1 = run_one([SW, "--maxhits", "1", "--minscore", "1", "--threads", "8", "--files", mpox_a, mpox_b])
    w2, c2, r2, sc2 = run_one([SW, "--maxhits", "1", "--minscore", "1", "--fastlocal", "--threads", "8", "--files", mpox_a, mpox_b])
    assert sc1 == sc2, f"Mpox fastlocal score mismatch: {sc1} vs {sc2}"
    print(f"  sw_default_8t        wall={w1:8.2f} ms  cpu={c1:7.1f} ms  rss={r1:7.1f} MiB  score={sc1}")
    print(f"  sw_fastlocal_8t      wall={w2:8.2f} ms  cpu={c2:7.1f} ms  rss={r2:7.1f} MiB  score={sc2}")

print("\n=== E. PREVIOUSLY SKIPPED 50 kb CELLS (5 SEEDS, +1/-2/-4/-1) ===")
for div in ["10", "rand"]:
    print(f"-- 50 kb {div} --")
    exp = []
    for s in range(1, 6):
        fa = os.path.join(WORK, "pairs", f"50000_{div}_s{s}.fa")
        _, _, _, sc = run_one([NW, "--scoreonly", "--threads", "8", "--file", fa])
        exp.append(sc)
    run_seeds("nw_wfa_1t", lambda s: [NW, "--wfa", "--threads", "1", "--printscores", "--file", os.path.join(WORK, "pairs", f"50000_{div}_s{s}.fa")], range(1, 6), exp)
    if div == "rand":
        run_seeds("wfa2_biwfa_1t", lambda s: [WFA2, "biwfa", os.path.join(WORK, "pairs", f"50000_rand_s{s}.fa")], range(1, 6), exp)

print("\n=== F. THREAD SCALING CURVE (1, 2, 4, 8, 16 THREADS, SEED 1, MEDIAN OF 3) ===")
for length, div in [(20000, "10"), (50000, "10"), (50000, "rand")]:
    fa = os.path.join(WORK, "pairs", f"{length}_{div}_s1.fa")
    print(f"-- Thread scaling {length//1000} kb {div} --")
    for mode, extra in [("full", []), ("linear", ["--linear"])]:
        base = None
        for thr in [1, 2, 4, 8, 16]:
            ws, cs, rs = [], [], []
            sc = None
            for _ in range(3):
                w, c, r, sc = run_one([NW, "--threads", str(thr), "--printscores", "--file", fa] + extra)
                ws.append(w); cs.append(c); rs.append(r)
            med_w = statistics.median(ws)
            if base is None: base = med_w
            print(f"  {mode:6s} {thr:2d}t  wall={med_w:8.2f} ms [{min(ws):7.2f}-{max(ws):7.2f}]  cpu={statistics.median(cs):7.1f} ms  rss={statistics.median(rs):7.1f} MiB  speedup={base/med_w:5.2f}x  score={sc}")

print("\n=== G. LONG-INDEL & UNEQUAL-LENGTH WORKLOAD (20 kb vs ~18.5 kb, 50-200 bp BLOCK INDELS, 3 SEEDS) ===")
bases = "ACGT"
for s in range(1, 4):
    rng = random.Random(1000 + s)
    a = [rng.choice(bases) for _ in range(20000)]
    b = list(a)
    # Apply 1% point substitutions
    for i in range(len(b)):
        if rng.random() < 0.01:
            b[i] = rng.choice([c for c in bases if c != b[i]])
    # Apply 12 deletions (50-200 bp) and 4 insertions (50-150 bp) -> unequal lengths ~18.5-19 kb
    events = []
    pos = 400
    for k in range(16):
        pos += rng.randint(600, 1100)
        if pos >= len(b) - 300: break
        if k % 4 == 0:
            ins_len = rng.randint(50, 150)
            ins = [rng.choice(bases) for _ in range(ins_len)]
            b[pos:pos] = ins
            pos += ins_len
        else:
            del_len = rng.randint(50, 200)
            del b[pos:pos + del_len]
    fa = os.path.join(WORK, "long_indel", f"long_indel_s{s}.fa")
    with open(fa, "w") as f:
        f.write(f">a_{s}\n{''.join(a)}\n>b_{s}\n{''.join(b)}\n")
    print(f"  seed {s}: len_a={len(a)}, len_b={len(b)}")

exp_li = []
for s in range(1, 4):
    fa = os.path.join(WORK, "long_indel", f"long_indel_s{s}.fa")
    _, _, _, sc = run_one([NW, "--scoreonly", "--threads", "8", "--file", fa])
    exp_li.append(sc)
run_seeds("nw_full_1t", lambda s: [NW, "--threads", "1", "--printscores", "--file", os.path.join(WORK, "long_indel", f"long_indel_s{s}.fa")], range(1, 4), exp_li)
run_seeds("nw_full_8t", lambda s: [NW, "--threads", "8", "--printscores", "--file", os.path.join(WORK, "long_indel", f"long_indel_s{s}.fa")], range(1, 4), exp_li)
run_seeds("nw_linear_1t", lambda s: [NW, "--linear", "--threads", "1", "--printscores", "--file", os.path.join(WORK, "long_indel", f"long_indel_s{s}.fa")], range(1, 4), exp_li)
run_seeds("nw_linear_8t", lambda s: [NW, "--linear", "--threads", "8", "--printscores", "--file", os.path.join(WORK, "long_indel", f"long_indel_s{s}.fa")], range(1, 4), exp_li)
run_seeds("nw_wfa_1t", lambda s: [NW, "--wfa", "--threads", "1", "--printscores", "--file", os.path.join(WORK, "long_indel", f"long_indel_s{s}.fa")], range(1, 4), exp_li)
run_seeds("ksw2_1t", lambda s: [KSW2, os.path.join(WORK, "long_indel", f"long_indel_s{s}.fa")], range(1, 4), exp_li)
run_seeds("wfa2_biwfa_1t", lambda s: [WFA2, "biwfa", os.path.join(WORK, "long_indel", f"long_indel_s{s}.fa")], range(1, 4), exp_li)
run_seeds("parasail_scan32_1t", lambda s: [PARA_DNA, "nw", os.path.join(WORK, "long_indel", f"long_indel_s{s}.fa")], range(1, 4), exp_li)
PY
