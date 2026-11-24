#!/usr/bin/env bash
# Journal-grade supplementary benchmark for Nucleics Align v1.3.0:
#   1. Builds pinned parasail v2.6.2 (600fb26, AVX2 16-bit striped and scan global traceback)
#   2. Benchmarks 5 real NCBI RefSeq protein pairs (BLOSUM62 -10/-1) with parasail
#      + independent BLOSUM62 CIGAR/alignment re-scoring across tools (Table 5)
#   3. Independently re-scores cross-tool 20 kb DNA alignments/CIGARs
#      (Nucleics Align full/--linear/--wfa, ksw2, WFA2 BiWFA)
#   4. Re-runs full 5 kb - 1 Mb scaling suite on v1.3.0 with score assertions
#      and matched WFA2 BiWFA comparisons at 100 kb, 300 kb, and 1 Mb (Table 3)
#
# Usage: bash benchmarks/journal_bench.sh [work_dir] [threads]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:-/tmp/nucleics_compare}"
THREADS="${2:-8}"
mkdir -p "$WORK/src" "$WORK/bin" "$WORK/scale"

NW="$ROOT/bin/needleman_wunsch"
SW="$ROOT/bin/smith_waterman"
ORIG_NW="$WORK/src/seq-align/bin/needleman_wunsch"
WFA2="$WORK/bin/run_wfa2"
KSW2="$WORK/bin/run_ksw2"

for req in "$WORK/bin/run_wfa2" "$WORK/src/ksw2/ksw2.h" "$WORK/src/WFA2-lib/lib/libwfa.a" "$WORK/proteins/h_tp53.fa" "$WORK/pairs/20000_1_s1.fa"; do
  if [ ! -e "$req" ]; then
    echo "Missing prerequisite $req — run benchmarks/compare_harness.sh and benchmarks/protein_bench.sh first with the same work_dir ($WORK)" >&2
    exit 1
  fi
done

echo "# host: $(uname -a)"
echo "# cpu:  $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')"
echo "# governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)"
echo "# smt:  $(cat /sys/devices/system/cpu/smt/active 2>/dev/null || echo unknown)"
echo "# gcc:  $(gcc --version | head -n1)"
echo "# python: $(python3 --version)"
echo "# nucleics-align: $(git -C "$ROOT" rev-parse --short HEAD)"
echo "# nw_sha256: $(sha256sum "$NW" | awk '{print $1}')"
echo "# sw_sha256: $(sha256sum "$SW" | awk '{print $1}')"
echo "# ksw2: $(git -C "$WORK/src/ksw2" rev-parse --short HEAD)"
echo "# WFA2-lib: $(git -C "$WORK/src/WFA2-lib" rev-parse --short HEAD)"

# 1. Build pinned parasail v2.6.2 (600fb26, AVX2 256-bit 16-bit striped NW + traceback)
if [ ! -d "$WORK/src/parasail" ]; then
  git clone --branch v2.6.2 https://github.com/jeffdaily/parasail "$WORK/src/parasail"
fi
git -C "$WORK/src/parasail" checkout -q v2.6.2
echo "# parasail: $(git -C "$WORK/src/parasail" rev-parse --short HEAD) (v2.6.2)"

cat > "$WORK/src/config.h" <<'EOF'
#define HAVE_POSIX_MEMALIGN 1
#define HAVE_AVX2 1
#define HAVE_SSE41 1
#define HAVE_SSE2 1
EOF

cat > "$WORK/src/run_parasail_prot.c" <<'EOF'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parasail.h"
#include "parasail/matrices/blosum62.h"

static char *read_fa(const char *path, int *out_len) {
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
  *out_len = (int)len;
  return buf;
}

int main(int argc, char **argv) {
  if (argc < 3) return 1;
  int print_cigar = (argc >= 4 && strcmp(argv[3], "--cigar") == 0);
  int la = 0, lb = 0;
  char *a = read_fa(argv[1], &la);
  char *b = read_fa(argv[2], &lb);
  /* Parasail open=11, extend=1 equals Gotoh gap_open=-10, gap_extend=-1 */
  parasail_result_t *res = parasail_nw_trace_striped_avx2_256_16(
      a, la, b, lb, 11, 1, &parasail_blosum62);
  if (!res || parasail_result_is_saturated(res)) {
    fprintf(stderr, "parasail saturated or failed\n");
    return 2;
  }
  parasail_cigar_t *cigar = parasail_result_get_cigar(res, a, la, b, lb, &parasail_blosum62);
  char *cstr = parasail_cigar_decode(cigar);
  printf("score: %d\n", res->score);
  if (print_cigar) printf("cigar: %s\n", cstr);
  free(cstr);
  parasail_cigar_free(cigar);
  parasail_result_free(res);
  free(a); free(b);
  return 0;
}
EOF

gcc -O3 -mavx2 -I "$WORK/src" \
  -I "$WORK/src/parasail" -I "$WORK/src/parasail/src" \
  "$WORK/src/run_parasail_prot.c" \
  "$WORK/src/parasail/src/nw_trace_striped_avx2_256_16.c" \
  "$WORK/src/parasail/src/memory.c" \
  "$WORK/src/parasail/src/memory_avx2.c" \
  "$WORK/src/parasail/src/cigar.c" \
  "$WORK/src/parasail/src/matrix_lookup.c" \
  -o "$WORK/bin/run_parasail_prot"

sed 's/parasail_nw_trace_striped_avx2_256_16/parasail_nw_trace_scan_avx2_256_16/g' \
  "$WORK/src/run_parasail_prot.c" > "$WORK/src/run_parasail_scan.c"

gcc -O3 -mavx2 -I "$WORK/src" \
  -I "$WORK/src/parasail" -I "$WORK/src/parasail/src" \
  "$WORK/src/run_parasail_scan.c" \
  "$WORK/src/parasail/src/nw_trace_scan_avx2_256_16.c" \
  "$WORK/src/parasail/src/memory.c" \
  "$WORK/src/parasail/src/memory_avx2.c" \
  "$WORK/src/parasail/src/cigar.c" \
  "$WORK/src/parasail/src/matrix_lookup.c" \
  -o "$WORK/bin/run_parasail_scan"

# Also build small CIGAR-printing wrappers for ksw2 and WFA2 BiWFA for independent re-scoring
cat > "$WORK/src/run_ksw2_cigar.c" <<'EOF'
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
  FILE *f = fopen(argv[1], "r"); int la=0, lb=0; uint8_t *a = read_fa_enc(f, &la), *b = read_fa_enc(f, &lb); fclose(f);
  int8_t mat[25]; for (int i=0;i<5;i++) for (int j=0;j<5;j++) mat[i*5+j] = (i==j && i<4) ? 1 : -2;
  ksw_extz_t ez; memset(&ez, 0, sizeof(ez));
  ksw_extz2_sse(NULL, la, a, lb, b, 5, mat, 4, 1, -1, -1, 0, 0, &ez);
  printf("score: %d\ncigar: ", ez.score);
  for (int i = 0; i < ez.n_cigar; i++) printf("%u%c", ez.cigar[i]>>4, "MIDNSHP=X"[ez.cigar[i]&0xf]);
  printf("\n"); free(ez.cigar); free(a); free(b); return 0;
}
EOF
gcc -O3 -msse4.1 "$WORK/src/run_ksw2_cigar.c" "$WORK/src/ksw2/ksw2_extz2_sse.c" \
  -I "$WORK/src/ksw2" -o "$WORK/bin/run_ksw2_cigar"

cat > "$WORK/src/run_wfa2_cigar.c" <<'EOF'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wavefront/wavefront_align.h"
static char *read_fa(FILE *f, int *out_len) {
  size_t cap = 4096, len = 0; char *buf = malloc(cap); char line[8192]; int started = 0;
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '>') { if (started) { fseek(f, -(long)strlen(line), SEEK_CUR); break; } started = 1; continue; }
    started = 1;
    for (char *p = line; *p; p++) { if (*p <= ' ') continue; if (len+2>cap){cap*=2;buf=realloc(buf,cap);} buf[len++]=*p; }
  }
  buf[len] = '\0'; *out_len = (int)len; return buf;
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); int la=0, lb=0; char *a = read_fa(f, &la), *b = read_fa(f, &lb); fclose(f);
  wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
  attr.distance_metric = gap_affine;
  attr.affine_penalties.match = 0; attr.affine_penalties.mismatch = 6;
  attr.affine_penalties.gap_opening = 8; attr.affine_penalties.gap_extension = 3;
  attr.memory_mode = wavefront_memory_ultralow;
  attr.heuristic.strategy = wf_heuristic_none;
  attr.alignment_scope = compute_alignment;
  wavefront_aligner_t *wf = wavefront_aligner_new(&attr);
  if (wavefront_align(wf, a, la, b, lb) != 0) return 2;
  int score = (la + lb + wf->cigar->score) / 2;
  printf("score: %d\ncigar: ", score);
  fwrite(wf->cigar->operations + wf->cigar->begin_offset, 1,
         wf->cigar->end_offset - wf->cigar->begin_offset, stdout);
  printf("\n"); wavefront_aligner_delete(wf); free(a); free(b); return 0;
}
EOF
gcc -O3 "$WORK/src/run_wfa2_cigar.c" -I "$WORK/src/WFA2-lib" "$WORK/src/WFA2-lib/lib/libwfa.a" \
  -lm -o "$WORK/bin/run_wfa2_cigar"

python3 - "$ROOT" "$WORK" "$THREADS" <<'PY'
import os, re, statistics, subprocess, sys, time

ROOT, WORK, THREADS = sys.argv[1], sys.argv[2], sys.argv[3]
NW = os.path.join(ROOT, "bin/needleman_wunsch")
SW = os.path.join(ROOT, "bin/smith_waterman")
PARASAIL = os.path.join(WORK, "bin/run_parasail_prot")
WFA2 = os.path.join(WORK, "bin/run_wfa2")
WFA2_CIG = os.path.join(WORK, "bin/run_wfa2_cigar")
KSW2_CIG = os.path.join(WORK, "bin/run_ksw2_cigar")

def load_blosum62(path):
    lines = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
    cols = lines[0].split()
    mat = {}
    for line in lines[1:]:
        parts = line.split()
        r = parts[0]
        for c, v in zip(cols, parts[1:]):
            mat[(r, c)] = int(v)
    return mat

BLOSUM62 = load_blosum62(os.path.join(ROOT, "scoring/BLOSUM62.txt"))

def read_seqs(path):
    seqs, cur = [], []
    for l in open(path):
        if l.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        else:
            cur.append(l.strip())
    if cur: seqs.append("".join(cur))
    return seqs

def rescore_aligned(a, b, sub_fn, go, ge):
    assert len(a) == len(b)
    sc, in_a, in_b = 0, False, False
    for x, y in zip(a, b):
        assert not (x == "-" and y == "-"), "Invalid double-gap column"
        if x == "-":
            sc += ge if in_a else (go + ge); in_a, in_b = True, False
        elif y == "-":
            sc += ge if in_b else (go + ge); in_b, in_a = True, False
        else:
            in_a = in_b = False; sc += sub_fn(x.upper(), y.upper())
    return sc

def rescore_cigar(seq_a, seq_b, cigar_str, sub_fn, go, ge, swap_id=False):
    i = j = sc = 0
    if re.match(r"^[MIDX=]+$", cigar_str):
        ops = [(1, ch) for ch in cigar_str]
    else:
        assert re.fullmatch(r"(?:[1-9]\d*[MIDX=])+", cigar_str), f"Malformed CIGAR: {cigar_str[:40]}"
        ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDX=])", cigar_str)]
    prev_op = None
    for n, op in ops:
        if op in ("M", "=", "X"):
            for _ in range(n):
                sc += sub_fn(seq_a[i].upper(), seq_b[j].upper())
                i += 1; j += 1
            prev_op = "M"
        elif (op == "I" and not swap_id) or (op == "D" and swap_id):
            sc += (ge * n) if prev_op == "I" else (go + ge * n)
            j += n; prev_op = "I"
        else:
            sc += (ge * n) if prev_op == "D" else (go + ge * n)
            i += n; prev_op = "D"
    assert i == len(seq_a) and j == len(seq_b), f"Span mismatch: ({i},{j}) vs ({len(seq_a)},{len(seq_b)})"
    return sc

def run_timed(cmd, reps=1):
    walls, cpus, rsss = [], [], []
    last_out = ""
    for _ in range(reps):
        t0 = time.perf_counter()
        p = subprocess.run(["/usr/bin/time", "-v"] + cmd, capture_output=True, text=True, check=True)
        dt = (time.perf_counter() - t0) * 1000.0
        rss = int(re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", p.stderr).group(1)) / 1024.0
        usr = float(re.search(r"User time \(seconds\):\s*([\d.]+)", p.stderr).group(1))
        sys_t = float(re.search(r"System time \(seconds\):\s*([\d.]+)", p.stderr).group(1))
        walls.append(dt); cpus.append((usr + sys_t) * 1000.0); rsss.append(rss)
        last_out = p.stdout
    m_sc = re.findall(r"(?:^|\s)score:\s*(-?\d+)", last_out)
    sc = int(m_sc[-1]) if m_sc else None
    return statistics.median(walls), min(walls), max(walls), statistics.median(cpus), statistics.median(rsss), sc, last_out

print("\n=== PARASAIL v2.6.2 (AVX2 16-bit striped, 1t) ON 5 REAL PROTEINS + INDEPENDENT RE-SCORING ===")
PROT_DIR = os.path.join(WORK, "proteins")
pairs = [
    ("Human/Mouse TP53", "h_tp53.fa", "m_tp53.fa"),
    ("Human POLR2A/Yeast Rpo21", "h_polr2a.fa", "y_rpo21.fa"),
    ("Human/Mouse BRCA1", "h_brca1.fa", "m_brca1.fa"),
    ("Human/Mouse MTOR", "h_mtor.fa", "m_mtor.fa"),
    ("Human/Mouse ATM", "h_atm.fa", "m_atm.fa"),
]
b62_fn = lambda x, y: BLOSUM62[(x, y)]
for name, f1, f2 in pairs:
    p1, p2 = os.path.join(PROT_DIR, f1), os.path.join(PROT_DIR, f2)
    sa, sb = read_seqs(p1)[0], read_seqs(p2)[0]
    med, lo, hi, cpu, rss, psc, _ = run_timed([PARASAIL, p1, p2], reps=3)
    # Re-score parasail CIGAR independently
    p_out = subprocess.check_output([PARASAIL, p1, p2, "--cigar"], text=True)
    cig = re.search(r"cigar:\s*(\S+)", p_out).group(1)
    # Parasail cigar: I consumes query (a), D consumes ref (b) -> swap_id=True
    rsc_p = rescore_cigar(sa, sb, cig, b62_fn, -10, -1, swap_id=True)
    assert rsc_p == psc, f"Parasail rescore mismatch {name}: {rsc_p} vs {psc}"
    # Re-score Nucleics Align full and --linear independently
    for mode, extra in [("full", []), ("linear", ["--linear"])]:
        out = subprocess.check_output([NW, "--threads", "1", "--scoring", "BLOSUM62", "--printscores", "--files", p1, p2] + extra, text=True)
        lines = [l for l in out.splitlines() if l.strip()]
        aa, ab = lines[0], lines[1]
        assert aa.replace("-", "") == sa and ab.replace("-", "") == sb
        rsc_n = rescore_aligned(aa, ab, b62_fn, -10, -1)
        assert rsc_n == psc, f"Nucleics {mode} rescore mismatch {name}: {rsc_n} vs {psc}"
    print(f"  {name:26s}  parasail_1t={med:7.2f} ms [{lo:6.2f}-{hi:6.2f}]  rss={rss:6.1f} MiB  score={psc} (parasail/full/linear CIGAR+span verified)")

print("\n=== INDEPENDENT CROSS-TOOL DNA CIGAR/ALIGNMENT RE-SCORING (20 kb seed 1, +1/-2/-4/-1) ===")
dna_fn = lambda x, y: 1 if x == y else -2
for div in ["1", "10", "rand"]:
    fa = os.path.join(WORK, "pairs", f"20000_{div}_s1.fa")
    sa, sb = read_seqs(fa)
    scores = {}
    for mode, extra in [("nw_full", []), ("nw_linear", ["--linear"]), ("nw_wfa", ["--wfa"])]:
        out = subprocess.check_output([NW, "--threads", "1", "--printscores", "--file", fa] + extra, text=True)
        lines = [l for l in out.splitlines() if l.strip()]
        aa, ab = lines[0], lines[1]
        assert aa.replace("-", "") == sa and ab.replace("-", "") == sb
        rep = int(re.search(r"score:\s*(-?\d+)", out).group(1))
        rsc = rescore_aligned(aa, ab, dna_fn, -4, -1)
        assert rsc == rep, f"{mode} rescore mismatch: {rsc} vs {rep}"
        scores[mode] = rsc
    # ksw2 CIGAR (ksw2: query=a, target=b; I consumes query a, D consumes target b -> swap_id=True)
    k_out = subprocess.check_output([KSW2_CIG, fa], text=True)
    k_rep = int(re.search(r"score:\s*(-?\d+)", k_out).group(1))
    k_cig = re.search(r"cigar:\s*(\S+)", k_out).group(1)
    k_rsc = rescore_cigar(sa, sb, k_cig, dna_fn, -4, -1, swap_id=True)
    assert k_rsc == k_rep
    scores["ksw2"] = k_rsc
    # WFA2 BiWFA CIGAR (pattern=a, text=b; I consumes text b, D consumes pattern a -> swap_id=False)
    w_out = subprocess.check_output([WFA2_CIG, fa], text=True)
    w_rep = int(re.search(r"score:\s*(-?\d+)", w_out).group(1))
    w_cig = re.search(r"cigar:\s*(\S+)", w_out).group(1)
    w_rsc = rescore_cigar(sa, sb, w_cig, dna_fn, -4, -1, swap_id=False)
    assert w_rsc == w_rep
    scores["biwfa"] = w_rsc
    assert len(set(scores.values())) == 1, f"Cross-tool mismatch {div}: {scores}"
    print(f"  20 kb {div:4s}  ALL 5 TOOLS RE-SCORED & SPAN-VERIFIED: score={k_rsc} ({list(scores.keys())})")

print("\n=== v1.3.0 SCALING TABLE 3 (5 kb - 1 Mb, 8 default threads, gen_pairs.py seed=11, score-asserted) ===")
SCALE_DIR = os.path.join(WORK, "scale")
GEN = os.path.join(ROOT, "benchmarks/gen_pairs.py")

cases = [
    ("5 kb",   "5000",    "1",    "5000_1_s11.fa",     3, True,  False),
    ("5 kb",   "5000",    "10",   "5000_10_s11.fa",    3, True,  False),
    ("5 kb",   "5000",    "rand", "5000_rand_s11.fa",  3, True,  False),
    ("20 kb",  "20000",   "1",    "20000_1_s11.fa",    3, True,  False),
    ("20 kb",  "20000",   "10",   "20000_10_s11.fa",   3, True,  False),
    ("20 kb",  "20000",   "rand", "20000_rand_s11.fa", 3, True,  False),
    ("50 kb",  "50000",   "1",    "50000_1_s11.fa",    3, True,  False),
    ("50 kb",  "50000",   "10",   "50000_10_s11.fa",   3, True,  False),
    ("50 kb",  "50000",   "rand", "50000_rand_s11.fa", 3, True,  False),
    ("100 kb", "100000",  "0.1",  "100000_0.1_s11.fa", 1, False, True),
    ("100 kb", "100000",  "1",    "100000_1_s11.fa",   1, False, True),
    ("100 kb", "100000",  "5",    "100000_5_s11.fa",   1, False, False),
    ("300 kb", "300000",  "0.1",  "300000_0.1_s11.fa", 1, False, True),
    ("300 kb", "300000",  "1",    "300000_1_s11.fa",   1, False, True),
    ("300 kb", "300000",  "5",    "300000_5_s11.fa",   1, False, False),
    ("1 Mb",   "1000000", "0.1",  "1000000_0.1_s11.fa",1, False, True),
    ("1 Mb",   "1000000", "1",    "1000000_1_s11.fa",  1, False, True),
    ("1 Mb",   "1000000", "5",    "1000000_5_s11.fa",  1, False, False),
]

import hashlib
for label, length, mut, fname, reps, do_loc, do_biwfa in cases:
    fa = os.path.join(SCALE_DIR, fname)
    subprocess.check_call(["python3", GEN, length, mut, "11", fa])
    h = hashlib.sha256(open(fa, "rb").read()).hexdigest()
    print(f"# sha256 {fname} = {h}")
    g_scores = {}
    row = [f"{label:6s} {mut:4s}"]
    for mname, cmd in [
        ("full",      [NW, "--threads", THREADS, "--printscores", "--file", fa]),
        ("wfa",       [NW, "--threads", THREADS, "--wfa", "--printscores", "--file", fa]),
        ("linear",    [NW, "--threads", THREADS, "--linear", "--printscores", "--file", fa]),
        ("scoreonly", [NW, "--threads", THREADS, "--scoreonly", "--file", fa]),
    ]:
        if length == "1000000" and mut == "5" and mname in ("full", "wfa"):
            row.append(f"{mname}=---"); continue
        med, lo, hi, _, rss, sc, _ = run_timed(cmd, reps=reps)
        assert sc is not None, f"Missing score for {mname} on {label} {mut}"
        g_scores[mname] = sc
        row.append(f"{mname}={lo:7.1f}/{rss:5.0f}")
    if do_loc:
        _, lo, _, _, rss, _, _ = run_timed([SW, "--threads", THREADS, "--maxhits", "1", "--minscore", "1", "--file", fa], reps=reps)
        row.append(f"loc1={lo:6.1f}/{rss:3.0f}")
    if do_biwfa:
        _, lo, _, _, rss, sc, _ = run_timed([WFA2, "biwfa", fa], reps=reps)
        g_scores["wfa2_biwfa"] = sc
        row.append(f"wfa2_biwfa={lo:7.1f}/{rss:4.0f}")
    assert len(set(g_scores.values())) == 1, f"Scale score mismatch {label} {mut}: {g_scores}"
    sc_val = next(iter(g_scores.values()))
    print("  " + " | ".join(row) + f" | score={sc_val}")
PY