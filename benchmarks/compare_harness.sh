#!/usr/bin/env bash
# Reproducible multi-seed comparator & fresh-process ablation harness for Nucleics Align.
# Builds exact pinned commits of seq-align (dc41988), ksw2 (289609b),
# WFA2-lib (v2.3.5 / cf3eb92), SSW (a66636b), and edlib (v1.2.7 / 4215d82),
# then benchmarks:
#   1. 9-configuration fresh-process ablation (20 kb, 5 seeds, wall ms + peak RSS MiB)
#   2. Global affine-gap comparators (20 kb & 50 kb, 5 seeds)
#   3. Unit-cost edit distance (--scoreonly vs edlib, 20 kb & 50 kb, 5 seeds)
#   4. Local alignment (--maxhits 1 vs SSW, 20 kb, 5 seeds)
#
# Usage: bash benchmarks/compare_harness.sh [work_dir] [threads] [all|ablation|global|edlib|local]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:-/tmp/nucleics_compare}"
THREADS="${2:-8}"
SECTION="${3:-all}"
mkdir -p "$WORK/src" "$WORK/bin" "$WORK/pairs"

NW="$ROOT/bin/needleman_wunsch"
SW="$ROOT/bin/smith_waterman"
. "$ROOT/benchmarks/bench_common.sh"
bench_require

SEQALIGN_REV="dc41988"
KSW2_REV="289609b"
WFA2_REV="cf3eb92"
SSW_REV="a66636b"
EDLIB_TAG="v1.2.7"

bench_header

# 1. Original seq-align (pinned commit dc41988)
if [ ! -d "$WORK/src/seq-align" ]; then
  git clone --recursive https://github.com/noporpoise/seq-align "$WORK/src/seq-align"
fi
git -C "$WORK/src/seq-align" checkout -q "$SEQALIGN_REV"
git -C "$WORK/src/seq-align" submodule update --init --recursive -q
CUR_SA="$(git -C "$WORK/src/seq-align" rev-parse HEAD)"
if [ ! -x "$WORK/src/seq-align/bin/needleman_wunsch" ] || [ "$(cat "$WORK/.stamp_seqalign" 2>/dev/null || true)" != "$CUR_SA" ]; then
  make -C "$WORK/src/seq-align" clean >/dev/null 2>&1 || true
  make -C "$WORK/src/seq-align" -j4 >/dev/null 2>&1
  echo "$CUR_SA" > "$WORK/.stamp_seqalign"
fi
echo "# seq-align: $(git -C "$WORK/src/seq-align" rev-parse --short HEAD)"

# 2. ksw2 (pinned commit 289609b)
if [ ! -d "$WORK/src/ksw2" ]; then
  git clone https://github.com/lh3/ksw2 "$WORK/src/ksw2"
fi
git -C "$WORK/src/ksw2" checkout -q "$KSW2_REV"
echo "# ksw2: $(git -C "$WORK/src/ksw2" rev-parse --short HEAD)"

cat > "$WORK/src/run_ksw2.c" <<'EOF'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "ksw2.h"

static uint8_t *read_fa_enc(FILE *f, int *out_len) {
  size_t cap = 4096, len = 0;
  uint8_t *buf = malloc(cap);
  char line[8192];
  int started = 0;
  static const uint8_t tbl[256] = {
    ['A']=0,['a']=0,['C']=1,['c']=1,['G']=2,['g']=2,['T']=3,['t']=3
  };
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '>') {
      if (started) { fseek(f, -(long)strlen(line), SEEK_CUR); break; }
      started = 1; continue;
    }
    started = 1;
    for (char *p = line; *p; p++) {
      if (*p <= ' ') continue;
      if (len + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
      buf[len++] = tbl[(unsigned char)*p];
    }
  }
  *out_len = (int)len;
  return buf;
}

int main(int argc, char **argv) {
  if (argc < 2) return 1;
  FILE *f = fopen(argv[1], "r");
  if (!f) { perror(argv[1]); return 1; }
  int la = 0, lb = 0;
  uint8_t *a = read_fa_enc(f, &la);
  uint8_t *b = read_fa_enc(f, &lb);
  fclose(f);
  int8_t mat[25];
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      mat[i*5+j] = (i == j && i < 4) ? 1 : -2;
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ez));
  ksw_extz2_sse(NULL, la, a, lb, b, 5, mat, 4, 1, -1, -1, 0, 0, &ez);
  printf("score: %d\ncigar_len: %d\n", ez.score, ez.n_cigar);
  free(ez.cigar);
  free(a); free(b);
  return 0;
}
EOF
gcc -O3 -msse4.1 "$WORK/src/run_ksw2.c" "$WORK/src/ksw2/ksw2_extz2_sse.c" \
  -I "$WORK/src/ksw2" -o "$WORK/bin/run_ksw2"

# 3. WFA2-lib (pinned commit cf3eb92 / v2.3.5)
if [ ! -d "$WORK/src/WFA2-lib" ]; then
  git clone --branch v2.3.5 https://github.com/smarco/WFA2-lib "$WORK/src/WFA2-lib"
fi
git -C "$WORK/src/WFA2-lib" checkout -q "$WFA2_REV"
CUR_WFA="$(git -C "$WORK/src/WFA2-lib" rev-parse HEAD)"
if [ ! -f "$WORK/src/WFA2-lib/lib/libwfa.a" ] || [ "$(cat "$WORK/.stamp_wfa2" 2>/dev/null || true)" != "${CUR_WFA}-O3" ]; then
  make -C "$WORK/src/WFA2-lib" clean >/dev/null 2>&1 || true
  make -C "$WORK/src/WFA2-lib" BUILD_EXAMPLES=0 BUILD_TOOLS=0 all -j4 >/dev/null 2>&1
  echo "${CUR_WFA}-O3" > "$WORK/.stamp_wfa2"
fi
echo "# WFA2-lib: $(git -C "$WORK/src/WFA2-lib" rev-parse --short HEAD) (v2.3.5)"

cat > "$WORK/src/run_wfa2.c" <<'EOF'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wavefront/wavefront_align.h"

static char *read_fa(FILE *f, int *out_len) {
  size_t cap = 4096, len = 0;
  char *buf = malloc(cap);
  char line[8192];
  int started = 0;
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '>') {
      if (started) { fseek(f, -(long)strlen(line), SEEK_CUR); break; }
      started = 1; continue;
    }
    started = 1;
    for (char *p = line; *p; p++) {
      if (*p <= ' ') continue;
      if (len + 2 > cap) { cap *= 2; buf = realloc(buf, cap); }
      buf[len++] = *p;
    }
  }
  buf[len] = '\0';
  *out_len = (int)len;
  return buf;
}

int main(int argc, char **argv) {
  if (argc < 3) return 1;
  int biwfa = (strcmp(argv[1], "biwfa") == 0);
  FILE *f = fopen(argv[2], "r");
  if (!f) { perror(argv[2]); return 1; }
  int la = 0, lb = 0;
  char *a = read_fa(f, &la);
  char *b = read_fa(f, &lb);
  fclose(f);

  wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
  attr.distance_metric = gap_affine;
  attr.affine_penalties.match = 0;
  attr.affine_penalties.mismatch = 6;
  attr.affine_penalties.gap_opening = 8;
  attr.affine_penalties.gap_extension = 3;
  attr.memory_mode = biwfa ? wavefront_memory_ultralow : wavefront_memory_high;
  attr.heuristic.strategy = wf_heuristic_none;
  attr.alignment_scope = compute_alignment;

  wavefront_aligner_t *wf = wavefront_aligner_new(&attr);
  int st = wavefront_align(wf, a, la, b, lb);
  if (st != 0) {
    fprintf(stderr, "WFA2 failed with status %d\n", st);
    wavefront_aligner_delete(wf);
    free(a); free(b);
    return 2;
  }
  int penalty = -wf->cigar->score;
  int score = (la + lb - penalty) / 2;
  printf("score: %d\ncigar_len: %d\n", score, wf->cigar->end_offset - wf->cigar->begin_offset);
  wavefront_aligner_delete(wf);
  free(a); free(b);
  return 0;
}
EOF
gcc -O3 "$WORK/src/run_wfa2.c" -I "$WORK/src/WFA2-lib" "$WORK/src/WFA2-lib/lib/libwfa.a" \
  -lm -o "$WORK/bin/run_wfa2"

# 4. SSW (pinned commit a66636b)
if [ ! -d "$WORK/src/SSW" ]; then
  git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library "$WORK/src/SSW"
fi
git -C "$WORK/src/SSW" checkout -q "$SSW_REV"
echo "# SSW: $(git -C "$WORK/src/SSW" rev-parse --short HEAD)"

cat > "$WORK/src/run_ssw.c" <<'EOF'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "ssw.h"

static int8_t *read_fa_enc(FILE *f, int *out_len) {
  size_t cap = 4096, len = 0;
  int8_t *buf = malloc(cap);
  char line[8192];
  int started = 0;
  static const int8_t tbl[256] = {
    ['A']=0,['a']=0,['C']=1,['c']=1,['G']=2,['g']=2,['T']=3,['t']=3
  };
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '>') {
      if (started) { fseek(f, -(long)strlen(line), SEEK_CUR); break; }
      started = 1; continue;
    }
    started = 1;
    for (char *p = line; *p; p++) {
      if (*p <= ' ') continue;
      if (len + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
      buf[len++] = tbl[(unsigned char)*p];
    }
  }
  *out_len = (int)len;
  return buf;
}

int main(int argc, char **argv) {
  if (argc < 2) return 1;
  FILE *f = fopen(argv[1], "r");
  if (!f) { perror(argv[1]); return 1; }
  int la = 0, lb = 0;
  int8_t *a = read_fa_enc(f, &la);
  int8_t *b = read_fa_enc(f, &lb);
  fclose(f);
  int8_t mat[25];
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      mat[i*5+j] = (i == j && i < 4) ? 2 : -2;
  s_profile *prof = ssw_init(a, la, mat, 5, 2);
  s_align *res = ssw_align(prof, b, lb, 3, 1, 1, 0, 0, 15);
  printf("score: %d\nssw_flag: %u\n", res->score1, (unsigned)res->flag);
  align_destroy(res);
  init_destroy(prof);
  free(a); free(b);
  return 0;
}
EOF
gcc -O3 "$WORK/src/run_ssw.c" "$WORK/src/SSW/src/ssw.c" -I "$WORK/src/SSW/src" -lm -o "$WORK/bin/run_ssw"

# 5. edlib (pinned tag v1.2.7, unit-cost global edit distance score-only)
if [ ! -d "$WORK/src/edlib" ]; then
  git clone --branch "$EDLIB_TAG" https://github.com/Martinsos/edlib "$WORK/src/edlib"
fi
git -C "$WORK/src/edlib" checkout -q "$EDLIB_TAG"
echo "# edlib: $(git -C "$WORK/src/edlib" rev-parse --short HEAD) ($EDLIB_TAG)"

cat > "$WORK/src/run_edlib.cpp" <<'EOF'
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "edlib.h"

static std::string read_fa(FILE *f) {
  std::string s;
  char line[8192];
  int started = 0;
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '>') {
      if (started) { fseek(f, -(long)strlen(line), SEEK_CUR); break; }
      started = 1; continue;
    }
    started = 1;
    for (char *p = line; *p; p++) {
      if (*p > ' ') s.push_back(*p);
    }
  }
  return s;
}

int main(int argc, char **argv) {
  if (argc < 2) return 1;
  FILE *f = fopen(argv[1], "r");
  if (!f) { perror(argv[1]); return 1; }
  std::string a = read_fa(f);
  std::string b = read_fa(f);
  fclose(f);
  EdlibAlignResult res = edlibAlign(a.c_str(), (int)a.size(), b.c_str(), (int)b.size(),
                                    edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
  if (res.status != EDLIB_STATUS_OK) return 2;
  /* Unit-cost NW score (0/-1/0/-1) is minus edit distance */
  printf("score: %d\n", -res.editDistance);
  edlibFreeAlignResult(res);
  return 0;
}
EOF
g++ -O3 -std=c++11 "$WORK/src/run_edlib.cpp" "$WORK/src/edlib/edlib/src/edlib.cpp" \
  -I "$WORK/src/edlib/edlib/include" -o "$WORK/bin/run_edlib"

# 6. Build fresh-process ablation driver
gcc -O3 "$ROOT/benchmarks/ablation.c" -I "$ROOT/include" -I "$ROOT/src/lib" \
  -L "$ROOT/src" -lalign -lpthread \
  -o "$WORK/bin/ablation"

python3 - "$ROOT" "$WORK" "$THREADS" "$SECTION" <<'PY'
import hashlib, os, re, statistics, subprocess, sys, time
from bench_common import sha256, timed

ROOT, WORK, THREADS, SECTION = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
NW = os.path.join(ROOT, "bin/needleman_wunsch")
SW = os.path.join(ROOT, "bin/smith_waterman")
ORIG_NW = os.path.join(WORK, "src/seq-align/bin/needleman_wunsch")
KSW2 = os.path.join(WORK, "bin/run_ksw2")
WFA2 = os.path.join(WORK, "bin/run_wfa2")
SSW = os.path.join(WORK, "bin/run_ssw")
EDLIB = os.path.join(WORK, "bin/run_edlib")
ABL = os.path.join(WORK, "bin/ablation")
GEN = os.path.join(ROOT, "benchmarks/gen_pairs.py")

def ensure_pairs(length, div, seeds):
    pairs = []
    for seed in seeds:
        fa = os.path.join(WORK, "pairs", f"{length}_{div}_s{seed}.fa")
        subprocess.check_call(["python3", GEN, str(length), div, str(seed), fa])
        pairs.append(fa)
    return pairs

def run_once(cmd, require_score=True):
    """(wall ms, CPU ms, RSS MiB, score, output) of one run"""
    return timed(cmd, require_score)

SEEDS = [1, 2, 3, 4, 5]

if SECTION in ("all", "ablation"):
    print("\n=== FRESH-PROCESS 9-CONFIGURATION ABLATION (20 kb, 5 seeds, lib wall ms + peak RSS MiB) ===")
    ABL_MODES = [
        "scalar_3mat_1t", "striped_32b_1t", "striped_diff8_1t",
        "banded_nowfa_1t", "banded_wfa_1t", "banded_nt",
        "linear_1t", "linear_nt", "wfa_1t"
    ]
    for div in ["1", "10"]:
        pairs = ensure_pairs(20000, div, SEEDS)
        print(f"-- 20 kb {div}% --")
        for s, fa in zip(SEEDS, pairs):
            print(f"# sha256 20000_{div}_s{s}.fa = {sha256(fa)}")
        ref_sc = {}
        for m in ABL_MODES:
            lib_ms, rsss = [], []
            for idx, fa in enumerate(pairs):
                _, _, rss, _, out = run_once([ABL, fa, m, THREADS], require_score=False)
                parts = out.strip().splitlines()[-1].split("\t")
                ms, sc = float(parts[1]), int(parts[2])
                if idx not in ref_sc:
                    ref_sc[idx] = sc
                else:
                    assert sc == ref_sc[idx], f"Ablation score mismatch {m} seed {idx+1}: {sc} vs {ref_sc[idx]}"
                lib_ms.append(ms); rsss.append(rss)
            raw = ",".join(f"{x:.2f}" for x in lib_ms)
            print(f"  {m:20s}  med={statistics.median(lib_ms):8.2f} ms [{min(lib_ms):7.2f}-{max(lib_ms):7.2f}]  rss={statistics.median(rsss):7.1f} MiB  seeds=[{raw}]")

if SECTION in ("all", "global"):
    print("\n=== GLOBAL AFFINE-GAP COMPARATORS (+1/-2/-4/-1, 5 seeds) ===")
    for length in [20000, 50000]:
        for div in ["1", "10", "rand"]:
            pairs = ensure_pairs(length, div, SEEDS)
            print(f"-- Global {length//1000} kb {div} --")
            for s, fa in zip(SEEDS, pairs):
                print(f"# sha256 {length}_{div}_s{s}.fa = {sha256(fa)}")
            ref = {}
            tools = [
                ("orig_1t",     lambda fa: [ORIG_NW, "--printscores", "--file", fa]),
                ("nw_1t",       lambda fa: [NW, "--threads", "1", "--printscores", "--file", fa]),
                (f"nw_{THREADS}t", lambda fa: [NW, "--threads", THREADS, "--printscores", "--file", fa]),
                ("nw_lin_1t",   lambda fa: [NW, "--threads", "1", "--linear", "--printscores", "--file", fa]),
                (f"nw_lin_{THREADS}t", lambda fa: [NW, "--threads", THREADS, "--linear", "--printscores", "--file", fa]),
                ("nw_wfa_1t",   lambda fa: [NW, "--threads", "1", "--wfa", "--printscores", "--file", fa]),
                ("nw_score_1t", lambda fa: [NW, "--threads", "1", "--scoreonly", "--file", fa]),
                ("ksw2_1t",     lambda fa: [KSW2, fa]),
                ("wfa2_high",   lambda fa: [WFA2, "high", fa]),
                ("biwfa",       lambda fa: [WFA2, "biwfa", fa]),
            ]
            for tname, fn in tools:
                if length == 50000 and tname == "orig_1t":
                    continue
                if length == 50000 and div == "rand" and tname in ("wfa2_high", "biwfa", "nw_wfa_1t"):
                    continue
                if length == 50000 and div == "10" and tname == "nw_wfa_1t":
                    continue
                use_pairs = pairs[:1] if (tname == "orig_1t" or (tname == "wfa2_high" and div == "rand")) else pairs
                walls, cpus, rsss = [], [], []
                for idx, fa in enumerate(use_pairs):
                    w, c, r, sc, _ = run_once(fn(fa), require_score=True)
                    if idx not in ref:
                        ref[idx] = sc
                    else:
                        assert sc == ref[idx], f"Global score mismatch {tname} {length} {div} seed {idx+1}: {sc} vs {ref[idx]}"
                    walls.append(w); cpus.append(c); rsss.append(r)
                raw = ",".join(f"{x:.2f}" for x in walls)
                print(f"  {tname:12s}  wall={statistics.median(walls):8.2f} ms [{min(walls):7.2f}-{max(walls):7.2f}]  cpu={statistics.median(cpus):8.1f} ms  rss={statistics.median(rsss):7.1f} MiB  seeds=[{raw}]")

if SECTION in ("all", "edlib"):
    print("\n=== UNIT-COST EDIT DISTANCE (--scoreonly vs edlib v1.2.7, 1 thread, 5 seeds) ===")
    for length, div in [(20000, "1"), (20000, "10"), (20000, "rand"),
                        (50000, "1"), (50000, "10"), (50000, "rand")]:
        pairs = ensure_pairs(length, div, SEEDS)
        print(f"-- Unit-cost {length//1000} kb {div} --")
        ref = {}
        for tname, fn in [
            ("nw_score_1t", lambda fa: [NW, "--threads", "1", "--scoreonly", "--match", "0", "--mismatch", "-1", "--gapopen", "0", "--gapextend", "-1", "--file", fa]),
            ("edlib_1t",    lambda fa: [EDLIB, fa]),
        ]:
            walls, rsss = [], []
            for idx, fa in enumerate(pairs):
                w, _, r, sc, _ = run_once(fn(fa), require_score=True)
                if idx not in ref:
                    ref[idx] = sc
                else:
                    assert sc == ref[idx], f"Unit-cost mismatch {tname} seed {idx+1}: {sc} vs {ref[idx]}"
                walls.append(w); rsss.append(r)
            raw = ",".join(f"{x:.2f}" for x in walls)
            print(f"  {tname:12s}  wall={statistics.median(walls):7.2f} ms [{min(walls):6.2f}-{max(walls):6.2f}]  rss={statistics.median(rsss):5.1f} MiB  seeds=[{raw}]")

if SECTION in ("all", "local"):
    print("\n=== LOCAL 20 kb (--maxhits 1 vs SSW, +2/-2/-2/-1, 5 seeds, per-seed scores & SSW flags) ===")
    for div in ["1", "10", "rand"]:
        pairs = ensure_pairs(20000, div, SEEDS)
        print(f"-- Local 20 kb {div} --")
        sw_scores = {}
        for tname, fn in [
            ("sw_1t", lambda fa: [SW, "--threads", "1", "--maxhits", "1", "--minscore", "1", "--file", fa]),
            (f"sw_{THREADS}t", lambda fa: [SW, "--threads", THREADS, "--maxhits", "1", "--minscore", "1", "--file", fa]),
            ("ssw", lambda fa: [SSW, fa]),
        ]:
            walls, cpus, rsss, scores, flags = [], [], [], [], []
            for idx, fa in enumerate(pairs):
                w, c, r, sc, out = run_once(fn(fa), require_score=True)
                if tname == "sw_1t":
                    sw_scores[idx] = sc
                elif tname == f"sw_{THREADS}t":
                    assert sc == sw_scores[idx], f"Local score mismatch {tname} seed {idx+1}"
                elif tname == "ssw":
                    mf = re.search(r"ssw_flag:\s*(\d+)", out)
                    flags.append(int(mf.group(1)) if mf else -1)
                    if div != "1":
                        assert sc == sw_scores[idx], f"SSW score mismatch on {div} seed {idx+1}: {sc} vs {sw_scores[idx]}"
                walls.append(w); cpus.append(c); rsss.append(r); scores.append(sc)
            extra = f"  flags={flags}" if flags else ""
            raw = ",".join(f"{x:.2f}" for x in walls)
            print(f"  {tname:12s}  wall={statistics.median(walls):8.2f} ms [{min(walls):7.2f}-{max(walls):7.2f}]  cpu={statistics.median(cpus):8.1f} ms  rss={statistics.median(rsss):7.1f} MiB  scores={scores}{extra}  seeds=[{raw}]")
PY
