#!/usr/bin/env bash
# Reproducible multi-pair protein benchmark (BLOSUM62, -10 / -1) across five
# real NCBI RefSeq orthologue pairs spanning ~50% to ~99% sequence identity
# and ~390 to ~3,060 amino acids. Compares original seq-align (dc41988) against
# Nucleics Align (full 1t/8t, --linear 1t/8t, --scoreonly, local --maxhits 1).
#
# Usage: bash benchmarks/protein_bench.sh [work_dir] [threads]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:-/tmp/nucleics_compare}"
THREADS="${2:-8}"
OUT="$WORK/proteins"
mkdir -p "$OUT"

NW="$ROOT/bin/needleman_wunsch"
SW="$ROOT/bin/smith_waterman"
ORIG_NW="$WORK/src/seq-align/bin/needleman_wunsch"
. "$ROOT/benchmarks/bench_common.sh"
bench_require
if [ ! -x "$ORIG_NW" ]; then
  echo "Run benchmarks/compare_harness.sh first (or build seq-align at $ORIG_NW)" >&2
  exit 1
fi

fetch_prot() { fetch_ncbi protein "$1" "$OUT/$2"; }

fetch_prot NP_000537.3 h_tp53.fa
fetch_prot NP_035770.2 m_tp53.fa
fetch_prot NP_000928.1 h_polr2a.fa
fetch_prot NP_010141.1 y_rpo21.fa
fetch_prot NP_009225.1 h_brca1.fa
fetch_prot NP_033894.3 m_brca1.fa
fetch_prot NP_004949.1 h_mtor.fa
fetch_prot NP_064393.2 m_mtor.fa
fetch_prot NP_000042.3 h_atm.fa
fetch_prot NP_031525.2 m_atm.fa

python3 - "$NW" "$SW" "$ORIG_NW" "$OUT" "$THREADS" <<'PY'
import hashlib, os, re, statistics, subprocess, sys, time
from bench_common import sha256, seq_len, median_summary, timed_reps

NW, SW, ORIG_NW, OUT, THREADS = sys.argv[1:6]


def run_best(cmd, reps=3):
    """(median wall, min, max, median RSS, score) of reps runs"""
    med, lo, hi, _, rss, sc, _ = median_summary(timed_reps(cmd, reps, require_score=True))
    return med, lo, hi, rss, sc

def identity_pct(fa1, fa2):
    p = subprocess.run([NW, "--threads", THREADS, "--scoring", "BLOSUM62", "--pretty", "--files", fa1, fa2],
                       capture_output=True, text=True, check=True)
    lines = [l for l in p.stdout.splitlines() if l.strip()]
    m = lines[1]
    return 100.0 * m.count("|") / len(m)

pairs = [
    ("Human/Mouse TP53", "h_tp53.fa", "m_tp53.fa"),
    ("Human/Yeast POLR2A", "h_polr2a.fa", "y_rpo21.fa"),
    ("Human/Mouse BRCA1", "h_brca1.fa", "m_brca1.fa"),
    ("Human/Mouse MTOR", "h_mtor.fa", "m_mtor.fa"),
    ("Human/Mouse ATM", "h_atm.fa", "m_atm.fa"),
]

print("\n=== PROTEIN BENCHMARK (BLOSUM62 -10/-1, median [min-max] of 3 runs) ===")
for name, f1, f2 in pairs:
    p1, p2 = os.path.join(OUT, f1), os.path.join(OUT, f2)
    la, lb = seq_len(p1), seq_len(p2)
    idp = identity_pct(p1, p2)
    print(f"\n-- {name} ({la} / {lb} aa, id={idp:.2f}%) --")
    print(f"# sha256: {f1}={sha256(p1)} {f2}={sha256(p2)}")
    g_scores = set()
    for mname, cmd in [
        ("seqalign_1t", [ORIG_NW, "--scoring", "BLOSUM62", "--printscores", "--files", p1, p2]),
        ("full_1t",     [NW, "--threads", "1", "--scoring", "BLOSUM62", "--printscores", "--files", p1, p2]),
        (f"full_{THREADS}t", [NW, "--threads", THREADS, "--scoring", "BLOSUM62", "--printscores", "--files", p1, p2]),
        ("lin_1t",      [NW, "--threads", "1", "--linear", "--scoring", "BLOSUM62", "--printscores", "--files", p1, p2]),
        (f"lin_{THREADS}t",  [NW, "--threads", THREADS, "--linear", "--scoring", "BLOSUM62", "--printscores", "--files", p1, p2]),
        ("scoreonly",   [NW, "--threads", THREADS, "--scoreonly", "--scoring", "BLOSUM62", "--files", p1, p2]),
        ("local1",      [SW, "--threads", THREADS, "--maxhits", "1", "--scoring", "BLOSUM62", "--files", p1, p2]),
    ]:
        med, lo, hi, rss, sc = run_best(cmd, reps=3)
        if mname != "local1":
            g_scores.add(sc)
        print(f"  {mname:12s}  wall={med:8.2f} ms [{lo:7.2f}-{hi:7.2f}]  rss={rss:6.1f} MiB  score={sc}")
    assert len(g_scores) == 1, f"Protein score mismatch for {name}: {g_scores}"
PY
