#!/usr/bin/env python3
"""Benchmark parasail scan (nw_trace_scan_avx2_256_16) and deep BiWFA CIGAR re-scoring."""
import re, statistics, subprocess, sys, time

from bench_common import (read_seqs, rescore_aligned, rescore_cigar,
                          median_summary, timed_reps)

# The work directory of compare_harness.sh and journal_bench.sh
WORK = sys.argv[1] if len(sys.argv) > 1 else "/tmp/nucleics_compare"
NW = "./bin/needleman_wunsch"
SCAN = WORK + "/bin/run_parasail_scan"
WFA2_CIG = WORK + "/bin/run_wfa2_cigar"

def run3(cmd):
    """(median wall, min, max, median RSS, score) of three runs"""
    med, lo, hi, _, rss, sc, _ = median_summary(timed_reps(cmd, 3))
    return med, lo, hi, rss, sc

print("=== PARASAIL v2.6.2 SCAN (nw_trace_scan_avx2_256_16, 1t) ===")
for name, f1, f2 in [
    ("TP53", "h_tp53.fa", "m_tp53.fa"),
    ("POLR2A/Rpo21", "h_polr2a.fa", "y_rpo21.fa"),
    ("BRCA1", "h_brca1.fa", "m_brca1.fa"),
    ("MTOR", "h_mtor.fa", "m_mtor.fa"),
    ("ATM", "h_atm.fa", "m_atm.fa"),
]:
    med, lo, hi, rss, sc = run3([SCAN, WORK + "/proteins/" + f1, WORK + "/proteins/" + f2])
    print(f"  {name:14s}  scan_1t={med:7.2f} ms [{lo:6.2f}-{hi:6.2f}]  rss={rss:5.1f} MiB  score={sc}")


# The DNA scoring of these runs: match 1, mismatch -2, gap -4 - 1 per letter
def dna(x, y):
    return 1 if x == y else -2


def rescore_aln(a, b):
    return rescore_aligned(a, b, dna, -4, -1)

def rescore_wfa_cig(sa, sb, cig):
    return rescore_cigar(sa, sb, cig, dna, -4, -1)

print("\n=== DEEP BiWFA INDEPENDENT CIGAR/ALIGNMENT RE-SCORING (100 kb 1% & 1 Mb 0.1%) ===")
for label, fa in [("100 kb 1%", WORK + "/scale/100000_1_s11.fa"), ("1 Mb 0.1%", WORK + "/scale/1000000_0.1_s11.fa")]:
    sa, sb = read_seqs(fa)
    out = subprocess.check_output([NW, "--threads", "1", "--wfa", "--printscores", "--file", fa], text=True)
    lines = [l for l in out.splitlines() if l.strip()]
    aa, ab = lines[0], lines[1]
    assert aa.replace("-", "") == sa and ab.replace("-", "") == sb
    rep_n = int(re.search(r"score:\s*(-?\d+)", out).group(1))
    rsc_n = rescore_aln(aa, ab)
    assert rsc_n == rep_n
    wout = subprocess.check_output([WFA2_CIG, fa], text=True)
    rep_w = int(re.search(r"score:\s*(-?\d+)", wout).group(1))
    cig_w = re.search(r"cigar:\s*(\S+)", wout).group(1)
    rsc_w = rescore_wfa_cig(sa, sb, cig_w)
    assert rsc_w == rep_w == rsc_n
    print(f"  {label:10s}  BiWFA & WFA2_BiWFA RE-SCORED & SPAN-VERIFIED: score={rsc_n}")
