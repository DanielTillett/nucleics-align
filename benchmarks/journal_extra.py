#!/usr/bin/env python3
"""Benchmark parasail scan (nw_trace_scan_avx2_256_16) and deep BiWFA CIGAR re-scoring."""
import re, statistics, subprocess, time

WORK = "/tmp/nucleics_compare"
NW = "./bin/needleman_wunsch"
SCAN = WORK + "/bin/run_parasail_scan"
WFA2_CIG = WORK + "/bin/run_wfa2_cigar"

def run3(cmd):
    ws, rs = [], []
    sc = None
    for _ in range(3):
        t0 = time.perf_counter()
        p = subprocess.run(["/usr/bin/time", "-v"] + cmd, capture_output=True, text=True, check=True)
        ws.append((time.perf_counter() - t0) * 1000.0)
        rs.append(int(re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", p.stderr).group(1)) / 1024.0)
        sc = int(re.findall(r"score:\s*(-?\d+)", p.stdout)[-1])
    return statistics.median(ws), min(ws), max(ws), statistics.median(rs), sc

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

def read_seqs(p):
    s, c = [], []
    for l in open(p):
        if l.startswith(">"):
            if c:
                s.append("".join(c))
                c = []
        else:
            c.append(l.strip())
    if c:
        s.append("".join(c))
    return s

def rescore_aln(a, b):
    assert len(a) == len(b)
    sc = 0
    ia = ib = False
    for x, y in zip(a, b):
        assert not (x == "-" and y == "-")
        if x == "-":
            sc += -1 if ia else -5
            ia, ib = True, False
        elif y == "-":
            sc += -1 if ib else -5
            ib, ia = True, False
        else:
            ia = ib = False
            sc += 1 if x.upper() == y.upper() else -2
    return sc

def rescore_wfa_cig(sa, sb, cig):
    i = j = sc = 0
    prev = None
    for op in cig:
        if op in ("M", "X", "="):
            sc += 1 if sa[i].upper() == sb[j].upper() else -2
            i += 1
            j += 1
            prev = "M"
        elif op == "I":
            sc += -1 if prev == "I" else -5
            j += 1
            prev = "I"
        else:
            sc += -1 if prev == "D" else -5
            i += 1
            prev = "D"
    assert i == len(sa) and j == len(sb)
    return sc

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
