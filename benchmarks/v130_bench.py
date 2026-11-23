#!/usr/bin/env python3
"""Targeted v1.3.0 benchmark on pop-os comparing v1.3.0 against v1.2.0 baselines:
1. BiWFA (--wfa, --wfabudget) on 20kb 1%, 50kb 1%, 100kb 1%, 300kb 1%, 1Mb 0.1%, 1Mb 1%
2. Banded Myers (--scoreonly) & Myers-Hirschberg (--linear) on 20kb/50kb unit-cost (0/-1/0/-1)
3. Opt-in Smith-Waterman --fastlocal --maxhits 1 on 20kb 1%, 10%, unrelated (1t & 8t)
4. 8-thread full matrix wall + CPU time on 50kb unrelated
5. Short protein (TP53 393 aa) wall time
"""
import os, re, statistics, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NW = os.path.join(ROOT, "bin", "needleman_wunsch")
SW = os.path.join(ROOT, "bin", "smith_waterman")
TMP = "/tmp/na_v130_bench"
os.makedirs(TMP, exist_ok=True)

def gen_pair(length, rate, seed, tag):
    import random
    fa1 = os.path.join(TMP, f"{tag}_a.fa")
    fa2 = os.path.join(TMP, f"{tag}_b.fa")
    if os.path.exists(fa1) and os.path.exists(fa2):
        return fa1, fa2
    rng = random.Random(seed)
    bases = "ACGT"
    seqa = [rng.choice(bases) for _ in range(length)]
    if rate is None:
        seqb = [rng.choice(bases) for _ in range(length)]
    else:
        seqb = []
        i = 0
        while i < len(seqa):
            if rng.random() < rate:
                r = rng.random()
                if r < 0.7:
                    seqb.append(rng.choice([b for b in bases if b != seqa[i]]))
                    i += 1
                elif r < 0.85:
                    i += 1  # deletion
                else:
                    seqb.append(rng.choice(bases))  # insertion
            else:
                seqb.append(seqa[i])
                i += 1
    with open(fa1, "w") as f:
        f.write(f">a\n{''.join(seqa)}\n")
    with open(fa2, "w") as f:
        f.write(f">b\n{''.join(seqb)}\n")
    return fa1, fa2

def run_timed(cmd, reps=3):
    walls, cpus, rsss, scores = [], [], [], []
    for _ in range(reps):
        t0 = time.perf_counter()
        p = subprocess.run(["/usr/bin/time", "-v"] + cmd,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        t1 = time.perf_counter()
        if p.returncode != 0:
            return None, None, None, f"ERR({p.returncode}): {p.stderr[-200:]}"
        walls.append((t1 - t0) * 1000.0)
        u = re.search(r"User time \(seconds\):\s*([0-9.]+)", p.stderr)
        s = re.search(r"System time \(seconds\):\s*([0-9.]+)", p.stderr)
        m = re.search(r"Maximum resident set size \(kbytes\):\s*([0-9]+)", p.stderr)
        cpus.append(((float(u.group(1)) if u else 0.0) + (float(s.group(1)) if s else 0.0)) * 1000.0)
        rsss.append((int(m.group(1)) if m else 0) / 1024.0)
        sc = re.search(r"score:\s*(-?\d+)", p.stdout)
        scores.append(int(sc.group(1)) if sc else None)
    return statistics.median(walls), statistics.median(cpus), statistics.median(rsss), scores[0]

print("=== 1. BiWFA (--wfa) ===")
for name, length, rate, extra in [
    ("20kb 1%",  20000,   0.01,  []),
    ("50kb 1%",  50000,   0.01,  []),
    ("100kb 1%", 100000,  0.01,  []),
    ("300kb 1%", 300000,  0.01,  ["--wfabudget", "1000000000"]),
    ("1Mb 0.1%", 1000000, 0.001, []),
]:
    a, b = gen_pair(length, rate, 11, name.replace(" ", "_").replace("%", "p"))
    w, c, r, sc = run_timed([NW, "--threads", "1", "--wfa", "--printscores"] + extra + ["--files", a, b], reps=3)
    print(f"  {name:10s}: wall={w:7.1f} ms  cpu={c:7.1f} ms  rss={r:6.1f} MiB  score={sc}")

print("\n=== 2. Unit-cost Myers (0/-1/0/-1) ===")
for name, length, rate in [
    ("20kb 1%",   20000, 0.01),
    ("20kb 10%",  20000, 0.10),
    ("20kb unrl", 20000, None),
    ("50kb 1%",   50000, 0.01),
    ("50kb 10%",  50000, 0.10),
    ("50kb unrl", 50000, None),
]:
    a, b = gen_pair(length, rate, 11, "u_" + name.replace(" ", "_").replace("%", "p"))
    w1, _, r1, sc1 = run_timed([NW, "--threads", "1", "--match", "0", "--mismatch", "-1",
                                "--gapopen", "0", "--gapextend", "-1", "--scoreonly", "--files", a, b], reps=3)
    w2, _, r2, sc2 = run_timed([NW, "--threads", "1", "--match", "0", "--mismatch", "-1",
                                "--gapopen", "0", "--gapextend", "-1", "--linear", "--printscores", "--files", a, b], reps=3)
    print(f"  {name:10s}: --scoreonly={w1:6.1f} ms ({r1:4.1f} MiB, sc={sc1})  --linear={w2:6.1f} ms ({r2:4.1f} MiB, sc={sc2})")

print("\n=== 3. Smith-Waterman --maxhits 1 --minscore 1 (default vs --fastlocal) ===")
for name, length, rate in [
    ("20kb 1%",   20000, 0.01),
    ("20kb 10%",  20000, 0.10),
    ("20kb unrl", 20000, None),
]:
    a, b = gen_pair(length, rate, 11, "sw_" + name.replace(" ", "_").replace("%", "p"))
    wd, _, rd, scd = run_timed([SW, "--threads", "1", "--maxhits", "1", "--minscore", "1", "--files", a, b], reps=3)
    wf, _, rf, scf = run_timed([SW, "--threads", "1", "--maxhits", "1", "--minscore", "1", "--fastlocal", "--files", a, b], reps=3)
    print(f"  {name:10s}: default_1t={wd:6.1f} ms ({rd:5.1f} MiB, sc={scd})  fastlocal_1t={wf:6.1f} ms ({rf:5.1f} MiB, sc={scf})")

print("\n=== 4. 8-thread full matrix (50kb unrelated) ===")
a, b = gen_pair(50000, None, 11, "u_50kb_unrl")
w8, c8, r8, sc8 = run_timed([NW, "--threads", "8", "--printscores", "--files", a, b], reps=3)
print(f"  50kb unrl 8t: wall={w8:6.1f} ms  cpu={c8:7.1f} ms  rss={r8:6.1f} MiB  score={sc8}")
