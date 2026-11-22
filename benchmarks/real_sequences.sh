#!/usr/bin/env bash
# Reproducible real-sequence benchmark for Table 1 of the Nucleics Align paper.
# Downloads seven public NCBI RefSeq/GenBank pairs, computes SHA-256 hashes,
# times needleman_wunsch (full, --wfa, --linear, --scoreonly) and
# smith_waterman (--maxhits 1) with pinned --threads 8, asserts identical
# optimal global scores across all global modes, and reports full-alignment % identity.
#
# Usage: bash benchmarks/real_sequences.sh [output_dir] [threads]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$ROOT/benchmarks/real_data}"
THREADS="${2:-8}"
mkdir -p "$OUT"

NW="$ROOT/bin/needleman_wunsch"
SW="$ROOT/bin/smith_waterman"
if [ ! -x "$NW" ] || [ ! -x "$SW" ]; then
  echo "Build the tools first with 'make' in $ROOT" >&2
  exit 1
fi
if [ ! -x /usr/bin/time ]; then
  echo "GNU time (/usr/bin/time) is required for CPU time and peak RSS" >&2
  exit 1
fi

echo "# host: $(uname -a)"
echo "# git:  $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "# threads: $THREADS"

fetch_acc() {
  local acc="$1" db="$2" dst="$OUT/$3"
  if [ ! -s "$dst" ]; then
    curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${db}&id=${acc}&rettype=fasta&retmode=text" -o "$dst"
    sleep 0.4
  fi
}

fetch_acc NC_012920.1 nuccore human_mt.fa
fetch_acc NC_011137.1 nuccore neanderthal_mt.fa
fetch_acc NC_001643.1 nuccore chimp_mt.fa
fetch_acc NC_045512.2 nuccore sars2_wuhan.fa
fetch_acc OL672836.1  nuccore sars2_omicron.fa
fetch_acc NC_004718.3 nuccore sars1.fa
fetch_acc NC_063383.1 nuccore mpox_ref.fa
fetch_acc ON563414.3  nuccore mpox_2022.fa
fetch_acc NC_000908.2 nuccore mgen_g37.fa
fetch_acc CP003773.1  nuccore mgen_m2288.fa
fetch_acc NP_009225.1 protein human_brca1.fa
fetch_acc NP_033894.3 protein mouse_brca1.fa

python3 - "$NW" "$SW" "$OUT" "$THREADS" <<'PY'
import hashlib, os, re, sys, subprocess, time

NW, SW, OUT, THREADS = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()[:12]

def seq_len(fa):
    with open(fa) as f:
        return sum(len(l.strip()) for l in f if not l.startswith(">"))

def run_cmd(cmd, reps=3):
    best = None
    for _ in range(reps):
        t0 = time.perf_counter()
        p = subprocess.run(["/usr/bin/time", "-v"] + cmd, capture_output=True, text=True)
        dt = (time.perf_counter() - t0) * 1000.0
        if p.returncode != 0:
            raise RuntimeError(f"Command failed ({p.returncode}): {cmd}\n{p.stderr}")
        m_rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", p.stderr)
        m_usr = re.search(r"User time \(seconds\):\s*([\d.]+)", p.stderr)
        m_sys = re.search(r"System time \(seconds\):\s*([\d.]+)", p.stderr)
        if not m_rss or not m_usr or not m_sys:
            raise RuntimeError(f"Failed to parse GNU time output:\n{p.stderr}")
        rss_mib = int(m_rss.group(1)) / 1024.0
        cpu_ms = (float(m_usr.group(1)) + float(m_sys.group(1))) * 1000.0
        m_sc = re.findall(r"(?:^|\s)score:\s*(-?\d+)", p.stdout)
        if not m_sc:
            raise RuntimeError(f"Failed to parse score from output for {cmd}")
        score = int(m_sc[-1])
        cand = (dt, cpu_ms, rss_mib, score)
        if best is None or dt < best[0]:
            best = cand
    return best

def identity_pct(fa1, fa2, extra):
    p = subprocess.run([NW, "--threads", THREADS, "--pretty", "--files", fa1, fa2] + extra,
                       capture_output=True, text=True, check=True)
    lines = [l for l in p.stdout.splitlines() if l.strip()]
    a, m, b = lines[0], lines[1], lines[2]
    mat = m.count("|")
    return 100.0 * mat / len(m), mat, m.count("*"), a.count("-") + b.count("-")

pairs = [
    ("Human vs Neanderthal mtDNA", "human_mt.fa", "neanderthal_mt.fa", [], 3),
    ("Human vs Chimp mtDNA", "human_mt.fa", "chimp_mt.fa", [], 3),
    ("SARS-CoV-2 Wuhan vs Omicron", "sars2_wuhan.fa", "sars2_omicron.fa", [], 3),
    ("SARS-CoV-2 vs SARS-CoV-1", "sars2_wuhan.fa", "sars1.fa", [], 3),
    ("Mpox Ref vs 2022 isolate", "mpox_ref.fa", "mpox_2022.fa", [], 3),
    ("M. genitalium G37 vs M2288", "mgen_g37.fa", "mgen_m2288.fa", [], 1),
    ("Human vs Mouse BRCA1 (BLOSUM62)", "human_brca1.fa", "mouse_brca1.fa", ["--scoring", "BLOSUM62"], 3),
]

print(f"{'Pair':34s} {'LenA':>7s} {'LenB':>7s} {'Mode':11s} {'Wall_ms':>9s} {'CPU_ms':>9s} {'RSS_MiB':>8s} {'Score':>10s}")
for name, f1, f2, extra, reps in pairs:
    p1, p2 = os.path.join(OUT, f1), os.path.join(OUT, f2)
    la, lb = seq_len(p1), seq_len(p2)
    print(f"# {name}: {f1}({sha256(p1)}) vs {f2}({sha256(p2)})")
    modes = [
        ("full",      [NW, "--threads", THREADS, "--printscores", "--files", p1, p2] + extra),
        ("wfa",       [NW, "--threads", THREADS, "--wfa", "--printscores", "--files", p1, p2] + extra),
        ("linear",    [NW, "--threads", THREADS, "--linear", "--printscores", "--files", p1, p2] + extra),
        ("scoreonly", [NW, "--threads", THREADS, "--scoreonly", "--files", p1, p2] + extra),
        ("local1",    [SW, "--threads", THREADS, "--maxhits", "1", "--files", p1, p2] + extra),
    ]
    global_scores = set()
    for mname, cmd in modes:
        if "M. genitalium" in name and mname == "local1":
            continue
        w, c, r, sc = run_cmd(cmd, reps=reps)
        if mname != "local1":
            global_scores.add(sc)
        print(f"{name:34s} {la:7d} {lb:7d} {mname:11s} {w:9.2f} {c:9.1f} {r:8.1f} {sc:10d}")
    assert len(global_scores) == 1, f"Score mismatch for {name}: {global_scores}"
    pct, mat, mis, gap = identity_pct(p1, p2, extra)
    print(f"  IDENTITY: {pct:.2f}% (matches={mat}, mismatches={mis}, gap_cols={gap})")
PY
