"""benchmarks/bench_common.py: what the benchmark scripts share.

Timing a command under GNU time (wall, CPU and peak memory, and the score
it prints), summarising repeated runs, file digests, FASTA reading, and
re-scoring alignments and CIGAR strings independently of the aligners.
The scripts here import it directly; the Python embedded in the shell
harnesses finds it through PYTHONPATH, which they set to this directory.
"""
import hashlib
import re
import statistics
import subprocess
import time
from collections import namedtuple

# One run of a command: wall and CPU time in ms, peak RSS in MiB, the last
# score it printed (None if none) and its standard output
Run = namedtuple("Run", "wall_ms cpu_ms rss_mib score stdout")

SCORE_RE = re.compile(r"(?:^|\s)score:\s*(-?\d+)")
_RSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
_USER_RE = re.compile(r"User time \(seconds\):\s*([\d.]+)")
_SYS_RE = re.compile(r"System time \(seconds\):\s*([\d.]+)")


def timed(cmd, require_score=False):
    """Run cmd once under /usr/bin/time -v. Raises RuntimeError if it fails,
    if GNU time's report cannot be read, or (require_score) if it prints no
    score."""
    t0 = time.perf_counter()
    p = subprocess.run(["/usr/bin/time", "-v"] + list(cmd),
                       capture_output=True, text=True)
    wall = (time.perf_counter() - t0) * 1000.0
    if p.returncode != 0:
        # The command's own messages come first; GNU time's report follows
        raise RuntimeError(f"Failed ({p.returncode}): {cmd}\n{p.stderr[:400]}")
    m_rss, m_usr, m_sys = (r.search(p.stderr) for r in (_RSS_RE, _USER_RE, _SYS_RE))
    if not m_rss or not m_usr or not m_sys:
        raise RuntimeError(f"Failed to parse GNU time output:\n{p.stderr}")
    scores = SCORE_RE.findall(p.stdout)
    if require_score and not scores:
        raise RuntimeError(f"Missing score in output for {cmd}:\n{p.stdout[:200]}")
    return Run(wall,
               (float(m_usr.group(1)) + float(m_sys.group(1))) * 1000.0,
               int(m_rss.group(1)) / 1024.0,
               int(scores[-1]) if scores else None,
               p.stdout)


def timed_reps(cmd, reps, require_score=False):
    """reps runs of cmd (see timed())."""
    return [timed(cmd, require_score) for _ in range(reps)]


def median_summary(runs):
    """(median wall, min wall, max wall, median CPU, median RSS, the last
    run's score, the last run's output)."""
    walls = [r.wall_ms for r in runs]
    return (statistics.median(walls), min(walls), max(walls),
            statistics.median(r.cpu_ms for r in runs),
            statistics.median(r.rss_mib for r in runs),
            runs[-1].score, runs[-1].stdout)


def best_run(runs):
    """The run with the least wall time (the first of equals)."""
    return min(runs, key=lambda r: r.wall_ms)


def sha256(path):
    """The SHA-256 digest of a file, in hex."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def seq_len(path):
    """Letters in a FASTA file (every line but the headers)."""
    with open(path) as f:
        return sum(len(line.strip()) for line in f if not line.startswith(">"))


def read_seqs(path):
    """The sequences of a FASTA file, in order."""
    seqs, cur = [], []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def rescore_aligned(a, b, sub_fn, go, ge):
    """The score of the aligned strings a and b ('-' for a gap): sub_fn(x, y)
    for each pair (upper case), and go + ge to open a gap, ge to extend it."""
    assert len(a) == len(b)
    sc, in_a, in_b = 0, False, False
    for x, y in zip(a, b):
        assert not (x == "-" and y == "-"), "Invalid double-gap column"
        if x == "-":
            sc += ge if in_a else (go + ge)
            in_a, in_b = True, False
        elif y == "-":
            sc += ge if in_b else (go + ge)
            in_b, in_a = True, False
        else:
            in_a = in_b = False
            sc += sub_fn(x.upper(), y.upper())
    return sc


def rescore_cigar(seq_a, seq_b, cigar_str, sub_fn, go, ge, swap_id=False):
    """The score of seq_a against seq_b by the CIGAR string (an operation
    per letter, or run-length), scored as rescore_aligned(). I consumes
    seq_b and D seq_a (swapped with swap_id); the whole of both must be
    consumed."""
    i = j = sc = 0
    if re.match(r"^[MIDX=]+$", cigar_str):
        ops = [(1, ch) for ch in cigar_str]
    else:
        assert re.fullmatch(r"(?:[1-9]\d*[MIDX=])+", cigar_str), \
            f"Malformed CIGAR: {cigar_str[:40]}"
        ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDX=])", cigar_str)]
    prev_op = None
    for n, op in ops:
        if op in ("M", "=", "X"):
            for _ in range(n):
                sc += sub_fn(seq_a[i].upper(), seq_b[j].upper())
                i += 1
                j += 1
            prev_op = "M"
        elif (op == "I" and not swap_id) or (op == "D" and swap_id):
            sc += (ge * n) if prev_op == "I" else (go + ge * n)
            j += n
            prev_op = "I"
        else:
            sc += (ge * n) if prev_op == "D" else (go + ge * n)
            i += n
            prev_op = "D"
    assert i == len(seq_a) and j == len(seq_b), \
        f"Span mismatch: ({i},{j}) vs ({len(seq_a)},{len(seq_b)})"
    return sc
