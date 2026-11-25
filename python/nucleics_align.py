"""Python interface to the Nucleics Align command line tools.

python/nucleics_align.py
url: https://github.com/DanielTillett/nucleics-align
license: BSD 3-Clause, see LICENSE

Runs bin/needleman_wunsch and bin/smith_waterman, writing the pairs to them
as FASTA on a background thread and reading each alignment back as it is
printed, so a batch of pairs runs on the tools' threads in one process.
Only the Python standard library is needed.

    from nucleics_align import NeedlemanWunsch, SmithWaterman

    nw = NeedlemanWunsch()
    aln = nw.align("ACGTTGCA", "ACGTGCA")
    print(aln.align_a, aln.markers, aln.align_b, aln.score, sep="\n")

    for aln in nw.align_pairs(pairs):      # pairs: (seq_a, seq_b) tuples
        ...
    print(nw.score("ACGTTGCA", "ACGTGCA"))  # score only, linear memory

    sw = SmithWaterman(max_hits=3)
    for hit in sw.align("ttacgtacgtaa", "acgtacgt"):
        print(hit.score, hit.pos_a, hit.align_a, hit.align_b)

The tools are found in bin_dir if given, else in $NUCLEICS_ALIGN_BIN, else
in ../bin next to this file (a built checkout), else on the PATH.
"""

import os
import queue
import re
import shutil
import subprocess
import threading
from dataclasses import dataclass
from typing import Iterable, Iterator, List, Tuple

__all__ = ["NeedlemanWunsch", "SmithWaterman", "GlobalAlignment", "LocalHit",
           "AlignmentError", "SCORING_SYSTEMS"]

# Built-in scoring systems (--scoring)
SCORING_SYSTEMS = ("PAM30", "PAM70", "BLOSUM80", "BLOSUM62",
                   "DNA_HYBRIDIZATION")


class AlignmentError(RuntimeError):
    """A tool refused the options or failed on a pair: the message is the
    tool's own (its standard error)."""


@dataclass
class GlobalAlignment:
    """One Needleman-Wunsch alignment. markers has '|' for a match, '*' for
    a mismatch and ' ' for a gap, one per column."""
    seq_a: str
    seq_b: str
    align_a: str
    align_b: str
    markers: str
    score: int


@dataclass
class LocalHit:
    """One Smith-Waterman hit: the aligned parts, with '-' for gaps, and
    where they lie in each sequence (0-based start, length in letters)."""
    align_a: str
    align_b: str
    markers: str
    score: int
    pos_a: int
    len_a: int
    pos_b: int
    len_b: int


def _find_tool(name, bin_dir):
    exe = name + (".exe" if os.name == "nt" else "")
    dirs = [bin_dir, os.environ.get("NUCLEICS_ALIGN_BIN"),
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         os.pardir, "bin")]
    for d in dirs:
        if d and os.path.isfile(os.path.join(d, exe)):
            return os.path.normpath(os.path.join(d, exe))
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(
        "cannot find %s: build the tools (make) or set NUCLEICS_ALIGN_BIN"
        % name)


# Letters the tools could read as something else: whitespace and control
# characters (line breaks end a sequence), and '>' starting a sequence (a
# FASTA header)
_BAD_LETTER = re.compile(r"[\x00-\x20\x7f]")


def _check_sequence(seq, which):
    if not isinstance(seq, str):
        raise TypeError("%s must be a str, not %s" % (which, type(seq).__name__))
    bad = _BAD_LETTER.search(seq)
    if bad:
        raise ValueError("%s has a whitespace or control character (%r) at "
                         "position %d" % (which, bad.group(), bad.start()))
    if seq.startswith(">"):
        raise ValueError("%s cannot start with '>'" % which)


def _int_option(name, value):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("%s must be an int" % name)
    return str(value)


class _Tool:
    """What both aligners share: scoring options, and running the tool on a
    stream of pairs."""

    _NAME = ""

    def __init__(self, *, match=None, mismatch=None, gap_open=None,
                 gap_extend=None, scoring=None, substitution_matrix=None,
                 substitution_pairs=None, wildcards=None,
                 case_sensitive=False, no_gaps_in_a=False, no_gaps_in_b=False,
                 no_gaps=False, no_gaps_in_longest=False, no_mismatches=False,
                 threads=None, bin_dir=None):
        self._exe = _find_tool(self._NAME, bin_dir)
        args = []
        for flag, value in (("--match", match), ("--mismatch", mismatch),
                            ("--gapopen", gap_open),
                            ("--gapextend", gap_extend)):
            if value is not None:
                args += [flag, _int_option(flag[2:], value)]
        if scoring is not None:
            if str(scoring).upper() not in SCORING_SYSTEMS:
                raise ValueError("unknown scoring %r (one of %s)"
                                 % (scoring, ", ".join(SCORING_SYSTEMS)))
            args += ["--scoring", str(scoring).upper()]
        for flag, path in (("--substitution_matrix", substitution_matrix),
                           ("--substitution_pairs", substitution_pairs)):
            if path is not None:
                path = os.fspath(path)
                if not os.path.isfile(path):
                    raise FileNotFoundError(path)
                args += [flag, path]
        for letter, score in (wildcards or {}).items():
            if not isinstance(letter, str) or len(letter) != 1:
                raise ValueError("wildcards are single letters, not %r"
                                 % (letter,))
            args += ["--wildcard", letter, _int_option("wildcard score", score)]
        if threads is not None:
            if _int_option("threads", threads) and not 1 <= threads <= 4096:
                raise ValueError("threads must be 1 to 4096")
            args += ["--threads", str(threads)]
        for flag, on in (("--case_sensitive", case_sensitive),
                         ("--nogapsin1", no_gaps_in_a),
                         ("--nogapsin2", no_gaps_in_b),
                         ("--nogaps", no_gaps),
                         ("--nogapsinlongest", no_gaps_in_longest),
                         ("--nomismatches", no_mismatches)):
            if on:
                args.append(flag)
        self._args = args

    @property
    def command(self) -> List[str]:
        """The command line the alignments run (the pairs go to its
        standard input)."""
        return [self._exe] + self._args + ["--file", "-"]

    def _stream(self, extra, pairs, skip_empty, read_one):
        """Run the tool with options extra on pairs, yielding read_one(lines,
        a, b) for each pair in order (or skip_empty's result for a pair with
        an empty sequence, which is not sent)."""
        proc = subprocess.Popen(
            [self._exe] + self._args + extra + ["--file", "-"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, encoding="latin-1")
        # In the order written: ("pair", a, b), ("empty", a, b), an error
        # from the input, then None. Bounded, so a caller that stops reading
        # results does not have the writer pull its pairs without end
        sent = queue.Queue(maxsize=1024)
        stop = threading.Event()  # the caller left the batch
        errors = []

        def put(item):
            """Queue item, giving up (False) once the caller has left."""
            while not stop.is_set():
                try:
                    sent.put(item, timeout=0.1)
                    return True
                except queue.Full:
                    pass
            return False

        def write():
            try:
                for i, pair in enumerate(pairs):
                    if stop.is_set():
                        return
                    try:
                        a, b = pair
                        _check_sequence(a, "seq_a of pair %d" % i)
                        _check_sequence(b, "seq_b of pair %d" % i)
                    except (TypeError, ValueError) as exc:
                        put(exc)
                        return
                    if skip_empty is not None and (a == "" or b == ""):
                        if not put(("empty", a, b)):
                            return
                        continue
                    try:
                        proc.stdin.write(">a%d\n%s\n>b%d\n%s\n" % (i, a, i, b))
                    except (BrokenPipeError, OSError):
                        return  # the tool stopped: its error is reported below
                    if not put(("pair", a, b)):
                        return
            except BaseException as exc:  # the caller's iterator failed
                put(exc)
            finally:
                try:
                    proc.stdin.close()
                except OSError:
                    pass
                put(None)

        def drain_errors():
            errors.append(proc.stderr.read())

        writer = threading.Thread(target=write, daemon=True)
        reader = threading.Thread(target=drain_errors, daemon=True)
        writer.start()
        reader.start()
        lines = iter(proc.stdout.readline, "")
        try:
            while True:
                item = sent.get()
                if item is None:
                    break
                if isinstance(item, BaseException):
                    raise item
                kind, a, b = item
                if kind == "empty":
                    yield skip_empty(a, b)
                    continue
                try:
                    yield read_one(lines, a, b)
                except _Truncated:
                    break
            rest = proc.stdout.read()
            proc.wait()
            reader.join()
            message = errors[0].strip() if errors else ""
            if proc.returncode != 0 or item is not None or rest.strip():
                raise AlignmentError(
                    message or "%s stopped (exit status %s)"
                    % (self._NAME, proc.returncode))
        finally:
            stop.set()
            if proc.poll() is None:
                proc.kill()
                proc.wait()
            # The tool has gone, so its pipes are at their end: the threads
            # finish, and the streams are closed once nothing reads them
            reader.join()
            writer.join(timeout=2.0)
            for stream in (proc.stdout, proc.stderr):
                stream.close()

    @staticmethod
    def _single(results):
        """The one result of a single-pair batch, once the batch has run to
        the end (so the tool's exit status is checked too)."""
        results = list(results)
        if len(results) != 1:
            raise AlignmentError("expected one result, got %d" % len(results))
        return results[0]


class _Truncated(Exception):
    """The tool's output ended part-way through a pair."""


def _next_line(lines):
    line = next(lines, None)
    if line is None:
        raise _Truncated()
    return line.rstrip("\n")


_SCORE_LINE = re.compile(r"^score: (-?\d+)$")


def _score_from(line):
    m = _SCORE_LINE.match(line)
    if m is None:
        raise AlignmentError("unexpected output line %r" % line)
    return int(m.group(1))


class NeedlemanWunsch(_Tool):
    """Optimal global alignment (bin/needleman_wunsch).

    Scoring, as the tool's options: match, mismatch, gap_open, gap_extend
    (a gap of n letters scores gap_open + n * gap_extend), scoring (one of
    SCORING_SYSTEMS), substitution_matrix or substitution_pairs (a file, as
    in scoring/), wildcards ({letter: score}, matching every letter),
    case_sensitive. free_start_gap and free_end_gap make gaps at the ends
    free. no_gaps_in_a, no_gaps_in_b, no_gaps, no_gaps_in_longest and
    no_mismatches restrict the alignment (experimental in the tool).

    linear aligns in linear memory (Hirschberg, Myers-Miller): the same
    score, but it may choose another of several equal-scoring alignments,
    and it refuses some options. linear_fallback does so only for pairs too
    big for the score matrices; linear_above (a length) for pairs whose
    longer sequence has at least that many letters. wfa tries the wavefront
    algorithm first on similar pairs, within wfa_budget diagonal updates
    (default: the tool's, min(2^28, 256 * (len_a + len_b))), falling back
    to the matrices: again the same score, perhaps another equal-scoring
    alignment. threads: threads within an alignment and for aligning pairs
    at once (default: half the CPUs).
    """

    _NAME = "needleman_wunsch"

    def __init__(self, *, free_start_gap=False, free_end_gap=False,
                 linear=False, linear_fallback=False, linear_above=None,
                 wfa=False, wfa_budget=None, **options):
        super().__init__(**options)
        for flag, on in (("--freestartgap", free_start_gap),
                         ("--freeendgap", free_end_gap),
                         ("--linear", linear),
                         ("--linearfallback", linear_fallback),
                         ("--wfa", wfa)):
            if on:
                self._args.append(flag)
        if linear_above is not None:
            if _int_option("linear_above", linear_above) and linear_above < 1:
                raise ValueError("linear_above must be at least 1")
            self._args += ["--linearabove", str(linear_above)]
        if wfa_budget is not None:
            if _int_option("wfa_budget", wfa_budget) and wfa_budget < 1:
                raise ValueError("wfa_budget must be at least 1")
            self._args += ["--wfabudget", str(wfa_budget)]

    def align(self, seq_a: str, seq_b: str) -> GlobalAlignment:
        """The best global alignment of seq_a and seq_b. For many pairs,
        align_pairs() runs them through one tool process."""
        return self._single(self.align_pairs([(seq_a, seq_b)]))

    def align_pairs(self, pairs: Iterable[Tuple[str, str]]
                    ) -> Iterator[GlobalAlignment]:
        """Alignments of each (seq_a, seq_b) in pairs, in order, as they
        finish. pairs may be any iterable, read as the tool needs them.
        Raises AlignmentError if the tool fails on a pair, after yielding
        the pairs before it."""
        def read_one(lines, a, b):
            align_a = _next_line(lines)
            markers = _next_line(lines)
            align_b = _next_line(lines)
            score = _score_from(_next_line(lines))
            _next_line(lines)  # blank line between pairs
            return GlobalAlignment(a, b, align_a, align_b, markers, score)
        return self._stream(["--pretty", "--printscores"], pairs, None,
                            read_one)

    def score(self, seq_a: str, seq_b: str) -> int:
        """The best global alignment's score, without the alignment: one
        pass in linear memory (--scoreonly). For many pairs, score_pairs()
        runs them through one tool process."""
        return self._single(self.score_pairs([(seq_a, seq_b)]))

    def score_pairs(self, pairs: Iterable[Tuple[str, str]]) -> Iterator[int]:
        """score() for each (seq_a, seq_b) in pairs, in order."""
        def read_one(lines, a, b):
            score = _score_from(_next_line(lines))
            _next_line(lines)
            return score
        return self._stream(["--scoreonly"], pairs, None, read_one)


_HIT_LINE = re.compile(r"^hit \d+\.\d+ score: (-?\d+)$")
_PART_LINE = re.compile(r"^  (.*)  \[pos: (\d+); len: (\d+)\]$")


class SmithWaterman(_Tool):
    """Optimal local alignments (bin/smith_waterman), best first.

    Scoring options as for NeedlemanWunsch (no free end gaps or linear
    memory). min_score: report hits scoring at least this (the tool's
    default is match * max(0.2 * shorter length, 2)); max_hits: at most this
    many hits per pair (default: all, which may be many for long sequences).
    fast_local (with max_hits=1) reconstructs the best hit of a long pair in
    little memory: the same score, perhaps another equal-scoring hit. A pair
    with an empty sequence has no hits.
    """

    _NAME = "smith_waterman"

    def __init__(self, *, min_score=None, max_hits=None, fast_local=False,
                 **options):
        super().__init__(**options)
        if min_score is not None:
            self._args += ["--minscore", _int_option("min_score", min_score)]
        if max_hits is not None:
            if _int_option("max_hits", max_hits) and max_hits < 1:
                raise ValueError("max_hits must be at least 1")
            self._args += ["--maxhits", str(max_hits)]
        if fast_local:
            self._args.append("--fastlocal")

    def align(self, seq_a: str, seq_b: str) -> List[LocalHit]:
        """The local alignments of seq_a and seq_b, best first. For many
        pairs, align_pairs() runs them through one tool process."""
        return self._single(self.align_pairs([(seq_a, seq_b)]))

    def align_pairs(self, pairs: Iterable[Tuple[str, str]]
                    ) -> Iterator[List[LocalHit]]:
        """The hits of each (seq_a, seq_b) in pairs, in order. Raises
        AlignmentError if the tool fails on a pair, after yielding the pairs
        before it."""
        def read_one(lines, a, b):
            line = _next_line(lines)
            if not line.startswith("== Alignment "):
                raise AlignmentError("unexpected output line %r" % line)
            _next_line(lines)
            hits = []
            while True:
                line = _next_line(lines)
                if line == "==":
                    return hits
                m = _HIT_LINE.match(line)
                part_a = _PART_LINE.match(_next_line(lines))
                markers = _next_line(lines)[2:]
                part_b = _PART_LINE.match(_next_line(lines))
                _next_line(lines)
                if m is None or part_a is None or part_b is None:
                    raise AlignmentError("unexpected hit output near %r"
                                         % line)
                hits.append(LocalHit(
                    part_a.group(1), part_b.group(1), markers,
                    int(m.group(1)), int(part_a.group(2)),
                    int(part_a.group(3)), int(part_b.group(2)),
                    int(part_b.group(3))))
        return self._stream(["--pretty"], pairs, lambda a, b: [], read_one)
