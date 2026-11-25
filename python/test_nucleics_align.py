"""Tests of python/nucleics_align.py against the built tools.

python/test_nucleics_align.py
url: https://github.com/DanielTillett/nucleics-align
license: BSD 3-Clause, see LICENSE

Run from the repository root after make: python3 -m unittest discover python
"""

import os
import random
import subprocess
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nucleics_align as na  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
INPUT = os.path.join(ROOT, "tests", "cli", "input")


def tool_output(name, *args, stdin=None):
    """What the tool prints itself, to compare with."""
    exe = os.path.join(ROOT, "bin", name)
    return subprocess.run([exe] + list(args), input=stdin, text=True,
                          capture_output=True).stdout


def random_pairs(n, seed=1, alphabet="ACGT"):
    rng = random.Random(seed)
    pairs = []
    for _ in range(n):
        a = "".join(rng.choice(alphabet) for _ in range(rng.randrange(1, 300)))
        b = "".join(c if rng.random() > 0.1 else rng.choice(alphabet)
                    for c in a if rng.random() > 0.03)
        pairs.append((a, b or "A"))
    return pairs


class NeedlemanWunschTest(unittest.TestCase):
    def test_one_pair(self):
        aln = na.NeedlemanWunsch().align("ACGT", "AGT")
        self.assertEqual((aln.align_a, aln.markers, aln.align_b, aln.score),
                         ("ACGT", "| ||", "A-GT", -2))

    def test_same_as_tool(self):
        # Every pair as the tool prints it with the same options
        pairs = random_pairs(60)
        fasta = "".join(">a\n%s\n>b\n%s\n" % p for p in pairs)
        for options, kwargs in (
                ([], {}),
                (["--linear"], {"linear": True}),
                (["--freestartgap", "--freeendgap"],
                 {"free_start_gap": True, "free_end_gap": True}),
                (["--match", "2", "--mismatch", "-3", "--gapopen", "-5",
                  "--gapextend", "-2"],
                 {"match": 2, "mismatch": -3, "gap_open": -5,
                  "gap_extend": -2})):
            want = tool_output("needleman_wunsch", "--pretty", "--printscores",
                               *options, "--file", "-", stdin=fasta)
            got = "".join("%s\n%s\n%s\nscore: %d\n\n" % (
                a.align_a, a.markers, a.align_b, a.score)
                for a in na.NeedlemanWunsch(**kwargs).align_pairs(pairs))
            self.assertEqual(got, want, options)

    def test_threads_keep_order(self):
        pairs = random_pairs(200, seed=2)
        one = list(na.NeedlemanWunsch(threads=1).align_pairs(pairs))
        many = list(na.NeedlemanWunsch(threads=4).align_pairs(iter(pairs)))
        self.assertEqual(one, many)
        self.assertEqual([(a.seq_a, a.seq_b) for a in many], pairs)

    def test_scores(self):
        pairs = random_pairs(40, seed=3)
        nw = na.NeedlemanWunsch()
        self.assertEqual(list(nw.score_pairs(pairs)),
                         [a.score for a in nw.align_pairs(pairs)])
        self.assertEqual(nw.score("ACGT", "AGT"), -2)

    def test_empty_sequences(self):
        nw = na.NeedlemanWunsch()
        got = list(nw.align_pairs([("", ""), ("ACGT", ""), ("", "AG"),
                                   ("ACGT", "AGT")]))
        self.assertEqual([a.score for a in got], [0, -8, -6, -2])
        self.assertEqual((got[1].align_a, got[1].align_b), ("ACGT", "----"))
        free_end = na.NeedlemanWunsch(free_end_gap=True)
        self.assertEqual([a.score for a in free_end.align_pairs(
            [("ACGT", ""), ("", "ACGT")])], [0, 0])

    def test_protein_scoring(self):
        a, b = "MKTAYIAKQRQISFVKSHFSRQ", "MKTAHIAKQRQISFVKSHFSRQ"
        builtin = na.NeedlemanWunsch(scoring="blosum62").align(a, b)
        from_file = na.NeedlemanWunsch(substitution_matrix=os.path.join(
            ROOT, "scoring", "BLOSUM62.txt")).align(a, b)
        self.assertEqual(builtin, from_file)

    def test_wildcards_and_case(self):
        self.assertEqual(na.NeedlemanWunsch(wildcards={"N": 1}).align(
            "ACNT", "ACGT").score, 4)
        self.assertEqual(na.NeedlemanWunsch(case_sensitive=True).align(
            "acgt", "ACGT").score,
            na.NeedlemanWunsch().align("acgt", "xxxx").score)

    def test_tool_error(self):
        # The pair before the failing one still comes out
        nw = na.NeedlemanWunsch(substitution_pairs=os.path.join(
            INPUT, "ac_pairs.txt"))
        it = nw.align_pairs([("ACCA", "ACA"), ("ACGT", "AGT")])
        self.assertEqual(next(it).score, -2)
        with self.assertRaisesRegex(na.AlignmentError, "Unknown character"):
            next(it)

    def test_linear_refused(self):
        nw = na.NeedlemanWunsch(linear=True, gap_open=-300000,
                                gap_extend=-300000)
        it = nw.align_pairs([("acgtacgt", "acgaacgt"), ("a" * 3000, "c" * 3000)])
        self.assertEqual(next(it).score, 5)
        with self.assertRaisesRegex(na.AlignmentError, "--linear"):
            next(it)

    def test_bad_input(self):
        nw = na.NeedlemanWunsch()
        with self.assertRaises(ValueError):
            nw.align("AC GT", "ACGT")
        with self.assertRaises(ValueError):
            nw.align(">ACGT", "ACGT")
        with self.assertRaises(ValueError):
            list(nw.align_pairs([("ACGT", "AGT"), ("A\nC", "AC")]))
        with self.assertRaises(ValueError):
            na.NeedlemanWunsch(scoring="NOPE")
        with self.assertRaises(TypeError):
            na.NeedlemanWunsch(match=1.5)
        with self.assertRaises(FileNotFoundError):
            na.NeedlemanWunsch(substitution_matrix="no/such/file")

    def test_stop_early(self):
        # Leaving a batch part-way stops the tool
        it = na.NeedlemanWunsch().align_pairs(random_pairs(500, seed=4))
        next(it)
        it.close()

    def test_single_pair_checks_exit_status(self):
        # align() finishes the tool's run: an exit status after the output
        # is still an error (a stand-in tool prints an alignment, then fails)
        import stat
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            fake = os.path.join(d, "needleman_wunsch")
            with open(fake, "w") as f:
                f.write("#!/bin/sh\nprintf 'ACGT\\n| ||\\nA-GT\\nscore: -2\\n\\n'\n"
                        "echo 'Error: late failure' >&2\nexit 3\n")
            os.chmod(fake, os.stat(fake).st_mode | stat.S_IXUSR)
            with self.assertRaisesRegex(na.AlignmentError, "late failure"):
                na.NeedlemanWunsch(bin_dir=d).align("ACGT", "AGT")

    def test_input_iterator_error(self):
        # An error from the caller's own iterator is the caller's error,
        # not a complete batch
        def pairs():
            yield ("ACGT", "AGT")
            raise OSError("input source failed")
        with self.assertRaisesRegex(OSError, "input source failed"):
            list(na.NeedlemanWunsch().align_pairs(pairs()))


class SmithWatermanTest(unittest.TestCase):
    def test_hits(self):
        hits = na.SmithWaterman(max_hits=2).align("acgtacgtacgt",
                                                 "ttacgtacgtaa")
        self.assertEqual(len(hits), 2)
        h = hits[0]
        self.assertEqual((h.score, h.align_a, h.markers, h.align_b, h.pos_a,
                          h.len_a, h.pos_b, h.len_b),
                         (18, "acgtacgta", "|||||||||", "acgtacgta", 0, 9,
                          2, 9))

    def test_same_as_tool(self):
        pairs = random_pairs(40, seed=5)
        fasta = "".join(">a\n%s\n>b\n%s\n" % p for p in pairs)
        want = tool_output("smith_waterman", "--pretty", "--maxhits", "3",
                           "--minscore", "5", "--file", "-", stdin=fasta)
        got = []
        for i, hits in enumerate(na.SmithWaterman(
                max_hits=3, min_score=5).align_pairs(pairs)):
            a, b = pairs[i]
            got.append("== Alignment %d lengths (%d, %d):\n\n"
                       % (i, len(a), len(b)))
            for j, h in enumerate(hits):
                got.append("hit %d.%d score: %d\n  %s  [pos: %d; len: %d]\n"
                           "  %s\n  %s  [pos: %d; len: %d]\n\n"
                           % (i, j, h.score, h.align_a, h.pos_a, h.len_a,
                              h.markers, h.align_b, h.pos_b, h.len_b))
            got.append("==\n")
        self.assertEqual("".join(got), want)

    def test_empty_and_no_hits(self):
        got = list(na.SmithWaterman().align_pairs(
            [("", "ACGT"), ("AAAA", "CCCC"), ("acgtacgt", "acgtacgt")]))
        self.assertEqual(got[0], [])
        self.assertEqual(got[1], [])
        self.assertEqual(got[2][0].score, 16)

    def test_stop_early_stops_reading_input(self):
        # Pairs with an empty sequence never reach the tool; closing the
        # batch must still stop the reader of the caller's pairs
        import time
        pulled = [0]

        def endless():
            while True:
                pulled[0] += 1
                yield ("", "ACGT")
        it = na.SmithWaterman().align_pairs(endless())
        self.assertEqual(next(it), [])
        it.close()
        time.sleep(0.3)
        seen = pulled[0]
        time.sleep(0.3)
        self.assertEqual(pulled[0], seen)
        self.assertLess(seen, 100000)


if __name__ == "__main__":
    unittest.main()
