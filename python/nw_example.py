#!/usr/bin/env python3
"""Global alignments with nucleics_align (as perl/nw_example.pl).

python/nw_example.py
url: https://github.com/DanielTillett/nucleics-align
license: BSD 3-Clause, see LICENSE

Usage: ./nw_example.py [seq1 seq2]
"""

import sys

from nucleics_align import NeedlemanWunsch


def main(argv):
    if len(argv) not in (1, 3):
        sys.exit("Usage: ./nw_example.py [seq1 seq2]")
    pairs = [tuple(argv[1:3])] if len(argv) == 3 else [
        ("asdfasdf", "asdfasdf"), ("dogg", "ggod"), ("asdfas", "fadsas"),
        ("", ""), ("asdf", "lkj")]

    nw = NeedlemanWunsch()
    # All the pairs go to one run of the tool, aligned on its threads
    for aln in nw.align_pairs(pairs):
        print("Doing '%s' vs '%s'" % (aln.seq_a, aln.seq_b))
        print(aln.align_a, aln.markers, aln.align_b, sep="\n")
        print("score: %d\n" % aln.score)

    # The score alone, in linear memory; and a linear-memory alignment
    a, b = pairs[-1]
    print("score only: %d" % nw.score(a, b))
    print("linear memory: %d" % NeedlemanWunsch(linear=True).align(a, b).score)


if __name__ == "__main__":
    main(sys.argv)
