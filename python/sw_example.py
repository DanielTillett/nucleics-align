#!/usr/bin/env python3
"""Local alignments with nucleics_align (as perl/sw_example.pl).

python/sw_example.py
url: https://github.com/DanielTillett/nucleics-align
license: BSD 3-Clause, see LICENSE

Usage: ./sw_example.py [seq1 seq2]
"""

import sys

from nucleics_align import SmithWaterman


def main(argv):
    if len(argv) not in (1, 3):
        sys.exit("Usage: ./sw_example.py [seq1 seq2]")
    pairs = [tuple(argv[1:3])] if len(argv) == 3 else [
        ("asdfasdf", "asdfasdf"), ("dogg", "ggod"), ("asdfas", "fadsas"),
        ("hello_henry", "wellohenry"), ("asdf", "lkj")]

    sw = SmithWaterman(min_score=2)
    for (a, b), hits in zip(pairs, sw.align_pairs(pairs)):
        print("Doing '%s' vs '%s'" % (a, b))
        if hits:
            best = hits[0]
            print("hit:", best.align_a, best.markers, best.align_b,
                  sep="\n  ")
        print("  (%d hits)" % len(hits))


if __name__ == "__main__":
    main(sys.argv)
