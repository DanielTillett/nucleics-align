seq-align
=========
Smith-Waterman and Needleman-Wunsch alignment in C
url: https://github.com/DanielTillett/seq-align
original author: Isaac Turner <turner.isaac@gmail.com>
license: BSD 3-Clause (see [LICENSE](LICENSE))
updated: 18 Nov 2026

[![CI](https://github.com/DanielTillett/seq-align/actions/workflows/ci.yml/badge.svg)](https://github.com/DanielTillett/seq-align/actions/workflows/ci.yml)

About
=====

C implementations of optimal local (Smith-Waterman) and global
(Needleman-Wunsch) alignment with affine gap penalties. The command line
tools `smith_waterman`, `needleman_wunsch` and `lcs` cover most uses, and the
library is easy to call from other programs (see `examples/`). Perl wrapper
modules are in `perl/`.

This is a fork of [noporpoise/seq-align](https://github.com/noporpoise/seq-align)
with bug fixes, a linear-space global aligner and a faster SIMD matrix fill.

Features:
* Align any pair of ASCII sequences (DNA, protein, words, etc.)
* Choose a scoring system (BLOSUM, PAM etc.) or load your own
* Wildcards that match any character with a given score
* Align with no mismatches (`--nomismatches`) or no gaps (`--nogaps`,
  `--nogapsin1`, `--nogapsin2`, `--nogapsinlongest`) in local and global
  alignment. Used together for local alignment, this lists common substrings
  in order of length/score
* Penalty-free gaps at the start or end of a global alignment
  (`--freestartgap` and `--freeendgap`)
* Global alignment in linear memory (`needleman_wunsch --linear`) for long
  sequences
* Striped SIMD matrix fill (Farrar's method) with output identical to the
  scalar code
* Read FASTA, FASTQ, SAM, BAM, plain (one sequence per line) and gzipped files
* Display in colour (`--colour`) and show context around local alignments
  (`--context <n>`)
* Case-insensitive matching by default; change with `--case_sensitive`

Build
-----

    $ git clone --recursive https://github.com/DanielTillett/seq-align
    $ cd seq-align
    $ make

Build options:

    $ make test                           # run the tests
    $ make clean && make SANITIZE=1 test  # tests with ASan and UBSan
    $ make NATIVE=1                       # tune for this CPU (e.g. AVX2)
    $ make DEBUG=1                        # no optimisation, debug symbols

`NATIVE=1` adds `-march=native`, which lets the SIMD fill use 8-lane AVX2 on
x86 (4-lane SSE2 or NEON otherwise). Binaries built this way may not run on
other CPUs.

To update:

    $ git pull
    $ git submodule update --init

Bundled dependencies (in `libs/`, under their own licences):

* string_buffer [https://github.com/noporpoise/string_buffer]
* seq_file [https://github.com/noporpoise/seq_file]
* sort_r [https://github.com/noporpoise/sort_r]

Using the library
-----------------

Link against `src/libalign.a` (and `libs/string_buffer/libstrbuf.a`, `-lz`).
The main entry points are:

* `needleman_wunsch_align()` (`src/needleman_wunsch.h`): global alignment
* `smith_waterman_align()` then `smith_waterman_fetch()` repeatedly
  (`src/smith_waterman.h`): local alignments, best first
* `hirschberg_align()` (`src/hirschberg.h`): global alignment in O(n) memory;
  same scores as `needleman_wunsch_align()`, but returns false for options it
  does not support (free end gaps, `no_gaps_in_*`, `no_mismatches`)
* `scoring_init()` and the `scoring_system_*()` functions
  (`src/alignment_scoring.h`): set up a scoring system

`alignment_striped_enabled` (`src/alignment.h`, default true) switches the
SIMD fill off, for testing or comparison. It gives exactly the same matrices
as the scalar fill.

See `examples/` for complete C and C++ programs.

Performance
-----------

Time per alignment of simulated sequences (about 10% substitutions and 2%
insertions/deletions), measured on an AMD Ryzen AI Max+ 395 under WSL2 with
gcc 15.2 at `-O3`. "Before" is the code at the start of this fork; every
version gives identical scores and alignments.

| Workload | Before | Now | Now, `NATIVE=1` |
|---|---|---|---|
| Needleman-Wunsch, 5 kb DNA | 0.252 s | 0.097 s | 0.063 s |
| Needleman-Wunsch, 2 kb protein (BLOSUM62) | 0.031 s | 0.011 s | 0.007 s |
| Linear space (`--linear`), 5 kb DNA | n/a | 0.077 s | 0.077 s |
| Smith-Waterman, 3 kb DNA, best 20 hits | 0.636 s | 0.157 s | 0.153 s |

Time and peak memory for one global alignment (`needleman_wunsch
--printscores`, whole program):

| Sequences | Full matrices | `--linear` |
|---|---|---|
| 5 kb | 0.24 s, 289 MB | 0.07 s, 2 MB |
| 20 kb | 4.8 s, 4.6 GB | 1.2 s, 2 MB |

Where the speed comes from:
* A substitution score table for the characters in the two sequences replaces
  a `scoring_lookup()` call per cell
* The matrix fill uses Farrar's striped SIMD layout (Farrar 2007,
  Bioinformatics 23:156-161), written with GCC/clang vector extensions so it
  builds on x86 and ARM. Options that change the recurrence for particular
  cells (free end gaps, `--nogaps..`, `--nomismatches`) use the scalar fill
* Smith-Waterman hits are ordered with a counting sort on score instead of a
  comparison sort of every cell. Tied hits are ordered by position in the
  first sequence, then the second, on every platform

Most of the remaining Smith-Waterman time is spent following hits back through
the matrices, which defines which hits are reported, so it is unchanged.

Examples
========

Basics:

    $ ./bin/needleman_wunsch CAGACGT CGATA
    C-AGACGT
    CGATA---

Print alignment scores:

    $ ./bin/needleman_wunsch --printscores CAGACGT CGATA
    C-AGACGT
    CGATA---
    score: -11

Read from file (dna.fa.gz):

    >seqA
    ACAATAGAC
    >seqB
    ACGAATAGAT
    >seqC
    ACGTGA
    CAGAT
    >seqD
    GTGGACG
    AGTA

    $ ./bin/needleman_wunsch --printscores --file dna.fa.gz
    AC-AATAGAC
    ACGAATAGAT
    score: -3

    ACGTGACAGAT
    GTGGACGAGTA
    score: -5


Reading from STDIN:

    $ gzip -d -c dna.fa.gz | ./needleman_wunsch --stdin
    AC-AATAGAC
    ACGAATAGAT

    ACGTGACAGAT
    GTGGACGAGTA

is the same as:

    $ gzip -d -c dna.fa.gz | ./needleman_wunsch --file -

Set different scoring systems:

    $ ./bin/needleman_wunsch --match 1 --mismatch 0 --gapopen -10 --gapextend 0 ACGTGCCCCACAGAT AGGTGGACGAGAT

Print the scoring matrices:

    $ ./bin/needleman_wunsch --printmatrices ACAGGT AAGGT
    seq_a: ACAGGT
    seq_b: AAGGT
    match_scores:
      0:    0 -2147483643 -2147483643 -2147483643 -2147483643 -2147483643 -2147483643
      1:  -2147483643   1  -7  -5  -9 -10 -11
      2:  -2147483643  -4  -1  -3  -7  -8  -9
      3:  -2147483643  -8  -6  -3  -2  -6 -10
      4:  -2147483643  -9  -7  -8  -2  -1  -8
      5:  -2147483643 -10  -8  -9 -10  -4   0
    gap_a_scores:
      0:    0 -2147483643 -2147483643 -2147483643 -2147483643 -2147483643 -2147483643
      1:   -5 -10 -11 -12 -13 -14 -15
      2:   -6  -4  -9 -10 -11 -12 -13
      3:   -7  -5  -6  -8 -12 -13 -14
      4:   -8  -6  -7  -8  -7 -11 -13
      5:   -9  -7  -8  -9  -7  -6 -11
    gap_b_scores:
      0:    0  -5  -6  -7  -8  -9 -10
      1:  -2147483643 -10  -4  -5  -6  -7  -8
      2:  -2147483643 -11  -9  -6  -7  -8  -9
      3:  -2147483643 -12 -10 -11  -8  -7  -8
      4:  -2147483643 -13 -11 -12 -13  -7  -6
      5:  -2147483643 -14 -12 -13 -14 -12  -9
    match: 1 mismatch: -2 gapopen: -4 gapexend: -1

    ACAGGT
    A-AGGT


Smith-Waterman
==============

    usage: ./bin/smith_waterman [OPTIONS] [seq1 seq2]
      Smith-Waterman optimal local alignment (maximises score).
      Takes a pair of sequences on the command line, or can read from a
      file and from sequence piped in.  Can read gzip files, FASTA and FASTQ.

      OPTIONS:
        --file <file>        Sequence file reading with gzip support - read two
                             sequences at a time and align them
        --files <f1> <f2>    Read one sequence from each file to align at one time
        --stdin              Read from STDIN (same as '--file -')

        --case_sensitive     Use case sensitive character comparison [default: off]

        --match <score>      [default: 2]
        --mismatch <score>   [default: -2]
        --gapopen <score>    [default: -2]
        --gapextend <score>  [default: -1]

        --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62|DNA_HYBRIDIZATION>
        --substitution_matrix <file>  see details for formatting
        --substitution_pairs <file>   see details for formatting

        --wildcard <w> <s>   Character <w> matches all characters with score <s>

        --minscore <score>   Minimum required score
                             [default: match * MAX(0.2 * length, 2)]
        --maxhits <hits>     Maximum number of results per alignment
                             [default: no limit]

        --context <n>        Print <n> bases of context
        --printseq           Print sequences before local alignments
        --printmatrices      Print dynamic programming matrices
        --printfasta         Print fasta header lines
        --pretty             Print with a descriptor line
        --colour             Print with colour

      Experimental Options:
        --nogapsin1          No gaps allowed within the first sequence
        --nogapsin2          No gaps allowed within the second sequence
        --nogaps             No gaps allowed in either sequence
        --nogapsinlongest    No gaps allowed within the longer sequence
        --nomismatches       No mismatches allowed

     DETAILS:
      * For help choosing scoring, see the README file.
      * Gap (of length N) penalty is: (open+N*extend)
      * To do alignment without affine gap penalty, set '--gapopen 0'.
      * Scoring files should be matrices, with entries separated by a single
        character or whitespace. See files in the 'scoring' directory for examples.

      https://github.com/DanielTillett/seq-align  (compiled: <date>)


Needleman-Wunsch
================

    usage: ./bin/needleman_wunsch [OPTIONS] [seq1 seq2]
      Needleman-Wunsch optimal global alignment (maximises score).
      Takes a pair of sequences on the command line, or can read from a
      file and from sequence piped in.  Can read gzip files, FASTA and FASTQ.

      OPTIONS:
        --file <file>        Sequence file reading with gzip support - read two
                             sequences at a time and align them
        --files <f1> <f2>    Read one sequence from each file to align at one time
        --stdin              Read from STDIN (same as '--file -')

        --case_sensitive     Use case sensitive character comparison [default: off]

        --match <score>      [default: 1]
        --mismatch <score>   [default: -2]
        --gapopen <score>    [default: -4]
        --gapextend <score>  [default: -1]

        --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62|DNA_HYBRIDIZATION>
        --substitution_matrix <file>  see details for formatting
        --substitution_pairs <file>   see details for formatting

        --wildcard <w> <s>   Character <w> matches all characters with score <s>


        --freestartgap       No penalty for gap at start of alignment
        --freeendgap         No penalty for gap at end of alignment

        --printscores        Print optimal alignment scores
        --zam                A funky type of output
        --linear             Align in linear memory (Hirschberg, Myers-Miller)
        --printmatrices      Print dynamic programming matrices
        --printfasta         Print fasta header lines
        --pretty             Print with a descriptor line
        --colour             Print with colour

      Experimental Options:
        --nogapsin1          No gaps allowed within the first sequence
        --nogapsin2          No gaps allowed within the second sequence
        --nogaps             No gaps allowed in either sequence
        --nogapsinlongest    No gaps allowed within the longer sequence
        --nomismatches       No mismatches allowed (cannot be used with --nogaps..)

     DETAILS:
      * For help choosing scoring, see the README file.
      * Gap (of length N) penalty is: (open+N*extend)
      * To do alignment without affine gap penalty, set '--gapopen 0'.
      * Scoring files should be matrices, with entries separated by a single
        character or whitespace. See files in the 'scoring' directory for examples.

      https://github.com/DanielTillett/seq-align  (compiled: <date>)

`--linear` cannot be combined with `--freestartgap`, `--freeendgap`,
`--nogaps..`, `--nomismatches` or `--printmatrices`. When several alignments
share the best score, it may return a different one from the default aligner.


Longest Common Substring
========================

    ./bin/lcs [options] <sequence>
      Print substrings in decreasing order of length


Scoring Penalties
-----------------

Proteins:

| Query Length  | Substitution Matrix | Gap Costs (open,extend) |
|---------------|---------------------|-------------------------|
| <35           | PAM-30              | (9,1)                   |
| 35-50         | PAM-70              | (10,1)                  |
| 50-85         | BLOSUM-80           | (10,1)                  |
| 85            | BLOSUM-62           | (10,1)                  |

[Table from: http://www.ncbi.nlm.nih.gov/blast/html/sub_matrix.html]

    gap (of length N) penalty: gap_open + N*gap_extend

NCBI BLAST Quote [from: http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml#Reward-penalty]:

    Many nucleotide searches use a simple scoring system that consists of a "reward"
    for a match and a "penalty" for a mismatch. The (absolute) reward/penalty ratio
    should be increased as one looks at more divergent sequences. A ratio of 0.33
    (1/-3) is appropriate for sequences that are about 99% conserved; a ratio of 0.5
    (1/-2) is best for sequences that are 95% conserved; a ratio of about one (1/-1)
    is best for sequences that are 75% conserved [1].

NCBI Gap (open, extend) values:

| open | extend |
|------|--------|
| -5   | -2     |
| -2   | -2     |
| -1   | -2     |
| -0   | -2     |
| -3   | -1     |
| -2   | -1     |
| -1   | -1     |

Our default (for now) are:  
gap_open/gap_extend: (-4,-1)  
match/mismatch: (1,-2)

 for help on choosing scoring matrices see:
 http://www.ebi.ac.uk/help/matrix.html

License
=======

seq-align is licensed under the BSD 3-Clause License; see [LICENSE](LICENSE).

It was originally written by Isaac Turner and released into the public domain
(Unlicense). The bundled libraries in `libs/` keep their own licences. The
substitution matrices in `scoring/` come from NCBI
(ftp://ftp.ncbi.nih.gov/blast/matrices/).


Development
===========

Bug reports and feature requests:
https://github.com/DanielTillett/seq-align/issues

Done:
* `--nogapsinlongest`: no gaps within the longer of the two sequences
  (the first sequence if they are the same length)
* Gotoh: affine gap penalties use Gotoh's three-matrix recurrence
  (match, gap in a, gap in b) for both Smith-Waterman and Needleman-Wunsch
  (Gotoh 1982, J Mol Biol 162:705-708)
* Hirschberg's linear-space algorithm with the affine-gap extension of
  Myers & Miller (1988, Comput Appl Biosci 4:11-17): `hirschberg_align()`
  and `needleman_wunsch --linear`
* Striped SIMD matrix fill (Farrar 2007, Bioinformatics 23:156-161)

Algorithms to investigate:
* SIMD for the linear-space aligner, which still uses a scalar loop
* Myers' bit-vector algorithm (1999, J ACM 46:395-415). Only computes unit
  cost edit distance, so it would be a separate fast path rather than a
  replacement for the scored aligners
* Combining Myers' and Hirschberg's algorithms for linear-space
  edit-distance alignment
