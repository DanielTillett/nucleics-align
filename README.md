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

`NATIVE=1` adds `-march=native`, which lets the SIMD code use 8-lane AVX2 on
x86 (4 lanes otherwise). Binaries built this way may not run on other CPUs.

Without it, x86-64 Linux builds with gcc or clang compile the SIMD code for
AVX2, SSE4.1 and plain SSE2, and pick one for the CPU when the program
starts. That gets most of the `NATIVE=1` speed in a portable binary (see
Performance). Plain SSE2 has no vector max instruction for 32-bit integers,
so on such CPUs the linear-space aligner uses scalar loops. Other platforms
use the compiler's target: ARM NEON (e.g. Apple silicon) gets all the SIMD
code, while x86 builds elsewhere (e.g. macOS on Intel) can add
`OPT="-O3 -msse4.1"` or `NATIVE=1`. ThreadSanitizer builds skip the run-time
selection, which it cannot handle.

To update:

    $ git pull
    $ git submodule update --init

Bundled dependencies (in `libs/`, under their own licences):

* string_buffer [https://github.com/noporpoise/string_buffer]
* seq_file [https://github.com/noporpoise/seq_file]
* sort_r [https://github.com/noporpoise/sort_r]

Using the library
-----------------

Link against `src/libalign.a` (and `libs/string_buffer/libstrbuf.a`, `-lz`,
`-lpthread`). The main entry points are:

* `needleman_wunsch_align()` (`src/needleman_wunsch.h`): global alignment
* `smith_waterman_align()` then `smith_waterman_fetch()` repeatedly
  (`src/smith_waterman.h`): local alignments, best first.
  `smith_waterman_set_min_score()` skips hits below a score: the hits
  returned are the same, but lower-scoring cells are never indexed
* `hirschberg_align()` (`src/hirschberg.h`): global alignment in O(n) memory;
  same scores as `needleman_wunsch_align()`, but returns false for options it
  does not support (a positive gap open, free end gaps, `no_gaps_in_*`,
  `no_mismatches`). Large sub-problems run their forward and reverse passes
  at the same time, each split into column blocks on several threads
  (tunable with `hirschberg_block_min_cells`, `hirschberg_block_min_cols`
  and `hirschberg_max_threads`). Passes skip cells that cannot be on a best
  path (`hirschberg_band_enabled`, default true), with the same result
* `hirschberg_score()` (`src/hirschberg.h`): only the global alignment score,
  in O(n) memory; same options as `hirschberg_align()`
* `scoring_init()` and the `scoring_system_*()` functions
  (`src/alignment_scoring.h`): set up a scoring system

`alignment_striped_enabled` (`src/alignment.h`, default true) switches SIMD
off in both the matrix fill and the linear-space aligner, for testing or
comparison; results are exactly the same either way. Large full-matrix fills
are split into column blocks on several threads, tunable with
`alignment_block_min_cells`, `alignment_block_min_cols` and
`alignment_max_threads`.

See `examples/` for complete C and C++ programs.

Performance
-----------

Time per alignment of simulated sequences (about 10% substitutions and 2%
insertions/deletions), measured on an AMD Ryzen AI Max+ 395 under WSL2 with
gcc 15.2 at `-O3`. "Before" is the code at the start of this fork; every
version gives identical scores and alignments.

| Workload | Before | Now | Now, `NATIVE=1` |
|---|---|---|---|
| Needleman-Wunsch, 5 kb DNA | 0.252 s | 0.046 s | 0.045 s |
| Needleman-Wunsch, 2 kb protein (BLOSUM62) | 0.031 s | 0.005 s | 0.005 s |
| Linear space (`--linear`), 5 kb DNA | n/a | 0.030 s | 0.026 s |
| Smith-Waterman, 3 kb DNA, best 20 hits | 0.636 s | 0.096 s | 0.094 s |

Timings vary by a few per cent between runs.

Time and peak memory for one global alignment (`needleman_wunsch
--printscores`, whole program):

| Sequences | Full matrices | `--linear` | `--scoreonly` |
|---|---|---|---|
| 5 kb | 0.03-0.15 s, 289 MB | 0.01 s, 2 MB | <0.01 s, 2 MB |
| 20 kb | 0.35-2.5 s, 4.6 GB | 0.10 s, 3 MB | 0.02 s, 3 MB |
| 50 kb | not run (needs ~30 GB) | 0.34 s, 4 MB | 0.09 s, 4 MB |

`--linear` and `--scoreonly` are fastest for similar sequences, where most
of the matrix can be skipped (see banding below). For 20 kb and 50 kb pairs
with substitutions and equal numbers of insertions and deletions:

| Divergence | 20 kb `--linear` | 20 kb `--scoreonly` | 50 kb `--linear` | 50 kb `--scoreonly` |
|---|---|---|---|---|
| 1% substitutions, 0.2% indels | 0.05 s | 0.01 s | 0.18 s | 0.05 s |
| 5%, 1% | 0.09 s | 0.03 s | 0.32 s | 0.06 s |
| 10%, 2% | 0.11 s | 0.02 s | 0.43 s | 0.12 s |
| 25%, 5% | 0.16 s | 0.04 s | 0.74 s | 0.20 s |
| unrelated | 0.13 s | 0.02 s | 0.60 s | 0.13 s |

Before banding the same runs took 0.11-0.22 s, 0.02-0.05 s, 0.53-1.06 s and
0.11-0.27 s, fastest for the unrelated pair.

Full-matrix times vary most: the first run after others have freed memory is
the slowest, while the kernel finds and clears huge pages.

Measured on 32 CPUs; large fills and passes use up to half of them. Threads
wait for each other by spinning, so CPU time can be several times the wall
time (50 kb unrelated sequences, `--linear`: 0.6 s wall, 6-8 s CPU).

The command line tools align one pair at a time. Splitting a batch across
processes scales well: 32 pairs ran 4x faster as 8 processes for
Needleman-Wunsch (5 kb) and 6x faster as 16 for Smith-Waterman (3 kb).

Where the speed comes from:
* A substitution score table for the characters in the two sequences replaces
  a `scoring_lookup()` call per cell
* The matrix fill uses Farrar's striped SIMD layout (Farrar 2007,
  Bioinformatics 23:156-161), written with GCC/clang vector extensions so it
  builds on x86 and ARM. Options that change the recurrence for particular
  cells (free end gaps, `--nogaps..`, `--nomismatches`) use the scalar fill
* The linear-space aligner's forward and reverse passes use the same striped
  layout, with Farrar's lazy F loop carrying insertion gaps across segment
  boundaries. Its rows match the scalar passes exactly, so the alignment is
  unchanged. Used where the CPU has a vector max instruction (see Build), on
  sub-problems large enough to repay building the score profile, with gap
  penalties of zero or less, and while every score stays far from 32-bit
  overflow
* In large passes the columns are split into blocks, one thread each, every
  thread one row behind the block on its left. Each block hands the next
  its last column's score and the insertion gap leaving it, so the result is
  the same as one thread. The forward and reverse passes also run at the
  same time (20 kb `--linear` 0.79 -> 0.19 s, `--scoreonly` 0.30 -> 0.04 s)
* The full-matrix striped fill is split into column blocks the same way,
  each block reading its left neighbour's last column straight from the
  matrices (20 kb Needleman-Wunsch 0.78 -> 0.35 s). Smith-Waterman gains less
  (20 kb 6.6 -> 6.2 s), as most of its time goes on finding hits
* On Linux the score matrices are allocated for transparent huge pages:
  page faults for a 5 kb x 5 kb alignment fall from ~72,000 to ~320, and
  fewer TLB misses speed up the fill itself (5 kb NW 0.060 -> 0.046 s)
* Exact banding in the linear-space aligner: a path through a cell far from
  the main diagonals needs many gap characters, which caps its score. Cells
  that cannot reach the best score are skipped, by running each pass in
  strips of rows over only the columns that might. The bound comes from the
  sub-problem's best score, which the split above it already found; the top
  level first tries a narrow band. If the best score found needs a wider band
  than was used, the pass is rerun wider (or whole), so the result is always
  the same. For similar sequences most cells are skipped: 50 kb at 1%
  divergence 0.9 -> 0.18 s and ~10 s -> 0.25 s of CPU; for unrelated ones
  the narrow first try costs a few milliseconds
* `--scoreonly` (`hirschberg_score()`): no alignment, no matrices. With a
  band, the forward and reverse passes run at the same time and meet in the
  middle row; otherwise one forward pass over every row
* On x86-64 Linux the SIMD code is compiled for several CPU levels and the
  best is picked when the program starts (see Build)
* Smith-Waterman hits are ordered with a counting sort on score instead of a
  comparison sort of every cell. Tied hits are ordered by position in the
  first sequence, then the second, on every platform
* Following Smith-Waterman hits back through the matrices uses the same
  substitution table and fixed gap penalties instead of a general step that
  re-checks every scoring option. Options that change penalties or forbid
  gaps for particular cells use the general step. A walk stops as soon as
  its next cell is already used by an earlier hit, before working out the
  step
* `smith_waterman` never indexes cells scoring below `--minscore`, and stops
  at `--maxhits` without searching for one more hit

Most of the remaining Smith-Waterman time is still spent following hits back
through the matrices. Each hit marks the cells it uses, which defines which
later hits are reported, so the walk itself is unchanged.

Comparison with other aligners
------------------------------

The same machine, with whole-program wall time (best of three runs) and peak
memory. Every tool used the same scoring (match +1, mismatch -2, gap of
length k: -(4+k)) and reported the same optimal score for every pair. The
sequences are random DNA: "1%" has 1% substitutions and 0.2% insertions and
deletions, "10%" has 10% and 2%, and "unrelated" is two independent random
sequences. seq-align uses up to 16 threads and the other tools one, so CPU
time is shown in brackets where it is well above the wall time.

The tools:
* seq-align before this fork: commit dc41988 (2017), the last upstream commit
* parasail 2.6.2: `nw_trace_striped_32`, `nw_striped_32` and
  `sw_trace_striped_32`
* ksw2 (289609b): `ksw_extz2_sse`, with no band
* WFA2-lib 2.3.5: exact (heuristics off), in the default memory mode and as
  BiWFA (ultralow memory)
* SSW 1.2.5

seq-align was run through `needleman_wunsch` and `smith_waterman`, and the
others through a small C program that reads the two sequences and aligns them
once.

Global alignment with the alignment returned, 20 kb:

| Tool | 1% | 10% | Unrelated | Peak memory |
|---|---|---|---|---|
| seq-align before this fork | 9.55 s | 9.75 s | 10.56 s | 4.6 GB |
| seq-align, full matrices | 0.37 s (2.1 s) | 0.43 s (2.3 s) | 0.41 s (2.2 s) | 4.6 GB |
| seq-align `--linear` | 0.05 s | 0.12 s | 0.12 s (0.73 s) | 7 MB |
| parasail | 1.17 s | 1.28 s | 0.76 s | 1.5 GB |
| ksw2 | 0.33 s | 0.35 s | 0.34 s | 524 MB |
| WFA2 | 0.01 s | 0.90 s | 23.3 s | 23 MB - 16.9 GB |
| WFA2 BiWFA | 0.01 s | 0.54 s | 10.3 s | 8 - 24 MB |

Global alignment, 50 kb (full matrices would need ~30 GB, so the original and
the full-matrix aligner were not run):

| Tool | 1% | 10% | Unrelated | Peak memory |
|---|---|---|---|---|
| seq-align `--linear` | 0.18 s (0.25 s) | 0.43 s (2.1 s) | 0.56 s (8.4 s) | 7 - 9 MB |
| parasail | 10.1 s | 8.74 s | 5.25 s | 9.5 GB |
| ksw2 | 1.77 s | 1.78 s | 1.71 s | 2.8 GB |
| WFA2 | 0.07 s | 5.23 s | out of memory (28 GB limit) | 87 MB - 6.2 GB |
| WFA2 BiWFA | 0.05 s | 3.07 s | 57.3 s | 7 - 44 MB |

Global alignment score only (every tool under 40 MB):

| Tool | 20 kb 1% | 20 kb 10% | 20 kb unrelated | 50 kb 1% | 50 kb 10% | 50 kb unrelated |
|---|---|---|---|---|---|---|
| seq-align `--scoreonly` | 0.01 s | 0.02 s (0.07 s) | 0.02 s (0.25 s) | 0.05 s (0.10 s) | 0.11 s (1.54 s) | 0.10 s (1.45 s) |
| parasail | 0.18 s | 0.17 s | 0.07 s | 1.14 s | 1.00 s | 0.45 s |
| ksw2 | 0.10 s | 0.10 s | 0.10 s | 0.63 s | 0.62 s | 0.61 s |
| WFA2 | <0.01 s | 0.46 s | 8.47 s | 0.03 s | 2.74 s | 49.5 s |

Local alignment, best hit only (`smith_waterman --maxhits 1`), 20 kb:

| Tool | 1% | 10% | Unrelated | Peak memory |
|---|---|---|---|---|
| seq-align before this fork | 71.1 s | 61.2 s | 34.5 s | 6.3 - 8.5 GB |
| seq-align | 3.03 s (4.6 s) | 3.32 s (4.9 s) | 0.65 s (2.1 s) | 4.7 - 5.6 GB |
| parasail | 0.90 s | 0.83 s | 0.72 s | 1.5 GB |
| SSW | 0.27 s | 0.22 s | 0.03 s | 7 MB |

For the unrelated pair seq-align reports no hit, as the best (score 14) is
below the default `--minscore`. At 5 kb every global aligner took 0.07 s or
less, except seq-align before this fork (0.53-0.63 s) and WFA2 on the
unrelated pair (0.65-0.96 s).

In short:
* Against the code before this fork, the same full-matrix alignment is
  18-30x faster, and `--linear` is up to ~190x faster on similar sequences
  while using megabytes instead of gigabytes
* For global alignment `--linear` has the lowest wall time of the exact
  aligners tested except WFA2 on near-identical sequences. WFA2's time grows
  with the number of differences, so it is far slower on divergent or
  unrelated pairs
* Per CPU core, ksw2 and parasail do more work than seq-align on divergent
  and unrelated sequences: seq-align's lead there comes from threads
* Local alignment is the weak spot: on related 20 kb pairs SSW finds the
  best hit 11-15x faster and parasail 3-4x faster. seq-align keeps the full
  matrices and indexes every cell above `--minscore` so that it can report
  further hits, which a best-hit-only search does not need

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
        --scoreonly          Print only the score (fast, linear memory)
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
`--nogaps..`, `--nomismatches`, `--printmatrices` or a positive `--gapopen`.
`--scoreonly` falls back to the full matrices for those scoring options.
When several alignments share the best score, `--linear` may return a
different one from the default aligner.


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
* Striped SIMD matrix fill (Farrar 2007, Bioinformatics 23:156-161), for
  the full matrices and the linear-space aligner

Algorithms to investigate:
* Myers' bit-vector algorithm (1999, J ACM 46:395-415). Only computes unit
  cost edit distance, so it would be a separate fast path rather than a
  replacement for the scored aligners
* Combining Myers' and Hirschberg's algorithms for linear-space
  edit-distance alignment
