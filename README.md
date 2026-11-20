Nucleics Align
==============
Smith-Waterman and Needleman-Wunsch alignment in C
url: https://github.com/DanielTillett/nucleics-align
original author: Isaac Turner <turner.isaac@gmail.com>
license: BSD 3-Clause (see [LICENSE](LICENSE))
updated: 19 Nov 2026

[![CI](https://github.com/DanielTillett/nucleics-align/actions/workflows/ci.yml/badge.svg)](https://github.com/DanielTillett/nucleics-align/actions/workflows/ci.yml)

About
=====

C implementations of optimal local (Smith-Waterman) and global
(Needleman-Wunsch) alignment with affine gap penalties. The command line
tools `smith_waterman`, `needleman_wunsch` and `lcs` cover most uses, and the
library is easy to call from other programs (see `examples/`). Perl wrapper
modules are in `perl/`.

Nucleics Align is a fork of Isaac Turner's
[seq-align](https://github.com/noporpoise/seq-align)
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

    $ git clone --recursive https://github.com/DanielTillett/nucleics-align
    $ cd nucleics-align
    $ make

Build options:

    $ make test                           # run the tests
    $ make SANITIZE=1 test                # tests with ASan and UBSan
    $ make NATIVE=1                       # tune for this CPU (e.g. AVX2)
    $ make DEBUG=1                        # no optimisation, debug symbols
    $ make LANES=16 test                  # tests with 4, 8 or 16 SIMD lanes

`make test` runs the library's unit tests, then `tests/cli/run.sh`, which
compares the tools' output, error messages and exit status with the
expected files in `tests/cli/expected`, running each batch case on one
thread and on four. The unit tests print their random seed; set
`SEQ_ALIGN_TEST_SEED` to repeat a run. Changing the compiler, flags or any
header rebuilds what depends on it, so a build never mixes objects compiled
different ways. `benchmarks/` has the scripts behind the timings below.

`NATIVE=1` adds `-march=native`, which lets the SIMD code use 16-lane
AVX-512 or 8-lane AVX2 on x86 (4 lanes otherwise). Binaries built this way
may not run on other CPUs.

Without it, x86-64 Linux builds with gcc or clang compile the SIMD code for
AVX2, SSE4.1 and plain SSE2, and pick one for the CPU when the program
starts. That gets most of the `NATIVE=1` speed in a portable binary (see
Performance). With gcc every version uses 8 lanes, run as two halves
without AVX2; clang does not allow that, so its versions use 4. AVX-512 is
only used with `NATIVE=1`: every version of a function shares one lane
count, and 16 lanes run 10x slower without AVX-512. Plain SSE2 has no vector
max instruction for 32-bit integers, so on such CPUs the linear-space
aligner uses scalar loops. Other platforms
use the compiler's target: ARM NEON (e.g. Apple silicon) gets all the SIMD
code, while x86 builds elsewhere (e.g. macOS on Intel) can add
`OPT="-O3 -msse4.1"` or `NATIVE=1`. ThreadSanitizer builds skip the run-time
selection, which it cannot handle.

`LANES=4`, `8` or `16` forces the lane count, so the tests can cover widths
the CPU would not pick (slower, but the same results). It is not for clang
on x86-64 Linux, whose run-time selection takes only 4 lanes. The CI runs
8 and 16 lanes on Arm.

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

* `needleman_wunsch_align()` (`src/needleman_wunsch.h`): global alignment.
  Set `compact` on the aligner to keep one traceback byte per cell instead
  of the three score matrices (12x less memory, same alignment); the
  matrices are then not available, e.g. to `alignment_print_matrices()`.
  `aligner_t` gained fields for this, so code built against an older
  `alignment.h` must be recompiled. Compact fills of similar sequences skip
  cells that cannot be on a best path (`alignment_band_enabled`, default
  true), with the same alignment
* `smith_waterman_align()` then `smith_waterman_fetch()` repeatedly
  (`src/smith_waterman.h`): local alignments, best first.
  `smith_waterman_set_min_score()` skips hits below a score: the hits
  returned are the same, but lower-scoring cells are never indexed.
  `smith_waterman_set_compact()` suits callers that mostly want the best
  hit: it is found from traceback bytes, and the score matrices are only
  filled if a second hit is fetched. Above `smith_waterman_rect_min_cells`
  (default 2^24), and where no gap scores above 0, only the rectangle of
  the best hit gets traceback bytes, unless the hit spans the sequences;
  above `smith_waterman_strip_min_cells` (default 2^28) the hit is followed
  back from rows kept every so often instead, in O(n sqrt(m)) memory. The
  hits are the same either way
* `hirschberg_align()` (`src/hirschberg.h`): global alignment in O(n) memory;
  same scores as `needleman_wunsch_align()`, but returns false for options it
  does not support (a positive gap open, free end gaps, `no_gaps_in_*`,
  `no_mismatches`). Large sub-problems run their forward and reverse passes
  at the same time, each split into column blocks on several threads
  (tunable with `hirschberg_block_min_cells`, `hirschberg_block_min_cols`
  and `hirschberg_max_threads`, shared by the two passes). Passes skip cells
  that cannot be on a best
  path (`hirschberg_band_enabled`, default true), and sub-problems reuse rows
  saved by the passes above them (`hirschberg_reuse_depth`, default 2; 0
  turns it off), with the same result
* `hirschberg_score()` (`src/hirschberg.h`): only the global alignment score,
  in O(n) memory; same options as `hirschberg_align()`
* `scoring_init()` and the `scoring_system_*()` functions
  (`src/alignment_scoring.h`): set up a scoring system

Settings that change speed but never a result (threads, SIMD, banding,
8- and 16-bit lanes, block sizes) are in an `alignment_options_t`
(`src/alignment.h`). Give an aligner its own with its `options` field
(Needleman-Wunsch) or `smith_waterman_set_options()`, and pass them to
`hirschberg_align_opts()`/`hirschberg_score_opts()`; threads with different
options can then align at the same time. Without options, the process-wide
defaults apply: `alignment_options_default()` reads them from the
`alignment_*`, `hirschberg_*` and `smith_waterman_rect_min_cells`
variables, which keep their old meaning (e.g. `alignment_striped_enabled =
false` switches SIMD off everywhere).

One setting can change a result, and is off by default: with
`linear_fallback` (`alignment_linear_fallback`), a Needleman-Wunsch
alignment too big for the score matrices or traceback bytes (too long, or
out of memory) is made in linear space instead, where that supports the
scoring. The score is the same, but it may choose another of several
equal-scoring alignments. `needleman_wunsch --linearfallback` turns it on.

Where several alignments score the same, which one is returned is fixed
and the same on every path of the full-matrix aligners (score matrices,
traceback bytes, bands, difference-encoded lanes, and bands or
Smith-Waterman hits kept in strips, which reproduce the same bytes):
* The traceback runs from the end of the alignment back to its start.
  Needleman-Wunsch starts in the matrix of the last cell that holds the
  best score, preferring GAP_A (a gap in the first sequence, from the cell
  above), then GAP_B (a gap in the second, from the cell on the left), then
  MATCH. Smith-Waterman starts in MATCH at the hit's last cell.
* At each step it moves to the first predecessor, in the same order GAP_A,
  GAP_B, MATCH, whose score plus the step's score gives the current score.
  So in GAP_A an extension wins over an opening, but in GAP_B an opening
  from GAP_A wins over an extension; from MATCH the diagonal cell's gap
  states win over its MATCH, which tends to put a gap as late in the
  alignment as it can go (`AAA` against `AA` gives `AA-`, `ACCA` against
  `ACA` gives `AC-A`). The `nw_ties` and `nw_free_gaps` command line
  fixtures pin such cases
* Smith-Waterman hits come out by score, highest first, then by lowest
  position in the first sequence, then in the second. A hit stops at a
  cell scoring 0, and one that reaches a cell an earlier hit used is
  dropped
* `--linear` (and `hirschberg_align()`) gives the same score but may pick
  another of the equal-scoring alignments: it splits each sub-problem at
  the first best crossing of its middle row (a crossing in MATCH or a gap
  before one inside a deletion), and in a one-row sub-problem puts the
  first sequence's letter at the first best column, or a deletion next to
  whichever neighbouring deletion it joins

The aligners above end the program on an error. The `*_checked()` versions
(`needleman_wunsch_align_checked()`, `smith_waterman_align_checked()`,
`smith_waterman_fetch_checked()`, `hirschberg_align_checked()`,
`hirschberg_score_checked()`) return it instead, as an `aln_status_t`.
Errors in the input are found before any work: a character pair the scoring
does not define (`ALN_ERR_SCORING`, naming the pair), sequences too long to
address (`ALN_ERR_TOO_LONG`), scores that could pass `INT_MAX`
(`ALN_ERR_OVERFLOW`) or, for linear space, options it refuses
(`ALN_ERR_UNSUPPORTED`). Running out of memory (`ALN_ERR_NOMEM`) is found as
it happens, on any of the alignment's threads: what was allocated is freed,
there is no alignment or hit to read, and the aligner can be used again.
Fetching a second Smith-Waterman hit can run out of memory too, as it fills
the score matrices if a compact alignment has not. Creating an aligner or a
result (`needleman_wunsch_new()`, `alignment_create()` and so on) still ends
the program if out of memory. `alignment_fprint_error()` prints the message
the aligners would have. Scores should lie within plus or minus
`SCORE_ABS_MAX` (`INT_MAX / 4`); the command line tools and the scoring file
loaders refuse others.

`alignment_colour_fprint_against()` and `alignment_fprint_spacer()` print to
a given stream instead of stdout. `scoring_make_case_sensitive()` switches a
scoring system to case-sensitive comparison, keeping its scores for upper
case too.

See `examples/` for complete C and C++ programs.

The command line tools use half the CPUs by default (`--threads <n>` to
change it). Pairs read from a file are aligned several at a time, on
separate threads, and printed in input order, so the output is the same as
one at a time. Each pair is printed as soon as those before it are. A pair
of more than 4 million cells runs on its own with all the threads. An error
in a pair (a character pair the scoring does not define, or scores big
enough to overflow) is reported after the output of the pairs before it,
and stops the program as before.

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
| 5 kb | 0.01 s, 26 MB | 0.01 s, 2 MB | <0.01 s, 2 MB |
| 20 kb | 0.03-0.08 s, 19-388 MB | 0.05-0.15 s, 7 MB | 0.01-0.03 s, 7 MB |
| 50 kb | 0.15-0.41 s, 92 MB-2.4 GB | 0.18-0.63 s, 7 MB | 0.05-0.18 s, 7 MB |

The full-matrix aligner keeps one traceback byte per cell, and for similar
sequences only near the main diagonals (see below); with `--printmatrices`
it keeps the score matrices, 12 bytes per cell (5 kb 289 MB, 20 kb 4.6 GB).
By divergence, full-matrix aligner, 50 kb: 1% 0.15 s and 92 MB, 5% 0.16 s
and 411 MB, 10% 0.18 s and 744 MB, 25% or unrelated 0.39-0.41 s and 2.4 GB
(the whole matrix); 20 kb: 1% 0.03 s and 19 MB, otherwise 0.07-0.08 s and
388 MB.

`--linear` and `--scoreonly` are fastest for similar sequences, where most
of the matrix can be skipped (see banding below). For 20 kb and 50 kb pairs
with substitutions and equal numbers of insertions and deletions:

| Divergence | 20 kb `--linear` | 20 kb `--scoreonly` | 50 kb `--linear` | 50 kb `--scoreonly` |
|---|---|---|---|---|
| 1% substitutions, 0.2% indels | 0.05 s | 0.01 s | 0.18 s | 0.05 s |
| 5%, 1% | 0.09 s | 0.03 s | 0.32 s | 0.06 s |
| 10%, 2% | 0.11 s | 0.03 s | 0.43 s | 0.10 s |
| 25%, 5% | 0.15 s | 0.04 s | 0.63 s | 0.18 s |
| unrelated | 0.08 s | 0.01 s | 0.33 s | 0.07 s |

Before banding the same runs took 0.11-0.22 s, 0.02-0.05 s, 0.53-1.06 s and
0.11-0.27 s, fastest for the unrelated pair.

Measured on 32 CPUs; large fills and passes use up to half of them. Threads
wait for each other by spinning, so CPU time can be several times the wall
time (50 kb unrelated sequences, `--linear`: 0.33 s wall, 2.1 s CPU).

Pairs from a file are aligned several at a time (see Using the library):
3,000 pairs of 0.2-2 kb took 0.50 s instead of 4.93 s (Needleman-Wunsch),
0.69 s instead of 7.48 s (`--linear`), 0.31 s instead of 2.54 s
(`--scoreonly`) and 0.63 s instead of 5.05 s (`smith_waterman --maxhits
1`), with the same output.

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
* The command line tools keep one traceback byte per cell instead of three
  4-byte scores (as ksw2 does). For each of the three matrices the byte
  records which matrix the traceback steps back to, found during the fill
  with the same tests in the same order as the usual traceback, so the
  alignment is identical. Most of the full-matrix time went on allocating
  and first touching the matrices: 20 kb Needleman-Wunsch 0.36 -> 0.08 s and
  4.6 GB -> 0.39 GB, and 50 kb (~30 GB before) now takes 0.4 s and 2.4 GB.
  Options the striped fill does not handle, and free start gaps, keep the
  score matrices
* `smith_waterman --maxhits 1` uses the same traceback bytes, with the fill
  also finding the best cell (highest score, then lowest position in the
  first sequence, then the second: the first hit of a full search). Only if
  a second hit is wanted are the score matrices filled and every cell
  indexed. 20 kb, best hit: 2.5-3.0 s -> 0.09-0.11 s, 5.6 GB -> 0.45 GB
* The traceback-byte fill uses the linear-space aligner's exact banding: a
  narrow band is filled first, in strips of rows, and widened (or the whole
  matrix filled) if the score found shows a best path could leave it. No
  cell outside a band that held can tie with a cell on a best path, so the
  traceback takes the same steps. Only the band's bytes are kept. 50 kb at
  1% divergence: 0.40 -> 0.15 s, 6.0 -> 0.14 s of CPU, 2.4 GB -> 92 MB.
  A band is used only where it beats the whole fill on its threads
* For a compact Smith-Waterman best hit of more than 2^24 cells, a
  score-only scan finds the best score and where it ends, and a scan of the
  reversed sequences up to there finds where every best path starts. Only
  that rectangle gets traceback bytes (as SSW does, but with the same hit
  as the full search). 100 kb sequences sharing a 2 kb region: 3.1 -> 0.34
  s and 11.6 GB -> 22 MB; 300 kb ones take 3.3 s and 66 MB (they did not fit
  in 40 GB before). The first scan gives up once the best rises at half the
  top substitution score per row from near the start: such a hit spans the
  sequences, and the whole fill is quicker. Not used where a gap scores
  above 0, as a hit can then start with a gap
* Smith-Waterman fills and scans use 16-bit lanes, twice as many per
  vector, where the scores allow (every cell lies between 0 and the best
  score so far), and fill each row in one pass: match, both gaps and the
  row's best cell together, with the lazy F loop adding what each lane
  passes to the next. A row whose best comes near the 16-bit limit, or a
  best rising fast enough to reach it, starts the fill again with 32-bit
  lanes. One thread, 20 kb: the score scan 0.125 -> 0.056 s (unrelated)
  and 0.29 -> 0.14 s (1%), the best hit with traceback bytes 0.33 -> 0.22 s
  and 0.49 -> 0.29 s; with `NATIVE=1` on AVX-512 the unrelated scan takes
  0.025 s (16 Gcells/s). 3,000 local best hits: 0.50 -> 0.33 s
* Smith-Waterman starts with 8-bit lanes, 4 times as many as 32-bit ones
  (as SSW does, but signed and with the same results), where every
  substitution score and gap penalty fits in a byte: then 16-bit lanes,
  then 32-bit, as the best grows. Lanes cannot hold a segment number here,
  so a row that reaches the best finds its first cell by searching the lane
  that holds it. Between unrelated sequences the best rises ever more
  slowly, so the check for a best that will outgrow the lanes projects from
  its rise over the last half of the rows, not from the start. One thread,
  20 kb unrelated: the score scan 0.057 -> 0.017 s (24 Gcells/s; 0.014 s
  with `NATIVE=1` on AVX-512). With the scan this cheap, compact best hits
  of more than 2^24 cells (was 2^30) start with it: 20 kb unrelated
  `smith_waterman --maxhits 1` 0.09 -> 0.01 s and 448 -> 5 MB, with related
  pairs as before
* Needleman-Wunsch traceback bytes are filled with difference-encoded lanes
  (after Suzuki & Kasahara 2018, as ksw2 does) where the whole matrix is
  filled: each cell is kept as its difference from the cell above, which
  stays small at any length, so 8-bit lanes (4 times as many as 32-bit
  ones) hold them with DNA scoring and 16-bit lanes with BLOSUM62. Every
  traceback byte comes from comparing values with the same reference, so
  the bytes, and the alignment, are those of the 32-bit fill. More lanes
  made the lazy F loop wrap more (gap_b carried across 4 times as many lane
  boundaries: 20 kb at 10%, 12 sweeps of each row instead of 3), so what
  each lane passes to the next is found first, by a scan over the lanes
  in absolute scores, and one pass then applies it. One thread, 20 kb fill:
  10% 0.40 -> 0.20 s, unrelated 0.28 -> 0.18 s. Bands are filled in strips
  of at least 64 rows (was 256), which makes a band 6-30% quicker and a
  first band that proves too narrow cheaper
* A fill across the whole width takes column blocks of at least 2,048
  columns (twice `block_min_cols`), so a 20 kb pair runs on about 10 threads
  instead of 16: the extra threads finished each row too soon after the
  block on their left and mostly waited. The same wall time with 40-45% less
  CPU (20 kb: Needleman-Wunsch at 10% 0.83 -> 0.48 s of CPU, the best local
  hit at 1% 0.92 -> 0.53 s). A band's strips keep 1,024, or narrow bands
  lose their threads (50 kb at 10%: 0.15 -> 0.29 s)
* Before banding, 128 words of 12 letters of the first sequence are looked
  for near where they would lie in the second, for a high guess at the best
  score. If even the band that guess needs could not pay, the band is not
  tried: 20 kb at 10-20% divergence and unrelated 0.048-0.053 -> 0.040-0.047
  s, 50 kb at 20% and unrelated about 10% quicker. The guess only decides
  whether to try a band, never the alignment
* Before that, with one match and one mismatch score and nothing else
  changed (no free end gaps, wildcards or gap rules), a score-only
  wavefront alignment (WFA, see Algorithms below) looks for the best score
  itself, for up to 16 diagonal updates per letter
  (`alignment_band_wfa_budget`). Its time grows with the square of the
  divergence, so it finishes only for close pairs, and then the band that
  score needs is filled once, instead of a 64-wide band first and wider
  ones after. Otherwise the guess and the narrow band go on as before. The
  band proof is unchanged, so the alignment is the same either way. 0.1%
  divergence: 50 kb 0.040 -> 0.018 s, 100 kb 0.098 -> 0.060 s, 300 kb
  0.51 -> 0.41 s; 20 kb at 0.5% 0.017 -> 0.011 s. Where it gives up it
  costs 0-5% (50 kb at 0.5%: 0.069 -> 0.072 s). Larger budgets let it
  finish at 1 Mb and 0.1%, but the exact band was then slower there (4.5
  -> 5.3 s), likely because the narrower band runs on fewer threads
* Once a band has a lower bound on the best score (the wavefront score, or
  the score of the alignment a band too narrow found), each strip covers
  only the diagonals that the cells above it can still reach while
  scoring that much: a cell's bound is its score plus a match for every
  letter left of the shorter remainder and a gap extension for the
  difference. A cell whose bound is below the best score is on no best
  path and cannot tie with one, so leaving it out changes no cell the
  traceback reads, and the alignment is the same. This is the
  computational volume of A*PA2 (see Algorithms below), there for edit
  distance, here with affine gaps and only between strips. On the 8-core
  machine below, best of 3: 1 Mb at 1% 13.6 -> 9.9 s, at 2% 26.4 -> 20.6
  s, at 5% 67.6 -> 55.0 s; 300 kb at 2% 2.71 -> 2.11 s; 100 kb at 1% 0.69
  -> 0.54 s; 50 kb at 0.5% 0.142 -> 0.100 s. 100 kb at 2% and 300 kb at
  0.5% were within 4% either way
* Long global alignments: a band whose traceback bytes would not fit in a
  quarter of the physical memory (`alignment_band_keep_max_cells`) keeps
  only one strip of rows' bytes at a time, and each strip's top row. The
  traceback fills a strip again from its top row when it reaches it, only
  up to the column it has reached, so the bytes and the alignment are the
  same. The band is filled about twice, so smaller bands keep every byte.
  1 Mb pairs at 1% divergence: 34.7 GB -> 350 MB (13.8 -> 18.2 s); at 2%
  and 5% they no longer run out of memory (see Long pairs below)
* Strip alignments: a compact best hit whose region (the rectangle, or the
  whole matrix where a hit spans the sequences) has more than 2^28 cells
  keeps no traceback bytes for it. The scan for the best cell also keeps
  every k-th row of the three matrices (k about sqrt(12 x rows)), and the
  hit is followed back strip by strip from the bottom, each strip's bytes
  filled again from the row kept above it and only up to the hit's column.
  The bytes are those of the whole fill, so the hit is too (a linear-space
  realignment would give the same score but could pick another of several
  equal-scoring alignments). 50 kb pairs sharing the whole sequence:
  1.27 -> 0.62 s and 2.9 GB -> 85 MB; 20 kb: 451 -> 23 MB, and slightly
  faster. Smaller regions keep the threaded fill of every byte, which is
  quicker there (10 kb: 0.03 s against 0.05 s)
* The linear-space aligner's forward and reverse passes share the threads
  rather than taking half the CPUs each, which oversubscribed the cores (50
  kb unrelated `--linear` 0.45 -> 0.33 s and 6.1 -> 2.1 s of CPU)
* The command line tools align pairs from a file several at a time, each
  on its own thread with its own aligner, and print them in input order
  (3,000 pairs of 0.2-2 kb: 8-11x faster)
* gcc builds with run-time selection use 8-lane vectors in every version,
  which also helps without AVX2: fewer segments per row, so the lazy F loop
  wraps less (20 kb unrelated `--scoreonly`, one thread, 0.20 -> 0.10 s).
  `NATIVE=1` on an AVX-512 CPU uses 16 lanes (the same pass 0.10 -> 0.055 s)
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
* Rows reused down the linear-space recursion: each split runs a forward
  pass over the top half of its rows and a reverse pass over the bottom
  half, and its first half then needs a forward pass from the same corner
  (its second half a reverse pass). A forward row depends only on the rows
  above it and the columns to its left, so the split's own row, cut to the
  half's columns, is exactly the row the half would make. A pass over all
  columns now also saves the rows the next two levels on its side will
  need (`hirschberg_reuse_depth`), and those levels skip that pass. Banded
  passes save nothing, as their rows are exact only along best paths. The
  alignments are the same. 400 kb pairs at 30% divergence 11.8 -> 9.5 s,
  unrelated 9.1 -> 7.5 s; 50 kb unrelated 0.31 -> 0.25 s; similar pairs,
  which are mostly banded, 1-5% faster. Peak memory is about the same (at
  most four more rows per level kept). Saving three levels was slower
  (10.2 s), as each pass is cut into more pieces
* The lazy F loop finds its carries first. Along a band, every cell right
  of the best path is best reached by an insertion gap from the left, so
  Farrar's loop carried each row's gap round all the lanes, one lane per
  sweep: about 8 sweeps of the segment per row against 1 for the fill
  itself. What lane l passes to lane l+1 is now found with a scan over the
  lanes (the better of its own and what it was passed, less a gap
  extension per cell crossed), as Block Aligner does with prefix scans
  (Liu & Steinegger 2023, Bioinformatics 39:btad487), and each lane takes
  its carry once. The values are the same, so are the alignments. 100 kb
  pairs: `--linear` 1% 0.50 -> 0.25 s and 5% 0.84 -> 0.38 s; `--scoreonly`
  1% 0.17 -> 0.08 s; the banded fill 1% 0.35 -> 0.25 s. Unrelated pairs
  and Smith-Waterman, where the carries rarely reach far, are unchanged
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

When more than one hit is wanted, most of the Smith-Waterman time is still
spent following hits back through the matrices. Each hit marks the cells it
uses, which defines which later hits are reported, so the walk itself is
unchanged.

Comparison with other aligners
------------------------------

The same machine, with whole-program wall time (best of three runs) and peak
memory. Every tool used the same scoring (match +1, mismatch -2, gap of
length k: -(4+k)) and reported the same optimal score for every pair. The
sequences are random DNA: "1%" has 1% substitutions and 0.2% insertions and
deletions, "10%" has 10% and 2%, and "unrelated" is two independent random
sequences. Nucleics Align uses up to 16 threads and the other tools one, so CPU
time is shown in brackets where it is well above the wall time.

The tools:
* seq-align before this fork: commit dc41988 (2017), the last upstream commit
* parasail 2.6.2: `nw_trace_striped_32`, `nw_striped_32` and
  `sw_trace_striped_32`
* ksw2 (289609b): `ksw_extz2_sse`, with no band
* WFA2-lib 2.3.5: exact (heuristics off), in the default memory mode and as
  BiWFA (ultralow memory)
* SSW 1.2.5

Nucleics Align was run through `needleman_wunsch` and `smith_waterman`, and the
others through a small C program that reads the two sequences and aligns them
once. The Nucleics Align rows were measured again after each change that
affects them (most recently the difference-encoded lanes, best of 5); the
other rows are from the earlier run, on the same sequences and machine.

Global alignment with the alignment returned, 20 kb:

| Tool | 1% | 10% | Unrelated | Peak memory |
|---|---|---|---|---|
| seq-align before this fork | 9.55 s | 9.75 s | 10.56 s | 4.6 GB |
| Nucleics Align, full matrix | 0.01 s | 0.05 s (0.8 s) | 0.05 s (0.7 s) | 13 - 388 MB |
| Nucleics Align `--linear` | 0.05 s | 0.12 s | 0.08 s (0.29 s) | 7 MB |
| parasail | 1.17 s | 1.28 s | 0.76 s | 1.5 GB |
| ksw2 | 0.33 s | 0.35 s | 0.34 s | 524 MB |
| WFA2 | 0.01 s | 0.90 s | 23.3 s | 23 MB - 16.9 GB |
| WFA2 BiWFA | 0.01 s | 0.54 s | 10.3 s | 8 - 24 MB |

Global alignment, 50 kb (the original would need ~30 GB, so was not run):

| Tool | 1% | 10% | Unrelated | Peak memory |
|---|---|---|---|---|
| Nucleics Align, full matrix | 0.10 s | 0.15 s (1.3 s) | 0.26 s (3.9 s) | 71 MB - 2.4 GB |
| Nucleics Align `--linear` | 0.18 s (0.23 s) | 0.41 s (1.2 s) | 0.31 s (1.8 s) | 7 MB |
| parasail | 10.1 s | 8.74 s | 5.25 s | 9.5 GB |
| ksw2 | 1.77 s | 1.78 s | 1.71 s | 2.8 GB |
| WFA2 | 0.07 s | 5.23 s | out of memory (28 GB limit) | 87 MB - 6.2 GB |
| WFA2 BiWFA | 0.05 s | 3.07 s | 57.3 s | 7 - 44 MB |

Global alignment score only (every tool under 40 MB):

| Tool | 20 kb 1% | 20 kb 10% | 20 kb unrelated | 50 kb 1% | 50 kb 10% | 50 kb unrelated |
|---|---|---|---|---|---|---|
| Nucleics Align `--scoreonly` | 0.01 s | 0.02 s (0.06 s) | 0.01 s (0.13 s) | 0.05 s (0.09 s) | 0.10 s (0.66 s) | 0.07 s (0.84 s) |
| parasail | 0.18 s | 0.17 s | 0.07 s | 1.14 s | 1.00 s | 0.45 s |
| ksw2 | 0.10 s | 0.10 s | 0.10 s | 0.63 s | 0.62 s | 0.61 s |
| WFA2 | <0.01 s | 0.46 s | 8.47 s | 0.03 s | 2.74 s | 49.5 s |

Long pairs, global alignment with the alignment returned
(`needleman_wunsch`), measured on another machine (AMD Ryzen 7 5700G, 8
cores, 16 threads, 62 GB, Linux), one run each, on 21 Nov 2026 (after the
lazy F scan and row reuse, see above):

| Pair | Band, every byte kept | Band kept in strips | `--linear` |
|---|---|---|---|
| 100 kb, 1% | 0.35 s, 268 MB | (kept) | 0.27 s, 6 MB |
| 300 kb, 5% | 2.77 s, 11.6 GB | (kept) | 2.47 s, 20 MB |
| 1 Mb, 0.1% | 3.46 s, 2.6 GB | (kept) | 2.86 s, 45 MB |
| 1 Mb, 0.5% | (in strips) | 5.36 s, 287 MB | 4.19 s, 47 MB |
| 1 Mb, 1% | (in strips) | 7.92 s, 319 MB | 6.33 s, 49 MB |
| 1 Mb, 2% | out of memory | 15.8 s, 467 MB | 10.6 s, 53 MB |
| 1 Mb, 5% | out of memory | 43.3 s, 1.1 GB | 23.6 s, 63 MB |

With 62 GB, bands of up to about 15 GB keep every byte; larger ones are
kept in strips. `--linear` is now the fastest on every pair here (1.1 to
1.8 times the speed of the band) and needs far less memory: the lazy F
scan sped up its passes more than the band fill (1 Mb at 1%: 19.2 -> 6.3 s
against 9.9 -> 7.9 s). Its alignment can differ where several share the
best score (see above). The score alone (`--scoreonly`) took 0.09-11.8 s.

Local alignment, best hit only (`smith_waterman --maxhits 1`), 20 kb:

| Tool | 1% | 10% | Unrelated | Peak memory |
|---|---|---|---|---|
| seq-align before this fork | 71.1 s | 61.2 s | 34.5 s | 6.3 - 8.5 GB |
| Nucleics Align | 0.12 s (0.9 s) | 0.11 s (0.8 s) | 0.01 s (0.14 s) | 5 - 23 MB |
| parasail | 0.90 s | 0.83 s | 0.72 s | 1.5 GB |
| SSW | 0.27 s | 0.22 s | 0.03 s | 7 MB |

For the unrelated pair Nucleics Align reports no hit, as the best (score 14) is
below the default `--minscore`. The Nucleics Align row was measured again
after the 8-bit lanes and strip alignments were added (best of 5; the
version before them took 0.13 s, 0.14 s and 0.09 s with 448-451 MB in the
same session). With `NATIVE=1` on AVX-512 the unrelated pair takes under
0.005 s. At 50 kb Nucleics Align's best hit takes 0.62 s (6.1 s CPU) and 85
MB when it spans both sequences (1.27 s and 2.9 GB before strip
alignments), and 0.06 s and 10 MB for the unrelated pair (the other tools
were not run). At 5 kb every global aligner
took 0.07 s or less, except seq-align before this fork (0.53-0.63 s) and
WFA2 on the unrelated pair (0.65-0.96 s).

In short:
* Against the code before this fork, the same full-matrix alignment is
  140-320x faster with 12-240x less memory, and the best local hit
  550-3,400x faster with over 250x less memory
* For global alignment Nucleics Align has the lowest wall time of the exact
  aligners tested, except WFA2 on near-identical sequences: 0.01-0.05 s for
  the full-matrix aligner at 20 kb against 0.33-0.35 s for ksw2. WFA2's time
  grows with the number of differences, so it is far slower on divergent or
  unrelated pairs. `--linear` needs only megabytes at any length
* On one core Nucleics Align is now faster than ksw2 too: `needleman_wunsch
  --threads 1` takes 0.11 s on the 20 kb 10% pair and 0.19 s on the
  unrelated one, against ksw2's 0.35 s and 0.34 s. The CPU times in brackets
  above come from the default 16 threads, which mostly wait for their
  neighbours (16 threads are only 4x faster than one); `--threads 4` takes
  0.07 s and 0.28 s of CPU on the unrelated pair. On similar pairs banding
  skips most cells (20 kb at 1%: 0.01 s of CPU)
* For the best local hit Nucleics Align is now faster than parasail and SSW
  in wall time, on unrelated pairs too (0.01 s and 5 MB against SSW's 0.03 s
  and 7 MB), and needs 23 MB where a 20 kb hit spans the sequences (SSW 7
  MB, parasail 1.5 GB). SSW still uses about 3x less CPU time on related
  pairs. On one core Nucleics Align's score scan of the unrelated pair takes
  0.017 s (0.014 s with `NATIVE=1` on AVX-512)

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

        --threads <n>        Threads to use, within an alignment and for aligning
                             pairs from a file at once [default: half the CPUs]

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

      https://github.com/DanielTillett/nucleics-align  (compiled: <date>)


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

        --threads <n>        Threads to use, within an alignment and for aligning
                             pairs from a file at once [default: half the CPUs]


        --freestartgap       No penalty for gap at start of alignment
        --freeendgap         No penalty for gap at end of alignment

        --printscores        Print optimal alignment scores
        --zam                A funky type of output
        --linear             Align in linear memory (Hirschberg, Myers-Miller)
        --scoreonly          Print only the score (fast, linear memory)
        --linearfallback     Align in linear memory where too big for the matrices
                             (same score, but may pick another equal-scoring
                             alignment)
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

      https://github.com/DanielTillett/nucleics-align  (compiled: <date>)

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

Nucleics Align is licensed under the BSD 3-Clause License; see [LICENSE](LICENSE).

It is derived from seq-align, originally written by Isaac Turner and released
into the public domain (Unlicense). The bundled libraries in `libs/` keep
their own licences. The substitution matrices in `scoring/` come from NCBI
(ftp://ftp.ncbi.nih.gov/blast/matrices/).


Changes that affect output
==========================

Output is otherwise the same as earlier versions, byte for byte. These fix
bugs:

* `--case_sensitive` now makes comparisons case sensitive with
  match/mismatch scoring and with substitution files. It used to change only
  how substitution files were loaded, so their upper-case entries were
  never found. Built-in scoring systems (`--scoring`) keep their scores for
  both cases; only a character against the other case of itself changes
* `--linear` scores sequences containing `-` or `\0` as the full-matrix
  aligner does (it used to read a `-` in the alignment as a gap)
* Empty number arguments (e.g. `--match ""`, once read as 0), scores beyond
  plus or minus 536870911, and alignments whose scores could overflow are
  refused with an error, instead of giving wrong scores
* Global alignments set the result's `pos_a`, `pos_b`, `len_a` and `len_b`
  (to the whole sequences), where they kept the values of an earlier local
  alignment
* The error for sequences too long to address reads `Error: Sequences too
  long to align`, without a source file and line in front
* Running out of memory while aligning prints `Error: Out of memory`,
  without a source file and line, and the tools stop the same way as for
  other errors (in batch mode, after the output of the pairs before)

Development
===========

Bug reports and feature requests:
https://github.com/DanielTillett/nucleics-align/issues

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
  the full matrices and the linear-space aligner, with the lazy F carries
  found by a scan over the lanes, as Block Aligner uses prefix scans (Liu &
  Steinegger 2023, Bioinformatics 39:btad487). Our own code
* Wavefront alignment (WFA: Marco-Sola et al. 2021, Bioinformatics
  37:456-463), score only, to size the band of a global alignment, with
  match scores converted to penalties as Eizenga & Paten describe (2022,
  bioRxiv 10.1101/2022.01.12.476087, a preprint). Our own code, written
  from the papers
* Computational volumes, as in A*PA2 (Groot Koerkamp 2024, WABI, LIPIcs
  312:17), adapted from edit distance to affine gaps, to narrow each strip
  of a band. Our own code, written from the paper

Algorithms to investigate:
* Myers' bit-vector algorithm (1999, J ACM 46:395-415). Only computes unit
  cost edit distance, so it would be a separate fast path rather than a
  replacement for the scored aligners
* Combining Myers' and Hirschberg's algorithms for linear-space
  edit-distance alignment
