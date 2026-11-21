# Nucleics Align

**Fast, exact pairwise sequence alignment in C:** optimal local
(Smith-Waterman) and global (Needleman-Wunsch) alignment with affine gaps,
from the command line, Python, Perl or C.

[![CI](https://github.com/DanielTillett/nucleics-align/actions/workflows/ci.yml/badge.svg)](https://github.com/DanielTillett/nucleics-align/actions/workflows/ci.yml)
[![Licence: BSD 3-Clause](https://img.shields.io/badge/licence-BSD%203--Clause-blue.svg)](LICENSE)
![C11](https://img.shields.io/badge/C-C11-555.svg)
![Python 3.7+](https://img.shields.io/badge/python-3.7%2B-3776ab.svg)
![Platforms: Linux | macOS](https://img.shields.io/badge/platforms-Linux%20%7C%20macOS-lightgrey.svg)

Nucleics Align is derived from Isaac Turner's
[seq-align](https://github.com/noporpoise/seq-align), with bug fixes, SIMD
and multi-threaded fills, exact banding and a linear-space global aligner.
Every speed-up returns the same optimal alignments as the plain dynamic
programming it replaces.

## Highlights

- **Exact.** Optimal alignments under your scoring, the same as a full
  dynamic-programming fill, with ties broken the same documented way on
  every code path
- **Fast.** A 20 kb global alignment with traceback takes 0.01-0.05 s, and
  the best local hit 0.01-0.12 s, where seq-align took 10-70 s
  ([benchmarks](docs/performance.md))
- **Long sequences.** `--linear` aligns megabase pairs in megabytes of
  memory: 1 Mb at 1% divergence in 6.3 s and 49 MB on an 8-core desktop
- **Any alphabet.** DNA, RNA, protein or any ASCII text, with built-in
  BLOSUM and PAM matrices, your own matrix files and wildcards
- **Batch friendly.** Reads FASTA, FASTQ, SAM, BAM, plain and gzipped files,
  aligning pairs on several threads with output in input order
- **Three ways in.** Command line tools, a Python module, Perl modules and a
  C library

## Contents

- [Quick start](#quick-start)
- [Building](#building)
- [Command line](#command-line)
- [Python](#python)
- [Perl](#perl)
- [C library](#c-library)
- [Performance](#performance)
- [Citing](#citing)
- [Contributing](#contributing)
- [Licence](#licence)

## Quick start

```console
$ git clone --recursive https://github.com/DanielTillett/nucleics-align
$ cd nucleics-align
$ make
$ ./bin/needleman_wunsch --printscores --pretty CAGACGT CGATA
C-AGACGT
| |*|
CGATA---
score: -11
```

## Building

You need a C11 compiler (gcc or clang), `make` and zlib (`zlib1g-dev` on
Debian and Ubuntu; included with macOS).

| Command | Does |
|---|---|
| `make` | Build the tools in `bin/`, the library `src/libalign.a` and the examples |
| `make test` | Run the unit tests and the command line fixtures |
| `make SANITIZE=1 test` | The same with AddressSanitizer and UBSan |
| `make NATIVE=1` | Tune for this CPU (e.g. AVX2 or AVX-512) |
| `make DEBUG=1` | No optimisation, with debug symbols |
| `make LANES=16 test` | Force 4, 8 or 16 SIMD lanes, to test other widths |

> [!TIP]
> Without `NATIVE=1`, x86-64 Linux builds already pick AVX2, SSE4.1 or plain
> SSE2 when the program starts, which gets most of the speed in a portable
> binary. `NATIVE=1` adds `-march=native` (16-lane AVX-512 where the CPU has
> it), and the binary may not run on other CPUs.

To update a checkout:

```console
$ git pull
$ git submodule update --init
```

<details>
<summary>Platform details</summary>

With gcc, every run-time version uses 8 lanes, run as two halves without
AVX2; clang does not allow that, so its versions use 4. AVX-512 is only used
with `NATIVE=1`: every version of a function shares one lane count, and 16
lanes run 10x slower without AVX-512. Plain SSE2 has no vector max
instruction for 32-bit integers, so on such CPUs the linear-space aligner
uses scalar loops. Other platforms use the compiler's target: Arm NEON
(e.g. Apple silicon) gets all the SIMD code, while x86 builds elsewhere
(e.g. macOS on Intel) can add `OPT="-O3 -msse4.1"` or `NATIVE=1`.
ThreadSanitizer builds skip the run-time selection, which it cannot handle.

`LANES=4`, `8` or `16` forces the lane count, so the tests cover widths the
CPU would not pick (slower, but the same results). It is not for clang on
x86-64 Linux, whose run-time selection takes only 4 lanes. CI runs 8 and 16
lanes on Arm.

Bundled dependencies, in `libs/` under their own licences:
[string_buffer](https://github.com/noporpoise/string_buffer),
[seq_file](https://github.com/noporpoise/seq_file) and
[sort_r](https://github.com/noporpoise/sort_r).

</details>

## Command line

```console
$ ./bin/smith_waterman --maxhits 1 --pretty ACGTTTGACGTACGTAGC TTTTACGTTTGACCTACGTAG
== Alignment 0 lengths (18, 21):

hit 0.0 score: 30
  ACGTTTGACGTACGTAG  [pos: 0; len: 17]
  |||||||||*|||||||
  ACGTTTGACCTACGTAG  [pos: 4; len: 17]

==
```

| Task | Command |
|---|---|
| Global alignment with scores | `needleman_wunsch --printscores seqA seqB` |
| Long sequences, linear memory | `needleman_wunsch --linear --printscores --file pairs.fa` |
| Score only (fastest) | `needleman_wunsch --scoreonly --file pairs.fa` |
| Best local hit | `smith_waterman --maxhits 1 --file pairs.fa.gz` |
| Protein | `needleman_wunsch --scoring BLOSUM62 --gapopen -10 --gapextend -1 ...` |
| Free end gaps (overlaps) | `needleman_wunsch --freestartgap --freeendgap ...` |
| Common substrings | `lcs ACGTACGTTACGACGT` |

The [command line guide](docs/command-line.md) has every option with
examples, and [Choosing scores](docs/scoring.md) has NCBI's guidance on
matrices and penalties.

> [!NOTE]
> `--linear` returns the same optimal score as the default aligner but,
> where several alignments share the best score, it may pick a different
> one.

## Python

[`python/nucleics_align.py`](python/nucleics_align.py) runs the tools from
Python 3.7 or later, with only the standard library. A batch of pairs goes
to one run of a tool and is aligned on its threads, in order.

```python
import sys
sys.path.insert(0, "python")
from nucleics_align import NeedlemanWunsch, SmithWaterman

nw = NeedlemanWunsch(gap_open=-5, threads=4)
aln = nw.align("ACGTTGCA", "ACGTGCA")
print(aln.align_a, aln.markers, aln.align_b, aln.score, sep="\n")

pairs = [("ACGT", "AGT"), ("GATTACA", "GCATGCA")]
for aln in nw.align_pairs(pairs):           # any iterable of pairs
    print(aln.score)
print(list(nw.score_pairs(pairs)))          # scores only, linear memory

for hit in SmithWaterman(max_hits=3).align("ttacgtacgtaa", "acgtacgt"):
    print(hit.score, hit.pos_a, hit.align_a, hit.align_b)
```

Every tool option has a keyword: `match`, `mismatch`, `gap_open`,
`gap_extend`, `scoring`, `substitution_matrix`, `substitution_pairs`,
`wildcards={"N": 0}`, `case_sensitive`, `free_start_gap`, `free_end_gap`,
`linear`, `linear_fallback`, the `no_gaps*` and `no_mismatches`
restrictions, `threads`, and for Smith-Waterman `min_score` and `max_hits`.

- A tool's error (a letter pair the scoring does not define, options linear
  space refuses) raises `AlignmentError` with the tool's message, after
  returning the pairs before it
- Sequences with whitespace, or starting with `>`, raise `ValueError`
- The tools are found in `bin/` beside `python/`, or set
  `NUCLEICS_ALIGN_BIN`
- Examples: [`nw_example.py`](python/nw_example.py),
  [`sw_example.py`](python/sw_example.py). Tests:
  `python3 -m unittest discover python`

## Perl

[`perl/NeedlemanWunsch.pm`](perl/NeedlemanWunsch.pm) and
[`perl/SmithWaterman.pm`](perl/SmithWaterman.pm) take the same options
under the tools' names:

```perl
use lib 'perl';
use NeedlemanWunsch;
use SmithWaterman;

my $nw = NeedlemanWunsch->new(gapopen => -5, linear => 1);
my $aln = $nw->do_alignment('ACGTTGCA', 'ACGTGCA');
$nw->print_alignment($aln);
print $nw->do_score('ACGTTGCA', 'ACGTGCA'), "\n";   # score only

my $sw = SmithWaterman->new(minscore => 5, maxhits => 3);
$sw->do_alignment('ttacgtacgtaa', 'acgtacgt');
$sw->print_hit($_) for $sw->get_all_hits();
```

They wait without a time limit unless `timeout` is set, number alignments
from 0, refuse unknown options and scoring names, and stop their tools when
the object goes. The Perl `SmithWaterman` reads hits one at a time, as the
tool's interactive mode gives them, where `minscore` defaults to 0 (every
hit); the Python module and the tool reading a file default to
`match * MAX(0.2 * length, 2)`. Tests: `perl perl/test.pl`.

## C library

Link against `src/libalign.a`, `libs/string_buffer/libstrbuf.a`, `-lpthread`
and `-lz`:

| Function | Header | Does |
|---|---|---|
| `needleman_wunsch_align()` | `needleman_wunsch.h` | Global alignment |
| `smith_waterman_align()`, `smith_waterman_fetch()` | `smith_waterman.h` | Local alignments, best first |
| `hirschberg_align()` | `hirschberg.h` | Global alignment in linear memory |
| `hirschberg_score()` | `hirschberg.h` | The global score alone, in linear memory |
| `scoring_init()`, `scoring_system_*()` | `alignment_scoring.h` | Scoring systems |

Each aligner has a `*_checked()` version that returns errors instead of
ending the program. The [library guide](docs/library.md) covers settings,
threads, tie-breaking and errors, and [`examples/`](examples) has complete C
and C++ programs.

## Performance

Global alignment of 20 kb random DNA with the alignment returned (whole
program wall time, best of three; AMD Ryzen AI Max+ 395; every tool reported
the same optimal score):

| Tool | 1% divergence | 10% | Unrelated | Peak memory |
|---|---|---|---|---|
| seq-align (original) | 9.55 s | 9.75 s | 10.56 s | 4.6 GB |
| **Nucleics Align** | **0.01 s** | **0.05 s** | **0.05 s** | 13 - 388 MB |
| **Nucleics Align `--linear`** | 0.05 s | 0.12 s | 0.08 s | **7 MB** |
| parasail | 1.17 s | 1.28 s | 0.76 s | 1.5 GB |
| ksw2 | 0.33 s | 0.35 s | 0.34 s | 524 MB |
| WFA2 | 0.01 s | 0.90 s | 23.3 s | 23 MB - 16.9 GB |

Nucleics Align used up to 16 threads and the others one; on one thread
Nucleics Align still takes 0.11-0.19 s against ksw2's 0.34-0.35 s.
[docs/performance.md](docs/performance.md) has the full comparison
(score-only, 50 kb, megabase pairs, best local hit), how each speed-up
works, and the scripts in [`benchmarks/`](benchmarks) reproduce the
timings.

## Citing

If you use Nucleics Align in published work, please cite it with the
details in [`CITATION.cff`](CITATION.cff) (GitHub's **Cite this
repository** button gives them in APA and BibTeX), and the papers behind
the methods you rely on, listed in [Algorithms and
references](docs/algorithms.md).

## Contributing

Bug reports and feature requests are welcome in the
[issue tracker](https://github.com/DanielTillett/nucleics-align/issues).
See [CONTRIBUTING.md](CONTRIBUTING.md) for building, testing and the
project's rule that speed-ups never change output, and
[CHANGELOG.md](CHANGELOG.md) for what has changed from seq-align.

## Licence

Nucleics Align is released under the [BSD 3-Clause License](LICENSE).

It is derived from seq-align, written by Isaac Turner and released into the
public domain (Unlicense). The bundled libraries in `libs/` keep their own
licences. The substitution matrices in `scoring/` come from NCBI
(<ftp://ftp.ncbi.nih.gov/blast/matrices/>).
