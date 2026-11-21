# Contributing to Nucleics Align

Bug reports, questions and pull requests are welcome. Please open an
[issue](https://github.com/DanielTillett/nucleics-align/issues) before a
large change, so we can agree the approach first.

## The one rule: speed-ups never change output

Every optimisation must give exactly the same scores and alignments as the
code it replaces, byte for byte, including which of several equal-scoring
alignments is returned (see [tie-breaking](docs/library.md#equal-scoring-alignments)).
A change that alters output on purpose is a bug fix, and goes in
[CHANGELOG.md](CHANGELOG.md) under *Changes to output*.

## Building and testing

```console
$ git clone --recursive https://github.com/DanielTillett/nucleics-align
$ cd nucleics-align
$ make test                      # unit tests, then the command line fixtures
$ make SANITIZE=1 test           # with AddressSanitizer and UBSan
$ make LANES=16 test             # other SIMD widths (4, 8 or 16)
$ python3 -m unittest discover python
$ perl perl/test.pl
```

- `make test` runs the library's unit tests (`bin/nucleics_align_tests`),
  then [`tests/cli/run.sh`](tests/cli/run.sh), which compares the tools'
  output, error messages and exit status with the files in
  `tests/cli/expected`, running each batch case on one thread and on four
- The unit tests print their random seed; set `NUCLEICS_ALIGN_TEST_SEED` to
  repeat a run
- Changing the compiler, flags or any header rebuilds what depends on it, so
  a build never mixes objects compiled different ways
- CI runs gcc and clang, with and without sanitisers, on Linux (x86 and Arm)
  and macOS, plus forced lane counts

## Pull requests

- Add a test that fails without your change: a unit test in
  [`src/tools/tests.c`](src/tools/tests.c), or a fixture in
  [`tests/cli/cases`](tests/cli/cases) for the tools' output
- For a speed-up, include timings before and after (best of several runs,
  with the machine), and check the output is unchanged against the previous
  build. [`benchmarks/`](benchmarks) has scripts for this
- Match the surrounding code's style and comment density
- Code from other projects needs its licence checked before it can be
  included; describe ideas taken from papers in
  [docs/algorithms.md](docs/algorithms.md) with the reference

## Documentation

User documentation is in [README.md](README.md) and [`docs/`](docs). It is
written in British English, with dates as DD MMM YYYY (e.g. 21 Nov 2026).
