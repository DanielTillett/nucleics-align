# Changelog

Nucleics Align is derived from Isaac Turner's
[seq-align](https://github.com/noporpoise/seq-align) (last commit dc41988,
2017). Output is the same as seq-align's, byte for byte, except for the
fixes listed under [Changes to output](#changes-to-output).

## Unreleased

### Added

- Striped SIMD fills (8-, 16- and 32-bit lanes, x86 and Arm) for both
  aligners, with run-time selection of AVX2, SSE4.1 or plain SSE2 on x86-64
  Linux
- Global alignment in linear memory (`needleman_wunsch --linear`,
  `--linearfallback`) and the score alone (`--scoreonly`)
- Exact banding, one traceback byte per cell, and bands kept in strips for
  long global alignments; strip alignments and rectangles for the best local
  hit
- Threads within an alignment and across pairs from a file (`--threads`),
  with output in input order
- `*_checked()` library functions that return errors instead of ending the
  program, and per-aligner settings (`alignment_options_t`)
- `--nogapsinlongest`
- A Python module, `python/nucleics_align.py`, with examples and tests
- Tests: unit tests, command line fixtures, and CI on Linux (x86 and Arm)
  and macOS

See [docs/performance.md](docs/performance.md) for what each change gained.

### Changed

- The public header is `nucleics_align.h` (was `seq_align.h`), with
  `NUCLEICS_ALIGN_*` macros; the test binary is `bin/nucleics_align_tests`
  and its seed variable `NUCLEICS_ALIGN_TEST_SEED`
- The Perl modules take every tool option, add
  `NeedlemanWunsch::do_score()` and `SmithWaterman::get_all_hits()`, wait
  without a time limit unless `timeout` is set, and refuse unknown options
  and scoring names
- Licensed under the BSD 3-Clause License (seq-align was released under the
  Unlicense)

### Fixed

- The Perl modules passed mixed-case options twice, accepted any scoring
  name, numbered alignments from -1, scored empty sequences with a formula
  that ignored free end gaps, and left stopped tools as zombie processes
- macOS builds failed to compile `alignment_cmdline.c`

### Changes to output

These fix bugs:

- `--case_sensitive` now makes comparisons case sensitive with
  match/mismatch scoring and with substitution files. It used to change only
  how substitution files were loaded, so their upper-case entries were
  never found. Built-in scoring systems (`--scoring`) keep their scores for
  both cases; only a character against the other case of itself changes
- `--linear` scores sequences containing `-` or `\0` as the full-matrix
  aligner does (it used to read a `-` in the alignment as a gap)
- Empty number arguments (e.g. `--match ""`, once read as 0), scores beyond
  plus or minus 536870911, and alignments whose scores could overflow are
  refused with an error, instead of giving wrong scores
- Global alignments set the result's `pos_a`, `pos_b`, `len_a` and `len_b`
  (to the whole sequences), where they kept the values of an earlier local
  alignment
- The error for sequences too long to address reads `Error: Sequences too
  long to align`, without a source file and line in front
- Running out of memory while aligning prints `Error: Out of memory`,
  without a source file and line, and the tools stop the same way as for
  other errors (in batch mode, after the output of the pairs before)
- `--freeendgap` makes the whole of a sequence aligned against an empty one
  free (score 0), as `--freestartgap` already did: that gap ends the
  alignment as well as starting it. It used to cost a gap
