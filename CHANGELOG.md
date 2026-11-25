# Changelog

Nucleics Align is derived from Isaac Turner's
[seq-align](https://github.com/noporpoise/seq-align) (last commit dc41988,
2017). Output is the same as seq-align's, byte for byte, except for the
fixes listed under [Changes to output](#changes-to-output).

## Unreleased

### Improved

- Band strips are narrowed further by seeds (`src/lib/alignment_seeds.c`):
  words of `a` with no exact copy in `b` each cost an error, and matched
  words on different diagonals cost gap letters, the gap-chaining seed
  heuristic of A\*PA (Groot Koerkamp & Ivanov 2024) with affine scores as
  penalties. Strips also get shorter as they get narrower. The alignment is
  unchanged; `alignment_seeds_enabled` and `alignment_seed_length` are the
  tuning and test hooks

## 1.3.0 - 24 Nov 2026

### Added

- `needleman_wunsch --wfabudget <n>` (`alignment_wfa_budget`): configurable
  diagonal-update cap for `--wfa` BiWFA, separate from `alignment_band_wfa_budget`
- `smith_waterman --fastlocal` (`smith_waterman_fast_rect`): opt-in best-hit
  reconstruction for single-hit compact runs (`--maxhits 1`) using forward and
  reverse $O(|B|)$-memory scans to locate a joint optimal start/end pair,
  followed by BiWFA or linear-space alignment and exact local score + strict
  positive-prefix verification (automatic fallback to exact Gotoh on mismatch;
  default Gotoh tie-breaking unchanged)
- Ukkonen exponential diagonal banding in `alignment_myers_score()` and
  `alignment_myers_score_banded()`, plus full-row `alignment_myers_row()` wired
  into single-threaded unit-cost Hirschberg (`--linear`)
- `alignment_sw_scan_hint()` to skip narrow 8/16-bit lane retries when a lower
  bound on the optimal score is already known

### Improved

- BiWFA (`src/alignment_wfa.c`): fixed `wf_compute_interior` pointer-offset UB
  and candidate bounds-check order, added reusable all-`WFA_NONE` null
  wavefront for empty sources, added leading/trailing `WFA_NONE` trimming
  (`wf_trim`), `gcd(px, po, pe)` penalty normalisation, 32-byte sentinel
  padding on forward and reversed buffers, runtime-dispatched 32-byte AVX2
  match extension, and empty-ring bailout (~2.0–2.3x faster on 20 kb–1 Mb DNA)
- Multi-threaded wave synchronisation (`alignment.c`, `hirschberg.c`): replaced
  unbounded `sched_yield()` spin-waits with bounded 512-`pause` spin +
  per-boundary condition-variable waits, cutting 8-thread full-matrix CPU time
  by ~25% while keeping wall time fast
- Small-input setup: cached online CPU count and physical RAM quarter, added
  early single-thread band rejection before WFA/k-mer probing, and made
  Hirschberg SIMD scratch allocation lazy

## 1.2.0 - 22 Nov 2026

### Added

- $O(s)$-space Bidirectional Wavefront Alignment (BiWFA) with decoupled
  AVX2/SSE4.1/SSE2 multi-versioned compute and extend loops for `--wfa`
  (`alignment_wfa_align()`), lifting the retained-wavefront memory ceiling so
  long similar pairs (50 kb--1 Mb at ~1% divergence) align in a few megabytes
  of RAM while returning the exact optimal score

### Improved

- Eliminated cache-line false sharing on multi-threaded block progress counters
  (`_Alignas(64)` in `alignment.c` and `hirschberg.c`), cached observed row
  progress, batched publication to every 16 rows, and added x86 `pause` hints
  to cut spin-wait CPU overhead on multi-threaded full-matrix fills
- `smith_waterman --maxhits 1` strip refills now reuse the precomputed
  substitution table and select 8- or 16-bit SIMD lanes whenever the known
  optimal score fits safely, preserving byte-identical traceback
- Capped traceback buffer growth above 1 MiB at 1.25x instead of power-of-two
  doubling to avoid over-allocating multi-gigabyte buffers

## 1.1.0 - 22 Nov 2026

### Added

- Unit-cost scores (match 0, mismatch and gap letter -1, no gap open
  penalty) on one thread, from `--scoreonly` or `hirschberg_score()`, use
  Myers' bit-vector algorithm where the band does not pay: the same score,
  2.5 times as fast for unrelated pairs
- `needleman_wunsch --linearabove <n>` (`linear_above` in
  `alignment_options_t`): pairs where either sequence has at least n
  letters are aligned in linear space, the faster and far smaller choice
  for long pairs, others with the matrices. Off by default, as linear space
  may choose another of several equal-scoring alignments
- `needleman_wunsch --wfa` (`wfa_align` in `alignment_options_t`,
  `alignment_wfa_align()`): aligns similar pairs by the wavefront algorithm
  where the scoring allows (one match and one mismatch score, with gap
  penalties), falling back to the matrices otherwise. Off by default, as WFA
  may choose another of several equal-scoring alignments

### Fixed

- Fixed memory-safety and edge-case bugs across WFA, Hirschberg, striped
  banding, Smith-Waterman fetch, and CLI parsing (including OOM handling in
  `alignment_wfa_align()`)

### Documentation

- Rewrote README and `docs/` guides with five SVG diagrams and updated
  high-resolution benchmarks up to 1 Mb (including `--wfa`)

## 1.0.0 - 21 Nov 2026

The first release of Nucleics Align.

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

- The version header is `nucleics_align.h` (was `seq_align.h`), with
  `NUCLEICS_ALIGN_*` macros; the test binary is `bin/nucleics_align_tests`
  and its seed variable `NUCLEICS_ALIGN_TEST_SEED`
- `aligner_t` gained fields (for one traceback byte per cell), so code
  built against seq-align's `alignment.h` must be recompiled
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
