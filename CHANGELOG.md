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
  tuning and test hooks. Band strips split into smaller thread blocks
  (`alignment_band_block_min_cols`, `alignment_band_block_min_cells`).
  1 Mb at 1%: 6.2 -> 3.1 s on 8 threads with a third of the CPU; 200 kb at
  1%: 0.62 -> 0.30 s and 657 -> 128 MB (same alignments)

### Fixed

- Builds without the SIMD kernels (a compiler without
  `__builtin_shufflevector`, or `-DALN_NO_SIMD`, new) compile and link
  again: `alignment_grow_capacity()` and the `ALN_TARGET_CLONES` fallback
  were declared only under `ALN_STRIPED`
- `hirschberg_simd_available()` takes the options to answer for (NULL: the
  defaults); it read the process-wide `alignment_striped_enabled` whatever
  the aligner's options said
- The unchecked aligners `needleman_wunsch_align()`,
  `hirschberg_align_opts()` and `hirschberg_score_opts()` are now the
  `*_checked()` forms plus the printed error and `exit()`, so they make
  the same checks: the linear-space pair gains the sequence-length guard,
  and Needleman-Wunsch reports a pair the scoring lacks with the aligners'
  message instead of `scoring_lookup()`'s. `smith_waterman_align()` keeps
  skipping the checks on purpose: the test suite compares the full,
  compact and rectangle fills through it under gaps scoring 0 or above,
  which the checked form refuses
- Perl: a tool still running at a read timeout is terminated before it is
  waited for (the modules could hang in `waitpid`), and a write to a tool
  that has gone fails quietly instead of `SIGPIPE` ending the caller
- `benchmarks/compare_harness.sh`, `journal_v130_extra.sh`, `ablation.c`
  and `percore.c` built against the `libs/` tree removed in v2.0
- The Makefile probed the compiler for `-falign-loops` once per recipe
  (`OPT` is now `:=`); `.gitignore` covers the `src/lib` and `src/cli`
  dependency files
- README: the 1 Mb 0.1% `--wfa` figure now matches `docs/performance.md`
  (160.7 ms / 11 MiB)

### Changed

- One parser for the scoring file formats, in the library
  (`scoring_parse_matrix()`, `scoring_parse_pairs()`,
  `scoring_parse_strerror()` in `alignment_scoring.h`): the public API's
  `na_scoring_parse_*()` and the tools' `--substitution_matrix` and
  `--substitution_pairs` read through it. The tools now read the whole
  (gzipped) file first, and their loads are transactional too; a matrix
  file with headings and no rows is refused ("No rows in matrix") instead
  of being accepted empty. The public API's error messages take the tools'
  wording. `align_scoring_load_*()` take the case setting from the scoring
- The banded compact fill (`alignment_fill_banded()`) runs as its phases:
  the cost model (`band_cost_init()`, `band_first_width()`), the plan of a
  band's strips (`band_plan()`), their room (`band_alloc()`) and the strip
  fill (`band_fill_strips()`), and `seeds_build()` as its (classes, seed
  length, the index of b, the seeds' copies, the chains), with one
  interpolation for sampled rows. Same alignments, timed the same
- The linear-space passes and the wavefront aligner pass argument
  structs (`hb_block_args_t`, `hb_pass_args_t`; `wf_dir_t`, `biwfa_ctx_t`,
  `biwfa_sub_t`, `biwfa_meet_t`) where 15 to 20 arguments went before; the
  kernels copy the fields to locals at entry. Same output; timed the same
  on 200 kb linear-space and banded fills and on wavefront alignment
- The tools' option parser is a table: one entry per option with its
  name, how many values it takes, which tool takes it and a handler, so
  the "only valid with" and "requires an argument" checks are made once,
  in one place. The same options, values, messages and exit statuses as
  before (every command line fixture is unchanged)
- The unit tests are one file per subject (`tests/test_*.c`, twelve
  suites) on a shared harness (`tests/test_common.h`): one random
  generator reseeded for every test from the run's seed and the test's
  name, so `NUCLEICS_ALIGN_TEST_SEED` repeats a run and a failing test
  repeats on its own; the library's tuning globals saved before and
  restored after each test, with a note for any a test left changed; one
  random sequence generator, mutator and alignment comparison in place of
  three or four each; a failed assertion names its test. The v1.3.0
  targeted tests (`tests/test_v130.c`) are folded into the suites. The
  unit tests no longer link the command line objects or zlib: the
  scoring-file loader and the option parsers' ranges are pinned by the
  command line fixtures instead, which gain cases for `--help`,
  `--fastlocal`, `--wfabudget`, `--threads 0` and out-of-range values
- `NA_VERSION_STR` is built from `NA_VERSION_MAJOR/MINOR/PATCH`, and
  `NA_VERSION_NUMBER` (0x020000) is new, for `#if` checks. The unused
  `NUCLEICS_ALIGN_VERSION` and `NUCLEICS_ALIGN_VERSION_STR` macros are gone,
  as are the unused `scoring_print()`, `alignment_colour_print_against()`
  and `alignment_print_spacer()` (the `f`-variants taking a `FILE *` stay)
- `scoring_try_lookup()` (new) returns false for a pair the scoring does not
  cover; `scoring_lookup()` wraps it and exits as before. The public API no
  longer reaches an `exit()` through `na_scoring_lookup()` or
  `na_score_alignment()`
- Internal tidying with no change in output: `sub_table_t` records its
  highest and lowest score (`smax`, `smin`) once, in place of five loops;
  `alignment_options_resolve()`, `alignment_lane_bits()`,
  `alignment_striped_applies()`, `band_strip_rows()`, `sw_trace_step()` and
  `sw_describe()` replace three to five copies each; the scoring-text
  parsers share one failure path; the three `na_*` aligner calls share one
  argument check

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
