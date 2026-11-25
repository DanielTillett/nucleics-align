# Performance

[Documentation index](README.md)

The point of forking `seq-align` was to make exact affine-gap alignment fast enough that you do not have to reach for a heuristic tool just to align two 50 kb or 1 Mb sequences. Every default optimisation here produces the same score and alignment as the original code (see [CHANGELOG.md](../CHANGELOG.md) for the few intentional bug fixes).

## How a global alignment is dispatched

When you run `needleman_wunsch`, the CLI dispatches a pair in this order (`src/cli/nw_cmdline.c`, `src/lib/hirschberg.c`, `src/lib/needleman_wunsch.c`):

```mermaid
flowchart TD
    A[Pair of sequences] --> S0{"--scoreonly?"}
    S0 -- yes --> HS0{"hirschberg_score_checked()"}
    HS0 -- ALN_ERR_UNSUPPORTED --> NW
    HS0 -- other error --> ERR["nw_error() (exit non-zero)"]
    HS0 -- ok --> MY{"Unit cost (0/-1/0/-1)?"}
    MY -- yes --> MYR["Ukkonen-banded Myers (64 cells/word)<br/>(MT striped fallback if divergent)"] --> DONE[Output]
    MY -- no --> HB{"Banded SIMD pass<br/>pays?"}
    HB -- yes --> DONE
    HB -- no --> FULLS["Full striped or scalar linear pass"] --> DONE

    S0 -- no --> L0{"--linear?"}
    L0 -- yes --> L1["hirschberg_align_checked()"]
    L1 -- ok --> DONE
    L1 -- error --> ERR["nw_error() (exit non-zero)"]

    L0 -- no --> NW["needleman_wunsch_align_checked()"]
    NW --> W0{"--wfa?"}
    W0 -- yes --> WA["Wavefront alignment (WFA)<br/>within cell budget"]
    WA -- finished --> DONE
    WA -- "unsupported / over budget" --> LA
    W0 -- no --> LA{"len >= --linearabove?"}
    LA -- yes --> LAH["hirschberg_align_checked()"]
    LAH -- ok --> DONE
    LAH -- "unsupported / error" --> C
    LA -- no --> C{Options handled by<br/>striped fill?}
    C -- no --> S[Scalar fill of 3 score matrices] --> T[Traceback]
    C -- yes --> W["Exact score from WFA (budgeted),<br/>or high guess from 12-mer matches"]
    W --> D{Does a band pay?}
    D -- yes --> E["Banded traceback-byte fill in strips<br/>(narrowed by A*PA2 volume and A*PA seeds)"]
    D -- no --> F["Whole traceback-byte fill<br/>(difference-encoded 8/16-bit lanes)"]
    E --> G{Band fits in 1/4 RAM?}
    G -- yes --> T
    G -- no --> K["Keep 1 strip at a time;<br/>refill strips during traceback"] --> T
    F --> T
    S -. "TOO_LONG / NOMEM<br/>with --linearfallback" .-> LFB["hirschberg_align_checked()"] --> DONE
    E -. "TOO_LONG / NOMEM<br/>with --linearfallback" .-> LFB
    F -. "TOO_LONG / NOMEM<br/>with --linearfallback" .-> LFB
    T --> DONE
```

<p align="center">
  <img src="figures/dp-banding-strips.svg" width="760" alt="Full DP matrix vs exact banding in strips vs linear space Hirschberg vs wavefront alignment">
</p>

## Benchmarks on Linux (`pop-os`, AMD Ryzen 7 5700G, v1.3.0)

Measured on an 8-core / 16-thread desktop (AMD Ryzen 7 5700G, 62 GiB RAM, Pop!_OS 24.04 Linux, GCC 13.3.0, default `make` build at `-O3`, default 8 threads, default DNA scoring `+1 / -2 / -4 / -1`, Nucleics Align v1.3.0 commit `2529b86`, archived in [`benchmarks/results/pop-os-v1.3.0-journal.txt`](../benchmarks/results/pop-os-v1.3.0-journal.txt)). Standard pairs (`5 kb – 50 kb`) are min wall / median RSS of 3 runs; large pairs (`100 kb – 1 Mb`) are 1 run each, timed with millisecond precision (`time.perf_counter()` + GNU `time -v` in binary MiB). **Bold** marks the fastest global mode returning a full alignment (excluding `--scoreonly` and `local --maxhits 1`).

<p align="center">
  <img src="figures/perf-pop-os.svg" width="760" alt="Global alignment wall time on pop-os (Ryzen 7 5700G) from 5 kb to 1 Mb">
</p>

### Standard lengths (5 kb – 50 kb, min of 3)

| Length | Div | `full` (wall / MiB) | `--wfa` (wall / MiB) | `--linear` (wall / MiB) | `--scoreonly` (wall / MiB) | `local --maxhits 1` (wall / MiB) |
|---|---|---:|---:|---:|---:|---:|
| 5 kb | 1% | 6.3 ms / 4 MiB | **2.8 ms** / 2 MiB | 7.1 ms / 2 MiB | 3.1 ms / 2 MiB | 23.5 ms / 30 MiB |
| 5 kb | 10% | 14.1 ms / 27 MiB | 23.2 ms / 27 MiB | **11.1 ms** / 3 MiB | 4.8 ms / 3 MiB | 21.6 ms / 31 MiB |
| 5 kb | rand | **12.4 ms** / 27 MiB | 21.0 ms / 27 MiB | 13.4 ms / 3 MiB | 4.9 ms / 2 MiB | 30.5 ms / 31 MiB |
| 20 kb | 1% | 25.6 ms / 16 MiB | **9.5 ms** / 3 MiB | 30.1 ms / 3 MiB | 8.7 ms / 3 MiB | 144.6 ms / 24 MiB |
| 20 kb | 10% | 71.0 ms / 387 MiB | 100.3 ms / 387 MiB | **70.1 ms** / 4 MiB | 15.6 ms / 4 MiB | 124.9 ms / 24 MiB |
| 20 kb | rand | **70.3 ms** / 387 MiB | 102.0 ms / 387 MiB | 76.4 ms / 4 MiB | 20.1 ms / 4 MiB | 153.0 ms / 24 MiB |
| 50 kb | 1% | 125.8 ms / 83 MiB | **39.4 ms** / 3 MiB | 105.4 ms / 4 MiB | 33.5 ms / 4 MiB | 618.8 ms / 87 MiB |
| 50 kb | 10% | 420.4 ms / 2,397 MiB | 498.1 ms / 2,397 MiB | **226.2 ms** / 7 MiB | 63.2 ms / 6 MiB | 550.1 ms / 86 MiB |
| 50 kb | rand | 402.1 ms / 2,395 MiB | 474.8 ms / 2,395 MiB | **277.9 ms** / 7 MiB | 96.8 ms / 6 MiB | 702.3 ms / 86 MiB |

### Large lengths (100 kb – 1 Mb, 1 run)

| Length | Div | `full` (wall / MiB) | `--wfa` (wall / MiB) | `--linear` (wall / MiB) | `--scoreonly` (wall / MiB) |
|---|---|---:|---:|---:|---:|
| 100 kb | 0.1% | 70.9 ms / 38 MiB | **9.9 ms** / 3 MiB | 211.8 ms / 7 MiB | 39.1 ms / 6 MiB |
| 100 kb | 1% | 476.7 ms / 308 MiB | **136.1 ms** / 4 MiB | 304.8 ms / 6 MiB | 108.3 ms / 6 MiB |
| 100 kb | 5% | **453.8 ms** / 1,474 MiB | 627.6 ms / 1,474 MiB | 482.0 ms / 9 MiB | 143.2 ms / 8 MiB |
| 300 kb | 0.1% | 374.1 ms / 282 MiB | **24.1 ms** / 5 MiB | 545.5 ms / 15 MiB | 126.3 ms / 15 MiB |
| 300 kb | 1% | **969.8 ms** / 2,693 MiB | 1,418.0 ms / 2,693 MiB | 1,170.1 ms / 16 MiB | 276.7 ms / 16 MiB |
| 300 kb | 5% | 3,562.8 ms / 206 MiB *(strips)* | 4,033.3 ms / 206 MiB | **2,870.9 ms** / 22 MiB | 1,073.0 ms / 20 MiB |
| 1 Mb | 0.1% | 4,041.0 ms / 3,125 MiB | **160.7 ms** / 11 MiB | 3,183.5 ms / 46 MiB | 1,072.8 ms / 44 MiB |
| 1 Mb | 1% | 9,349.3 ms / 342 MiB *(strips)* | 10,220.5 ms / 342 MiB | **7,293.5 ms** / 50 MiB | 2,658.4 ms / 47 MiB |
| 1 Mb | 5% | *(not run)* | *(not run)* | **95,338.5 ms** / 110 MiB | 58,836.2 ms / 94 MiB |

With v1.3.0 BiWFA (wavefront trimming, reusable null fronts, and 32-byte AVX2 match extension), `--wfa` aligns **1 Mb at 0.1% divergence in 160.7 ms and 11 MiB** (25.1x faster and 284x smaller than `full`, 19.8x faster than `--linear`, and 6.7x faster than `--scoreonly` while returning the full alignment), and **100 kb at 1% divergence in 136.1 ms and 4 MiB** (3.5x faster and 77x smaller than `full`). At 300 kb 1% and 1 Mb 1%, `--wfa` hits its default work budget (`min(2^28, 256 * (len_a + len_b))`, adjustable via `--wfabudget <N>`) and falls back to the matrix aligner.

For longer or more divergent pairs, `full` and `--linear` trade the lead (`full` is faster at 100 kb 5% and 300 kb 1%; `--linear` is faster at 50 kb 10%/rand, 300 kb 5%, and 1 Mb 1%), while `--linear` keeps peak RAM at **46–110 MiB** even at 1 Mb and 5% divergence. At 1 Mb and 1%, `full` automatically switches to strip-kept bands, dropping peak RAM to **342 MiB** while keeping fixed tie-breaking.

## Other workloads (protein, batching, 1 Mb strips)

- **Protein (2 kb, BLOSUM62, Ryzen AI Max+ 395):** `0.031 s` in original `seq-align` down to `0.005 s` now (same with `NATIVE=1`).
- **Batch threading (3,000 DNA pairs, 0.2–2 kb, Ryzen AI Max+ 395):** `4.93 s → 0.50 s` (`needleman_wunsch`), `7.48 s → 0.69 s` (`--linear`), `2.54 s → 0.31 s` (`--scoreonly`), and `5.05 s → 0.63 s` (`smith_waterman --maxhits 1`), with byte-identical output.
- **1 Mb strip-kept bands on `pop-os` (62 GB RAM, separate strip-memory run on 21 Nov 2026):** once a band would exceed a quarter of RAM (~15 GB here), `full` switches to keeping one strip at a time. That earlier run timed `full` (strips) against `--linear` across four divergences (which is why `1 Mb 5% full` was not repeated in the `--wfa` table above):
  - `0.5%` divergence: `5.36 s / 287 MB` (`--linear`: `4.19 s / 47 MB`)
  - `1%` divergence: `7.92 s / 319 MB` (`--linear`: `6.33 s / 49 MB`)
  - `2%` divergence: `15.8 s / 467 MB` (`--linear`: `10.6 s / 53 MB`)
  - `5%` divergence: `43.3 s / 1.1 GB` (`--linear`: `23.6 s / 63 MB`)

## Comparison with other exact aligners

Measured on an AMD Ryzen AI Max+ 395 (16 cores / 32 threads, WSL2, gcc 15.2 `-O3`), best of 3 runs, same scoring (`+1 / -2 / -4 / -1`), every tool returning the same optimal score. Nucleics Align uses up to 16 threads; the other tools are single-threaded. Because threads wait on neighbouring blocks by spinning, CPU time (shown in brackets where well above wall time) can be several times wall time.

### Global alignment, 20 kb

| Tool | 1% | 10% | Unrelated | Peak memory |
|---|---|---|---|---|
| seq-align (original) | 9.55 s | 9.75 s | 10.56 s | 4.6 GB |
| Nucleics Align, full matrix | 0.01 s | 0.06 s (0.57 s) | 0.06 s (0.58 s) | 11 – 386 MB |
| Nucleics Align `--linear` | 0.02 s | 0.04 s (0.07 s) | 0.04 s (0.12 s) | 2 MB |
| parasail 2.6.2 | 1.17 s | 1.28 s | 0.76 s | 1.5 GB |
| ksw2 | 0.33 s | 0.35 s | 0.34 s | 524 MB |
| WFA2 2.3.5 | 0.01 s | 0.90 s | 23.3 s | 23 MB – 16.9 GB |
| WFA2 BiWFA | 0.01 s | 0.54 s | 10.3 s | 8 – 24 MB |

Even restricted to a single thread (`--threads 1`), Nucleics Align takes **0.11 s** on the 20 kb 10% pair and **0.19 s** on the unrelated pair (against ksw2's 0.35 s and 0.34 s).

### Global alignment, 50 kb

Original `seq-align` would need ~30 GB here, so it was not run.

| Tool | 1% | 10% | Unrelated | Peak memory |
|---|---|---|---|---|
| Nucleics Align, full matrix | 0.06 s | 0.14 s (0.97 s) | 0.29 s (4.4 s) | 54 MB – 2.4 GB |
| Nucleics Align `--linear` | 0.07 s (0.09 s) | 0.17 s (0.46 s) | 0.18 s (1.0 s) | 3 – 6 MB |
| parasail 2.6.2 | 10.1 s | 8.74 s | 5.25 s | 9.5 GB |
| ksw2 | 1.77 s | 1.78 s | 1.71 s | 2.8 GB |
| WFA2 2.3.5 | 0.07 s | 5.23 s | OOM (28 GB) | 87 MB – 6.2 GB |
| WFA2 BiWFA | 0.05 s | 3.07 s | 57.3 s | 7 – 44 MB |

### Global alignment, score only

| Tool | 20 kb 1% | 20 kb 10% | 20 kb rand | 50 kb 1% | 50 kb 10% | 50 kb rand |
|---|---|---|---|---|---|---|
| Nucleics Align `--scoreonly` | <0.01 s | 0.01 s (0.03 s) | 0.01 s (0.07 s) | 0.02 s (0.03 s) | 0.03 s (0.25 s) | 0.06 s (0.67 s) |
| parasail | 0.18 s | 0.17 s | 0.07 s | 1.14 s | 1.00 s | 0.45 s |
| ksw2 | 0.10 s | 0.10 s | 0.10 s | 0.63 s | 0.62 s | 0.61 s |
| WFA2 | <0.01 s | 0.46 s | 8.47 s | 0.03 s | 2.74 s | 49.5 s |

### Best local hit, 20 kb (`smith_waterman --maxhits 1`)

| Tool | 1% | 10% | Unrelated | Peak memory |
|---|---|---|---|---|
| seq-align (original) | 71.1 s | 61.2 s | 34.5 s | 6.3 – 8.5 GB |
| Nucleics Align | 0.11 s (0.52 s) | 0.12 s (0.56 s) | 0.01 s (0.11 s) | 3 – 21 MB |
| parasail | 0.90 s | 0.83 s | 0.72 s | 1.5 GB |
| SSW 1.2.5 | 0.27 s | 0.22 s | 0.03 s | 7 MB |

## Where the speed comes from

In rough order of what made the biggest difference:

1. **One traceback byte per cell (12x less RAM):** Instead of storing three 4-byte score matrices, compact mode packs the three 2-bit predecessor directions into one byte per cell. Most of the old wall time was actually page-faulting gigabytes of RAM.
2. **Striped SIMD fills (Farrar 2007) + prefix-scan lazy F:** Written with portable GCC/clang vector extensions (AVX2/SSE4.1/SSE2 multi-versioned on x86-64 Linux, AVX-512 with `NATIVE=1`, NEON on Arm64). Precomputing inter-lane carries with a single scan cuts lazy-F sweeps along bands from ~8 per row down to 1.
3. **8- and 16-bit narrow / difference-encoded lanes:** Smith-Waterman starts in 8-bit lanes (4x width) and widens to 16- or 32-bit only if scores outgrow the lane. Needleman-Wunsch full traceback fills store differences from the cell above (Suzuki & Kasahara 2018), keeping values in 8 bits for DNA and 16 bits for BLOSUM62 at any sequence length. Linear-space passes use 16-bit relative blocks.
4. **Exact banding, WFA sizing, and A\*PA strip narrowing:** Before filling a global matrix, a budgeted WFA score pass (or 12-mer sample guess) estimates the band width. Each strip is narrowed using the exact upper bound from the row above, tightened by seeds: words of `a` with no exact copy in `b` must each hold an error, and matched words on different diagonals cost gap letters (the gap-chaining seed heuristic of A\*PA), so a strip covers only the diagonals a best path can use, and strips get shorter as they get narrower. If the band is too big for a quarter of RAM, only one strip's bytes are kept at a time and refilled on traceback.
5. **Bidirectional wavefront alignment (`--wfa`, `--wfabudget`):** On close pairs, trims empty wavefront diagonals, uses 32-byte AVX2 match extension (8-byte fallback), and finds breakpoints in $O(s)$ wavefront memory.
6. **Linear-space global alignment (`--linear`) with row reuse:** Hirschberg / Myers-Miller with exact banded passes, 16-bit blocks, two-level intermediate row reuse, and Myers bit-vector split rows for single-threaded unit-cost pairs.
7. **Ukkonen-banded Myers' bit-vector edit distance:** Under unit costs (`0 / -1 / 0 / -1`), `--scoreonly` computes 64 cells per CPU word inside an exponential Ukkonen diagonal band (matching `edlib` at 20 kb 1% in 3.4 ms).
8. **Smith-Waterman rectangles, strip checkpoints, & `--fastlocal`:** Large best-hit searches scan forward and reverse in $O(n)$ memory with lane-width hints to find the exact rectangle of the best hit; if that rectangle still exceeds $2^{28}$ cells, periodic row checkpoints refill one strip at a time. Opt-in `--fastlocal` reconstructs the sub-rectangle with BiWFA or Hirschberg and verifies exact local optimality.
9. **Column-block & batch threading:** Large single alignments split columns across threads in a pipelined wavefront with bounded 512-`pause` spins and per-boundary condition-variable waits; file batches align independent pairs in parallel and emit them in strict input order.

To reproduce the `5 kb – 50 kb` Nucleics Align timings (best of 3 with GNU `time` at 10 ms resolution, without `--wfa`), see [`benchmarks/README.md`](../benchmarks/README.md) and [`benchmarks/run.sh`](../benchmarks/run.sh). The `pop-os` `5 kb – 1 Mb` millisecond tables (including `--wfa`) were timed directly on simulated DNA pairs (`seed=11`, substitutions + indels) using Python's `time.perf_counter()` and GNU `time -v`.
