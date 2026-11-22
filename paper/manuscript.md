---
title: "Nucleics Align: fast, exact affine-gap pairwise sequence alignment for DNA and protein sequences up to megabase scale"
author: "Daniel Tillett^[Affiliation placeholder, City, Country; email@example.org; ORCID: 0000-0000-0000-0000]"
date: "22 Nov 2026"
bibliography: references.bib
---

## Abstract

**Background:** Exact pairwise sequence alignment with affine gap penalties (Needleman--Wunsch, Smith--Waterman, Gotoh) remains necessary whenever heuristic mappers may miss optimal alignments or when reproducible scores and fixed tie-breaking are required. Classic scalar tools such as Isaac Turner's `seq-align` are simple and exact, but aligning a 20 kb DNA pair takes ~11 s and 4.5 GiB of RAM, and a 50 kb pair would require ~28 GiB. Modern libraries (`parasail`, `ksw2`, `SSW`, `WFA2`, `SeqAn`, `EMBOSS`) each address parts of this space, but standalone C tools combining general substitution matrices, single-pair multi-threading, linear-space global alignment at high divergence, and budgeted wavefront alignment remain scarce.

**Results:** I present Nucleics Align (v1.1.0), a standalone C11 library and CLI suite derived from `seq-align`. It combines 1-byte-per-cell traceback storage, portable striped SIMD kernels (SSE2, SSE4.1, AVX2, AVX-512, Arm64 NEON) with prefix-scan carries, 8- and 16-bit difference-encoded lanes, formally certified diagonal banding with strip pruning and strip-kept memory, linear-space Hirschberg/Myers-Miller alignment with row reuse (`--linear`), budgeted wavefront alignment (`--wfa`), and Myers' bit-vector edit distance for unit-cost `--scoreonly`. Every default optimisation reproduces the exact optimal score and tie-breaking of this codebase's scalar Gotoh baseline (optional `--linear` and `--wfa` guarantee the identical optimal score but may select another equal-scoring alignment), verified across ~2.0 million unit-test assertion checks and 63 CLI fixtures (75 runs). On native Linux (8-core AMD Ryzen 7 5700G), Nucleics Align is 45–412x faster than original `seq-align` (1-seed baseline vs 5-seed median) on 20 kb global DNA on 1 thread (26.5–248.6 ms vs 10.9–12.1 s) and 147–400x faster on 8 threads (27.3–76.3 ms). On matched 1-thread 20 kb and 50 kb pairs, `--linear` is 2.5–14.9x faster than `ksw2` across all six conditions using 2.9–6.5 MiB (138–690x less RAM), and on 8 threads reduces wall time further while using less total CPU time than single-threaded `ksw2`. A nine-configuration ablation isolates each layer. On real RefSeq/GenBank pairs (8 threads), `--wfa` aligns two 197 kb Mpox genomes (99.92% identity) in 7.6 ms / 5.9 MiB, `--linear` aligns SARS-CoV-2 vs SARS-CoV-1 (29.9 kb, 79.54% identity) in 126.8 ms / 4.1 MiB, and `--linear` aligns two 580 kb *Mycoplasma genitalium* chromosomes (99.09% identity) in 2.28 s / 29.3 MiB. On simulated 1 Mb DNA, `--wfa` completes at 0.1% nominal mutation in 315.8 ms (312 MiB) and `--linear` at 5% in 24.3 s (64 MiB).

**Conclusions:** Nucleics Align makes exact affine-gap global alignment practical up to megabase scale and accelerates local and protein alignment on ordinary desktop hardware, without external alignment dependencies.

**Availability:** Source code, C library, CLI tools, Python 3 and Perl interfaces, and reproducible benchmark harnesses are released under the BSD 3-Clause licence at https://github.com/DanielTillett/nucleics-align.

---

## 1. Introduction

Exact pairwise sequence alignment by dynamic programming has been a core method in bioinformatics since the 1970s. Needleman and Wunsch introduced optimal global alignment [@needleman1970general], Smith and Waterman formulated optimal local alignment [@smith1981identification], and Gotoh reduced affine gap-penalty alignment from cubic to quadratic time using three coupled dynamic programming matrices [@gotoh1982improved]. Under an affine gap model, a contiguous insertion or deletion of $n$ residues scores $g_{\text{open}} + n \cdot g_{\text{extend}}$, penalising gap initiation more heavily than extension so that one contiguous indel is preferred over scattered single-base gaps.

In high-throughput genomics, heuristic seed-and-extend tools such as BLAST/BLAST+ [@altschul1990blast; @camacho2009blastplus] and minimap2 [@li2018minimap2] dominate database search and read mapping because $O(mn)$ dynamic programming is prohibitive across millions of reads or gigabase references. Exact pairwise alignment remains necessary, however, when comparing viral and organelle genomes, validating variants near complex indels, scoring homologous proteins with BLOSUM [@henikoff1992amino] or PAM [@dayhoff1978model] matrices, or establishing ground-truth reference alignments when testing heuristic tools.

Over the past five decades, four main families of exact speed-ups have been developed:

1. **SIMD vectorisation and multi-threading.** Early implementations exploited multimedia instructions along anti-diagonals [@wozniak1997using] or query profiles [@rognes2000sixfold]. Farrar's striped SIMD layout [@farrar2007striped] moved the vertical-gap dependency out of the inner loop into a lazy-F correction pass and underpins `SSW` [@zhao2013ssw; @sswrepo] and `parasail` [@daily2016parasail]. Prefix-scan carry propagation was introduced for GPU Smith-Waterman by Khajeh-Saeed et al. [@khajehsaeed2010acceleration], implemented in `parasail`'s scan kernels [@daily2016parasail], and used inside adaptive SIMD blocks by Block Aligner [@liu2023block]. Suzuki and Kasahara [@suzuki2018introducing] showed that storing vertical and horizontal score differences keeps lane values within 8 or 16 bits regardless of sequence length, an approach used in `ksw2` [@li2018minimap2; @ksw2repo]. The C++ `SeqAn` library [@rahn2018seqan] combines vectorised DP with multi-threading across and within sequence pairs.
2. **Linear-space divide-and-conquer.** Hirschberg [@hirschberg1975linear] and Myers and Miller [@myers1988optimal] reduced the working memory of exact global alignment from $O(mn)$ to linear in sequence length at the cost of roughly twice the arithmetic work, later combined with diagonal banding by Chao et al. [@chao1992aligning]. Many biologists encounter this algorithm via EMBOSS `stretcher` (alongside full-matrix `needle`) [@rice2000emboss].
3. **Output-sensitive and volume-pruning methods.** Ukkonen [@ukkonen1985algorithms] showed that edit distance can be computed in $O(s \cdot \min(m, n))$ time for edit distance $s$, and Spouge [@spouge1989speeding] formalised computational volumes for pruning unreachable lattice states. For affine gaps, the Wavefront Alignment (WFA) algorithm [@marcosola2021fast; @wfa2lib] runs in $O(ns)$ worst-case time and $O(s^2)$ space; Eizenga and Paten (preprint) [@eizenga2022improving] generalised its scoring and analysed an $O(n + s^2)$ expected-time variant under uniform random errors; bidirectional WFA (`BiWFA`) reduces space to $O(s)$ [@marcosola2023optimal]; and `Singletrack` [@lopezvillellas2026singletrack] cuts affine backtrace memory 3x inside `WFA2-lib`. Groot Koerkamp's `A*PA2` [@grootkoerkamp2024apa2] applies admissible heuristics and strip-based computational volumes to unit-cost edit distance.
4. **Bit-vector parallelism.** For unit-cost Levenshtein edit distance, Myers' bit-vector algorithm [@myers1999fast] computes 64 DP cells per machine word using a small constant number of bitwise operations and underlies `edlib` [@sosic2017edlib].

Each existing tool occupies a distinct niche. `edlib` [@sosic2017edlib] is limited to unit-cost edit distance without affine gaps or substitution matrices. `WFA2` and `BiWFA` [@marcosola2021fast; @marcosola2023optimal] are fast on closely related pairs, but their work grows with penalty score $s$ so they slow sharply on divergent or unrelated 20–50 kb pairs, and they do not support arbitrary protein substitution matrices. `parasail` [@daily2016parasail] and `ksw2` [@li2018minimap2] support general scoring well, but run single-threaded per pair and consume gigabytes of RAM on 50 kb global alignments. `SSW` [@zhao2013ssw] is fast for local alignment but does not provide global affine-gap alignment. `EMBOSS` `needle` and `stretcher` [@rice2000emboss] are scalar implementations, while `SeqAn` [@rahn2018seqan] is a C++ template library rather than a lightweight standalone C tool.

I started Nucleics Align by forking Isaac Turner's `seq-align` [@turner2017seqalign] for its exact Gotoh formulation, fixed tie-breaking, and minimal C dependency footprint. Its limitation was resource use: on a modern desktop, aligning a single 20 kb DNA pair took ~11 s and 4.5 GiB of RAM, and a 50 kb pair would have required ~28 GiB. Nucleics Align keeps `seq-align`'s exact mathematics and scoring flexibility while integrating all four families of exact speed-ups into one portable C11 library and CLI suite. The engineering rule throughout has been strict: every default optimisation must produce the exact same score and alignment, byte for byte, as this codebase's scalar Gotoh baseline (the few intentional bug fixes relative to the original 2017 `seq-align` are listed in `CHANGELOG.md`).

---

## 2. Implementation

Nucleics Align is written in portable C11 with no external alignment library dependencies (bundling only Isaac Turner's `string_buffer`, `seq_file`, and `sort_r`, plus `zlib` and POSIX threads). It builds on Linux, macOS, and Windows under WSL (native Windows compilers are not supported), and provides three CLI tools (`needleman_wunsch`, `smith_waterman`, `lcs`), a static C library (`libalign.a`), a zero-dependency Python 3 wrapper, and Perl modules.

### 2.1 Execution dispatch and modes

Global alignment in `needleman_wunsch` dispatches according to the requested mode and scoring options:

1. **Score-only (`--scoreonly`):** Invokes the checked linear-space score pass. Where scoring options are supported, it first tries a banded striped SIMD pass; if banding does not pay, it attempts Myers' bit-vector algorithm [@myers1999fast] under unit-cost scoring (`0 / -1 / 0 / -1`) on one thread (or when SIMD is unavailable), falling back to a full striped or scalar linear pass. On unsupported scoring (such as free terminal gaps or positive gap-open penalties), it falls back to the full global aligner; other errors terminate with an error status.
2. **Explicit linear space (`--linear`):** Calls the checked Hirschberg/Myers-Miller aligner directly and exits non-zero on unsupported options or errors.
3. **Default global alignment:** When `--wfa` is enabled, it first attempts budgeted wavefront alignment; if `--linearabove <n>` is set and the longer sequence has at least $n$ letters, it attempts linear-space alignment. If either is unsupported or does not complete, execution falls back to the matrix aligner. The matrix aligner uses a compact 1-byte-per-cell traceback fill whenever the scoring options are supported by the striped kernels: it first attempts a short budgeted score-only WFA pass (default budget 16 diagonal updates per letter of `seq_a`) or 12-mer estimate to size an exact band (Section 2.4), falling back to a full difference-encoded matrix fill if banding does not pay, and to the scalar three-matrix fill if scoring options are unsupported by the striped kernels. If `--linearfallback` is enabled and matrix allocation fails or overflows `size_t`, it falls back to linear space.

### 2.2 Compact traceback storage and deterministic tie-breaking

The baseline Gotoh implementation stores three 32-bit score matrices (`MATCH`, `GAP_A`, `GAP_B`), requiring 12 bytes per cell (~4.47 GiB for $20{,}000 \times 20{,}000$ cells). In Nucleics Align, compact mode (default in the CLI) packs the three 2-bit predecessor state indicators into one byte per cell while keeping only two score rows in RAM during the forward fill, cutting matrix memory 12-fold.

Whenever multiple paths tie for the optimal score, traceback in every full-matrix path (scalar matrices, compact traceback bytes, exact bands, strip-kept bands, and difference-encoded lanes) follows the same deterministic rule. In Gotoh traceback, the current matrix state fixes the spatial step (`GAP_A` steps to $(x, y-1)$, `GAP_B` to $(x-1, y)$, `MATCH` to $(x-1, y-1)$), and at that predecessor cell the next state is chosen in strict priority order `GAP_A`, `GAP_B`, `MATCH`. Global alignment starts at $(|A|, |B|)$ using the same state priority; Smith-Waterman always starts in `MATCH` at the highest-scoring cell (breaking ties by lowest end coordinate $x$, then $y$). This priority tends to push gaps toward the right end of runs (e.g., `AAA` vs `AA` gives `AA-`, and `ACCA` vs `ACA` gives `AC-A`). By contrast, `--linear` and `--wfa` always return the identical optimal score but may pick a different equal-scoring alignment due to midpoint or wavefront tie-breaking.

### 2.3 Portable striped SIMD, prefix-scan carries, and narrow/difference-encoded lanes

All vectorised kernels use portable GCC/clang vector extensions rather than architecture-specific intrinsics. On x86-64 Linux, GCC/clang function multi-versioning (`target_clones` for `avx2`, `sse4.1`, and `default` SSE2) compiles all three variants into one portable binary (8-lane 32-bit vectors under GCC; 4-lane under Clang dispatch, as Clang rejects 8-lane vectors inside `target_clones` multi-versioned functions). Building with `NATIVE=1` (`-march=native`) enables 16 × 32-bit AVX-512 lanes where present, and Arm64 builds (including Apple silicon) use NEON automatically. On non-Linux x86-64 builds without `NATIVE=1`, multi-versioning is disabled and the build uses baseline SSE2 (with linear-space passes falling back to scalar code).

The striped fill follows Farrar's layout [@farrar2007striped], precomputing a per-query substitution profile for the active alphabet. Drawing on prefix-scan carry propagation [@khajehsaeed2010acceleration; @daily2016parasail; @liu2023block], Nucleics Align precomputes inter-lane vertical-gap carries in a single scan rather than iterating multi-sweep lazy-F loops across banded rows. To increase lane density without losing exactness:
- **Global full traceback fills** store vertical and horizontal score differences [@suzuki2018introducing], keeping values in 8-bit lanes when the per-step score range fits (including standard DNA scoring) and 16-bit lanes for wider matrices such as BLOSUM62 at any sequence length.
- **Local fills** begin in 8-bit lanes (four times as many lanes as 32-bit vectors) and widen to 16- or 32-bit lanes only if a score threatens to saturate.
- **Linear-space passes** use 16-bit relative score blocks.

### 2.4 Exact certified banding, strip pruning, and strip-kept memory

Let $m = |A|$, $n = |B|$, $D = m - n$, $s_{\max} = \max\big(0, \max_{a,b} S(a,b)\big)$, and assume non-positive gap penalties $g_{\text{open}} \le 0$ and $g_{\text{extend}} \le 0$ with $s_{\max} - 2 g_{\text{extend}} > 0$. Before filling a banded global matrix, Nucleics Align first attempts a budgeted score-only WFA pass [@marcosola2021fast] (default 16 diagonal updates per letter of $A$); when it finishes, it yields the exact optimal score. Otherwise, a fast 12-mer sampling routine produces a deliberately high heuristic score estimate used solely to decide whether a band is worth trying; the final alignment never depends on that estimate.

Exactness is certified mathematically from the score achieved inside the closed band $B(b) = \{(x, y) : \min(0, D) - b \le x - y \le \max(0, D) + b\}$. Let $S_{\text{out}}(b)$ denote the maximum score of any global path that leaves $B(b)$ (or $-\infty$ if no such path exists). Any global path visiting a cell outside $B(b)$ requires at least $I = |D| + 2(b + 1)$ gap characters and at most $(m + n - I)/2 = \min(m, n) - (b + 1)$ aligned residue pairs; because $g_{\text{open}} \le 0$ and $s_{\max} \ge 0$, every path leaving $B(b)$ satisfies:

$$S_{\text{out}}(b) \le s_{\max}\big(\min(m, n) - (b + 1)\big) + g_{\text{extend}}\big(|D| + 2(b + 1)\big)$$

Consequently, once a banded fill finds a feasible global alignment of score $L$, every path leaving $B(b)$ satisfies the strict inequality $S_{\text{out}}(b) < L$ (excluding both better and tied optimal paths outside $B(b)$) whenever:

$$b \ge \left\lfloor \frac{s_{\max}\min(m, n) + g_{\text{extend}}|D| - L}{s_{\max} - 2 g_{\text{extend}}} \right\rfloor$$

Excluded predecessor cells outside $B(b)$ are set to the DP floor so no cell score inside $B(b)$ is ever overestimated. If this certificate fails, the band is widened and re-run until the inequality holds, or abandoned for a full fill (whose exactness does not depend on the band certificate), certifying that every optimal and tied-optimal path stays within the closed band $B(b)$ and preserving byte-identical full-matrix traceback.

Inspired by `A*PA2` [@grootkoerkamp2024apa2] and Spouge's computational volumes [@spouge1989speeding], the band is filled in horizontal strips. Once a feasible lower bound $L$ is known (from the exact WFA score or from a prior banded pass) and $g_{\text{extend}} < 0$, for each cell $(x, y)$ on the bottom row of a strip with best state score $v = \max(M[x,y], G_A[x,y], G_B[x,y])$ and remaining lengths $r_i = m - x$, $r_j = n - y$, any completion to $(m, n)$ scores at most:

$$U(x, y) = v + s_{\max}\min(r_i, r_j) + g_{\text{extend}}|r_i - r_j|$$

Because $s_{\max} \ge 0$, $U(x, y)$ remains a valid upper bound even when all substitution scores are negative. Cells with $U(x, y) < L$ cannot lie on any optimal or tied-optimal path and are pruned; each surviving cell on diagonal $d = x - y$ can reach at most $k = \lfloor (U(x, y) - L) / (-2 g_{\text{extend}}) \rfloor$ diagonals beyond $[\min(d, D), \max(d, D)]$. When the estimated band footprint exceeds one quarter of physical RAM (falling back to 1 GiB when physical memory cannot be queried), the aligner switches automatically to strip-kept mode, retaining only strip boundary rows during the forward pass and refilling one strip's traceback bytes at a time during traceback. For a 1 Mb, 1% nominal mutation global alignment on a 62 GiB host, this avoids an estimated ~35 GiB band allocation and runs in 319 MiB while preserving full-matrix tie-breaking.

### 2.5 Budgeted wavefront alignment (`--wfa`)

When `--wfa` is requested, Nucleics Align converts affine match/mismatch/gap scoring into non-negative penalty costs [@eizenga2022improving] and advances wavefronts by penalty score [@marcosola2021fast], skipping matching runs 8 bytes at a time via 64-bit word comparisons and tracing back through retained wavefronts. Unlike `BiWFA`'s $O(s)$ space [@marcosola2023optimal], retained-wavefront traceback uses $O(s^2)$ space and $O(ns)$ worst-case time (with expected $O(n + s^2)$ scaling under uniform random errors [@eizenga2022improving]). The global aligner therefore caps WFA at $\min(2^{25},\; 64 \cdot (|A| + |B|))$ wavefront cells; if the scoring does not reduce to one match score, one mismatch score, and standard affine gap penalties, or if the budget is exhausted, execution falls back to the matrix aligner.

### 2.6 Linear-space global alignment (`--linear`) and Myers bit-vector score-only

The `--linear` aligner implements Hirschberg's divide-and-conquer algorithm [@hirschberg1975linear] extended to Gotoh's affine gap model [@myers1988optimal], combined with certified banded striped passes [@chao1992aligning; @farrar2007striped], 16-bit relative blocks, and two-level intermediate row reuse so child sub-problems skip redundant passes. Working score-row memory is $O(|B|)$ (with constant-factor overhead for row reuse and per-thread buffers), plus $O(|A| + |B|)$ for the returned alignment strings.

For `--scoreonly`, when scoring is unit-cost edit distance (`0 / -1 / 0 / -1`), banding does not pay, and execution is single-threaded (or SIMD is unavailable), Nucleics Align dispatches to Myers' bit-vector algorithm [@myers1999fast], computing 64 cells per machine word using a small constant number of bitwise operations.

### 2.7 Best-hit local alignment (`smith_waterman --maxhits 1`)

When only the single highest-scoring local hit is requested (`--maxhits 1`), Nucleics Align scans forward in $O(|B|)$ memory to locate the optimal end cell and score, then scans in reverse from that cell in $O(|B|)$ memory to locate the start cell. Only the resulting sub-rectangle is filled for traceback, with periodic row checkpoints if the rectangle exceeds $2^{28}$ cells.

### 2.8 Multi-threading, error handling, and verification

By default, the CLI tools use half of the host's logical CPUs (`--threads <n>` overrides this). Large single alignments (exceeding $4 \times 10^6$ cells and 1,024 columns per thread block) split columns into pipelined wavefront blocks across threads, while file batches (`--file`, `--files`) align independent pairs concurrently and emit results in strict input order so output is byte-identical to `--threads 1`.

In the C library, every major aligner provides a checked entry point returning an error status (or `1`/`0` with an `aln_error_t` parameter for `smith_waterman_fetch_checked()`) instead of calling `exit()`. Continuous integration runs on Linux (x86-64 and Arm64) and macOS with GCC and Clang, including AddressSanitizer and UndefinedBehaviorSanitizer, ~2.0 million unit-test assertion checks (many on randomised inputs comparing every optimised path against the scalar Gotoh baseline), and 63 CLI regression fixtures (75 test runs, as 12 batch fixtures run at both 1 and 4 threads).

---

## 3. Results

All main benchmarks were run on a native Linux desktop (AMD Ryzen 7 5700G, 8 cores / 16 logical CPUs, 62 GiB RAM, Pop!_OS 24.04 Linux, GCC 13.3.0, default `make` `-O3` build). CLI timings include process startup, input parsing, alignment, and output writing (so sub-10 ms runs include ~2–4 ms of process startup overhead); the C library ablation (Table 2) times the alignment call directly via `clock_gettime(CLOCK_MONOTONIC)`. Memory figures report peak resident set size (RSS) from GNU `time -v` in binary MiB/GiB (`1 MiB = 1024^2` bytes; `1 GiB = 1024 MiB`), and CPU times come from GNU `time -v` (10 ms resolution, where `0.00 s` is displayed as `<10 ms`). Unless noted, global DNA runs used `+1` match, `-2` mismatch, `-4` gap open, `-1` gap extend; local DNA runs used `+2 / -2 / -2 / -1`; and protein runs used `--scoring BLOSUM62` (`-10 / -1`). All global modes and comparators in Tables 1, 2, and 4 were programmatically asserted to return the identical optimal score on every completed pair.

### 3.1 Real biological sequences

Seven public NCBI RefSeq/GenBank pairs spanning mitochondrial genomes, coronaviruses, orthopoxviruses, whole bacterial chromosomes, and mammalian proteins were timed with 8 default threads (best of 3 runs for pairs $\le 200$ kb, 1 run for *M. genitalium*; reproducible via `benchmarks/real_sequences.sh`; Table 1). Circular mitochondrial and bacterial genomes were aligned as deposited in RefSeq without origin rotation; the 580 kb *M. genitalium* pair serves as a whole-chromosome scaling stress test rather than a rearrangement-aware comparison. Reported sequence identities belong to the `full` traceback and count case-insensitive equal symbols divided by total alignment columns (so `N/N` counts as an identity match under `--pretty`, whereas under default DNA alignment scoring ambiguous bases such as `N` are scored as mismatches).

**Table 1. Performance on real biological sequence pairs (`pop-os`, AMD Ryzen 7 5700G, 8 default threads, v1.1.0; wall ms / peak MiB).** DNA global modes used `+1/-2/-4/-1`; DNA `local --maxhits 1` used `+2/-2/-2/-1`; BRCA1 used `BLOSUM62` (`-10/-1`). Bold marks the fastest global mode returning a full alignment. `†` = `--wfa` exceeded its cell budget and fell back to matrices; `‡` = `--wfa` unsupported for substitution matrix and fell back; `—` = not run.

| Pair (NCBI Accessions) | Lengths | Identity | `full` (ms / MiB) | `--wfa` (ms / MiB) | `--linear` (ms / MiB) | `--scoreonly` (ms / MiB) | `local --maxhits 1` (ms / MiB) | Optimal global score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Human/Neanderthal mtDNA (`NC_012920.1`, `NC_011137.1`) | 16.6 / 16.6 kb | 98.72% | 18.2 / 11.5 | **8.5** / 8.6 | 23.6 / 3.0 | 6.8 / 3.1 | 118.3 / 20.0 | 15,928 |
| Human/Chimp mtDNA (`NC_012920.1`, `NC_001643.1`) | 16.6 / 16.6 kb | 85.19% | **52.7** / 266.7 | 75.5† / 267.2 | 53.2 / 4.2 | 20.7 / 3.7 | 113.4 / 19.7 | 10,632 |
| SARS-CoV-2 Wuhan/Omicron (`NC_045512.2`, `OL672836.1`) | 29.9 / 29.7 kb | 99.03% | 23.1 / 16.1 | **7.7** / 7.2 | 47.0 / 3.6 | 13.2 / 3.4 | 259.6 / 42.0 | 29,240 |
| SARS-CoV-2/SARS-CoV-1 (`NC_045512.2`, `NC_004718.3`) | 29.9 / 29.8 kb | 79.54% | 168.9 / 856.6 | 210.5† / 857.0 | **126.8** / 4.1 | 32.1 / 4.3 | 222.4 / 41.9 | 11,467 |
| Mpox Ref/2022 isolate (`NC_063383.1`, `ON563414.3`) | 197.2 / 197.2 kb | 99.92% | 91.4 / 63.7 | **7.6** / 5.9 | 286.6 / 11.1 | 50.0 / 10.4 | 9,591.1 / 625.9 | 196,823 |
| *M. genitalium* G37/M2288 (`NC_000908.2`, `CP003773.1`) | 580.1 / 579.6 kb | 99.09% | **2,100.2** / 4,430.7 | 2,429.6† / 4,430.9 | 2,281.4 / 29.3 | 577.6 / 27.4 | — | 564,775 |
| Human/Mouse BRCA1, BLOSUM62 (`NP_009225.1`, `NP_033894.3`) | 1,863 / 1,812 aa | 55.76% | **6.3** / 6.4 | 6.6‡ / 6.5 | 6.4 / 2.5 | 4.1 / 2.4 | 7.5 / 7.2 | 4,827 |

On closely related viral genomes, `--wfa` aligns two 197 kb Mpox genomes (99.92% identity) in 7.6 ms and 5.9 MiB (12.0x faster and 10.8x smaller than `full`, and 6.6x faster than `--scoreonly`), and aligns SARS-CoV-2 Wuhan vs Omicron (99.03% identity) in 7.7 ms (3.0x faster than `full`). On divergent SARS-CoV-2 vs SARS-CoV-1 (79.54% identity), `--linear` is faster than `full` (126.8 ms vs 168.9 ms) while using 209x less RAM (4.1 MiB vs 856.6 MiB). On Human vs Chimp mtDNA (85.19% identity), `--linear` matches `full` wall time (53.2 ms vs 52.7 ms) while cutting RAM 63.5-fold (4.2 MiB vs 266.7 MiB). Two 580 kb *M. genitalium* chromosomes (99.09% identity, 575,139 matches, 4,101 mismatches, 1,154 gap columns) align in 2.10 s with `full` (4.33 GiB) or 2.28 s with `--linear` using 29.3 MiB (151x less RAM). On Human vs Mouse BRCA1 under BLOSUM62 (`-10 / -1`), `full` and `--linear` finish in 6.3–6.4 ms, and `--wfa` falls back to `full`.

### 3.2 Nine-configuration ablation

To quantify each configuration independently, `benchmarks/ablation.c` timed nine configurations of the C11 library API (rows 1–6 add successive optimisations to the matrix kernel; rows 7–9 are alternative linear-space and wavefront modes) across five independently generated 20 kb DNA pairs (`benchmarks/gen_pairs.py`, seeds 1–5, substitution probability $p \in \{1\%, 10\%\}$ with insertion and deletion each at probability $p/5$; Table 2). Every configuration returned the exact same optimal score on every pair.

**Table 2. Nine-configuration ablation on 20 kb DNA (`pop-os`, AMD Ryzen 7 5700G, 5 independent seeds, C11 library API, median [min–max] wall ms).** `†` = `--wfa` hit its cell budget and fell back.

| Configuration | 20 kb, 1% sub. (ms) | 20 kb, 10% sub. (ms) |
|---|---:|---:|
| 1. Scalar 3-matrix baseline (1 thr.) | 2,743.2 [2,737.8–2,756.3] | 2,858.1 [2,857.4–2,868.9] |
| 2. + Compact 1-byte traceback (scalar, 1 thr.) | 2,752.1 [2,738.0–2,756.8] | 2,847.1 [2,843.0–2,851.1] |
| 3. + Striped 32-bit SIMD (1 thr.) | 440.3 [438.8–443.9] | 428.0 [423.3–430.6] |
| 4. + Difference-encoded 8-bit lanes (1 thr.) | 292.9 [290.3–295.3] | 289.9 [288.8–292.8] |
| 5. + Certified exact banding (1 thr.) | 23.1 [21.9–24.5] | 297.7 [297.4–298.6] |
| 6. + Multi-threading (8 thr.) | 23.0 [21.3–23.4] | 88.3 [82.1–91.9] |
| 7. Linear space `--linear` (1 thr.) | 38.7 [37.6–38.9] | 165.2 [163.1–179.7] |
| 8. Linear space `--linear` (8 thr.) | 29.2 [28.5–29.6] | 77.1 [74.9–83.6] |
| 9. Wavefront `--wfa` (1 thr.) | 16.0 [11.3–16.8] | 382.7† [375.3–403.4] |

Nucleics Align's scalar 3-matrix baseline (2.74 s) is already ~4x faster than original `seq-align` (10.9 s in Table 4) because it replaces per-cell `scoring_lookup()` calls with a precomputed substitution table (`sub_table_t`) and avoids cold 4.5 GiB `calloc` overhead. Compact 1-byte traceback cuts matrix memory 12x (4.47 GiB to 386 MiB) with equal scalar speed. Striped 32-bit SIMD yields a 6.2–6.7x speedup, and 8-bit difference encoding adds a further 1.5x (290–293 ms). At 1% substitution, certified exact banding gives a 12.7x jump (292.9 ms to 23.1 ms) and `--wfa` reaches 16.0 ms. At 10% substitution, enabling banding adds ~8 ms of overhead over the unbanded 8-bit fill (297.7 vs 289.9 ms), 8-thread parallelism adds 3.4x (88.3 ms), and `--linear` reaches 165.2 ms on 1 thread and 77.1 ms on 8 threads. (Note that Table 2's 10% times slightly exceed Table 4's CLI times because `ablation.c` allocates and frees each aligner inside the timed call immediately after the 4.5 GiB scalar run without warm-up.)

### 3.3 Scaling from 5 kb to 1 Mb on simulated DNA

Simulated DNA pairs from 5 kb to 1 Mb were generated with `benchmarks/make_pairs.sh` (`seed=11`; each mutated position draws uniformly from all four bases so the effective substitution rate is ~0.75x the nominal mutation setting, with indels at one fifth of the nominal setting; `awk`'s `srand`/`rand` differs between `mawk` and `gawk`) and timed on `pop-os` (8 default threads; Table 3).

**Table 3. Scaling on simulated DNA from 5 kb to 1 Mb (`pop-os`, AMD Ryzen 7 5700G, 8 default threads; wall ms / peak MiB).** Standard lengths (5–50 kb): best of 3 runs; large lengths (100 kb – 1 Mb): 1 run. DNA global modes used `+1/-2/-4/-1`; `local --maxhits 1` used `+2/-2/-2/-1`. Bold marks the fastest global mode returning a full alignment. `†` = `--wfa` hit its cell budget and fell back; `—` = not run.

| Length | Nominal mut. | `full` (ms / MiB) | `--wfa` (ms / MiB) | `--linear` (ms / MiB) | `--scoreonly` (ms / MiB) | `local --maxhits 1` (ms / MiB) |
|---|---|---:|---:|---:|---:|---:|
| 5 kb | 1% | 3.4 / 3 | **1.8** / 3 | 6.5 / 3 | 2.6 / 3 | 21.1 / 31 |
| 5 kb | 10% | 11.7 / 9 | 17.4† / 10 | **9.7** / 3 | 4.0 / 3 | 19.9 / 31 |
| 5 kb | unrelated | 13.3 / 27 | 18.8† / 27 | **12.7** / 3 | 4.4 / 2 | 7.0 / 3 |
| 20 kb | 1% | 22.7 / 16 | **13.1** / 14 | 29.0 / 3 | 8.1 / 3 | 149.7 / 25 |
| 20 kb | 10% | 73.8 / 388 | 102.0† / 388 | **62.6** / 4 | 14.5 / 4 | 135.0 / 24 |
| 20 kb | unrelated | **71.3** / 387 | 100.2† / 388 | 73.4 / 4 | 20.2 / 4 | 23.9 / 4 |
| 50 kb | 1% | 101.2 / 71 | 166.4† / 76 | **97.0** / 4 | 28.0 / 4 | 577.4 / 87 |
| 50 kb | 10% | **165.7** / 584 | 244.1† / 584 | 216.8 / 6 | 60.7 / 6 | 546.4 / 87 |
| 50 kb | unrelated | 417.1 / 2,395 | 489.2† / 2,396 | **276.4** / 8 | 97.6 / 6 | 129.7 / 7 |
| 100 kb | 0.1% | 44.1 / 34 | **6.1** / 6 | 143.2 / 7 | 26.6 / 6 | — |
| 100 kb | 1% | 350.6 / 268 | 488.4† / 268 | **274.7** / 7 | 94.6 / 6 | — |
| 100 kb | 5% | **409.3** / 1,234 | 567.3† / 1,234 | 515.4 / 8 | 122.6 / 7 | — |
| 300 kb | 0.1% | 337.9 / 256 | **31.4** / 32 | 523.9 / 15 | 117.5 / 15 | — |
| 300 kb | 1% | **876.1** / 2,376 | 1,246.2† / 2,376 | 1,054.1 / 16 | 246.8 / 16 | — |
| 300 kb | 5% | 2,731.3 / 11,625 | 3,185.1† / 11,625 | **2,522.6** / 20 | 944.7 / 19 | — |
| 1 Mb | 0.1% | 3,480.8 / 2,641 | **315.8** / 312 | 2,910.7 / 46 | 935.2 / 43 | — |
| 1 Mb | 1% | 8,103.1 / 319 *(strips)* | 8,517.2† / 392 | **6,402.6** / 49 | 2,256.7 / 46 | — |
| 1 Mb | 5% | — | — | **24,257.5** / 64 | 11,955.8 / 61 | — |

At 0.1% nominal mutation, `--wfa` aligns 1 Mb in 315.8 ms and 312 MiB (11.0x faster than `full`, 9.2x faster than `--linear`, and 3.0x faster than `--scoreonly`). At 1%, `--wfa` wins at 5 kb and 20 kb, but from 50 kb upward it exhausts its cell budget and falls back to the matrix aligner. Across longer divergent pairs, `--linear` keeps peak RSS at 46–64 MiB even at 1 Mb and 5%. In a separate earlier strip-memory run on the same host (21 Nov 2026), `full` strip-kept bands completed 1 Mb pairs in 5.36 s / 287 MiB (0.5%), 7.92 s / 319 MiB (1%), 15.8 s / 467 MiB (2%), and 43.3 s / 1.1 GiB (5%) while preserving full-matrix tie-breaking.

### 3.4 Reproducible 5-seed comparison with existing exact aligners

Using the self-contained harness `benchmarks/compare_harness.sh` on native Linux `pop-os`, Nucleics Align (`d434eb2`, default `make` `-O3`) was compared against pinned builds of original `seq-align` (`dc41988`, default `make`) [@turner2017seqalign], `ksw2` (`289609b`, `gcc -O3 -msse4.1`, `ksw_extz2_sse` with `w=-1`, `zdrop=-1`, `flag=0`, returning score + CIGAR) [@li2018minimap2; @ksw2repo], `WFA2-lib` (`v2.3.5` / `cf3eb92`, `gcc -O3`, exact high-memory and ultralow `BiWFA` modes with `wf_heuristic_none`, `compute_alignment` CIGAR, affine penalties `match=0, mismatch=6, gap_open=8, gap_extend=3` mapped back via `score = (|A|+|B|-penalty)/2`) [@marcosola2021fast; @marcosola2023optimal; @wfa2lib], and `SSW` (`a66636b`, `gcc -O3`, `ssw_align` with `+2/-2/-2/-1` local scoring, score + CIGAR) [@zhao2013ssw; @sswrepo]. All tools were timed across five independently generated deterministic DNA pairs per condition (`benchmarks/gen_pairs.py`, seeds 1–5; 1 seed for the >9 s `seq-align` and `WFA2 high` unrelated baselines; `WFA2` on 50 kb unrelated was skipped in the harness). Every global tool returned the exact same optimal score on every completed pair (Table 4).

**Table 4. Reproducible 5-seed comparator benchmark (`pop-os`, AMD Ryzen 7 5700G, native Linux; median wall ms [min–max] and median CPU ms in parentheses).** Bold marks the lowest wall time per column within each block. `†` = `--wfa` hit its cell budget and fell back.

| Benchmark & tool | 1% sub. (ms) | 10% sub. (ms) | Unrelated (ms) | Median RSS (1% / 10% / rand) |
|---|---|---|---|---|
| **Global, 20 kb** | | | | |
| `seq-align` (`dc41988`, 1 thr.) | 10,912 ms | 11,239 ms | 12,147 ms | 4,577 / 4,586 / 4,580 MiB |
| Nucleics Align (default, 1 thr.) | 26.5 [25.3–27.3] (20) | 248.6 [246.1–250.7] (240) | 222.3 [219.0–224.9] (210) | 17.5 / 385.8 / 385.8 MiB |
| Nucleics Align `--linear` (1 thr.) | 42.1 [41.8–43.8] (30) | 147.0 [144.0–151.1] (140) | 170.1 [168.0–172.4] (160) | 2.9 / 3.2 / 3.8 MiB |
| Nucleics Align `--wfa` (1 thr.) | 19.3 [17.8–19.4] (10) | 273.6† [271.9–275.7] (260) | 247.8† [247.1–252.4] (240) | 19.9 / 386.1 / 385.8 MiB |
| Nucleics Align (default, 8 thr.) | 27.3 [24.9–27.8] (20) | 76.3 [72.6–78.6] (540) | **76.2** [73.8–83.0] (550) | 17.5 / 387.2 / 387.2 MiB |
| Nucleics Align `--linear` (8 thr.) | 30.8 [30.5–32.4] (40) | **73.0** [71.8–73.6] (160) | 78.0 [77.5–81.1] (220) | 3.2 / 3.8 / 4.3 MiB |
| `ksw2` (`289609b`, 1 thr.) | 424.4 [419.0–434.7] (410) | 427.7 [423.5–441.6] (420) | 433.7 [429.9–446.0] (430) | 523.4 / 523.4 / 523.4 MiB |
| `WFA2` (`v2.3.5` high, 1 thr.) | 17.8 [16.4–19.1] (10) | 888.8 [855.5–968.3] (870) | 9,491 (9,480) | 29.2 / 1,491 / 17,029 MiB |
| `WFA2 BiWFA` (`v2.3.5`, 1 thr.) | **6.3** [6.2–7.2] (<10) | 156.1 [152.1–170.8] (150) | 1,653.2 [1,631–1,682] (1,640) | 3.9 / 9.7 / 24.0 MiB |
| **Global, 50 kb** | | | | |
| Nucleics Align (default, 1 thr.) | 120.3 [111.9–123.2] (110) | 1,491.2 [1,488–1,498] (1,480) | 1,338.6 [1,337–1,353] (1,330) | 83.4 / 2,390 / 2,392 MiB |
| Nucleics Align `--linear` (1 thr.) | 159.2 [154.5–164.3] (150) | 833.8 [829.6–852.5] (820) | 948.7 [944.5–1,016] (940) | 4.0 / 4.9 / 6.5 MiB |
| Nucleics Align (default, 8 thr.) | 118.8 [111.7–121.4] (110) | 450.3 [439.4–460.0] (3,330) | 430.3 [420.9–437.5] (3,300) | 83.4 / 2,393 / 2,395 MiB |
| Nucleics Align `--linear` (8 thr.) | 106.0 [104.9–110.5] (160) | **229.0** [219.9–241.5] (860) | **279.1** [277.0–291.8] (1,170) | 4.4 / 6.3 / 7.7 MiB |
| `ksw2` (`289609b`, 1 thr.) | 2,376.8 [2,366–2,389] (2,360) | 2,365.0 [2,350–2,389] (2,360) | 2,358.9 [2,351–2,367] (2,350) | 2,762 / 2,759 / 2,761 MiB |
| `WFA2` (`v2.3.5` high, 1 thr.) | 71.9 [65.5–81.0] (60) | 5,087.3 [5,044–5,337] (5,080) | *(not run)* | 120.9 / 9,226 / — MiB |
| `WFA2 BiWFA` (`v2.3.5`, 1 thr.) | **18.6** [17.4–19.5] (10) | 866.5 [855.0–899.9] (860) | *(not run)* | 5.5 / 18.5 / — MiB |
| **Best local hit, 20 kb (`--maxhits 1`)** | | | | |
| Nucleics Align (1 thr.) | 507.6 [499.8–511.1] (490) | 409.4 [405.7–411.5] (390) | 153.7 [152.5–156.6] (150) | 22.1 / 21.9 / 2.9 MiB |
| Nucleics Align (8 thr.) | **147.4** [143.0–155.0] (680) | 137.9 [134.2–144.2] (590) | **24.3** [24.1–24.9] (170) | 24.7 / 24.6 / 4.4 MiB |
| `SSW` (`a66636b`, 1 thr.) | 2,147.8* [196–2,220] (2,140) | **128.5** [126.0–154.3] (120) | 519.3 [271–796] (510) | 884.5 / 7.7 / 95.4 MiB |

In matched single-thread runs, Nucleics Align `--linear` is 2.5–14.9x faster than `ksw2` across all six 20 kb and 50 kb conditions (using 138–690x less RAM), beats `WFA2 BiWFA` at 10% substitution (147.0 vs 156.1 ms at 20 kb; 833.8 vs 866.5 ms at 50 kb) and by 9.7x on 20 kb unrelated pairs (170.1 vs 1,653.2 ms), while `BiWFA` remains fastest at 1% substitution (6.3 ms at 20 kb; 18.6 ms at 50 kb). On 8 threads, Nucleics Align `--linear` drops 50 kb 10% and unrelated pairs to 229.0 ms and 279.1 ms (6.3–7.7 MiB). On 20 kb local alignment (`--maxhits 1`), `SSW` is faster than 1-thread Nucleics Align at 10% substitution (128.5 vs 409.4 ms), but at 1% substitution (`*`) the optimal local score (~38,000–39,000) exceeds SSW's signed 16-bit limit (`INT16_MAX = 32,767`; SSW has no 32-bit path, so scores saturate and exact traceback is not guaranteed), and SSW is slower on unrelated pairs (519.3 ms vs 153.7 ms 1-thread / 24.3 ms 8-thread Nucleics Align). In a separate batch benchmark on a 16-core host (Ryzen AI Max+ 395), 3,000 DNA pairs (0.2–2 kb) ran in 0.50 s vs 4.93 s (`needleman_wunsch`), 0.69 s vs 7.48 s (`--linear`), 0.31 s vs 2.54 s (`--scoreonly`), and 0.63 s vs 5.05 s (`smith_waterman --maxhits 1`) comparing 16 threads against 1 thread, and a 2 kb BLOSUM62 protein pair dropped from 31 ms in original `seq-align` to 5 ms.

---

## 4. Discussion

Exact pairwise alignment has often forced a choice between general-purpose DP libraries that are predictable but slow and memory-hungry on long sequences, and output-sensitive wavefront aligners that are fast on near-identical inputs but degrade sharply on divergent ones. Nucleics Align combines compact traceback storage, difference-encoded striped SIMD with prefix-scan carry propagation, formally certified strip banding, Hirschberg/Myers-Miller linear space, and budgeted WFA inside one portable C11 codebase, giving predictable performance from kilobase proteins and viral genomes up to megabase DNA pairs.

Nucleics Align is not a replacement for seed-and-extend read mappers such as `minimap2` [@li2018minimap2] or database search tools such as `BLAST+` [@camacho2009blastplus], nor for `edlib` [@sosic2017edlib] when only unit-cost edit distance is needed. It targets workflows that require exact affine-gap global, overlap, or local alignment — including custom substitution matrices, wildcards, and fixed tie-breaking — from both the command line and C, Python, or Perl. Recent advances such as `Singletrack` [@lopezvillellas2026singletrack] inside `WFA2-lib` reduce wavefront backtrace memory 3x and remain complementary for wavefront-dominated workloads.

### Limitations and trade-offs

1. **Spin-waiting CPU overhead:** Multi-threaded single-pair fills synchronise column blocks by spinning, so total CPU time on 8 threads can exceed single-threaded comparators on full-matrix fills (e.g., 3.30 s CPU vs 2.36 s for `ksw2` on 50 kb unrelated `full`, despite 430 ms vs 2,359 ms wall time). By contrast, `--linear` on 8 threads uses *less* total CPU time than `ksw2` (0.86–1.17 s vs 2.36 s), and `--threads 1` avoids spin overhead altogether on busy multi-job servers.
2. **Local alignment remains $O(mn)$ in time:** Because `smith_waterman --maxhits 1` uses neither banding nor WFA, long similar local pairs remain quadratic in time (e.g., 9.6 s and 626 MiB on the 197 kb Mpox pair vs 7.6 ms for global `--wfa`).
3. **WFA crossover and scoring scope:** `--wfa` applies only when scoring reduces to one match score, one mismatch score, and standard affine gap penalties, and on our tested pairs at $\ge 1\%$ divergence from 50 kb upward it exhausts its cell budget and falls back to the matrix aligner. It is therefore opt-in.
4. **Memory vs tie-breaking on large pairs:** Full compact alignment on 580 kb *M. genitalium* still uses 4.33 GiB, and strip-kept mode on 1 Mb 1% pairs incurs refill overhead (8.10 s vs 6.40 s for `--linear` in Table 3; 7.92 s vs 6.33 s in the separate strip run). `--linear` and `--wfa` avoid this memory cost while returning the identical optimal score, but may select a different equal-scoring alignment when ties exist.
5. **Single-piece affine gaps and platform dispatch:** Nucleics Align implements Gotoh's single-piece affine gap model rather than the dual-affine convex gap cost used by `minimap2`/`ksw2` [@li2018minimap2]. Automatic AVX2/SSE4.1/SSE2 multi-versioning applies on x86-64 Linux; macOS and non-Linux x86-64 builds require `NATIVE=1` for AVX2/AVX-512.

### Future work

Two planned extensions build directly on the current architecture: combining Myers' bit-vector passes [@myers1999fast] with Hirschberg's divide-and-conquer [@hirschberg1975linear] so unit-cost `--linear` gains 64-cell-per-word throughput for full alignments, and exposing `--wfa` and `--linearabove` in the Python wrapper.

---

## 5. Availability, acknowledgements, and declarations

- **Code and data availability:** Nucleics Align is released under the BSD 3-Clause licence at https://github.com/DanielTillett/nucleics-align (citation metadata in `CITATION.cff`; a versioned release tag and Zenodo archive DOI will be minted when this manuscript is posted). All real biological sequences in Table 1 are publicly available from NCBI RefSeq/GenBank under the listed accessions and can be reproduced with `benchmarks/real_sequences.sh`; the ablation and 5-seed comparator suite are reproduced by `benchmarks/ablation.c`, `benchmarks/gen_pairs.py`, and `benchmarks/compare_harness.sh`.
- **Acknowledgements:** Nucleics Align is derived from Isaac Turner's `seq-align` (Unlicense) and bundles his `string_buffer` (CC0 1.0), `seq_file` (public domain), and `sort_r` (public domain) libraries. Substitution matrices in `scoring/` are from the NCBI BLAST matrix archive.
- **Funding:** [Placeholder: None / Funding source.]
- **Competing interests:** Daniel Tillett is the founder and director of Nucleics Pty Ltd; Nucleics Align is released under the open-source BSD 3-Clause licence with no commercial restrictions or patents.
- **Declaration of AI-assisted technologies:** Large language models (Claude, Codex, and Gemini) were used during software optimisation, code auditing, testing, and manuscript drafting. All algorithms, benchmark measurements, test verification, and scientific claims were reviewed and verified by the author, who takes full responsibility for this work.

---

## References

*(Rendered from `paper/references.bib` via BibTeX / Pandoc `--citeproc`.)*
