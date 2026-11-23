---
title: "Nucleics Align: fast, exact affine-gap pairwise sequence alignment for DNA and protein sequences up to megabase scale"
author: "Daniel Tillett^[Affiliation placeholder, City, Country; email@example.org; ORCID: 0000-0000-0000-0000]"
date: "23 Nov 2026"
bibliography: references.bib
---

## Abstract

**Background:** Exact pairwise sequence alignment with affine gap penalties (Needleman--Wunsch, Smith--Waterman, Gotoh) remains necessary whenever heuristic mappers may miss optimal alignments or when reproducible scores and fixed tie-breaking are required. Classic scalar tools such as Isaac Turner's `seq-align` are simple and exact, but aligning a 20 kb DNA pair takes ~11 s and 4.5 GiB of RAM, and a 50 kb pair has an estimated ~28 GiB matrix footprint. Modern libraries (`parasail`, `ksw2`, `SSW`, `WFA2`, `edlib`, `SeqAn`, `EMBOSS`) each address parts of this space, but standalone C tools combining general substitution matrices, single-pair multi-threading, linear-space global alignment at high divergence, and $O(s)$-space bidirectional wavefront alignment remain scarce.

**Results:** I present Nucleics Align (v1.2.0), a standalone C11 library and CLI suite derived from `seq-align`. It combines 1-byte-per-cell traceback storage, portable striped SIMD kernels (SSE2, SSE4.1, AVX2, AVX-512, Arm64 NEON) with prefix-scan carries and cache-line-padded multi-threading, 8- and 16-bit difference-encoded lanes, formally certified diagonal banding with strip pruning and strip-kept memory, linear-space Hirschberg/Myers-Miller alignment with row reuse (`--linear`), bidirectional wavefront alignment with $O(s)$ wavefront state (`--wfa`), and Myers' bit-vector edit distance for unit-cost `--scoreonly`. Every default optimisation reproduces the exact optimal score and tie-breaking of this codebase's scalar Gotoh baseline (optional `--linear` and `--wfa` guarantee the identical optimal score), verified across ~2.0 million unit-test checks (including ASan/UBSan) and 63 CLI fixtures (75 runs). On native Linux (8-core AMD Ryzen 7 5700G), Nucleics Align is 45–424x faster than original `seq-align` (1-seed baseline vs 5-seed median) on 20 kb global DNA on 1 thread (25.9–248.2 ms vs 11.0–12.1 s) and 153–414x faster on 8 threads (26.5–73.7 ms). On 1 thread, `--wfa` aligns 20 kb and 50 kb 1% pairs in 11.3 ms (2.8 MiB) and 47.6 ms (3.3 MiB), 1.6x faster and 10–37x smaller than `WFA2` high-memory mode (while `WFA2 BiWFA` is 1.8–2.6x faster in wall time at 4.0–5.4 MiB); `--linear` is 2.3–14.8x faster than `ksw2` across all six 20–50 kb conditions using 2.8–6.5 MiB (138–690x less RAM). A fresh-process nine-configuration ablation quantifies each layer's runtime and peak RAM. On real RefSeq/GenBank pairs (8 threads), `--wfa` aligns 197 kb Mpox genomes (99.92% identity) in 8.1 ms / 4.3 MiB, `--linear` aligns SARS-CoV-2 vs SARS-CoV-1 (29.9 kb, 79.54% identity) in 127.0 ms / 4.5 MiB and 580 kb *Mycoplasma genitalium* chromosomes in 2.28 s / 29.2 MiB, and five BLOSUM62 protein orthologues (390–3,066 aa) run 1.3–3.3x faster on 1 thread (up to 4.2x on 8 threads) than `parasail` AVX2 and up to 22.2x (1 thread) or 28.4x (8 threads) faster than `seq-align`. On simulated 1 Mb DNA, `--wfa` aligns 0.1% substitution pairs in 186.5 ms (11 MiB) and `--linear` aligns 1% in 7.19 s (50 MiB).

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
3. **Default global alignment:** When `--wfa` is enabled, it first attempts budgeted bidirectional wavefront alignment; if `--linearabove <n>` is set and the longer sequence has at least $n$ letters, it attempts linear-space alignment. If either is unsupported or does not complete, execution falls back to the matrix aligner. The matrix aligner uses a compact 1-byte-per-cell traceback fill whenever `opts.striped` is enabled and the scoring options are supported by the striped kernels: it first attempts a short budgeted score-only WFA pass (default budget 16 diagonal updates per letter of `seq_a`) or 12-mer estimate to size an exact band (Section 2.4), falling back to a full difference-encoded matrix fill if banding does not pay, and to the scalar three-matrix fill if `opts.striped` is disabled or scoring options are unsupported by the striped kernels. If `--linearfallback` is enabled and matrix allocation fails or overflows `size_t`, it falls back to linear space.

### 2.2 Compact traceback storage and deterministic tie-breaking

The baseline Gotoh implementation stores three 32-bit score matrices (`MATCH`, `GAP_A`, `GAP_B`), requiring 12 bytes per cell (~4.47 GiB for $20{,}000 \times 20{,}000$ cells). In Nucleics Align, compact mode (default in the CLI whenever `opts.striped` is enabled) packs the three 2-bit predecessor state indicators into one byte per cell while keeping only two score rows in RAM during the forward fill, cutting matrix memory 12-fold. Above 1 MiB, traceback buffers grow by 1.25x rather than power-of-two doubling to avoid over-allocating multi-gigabyte matrices.

Whenever multiple paths tie for the optimal score, traceback in every full-matrix path (scalar matrices, compact traceback bytes, exact bands, strip-kept bands, and difference-encoded lanes) follows the same deterministic rule. In Gotoh traceback, the current matrix state fixes the spatial step (`GAP_A` steps to $(x, y-1)$, `GAP_B` to $(x-1, y)$, `MATCH` to $(x-1, y-1)$), and at that predecessor cell the next state is chosen in strict priority order `GAP_A`, `GAP_B`, `MATCH`. Global alignment starts at $(|A|, |B|)$ using the same state priority; Smith-Waterman always starts in `MATCH` at the highest-scoring cell (breaking ties by lowest end coordinate $x$, then $y$). This priority tends to push gaps toward the right end of runs (e.g., `AAA` vs `AA` gives `AA-`, and `ACCA` vs `ACA` gives `AC-A`). By contrast, `--linear` and `--wfa` always return the identical optimal score but may pick a different equal-scoring alignment due to midpoint or wavefront tie-breaking.

### 2.3 Portable striped SIMD, prefix-scan carries, and narrow/difference-encoded lanes

All vectorised kernels use portable GCC/clang vector extensions rather than architecture-specific intrinsics. On x86-64 Linux, GCC/clang function multi-versioning (`target_clones` for `avx2`, `sse4.1`, and `default` SSE2) compiles all three variants into one portable binary (8-lane 32-bit vectors under GCC; 4-lane under Clang dispatch, as Clang rejects 8-lane vectors inside `target_clones` multi-versioned functions). Building with `NATIVE=1` (`-march=native`) enables 16 × 32-bit AVX-512 lanes where present, and Arm64 builds (including Apple silicon) use NEON automatically. On non-Linux x86-64 builds without `NATIVE=1`, multi-versioning is disabled and the build uses baseline SSE2 (with linear-space passes falling back to scalar code).

The striped fill follows Farrar's layout [@farrar2007striped], precomputing a per-query substitution profile for the active alphabet. Drawing on prefix-scan carry propagation [@khajehsaeed2010acceleration; @daily2016parasail; @liu2023block], Nucleics Align precomputes inter-lane vertical-gap carries in a single scan rather than iterating multi-sweep lazy-F loops across banded rows. To increase lane density without losing exactness:
- **Global full traceback fills** store vertical and horizontal score differences [@suzuki2018introducing], keeping values in 8-bit lanes when the per-step score range fits (including standard DNA scoring) and 16-bit lanes for wider matrices such as BLOSUM62 at any sequence length.
- **Local fills** begin in 8-bit lanes (four times as many lanes as 32-bit vectors) and widen to 16- or 32-bit lanes only if a score threatens to saturate.
- **Linear-space passes** use 16-bit relative score blocks.

### 2.4 Exact certified banding, strip pruning, and strip-kept memory

Let $m = |A|$, $n = |B|$, $D = m - n$, $s_{\max} = \max\big(0, \max_{a,b} S(a,b)\big)$, and assume standard affine gap scoring $g_{\text{open}} + \ell \cdot g_{\text{extend}}$ with fully penalised terminal gaps, non-positive penalties $g_{\text{open}} \le 0$ and $g_{\text{extend}} \le 0$, and $s_{\max} - 2 g_{\text{extend}} > 0$. Before filling a banded global matrix, Nucleics Align first attempts a budgeted score-only WFA pass [@marcosola2021fast] (default 16 diagonal updates per letter of $A$); when it finishes, it yields the exact optimal score. Otherwise, a fast 12-mer sampling routine produces a deliberately high heuristic score estimate used solely to decide whether a band is worth trying; the final alignment never depends on that estimate.

Exactness is certified mathematically from the score achieved inside the closed band $B(b) = \{(x, y) : \min(0, D) - b \le x - y \le \max(0, D) + b\}$. Let $S_{\text{out}}(b)$ denote the maximum score of any global path that leaves $B(b)$ (or $-\infty$ if no such path exists). Any global path visiting a cell outside $B(b)$ requires at least $I = |D| + 2(b + 1)$ gap characters and at most $(m + n - I)/2 = \min(m, n) - (b + 1)$ aligned residue pairs; because $g_{\text{open}} \le 0$ and $s_{\max} \ge 0$, every path leaving $B(b)$ satisfies:

$$S_{\text{out}}(b) \le s_{\max}\big(\min(m, n) - (b + 1)\big) + g_{\text{extend}}\big(|D| + 2(b + 1)\big)$$

Consequently, once a banded fill finds a feasible global alignment of score $L$, every path leaving $B(b)$ satisfies the strict inequality $S_{\text{out}}(b) < L$ (excluding both better and tied optimal paths outside $B(b)$) whenever:

$$b \ge \left\lfloor \frac{s_{\max}\min(m, n) + g_{\text{extend}}|D| - L}{s_{\max} - 2 g_{\text{extend}}} \right\rfloor$$

Excluded predecessor cells outside $B(b)$ are set to the DP floor so no cell score inside $B(b)$ is ever overestimated. When this inequality holds, every optimal and tied-optimal path stays within the closed band $B(b)$, preserving byte-identical full-matrix traceback. If the certificate fails, the band is widened and re-run until the inequality holds, or abandoned for a full matrix fill (which is exact unconditionally).

Inspired by `A*PA2` [@grootkoerkamp2024apa2] and Spouge's computational volumes [@spouge1989speeding], the band is filled in horizontal strips. Once a feasible lower bound $L$ is known (from the exact WFA score or from a prior banded pass) and $g_{\text{extend}} < 0$, for each cell $(x, y)$ on the bottom row of a strip with best state score $v = \max(M[x,y], G_A[x,y], G_B[x,y])$ and remaining lengths $r_i = m - x$, $r_j = n - y$, any completion to $(m, n)$ scores at most:

$$U(x, y) = v + s_{\max}\min(r_i, r_j) + g_{\text{extend}}|r_i - r_j|$$

Because $s_{\max} \ge 0$, $U(x, y)$ remains a valid upper bound even when all substitution scores are negative. Cells with $U(x, y) < L$ cannot lie on any optimal or tied-optimal path and are pruned; each surviving cell on diagonal $d = x - y$ can reach at most $k = \lfloor (U(x, y) - L) / (-2 g_{\text{extend}}) \rfloor$ diagonals beyond $[\min(d, D), \max(d, D)]$. When the estimated band footprint exceeds one quarter of physical RAM (falling back to 1 GiB when physical memory cannot be queried), the aligner switches automatically to strip-kept mode, retaining only strip boundary rows during the forward pass and refilling one strip's traceback bytes at a time during traceback. For a 1 Mb, 1% nominal mutation global alignment on a 62 GiB host, this avoids an estimated ~35 GiB band allocation and runs in 319 MiB while preserving full-matrix tie-breaking.

### 2.5 Budgeted bidirectional wavefront alignment (`--wfa`)

When `--wfa` is requested, Nucleics Align converts affine match/mismatch/gap scoring into non-negative penalty costs [@eizenga2022improving] and advances wavefronts by penalty score [@marcosola2021fast]. For short pairs ($|A| + |B| \le 2{,}048$) or non-progressing boundary splits, it traces back directly through retained wavefronts. For larger pairs, v1.2.0 implements bidirectional wavefront alignment (`BiWFA`) [@marcosola2023optimal] with $O(s)$ auxiliary wavefront state (in addition to $O(|A| + |B|)$ storage for the recoded forward/reversed sequences and output buffers): forward and reverse ring-buffer wavefronts advance alternately until they meet across match or affine gap states (accounting for the shared gap-open charge on $I$--$I$ and $D$--$D$ breakpoints), recursively splitting the problem down to the retained-wavefront base case. Each wavefront step decouples its branch-free interior compute loop (compiled with `target_clones` for AVX2, SSE4.1, and SSE2) from its 8-byte word-comparison match-extension loop, and skips overlap checks until the sum of forward and reverse furthest offsets reaches $|B|$. To fall back quickly on divergent pairs, the global aligner caps total diagonal updates (across both bidirectional and base-case steps) at $\min(2^{28},\; 256 \cdot (|A| + |B|))$; if scoring does not reduce to one match score, one mismatch score, and standard affine gap penalties, or if the work budget is exhausted, execution falls back to the matrix aligner.

### 2.6 Linear-space global alignment (`--linear`) and Myers bit-vector score-only

The `--linear` aligner implements Hirschberg's divide-and-conquer algorithm [@hirschberg1975linear] extended to Gotoh's affine gap model [@myers1988optimal], combined with certified banded striped passes [@chao1992aligning; @farrar2007striped], 16-bit relative blocks, and two-level intermediate row reuse so child sub-problems skip redundant passes. Working score-row memory is $O(|B|)$ (with constant-factor overhead for row reuse and per-thread buffers), plus $O(|A| + |B|)$ for the returned alignment strings.

For `--scoreonly`, when scoring is unit-cost edit distance (`0 / -1 / 0 / -1`), banding does not pay, and execution is single-threaded (or SIMD is unavailable), Nucleics Align dispatches to Myers' bit-vector algorithm [@myers1999fast], computing 64 cells per machine word using a small constant number of bitwise operations.

### 2.7 Best-hit local alignment (`smith_waterman --maxhits 1`)

When only the single highest-scoring local hit is requested (`--maxhits 1`), Nucleics Align scans forward in $O(|B|)$ memory to locate the optimal end cell and score, then scans in reverse from that cell in $O(|B|)$ memory to locate the start cell. Only the resulting sub-rectangle is filled for traceback, with periodic row checkpoints if the rectangle exceeds $2^{28}$ cells. During checkpointed strip refills, v1.2.0 reuses the precomputed substitution table and selects 8- or 16-bit SIMD lanes whenever the known optimal score fits safely, preserving byte-identical traceback.

### 2.8 Multi-threading, error handling, and verification

By default, the CLI tools use half of the host's logical CPUs (`--threads <n>` overrides this). Large single alignments (exceeding $4 \times 10^6$ cells and 1,024 columns per thread block) split columns into pipelined wavefront blocks across threads, while file batches (`--file`, `--files`) align independent pairs concurrently and emit results in strict input order so output is byte-identical to `--threads 1`. Per-block progress counters are 64-byte aligned (`_Alignas(64)`) to prevent cache-line false sharing, waiters cache the last observed row to avoid redundant atomic loads, and progress is published every 16 rows with x86 `pause` hints during short spins.

In the C library, every major aligner provides a checked entry point returning an error status (or `1`/`0` with an `aln_error_t` parameter for `smith_waterman_fetch_checked()`) instead of calling `exit()`. Continuous integration runs on Linux (x86-64 and Arm64) and macOS with GCC and Clang, including AddressSanitizer and UndefinedBehaviorSanitizer (`SANITIZE=1 make test`), ~2.0 million unit-test assertion checks (many on randomised inputs comparing every optimised path against the scalar Gotoh baseline), and 63 CLI regression fixtures (75 test runs, as 12 batch fixtures run at both 1 and 4 threads).

---

## 3. Results

All main benchmarks were run on a native Linux desktop (AMD Ryzen 7 5700G, 8 cores / 16 logical CPUs, 62 GiB RAM, Pop!_OS 24.04 Linux, GCC 13.3.0, default `make` `-O3` build; CPU frequency governor `powersave`, SMT active, no manual core pinning). The archived benchmark logs `benchmarks/results/pop-os-v1.2.0.txt` and `benchmarks/results/pop-os-v1.2.0-journal.txt` (Nucleics Align v1.2.0, commit `aaf0754`) record all measurements for Tables 1–5, including per-seed wall times, median CPU times and peak RSS, local per-seed scores and SSW flags, independent cross-tool CIGAR/alignment re-scoring checks, full SHA-256 digests for all synthetic/protein FASTAs and Nucleics Align CLI binaries, and 12-character SHA-256 prefixes for Table 1 FASTAs. CLI timings include process startup, input parsing, alignment, and output writing (so sub-10 ms runs include ~2–4 ms of process startup overhead); the C library ablation (Table 2) times the alignment call directly via `clock_gettime(CLOCK_MONOTONIC)`. Memory figures report peak resident set size (RSS) from GNU `time -v` in binary MiB/GiB (`1 MiB = 1024^2` bytes; `1 GiB = 1024 MiB`), and CPU times come from GNU `time -v` (10 ms resolution, where `0.00 s` is displayed as `<10 ms`). Unless noted, global DNA runs used `+1` match, `-2` mismatch, `-4` gap open, `-1` gap extend; local DNA runs used `+2 / -2 / -2 / -1`; and protein runs used `--scoring BLOSUM62` (`-10 / -1`). All global modes and comparators in Tables 1, 2, 4, and 5 were programmatically asserted to return the identical optimal score on every completed pair.

### 3.1 Real biological sequences

Seven public NCBI RefSeq/GenBank pairs spanning mitochondrial genomes, coronaviruses, orthopoxviruses, whole bacterial chromosomes, and mammalian proteins were timed with 8 default threads on v1.2.0 (best of 3 runs for pairs $\le 200$ kb, 1 run for *M. genitalium*; reproducible via `benchmarks/real_sequences.sh`; Table 1). Circular mitochondrial and bacterial genomes were aligned as deposited in RefSeq without origin rotation; the 580 kb *M. genitalium* pair serves as a whole-chromosome scaling stress test rather than a rearrangement-aware comparison. Reported sequence identities belong to the `full` traceback and count case-insensitive equal symbols divided by total alignment columns (so `N/N` counts as an identity match under `--pretty`, whereas under default DNA alignment scoring ambiguous bases such as `N` are scored as mismatches).

**Table 1. Performance on real biological sequence pairs (`pop-os`, AMD Ryzen 7 5700G, 8 default threads, v1.2.0 `aaf0754`; wall ms / peak MiB).** DNA global modes used `+1/-2/-4/-1`; DNA `local --maxhits 1` used `+2/-2/-2/-1`; BRCA1 used `BLOSUM62` (`-10/-1`). Bold marks the fastest global mode returning a full alignment. `†` = `--wfa` exceeded its work budget and fell back to matrices; `‡` = `--wfa` unsupported for substitution matrix and fell back; `—` = not run.

| Pair (NCBI Accessions) | Lengths | Identity | `full` (ms / MiB) | `--wfa` (ms / MiB) | `--linear` (ms / MiB) | `--scoreonly` (ms / MiB) | `local --maxhits 1` (ms / MiB) | Optimal global score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Human/Neanderthal mtDNA (`NC_012920.1`, `NC_011137.1`) | 16.6 / 16.6 kb | 98.72% | 17.6 / 11.5 | **5.9** / 2.9 | 23.0 / 3.1 | 6.9 / 3.1 | 89.3 / 19.2 | 15,928 |
| Human/Chimp mtDNA (`NC_012920.1`, `NC_001643.1`) | 16.6 / 16.6 kb | 85.19% | **50.0** / 266.4 | 81.7† / 266.8 | 55.9 / 4.1 | 20.2 / 3.7 | 85.6 / 19.0 | 10,632 |
| SARS-CoV-2 Wuhan/Omicron (`NC_045512.2`, `OL672836.1`) | 29.9 / 29.7 kb | 99.03% | 22.2 / 16.2 | **5.7** / 3.0 | 46.6 / 3.6 | 12.8 / 3.4 | 246.5 / 41.9 | 29,240 |
| SARS-CoV-2/SARS-CoV-1 (`NC_045512.2`, `NC_004718.3`) | 29.9 / 29.8 kb | 79.54% | 160.9 / 856.9 | 220.4† / 856.9 | **127.0** / 4.5 | 31.9 / 4.0 | 226.7 / 42.0 | 11,467 |
| Mpox Ref/2022 isolate (`NC_063383.1`, `ON563414.3`) | 197.2 / 197.2 kb | 99.92% | 91.0 / 63.6 | **8.1** / 4.3 | 287.1 / 10.9 | 49.9 / 10.4 | 9,376.7 / 625.8 | 196,823 |
| *M. genitalium* G37/M2288 (`NC_000908.2`, `CP003773.1`) | 580.1 / 579.6 kb | 99.09% | **2,123.3** / 4,430.9 | 3,061.6† / 4,431.4 | 2,276.2 / 29.2 | 574.0 / 27.3 | — | 564,775 |
| Human/Mouse BRCA1, BLOSUM62 (`NP_009225.1`, `NP_033894.3`) | 1,863 / 1,812 aa | 55.76% | 7.3 / 6.4 | 6.9‡ / 6.4 | **6.5** / 2.6 | 4.1 / 2.5 | 6.9 / 6.9 | 4,827 |

On closely related viral genomes, `--wfa` (BiWFA) aligns two 197 kb Mpox genomes (99.92% identity) in 8.1 ms and 4.3 MiB (11.3x faster and 14.8x smaller than `full`, and 6.2x faster than `--scoreonly`), SARS-CoV-2 Wuhan vs Omicron (99.03% identity) in 5.7 ms and 3.0 MiB (3.9x faster than `full`), and Human vs Neanderthal mtDNA (98.72% identity) in 5.9 ms and 2.9 MiB. On divergent SARS-CoV-2 vs SARS-CoV-1 (79.54% identity), `--linear` is faster than `full` (127.0 ms vs 160.9 ms) while using 190x less RAM (4.5 MiB vs 856.9 MiB). On Human vs Chimp mtDNA (85.19% identity), `--linear` runs in 55.9 ms vs 50.0 ms for `full` while cutting RAM 65.0-fold (4.1 MiB vs 266.4 MiB). Two 580 kb *M. genitalium* chromosomes (99.09% identity, 575,139 matches, 4,101 mismatches, 1,154 gap columns) align in 2.12 s with `full` (4.33 GiB) or 2.28 s with `--linear` using 29.2 MiB (152x less RAM). On Human vs Mouse BRCA1 under BLOSUM62 (`-10 / -1`), `--linear` finishes in 6.5 ms (2.6 MiB) and `full` in 7.3 ms (6.4 MiB).

### 3.2 Fresh-process nine-configuration ablation

To quantify each layer's runtime and peak RAM without cross-run allocator effects, `benchmarks/ablation.c` timed nine configurations of the C11 library API, each in its own fresh process with GNU `time -v` RSS measurement, across five independently generated 20 kb DNA pairs (`benchmarks/gen_pairs.py`, seeds 1–5, substitution probability $p \in \{1\%, 10\%\}$ with single-base insertions and deletions each at probability $p/5$; Table 2). Rows 1–6 add successive optimisations to the matrix kernel (compact 1-byte traceback requires `opts.striped = true`, so row 2 introduces both together; rows 4 and 5 separate 12-mer band sizing from budgeted score-only WFA sizing); rows 7–9 are alternative linear-space and bidirectional wavefront modes. Every configuration returned the exact same optimal score on every pair.

**Table 2. Fresh-process nine-configuration ablation on 20 kb DNA (`pop-os`, AMD Ryzen 7 5700G, v1.2.0 `aaf0754`, 5 independent seeds, C11 library API; median [min–max] wall ms and median peak RSS in MiB).** At 10% substitution, rows 4 and 5 skip banding via the cost check and execute a full 8-bit difference fill plus the sizing probe (`385.5–385.8 MiB`). `†` = `--wfa` hit its work budget and fell back to matrices.

| Configuration | 20 kb, 1% sub. (ms / MiB) | 20 kb, 10% sub. (ms / MiB) |
|---|---:|---:|
| 1. Scalar 3-matrix baseline (1 thr.) | 2,730.3 [2,724.6–2,734.9] / 4,580.2 | 2,838.4 [2,826.8–2,865.3] / 4,580.2 |
| 2. + Compact 1-byte & striped 32-bit SIMD (1 thr.) | 443.0 [441.1–459.1] / 385.7 | 431.0 [425.3–443.5] / 385.8 |
| 3. + Difference-encoded 8-bit lanes (1 thr.) | 237.1 [236.2–243.1] / 385.4 | 237.7 [233.9–242.2] / 385.4 |
| 4. + Certified banding, 12-mer sizing only (1 thr.) | 22.8 [20.4–24.0] / 17.3 | 242.9 [239.8–247.4] / 385.5 |
| 5. + Certified banding, WFA sizing enabled (12-mer fallback, 1 thr.) | 24.2 [22.9–24.9] / 17.6 | 243.3 [241.5–248.2] / 385.8 |
| 6. + Multi-threading (8 thr.) | 24.2 [22.5–24.6] / 17.6 | 69.6 [65.9–75.1] / 387.4 |
| 7. Linear space `--linear` (1 thr.) | 39.8 [39.2–40.1] / 2.8 | 177.1 [151.9–183.8] / 3.2 |
| 8. Linear space `--linear` (8 thr.) | 28.7 [28.2–28.9] / 3.2 | 64.1 [63.0–69.6] / 3.6 |
| 9. Bidirectional wavefront `--wfa` (1 thr.) | 8.6 [8.0–10.0] / 2.9 | 278.5† [275.6–280.0] / 386.0 |

Nucleics Align's scalar 3-matrix baseline (2.73–2.84 s, 4,580 MiB) is already ~4x faster than original `seq-align` (11.0–11.3 s in Table 4), likely reflecting replacement of per-cell `scoring_lookup()` calls with a precomputed `sub_table_t` and allocation differences (not isolated separately here). Compact 1-byte traceback combined with striped 32-bit SIMD cuts peak RAM 11.9-fold (4,580 to 386 MiB) and speeds up the fill 6.2–6.6x (431–443 ms). Difference-encoded 8-bit lanes add a further 1.81–1.87x speedup (237–238 ms). At 1% substitution, certified exact banding yields a 10.4x jump (237.1 to 22.8 ms, 17.3 MiB) and BiWFA `--wfa` reaches 8.6 ms in 2.9 MiB. At 10% substitution, the band cost check skips banding and adds ~5–6 ms of probe overhead over the unbanded 8-bit fill (243.3 vs 237.7 ms, matching Table 4's 248.2 ms CLI time once ~5 ms of process startup is included); 8-thread parallelism adds 3.5x (69.6 ms), and `--linear` reaches 177.1 ms on 1 thread and 64.1 ms on 8 threads in 3.2–3.6 MiB.

### 3.3 Scaling from 5 kb to 1 Mb on simulated DNA

Simulated DNA pairs from 5 kb to 1 Mb were generated deterministically with `benchmarks/gen_pairs.py` (`seed=11`, substitution probability $p \in \{0.1\%, 1\%, 5\%, 10\%\}$ with single-base insertions and deletions each at probability $p/5$, or unrelated pairs) and timed on `pop-os` (v1.2.0 `aaf0754`, 8 default threads; reproducible via `benchmarks/journal_bench.sh`; Table 3). Every global mode — including matched single-threaded `WFA2 BiWFA` (`v2.3.5`) runs at 100 kb, 300 kb, and 1 Mb (0.1% and 1%) — was programmatically asserted to return the identical optimal score on every completed pair.

**Table 3. Scaling on simulated DNA from 5 kb to 1 Mb (`pop-os`, AMD Ryzen 7 5700G, 8 default threads, v1.2.0 `aaf0754`; wall ms / peak MiB).** Standard lengths (5–50 kb): best of 3 runs; large lengths (100 kb – 1 Mb): 1 run. DNA global modes used `+1/-2/-4/-1`; `local --maxhits 1 --minscore 1` used `+2/-2/-2/-1`. Bold marks the fastest Nucleics Align global mode returning a full alignment. `†` = `--wfa` hit its work budget and fell back to matrices; `—` = not run.

| Length | Sub. rate | `full` (ms / MiB) | `--wfa` (ms / MiB) | `--linear` (ms / MiB) | `--scoreonly` (ms / MiB) | `local1` / `WFA2 BiWFA` (1 thr.) |
|---|---|---:|---:|---:|---:|---:|
| 5 kb | 1% | 6.2 / 4 | **2.9** / 3 | 7.1 / 2 | 3.2 / 2 | 21.5 / 31 *(local1)* |
| 5 kb | 10% | 14.3 / 27 | 25.0† / 27 | **10.8** / 3 | 4.8 / 2 | 20.0 / 31 *(local1)* |
| 5 kb | unrelated | **12.6** / 27 | 23.0† / 27 | 12.9 / 3 | 5.0 / 2 | 29.3 / 30 *(local1)* |
| 20 kb | 1% | 26.3 / 16 | **11.2** / 3 | 30.2 / 3 | 8.6 / 3 | 137.0 / 24 *(local1)* |
| 20 kb | 10% | **69.1** / 387 | 109.6† / 388 | **69.1** / 4 | 15.3 / 4 | 103.6 / 24 *(local1)* |
| 20 kb | unrelated | **68.7** / 387 | 107.9† / 387 | 75.4 / 4 | 20.2 / 4 | 133.4 / 24 *(local1)* |
| 50 kb | 1% | 125.2 / 83 | **47.5** / 3 | 105.7 / 4 | 33.3 / 4 | 576.7 / 87 *(local1)* |
| 50 kb | 10% | 421.3 / 2,397 | 523.5† / 2,398 | **218.2** / 6 | 62.6 / 6 | 543.8 / 87 *(local1)* |
| 50 kb | unrelated | 405.3 / 2,395 | 501.1† / 2,395 | **274.2** / 8 | 97.3 / 6 | 617.0 / 86 *(local1)* |
| 100 kb | 0.1% | 71.2 / 38 | **10.7** / 3 | 213.8 / 6 | 40.5 / 6 | 6.4 / 5 *(BiWFA)* |
| 100 kb | 1% | 496.8 / 308 | **163.4** / 4 | 303.4 / 7 | 108.5 / 6 | 56.1 / 8 *(BiWFA)* |
| 100 kb | 5% | **439.5** / 1,474 | 692.5† / 1,474 | 491.6 / 9 | 143.7 / 8 | — |
| 300 kb | 0.1% | 369.1 / 282 | **27.3** / 5 | 545.9 / 15 | 127.7 / 15 | 12.2 / 7 *(BiWFA)* |
| 300 kb | 1% | **904.9** / 2,693 | 1,482.0† / 2,693 | 1,157.2 / 17 | 271.2 / 16 | 384.8 / 16 *(BiWFA)* |
| 300 kb | 5% | 3,470.7 / 206 *(strips)* | 4,036.2† / 206 | **2,850.9** / 22 | 1,056.5 / 20 | — |
| 1 Mb | 0.1% | 4,034.8 / 3,125 | **186.5** / 11 | 3,183.4 / 46 | 1,079.7 / 43 | 74.4 / 16 *(BiWFA)* |
| 1 Mb | 1% | 9,232.1 / 342 *(strips)* | 9,993.4† / 342 | **7,193.2** / 50 | 2,638.4 / 47 | 4,199.2 / 44 *(BiWFA)* |
| 1 Mb | 5% | — | — | **94,052.5** / 110 | 59,388.0 / 94 | — |

With v1.2.0 BiWFA, `--wfa` aligns 50 kb 1% in 47.5 ms / 3 MiB, 100 kb 1% in 163.4 ms / 4 MiB (3.0x faster and 77x smaller than `full`), and 1 Mb 0.1% in 186.5 ms / 11 MiB (21.6x faster and 284x smaller than `full`, and 17.1x faster than `--linear`). Across 100 kb – 1 Mb low-divergence pairs, `WFA2 BiWFA` remains 1.7–2.9x faster in wall time (e.g., 74.4 ms vs 186.5 ms at 1 Mb 0.1%) while Nucleics Align `--wfa` uses less peak RSS (3–11 MiB vs 5–16 MiB). On longer divergent pairs, `full` automatically activates strip-kept mode (206 MiB at 300 kb 5%; 342 MiB at 1 Mb 1%) while preserving full-matrix tie-breaking, and `--linear` keeps peak RSS at 46–110 MiB even at 1 Mb and 5% true substitution rate (94.1 s).

### 3.4 Reproducible 5-seed comparison with existing exact aligners

Using the self-contained harness `benchmarks/compare_harness.sh` on native Linux `pop-os`, Nucleics Align (v1.2.0 `aaf0754`, default `make` `-O3`, portable AVX2/SSE4.1/SSE2 dispatch) was compared against clean pinned builds of original `seq-align` (`dc41988`, default `make`) [@turner2017seqalign], `ksw2` (`289609b`, `gcc -O3 -msse4.1`, 128-bit SSE4.1 `ksw_extz2_sse` with `w=-1`, `zdrop=-1`, `flag=0` for unbanded exact score + CIGAR) [@li2018minimap2; @ksw2repo], `WFA2-lib` (`v2.3.5` / `cf3eb92`, built via `make BUILD_EXAMPLES=0 BUILD_TOOLS=0 all` with `-O3 -march=native`, exact high-memory and ultralow `BiWFA` modes with `wf_heuristic_none`, `compute_alignment` CIGAR, affine penalties `match=0, mismatch=6, gap_open=8, gap_extend=3` mapped back via `score = (|A|+|B|-penalty)/2`) [@marcosola2021fast; @marcosola2023optimal; @wfa2lib], `edlib` (`v1.2.7` / `ec2310e`, `g++ -O3`, `EDLIB_MODE_NW`, `EDLIB_TASK_DISTANCE`, unit-cost global edit distance) [@sosic2017edlib], and `SSW` (`a66636b`, `gcc -O3`, `ssw_align` with `+2/-2/-2/-1` local scoring, score + CIGAR) [@zhao2013ssw; @sswrepo]. All tools were timed across five independently generated deterministic DNA pairs per condition (`benchmarks/gen_pairs.py`, seeds 1–5, single-base indels at $p/5$; 1 seed for original `seq-align` in all three 20 kb conditions and for `WFA2 high` on 20 kb unrelated; `WFA2` on 50 kb unrelated was skipped). Every global and unit-cost tool returned the exact same optimal score on every completed pair. For local alignment, Nucleics Align used `--minscore 1` so it reported the optimal hit on every seed, matching `SSW` identically on all 10% and unrelated seeds (while `SSW` saturated at 32,767 on all five 1% seeds; Table 4).

**Table 4. Reproducible 5-seed comparator benchmark (`pop-os`, AMD Ryzen 7 5700G, native Linux, v1.2.0 `aaf0754`; median wall ms [min–max] and median CPU ms in parentheses).** Bold marks the lowest wall time per column within each block. `†` = `--wfa` hit its work budget and fell back. `*` = `SSW` saturated at `INT16_MAX = 32,767` on all 5 seeds (true optimal scores 38,776–38,940; `ssw_flag = [1, 1, 0, 1, 0]`).

| Benchmark & tool | 1% sub. (ms) | 10% sub. (ms) | Unrelated (ms) | Median RSS (1% / 10% / rand) |
|---|---|---|---|---|
| **Global affine, 20 kb (`+1/-2/-4/-1`)** | | | | |
| `seq-align` (`dc41988`, 1 thr.) | 10,976 (10,960) | 11,254 (11,240) | 12,101 (12,080) | 4,577 / 4,586 / 4,580 MiB |
| Nucleics Align (default, 1 thr.) | 25.9 [24.6–27.0] (20) | 248.2 [246.3–249.1] (240) | 220.2 [219.4–223.9] (210) | 17.6 / 385.8 / 385.7 MiB |
| Nucleics Align `--linear` (1 thr.) | 42.0 [41.0–42.4] (30) | 172.7 [149.3–178.0] (160) | 191.4 [171.8–202.4] (180) | 2.8 / 3.2 / 3.8 MiB |
| Nucleics Align `--wfa` (1 thr.) | 11.3 [10.9–12.1] (<10) | 283.0† [281.0–285.7] (270) | 257.7† [256.2–260.3] (240) | 2.8 / 385.8 / 385.8 MiB |
| Nucleics Align (default, 8 thr.) | 26.5 [24.8–27.2] (20) | 73.7 [72.6–79.5] (520) | 73.6 [69.6–78.6] (520) | 17.6 / 387.4 / 387.3 MiB |
| Nucleics Align `--linear` (8 thr.) | 30.6 [29.9–31.4] (40) | **67.3** [65.3–68.0] (140) | **71.8** [70.4–74.7] (200) | 3.2 / 3.7 / 4.4 MiB |
| `ksw2` (`289609b`, SSE4.1, 1 thr.) | 415.7 [413.0–422.0] (410) | 422.9 [415.3–425.5] (410) | 432.2 [427.1–432.7] (420) | 523.2 / 523.2 / 523.4 MiB |
| `WFA2` (`v2.3.5` high, 1 thr.) | 18.4 [17.7–19.7] (10) | 865.6 [830.5–957.9] (850) | 9,555 (9,540) | 28.9 / 1,491 / 17,029 MiB |
| `WFA2 BiWFA` (`v2.3.5`, 1 thr.) | **6.3** [6.1–7.9] (<10) | 156.5 [152.1–170.1] (140) | 1,655.9 [1,634–1,673] (1,640) | 4.0 / 9.7 / 24.1 MiB |
| **Global affine, 50 kb (`+1/-2/-4/-1`)** | | | | |
| Nucleics Align (default, 1 thr.) | 117.3 [110.4–121.6] (100) | 1,486.4 [1,479–1,517] (1,480) | 1,336.9 [1,329–1,346] (1,330) | 83.4 / 2,390 / 2,392 MiB |
| Nucleics Align `--linear` (1 thr.) | 158.9 [154.5–163.5] (150) | 833.6 [825.8–910.7] (820) | 948.0 [939.1–1,062] (940) | 4.0 / 5.0 / 6.5 MiB |
| Nucleics Align `--wfa` (1 thr.) | 47.6 [43.3–53.3] (40) | *(not run)* | *(not run)* | 3.3 / — / — MiB |
| Nucleics Align (default, 8 thr.) | 117.0 [111.9–123.9] (100) | 423.8 [411.3–436.0] (3,120) | 413.5 [394.5–415.4] (3,140) | 83.5 / 2,394 / 2,395 MiB |
| Nucleics Align `--linear` (8 thr.) | 105.3 [102.5–107.1] (160) | **218.7** [214.2–221.4] (790) | **277.9** [276.8–292.7] (1,170) | 4.5 / 6.2 / 7.4 MiB |
| `ksw2` (`289609b`, SSE4.1, 1 thr.) | 2,357.3 [2,349–2,363] (2,340) | 2,342.9 [2,340–2,360] (2,330) | 2,355.3 [2,350–2,365] (2,340) | 2,762 / 2,760 / 2,761 MiB |
| `WFA2` (`v2.3.5` high, 1 thr.) | 74.2 [65.9–80.7] (60) | 5,102.6 [5,018–5,334] (5,090) | *(not run)* | 121.0 / 9,226 / — MiB |
| `WFA2 BiWFA` (`v2.3.5`, 1 thr.) | **18.6** [17.4–19.3] (10) | 857.9 [851.9–895.9] (850) | *(not run)* | 5.4 / 18.5 / — MiB |
| **Unit-cost edit distance, 20 kb / 50 kb (`0/-1/0/-1`, score-only, 1 thr.)** | | | | |
| Nucleics Align `--scoreonly` (20 kb) | 8.3 [8.0–8.3] | 20.8 [20.4–20.8] | 28.1 [28.0–29.3] | 2.7 / 2.8 / 2.8 MiB |
| `edlib` (`v1.2.7`, 20 kb) | **3.2** [3.2–3.5] | **6.7** [6.2–8.1] | **17.8** [17.5–17.9] | 3.9 / 3.9 / 4.0 MiB |
| Nucleics Align `--scoreonly` (50 kb) | 22.4 [21.5–22.5] | 101.1 [100.4–104.4] | 152.7 [151.9–154.8] | 3.8 / 4.1 / 3.8 MiB |
| `edlib` (`v1.2.7`, 50 kb) | **6.1** [5.8–6.9] | **18.2** [17.7–18.7] | **72.5** [72.4–74.2] | 4.1 / 4.1 / 4.1 MiB |
| **Best local hit, 20 kb (`--maxhits 1 --minscore 1`, `+2/-2/-2/-1`)** | | | | |
| Nucleics Align (1 thr.) | 487.1 [481.7–489.3] (480) | 321.8 [319.2–325.8] (310) | 575.2 [572.4–585.4] (570) | 22.1 / 21.6 / 21.2 MiB |
| Nucleics Align (8 thr.) | **138.8** [137.5–145.4] (650) | **105.9** [103.6–108.9] (480) | **131.7** [131.1–142.9] (740) | 24.6 / 23.9 / 24.2 MiB |
| `SSW` (`a66636b`, 1 thr.) | 2,182.3* [195.5–2,251.4] (2,170) | 126.2 [124.3–141.5] (120) | 523.2 [271.7–800.6] (510) | 884.5 / 7.8 / 95.6 MiB |

In matched single-thread runs, Nucleics Align `--linear` is 2.3–14.8x faster than `ksw2` (note `ksw2` runs 128-bit SSE4.1 unbanded for exactness, whereas Nucleics Align dispatches 256-bit AVX2) across all six 20 kb and 50 kb conditions using 138–690x less RAM. With v1.2.0 BiWFA, single-threaded `--wfa` aligns 20 kb 1% in 11.3 ms (2.8 MiB) and 50 kb 1% in 47.6 ms (3.3 MiB), beating `WFA2` high-memory mode by 1.6x while using 10–37x less RAM (and less RAM than `WFA2 BiWFA` at 4.0–5.4 MiB), while `WFA2 BiWFA` remains 1.8–2.6x faster in wall time at 1% substitution (6.3 ms and 18.6 ms). At 10% and unrelated divergence on 1 thread, `WFA2 BiWFA` is slightly faster than `--linear` at 20 kb 10% (156.5 vs 172.7 ms), the two match within overlapping ranges at 50 kb 10% (833.6 [825.8–910.7] vs 857.9 [851.9–895.9] ms), and `--linear` is 8.7x faster on 20 kb unrelated (191.4 vs 1,655.9 ms); on 8 threads, `--linear` drops 50 kb 10% and unrelated pairs to 218.7 ms and 277.9 ms (6.2–7.4 MiB).

On unit-cost edit distance (`--scoreonly` vs `edlib v1.2.7`), both tools agreed on every score across all 30 pairs, with `edlib` 1.6–5.5x faster (3.2–72.5 ms vs 8.3–152.7 ms) because `edlib` combines Myers' bit-vector algorithm with Ukkonen banding. To verify traceback validity beyond reported scores, `benchmarks/journal_bench.sh` independently re-scored the returned alignments or CIGAR strings and verified exact input span coverage on 20 kb (`1%`, `10%`, and `rand`, seed 1) across all five global tools (`full`, `--linear`, `--wfa`, `ksw2`, and `WFA2 BiWFA`), confirming exact agreement (`18,934`, `10,215`, and `-13,508`). On 20 kb local alignment (`--maxhits 1 --minscore 1` for Nucleics Align), all three local configurations returned identical scores on 10% (28,686–29,521) and unrelated pairs (3,241–3,332) with `ssw_flag = 0` on every seed; `SSW` is faster than 1-thread Nucleics Align at 10% (126.2 vs 321.8 ms) and similar on unrelated pairs (523.2 vs 575.2 ms), while 8-thread Nucleics Align is fastest across all three conditions (105.9–138.8 ms). At 1% substitution (`*`), the true optimal local score (38,776–38,940) exceeds SSW's signed 16-bit limit (`INT16_MAX = 32,767`; SSW has no 32-bit path, so all five seeds saturate at 32,767 and three of five seeds return `ssw_flag = 1`, which SSW defines as failure of banded traceback); these saturated 1% SSW runs are therefore reported as a failure-mode observation rather than an exact-alignment timing comparison.

### 3.5 Multi-pair real protein benchmark (BLOSUM62)

Because wavefront aligners such as `WFA2` and edit-distance aligners such as `edlib` do not support arbitrary protein substitution matrices, five real NCBI RefSeq orthologue pairs spanning 390 to 3,066 amino acids and 49.1% to 98.9% sequence identity (column identity of the default global traceback) were benchmarked under `BLOSUM62` (`-10 / -1`) comparing original `seq-align` (`dc41988`) and pinned `parasail` (`v2.6.2` / `600fb26`, `gcc -O3 -mavx2`, `parasail_nw_trace_striped_avx2_256_16` with `open=11, extend=1` matching Gotoh `-10/-1`, full CIGAR traceback) [@daily2016parasail] against Nucleics Align (`benchmarks/protein_bench.sh` and `benchmarks/journal_bench.sh`, median of 3 runs; Table 5). Every global tool returned the identical optimal score on every pair, and independent BLOSUM62 re-scoring of the returned alignments/CIGARs for `parasail`, `full`, and `--linear` confirmed exact scores and input span coverage across all five pairs.

**Table 5. Real protein orthologue benchmark under BLOSUM62 (`-10/-1`, `pop-os`, AMD Ryzen 7 5700G, v1.2.0 `aaf0754`; median wall ms / peak MiB over 3 runs).** Bold marks the fastest global mode returning a full alignment.

| Protein pair (NCBI Accessions) | Len. (aa) | Id. | `seq-align` (1t) | `parasail` (1t) | `full` (1t) | `full` (8t) | `--linear` (1t) | `--linear` (8t) | `--scoreonly` | `local1` | Score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Human/Mouse TP53 (`NP_000537.3`, `NP_035770.2`) | 393/390 | 76.77% | 5.8 / 3.5 | 2.5 / 2.0 | 2.0 / 2.3 | 2.2 / 2.4 | **2.0** / 2.3 | 2.1 / 2.3 | 1.9 / 2.4 | 2.3 / 2.6 | 1,540 |
| Human POLR2A / Yeast Rpo21 (`NP_000928.1`, `NP_010141.1`) | 1,970/1,733 | 49.09% | 85.6 / 40.6 | 14.1 / 8.2 | 6.7 / 6.3 | 6.8 / 6.5 | **5.7** / 2.5 | 6.1 / 2.6 | 4.1 / 2.5 | 7.0 / 6.9 | 4,630 |
| Human/Mouse BRCA1 (`NP_009225.1`, `NP_033894.3`) | 1,863/1,812 | 55.76% | 84.2 / 40.3 | 12.7 / 8.2 | 6.7 / 6.3 | 6.9 / 6.4 | **5.8** / 2.4 | 6.0 / 2.5 | 3.9 / 2.5 | 7.1 / 7.0 | 4,827 |
| Human/Mouse MTOR (`NP_004949.1`, `NP_064393.2`) | 2,549/2,549 | 98.94% | 153.6 / 76.0 | 25.2 / 14.1 | 10.3 / 10.4 | 10.1 / 10.5 | 7.8 / 2.5 | **6.5** / 2.5 | 4.0 / 2.6 | 11.4 / 11.5 | 13,042 |
| Human/Mouse ATM (`NP_000042.3`, `NP_031525.2`) | 3,056/3,066 | 83.84% | 230.1 / 108.9 | 33.9 / 19.7 | 14.3 / 12.5 | 14.0 / 12.5 | 10.4 / 2.5 | **8.1** / 2.7 | 4.7 / 2.5 | 15.2 / 14.6 | 13,333 |

Across all five proteins, Nucleics Align beats `parasail` AVX2 on 1 thread by 1.3–2.4x in `full` mode (2.0–14.3 ms vs 2.5–33.9 ms) and by 1.3–3.3x in `--linear` mode (2.0–10.4 ms, using 3.3–7.9x less RAM on the four large proteins); on 8 threads `--linear` reaches 6.5 ms on MTOR (3.9x faster than `parasail`) and 8.1 ms on ATM (4.2x faster than `parasail`, 28.4x faster and 40x smaller than original `seq-align`).

---

## 4. Discussion

Exact pairwise alignment has often forced a choice between general-purpose DP libraries that are predictable but slow and memory-hungry on long sequences, and output-sensitive wavefront aligners that are fast on near-identical inputs but degrade sharply on divergent ones. Nucleics Align combines compact traceback storage, difference-encoded striped SIMD with prefix-scan carry propagation, formally certified strip banding, Hirschberg/Myers-Miller linear space, and $O(s)$-space bidirectional WFA inside one portable C11 codebase, giving predictable performance from kilobase proteins and viral genomes up to megabase DNA pairs.

### Novelty and engineering contributions

Nucleics Align builds on established algorithmic foundations [@gotoh1982improved; @hirschberg1975linear; @myers1988optimal; @farrar2007striped; @suzuki2018introducing; @marcosola2021fast; @marcosola2023optimal; @myers1999fast; @grootkoerkamp2024apa2]. Its primary contributions are engineering and systems-level: (1) a closed-band exactness certificate (Section 2.4) that excludes tied optimal paths outside the band, combined with strip-level upper-bound pruning and automatic strip-kept memory fallback so full-matrix tie-breaking is preserved byte-for-byte even when a band would otherwise exceed RAM; and (2) the unified integration of all four exact speed-up families — with portable multi-versioned SIMD, cache-line-padded single-pair multi-threading, and full general substitution-matrix support — into a zero-external-alignment-dependency C11 library and CLI suite.

Nucleics Align is not a replacement for seed-and-extend read mappers such as `minimap2` [@li2018minimap2] or database search tools such as `BLAST+` [@camacho2009blastplus], nor for `edlib` [@sosic2017edlib] when only unit-cost edit distance is needed (as confirmed by Table 4). It targets workflows that require exact affine-gap global, overlap, or local alignment — including custom substitution matrices, wildcards, and fixed tie-breaking — from both the command line and C, Python, or Perl. Recent advances such as `Singletrack` [@lopezvillellas2026singletrack] inside `WFA2-lib` reduce affine backtrace storage 3x in its demonstrated configuration and remain complementary for wavefront-dominated workloads.

### Limitations and trade-offs

1. **Spin-waiting CPU overhead:** Multi-threaded single-pair fills synchronise column blocks by spinning (even after v1.2.0's 64-byte padding and 16-row batched publication), so total CPU time on 8 threads can still exceed single-threaded comparators on full-matrix fills (e.g., 3.14 s CPU vs 2.34 s for `ksw2` on 50 kb unrelated `full`, despite 414 ms vs 2,355 ms wall time). By contrast, `--linear` on 8 threads uses *less* total CPU time than `ksw2` (0.79–1.17 s vs 2.33–2.34 s), and `--threads 1` avoids spin overhead altogether on busy multi-job servers.
2. **Local alignment remains $O(mn)$ in time:** Because `smith_waterman --maxhits 1` uses neither banding nor WFA, long similar local pairs remain quadratic in time (e.g., 9.4 s and 626 MiB on the 197 kb Mpox pair vs 8.1 ms for global `--wfa`).
3. **WFA crossover and scoring scope:** `--wfa` applies only when scoring reduces to one match score, one mismatch score, and standard affine gap penalties. Although v1.2.0 BiWFA removes the retained-wavefront memory ceiling and aligns 50 kb 1% pairs in 47.6 ms / 3.3 MiB, on divergent ($\ge 10\%$) or unrelated pairs its work budget is exhausted and it falls back to the matrix aligner; it therefore remains opt-in.
4. **Memory vs tie-breaking on large pairs:** Full compact alignment on 580 kb *M. genitalium* still uses 4.33 GiB, and strip-kept mode on 1 Mb 1% pairs incurs refill overhead (9.23 s vs 7.19 s for `--linear` in Table 3). `--linear` and `--wfa` avoid this memory cost while returning the identical optimal score, but may select a different equal-scoring alignment when ties exist.
5. **Single-piece affine gaps, synthetic indel model, and platform dispatch:** Nucleics Align implements Gotoh's single-piece affine gap model rather than the dual-affine convex gap cost used by `minimap2`/`ksw2` [@li2018minimap2]. Synthetic benchmarks in Tables 2–4 model single-base indels rather than long block indels (which are represented in the real viral and protein pairs in Tables 1 and 5). Automatic AVX2/SSE4.1/SSE2 multi-versioning applies on x86-64 Linux; macOS and non-Linux x86-64 builds require `NATIVE=1` for AVX2/AVX-512.

### Future work

Two planned extensions build directly on the current architecture: combining Myers' bit-vector passes [@myers1999fast] with Hirschberg's divide-and-conquer [@hirschberg1975linear] so unit-cost `--linear` gains 64-cell-per-word throughput for full alignments, and exposing `--wfa` and `--linearabove` in the Python wrapper.

---

## 5. Availability, acknowledgements, and declarations

- **Code and data availability:** Nucleics Align is released under the BSD 3-Clause licence at https://github.com/DanielTillett/nucleics-align (citation metadata in `CITATION.cff`; a versioned release tag and Zenodo archive DOI will be minted when this manuscript is posted). All real biological sequences in Tables 1 and 5 are publicly available from NCBI RefSeq/GenBank under the listed accessions and can be reproduced with `benchmarks/real_sequences.sh` and `benchmarks/protein_bench.sh`; the ablation and 5-seed comparator suite are reproduced by `benchmarks/ablation.c`, `benchmarks/gen_pairs.py`, and `benchmarks/compare_harness.sh`, with raw logs archived in `benchmarks/results/pop-os-v1.2.0.txt` and `benchmarks/results/pop-os-v1.2.0-journal.txt` (reproduced by `benchmarks/compare_harness.sh`, `benchmarks/protein_bench.sh`, and `benchmarks/journal_bench.sh`).
- **Acknowledgements:** Nucleics Align is derived from Isaac Turner's `seq-align` (Unlicense) and bundles his `string_buffer` (CC0 1.0), `seq_file` (public domain), and `sort_r` (public domain) libraries. Substitution matrices in `scoring/` are from the NCBI BLAST matrix archive.
- **Funding:** [Placeholder: None / Funding source.]
- **Competing interests:** Daniel Tillett is the founder and director of Nucleics Pty Ltd; Nucleics Align is released under the open-source BSD 3-Clause licence with no commercial restrictions or patents.
- **Declaration of AI-assisted technologies:** Large language models (Claude, Codex, and Gemini) were used during software optimisation, code auditing, testing, and manuscript drafting. All algorithms, benchmark measurements, test verification, and scientific claims were reviewed and verified by the author, who takes full responsibility for this work.

---

## References

*(Rendered from `paper/references.bib` via BibTeX / Pandoc `--citeproc`.)*
