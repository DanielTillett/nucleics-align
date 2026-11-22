# Benchmarks

Scripts for reproducing the timings in [docs/performance.md](../docs/performance.md) and the manuscript ([paper/manuscript.md](../paper/manuscript.md)). Build the tools first with `make` in the repository root; GNU `time` (`/usr/bin/time`) is used for CPU time and peak RSS.

- **`benchmarks/run.sh [lengths]`** — runs `needleman_wunsch` (`full`, `--linear`, `--scoreonly`) and `smith_waterman --maxhits 1` on simulated DNA pairs from `benchmarks/make_pairs.sh` (best of 3 runs).
- **`benchmarks/real_sequences.sh [output_dir] [threads]`** — downloads the seven public NCBI RefSeq/GenBank pairs in Table 1, logs SHA-256 digests, times all modes with pinned threads (default 8), verifies identical global scores across modes, and computes full-alignment identity.
- **`benchmarks/gen_pairs.py <length> <div|rand> <seed> [out.fa]`** — deterministic, platform-independent Python 3 DNA pair generator (MT19937 PRNG).
- **`benchmarks/ablation.c`** — C11 library-level ablation driver timing each optimisation layer (scalar 3-matrix, compact 1-byte, striped 32-bit, difference-encoded 8-bit, certified banding, multi-threading, linear space, and WFA) and asserting identical optimal scores.
- **`benchmarks/compare_harness.sh [work_dir] [threads]`** — builds pinned revisions of original `seq-align` (`dc41988`), `ksw2` (`289609b`), `WFA2-lib` (`v2.3.5` / `cf3eb92`), and `SSW` (`a66636b`), then runs the 5-seed ablation and comparator suite (20 kb and 50 kb global, 20 kb local) with median and min–max wall time, CPU time, peak RSS, and global score verification.
