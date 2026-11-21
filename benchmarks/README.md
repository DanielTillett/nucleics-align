# Benchmarks

Scripts to reproduce the timings in [docs/performance.md](../docs/performance.md)
for Nucleics Align itself. Comparisons with other aligners (parasail, ksw2,
WFA2, SSW) need those tools built separately and are not included.

| Script | Does |
|---|---|
| [`make_pairs.sh`](make_pairs.sh) `<length> <divergence> <seed>` | A pair of random DNA sequences, the same for the same arguments. Divergence `1` has 1% substitutions and 0.2% insertions and deletions, `10` has 10% and 2%, and `rand` gives two unrelated sequences |
| [`run.sh`](run.sh) `[lengths]` | Times `needleman_wunsch` (full matrix, `--linear`, `--scoreonly`) and `smith_waterman --maxhits 1` on those pairs (seed 11), best of three, and prints a tab-separated table of wall time, CPU time, peak memory and score. Needs GNU time (`/usr/bin/time -f`) |
| [`percore.c`](percore.c) | Library speed on one thread and on the default threads for one pair, as cells per CPU second. Build instructions are at the top of the file |

From the repository root after `make`:

```console
$ sh benchmarks/run.sh "5000 20000"
$ sh benchmarks/make_pairs.sh 20000 rand 11 > pair.txt
$ cc -O2 -I src -I libs benchmarks/percore.c -L src -lalign -lpthread -o percore
$ ./percore pair.txt
```

> [!NOTE]
> Timings depend on the CPU, the number of cores and the build (`make
> NATIVE=1` tunes for this CPU, e.g. 16-lane AVX-512 vectors). The published
> numbers come from an AMD Ryzen AI Max+ 395 (16 cores, 32 threads) under
> WSL, and an AMD Ryzen 7 5700G (8 cores) under Linux for the megabase pairs.
