Benchmarks
==========

Scripts to reproduce the timings in the main README for Nucleics Align
itself. Comparisons with other aligners (parasail, ksw2, WFA2, SSW) need
those tools built separately and are not included.

* `make_pairs.sh <length> <divergence> <seed>`: a pair of random DNA
  sequences, the same for the same arguments. Divergence `1` has 1%
  substitutions and 0.2% insertions and deletions, `10` has 10% and 2%,
  and `rand` gives two unrelated sequences.
* `run.sh [lengths]`: times `needleman_wunsch` (full matrix, `--linear`,
  `--scoreonly`) and `smith_waterman --maxhits 1` on those pairs (seed 11),
  best of three, and prints a tab-separated table with wall time, CPU time,
  peak memory and score. Needs GNU time (`/usr/bin/time -f`).
* `percore.c`: library speed on one thread and on the default threads for
  one pair, as cells per CPU second. Build instructions are at the top of
  the file.

Usage, from the repository root after `make`:

    $ sh benchmarks/run.sh "5000 20000"
    $ sh benchmarks/make_pairs.sh 20000 rand 11 > pair.txt
    $ cc -O2 -I src -I libs benchmarks/percore.c -L src -lalign -lpthread -o percore
    $ ./percore pair.txt

Timings depend on the CPU, the number of cores and the build (`make
NATIVE=1` tunes for this CPU, e.g. 16-lane AVX-512 vectors). The README's
numbers come from an AMD Ryzen AI Max+ 395 (16 cores, 32 threads) under WSL.
