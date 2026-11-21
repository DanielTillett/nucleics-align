# Nucleics Align documentation

If you are new here, start with the [README](../README.md) for building and quick examples, then jump into whichever guide matches what you are doing:

| Guide | What you will find |
|---|---|
| [Command line tools](command-line.md) | `needleman_wunsch`, `smith_waterman` and `lcs`: worked examples, batching, and every flag |
| [Choosing scores](scoring.md) | How affine gap penalties work, and sensible DNA/RNA and protein defaults |
| [Using the C library](library.md) | Function signatures, `alignment_options_t`, tie-breaking rules and error handling |
| [Performance](performance.md) | Benchmarks from 5 kb to 1 Mb (including `--wfa`), comparison with other aligners, and how each speed-up works |
| [Algorithms and references](algorithms.md) | The papers behind the algorithms, with DOIs |
| [Examples](../examples/README.md) | Small, complete C and C++ programs linking `libalign.a` |
| [Benchmarks](../benchmarks/README.md) | Scripts to reproduce the timings on your own machine |

For what has changed since Isaac Turner's `seq-align`, see [CHANGELOG.md](../CHANGELOG.md).
