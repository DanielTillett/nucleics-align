# Examples

C and C++ programs that call the Nucleics Align library, as a starting
point for your own:

| File | Shows |
|---|---|
| [`nw_example.c`](nw_example.c) | Global alignment (Needleman-Wunsch) with free start and end gaps and custom substitution scores |
| [`sw_example.c`](sw_example.c) | Local alignments (Smith-Waterman) with custom substitution scores, fetching hits best first |
| [`nw_example.cpp`](nw_example.cpp) | Global alignment from C++ with no gaps or mismatches allowed, so it needs sequences that match exactly |
| [`sw_example.cpp`](sw_example.cpp) | Local alignment from C++ with no gaps or mismatches allowed: the exact common substrings |

Build them after building the library (`make` at the top level builds them
too):

```console
$ make all
$ ./nw_example ACGTTGCA ACGTGCA
$ ./sw_example ttacgtacgtaa acgtacgt
```

The [library guide](../docs/library.md) describes the API, settings and
error handling. Python and Perl examples are in [`python/`](../python) and
[`perl/`](../perl).

Originally written by Isaac Turner <turner.isaac@gmail.com>.
