# Examples

C and C++ programs that call the Nucleics Align library, as a starting
point for your own:

| File | Shows |
|---|---|
| [`nw_example.c`](nw_example.c), [`nw_example.cpp`](nw_example.cpp) | Global alignment (Needleman-Wunsch) with a custom substitution score |
| [`sw_example.c`](sw_example.c), [`sw_example.cpp`](sw_example.cpp) | Local alignments (Smith-Waterman), fetching hits best first |

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
