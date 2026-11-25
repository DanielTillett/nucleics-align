# Using the C library

[Documentation index](README.md)

The public API is one header, `include/nucleics_align.h`, with the `na_*` functions. It is C11 and works from C++ unchanged. After `make`, link against `src/libalign.a` and pthreads (the library needs no zlib):

```sh
cc -O3 myprog.c -I include -L src -lalign -lpthread
c++ -O3 myprog.cpp -I include -L src -lalign -lpthread
```

## Minimal example

```c
#include <stdio.h>
#include <string.h>
#include "nucleics_align.h"

int main(void) {
  na_scoring_t *sc = na_scoring_new(1, -2, -4, -1); // match, mismatch, gap open, gap extend
  na_aligner_t *al = na_aligner_new();
  na_alignment_t *res = na_alignment_new();
  if(!sc || !al || !res) return 1;

  const char *a = "GATTACAGATTACA", *b = "GATCACAGTTACA";
  if(na_global_align(al, a, strlen(a), b, strlen(b), sc, res) != NA_OK) {
    char msg[128];
    na_error_format(na_aligner_last_error(al), msg, sizeof(msg));
    fprintf(stderr, "%s\n", msg);
    return 1;
  }
  printf("%s\n%s\nscore: %d\n", na_alignment_a(res), na_alignment_b(res),
         na_alignment_score(res));

  na_alignment_free(res);
  na_aligner_free(al);
  na_scoring_free(sc);
  return 0;
}
```

Complete C and C++ programs live in [`examples/`](../examples/README.md).

## The three objects

| Object | Made by | Holds | Threads |
|---|---|---|---|
| `na_scoring_t` | `na_scoring_new()`, `na_scoring_new_preset()` | match/mismatch, gap penalties, pairs, wildcards, flags | share one read-only across threads once set up |
| `na_aligner_t` | `na_aligner_new()` | an aligner's working memory, settings and last error | one per thread |
| `na_alignment_t` | `na_alignment_new()` | one result: the two aligned strings, score and span | one per thread, reused across calls |

All three are opaque; free them with the matching `na_*_free()`, which accept NULL. Reuse one aligner and one result buffer for a loop of pairs: their buffers grow on demand and stay allocated.

## Scoring

`na_scoring_new(match, mismatch, gap_open, gap_extend)`: a gap of k letters scores `gap_open + k * gap_extend`. `na_scoring_new_preset(name)` takes `"PAM30"`, `"PAM70"`, `"BLOSUM80"`, `"BLOSUM62"`, `"DNA_HYBRIDIZATION"` (the list from `na_scoring_preset_names()`), and also `"DEFAULT"` or `"DNA_DEFAULT"` for the tools' default DNA scoring; it returns NULL for any other name. Scores must lie within `-INT_MAX/4 .. INT_MAX/4` or the setters return `NA_ERR_INVALID_ARG`.

- `na_scoring_set_flags(sc, flags)` **replaces** every flag with the set given (so pass the whole set each time). `NA_CASE_SENSITIVE` is one-way: once set, a call without it is refused with `NA_ERR_INVALID_ARG` and changes nothing. Pairs and wildcards added before it are copied to every mix of cases; those added after are exact.
- `na_scoring_set_pair(sc, a, b, score)` and `na_scoring_set_wildcard(sc, c, score)` add scores for a character pair, or a character that matches everything. Two wildcards against each other score the lower of the two.
- `na_scoring_parse_matrix()` and `na_scoring_parse_pairs()` read the file formats in [`scoring/`](../scoring) from memory; on error nothing changes and `err->line` says where.
- `na_scoring_lookup()` answers `NA_ERR_SCORING_PAIR` for a pair the scoring does not cover.

Flags for global alignment: `NA_FREE_START_GAP` and `NA_FREE_END_GAP` make a gap free while it runs along the first row or column (before any letter of the other sequence) or the last (after all of them). `NA_NO_GAPS_A`, `NA_NO_GAPS_B`, `NA_NO_GAPS_LONGEST` and `NA_NO_MISMATCHES` restrict the alignment; `NA_NO_MISMATCHES` with any `NA_NO_GAPS_*` is refused by `na_global_align()` (such a pair may have no alignment at all).

## Global alignment

`na_global_align(al, a, len_a, b, len_b, sc, out)` writes the best global alignment of `a` and `b` into `out`. Sequences are byte strings of the given lengths: empty sequences and embedded NUL bytes are fine. On any error `out` is cleared (length 0, score 0) and the status is also kept in `na_aligner_last_error(al)`.

`na_aligner_set_method()` picks how:

| Method | What runs | Settings it honours |
|---|---|---|
| `NA_METHOD_AUTO` (default) | exact DP with one traceback byte per cell, banded and SIMD where it applies; the documented tie-breaking | `linear_above`, `linear_fallback`, `threads` |
| `NA_METHOD_FULL` | the three full score matrices (12 bytes per cell); needed for `na_debug_fprint_matrices()` | `threads` |
| `NA_METHOD_LINEAR` | Hirschberg / Myers-Miller in linear memory; refuses the free-gap, no-gap and no-mismatch flags and a positive gap open (`NA_ERR_UNSUPPORTED`) | `threads` |
| `NA_METHOD_WFA` | bidirectional wavefront alignment within `wfa_budget` diagonal updates, otherwise as AUTO | `wfa_budget`, then as AUTO |

Every method returns the same optimal score. AUTO and FULL return the same alignment, chosen by the tie-breaking rules below; LINEAR and WFA, and AUTO once `na_aligner_set_linear_above()` or `na_aligner_set_linear_fallback()` has sent a pair to linear memory, may return another alignment with that score.

`na_global_score(al, a, len_a, b, len_b, sc, &score)` gives the score alone in linear memory (unit-cost scoring uses Myers' bit-vector algorithm). It ignores the method; for flags linear memory refuses, it aligns with the compact matrices and reports their score.

## Local alignment

`na_local_align(al, a, len_a, b, len_b, sc)` finds the local alignments; `na_local_next(al, out)` then fetches them best first, returning `NA_OK` with a hit in `out`, `NA_DONE` after the last, or an error. Hits never share a cell with an earlier one, so each is a distinct alignment. Only hits scoring at least `na_aligner_set_min_local_score()` are fetched: the default is 1 (every positive hit), and values below 1 are clamped to 1. The command line tools instead default to `match * max(0.2 * shorter length, 2)`, which keeps the number of hits on long sequences sensible; set it yourself for long pairs.

`a`, `b` and `sc` must stay valid and unchanged until the next `na_local_align()`, `na_global_align()`, `na_global_score()` or `na_aligner_free()`: any of those ends the iteration (`na_local_next()` then returns `NA_DONE`). The first hit of a long pair is found with little memory; the second fills the score matrices (12 bytes per cell).

`na_aligner_set_fast_local(al, true)` reconstructs the first hit of a long pair with wavefront or linear-memory alignment, verified against the exact score. The score is the same, but among equal-scoring paths another may be chosen, and because the cells it marks as used differ, later hits can differ too.

## Results

`na_alignment_a()` and `na_alignment_b()` return the aligned strings with `'-'` for gaps, `na_alignment_length()` their length, `na_alignment_score()` the score, and `na_alignment_span()` the 0-based start and the letters consumed in each sequence (for a global alignment `0, len_a, 0, len_b`). The strings belong to the result object: they are valid until the next call that writes to it or frees it. Copy them if you need them longer. Passing NULL to these accessors gives `""`, 0 or no span.

Do not pass a result's own strings back in as input to a call that writes to the same result.

`na_score_alignment(res_a, res_b, len, sc, &score)` scores a pair of aligned strings as the aligner would, including the free start/end gap flags; a column of two gaps is `NA_ERR_INVALID_ARG`, an unscored pair `NA_ERR_SCORING_PAIR` (so is a mismatch under `NA_NO_MISMATCHES`).

## Settings

| Setter | Meaning | Default |
|---|---|---|
| `na_aligner_set_threads(al, n)` | threads within one alignment (0: half the CPUs; at most 4096) | 0 |
| `na_aligner_set_method(al, m)` | see above | `NA_METHOD_AUTO` |
| `na_aligner_set_linear_above(al, n)` | under AUTO, align pairs whose longer sequence has at least n letters in linear memory (0: never) | 0 |
| `na_aligner_set_linear_fallback(al, on)` | under AUTO, align in linear memory when the matrices would not fit or ran out of memory | off |
| `na_aligner_set_wfa_budget(al, n)` | diagonal updates `NA_METHOD_WFA` may spend before giving the pair to AUTO (0: `min(2^28, 256 * (len_a + len_b))`) | 0 |
| `na_aligner_set_min_local_score(al, s)` | lowest local hit score fetched (clamped to at least 1) | 1 |
| `na_aligner_set_fast_local(al, on)` | see local alignment | off |

Settings apply to the calls after them.

## Errors

Every fallible call returns an `na_status_t`; the aligner's calls also store it, with any detail, in `na_aligner_last_error(al)`, which stays valid until the aligner's next call. `na_strerror()` names a status and `na_error_format()` writes the detail (the offending character pair, or the line of a parse error) into a buffer.

| Status | Meaning |
|---|---|
| `NA_OK` | Success |
| `NA_DONE` | `na_local_next()`: no more hits |
| `NA_ERR_INVALID_ARG` | NULL where a pointer is needed, an unknown flag, a score out of range, or a change the scoring refuses |
| `NA_ERR_SCORING_PAIR` | A character pair in the sequences has no score (`err->a`, `err->b`) |
| `NA_ERR_TOO_LONG` | The score matrices would not fit in memory addresses |
| `NA_ERR_OVERFLOW` | Scores could pass `INT_MAX`, or reach the matrices' floor |
| `NA_ERR_UNSUPPORTED` | The method refuses these flags or penalties (linear memory some; Smith-Waterman refuses gaps that score above 0) |
| `NA_ERR_NOMEM` | Out of memory; the aligner stays usable |
| `NA_ERR_PARSE` | A scoring file's text is malformed (`err->line`) |

`na_debug_fprint_matrices(al, fp)` prints the three score matrices of the last `NA_METHOD_FULL` alignment, for checking a scoring rule by hand; it reads the sequences and scoring given to that call, which must still be valid. It prints nothing after any other method.

## Equal-scoring alignments

When two alignments share the best score, the choice under `NA_METHOD_AUTO` and `NA_METHOD_FULL` (and the command line tools) is deterministic and identical across every full-matrix path (score matrices, traceback bytes, exact bands, strip-kept bands and difference-encoded lanes):

<p align="center">
  <img src="figures/tie-breaking.svg" width="760" alt="Gotoh traceback tie-breaking: current state fixes the spatial step; GAP_A before GAP_B before MATCH picks the predecessor state, shifting gaps right">
</p>

- **Starting cell:** Needleman-Wunsch starts at the bottom-right cell `(len_a, len_b)` in the best matrix in the order `GAP_A`, `GAP_B`, `MATCH`. Smith-Waterman starts in `MATCH` at the highest-scoring cell (breaking ties by lowest end coordinate `x` in `seq_a`, then `y` in `seq_b`).
- **Each step back:** The current matrix fixes the spatial step (`GAP_A` moves to `(x, y-1)`, `GAP_B` to `(x-1, y)`, `MATCH` to `(x-1, y-1)`). At that predecessor cell, traceback picks the first matrix in the order `GAP_A`, `GAP_B`, `MATCH` whose score plus the step's penalty equals the current score. That tends to push gaps toward the right end of runs (for example, `AAA` vs `AA` gives `AA-`, and `ACCA` vs `ACA` gives `AC-A`).
- **Local hits:** Candidate hits are sorted by score (highest first), then by the hit's end cell (lowest `x`, then lowest `y`). The returned start positions are computed during traceback. A hit stops at score `0`, and any candidate path that touches a cell already used by an earlier hit is discarded.
- **Linear memory and WFA:** `NA_METHOD_LINEAR` and `NA_METHOD_WFA` always find the same optimal score, but use their own midpoint or wavefront tie-breaking.

## Internal API (`src/lib`)

The `na_*` functions wrap the aligners in `src/lib`, which the command line tools and the test suite use directly. They are not a stable interface, but they expose a few things the public API does not (the raw `aligner_t` matrices, direct wavefront and Myers calls, the tuning globals). Compile against them with `-I include -I src/lib`.

| Function | Header | What it does |
|---|---|---|
| `needleman_wunsch_align_checked()` | `needleman_wunsch.h` | Global alignment (`nw->compact = true` for one traceback byte per cell) |
| `smith_waterman_align_checked()`, `smith_waterman_fetch_checked()` | `smith_waterman.h` | Local alignments, best first |
| `hirschberg_align_checked()`, `hirschberg_score_checked()` | `hirschberg.h` | Linear-memory global alignment and score |
| `alignment_wfa_align()`, `alignment_wfa_score()` | `alignment.h` | Wavefront alignment and score, given a `sub_table_t` built by `sub_table_build()` |
| `alignment_myers_score()` | `alignment.h` | Bit-vector unit-cost edit distance |
| `scoring_init()`, `scoring_system_by_name()`, `scoring_add_mutation()`, `scoring_add_wildcard()` | `alignment_scoring.h` | Scoring setup |

Every aligner has an unchecked form (`needleman_wunsch_align()`, `smith_waterman_align()`, `hirschberg_align()`) that prints the error and calls `exit()` instead of returning a status, and a `*_checked()` form returning an `aln_status_t` (`ALN_OK`, `ALN_ERR_SCORING`, `ALN_ERR_TOO_LONG`, `ALN_ERR_OVERFLOW`, `ALN_ERR_UNSUPPORTED`, `ALN_ERR_NOMEM`, which `na_status_t` mirrors) with an optional `aln_error_t` that `alignment_fprint_error()` formats. `smith_waterman_fetch_checked()` returns 1 for a hit and 0 both at the end and on error: check `err->status`.

Speed and threading settings that do not change the result live in `alignment_options_t` (`alignment_options_default()`, then `nw->options = &opts`, `smith_waterman_set_options(sw, &opts)`, or the `hirschberg_*_opts()` calls). `linear_fallback`, `linear_above` and `wfa_align` are the three that may change which equal-scoring alignment comes back.

`alignment_wfa_align()` takes a budget in diagonal updates and returns `ALN_ERR_UNSUPPORTED` when the scoring is not one match score, one mismatch score and gap penalties, or the budget runs out, and `ALN_ERR_NOMEM` when memory does. Direct callers get no fallback; `needleman_wunsch_align_checked()` with `opts.wfa_align` falls back to the matrices itself.

`alignment_colour_fprint_against()` and `alignment_fprint_spacer()` (`alignment.h`) write coloured alignments and `|`/`*` marker lines to a `FILE *`.
