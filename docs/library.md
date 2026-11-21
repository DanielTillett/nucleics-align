# Using the C library

[Documentation index](README.md)

The aligners are a static C library, `src/libalign.a`, with headers in
`src/`. Complete C and C++ programs are in [`examples/`](../examples).

```c
#include <stdio.h>
#include "needleman_wunsch.h"

int main(void)
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);

  scoring_t scoring;
  scoring_init(&scoring, 1, -2, -4, -1,   // match, mismatch, gap open, extend
               false, false,              // free start, end gaps
               false, false, false,       // no gaps in a, in b, no mismatches
               false);                    // case sensitive

  needleman_wunsch_align("ACGTTGCA", "ACGTGCA", &scoring, nw, result);
  printf("%s\n%s\nscore: %i\n", result->result_a, result->result_b,
         result->score);

  alignment_free(result);
  needleman_wunsch_free(nw);
  return 0;
}
```

Build it from the repository root after `make`:

```console
$ cc -std=c11 -I src -I libs example.c -L src -L libs/string_buffer \
    -lalign -lstrbuf -lpthread -lz -o example
```

## Entry points

* `needleman_wunsch_align()` (`src/needleman_wunsch.h`): global alignment.
  Set `compact` on the aligner to keep one traceback byte per cell instead
  of the three score matrices (12x less memory, same alignment); the
  matrices are then not available, e.g. to `alignment_print_matrices()`.
  Compact fills of similar sequences skip
  cells that cannot be on a best path (`alignment_band_enabled`, default
  true), with the same alignment
* `smith_waterman_align()` then `smith_waterman_fetch()` repeatedly
  (`src/smith_waterman.h`): local alignments, best first.
  `smith_waterman_set_min_score()` skips hits below a score: the hits
  returned are the same, but lower-scoring cells are never indexed.
  `smith_waterman_set_compact()` suits callers that mostly want the best
  hit: it is found from traceback bytes, and the score matrices are only
  filled if a second hit is fetched. Above `smith_waterman_rect_min_cells`
  (default 2^24), and where no gap scores above 0, only the rectangle of
  the best hit gets traceback bytes, unless the hit spans the sequences;
  above `smith_waterman_strip_min_cells` (default 2^28) the hit is followed
  back from rows kept every so often instead, in O(n sqrt(m)) memory. The
  hits are the same either way
* `hirschberg_align()` (`src/hirschberg.h`): global alignment in O(n) memory;
  same scores as `needleman_wunsch_align()`, but returns false for options it
  does not support (a positive gap open, free start or end gaps,
  `no_gaps_in_*`, `no_mismatches`, or penalties large enough to risk
  overflow). Large sub-problems run their forward and reverse passes
  at the same time, each split into column blocks on several threads
  (tunable with `hirschberg_block_min_cells`, `hirschberg_block_min_cols`
  and `hirschberg_max_threads`, shared by the two passes). Passes skip cells
  that cannot be on a best
  path (`hirschberg_band_enabled`, default true), and sub-problems reuse rows
  saved by the passes above them (`hirschberg_reuse_depth`, default 2; 0
  turns it off), with the same result
* `hirschberg_score()` (`src/hirschberg.h`): only the global alignment score,
  in O(n) memory; same options as `hirschberg_align()`
* `scoring_init()` and the `scoring_system_*()` functions
  (`src/alignment_scoring.h`): set up a scoring system;
  `scoring_make_case_sensitive()` switches one to case-sensitive
  comparison, keeping its scores for upper case too

## Settings

Settings that change speed but never a result (threads, SIMD, banding,
8- and 16-bit lanes, block sizes) are in an `alignment_options_t`
(`src/alignment.h`). Give an aligner its own with its `options` field
(Needleman-Wunsch) or `smith_waterman_set_options()`, and pass them to
`hirschberg_align_opts()`/`hirschberg_score_opts()`; threads with different
options can then align at the same time. Without options, the process-wide
defaults apply: `alignment_options_default()` reads them from the
`alignment_*`, `hirschberg_*`, `smith_waterman_rect_min_cells` and
`smith_waterman_strip_min_cells`
variables, which keep their old meaning (e.g. `alignment_striped_enabled =
false` switches SIMD off everywhere).

> [!IMPORTANT]
> Two settings can change a result, and both are off by default. With
> `linear_fallback` (`alignment_linear_fallback`), a Needleman-Wunsch
> alignment too big for the score matrices or traceback bytes (too long, or
> out of memory) is made in linear space instead, where that supports the
> scoring; with `linear_above` (`alignment_linear_above`, a length), so is
> any pair whose longer sequence has at least that many letters. The score
> is the same, but linear space may choose another of several equal-scoring
> alignments. `needleman_wunsch --linearfallback` and `--linearabove <n>`
> turn them on.

## Equal-scoring alignments

Where several alignments score the same, which one is returned is fixed
and the same on every path of the full-matrix aligners (score matrices,
traceback bytes, bands, difference-encoded lanes, and bands or
Smith-Waterman hits kept in strips, which reproduce the same bytes):
* The traceback runs from the end of the alignment back to its start.
  Needleman-Wunsch starts in the matrix of the last cell that holds the
  best score, preferring GAP_A (a gap in the first sequence, from the cell
  above), then GAP_B (a gap in the second, from the cell on the left), then
  MATCH. Smith-Waterman starts in MATCH at the hit's last cell.
* At each step it moves to the first predecessor, in the same order GAP_A,
  GAP_B, MATCH, whose score plus the step's score gives the current score.
  So in GAP_A an extension wins over an opening, but in GAP_B an opening
  from GAP_A wins over an extension; from MATCH the diagonal cell's gap
  states win over its MATCH, which tends to put a gap as late in the
  alignment as it can go (`AAA` against `AA` gives `AA-`, `ACCA` against
  `ACA` gives `AC-A`). The `nw_ties` and `nw_free_gaps` command line
  fixtures pin such cases
* Smith-Waterman hits come out by score, highest first, then by lowest
  position in the first sequence, then in the second. A hit stops at a
  cell scoring 0, and one that reaches a cell an earlier hit used is
  dropped
* `--linear` (and `hirschberg_align()`) gives the same score but may pick
  another of the equal-scoring alignments: it splits each sub-problem at
  the first best crossing of its middle row (a crossing in MATCH or a gap
  before one inside a deletion), and in a one-row sub-problem puts the
  first sequence's letter at the first best column, or a deletion next to
  whichever neighbouring deletion it joins

## Errors

The aligners above end the program on an error. The `*_checked()` versions
(`needleman_wunsch_align_checked()`, `smith_waterman_align_checked()`,
`smith_waterman_fetch_checked()`, `hirschberg_align_checked()`,
`hirschberg_score_checked()`) return it instead, as an `aln_status_t`, and
fill an `aln_error_t` if you pass one (`alignment_fprint_error()` prints
the message the aligners would have):

| Status | Cause |
|---|---|
| `ALN_OK` | Success |
| `ALN_ERR_SCORING` | A character pair has no score (`err->a`, `err->b` name it) |
| `ALN_ERR_TOO_LONG` | The sequences are too long to address the matrices |
| `ALN_ERR_OVERFLOW` | Scores could pass `INT_MAX` |
| `ALN_ERR_UNSUPPORTED` | Linear space refuses these options, or penalties this large |
| `ALN_ERR_NOMEM` | Out of memory, while aligning or when fetching a second Smith-Waterman hit |

Errors in the input are found before any work. Running out of memory is
found as it happens, on any of the alignment's threads: what was allocated
is freed, there is no alignment or hit to read, and the aligner can be used
again. A second Smith-Waterman hit can run out of memory too, as it fills
the score matrices if a compact alignment has not. Creating an aligner or a
result (`needleman_wunsch_new()`, `alignment_create()` and so on) still ends
the program if out of memory. Scores should lie within plus or minus
`SCORE_ABS_MAX` (`INT_MAX / 4`); the command line tools and the scoring file
loaders refuse others.

## Printing

`alignment_colour_fprint_against()` and `alignment_fprint_spacer()` print to
a given stream instead of stdout.
