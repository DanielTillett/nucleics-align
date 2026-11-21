# Algorithms and references

## Used here

| Method | Where | Reference |
|---|---|---|
| Affine gaps with three matrices (match, gap in a, gap in b) | Both aligners | Gotoh 1982, *J Mol Biol* 162:705-708 |
| Linear-space alignment with affine gaps | `needleman_wunsch --linear`, `hirschberg_align()` | Hirschberg's algorithm with the extension of Myers & Miller 1988, *Comput Appl Biosci* 4:11-17 |
| Striped SIMD fill | Full matrices and linear-space passes | Farrar 2007, *Bioinformatics* 23:156-161 |
| Lazy F carries found by a scan over the lanes | Striped fills | The idea of Block Aligner's prefix scans, Liu & Steinegger 2023, *Bioinformatics* 39:btad487 ([doi:10.1093/bioinformatics/btad487](https://doi.org/10.1093/bioinformatics/btad487)) |
| Difference-encoded lanes | Needleman-Wunsch traceback bytes | Suzuki & Kasahara 2018, *BMC Bioinformatics* 19:45, as ksw2 uses them |
| Wavefront alignment (WFA), score only, to size a band | Banded global alignment | Marco-Sola et al. 2021, *Bioinformatics* 37:456-463; match scores converted to penalties as Eizenga & Paten describe (2022, bioRxiv 10.1101/2022.01.12.476087, a preprint) |
| Computational volumes, adapted from edit distance to affine gaps | Narrowing each strip of a band | A\*PA2, Groot Koerkamp 2024, WABI, *LIPIcs* 312:17 |

The code is this project's own, written from the papers.

`--nogapsinlongest` (no gaps within the longer of the two sequences, or the
first if they are the same length) is this project's addition.

## To investigate

- Myers' bit-vector algorithm (1999, *J ACM* 46:395-415). It only computes
  unit-cost edit distance, so it would be a separate fast path rather than a
  replacement for the scored aligners
- Combining Myers' and Hirschberg's algorithms for linear-space
  edit-distance alignment
