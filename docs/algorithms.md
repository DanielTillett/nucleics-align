# Algorithms and references

[Documentation index](README.md)

## Foundations

| Method | Reference |
|---|---|
| Global alignment | Needleman & Wunsch 1970, *J Mol Biol* 48:443-453 ([doi:10.1016/0022-2836(70)90057-4](https://doi.org/10.1016/0022-2836(70)90057-4)) |
| Local alignment | Smith & Waterman 1981, *J Mol Biol* 147:195-197 ([doi:10.1016/0022-2836(81)90087-5](https://doi.org/10.1016/0022-2836(81)90087-5)) |
| Affine gaps with three matrices (match, gap in a, gap in b) | Gotoh 1982, *J Mol Biol* 162:705-708 ([doi:10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836(82)90398-9)) |
| Linear-space alignment | Hirschberg 1975, *Commun ACM* 18:341-343 ([doi:10.1145/360825.360861](https://doi.org/10.1145/360825.360861)) |
| Linear space with affine gaps | Myers & Miller 1988, *Comput Appl Biosci* 4:11-17 ([doi:10.1093/bioinformatics/4.1.11](https://doi.org/10.1093/bioinformatics/4.1.11)) |

## Speed-ups used here

| Method | Where | Reference |
|---|---|---|
| Striped SIMD fill | Full matrices and linear-space passes | Farrar 2007, *Bioinformatics* 23:156-161 ([doi:10.1093/bioinformatics/btl582](https://doi.org/10.1093/bioinformatics/btl582)) |
| Lazy F carries found by a scan over the lanes | Striped fills | The idea of Block Aligner's prefix scans: Liu & Steinegger 2023, *Bioinformatics* 39:btad487 ([doi:10.1093/bioinformatics/btad487](https://doi.org/10.1093/bioinformatics/btad487)) |
| Difference-encoded lanes | Needleman-Wunsch traceback bytes | Suzuki & Kasahara 2018, *BMC Bioinformatics* 19(Suppl 1):45 ([doi:10.1186/s12859-018-2014-8](https://doi.org/10.1186/s12859-018-2014-8)), as ksw2 uses them |
| Wavefront alignment (WFA), score only, to size a band | Banded global alignment | Marco-Sola et al. 2021, *Bioinformatics* 37:456-463 ([doi:10.1093/bioinformatics/btaa777](https://doi.org/10.1093/bioinformatics/btaa777)); match scores converted to penalties as Eizenga & Paten describe (2022, bioRxiv, [doi:10.1101/2022.01.12.476087](https://doi.org/10.1101/2022.01.12.476087), a preprint) |
| Computational volumes, adapted from edit distance to affine gaps | Narrowing each strip of a band | A\*PA2: Groot Koerkamp 2024, WABI, *LIPIcs* 312:17 ([doi:10.4230/LIPIcs.WABI.2024.17](https://doi.org/10.4230/LIPIcs.WABI.2024.17)) |

The code is this project's own, written from the papers.

`--nogapsinlongest` (no gaps within the longer of the two sequences, or the
first if they are the same length) is this project's addition.

## To investigate

- Myers' bit-vector algorithm (1999, *J ACM* 46:395-415,
  [doi:10.1145/316542.316550](https://doi.org/10.1145/316542.316550)). It only computes
  unit-cost edit distance, so it would be a separate fast path rather than a
  replacement for the scored aligners
- Combining Myers' and Hirschberg's algorithms for linear-space
  edit-distance alignment
