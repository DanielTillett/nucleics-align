# Choosing scores

[Documentation index](README.md)

A gap of n letters scores `gap_open + n * gap_extend` (both 0 or negative).
The tools' defaults are:

| Tool | Match | Mismatch | Gap open | Gap extend |
|---|---|---|---|---|
| `needleman_wunsch` | 1 | -2 | -4 | -1 |
| `smith_waterman` | 2 | -2 | -2 | -1 |

## DNA and RNA

NCBI's BLAST documentation suggests raising the size of the mismatch
penalty relative to the match reward for closer sequences:

> Many nucleotide searches use a simple scoring system that consists of a
> "reward" for a match and a "penalty" for a mismatch. The (absolute)
> reward/penalty ratio should be increased as one looks at more divergent
> sequences. A ratio of 0.33 (1/-3) is appropriate for sequences that are
> about 99% conserved; a ratio of 0.5 (1/-2) is best for sequences that are
> 95% conserved; a ratio of about one (1/-1) is best for sequences that are
> 75% conserved.

From the NCBI BLAST help pages
(<http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml#Reward-penalty>).

Gap penalties (open, extend) that NCBI offered alongside these:

| Open | Extend |
|---|---|
| -5 | -2 |
| -2 | -2 |
| -1 | -2 |
| 0 | -2 |
| -3 | -1 |
| -2 | -1 |
| -1 | -1 |

`--scoring DNA_HYBRIDIZATION` is a built-in alternative, and `--wildcard N
0` lets `N` match anything.

## Proteins

NCBI's suggested substitution matrix and gap penalties (open, extend) by
query length:

| Query length | Matrix | Gap penalties |
|---|---|---|
| < 35 | PAM30 | (9, 1) |
| 35-50 | PAM70 | (10, 1) |
| 50-85 | BLOSUM80 | (10, 1) |
| > 85 | BLOSUM62 | (10, 1) |

From <http://www.ncbi.nlm.nih.gov/blast/html/sub_matrix.html>. In this
project's sign convention, (10, 1) is `--gapopen -10 --gapextend -1`. The
four matrices are built in, with those gap penalties already set
(`--scoring BLOSUM62` gives -10, -1; `--gapopen` and `--gapextend` override
them), and many more are in
[`scoring/`](../scoring) for `--substitution_matrix`, from the
[NCBI BLAST matrix archive](https://ftp.ncbi.nlm.nih.gov/blast/matrices/).

The EBI's notes on choosing a matrix: <http://www.ebi.ac.uk/help/matrix.html>.

> [!NOTE]
> These NCBI and EBI pages are old and may have moved; the advice is
> quoted as it was.
