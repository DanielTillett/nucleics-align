# Command line tools

`make` builds three tools in `bin/`:

| Tool | Does |
|---|---|
| `needleman_wunsch` | Optimal global alignment: whole sequences, end to end |
| `smith_waterman` | Optimal local alignments: the best-scoring regions, best first |
| `lcs` | Longest common substrings of a sequence with itself |

Each takes a pair of sequences on the command line, or reads pairs from
FASTA, FASTQ, SAM, BAM, plain (one sequence per line) or gzipped files, or
from standard input. Outputs below are from the current build.

- [Global alignment](#global-alignment)
- [Local alignment](#local-alignment)
- [Reading files](#reading-files)
- [Scoring](#scoring)
- [Longest common substrings](#longest-common-substrings)
- [All options](#all-options)

## Global alignment

```console
$ ./bin/needleman_wunsch CAGACGT CGATA
C-AGACGT
CGATA---
```

`--printscores` adds the score, and `--pretty` a line marking matches (`|`),
mismatches (`*`) and gaps (space):

```console
$ ./bin/needleman_wunsch --printscores --pretty CAGACGT CGATA
C-AGACGT
| |*|
CGATA---
score: -11
```

In linear memory, for long sequences (`--linear`), or the score alone
(`--scoreonly`, the fastest):

```console
$ ./bin/needleman_wunsch --linear --printscores CAGACGT CGATA
CAGACGT
C-GA-TA
score: -11

$ ./bin/needleman_wunsch --scoreonly CAGACGT CGATA
score: -11
```

> [!NOTE]
> `--linear` gives the same score as the default aligner but, where several
> alignments share the best score, it may return a different one, as above.
> It cannot be combined with `--freestartgap`, `--freeendgap`, `--nogaps..`,
> `--nomismatches`, `--printmatrices` or a positive `--gapopen`;
> `--scoreonly` uses the full matrices for those. `--linearfallback` uses
> linear memory only for pairs too big for the matrices.

`--printmatrices` prints the three dynamic programming matrices (match, gap
in the first sequence, gap in the second) before the alignment:

<details>
<summary>Example output</summary>

```console
$ ./bin/needleman_wunsch --printmatrices ACAGGT AAGGT
seq_a: ACAGGT
seq_b: AAGGT
match_scores:
  0:	  0	-2147483643	-2147483643	-2147483643	-2147483643	-2147483643	-2147483643
  1:	-2147483643	  1	 -7	 -5	 -9	-10	-11
  2:	-2147483643	 -4	 -1	 -3	 -7	 -8	 -9
  3:	-2147483643	 -8	 -6	 -3	 -2	 -6	-10
  4:	-2147483643	 -9	 -7	 -8	 -2	 -1	 -8
  5:	-2147483643	-10	 -8	 -9	-10	 -4	  0
gap_a_scores:
  0:	  0	-2147483643	-2147483643	-2147483643	-2147483643	-2147483643	-2147483643
  1:	 -5	-10	-11	-12	-13	-14	-15
  2:	 -6	 -4	 -9	-10	-11	-12	-13
  3:	 -7	 -5	 -6	 -8	-12	-13	-14
  4:	 -8	 -6	 -7	 -8	 -7	-11	-13
  5:	 -9	 -7	 -8	 -9	 -7	 -6	-11
gap_b_scores:
  0:	  0	 -5	 -6	 -7	 -8	 -9	-10
  1:	-2147483643	-10	 -4	 -5	 -6	 -7	 -8
  2:	-2147483643	-11	 -9	 -6	 -7	 -8	 -9
  3:	-2147483643	-12	-10	-11	 -8	 -7	 -8
  4:	-2147483643	-13	-11	-12	-13	 -7	 -6
  5:	-2147483643	-14	-12	-13	-14	-12	 -9
match: 1 mismatch: -2 gapopen: -4 gapexend: -1

ACAGGT
A-AGGT
```

</details>

## Local alignment

`smith_waterman` lists the hits best first, each with where it lies in each
sequence (0-based start and length):

```console
$ ./bin/smith_waterman --pretty \
    acgtttgacgtacgtagctagctgatcgtagctagctagctgatcgatcgtagctagc \
    ttttttacgtttgacgtacgtagcccccgctgatcgtagctagctagctgatcgatcgtaggggggg
== Alignment 0 lengths (58, 67):

hit 0.0 score: 94
  acgtttgacgtacgtagcta--gctgatcgtagctagctagctgatcgatcgtag  [pos: 0; len: 53]
  ||||||||||||||||||**  |||||||||||||||||||||||||||||||||
  acgtttgacgtacgtagcccccgctgatcgtagctagctagctgatcgatcgtag  [pos: 6; len: 55]

hit 0.1 score: 35
  tga-cgtacg-tagctagctgatcg  [pos: 5; len: 23]
  ||| |||| | ||||||||||||||
  tgatcgta-gctagctagctgatcg  [pos: 30; len: 24]

==
```

`--maxhits <n>` stops after n hits (`--maxhits 1` is much faster for the
best hit alone), `--minscore <s>` drops hits scoring below s, and
`--context <n>` shows n letters either side (marked `.`):

```console
$ ./bin/smith_waterman --maxhits 1 --context 4 --pretty \
    acgtttgacgtacgtagctagctgatcgtagctagctagctgatcgatcgtagctagc \
    ttttttacgtttgacgtacgtagcccccgctgatcgtagctagctagctgatcgatcgtaggggggg
== Alignment 0 lengths (58, 67):

hit 0.0 score: 94
      acgtttgacgtacgtagcta--gctgatcgtagctagctagctgatcgatcgtagctag  [pos: 0; len: 53]
      ||||||||||||||||||**  |||||||||||||||||||||||||||||||||....
  ttttacgtttgacgtacgtagcccccgctgatcgtagctagctagctgatcgatcgtaggggg  [pos: 6; len: 55]

==
```

## Reading files

`--file` reads two sequences at a time, here from a gzipped FASTA file
`dna.fa.gz` holding:

```text
>seqA
ACAATAGAC
>seqB
ACGAATAGAT
>seqC
ACGTGA
CAGAT
>seqD
GTGGACG
AGTA
```

```console
$ ./bin/needleman_wunsch --printscores --file dna.fa.gz
AC-AATAGAC
ACGAATAGAT
score: 1

ACGTGAC-AGAT
GTG-GACGAGTA
score: -12
```

`--file -` reads standard input, and `--files <f1> <f2>` takes one sequence
from each file. Pairs from a file are aligned several at a time, on
`--threads <n>` threads (default: half the CPUs), and printed in input
order, so the output is the same as one at a time:

```console
$ gzip -d -c dna.fa.gz | ./bin/needleman_wunsch --file -
AC-AATAGAC
ACGAATAGAT

ACGTGAC-AGAT
GTG-GACGAGTA
```

`--stdin` also reads standard input, but one pair at a time as it arrives;
`smith_waterman --stdin` then waits for `h` (next hit) or `a` (next
alignment) after each hit, as the [Perl module](../README.md#perl) uses it.

## Scoring

A gap of n letters scores `gapopen + n * gapextend`. The defaults are match
1, mismatch -2, gap open -4, gap extend -1 for `needleman_wunsch`, and match
2, mismatch -2, gap open -2, gap extend -1 for `smith_waterman`.

```console
$ ./bin/needleman_wunsch --match 1 --mismatch 0 --gapopen -10 --gapextend 0 \
    --printscores ACGTGCCCCACAGAT AGGTGGACGAGAT
ACGTGCCCCACAGAT
AGGTGGACGA--GAT
score: -1

$ ./bin/needleman_wunsch --scoring BLOSUM62 --printscores --pretty \
    MKTAYIAKQRQISFVKSHFSRQ MKTAHIAKQRISFVKSHFSRQ
MKTAYIAKQRQISFVKSHFSRQ
||||*||||| |||||||||||
MKTAHIAKQR-ISFVKSHFSRQ
score: 88
```

- `--scoring` picks a built-in system: `PAM30`, `PAM70`, `BLOSUM80`,
  `BLOSUM62` or `DNA_HYBRIDIZATION`
- `--substitution_matrix <file>` loads a matrix and `--substitution_pairs
  <file>` a list of pairs; see [`scoring/`](../scoring) for examples
- `--wildcard <w> <s>`: letter w matches every letter with score s (e.g.
  `--wildcard N 0`)
- `--case_sensitive` compares letters case sensitively (default: off)
- `--freestartgap`, `--freeendgap` (global only) make gaps at the ends free
- `--nogaps`, `--nogapsin1`, `--nogapsin2`, `--nogapsinlongest` forbid gaps
  within a sequence, and `--nomismatches` forbids mismatches. Used together
  for local alignment, they list the common substrings in order of length

[Choosing scores](scoring.md) has guidance on matrices and penalties.

## Longest common substrings

```console
$ ./bin/lcs ACGTACGTTACGACGT
ACGT [0,4]
ACGT [0,12]
TACG [3,8]
ACGT [4,12]
ACG [0,9]
ACG [9,12]
T [7,8]
T [8,15]
```

## All options

<details>
<summary><code>needleman_wunsch</code></summary>

```text
usage: ./bin/needleman_wunsch [OPTIONS] [seq1 seq2]
  Needleman-Wunsch optimal global alignment (maximises score).
  Takes a pair of sequences on the command line, or can read from a
  file and from sequence piped in.  Can read gzip files, FASTA and FASTQ.

  OPTIONS:
    --file <file>        Sequence file reading with gzip support - read two
                         sequences at a time and align them
    --files <f1> <f2>    Read one sequence from each file to align at one time
    --stdin              Read from STDIN (same as '--file -')

    --case_sensitive     Use case sensitive character comparison [default: off]

    --match <score>      [default: 1]
    --mismatch <score>   [default: -2]
    --gapopen <score>    [default: -4]
    --gapextend <score>  [default: -1]

    --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62|DNA_HYBRIDIZATION>
    --substitution_matrix <file>  see details for formatting
    --substitution_pairs <file>   see details for formatting

    --wildcard <w> <s>   Character <w> matches all characters with score <s>

    --threads <n>        Threads to use, within an alignment and for aligning
                         pairs from a file at once [default: half the CPUs]


    --freestartgap       No penalty for gap at start of alignment
    --freeendgap         No penalty for gap at end of alignment

    --printscores        Print optimal alignment scores
    --zam                A funky type of output
    --linear             Align in linear memory (Hirschberg, Myers-Miller)
    --scoreonly          Print only the score (fast, linear memory)
    --linearfallback     Align in linear memory where too big for the matrices
                         (same score, but may pick another equal-scoring
                         alignment)
    --printmatrices      Print dynamic programming matrices
    --printfasta         Print fasta header lines
    --pretty             Print with a descriptor line
    --colour             Print with colour

  Experimental Options:
    --nogapsin1          No gaps allowed within the first sequence
    --nogapsin2          No gaps allowed within the second sequence
    --nogaps             No gaps allowed in either sequence
    --nogapsinlongest    No gaps allowed within the longer sequence
    --nomismatches       No mismatches allowed (cannot be used with --nogaps..)

 DETAILS:
  * For help choosing scoring, see the README file.
  * Gap (of length N) penalty is: (open+N*extend)
  * To do alignment without affine gap penalty, set '--gapopen 0'.
  * Scoring files should be matrices, with entries separated by a single
    character or whitespace. See files in the 'scoring' directory for examples.

  https://github.com/DanielTillett/nucleics-align  (compiled: <date>)
```

</details>

<details>
<summary><code>smith_waterman</code></summary>

```text
usage: ./bin/smith_waterman [OPTIONS] [seq1 seq2]
  Smith-Waterman optimal local alignment (maximises score).
  Takes a pair of sequences on the command line, or can read from a
  file and from sequence piped in.  Can read gzip files, FASTA and FASTQ.

  OPTIONS:
    --file <file>        Sequence file reading with gzip support - read two
                         sequences at a time and align them
    --files <f1> <f2>    Read one sequence from each file to align at one time
    --stdin              Read from STDIN (same as '--file -')

    --case_sensitive     Use case sensitive character comparison [default: off]

    --match <score>      [default: 2]
    --mismatch <score>   [default: -2]
    --gapopen <score>    [default: -2]
    --gapextend <score>  [default: -1]

    --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62|DNA_HYBRIDIZATION>
    --substitution_matrix <file>  see details for formatting
    --substitution_pairs <file>   see details for formatting

    --wildcard <w> <s>   Character <w> matches all characters with score <s>

    --threads <n>        Threads to use, within an alignment and for aligning
                         pairs from a file at once [default: half the CPUs]

    --minscore <score>   Minimum required score
                         [default: match * MAX(0.2 * length, 2)]
    --maxhits <hits>     Maximum number of results per alignment
                         [default: no limit]

    --context <n>        Print <n> bases of context
    --printseq           Print sequences before local alignments
    --printmatrices      Print dynamic programming matrices
    --printfasta         Print fasta header lines
    --pretty             Print with a descriptor line
    --colour             Print with colour

  Experimental Options:
    --nogapsin1          No gaps allowed within the first sequence
    --nogapsin2          No gaps allowed within the second sequence
    --nogaps             No gaps allowed in either sequence
    --nogapsinlongest    No gaps allowed within the longer sequence
    --nomismatches       No mismatches allowed

 DETAILS:
  * For help choosing scoring, see the README file.
  * Gap (of length N) penalty is: (open+N*extend)
  * To do alignment without affine gap penalty, set '--gapopen 0'.
  * Scoring files should be matrices, with entries separated by a single
    character or whitespace. See files in the 'scoring' directory for examples.

  https://github.com/DanielTillett/nucleics-align  (compiled: <date>)
```

</details>

<details>
<summary><code>lcs</code></summary>

```text
./bin/lcs [options] <sequence>
  Print substrings in decreasing order of length
```

</details>
