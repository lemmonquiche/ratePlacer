# File formats

All files are whitespace-delimited plain text. In the examples, `→` marks a tab, but any
whitespace works. The files in [`example/`](../example) are complete, working examples
of every input.

**Reference inputs** (one set per tree, in the `-p` directory):

- [`<tree>_reference.txt`](#tree_referencetxt): reference alignment and tree
- [`<tree>_parameter.txt`](#tree_parametertxt): GTR+Γ parameters
- [`<tree>_likelihood.txt`](#tree_likelihoodtxt-gfl-output): written by `gfl`

**Sample inputs:**

- [Assignment file (`-a`)](#assignment-file--a)
- [Alignment file (`-q`)](#alignment-file--q)
- [Error profile (`-e`)](#error-profile--e)

**Outputs:**

- [Outputs](#outputs)

**Node numbering:**

- [Node numbering](#node-numbering)

Trees are numbered `0` … `N−1`. The number is the filename prefix, e.g. `0_reference.txt`
or `1_parameter.txt`.

---

Both reference files are written by
[`scripts/make_reference.R`](../scripts/make_reference.R)
(see the [README](../README.md#1-build-the-reference-database-once-per-reference-set)).
The formats below matter if you build the files yourself.

## `<tree>_reference.txt`

The reference alignment, followed by the tree.

```
12 5000
CTAATCGATACTCGCAAATAAACTACACGGCCAAAATGCT...
TTAATCGATACTCGCAAATAAACTACACGGCCGAAATGCT...
...                                          (12 sequence lines in total)
(((3:0.0001778704,7:0.0001778704):0.0031334660,12:0.0033113364):0.0609250274,(...));
```

| Line(s) | Content |
|---------|---------|
| 1 | Number of sequences *n* (at least 3) and alignment length *L*. |
| 2 … *n*+1 | One aligned sequence per line, each *L* characters long, without names. Allowed characters: `A C G T` (any case). `N`, `-` and `~` mean missing data. |
| last | The tree in Newick format, on a single line, ending with `;` and a newline. |

Requirements for the tree:

- Leaves are labelled with integers `1` … *n*. Leaf *k* is the *k*-th sequence line.
- The tree is **rooted**, **strictly binary** and **ultrametric** (gfl allows a tolerance of
  about 1e-5).
- Every branch has a length, in expected substitutions per site.
- There are no internal node labels or support values. For example, `)0.95:0.01` is
  not allowed.
- The file ends with a newline after the tree. ratePlacer locates the tree by reading the
  last line of the file.

## `<tree>_parameter.txt`

```
0.5
0.3 0.2 0.2 0.3
1.0 4.0 1.0 1.0 4.0 1.0
```

| Line | Content |
|------|---------|
| 1 | Gamma shape parameter α for rate variation across sites. 4 discrete categories are used. |
| 2 | Base frequencies π<sub>A</sub> π<sub>C</sub> π<sub>G</sub> π<sub>T</sub>. They must sum to 1 within 3e-7. |
| 3 | GTR exchangeabilities in the order **AC AG AT CG CT GT**. |

The values are read as a sequence of 11 numbers, so line breaks don't matter. The rate matrix is normalized
to one expected substitution per unit time, so only the relative values of the
exchangeabilities matter. With phangorn, these are `fit$shape`, `fit$bf` and `fit$Q` of a
`pml` fit, which uses the same AC, AG, AT, CG, CT, GT order.

## `<tree>_likelihood.txt` (gfl output)

`gfl` writes this file and ratePlacer reads it. You do not need to edit it. It contains:

1. A header: `n L 4 P`, where *P* is the number of unique site patterns.
2. The relative rates of the 4 gamma categories.
3. 2n − 1 lines of `node: age branch_length`, sorted by the age of the branch's top.
4. The tree height.
5. The pattern index of each alignment column (*L* values).
6. For each gamma category (`C1` … `C4`) and each site pattern (`S1:` …), the
   log fractional likelihoods for every node.

The file must match the `<tree>_reference.txt` it was built from. Rerun `gfl` if the
reference changes.

---

## Node numbering

The assignment file refers to tree nodes by number. Each node *v* stands for the branch
**above** it, from *v* to its parent. A tree with *n* leaves has 2n − 1 nodes, numbered
from 1 as follows:

- **Leaves:** leaf label *k* → node *k* (1 … *n*).
- **Internal nodes:** take the Newick string and number its commas 1 … n−1 from left to
  right. Every internal node has exactly one comma that separates its two children. If
  that comma is the *c*-th comma, the node number is *n* + *c* (*n*+1 … 2n−1).

For example, with *n* = 4:

```
((1:0.1,2:0.1):0.2,(3:0.15,4:0.15):0.15);
   ,1          ,2      ,3
```

| Node | Branch above |
|------|--------------|
| 1, 2, 3, 4 | the terminal branches of leaves 1–4 |
| 5 (comma 1) | the branch above the clade (1,2) |
| 6 (comma 2) | the root. It has no branch above it. |
| 7 (comma 3) | the branch above the clade (3,4) |

Reads assigned to the root must be moved onto a branch with `-k 1` or `-k 2` (see
[usage.md](usage.md#read-placement--k)).

Two other places use **0-based** node numbers, i.e. one less than above: the
`--toss-branches` file, and the node column of the mode 7 output.

---

## Assignment file (`-a`)

```
600
0→14
0→16
0→16
...
```

| Line(s) | Content |
|---------|---------|
| 1 | Number of reads *R*. |
| 2 … *R*+1 | One line per read: the tree number (0-based) and the node number (1-based, see [Node numbering](#node-numbering)). |

The reads must be in the same order as in the alignment file and the error profile.

## Alignment file (`-q`)

```
600
62→3948→TACTTAGGTGAACGTACTCATAGTGCGTATTTAATTCAACAGACTCTTGTCCTGGGCGGGGG
64→2670→GGCGTGCAGTTCATTACATCCAATGCGCGCAACCAGTAAGAGAGTTGGGGAATCAACATGGTAT
...
```

| Line(s) | Content |
|---------|---------|
| 1 | Number of reads *R*. It must match the assignment file. |
| 2 … *R*+1 | One line per read: the length ℓ, the start position, and the sequence. |

- **Start** is the 0-based column of the reference alignment (of the read's assigned
  tree) where the read begins.
- **Sequence** is the read as aligned to the reference alignment. It has exactly ℓ
  characters, one per reference column from `start` to `start + ℓ − 1`. Deletions relative
  to the reference are written as `-`. Insertions relative to the reference cannot be
  represented, so remove them. Allowed characters are `A C G T`, plus `N`, `-` and `~` for
  missing data.

## Error profile (`-e`)

The error profile gives, for each read base, the probability of the observed base given
each possible true base. This is where sequencing error and post-mortem damage enter the
model: the per-site probability *p<sub>s</sub>(a)* in Eq. 1 of the paper.

```
0→0→-4.605170→-4.605170→-4.605170→-0.030459
0→1→-0.030459→-4.605170→-4.605170→-4.605170
...
```

Each line is one read position:

| Column | Content |
|--------|---------|
| 1 | Read index (0-based, in the order of the alignment file). |
| 2 | Position within the read's aligned sequence (0-based, counting `-` characters). |
| 3–6 | Natural-log likelihoods log P(observed base \| true base = A), (… = C), (… = G), (… = T). |

- **Positions that are not listed** get a flat 1% error: log(0.99) if the true base equals
  the observed base, and log(0.01/3) otherwise. An **empty file** therefore applies a 1%
  error to every base.
- Entries that point at a gap or missing base are ignored, with a warning.
- Lines can be in any order and must be shorter than 256 characters.

To build these profiles from base quality scores and an estimated damage matrix, see
Section 4.3 of the paper and the scripts in `scripts/`.

---

## Outputs

| File | Written | Content |
|------|---------|---------|
| `PREFIX.summary` | always | `key=value` lines: the settings used, then the results. The results depend on the mode; see [usage.md](usage.md#analysis-modes). |
| `PREFIX.log` | always | A one-line summary of the estimate. With `-A 1`, also an estimate for each tree. With `-D`, also the placement diagnostics. |
| `PREFIX.reads` | `-R` | The merged sequences, one header line plus the sequence per branch. |
| `PREFIX.tossed` | `-T` | The discarded reads in FASTA format (`>read_i tree=… node=… start=…`). |

Example `PREFIX.summary` from mode 3:

```
assign_file=example/sample.assign
...
num_threads=1
toss_root_reads=0
debug_placements=0
input_reads=3
reads_used=3
estimated_age=0.0077983155228159
likelihood=19278.1937129003454174
CI_lower=0.0031770157624604
CI_upper=0.0124196152831714
```

| Key | Meaning |
|-----|---------|
| `input_reads` | Number of sequences after placement, merging and tossing. |
| `reads_used` | Number of sequences left after mode 3 drops incompatible reads. |
| `estimated_age` | Maximum-likelihood sample age, in substitutions per site. |
| `likelihood` | **Negative** log-likelihood at the estimate. |
| `CI_lower`, `CI_upper` | 95% confidence interval from the observed Fisher information. |
