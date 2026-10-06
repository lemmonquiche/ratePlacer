# Usage

- [gfl](#gfl)
- [ratePlacer](#rateplacer)
  - [Options](#options)
  - [Analysis modes](#analysis-modes)
  - [Read placement (`-k`)](#read-placement--k)
  - [Merging reads per branch (`-r`, `-C`)](#merging-reads-per-branch--r--c)
  - [Discarding reads near the root (`-T`, `-b`)](#discarding-reads-near-the-root--t--b)
  - [Multiple reference trees](#multiple-reference-trees)
  - [Threads and performance](#threads-and-performance)
- [Settings used in the paper](#settings-used-in-the-paper)
- [Troubleshooting](#troubleshooting)

File formats are described in [file_formats.md](file_formats.md).

---

## gfl

```
gfl INPUT_DIR OUTPUT_DIR NUM_TREES
```

| Argument     | Description |
|--------------|-------------|
| `INPUT_DIR`  | Directory containing `<tree>_reference.txt` and `<tree>_parameter.txt` for trees `0` … `NUM_TREES-1`. |
| `OUTPUT_DIR` | Directory to write `<tree>_likelihood.txt` to. It is created if it does not exist. |
| `NUM_TREES`  | Number of reference trees. |

For each tree, `gfl` does the following:

1. Reads the alignment, tree and GTR+Γ parameters.
2. Discretizes the gamma distribution into 4 rate categories.
3. Collapses identical alignment columns (site patterns).
4. Computes the fractional likelihoods for every branch. A postorder pass gives the
   likelihood of the subtree below each branch, and a preorder pass gives the likelihood of
   the rest of the tree.

`gfl` exits with `Tree not ultrametric` if the root-to-tip distances differ by more than
about 1e-5.

You only need to run `gfl` once per reference set. Its output is reused for every sample.
The output is large: it holds (2n − 1) nodes × site patterns × 4 categories × up to 8
values per tree.

---

## ratePlacer

```
ratePlacer -a ASSIGN_FILE -q QUERY_FILE -e ERROR_FILE -l LIK_DIR -p PARAM_DIR \
           -n NUM_TREES -m MODE -A ALL_TREES -c COMPARE_AGE -r MERGE_READS \
           -k REASSIGN -o OUT_PREFIX [options]
```

### Options

All options in the first two tables are required. Run `ratePlacer -h` to print this list.

**Inputs and outputs**

| Option | Description |
|--------|-------------|
| `-a`, `--assign-file FILE` | Read assignments: the tree and node for each read. |
| `-q`, `--query-file FILE`  | Reads aligned to the reference coordinates. |
| `-e`, `--error-file FILE`  | Per-base nucleotide log-likelihoods (error/damage profile). Can be empty. |
| `-l`, `--lik-dir DIR`      | `gfl` output directory (`<tree>_likelihood.txt`). |
| `-p`, `--param-dir DIR`    | Reference directory (`<tree>_reference.txt`, `<tree>_parameter.txt`). |
| `-n`, `--num-trees N`      | Number of reference trees. |
| `-o`, `--out-prefix PREFIX`| Output prefix. ratePlacer writes `PREFIX.summary` and `PREFIX.log`. The directory must already exist. |

`-a`, `-q` and `-e` can be gzipped (`.gz`). The files in `-l` and `-p` must be uncompressed.

**Analysis settings**

| Option | Values | Description |
|--------|--------|-------------|
| `-m`, `--mode`        | 0–9 | Analysis to run. See [Analysis modes](#analysis-modes). Use `3` to estimate a sample age. |
| `-A`, `--all-trees`   | 0/1 | In modes 3, 6 and 9, `1` also estimates a separate age for each tree and writes it to the log. The joint estimate over all trees is always computed. |
| `-c`, `--compare-age` | number | Age tested by the likelihood ratio test in mode 2. Ignored in other modes, so use `0`. |
| `-r`, `--merge-reads` | 0/1 | `1` merges reads assigned to the same branch into one sequence. See [Merging](#merging-reads-per-branch--r--c). |
| `-k`, `--reassign`    | 0/1/2 | How to place reads. See [Read placement](#read-placement--k). |

**Optional**

| Option | Default | Description |
|--------|---------|-------------|
| `-C`, `--coverage T` | `0.05f` | With `-r 1`, drops merged sequences that cover fewer reference sites than `T`. Write `T` as a fraction of the reference length (`0.05f`) or as a number of bases (`500b`). |
| `-t`, `--threads N`  | 1 | Number of OpenMP threads. |
| `-R`, `--write-reads`| off | With `-r 1`, writes the merged sequences to `PREFIX.reads`. |
| `-T`, `--toss-root-reads` | off | Discards reads placed on the two branches directly below the root and writes them to `PREFIX.tossed`. |
| `-b`, `--toss-branches FILE` | — | With `-T`, discards reads on the branches listed in `FILE` instead of the root branches. |
| `-D`, `--debug-placements` | off | After merging, checks each merged sequence against neighbouring branches and writes the results to the log. |
| `-h`, `--help` | | Prints usage. |

### Analysis modes

Mode 3 is the main analysis. It estimates the age of the whole sample. The other modes
are diagnostics, or analyses of single reads.

| Mode | Analysis | Output in `PREFIX.summary` |
|------|----------|----------------------------|
| **3** | **Sample age.** Finds the maximum-likelihood age shared by all reads, with a 95% confidence interval from the observed Fisher information. | `reads_used`, `estimated_age`, `likelihood`, `CI_lower`, `CI_upper` |
| 6 | Like mode 3, but keeps every read. The age is therefore bounded by the youngest branch top among the reads. | as mode 3 |
| 9 | Like mode 3, but moves incompatible reads to their parent branch instead of dropping them. | as mode 3 |
| 1 | Age and placement of each read, estimated separately. | One block per read |
| 0 | Like mode 1, plus the per-site likelihood scores for each read. | One block per read |
| 2 | Per-read likelihood ratio test of the read's best age against `-c COMPARE_AGE`. | One block per read |
| 4 | Joint likelihood of the sample evaluated on a grid of ages, as a likelihood profile. | One block per age |
| 5 | Like mode 4, for ages up to 0.005 only. | One block per age |
| 7 | Per-read likelihood at points along the assigned branch. | CSV: `read,tree,node,position,likelihood` (node is 0-based) |
| 8 | Per-read likelihood surface over age × placement. | CSV: `read,alpha,placement,likelihood,age` |

**How mode 3 handles reads.** A read cannot be older than the top (parent node) of the
branch it is placed on. Before optimizing, mode 3 does a coarse search over ages. It drops
reads on branches younger than the candidate age when doing so improves the fit. The
number of reads that remain is reported as `reads_used`. ratePlacer then optimizes the age
with Brent's method. In each step, it optimizes each read's position along its branch
given the candidate age.

**Likelihood values.** The `likelihood` reported in the summary and log is the
**negative** log-likelihood, so smaller is better.

### Read placement (`-k`)

The assignment file gives each read a starting node. `-k` controls what ratePlacer does
with it:

| `-k` | Behaviour | Use when |
|------|-----------|----------|
| 0 | Uses the given node as is. | The branch is known, e.g. in simulations with true placements. |
| 1 | Greedy maximum-likelihood search over branches, starting at the given node and moving up or down the tree until the likelihood stops improving. | You supply a starting node directly, without tronko. |
| 2 | Compares the given branch with the branches immediately around it and keeps the best. | Starting nodes come from tronko-assign. |

A read is placed on the branch *above* its node, which runs from the node to its parent.
The root has no branch above it. Reads assigned to the root therefore need `-k 1` or
`-k 2`, which move them onto one of the root's child branches. With `-k 0`, a root
assignment produces a likelihood of 1e9.

### Merging reads per branch (`-r`, `-C`)

With `-r 1`, ratePlacer merges all reads placed on the same branch of the same tree into
a single sequence that spans the reference positions they cover. This happens after
placement (`-k`). Where reads overlap, ratePlacer combines their per-base nucleotide
likelihoods. Merged sequences that cover fewer sites than the `-C` threshold are dropped.
The paper used merging with `-C 500b` in the simulations and the default `-C` for Kap København.

Merging greatly reduces run time and improves accuracy. It avoids optimizing a separate
placement for each short read, which carries little information on its own. After
merging, `reads_used` counts merged sequences, not raw reads.

### Discarding reads near the root (`-T`, `-b`)

Reads near the root are often uninformative or misassigned. `-T` discards reads placed on
the two branches directly below the root of each tree. To discard reads from other
branches instead, list them in a file and pass it with `-b`. Each line holds a tree
number and a 0-based node number (one less than the node number used in the assignment
file):

```
0 14
0 17
```

### Multiple reference trees

A sample can be dated against several reference trees at once, e.g. one tree per taxon or
genomic region. Number the trees `0` … `N−1`, run `gfl` on all of them, and give each read
the tree it belongs to in the assignment file. Mode 3 estimates one age shared across all
trees. Add `-A 1` to also log an age estimated from each tree on its own. Trees with no
assigned reads are not loaded.

### Threads and performance

Use `-t N` to run on N threads. This requires an OpenMP build; see the README. The input
reading, read placement (`-k 1`/`-k 2`) and likelihood steps run in parallel.

Memory is dominated by the fractional likelihoods of the trees in use. Placement with
`-k 1` or `-k 2` is the slowest step when there are many unmerged reads.

---

## Settings used in the paper

| Setting | Value |
|---------|-------|
| Mode | `-m 3` (sample age estimation) |
| Merging | `-r 1`, with `-C 500b` in the simulations and the default `-C 0.05f` for Kap København |
| Placement | tronko-assign, then ML over the surrounding branches. Or supply a node directly and use ratePlacer's ML search. |
| Damage | Nucleotide likelihoods from the damage model (Section 4.3 of the paper). Run once with the high-α profile and once with the low-α profile, and report both ages. |

---

## Troubleshooting

| Message | Cause |
|---------|-------|
| `Tree not ultrametric` (gfl) | The root-to-tip distances differ. Make the tree ultrametric before running gfl. |
| `Nucleotide frequencies not properly scaled` | The base frequencies in `<tree>_parameter.txt` must sum to 1 (within 3e-7). |
| `Number of sequences and sequence lengths do no match in input files` | `<tree>_reference.txt` changed after `gfl` was run. Rerun `gfl`. |
| `Different number of query sequences found in query data file and in assignment file` | The first lines of `-a` and `-q` disagree, or the files have different numbers of reads. |
| `Error reading assignments for read …` | A node number is out of range for its tree. Nodes are 1 … 2n−1. See [file_formats.md](file_formats.md#node-numbering). |
| `BAD BASE (x)` | A sequence contains a character other than `ACGTN-~`. |
| `Warning, error profile includes positions not in query alignment` | An error profile entry points at a gap or missing base. The entry is ignored. |
| `Warning! Age estimate was near boundary` | The estimate is close to the oldest age allowed by the remaining reads. Check the read placements. |
| `estimated_age` ≈ 3e-8 | ratePlacer's lower bound. In our simulations, this indicated uncorrected damage. See the README. |
| `Cannot open summary file` | The directory in `-o` does not exist. |
