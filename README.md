# ratePlacer

ratePlacer estimates the age of ancient DNA samples by maximum likelihood, using a
molecular clock on a reference phylogeny. It was built for ancient environmental DNA
(aeDNA): it combines many short, damaged reads from one or more taxa into a single sample
age, without building a consensus sequence or assuming a single source population.

Each read is placed on a branch of a fixed, ultrametric reference tree. The sample age is
the age at which the placed reads are jointly most likely. ratePlacer reports the
maximum-likelihood age and a 95% confidence interval. Per-base nucleotide likelihoods let
you account for sequencing error and post-mortem damage.

The repository contains two programs:

| Program      | What it does |
|--------------|--------------|
| `gfl`        | "Get fractional likelihoods". Precomputes, once per reference tree, the fractional likelihoods that ratePlacer needs. |
| `ratePlacer` | Places reads on the reference tree and estimates the sample age. |

> **Paper:** Lemmon-Kishi M, Pipes L, De Sanctis B, Nielsen R. *Molecular Clock Dating of
> Ancient Environmental DNA Reveals Damage Beyond Deamination.* bioRxiv (2026).
> https://doi.org/10.64898/2026.07.03.735781

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Workflow](#workflow)
- [Interpreting results](#interpreting-results)
- [Documentation](#documentation)
- [Repository layout](#repository-layout)
- [Citation](#citation)
- [License](#license)

## Installation

Requirements: a C compiler (gcc or clang), `make`, and optionally OpenMP for
multithreading. `gzip` must be on your `PATH` if you pass gzipped inputs. Building a
reference database with the provided R script also needs R with `phangorn` and
Bioconductor's `Biostrings` (see [Workflow](#1-build-the-reference-database-once-per-reference-set)).
[`environment.yml`](environment.yml) installs these and the other pipeline dependencies
with conda.

```bash
git clone https://github.com/lemmonquiche/ratePlacer.git
cd ratePlacer
make
```

This builds `bin/gfl` and `bin/ratePlacer`. To put them on your `PATH`:

```bash
make install                    # installs to /usr/local/bin
make install PREFIX=$HOME/.local  # or anywhere else
```

**Linux:** gcc supports OpenMP out of the box.

**macOS:** Apple clang needs Homebrew's OpenMP runtime for multithreading. Without it,
`make` still builds a single-threaded ratePlacer and prints a note.

```bash
brew install libomp
make clean && make
```

Other build options:

| Command | Effect |
|---------|--------|
| `make OPENMP=0` | Build without OpenMP (single-threaded). |
| `make NATIVE=0` | Leave out `-march=native`, e.g. for binaries that run on a cluster with different CPUs. |
| `make CC=gcc-14` | Use a different compiler. |
| `make debug` | Build `bin/*_debug` with AddressSanitizer and UBSan. |
| `make help` | List all targets and options. |

## Quick start

The `example/` directory holds a small simulated dataset: 12 reference taxa, and 600
reads from three ancient lineages that share a true age of 0.01 substitutions/site. To
run it:

```bash
make example
```

This runs `gfl` and then `ratePlacer` (see [`example/run_example.sh`](example/run_example.sh))
and prints the estimate, which should be close to 0.01:

```
reads_used=3
estimated_age=0.0077983155228159
CI_lower=0.0031770157624604
CI_upper=0.0124196152831714
```

`reads_used=3` is expected. With `-r 1`, ratePlacer merges the reads on each branch into
one sequence, and these reads fall on three branches.

## Workflow

```
 reference sequences ──► ultrametric tree + GTR+Γ fit ──► <tree>_reference.txt
                                                           <tree>_parameter.txt
                                                                   │
                                                                   ▼
                                                     gfl ──► <tree>_likelihood.txt
                                                                   │
 sample reads ──► align to reference ──► assignments (.assign)     │
                  + assign to node        aligned reads (.align) ──┼──► ratePlacer ──► age + 95% CI
                                          damage/error profile ────┘
```

### 1. Build the reference database (once per reference set)

Each reference tree is numbered 0, 1, …, N−1 and needs two files in one directory:

- `<tree>_reference.txt`: the reference alignment, followed by the tree in Newick
  format. The tree must be rooted, binary and **ultrametric**, with branch lengths in
  substitutions per site.
- `<tree>_parameter.txt`: GTR+Γ parameters (gamma shape, base frequencies, exchangeabilities).

[`scripts/make_reference.R`](scripts/make_reference.R)
creates both files from a reference alignment and a maximum-likelihood tree. In the paper,
the tree came from RAxML under GTRGAMMA.

```bash
Rscript scripts/make_reference.R ALIGNMENT.fasta TREE.newick TREE_NUM OUT_DIR [phy]
```

| Argument | Description |
|----------|-------------|
| `ALIGNMENT.fasta` | Reference multiple sequence alignment. A sequence's name is its FASTA header up to the first space, and it must match a tip label in the tree. |
| `TREE.newick` | Tree for the reference sequences. It can be unrooted; the script roots it at the midpoint. |
| `TREE_NUM` | The number of this tree (0, 1, …), used as the prefix of the output files. |
| `OUT_DIR` | Directory to write to. It is created if it doesn't exist. Use the same directory for all trees; it becomes ratePlacer's `-p` directory. |
| `phy` | Optional. Also writes the alignment in PHYLIP format. |

The script does the following:

1. Roots the tree at its midpoint and drops alignment sequences that are not in the tree.
   It stops with an error if a tip has no sequence or if two sequences share a name.
2. Makes the tree ultrametric with phangorn (`minEdge(tree, tau = 1e-5, enforce_ultrametric = TRUE)`).
3. Fits a GTR+Γ model with 4 rate categories (`pml`, then `optim.pml`). This keeps the
   tree topology and its ultrametric shape fixed, and re-estimates the branch lengths,
   rates, gamma shape and base frequencies.
4. Renumbers the tips 1 … n in the order they appear in the Newick string, and writes the
   files below.

| Output | Content |
|--------|---------|
| `<TREE_NUM>_reference.txt` | Input for gfl and ratePlacer: the alignment plus the ultrametric tree with numbered tips. |
| `<TREE_NUM>_parameter.txt` | Input for gfl and ratePlacer: the GTR+Γ parameters. |
| `<TREE_NUM>_pared.fa` | The alignment in tip-number order: the *k*-th sequence is leaf *k*. Use it to map leaf numbers back to sequence names, and as the reference alignment when you align reads. |
| `<TREE_NUM>_ultrametric.tree` | The fitted ultrametric tree with the original tip names. |
| `<TREE_NUM>_pared.phy` | Only with `phy`: the alignment in PHYLIP format. |

Requirements and notes:

- **Requirements:** R with `phangorn` and Bioconductor's `Biostrings`. `seqmagick` is
  also needed for the `phy` option.
- **Several trees:** run the script once per tree, with the same `OUT_DIR` and
  `TREE_NUM` = 0, 1, …. Then run `gfl` on that directory.
- **Gamma categories:** the script fits 4 gamma rate categories (`k = 4`). gfl and
  ratePlacer require exactly 4.

See [docs/file_formats.md](docs/file_formats.md) for the exact formats if you build these
files another way.

### 2. Precompute fractional likelihoods with `gfl`

```bash
bin/gfl reference_dir likelihood_dir NUM_TREES
```

`gfl` writes `<tree>_likelihood.txt` for each tree. It only needs to be rerun if the
reference tree or parameters change.

### 3. Prepare the sample inputs

ratePlacer takes three files describing the sample's reads:

- **Assignment file** (`-a`): the tree and node that each read is assigned to. Nodes can
  come from [tronko](https://github.com/lpipes/tronko)-assign, or you can supply them
  directly and let ratePlacer search for the best branch (`-k`).
- **Alignment file** (`-q`): each read aligned to its tree's reference alignment
  coordinates.
- **Error profile** (`-e`): per-base nucleotide log-likelihoods that encode sequencing
  error and DNA damage. The file can be empty, in which case ratePlacer assumes a flat 1%
  error rate.

The scripts in `scripts/` generate these files from mapped reads, including the damage
estimation. [docs/pipeline.md](docs/pipeline.md) walks through the full pipeline used in
the paper, from remapped reads to age estimates.

### 4. Estimate the sample age

```bash
bin/ratePlacer -a sample.assign -q sample.align -e sample.error \
    -l likelihood_dir -p reference_dir -n NUM_TREES \
    -m 3 -A 0 -c 0 -r 1 -k 1 -t 8 \
    -o results/sample
```

This writes `results/sample.summary` (the settings, `estimated_age`, `CI_lower` and
`CI_upper`) and `results/sample.log`. Run `bin/ratePlacer -h` for all options, or see
[docs/usage.md](docs/usage.md).

## Interpreting results

- **Units.** Ages are in expected substitutions per site, the same units as the reference
  tree's branch lengths. To convert to years, scale by a calibration. For example, if the
  ultrametric tree has root depth *D* (substitutions/site) and the root divergence is
  *T* years, then age in years = `estimated_age × T / D`.
- **What the confidence interval covers.** The 95% interval comes from the observed Fisher
  information. It conditions on the fixed reference tree, divergence calibration and
  substitution model, so it is narrower than a Bayesian HPD that integrates over these.
- **An age near zero (~3×10⁻⁸)** is the smallest value ratePlacer can return. In our
  simulations, this happened when damage was not corrected: the damage-induced mismatches
  lengthen the ancient branches. Check the damage correction before trusting such a
  result.
- **Damage correction changes the age.** Under-correcting biases ages younger, and
  over-correcting biases them older. We recommend running ratePlacer with both the high-α
  and low-α damage profiles and reporting the range.
- **Data requirements.** As a rough guide, older samples can be dated from about 1,000
  reads and younger samples need about 10,000. Damage correction raises these numbers.

ratePlacer assumes a strict molecular clock, an ultrametric reference tree with fixed
branch lengths, a single reference topology, and that all reads in a sample share one age.
See the Discussion of the paper for these and other limitations.

## Documentation

- [docs/pipeline.md](docs/pipeline.md): the paper's end-to-end analysis pipeline, from
  remapped reads through damage estimation to ratePlacer. The conda environment for its
  scripts is in [`environment.yml`](environment.yml).
- [ancientDNASymmetry](https://github.com/lemmonquiche/ancientDNASymmetry): the
  substitution symmetry and observed-to-expected analyses from the paper (Tables 1–3,
  Figs. 2–3).
- [docs/usage.md](docs/usage.md): command-line reference for `gfl` and `ratePlacer`,
  including analysis modes, read placement, merging, threading and troubleshooting.
- [docs/file_formats.md](docs/file_formats.md): input and output file formats, and how
  tree nodes are numbered.

## Repository layout

```
src/             C source for ratePlacer (RatePlacer.c, opt.h, minfunc.c) and gfl (get_frac_like.c)
example/         Simulated example dataset, its generator and run script
docs/            Usage, file format and pipeline documentation
scripts/         Pipeline scripts: reference building, read processing, damage estimation
                 and ratePlacer input generation (see docs/pipeline.md)
simulation/      Simulation pipeline from the paper
environment.yml  Conda environment for the pipeline scripts
Makefile
```

## Citation

If you use ratePlacer, please cite:

```bibtex
@article{lemmonkishi2026ratePlacer,
  title   = {Molecular Clock Dating of Ancient Environmental DNA Reveals Damage Beyond Deamination},
  author  = {Lemmon-Kishi, Maya and Pipes, Lenore and De Sanctis, Bianca and Nielsen, Rasmus},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.07.03.735781}
}
```

If you use tronko for read assignment, please also cite Pipes & Nielsen (2024).

## License

ratePlacer is released under the GNU General Public License v3.0; see [LICENSE](LICENSE).

It includes code by others:

- `src/minfunc.c`: Brent's minimization routines (Brent 1973), translated to C by
  David L. Swofford.
- The gamma-distribution and eigen-decomposition routines in `src/RatePlacer.c` and
  `src/get_frac_like.c`, written by Ziheng Yang.

## Contact

Questions and bug reports: please open a GitHub issue, or email Maya Lemmon-Kishi
(maya_lemmon-kishi@berkeley.edu).
