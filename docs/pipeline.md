# Analysis pipeline

This is the pipeline used to date the nuclear reads from Kap København in the paper
(Section 4.5). It starts from reads that have already been assigned to one taxon, and
ends with ratePlacer age estimates that bracket the effect of damage correction. Read
placement is done by ratePlacer's own maximum-likelihood search, without tronko.

- [Overview](#overview)
- [Requirements](#requirements)
- [Before you start](#before-you-start)
  - [Build the nuclear reference MSA](#build-the-nuclear-reference-msa)
- [1. Remap reads to the compartment database](#1-remap-reads-to-the-compartment-database)
- [2. Split reads by genomic compartment](#2-split-reads-by-genomic-compartment-filter_compartmentspy)
- [3. Trim the reference MSA and convert reads to MSA coordinates](#3-trim-the-reference-msa-and-convert-reads-to-msa-coordinates-filter_msa_by_sampy)
- [4. Build the reference tree](#4-build-the-reference-tree)
- [5. Make the ratePlacer reference files](#5-make-the-rateplacer-reference-files-make_referencer)
- [6. Precompute fractional likelihoods](#6-precompute-fractional-likelihoods-gfl)
- [7. Count mismatches for damage estimation](#7-count-mismatches-for-damage-estimation-generate_matricespy)
- [8. Estimate damage](#8-estimate-damage-run_analysispy)
- [9. Generate the ratePlacer inputs](#9-generate-the-rateplacer-inputs-generaterateplacerinputs_msareads_damagematrixpy)
- [10. Estimate the sample age](#10-estimate-the-sample-age-rateplacer)
- [Variations used in the paper](#variations-used-in-the-paper)
- [Scripts reference](#scripts-reference)

## Overview

```
taxon reads (bamdam extract)
   │
   ▼
1. bowtie2 remap ──► full_database.sam
   │
   ▼
2. filter_compartments.py ──► <compartment>.sam
   │                               │
   ▼                               ▼
3. filter_msa_by_sam.py       7. generate_matrices.py ──► A, C, D, ref_counts (.npy)
   │                     │         (uses 0_parameter.txt)            │
   ▼                     ▼                                           ▼
trimmed MSA         reads.msa.fq                   8. run_analysis.py ──► B_lo.npy, B_hi.npy
   │                     │                                           │
   ▼                     └───────────────────┐                       │
4. RAxML ──► tree                            ▼                       │
   │                     9. generateRatePlacerInputs_msaReads_damageMatrix.py
   ▼                        (once with B_lo, once with B_hi) ◄───────┘
5. make_reference.R ──► 0_reference.txt, 0_parameter.txt
   │                                         │
   ▼                                         ▼
6. gfl ──► 0_likelihood.txt ──────► 10. ratePlacer (lo and hi) ──► sample ages
```

The commands below use one taxon, the nuclear compartment and one reference tree
(tree `0`), with file names chosen for illustration.

## Requirements

Everything except Cactus, `gfl` and `ratePlacer` is in the conda environment in
[`environment.yml`](../environment.yml):

```bash
conda env create -f environment.yml
conda activate rateplacer
```

| Tool | Used in |
|------|---------|
| [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) (progressive Cactus and `cactus-hal2maf`; 9.1.2 in the paper) | reference MSA |
| [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/) | step 1 |
| Python ≥ 3.10 with `pysam`, `numpy`, `scipy`, `matplotlib` | reference MSA, steps 2, 3, 7, 8, 9 |
| [RAxML](https://github.com/stamatak/standard-RAxML) (8.2.12 in the paper) | step 4 |
| R with `phangorn` and `Biostrings` | step 5 |
| `gfl`, `ratePlacer` (this repository) | steps 6, 10 |

## Before you start

The pipeline expects the following, which are produced upstream:

- **Reads for one taxon.** These come from competitive mapping and taxonomic assignment.
  Section 4.5.1 of the paper gives the commands and parameters. In summary:
  - reads were quality-filtered and deduplicated with fastp and SGA;
  - they were mapped with `bowtie2 -k 100` to eukaryotic and microbial databases;
  - they were assigned taxonomically with ngsLCA (`-simscorelow 0.95 -simscorehigh 1.0`);
  - damage was summarized with `bamdam shrink` and `bamdam compute`.

  A taxon was classified as ancient if it had mean damage > 0.4, DUST < 7 and more than
  10,000 reads. The taxon's reads were then extracted with
  `bamdam extract --only-top-alignment`.

  > **To be added:** the script for classifying taxa as ancient or modern.

- **A reference MSA** for the compartment being dated, containing the focal genus and its
  outgroup. See [Build the nuclear reference MSA](#build-the-nuclear-reference-msa).

- **A compartment database.** This is the Bowtie 2 index that reads are remapped to. It
  contains the nuclear, mitochondrial and chloroplast references of the focal genus and
  its outgroup. Which genomes to include is up to you; the paper used assemblies from
  NCBI RefSeq and Phytozome (Supplementary S2.4). Two files describe it:

  `compartments.csv` maps every reference sequence name in the index to its compartment:

  ```
  reference,compartment
  Betula_pendula,nuclear
  Betula_pendula_chloroplast,chloroplast
  Betula_pendula_mito,mitochondria
  ...
  ```

  `outgroup.txt` lists the outgroup reference names, one per line.

  Steps 3 and 9 map reads onto the MSA by reference name. So for the compartment you
  analyse, the sequences in the index must be the MSA's sequences with the gaps removed,
  under the same names. For the nuclear compartment, these are the stitched genome rows
  from `maf_to_fasta.py`, not the original chromosomes:

  ```bash
  sed '/^>/!s/-//g' nuclear_msa.fasta > nuclear_ungapped.fasta
  ```

### Build the nuclear reference MSA

The nuclear genomes of the focal genus and its outgroup are aligned with progressive
Cactus. The alignment is exported as a MAF anchored on one genome from the focal genus,
and then stitched into a single FASTA alignment. The examples below use Betula with Alnus
as the outgroup.

**1. Align the genomes with Cactus.** The seq file lists one `name path` pair per genome.
A first line with a Newick guide tree is optional. Genome names must not contain `.`,
because step 3 takes the genome name to be everything before the first `.`.

```
Betula_pendula   /path/to/Betula_pendula.fa
Betula_nana      /path/to/Betula_nana.fa
Alnus_glutinosa  /path/to/Alnus_glutinosa.fa
```

```bash
cactus ./jobstore betula_alnus_seq_file.txt betula_alnus.hal \
    --maxCores 20 --maxMemory 500Gi --consMemory 500Gi \
    --batchSystem single_machine --workDir ./temp/
```

**2. Export a MAF anchored on a focal-genus genome.** `--refGenome` sets the backbone
genome. The paper used a genome from the focal genus, because mapping rates varied with
the choice of reference.

```bash
cactus-hal2maf ./jobstore_maf betula_alnus.hal betula_alnus.maf.gz \
    --refGenome Betula_pendula --noAncestors --outType single \
    --chunkSize 500000 --batchCores 40 --batchMemory 750Gi
```

**3. Stitch the MAF into one FASTA alignment** (`maf_to_fasta.py`).

```bash
python scripts/maf_to_fasta.py betula_alnus.maf.gz nuclear_msa.fasta
```

The script concatenates the alignment blocks in MAF order. It writes one row per genome,
named by genome. Gaps fill the blocks where a genome is absent. Then it removes the
columns that contain no A, C, G or T in any genome. If a genome has more than one
sequence in a block, it keeps the first and prints a warning. The input can be gzipped.

Populus and Salix share one Cactus alignment in the paper, but are dated in separate
ratePlacer runs, each with its own `--refGenome`.

## 1. Remap reads to the compartment database

Build a Bowtie 2 index of the compartment database. For the nuclear compartment, the
references are the aligned references (the gap-stripped MSA rows, see
[Before you start](#before-you-start)), not the original chromosomes. The other
compartments and the outgroup are included alongside them.

```bash
bowtie2-build compartment_refs.fasta full.bowtie
```

Then remap the taxon's reads:

```bash
bowtie2 -p 8 -k 100 -N 0 -x full.bowtie -U taxon_reads.fq --no-unal -S full_database.sam
```

- `-k 100` reports up to 100 alignments per read. Step 2 needs these to detect
  multimapping and hits to the outgroup.
- `-N 0` allows no mismatches in the seed (the Bowtie 2 default). Use 1 to allow for a seed mismatch.
- `--no-unal` leaves unaligned reads out of the SAM.

The reads are filtered in step 2.

## 2. Split reads by genomic compartment (`filter_compartments.py`)

```bash
python scripts/filter_compartments.py \
    --input full_database.sam --compartments compartments.csv --outgroup outgroup.txt \
    --outdir compartments/ --prefix sample_
```

The script groups the alignments by read and **discards** a read if any of the
following holds:

- none of its alignments are mapped;
- any alignment hits an outgroup reference;
- any alignment hits a reference that isn't in `compartments.csv`;
- it maps more than once to the same reference;
- its alignments span more than one compartment;
- it doesn't have exactly one primary alignment;
- the primary alignment's score is below `--min-score` (default `-24`).

For each kept read, the script writes the primary alignment to
`compartments/sample_<compartment>.sam`. It also writes the alignment score distribution
to `compartments/sample_<compartment>_AS_distribution.tsv`. The input does not need to be
sorted; it is name-sorted to a temporary file if needed.

## 3. Trim the reference MSA and convert reads to MSA coordinates (`filter_msa_by_sam.py`)

```bash
python scripts/filter_msa_by_sam.py \
    -s compartments/sample_nuclear.sam -m nuclear_msa.fa \
    -o nuclear_trimmed.fa --extract-reads sample_nuclear.msa.fq -v
```

The script writes two files:

- **`nuclear_trimmed.fa`**: the MSA reduced to the columns covered by at least one read,
  with all-gap columns removed. This is the reference alignment for all later steps.
- **`sample_nuclear.msa.fq`**: the reads converted to `nuclear_trimmed.fa` coordinates.
  Deletions relative to the reference become `-` (quality `!`) and insertions are
  removed. Each header records the read's position:
  `@read start=<column> len=<length> ref=<reference> rc=<0|1>`.

To date several samples against one reference, pass all their SAM files to `-s` and one
FASTQ per SAM to `--extract-reads`. The trimmed MSA then covers the union of their reads.

Other options:

- `--trim`, `--trim-5`, `--trim-3`: trim read ends before conversion.
- `--coverage-stats`, `--depth-file`: report read depth.
- `--validate-mismatches FILE`: check the converted reads against the MD tags.

## 4. Build the reference tree

Estimate a maximum-likelihood tree from the trimmed MSA. The paper used RAxML 8.2.12
under GTRGAMMA:

```bash
raxmlHPC-PTHREADS -s nuclear_trimmed.fa -n nuclear -m GTRGAMMA -p 12345 -T 10
```

The tree is `RAxML_bestTree.nuclear`.

## 5. Make the ratePlacer reference files (`make_reference.R`)

```bash
Rscript scripts/make_reference.R nuclear_trimmed.fa RAxML_bestTree.nuclear 0 reference/
```

The script roots the tree at its midpoint, makes it ultrametric, and fits GTR+Γ with 4
rate categories. It writes:

- `reference/0_reference.txt` and `reference/0_parameter.txt`, which ratePlacer reads;
- `reference/0_pared.fa`, which maps leaf numbers back to sequence names.

See the [README](../README.md#1-build-the-reference-database-once-per-reference-set) for
details. `0_parameter.txt` is also the GTR input for damage estimation in step 7.

## 6. Precompute fractional likelihoods (`gfl`)

```bash
bin/gfl reference/ likelihoods/ 1
```

## 7. Count mismatches for damage estimation (`generate_matrices.py`)

```bash
python scripts/generate_matrices.py compartments/sample_nuclear.sam \
    --params reference/0_parameter.txt --n-positions 10 \
    --output-dir damage/ --prefix sample_
```

The script uses each read's primary alignment against the compartment references. The
references come from the MD tags, so the SAM must have MD tags; Bowtie 2 adds them. The
counts are taken in the orientation of the original molecule. The script writes the
following to `damage/`:

| File | Content |
|------|---------|
| `sample_A.npy` | The GTR rate matrix Q from `0_parameter.txt`. A(α) = exp(αQ). |
| `sample_C.npy` | The sequencing error matrix, from the mean base qualities. |
| `sample_D.npy` | The observed mismatch matrices, shape (2N+1, 4, 4). They cover 5′ positions 1…N, one middle bin, and 3′ positions N…1. |
| `sample_ref_counts.npy` | Counts of each reference base in each position bin, shape (2N+1, 4). |
| `sample_ref_freq.npy` | The same as frequencies. |

`--n-positions 10` gives the paper's 21 bins.

**Note:** `--prefix` must end in `_` (`sample_`), because step 8 looks for
`<prefix>_A.npy` and is given the prefix without the underscore (`sample`).

## 8. Estimate damage (`run_analysis.py`)

```bash
python scripts/run_analysis.py --prefixes sample --data-dir damage/ \
    --out-dir damage/results/ --n-ref sample damage/sample_ref_counts.npy
```

This fits the damage model of Section 4.3 of the paper. It first finds the feasible
interval of α, then profiles the likelihood over α to get a 95% confidence interval.
Finally it evaluates the damage matrices B at the two ends of that interval. It writes
the following to `damage/results/`:

| File | Content |
|------|---------|
| `sample_B_lo.npy`, `sample_B_hi.npy` | Damage matrices at the lower and upper CI bounds of α, shape (2N+1, 4, 4). These are the inputs to step 9. |
| `sample_alpha.npz` | The feasible interval, the MLE and CI of α, the likelihood profile, and diagnostics: `alpha_max`, `ci_hi_censored`, and `B_lo_method`. |
| `sample_profile_likelihood.png`, `sample_damage_profiles.png`, `sample_feasibility_scan.png` | Diagnostic plots. |

`--n-ref` uses the reference base counts from step 7 as the multinomial sample sizes. The
alternative, `--n-reads N --read-length L`, derives the counts from the number of reads
and the read length.

α is searched up to α<sub>max</sub> = 1/max<sub>i</sub>(−Q<sub>ii</sub>), the largest valid
value (Section 4.3.1). `--alpha-upper` can lower this limit but not raise it.

Two diagnostics are worth checking:

- **Censored upper bound.** If the likelihood stays within 1.92 of its maximum all the way
  to the upper limit, the script prints a `WARNING` and sets `ci_hi_censored` in the npz.
  `B_hi` is then evaluated at the limit, not at a likelihood-based bound.
- **How `B_lo` was solved.** `B_lo_method` records whether `B_lo` was computed directly
  as A⁻¹DC⁻¹ (`analytic`), which is exact when the lower bound lies in the feasible
  interval. Otherwise it is fitted with the constrained solver (`constrained`), the same
  way as `B_hi`.

`damage_estimator.py` and `profile_likelihood.py` are modules used by this script and are
not run directly.

## 9. Generate the ratePlacer inputs (`generateRatePlacerInputs_msaReads_damageMatrix.py`)

Run this once for each damage estimate:

```bash
for bound in lo hi; do
  python scripts/generateRatePlacerInputs_msaReads_damageMatrix.py \
      -i sample_nuclear.msa.fq -m nuclear_trimmed.fa -n NODE --doublestrand \
      --damage-matrix damage/results/sample_B_${bound}.npy \
      -b sample.align -c sample.assign -e sample.${bound}.error -o sample.map.txt
done
```

| Option | Description |
|--------|-------------|
| `-i` | Reads in MSA coordinates, from step 3. |
| `-m` | The **trimmed** MSA from step 3. The read coordinates refer to it. |
| `-n` | The node every read is assigned to on tree 0. ratePlacer's `-k 1` search then places each read starting from this node. The node number depends on the taxon's tree; see [node numbering](file_formats.md#node-numbering). |
| `--doublestrand` / `--singlestrand` | Library type (required). |
| `--damage-matrix` | Damage matrices from step 8. |
| `-b`, `-c`, `-e` | Output alignment (`-q`), assignment (`-a`) and error profile (`-e`) files for ratePlacer. |
| `-o` | Output table mapping read names to ratePlacer read indices. Give it a `.txt` name, because the removed-reads report is named after it. |

For each base, the nucleotide likelihoods combine the read's base quality with the damage
matrix for the base's position in the molecule (Section 4.3.3 of the paper).

Reads with more than `--max-mismatches` mismatches to their reference (default 5) are
removed. So are reads that lack information, e.g. a quality string. The paper used the
default, and in the Kap København data no reads were removed by the mismatch filter. The alignment and
assignment files are the same in both runs; only the error profile differs.

Without `--damage-matrix`, the error profiles come from base quality alone. Combined with
`--damage-masking`, `--mask-transitions` or `--trim-ends`, this reproduces the paper's
"age estimates without the damage model".

## 10. Estimate the sample age (`ratePlacer`)

```bash
for bound in lo hi; do
  bin/ratePlacer -a sample.assign -q sample.align -e sample.${bound}.error \
      -l likelihoods/ -p reference/ -n 1 \
      -m 3 -A 0 -c 0 -k 1 -r 1 -t 16 \
      -o results/sample.${bound}
done
```

Merged sequences covering less than 5% of the reference are dropped (the `-C` default,
which the paper used for Kap København; the simulations used `-C 500b`). Report the age
range spanned by the `lo` and `hi` estimates. To convert to years, use a
divergence time for the reference tree (see the
[README](../README.md#interpreting-results)).

---

## Variations used in the paper

### Dating a consensus sequence (Betula chloroplast)

To compare with the BEAST analysis of Kjær et al. (2022), we took the reference
sequences and the Kap København consensus from their BEAUti XML. The references went
through steps 4–6. The consensus was treated as a single read spanning the alignment.
For that kind of input:

- the assignment and alignment files hold one entry, the consensus, with start `0`;
- the error file is empty, because a consensus from sites with ≥ 20× depth has no damage
  profile;
- run with `-r 0`, as there is nothing to merge.

The paper's run used an earlier ratePlacer without an error model. The current version
applies a flat 1% error when the error file is empty, so results can differ slightly.

### Age estimates without the damage model

Before the damage model, error profiles came from base quality alone, and damage was
handled by masking or trimming the read ends. With the current step 9 script, leave out
`--damage-matrix` and add one of these option sets:

| Treatment | Options |
|-----------|---------|
| Mask possible deamination in the first and last 7 bases | `--damage-masking ends --damage-masking-bases 7 --damage-types deamination` |
| Mask all transitions | `--mask-transitions` |
| Trim 7 bases from each end | `--trim-ends 7` |

These runs used an earlier version of the scripts, so they may not reproduce the paper's
numbers exactly. As the paper shows, none of these treatments corrects damage adequately.

## Scripts reference

| Script | Step | Role |
|--------|------|------|
| `maf_to_fasta.py` | reference MSA | Stitches a Cactus MAF into one FASTA alignment with one row per genome. |
| `filter_compartments.py` | 2 | Splits remapped reads by compartment and removes outgroup, multimapping and low-scoring reads. |
| `filter_msa_by_sam.py` | 3 | Trims the MSA to read-covered columns and converts reads to MSA coordinates. |
| `make_reference.R` | 5 | Builds the ultrametric tree, GTR+Γ parameters and ratePlacer reference files. |
| `generate_matrices.py` | 7 | Builds the A, C and D matrices for damage estimation. |
| `run_analysis.py` | 8 | Estimates damage (B at both ends of the α interval). |
| `damage_estimator.py`, `profile_likelihood.py` | 8 | Modules imported by `run_analysis.py`. |
| `generateRatePlacerInputs_msaReads_damageMatrix.py` | 9 | Writes the ratePlacer alignment, assignment and error profile files. |
