#!/usr/bin/env python3
"""
Generate matrices A, C, and D for ancient DNA damage analysis.

The framework models: D = A(α) · B · C, where:
  A(α) = expm(α·Q)  — evolutionary substitution (this script outputs Q)
  B               — DNA damage (estimated downstream from D, A, C)
  C               — sequencing error
  D               — observed positional mismatch matrix
All are row-stochastic, M[from, to], so the product reads left to right in the physical
order: reference → evolution (A) → damage (B) → sequencing error (C) → read.

Outputs:
  A.npy        : normalized GTR rate matrix Q (shape 4×4); A(α) = expm(α·Q)
  C.npy        : sequencing error matrix (shape 4×4)
  D.npy        : positional mismatch matrix (shape 2N+1 × 4 × 4)
  ref_freq.npy   : positional reference base frequencies (shape 2N+1 × 4)
  ref_counts.npy : positional reference base raw counts (shape 2N+1 × 4, int64)
                 same positional indexing as D:
                   indices 0..N-1  = 5' positions 1..N
                   index  N        = middle (averaged)
                   indices N+1..2N = 3' positions N..1 (2N = 3'-most base)

All matrices use row/col order A/C/G/T.
"""

import argparse
import os
import sys

import numpy as np
import pysam


BASES = 'ACGT'
BASE_IDX = {b: i for i, b in enumerate(BASES)}
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Exchange rate pair indices — order in rate_matrix: [AC, AG, AT, CG, CT, GT]
PAIR_TO_IDX = {
    ('A', 'C'): 0, ('C', 'A'): 0,
    ('A', 'G'): 1, ('G', 'A'): 1,
    ('A', 'T'): 2, ('T', 'A'): 2,
    ('C', 'G'): 3, ('G', 'C'): 3,
    ('C', 'T'): 4, ('T', 'C'): 4,
    ('G', 'T'): 5, ('T', 'G'): 5,
}


def complement(base):
    return COMP.get(base, 'N')


def parse_gtr_parameters(path):
    """
    Parse GTR parameter file (<tree>_parameter.txt from make_reference.R).
    Read as a sequence of 11 numbers, as ratePlacer does, so line breaks don't matter:
      gamma shape (ignored), base frequencies A C G T, exchange rates AC AG AT CG CT GT
    Returns dict with 'base_frequencies' and 'rate_matrix' (list of 6 floats).
    """
    with open(path) as f:
        values = [float(x) for x in f.read().split()]
    if len(values) != 11:
        sys.exit(f"Error: expected 11 numbers in GTR parameter file {path}, found {len(values)}")

    base_frequencies = dict(zip(BASES, values[1:5]))
    rate_matrix = values[5:11]

    return {'base_frequencies': base_frequencies, 'rate_matrix': rate_matrix}


def build_A_matrix(gtr_params):
    """
    Build the normalized GTR rate matrix Q.
    Off-diagonal: Q[i][j] = exchange_rate(i,j) × π[j]
    Diagonal:     Q[i][i] = -sum_{j≠i} Q[i][j]
    Normalized so that -sum_i π[i]·Q[i][i] = 1 (expected substitution rate = 1).
    """
    pi = gtr_params['base_frequencies']
    rates = gtr_params['rate_matrix']

    Q = np.zeros((4, 4))
    for i, bi in enumerate(BASES):
        for j, bj in enumerate(BASES):
            if i != j:
                Q[i, j] = rates[PAIR_TO_IDX[(bi, bj)]] * pi[bj]
        Q[i, i] = -Q[i, :].sum()

    pi_arr = np.array([pi[b] for b in BASES])
    norm_factor = -np.dot(pi_arr, np.diag(Q))
    if norm_factor > 0:
        Q /= norm_factor

    return Q


def _row_normalize(counts):
    """Row-normalize a 4×4 count matrix. Zero rows become uniform (0.25)."""
    mat = counts.astype(np.float64)
    for i in range(4):
        s = mat[i].sum()
        if s > 0:
            mat[i] /= s
        else:
            mat[i] = 0.25
    return mat


def scan_alignment(input_file, n_positions, min_baseq=0, min_mapq=0):
    """
    Single pass over SAM/BAM to collect data for C and D.

    For each aligned base (MD-tag derived, no FASTA needed):
      - Strand-correct to original molecule orientation (reverse reads:
        complement both ref and read bases, flip position from 5' end).
      - Route (ref_base, read_base) to the appropriate D count array.
      - For C, bucket the Phred quality score by the molecule-orientation base.

    Returns:
      quality_scores  : dict[base] -> list of Phred scores (molecule orientation)
      five_counts     : ndarray (n_positions, 4, 4) — 5' end raw counts
      three_counts    : ndarray (n_positions, 4, 4) — 3' end raw counts
                        index 0 = 3'-most position
      middle_counts   : ndarray (4, 4) — middle raw counts
      reads_processed : int — total mapped primary reads
      total_read_length : int — sum of query_length across all processed reads
    """
    quality_scores = {b: [] for b in BASES}
    five_counts  = np.zeros((n_positions, 4, 4), dtype=np.float64)
    three_counts = np.zeros((n_positions, 4, 4), dtype=np.float64)
    middle_counts = np.zeros((4, 4), dtype=np.float64)

    aln = pysam.AlignmentFile(input_file, "r")
    reads_processed = 0
    reads_reverse = 0
    bases_processed = 0
    total_read_length = 0
    md_error_shown = False

    print("Processing alignment file...", file=sys.stderr)

    for read in aln:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue

        reads_processed += 1
        if reads_processed % 100_000 == 0:
            print(f"  {reads_processed:,} reads  {bases_processed:,} bases...",
                  file=sys.stderr)

        is_reverse = read.is_reverse
        if is_reverse:
            reads_reverse += 1
        total_read_length += read.query_length

        qual_array = read.query_qualities
        read_seq = read.query_sequence
        if read_seq is None:
            continue
        read_seq = read_seq.upper()
        read_length = read.query_length

        for query_pos, _ref_pos, ref_base in read.get_aligned_pairs(
                matches_only=True, with_seq=True):

            if ref_base is None:
                if not md_error_shown:
                    print(
                        "\nERROR: MD tag missing or incomplete — ref_base is None.\n"
                        "Add MD tags with:\n"
                        "  samtools calmd -b in.bam ref.fa > out.bam",
                        file=sys.stderr,
                    )
                    aln.close()
                    sys.exit(1)
                continue

            ref_base  = ref_base.upper()
            read_base = read_seq[query_pos]

            if ref_base not in BASES or read_base not in BASES:
                continue

            if qual_array is not None and qual_array[query_pos] < min_baseq:
                continue

            bases_processed += 1

            # -- C matrix: quality score bucketed to molecule-orientation base --
            mol_read_base = complement(read_base) if is_reverse else read_base
            if qual_array is not None:
                quality_scores[mol_read_base].append(qual_array[query_pos])

            # -- D matrix: molecule-orientation bases and 5'-anchored position --
            if is_reverse:
                mol_ref  = complement(ref_base)
                mol_read = complement(read_base)
                mol_pos  = read_length - 1 - query_pos
            else:
                mol_ref  = ref_base
                mol_read = read_base
                mol_pos  = query_pos

            ri = BASE_IDX[mol_ref]
            oi = BASE_IDX[mol_read]

            dist_from_3prime = read_length - 1 - mol_pos
            if mol_pos < n_positions:
                five_counts[mol_pos, ri, oi] += 1
            elif dist_from_3prime < n_positions:
                three_counts[dist_from_3prime, ri, oi] += 1
            else:
                middle_counts[ri, oi] += 1

    aln.close()
    print(
        f"\nProcessed {reads_processed:,} reads ({reads_reverse:,} reverse), "
        f"{bases_processed:,} bases",
        file=sys.stderr,
    )
    return quality_scores, five_counts, three_counts, middle_counts, reads_processed, total_read_length


def build_C_matrix(quality_scores):
    """
    Build sequencing error matrix C from Phred quality scores.
      C[i][i] = 1 - avg_err[i]
      C[i][j] = avg_err[i] / 3  for j ≠ i
    where avg_err[i] = mean(10^(-Q/10)) over all bases of molecule type i.
    """
    C = np.zeros((4, 4))
    for i, base in enumerate(BASES):
        scores = quality_scores[base]
        if len(scores) == 0:
            C[i, i] = 1.0
        else:
            avg_err = np.mean(10.0 ** (-np.array(scores, dtype=np.float64) / 10.0))
            C[i, i] = 1.0 - avg_err
            for j in range(4):
                if j != i:
                    C[i, j] = avg_err / 3.0
    return C


def build_ref_counts_matrix(five_counts, three_counts, middle_counts, n_positions):
    """
    Build positional reference base count array of shape (2N+1, 4).
    Each row is the raw count of [A, C, G, T] reference bases at that position,
    using the same positional indexing as D.
    """
    n = n_positions
    ref_counts = np.zeros((2 * n + 1, 4), dtype=np.int64)

    for p in range(n):
        ref_counts[p]                  = five_counts[p].sum(axis=1)
        ref_counts[n + 1 + (n - 1 - p)] = three_counts[p].sum(axis=1)

    ref_counts[n] = middle_counts.sum(axis=1)
    return ref_counts


def build_ref_freq_matrix(five_counts, three_counts, middle_counts, n_positions):
    """
    Build positional reference base frequency array of shape (2N+1, 4).
    Each row is the frequency of [A, C, G, T] among reference bases at that position,
    using the same positional indexing as D.
    """
    n = n_positions
    ref_freq = np.zeros((2 * n + 1, 4))

    def _normalize_row(counts_4x4):
        row = counts_4x4.sum(axis=1)  # sum over observed bases → counts per ref base
        total = row.sum()
        return row / total if total > 0 else np.full(4, 0.25)

    for p in range(n):
        ref_freq[p]             = _normalize_row(five_counts[p])
        ref_freq[n + 1 + (n - 1 - p)] = _normalize_row(three_counts[p])

    ref_freq[n] = _normalize_row(middle_counts)
    return ref_freq


def build_D_matrix(five_counts, three_counts, middle_counts, n_positions):
    """
    Assemble the (2N+1, 4, 4) normalized D array.
      D[0..N-1]    : 5' positions 1..N   (D[0] = 5'-most base)
      D[N]         : middle
      D[N+1..2N]   : 3' positions N..1  (D[2N] = 3'-most base)
    """
    n = n_positions
    D = np.zeros((2 * n + 1, 4, 4))

    for p in range(n):
        D[p] = _row_normalize(five_counts[p])

    D[n] = _row_normalize(middle_counts)

    # three_counts[0] = 3'-most base → D[2N]
    # three_counts[N-1] = position N from 3' end → D[N+1]
    for p in range(n):
        D[n + 1 + (n - 1 - p)] = _row_normalize(three_counts[p])

    return D


def _print_matrix(label, mat):
    """Print a labeled 4×4 matrix with A/C/G/T row and column headers."""
    print(f"\n{label}")
    print(f"{'':>6}" + "".join(f"{c:>10}" for c in BASES))
    for i, row_base in enumerate(BASES):
        row = f"{row_base:>6}" + "".join(f"{mat[i, j]:>10.6f}" for j in range(4))
        print(row)


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Generate matrices A (GTR rate Q), C (sequencing error), and '
            'D (positional mismatch) for ancient DNA damage analysis.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('input_file',
                        help='Input SAM or BAM file (must have MD tags)')
    parser.add_argument('--params', metavar='FILE', default=None,
                        help='GTR parameter file '
                             '(default: 0_parameter.txt one directory above this script)')
    parser.add_argument('--n-positions', type=int, default=10, metavar='INT',
                        help='Positions at each end for D matrix (default: 10)')
    parser.add_argument('--output-dir', metavar='DIR', default='.',
                        help='Directory for .npy output files (default: .)')
    parser.add_argument('--prefix', metavar='STR', default='',
                        help='Prefix for output filenames, e.g. "ds_" → ds_A.npy (default: none)')
    parser.add_argument('--min-baseq', type=int, default=0, metavar='INT',
                        help='Minimum base quality to include (default: 0)')
    parser.add_argument('--min-mapq', type=int, default=0, metavar='INT',
                        help='Minimum mapping quality to include (default: 0)')

    args = parser.parse_args()

    if args.params is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        args.params = os.path.join(script_dir, '..', '0_parameter.txt')

    # --- Matrix A ---
    print(f"Loading GTR parameters from: {args.params}", file=sys.stderr)
    gtr_params = parse_gtr_parameters(args.params)
    bf = gtr_params['base_frequencies']
    rm = gtr_params['rate_matrix']
    print(f"  Base frequencies: A={bf['A']:.6f} C={bf['C']:.6f} "
          f"G={bf['G']:.6f} T={bf['T']:.6f}", file=sys.stderr)
    print(f"  Rate matrix (AC AG AT CG CT GT): "
          + " ".join(f"{r:.6f}" for r in rm), file=sys.stderr)
    A = build_A_matrix(gtr_params)

    # --- Scan alignment (C and D data) ---
    quality_scores, five_counts, three_counts, middle_counts, reads_processed, total_read_length = scan_alignment(
        args.input_file,
        n_positions=args.n_positions,
        min_baseq=args.min_baseq,
        min_mapq=args.min_mapq,
    )

    # --- Matrix C ---
    C = build_C_matrix(quality_scores)

    # --- Matrix D and reference base frequencies ---
    D = build_D_matrix(five_counts, three_counts, middle_counts, args.n_positions)
    ref_freq   = build_ref_freq_matrix(five_counts, three_counts, middle_counts, args.n_positions)
    ref_counts = build_ref_counts_matrix(five_counts, three_counts, middle_counts, args.n_positions)

    # --- Print read summary ---
    avg_len = total_read_length / reads_processed if reads_processed > 0 else 0.0
    print("\nRead summary")
    print(f"  Total reads  : {reads_processed:,}")
    print(f"  Average length : {avg_len:.1f} bp")

    # --- Print matrices ---
    n = args.n_positions
    _print_matrix("A  (normalized GTR rate matrix Q;  A(α) = expm(α·Q))", A)
    _print_matrix("C  (sequencing error matrix)", C)
    for p in range(n):
        _print_matrix(f"D  5' position {p + 1}", D[p])
    _print_matrix("D  Middle", D[n])
    for p in range(n):
        dist = n - p  # n down to 1
        _print_matrix(f"D  3' position {dist}", D[n + 1 + p])

    print(f"\nRef base frequencies  (shape {ref_freq.shape},  same positional indexing as D)")
    print(f"{'':>12}" + "".join(f"{b:>10}" for b in BASES))
    for p in range(n):
        print(f"{'5pr pos '+str(p+1):>12}" + "".join(f"{ref_freq[p, j]:>10.6f}" for j in range(4)))
    print(f"{'Middle':>12}" + "".join(f"{ref_freq[n, j]:>10.6f}" for j in range(4)))
    for p in range(n):
        dist = n - p
        print(f"{'3pr pos '+str(dist):>12}" + "".join(f"{ref_freq[n + 1 + p, j]:>10.6f}" for j in range(4)))

    print(f"\nRef base counts  (shape {ref_counts.shape},  same positional indexing as D)")
    print(f"{'':>12}" + "".join(f"{b:>12}" for b in BASES))
    for p in range(n):
        print(f"{'5pr pos '+str(p+1):>12}" + "".join(f"{ref_counts[p, j]:>12,}" for j in range(4)))
    print(f"{'Middle':>12}" + "".join(f"{ref_counts[n, j]:>12,}" for j in range(4)))
    for p in range(n):
        dist = n - p
        print(f"{'3pr pos '+str(dist):>12}" + "".join(f"{ref_counts[n + 1 + p, j]:>12,}" for j in range(4)))

    # --- Save ---
    os.makedirs(args.output_dir, exist_ok=True)
    a_path   = os.path.join(args.output_dir, f'{args.prefix}A.npy')
    c_path   = os.path.join(args.output_dir, f'{args.prefix}C.npy')
    d_path   = os.path.join(args.output_dir, f'{args.prefix}D.npy')
    ref_path    = os.path.join(args.output_dir, f'{args.prefix}ref_freq.npy')
    counts_path = os.path.join(args.output_dir, f'{args.prefix}ref_counts.npy')
    np.save(a_path, A)
    np.save(c_path, C)
    np.save(d_path, D)
    np.save(ref_path, ref_freq)
    np.save(counts_path, ref_counts)
    print(f"\nSaved A.npy          →  {a_path}", file=sys.stderr)
    print(f"Saved C.npy          →  {c_path}", file=sys.stderr)
    print(f"Saved D.npy          →  {d_path}  shape: {D.shape}", file=sys.stderr)
    print(f"Saved ref_freq.npy   →  {ref_path}  shape: {ref_freq.shape}", file=sys.stderr)
    print(f"Saved ref_counts.npy →  {counts_path}  shape: {ref_counts.shape}", file=sys.stderr)
    print(
        f"  D index guide: 0..{n-1} = 5' pos 1..{n}, "
        f"{n} = middle, "
        f"{n+1}..{2*n} = 3' pos {n}..1",
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
