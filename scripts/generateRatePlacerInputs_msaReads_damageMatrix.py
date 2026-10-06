import argparse
import gzip
import math
import os

parser = argparse.ArgumentParser(description='Generate RatePlacer inputs from MSA-aligned reads with deamination masking and error profiles')
parser.add_argument("-o", "--outMap", required=True, help='Output read mapping file')
parser.add_argument("-i", "--inputReads", required=True, help='Input FASTQ file from filter_msa_by_sam.py --extract-reads (header: start, len, ref)')
parser.add_argument("-n", "--node", type=int, required=True, help='ratePlacer node number for assignment (all reads assigned to this node, tree=0)')
parser.add_argument("-sam", "--samFile", required=False, help='Optional SAM file for reverse complement flags (fallback if FASTQ headers lack rc= field)')
parser.add_argument("-f", "--taxaFilter", required=False, help='Optional taxa filter file - subset of taxa file containing acceptable reference IDs')
parser.add_argument("-b", "--outAlign", required=True, help='Output alignment file')
parser.add_argument("-c", "--outAssign", required=True, help='Output assignment file')
parser.add_argument("-e", "--outError", required=True, help='Output error profile file')
parser.add_argument("-m", "--refMSA", required=True, help='Reference MSA file')
parser.add_argument("-g", "--gtrParams", required=False, help='GTR model parameter file containing base frequencies and rate matrix parameters')
parser.add_argument("--error-fudge-factor", type=float, default=0.0,
                   help='Optional fudge/error factor (e) added to error rate calculation: error_rate = ref_norm - expected_rate × correction_factor + e. Default is 0.0 (no adjustment). Use to account for potential unknown errors in the correction factor.')
parser.add_argument("--positional-matrix-bases", type=int, default=15,
                   help='Number of bases from 5\' and 3\' ends to use separate positional rate matrices (default: 15). Set to 0 to use single global matrix.')
parser.add_argument("--quality-score-only", action='store_true',
                   help='Use only base quality scores for error profiles, skip phylogenetic error rate analysis')
parser.add_argument("--encoding", choices=['phred33', 'phred64'], default='phred33',
                   help='Quality score encoding (default: phred33)')
parser.add_argument("--error-profile-bases", type=int, default=None,
                   help='Only generate error profiles for first and last N bases of each read (default: all bases)')
parser.add_argument("--no-error-profile", action='store_true',
                   help='Disable error profile generation')
parser.add_argument("--mask-transitions", action='store_true',
                   help='Mask all transition substitutions (A<->G, C<->T), keeping only transversions. Applied after trimming')
parser.add_argument("--max-mismatches", type=int, default=5,
                   help='Maximum number of allowed mismatches to reference after masking. Applied after trimming (default: 5)')
parser.add_argument("--substitution-analysis", action='store_true',
                   help='Perform substitution analysis comparing filtered reads to reference sequences')
parser.add_argument("--output-mismatched-fasta", type=str, default=None,
                   help='Output file path to write reads with mismatches in FASTA format after masking')
parser.add_argument("--damage-masking", choices=['none', 'ends', 'all'], default='none',
                   help='Mask damage patterns: "none" = no masking, "ends" = mask at sequence ends, "all" = mask across full sequence. Applied after trimming (default: none)')
parser.add_argument("--damage-masking-bases", type=int, default=15,
                   help='Number of bases from each end to mask when using --damage-masking ends. Applied after trimming (default: 15)')
parser.add_argument("--damage-types", choices=['deamination', 'artifacts', 'both'], default='deamination',
                   help='Types of damage to mask: "deamination" = C→T and G→A only (default), "artifacts" = G→T and C→A only, "both" = all four patterns')
parser.add_argument("--debug", action='store_true',
                   help='Enable debug output messages (default: disabled)')
parser.add_argument("--show-read-examples", action='store_true',
                   help='Show detailed examples of forward and reverse complement reads for analysis (default: disabled)')
parser.add_argument("--trim-ends", type=int, default=0,
                   help='Number of bases to trim from both 5\' and 3\' ends of reads (default: 0, no trimming)')
parser.add_argument("--damage-matrix", type=str, default=None,
                   help='Path to pre-computed .npy damage probability matrix (shape: 2N+1, 4, 4). '
                        'N is auto-detected from shape. Replaces on-the-fly obs/exp damage calculation. '
                        'Skips --substitution-analysis / --gtrParams damage estimation.')

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--singlestrand', action='store_true', help='Single-stranded DNA (mask T->C damage)')
group.add_argument('--doublestrand', action='store_true', help='Double-stranded DNA (mask A->G damage)')

args = parser.parse_args()

def parse_gtr_parameters(gtr_file_path):
    """
    Parse GTR model parameters from parameter file.

    Args:
        gtr_file_path: Path to the GTR parameter file

    Returns:
        dict: Dictionary containing 'base_frequencies' and 'rate_matrix' parameters
              Returns None if file not provided or parsing fails
    """
    if not gtr_file_path:
        return None

    try:
        # Read as a sequence of numbers (as ratePlacer does), so line breaks don't matter:
        # gamma shape, base frequencies A C G T, exchange rates AC AG AT CG CT GT
        with open(gtr_file_path, 'r') as f:
            values = [float(x) for x in f.read().split()]
        if len(values) != 11:
            raise ValueError(f"expected 11 numbers, found {len(values)}")

        base_frequencies = {
            'A': values[1],
            'C': values[2],
            'G': values[3],
            'T': values[4]
        }
        rate_matrix = values[5:11]

        if args.debug:
            print(f"\nGTR Parameters loaded from {gtr_file_path}:")
            print(f"Base frequencies: {base_frequencies}")
            print(f"Rate matrix (AC AG AT CG CT GT): {rate_matrix}")

        return {
            'base_frequencies': base_frequencies,
            'rate_matrix': rate_matrix
        }

    except (FileNotFoundError, ValueError) as e:
        raise SystemExit(f"Error parsing GTR parameter file {gtr_file_path}: {e}")

def load_damage_matrix_from_npy(path):
    """
    Load a pre-computed damage probability matrix from a .npy file.

    Expected shape: (2N+1, 4, 4) where N is the number of end positions covered.
    Axis 0: positional index (0..N-1 = 5' positions, N = middle, N+1..2N = 3' positions)
    Axis 1: reference/latent base row  (A=0, C=1, G=2, T=3)
    Axis 2: damaged/intermediate base column (A=0, C=1, G=2, T=3)

    Returns:
        (damage_npy, N) — the loaded array and the number of end positions
    """
    import numpy as np
    mat = np.load(path)
    if mat.ndim != 3 or mat.shape[1] != 4 or mat.shape[2] != 4:
        raise ValueError(f"Damage matrix must have shape (2N+1, 4, 4), got {mat.shape}")
    total_pos = mat.shape[0]
    if total_pos % 2 == 0:
        raise ValueError(f"Damage matrix first dimension must be odd (2N+1), got {total_pos}")
    n = (total_pos - 1) // 2
    print(f"Loaded damage matrix from {path}: shape={mat.shape}, N={n} end positions")
    return mat, n

def build_positional_matrices_from_npy(damage_npy, n):
    """
    Convert a (2N+1, 4, 4) damage matrix to the positional_matrices_arrays dict
    expected by generate_error_profile().

    Positional index mapping (matches generate_matrices.py D-matrix ordering):
      damage_npy[i]       -> ('5prime', i)     for i in 0..N-1
                             (i=0 is the 5'-most base of the molecule)
      damage_npy[N]       -> 'middle'
      damage_npy[2N - j]  -> ('3prime', j)     for j in 0..N-1
                             (j=0 is the 3'-most base; j=N-1 is Nth from 3')

    Args:
        damage_npy: numpy array of shape (2N+1, 4, 4)
        n: number of end positions (N)

    Returns:
        Dict with numpy array values keyed by ('5prime', i), 'middle', ('3prime', j)
    """
    matrices = {}
    for i in range(n):
        matrices[('5prime', i)] = damage_npy[i]
    matrices['middle'] = damage_npy[n]
    for j in range(n):
        matrices[('3prime', j)] = damage_npy[2 * n - j]
    return matrices

def get_gtr_rate_for_substitution(sub_type, gtr_params):
    """
    Get the GTR substitution rate for a specific substitution type.

    In GTR model: rate(i→j) = rate_matrix[ij] × freq[j]

    Args:
        sub_type: Substitution type (e.g., 'A→C')
        gtr_params: GTR parameters dictionary

    Returns:
        GTR substitution rate for the substitution, or None if not available
    """
    if not gtr_params:
        return None

    # GTR rate matrix order: [AC, AG, AT, CG, CT, GT]
    rate_matrix = gtr_params['rate_matrix']
    base_frequencies = gtr_params['base_frequencies']

    # Map substitution types to GTR rate calculations
    # GTR substitution rate = rate_matrix[ij] × freq[target_base]
    substitution_calculations = {
        'A→C': rate_matrix[0] * base_frequencies['C'],  # AC rate × freq(C)
        'C→A': rate_matrix[0] * base_frequencies['A'],  # CA rate × freq(A)
        'A→G': rate_matrix[1] * base_frequencies['G'],  # AG rate × freq(G)
        'G→A': rate_matrix[1] * base_frequencies['A'],  # GA rate × freq(A)
        'A→T': rate_matrix[2] * base_frequencies['T'],  # AT rate × freq(T)
        'T→A': rate_matrix[2] * base_frequencies['A'],  # TA rate × freq(A)
        'C→G': rate_matrix[3] * base_frequencies['G'],  # CG rate × freq(G)
        'G→C': rate_matrix[3] * base_frequencies['C'],  # GC rate × freq(C)
        'C→T': rate_matrix[4] * base_frequencies['T'],  # CT rate × freq(T)
        'T→C': rate_matrix[4] * base_frequencies['C'],  # TC rate × freq(C)
        'G→T': rate_matrix[5] * base_frequencies['T'],  # GT rate × freq(T)
        'T→G': rate_matrix[5] * base_frequencies['G'],  # TG rate × freq(G)
    }

    return substitution_calculations.get(sub_type, None)

def get_expected_gtr_rate_from_data(sub_type, gtr_params, ref_base_counts):
    """
    Calculate expected GTR substitution rate using data-based base frequencies.

    In GTR model: rate(i→j) = rate_matrix[ij] × freq[j]
    Here freq[j] is calculated from the actual reference base counts in the data.

    Args:
        sub_type: Substitution type (e.g., 'A→C')
        gtr_params: GTR parameters dictionary (for rate matrix)
        ref_base_counts: Dictionary of reference base counts from the data

    Returns:
        Expected GTR substitution rate, or None if not available
    """
    if not gtr_params or not ref_base_counts:
        return None

    # Calculate base frequencies from data
    total_bases = sum(ref_base_counts.values())
    if total_bases == 0:
        print(f"DEBUG: get_expected_gtr_rate_from_data - total_bases is 0 for {sub_type}")
        return None

    data_base_frequencies = {
        base: count / total_bases
        for base, count in ref_base_counts.items()
    }

    # GTR rate matrix order: [AC, AG, AT, CG, CT, GT]
    rate_matrix = gtr_params['rate_matrix']

    # Map substitution types to GTR rate calculations using data-based frequencies
    # GTR substitution rate = rate_matrix[ij] × data_freq[target_base]
    target_base = sub_type.split('→')[1]

    substitution_to_matrix_index = {
        'A→C': 0, 'C→A': 0,  # AC
        'A→G': 1, 'G→A': 1,  # AG
        'A→T': 2, 'T→A': 2,  # AT
        'C→G': 3, 'G→C': 3,  # CG
        'C→T': 4, 'T→C': 4,  # CT
        'G→T': 5, 'T→G': 5,  # GT
    }

    matrix_idx = substitution_to_matrix_index.get(sub_type)
    if matrix_idx is None:
        print(f"DEBUG: get_expected_gtr_rate_from_data - no matrix index for {sub_type}")
        return None

    target_freq = data_base_frequencies.get(target_base, 0)
    expected_rate = rate_matrix[matrix_idx] * target_freq

    # # Debug first few calculations
    # import random
    # if random.random() < 0.1:  # Print ~10% of calculations
    #     print(f"DEBUG: {sub_type}: matrix[{matrix_idx}]={rate_matrix[matrix_idx]:.6f} × freq[{target_base}]={target_freq:.6f} = {expected_rate:.6f}")

    return expected_rate

def print_gtr_substitution_matrix(gtr_params):
    """
    Print the full GTR substitution rate matrix for verification.

    Args:
        gtr_params: GTR parameters dictionary
    """
    if not gtr_params:
        print("No GTR parameters available to display matrix")
        return

    print("\nGTR SUBSTITUTION RATE MATRIX")
    print("=" * 60)
    print("Formula: rate(i→j) = rate_matrix[ij] × freq[j]")
    print()

    # Print base frequencies
    base_frequencies = gtr_params['base_frequencies']
    print("Base frequencies:")
    for base in ['A', 'C', 'G', 'T']:
        print(f"  freq[{base}] = {base_frequencies[base]:.6f}")
    print()

    # Print rate matrix parameters
    rate_matrix = gtr_params['rate_matrix']
    print("Rate matrix parameters:")
    print(f"  AC/CA = {rate_matrix[0]:.6f}")
    print(f"  AG/GA = {rate_matrix[1]:.6f}")
    print(f"  AT/TA = {rate_matrix[2]:.6f}")
    print(f"  CG/GC = {rate_matrix[3]:.6f}")
    print(f"  CT/TC = {rate_matrix[4]:.6f}")
    print(f"  GT/TG = {rate_matrix[5]:.6f}")
    print()

    # Calculate and print full substitution matrix
    bases = ['A', 'C', 'G', 'T']
    substitution_types = [
        'A→C', 'A→G', 'A→T',
        'C→A', 'C→G', 'C→T',
        'G→A', 'G→C', 'G→T',
        'T→A', 'T→C', 'T→G'
    ]

    print("Substitution rates:")
    print("-" * 40)
    for sub_type in substitution_types:
        rate = get_gtr_rate_for_substitution(sub_type, gtr_params)
        print(f"  {sub_type}: {rate:.6f}")
    print()

    # Print as matrix format
    print("Matrix format (from → to):")
    print("-" * 40)
    print(f"{'':>8} {'A':>10} {'C':>10} {'G':>10} {'T':>10}")
    print("-" * 48)

    for i, from_base in enumerate(bases):
        row = f"{from_base:>8}"
        for j, to_base in enumerate(bases):
            if from_base == to_base:
                row += f"{'—':>10}"  # Diagonal (no substitution)
            else:
                sub_type = f"{from_base}→{to_base}"
                rate = get_gtr_rate_for_substitution(sub_type, gtr_params)
                row += f"{rate:>10.6f}"
        print(row)
    print()

def calculate_error_rate_matrix(dataset, gtr_params, fudge_factor=0.0):
    """
    Calculate error rate matrix using the substitution with smallest obs/exp ratio as baseline.

    Formula: error_rate[ref→obs] = ref_norm[ref→obs] - expected_rate[ref→obs] × (ref_norm_small / expected_rate_small) + fudge_factor

    The obs/exp ratio is calculated as:
    - obs = ref_norm (observed substitution rate normalized by reference base count)
    - exp = GTR rate calculated using rate matrix and data-based base frequencies

    Args:
        dataset: Substitution data for analysis
        gtr_params: GTR parameters dictionary
        fudge_factor: Optional error adjustment factor (default: 0.0)

    Returns:
        Dictionary containing error rates and baseline information
    """
    if not gtr_params:
        return None

    # Calculate obs/exp ratios for all substitutions
    obs_exp_ratios = {}
    ref_norms = {}
    expected_rates = {}

    # Debug: print reference base counts
    print(f"DEBUG: Reference base counts for obs/exp calculation:")
    for base in ['A', 'C', 'G', 'T']:
        count = dataset['ref_base_counts'].get(base, 0)
        total = sum(dataset['ref_base_counts'].values())
        freq = count / total if total > 0 else 0
        print(f"  {base}: {count} ({freq:.4f})")

    for sub_type, count in dataset['substitutions'].items():
        ref_base = sub_type.split('→')[0]
        ref_count = dataset['ref_base_counts'].get(ref_base, 0)

        if ref_count > 0:
            # Observed rate (normalized by reference base count)
            ref_norm = count / ref_count
            # Expected rate using GTR matrix with data-based base frequencies
            expected_rate = get_expected_gtr_rate_from_data(sub_type, gtr_params, dataset['ref_base_counts'])

            if expected_rate and expected_rate > 0:
                obs_exp_ratio = ref_norm / expected_rate
                obs_exp_ratios[sub_type] = obs_exp_ratio
                ref_norms[sub_type] = ref_norm
                expected_rates[sub_type] = expected_rate
            else:
                print(f"DEBUG: Skipping {sub_type} - expected_rate is {expected_rate}")

    print(f"DEBUG: Total obs_exp_ratios calculated: {len(obs_exp_ratios)} out of {len(dataset['substitutions'])} substitutions")

    if not obs_exp_ratios:
        print("DEBUG: No valid obs/exp ratios calculated - returning None from calculate_error_rate_matrix")
        return None

    # Find substitution with smallest obs/exp ratio
    smallest_sub = min(obs_exp_ratios, key=obs_exp_ratios.get)
    smallest_ratio = obs_exp_ratios[smallest_sub]
    ref_norm_small = ref_norms[smallest_sub]
    expected_rate_small = expected_rates[smallest_sub]

    # Calculate correction factor
    correction_factor = ref_norm_small / expected_rate_small

    # Calculate error rates for all substitutions
    error_rates = {}
    for sub_type in obs_exp_ratios:
        ref_norm = ref_norms[sub_type]
        expected_rate = expected_rates[sub_type]

        # Formula: error_rate = ref_norm - expected_rate × correction_factor + fudge_factor
        error_rate = ref_norm - (expected_rate * correction_factor) + fudge_factor
        error_rates[sub_type] = error_rate

    return {
        'smallest_substitution': smallest_sub,
        'smallest_ratio': smallest_ratio,
        'ref_norm_small': ref_norm_small,
        'expected_rate_small': expected_rate_small,
        'correction_factor': correction_factor,
        'fudge_factor': fudge_factor,
        'error_rates': error_rates,
        'ref_norms': ref_norms,
        'expected_rates': expected_rates,
        'obs_exp_ratios': obs_exp_ratios
    }

### TO DO: IS THIS THE EXACT SAME AS calculate_error_rate_matrix?
def calculate_error_rate_matrix_from_data(data, gtr_params, fudge_factor=0.0):
    """
    Calculate error rate matrix from positional substitution data using smallest obs/exp ratio.

    Formula: error_rate[ref→obs] = ref_norm[ref→obs] - expected_rate[ref→obs] × (ref_norm_small / expected_rate_small) + fudge_factor

    The obs/exp ratio is calculated as:
    - obs = ref_norm (observed substitution rate normalized by reference base count)
    - exp = GTR rate calculated using rate matrix and data-based base frequencies

    Args:
        data: Dictionary with 'substitutions' (Counter) and 'ref_base_counts' (Counter)
        gtr_params: GTR parameters from file
        fudge_factor: Optional error adjustment factor (default: 0.0)

    Returns:
        Dictionary with 'error_rates' key containing error rates for each substitution type
    """
    if not gtr_params:
        return None

    # Calculate obs/exp ratios for all substitutions
    obs_exp_ratios = {}
    ref_norms = {}
    expected_rates = {}

    for sub_type, count in data['substitutions'].items():
        ref_base = sub_type.split('→')[0]
        ref_count = data['ref_base_counts'].get(ref_base, 0)

        if ref_count > 0:
            # Observed rate (normalized by reference base count)
            ref_norm = count / ref_count
            # Expected rate using GTR matrix with data-based base frequencies
            expected_rate = get_expected_gtr_rate_from_data(sub_type, gtr_params, data['ref_base_counts'])

            if expected_rate and expected_rate > 0:
                obs_exp_ratio = ref_norm / expected_rate
                obs_exp_ratios[sub_type] = obs_exp_ratio
                ref_norms[sub_type] = ref_norm
                expected_rates[sub_type] = expected_rate

    if not obs_exp_ratios:
        return None

    # Find substitution with smallest obs/exp ratio
    smallest_sub = min(obs_exp_ratios, key=obs_exp_ratios.get)
    smallest_ratio = obs_exp_ratios[smallest_sub]
    ref_norm_small = ref_norms[smallest_sub]
    expected_rate_small = expected_rates[smallest_sub]

    # Calculate correction factor
    correction_factor = ref_norm_small / expected_rate_small

    # Calculate error rates for all substitutions
    error_rates = {}
    for sub_type in obs_exp_ratios:
        ref_norm = ref_norms[sub_type]
        expected_rate = expected_rates[sub_type]

        # Formula: error_rate = ref_norm - expected_rate × correction_factor + fudge_factor
        error_rate = ref_norm - (expected_rate * correction_factor) + fudge_factor
        error_rates[sub_type] = error_rate

    return {
        'smallest_substitution': smallest_sub,
        'smallest_ratio': smallest_ratio,
        'ref_norm_small': ref_norm_small,
        'expected_rate_small': expected_rate_small,
        'correction_factor': correction_factor,
        'fudge_factor': fudge_factor,
        'error_rates': error_rates,
        'ref_norms': ref_norms,
        'expected_rates': expected_rates,
        'obs_exp_ratios': obs_exp_ratios
    }

def create_error_rate_matrix(error_matrix):
    """
    Create a proper rate matrix from error rates where each row sums to zero.

    The diagonal elements are set to -1 × (sum of off-diagonal elements in that row)

    Args:
        error_matrix: Dictionary containing error rates from calculate_error_rate_matrix

    Returns:
        Dictionary containing the rate matrix
    """
    if not error_matrix:
        return None

    bases = ['A', 'C', 'G', 'T']
    rate_matrix = {}

    # Initialize matrix with error rates for off-diagonal elements
    for from_base in bases:
        for to_base in bases:
            if from_base != to_base:
                sub_type = f"{from_base}→{to_base}"
                if sub_type in error_matrix['error_rates']:
                    rate_matrix[sub_type] = error_matrix['error_rates'][sub_type]
                else:
                    rate_matrix[sub_type] = 0.0

    # Calculate diagonal elements to make each row sum to zero
    for from_base in bases:
        # Sum all off-diagonal elements in this row
        row_sum = 0.0
        for to_base in bases:
            if from_base != to_base:
                sub_type = f"{from_base}→{to_base}"
                row_sum += rate_matrix[sub_type]

        # Set diagonal element to negative sum of off-diagonal elements
        diagonal_key = f"{from_base}→{from_base}"
        rate_matrix[diagonal_key] = -row_sum

    # Normalize the matrix so the mean substitution rate equals 1.
    # Using uniform stationary frequencies (1/4 each): mu = -(1/4) * sum_i(Q_ii)
    mu = -sum(rate_matrix[f"{b}→{b}"] for b in bases) / len(bases)
    if mu > 0:
        for key in rate_matrix:
            rate_matrix[key] /= mu

    return rate_matrix

def create_error_prob_matrix(error_matrix):
    """
    Create a probability matrix from error rates where each row sums to one.

    Args:
        error_matrix: Dictionary containing error rates from calculate_error_rate_matrix

    Returns:
        Dictionary containing the probability matrix
    """
    if not error_matrix:
        return None

    bases = ['A', 'C', 'G', 'T']
    rate_matrix = {}

    # Initialize matrix with error rates for off-diagonal elements.
    # Clamp to [0, 1]: negative rates mean "not distinguishable from background"
    # (no damage correction for that type). Allowing negatives pushes the diagonal
    # above 1, which produces p_x_given_y > 1 and positive log likelihoods.
    for from_base in bases:
        for to_base in bases:
            if from_base != to_base:
                sub_type = f"{from_base}→{to_base}"
                if sub_type in error_matrix['error_rates']:
                    rate_matrix[sub_type] = max(0.0, error_matrix['error_rates'][sub_type])
                else:
                    rate_matrix[sub_type] = 0.0

    # Calculate diagonal elements to make each row sum to one.
    for from_base in bases:
        row_sum = 1.0
        for to_base in bases:
            if from_base != to_base:
                sub_type = f"{from_base}→{to_base}"
                row_sum -= rate_matrix[sub_type]
        rate_matrix[f"{from_base}→{from_base}"] = row_sum

        # Set diagonal element to negative sum of off-diagonal elements
        #diagonal_key = f"{from_base}→{from_base}"
        #rate_matrix[diagonal_key] = -row_sum

    # Normalize the matrix so the mean substitution rate equals 1.
    # Using uniform stationary frequencies (1/4 each): mu = -(1/4) * sum_i(Q_ii)
    #mu = -sum(rate_matrix[f"{b}→{b}"] for b in bases) / len(bases)
    #if mu > 0:
    #    for key in rate_matrix:
    #        rate_matrix[key] /= mu

    return rate_matrix


def matrix_exponential(rate_matrix, time=1.0):
    """
    Calculate the matrix exponential exp(rate_matrix * time) to get probabilities.

    Uses numpy's matrix exponentiation if available, otherwise scipy.

    Args:
        rate_matrix: Dictionary containing rate matrix elements
        time: Time parameter (default: 1.0)

    Returns:
        Dictionary containing probability matrix elements
    """
    try:
        import numpy as np
        from scipy.linalg import expm
    except ImportError:
        print("Error: numpy and scipy are required for matrix exponentiation")
        return None

    bases = ['A', 'C', 'G', 'T']

    # Convert rate_matrix dictionary to numpy array
    matrix = np.zeros((4, 4))
    base_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    for i, from_base in enumerate(bases):
        for j, to_base in enumerate(bases):
            rate_key = f"{from_base}→{to_base}"
            if rate_key in rate_matrix:
                matrix[i, j] = rate_matrix[rate_key]

    # Calculate matrix exponential
    prob_matrix = expm(matrix * time)

    # Convert back to dictionary format
    prob_dict = {}
    for i, from_base in enumerate(bases):
        for j, to_base in enumerate(bases):
            prob_key = f"{from_base}→{to_base}"
            prob_dict[prob_key] = prob_matrix[i, j]

    return prob_dict

def quality_to_error_prob(quality_char: str, encoding: str = 'phred33') -> float:
    """
    Convert quality character to error probability.

    Args:
        quality_char: Single quality character
        encoding: Quality score encoding ('phred33' or 'phred64')

    Returns:
        Error probability (0.0 to 1.0)
    """
    if encoding == 'phred33':
        q_score = ord(quality_char) - 33
    elif encoding == 'phred64':
        q_score = ord(quality_char) - 64
    else:
        raise ValueError(f"Unsupported encoding: {encoding}")

    # Convert quality score to error probability
    error_prob = 10 ** (-q_score / 10)
    return error_prob

def apply_transition_masking(fastaSeq, sequence, ref_msa):
    """
    Apply transition masking to the sequence, keeping only transversions.

    Transitions are A<->G and C<->T substitutions.
    Transversions are A<->C, A<->T, G<->C, G<->T substitutions.

    Args:
        fastaSeq: List of bases in the sequence
        sequence: Alignment positions
        ref_msa: Reference MSA sequence

    Returns:
        Modified fastaSeq with transition substitutions masked
    """
    if not args.mask_transitions:
        return fastaSeq

    # Define transition pairs
    transitions = {
        ('A', 'G'), ('G', 'A'),  # A<->G
        ('C', 'T'), ('T', 'C')   # C<->T
    }

    for pos, ref_index in enumerate(sequence):
        if ref_index == -1:  # Skip insertions
            continue

        if ref_index < len(ref_msa):
            ref_base = ref_msa[ref_index].upper()
            read_base = fastaSeq[pos].upper()

            # Check if this is a transition substitution
            if (ref_base, read_base) in transitions:
                fastaSeq[pos] = "-"  # Mask the transition

    return fastaSeq


def apply_damage_masking(fastaSeq, sequence, ref_msa, damage_masking_type, damage_masking_bases=15, strand_type='doublestrand', is_reverse_complement=False, damage_types='both'):
    """
    Apply damage masking to mask common damage patterns and sequencing artifacts.

    Masks deamination (strand-specific) and/or sequencing artifacts (G→T for single-strand,
    G→T and C→A for double-strand) either across the full sequence or at the sequence ends.

    Args:
        fastaSeq: List of bases in the sequence
        sequence: Alignment positions
        ref_msa: Reference MSA sequence
        damage_masking_type: 'all' for full sequence, 'ends' for sequence ends, or None to skip
        damage_masking_bases: Number of bases from each end to mask (only used when type='ends')
        strand_type: 'singlestrand' or 'doublestrand' to determine deamination and artifact masking
        is_reverse_complement: Whether the read is reverse complement
        damage_types: 'deamination' for C→T/G→A only, 'artifacts' for G→T/C→A only, 'both' for all

    Returns:
        Modified fastaSeq with damage patterns masked
    """
    if not damage_masking_type or damage_masking_type == 'none':
        return fastaSeq

    # Define damage patterns based on strand type and damage_types selection
    damage_patterns = set()

    if strand_type == 'singlestrand':
        # Single-stranded: C→T for forward, G→A for reverse complement
        if is_reverse_complement:
            if damage_types in ('deamination', 'both'):
                damage_patterns.add(('G', 'A'))  # G→A deamination (RC)
            if damage_types in ('artifacts', 'both'):
                damage_patterns.add(('C', 'A'))  # G→T sequencing artifact (complement)
        else:
            if damage_types in ('deamination', 'both'):
                damage_patterns.add(('C', 'T'))  # C→T deamination (forward)
            if damage_types in ('artifacts', 'both'):
                damage_patterns.add(('G', 'T'))  # G→T sequencing artifact
    else:  # doublestrand
        # Double-stranded: both C→T and G→A deamination
        if damage_types in ('deamination', 'both'):
            damage_patterns.add(('C', 'T'))  # C→T deamination
            damage_patterns.add(('G', 'A'))  # G→A deamination
        if damage_types in ('artifacts', 'both'):
            damage_patterns.add(('G', 'T'))  # G→T sequencing artifact
            damage_patterns.add(('C', 'A'))  # C→A sequencing artifact (complement of G→T)

    # Build mapping from gapped position to ungapped position for "ends" detection
    gapped_to_ungapped = {}
    ungapped_pos = 0
    for gapped_pos, base in enumerate(fastaSeq):
        if base != '-':
            gapped_to_ungapped[gapped_pos] = ungapped_pos
            ungapped_pos += 1
    ungapped_read_length = ungapped_pos

    for pos, ref_index in enumerate(sequence):
        if ref_index == -1:  # Skip insertions
            continue

        if ref_index < len(ref_msa):
            ref_base = ref_msa[ref_index].upper()
            read_base = fastaSeq[pos].upper()

            # Skip gap positions
            if read_base == '-':
                continue

            # Check if this position should be considered for masking
            should_mask = False

            if damage_masking_type == 'all':
                # Mask across the full sequence
                should_mask = True
            elif damage_masking_type == 'ends':
                # Mask only at the ends of the sequence (using ungapped positions)
                ungapped_pos = gapped_to_ungapped.get(pos)
                if ungapped_pos is not None:
                    if ungapped_pos < damage_masking_bases or ungapped_pos >= ungapped_read_length - damage_masking_bases:
                        should_mask = True

            # Apply masking if this position should be masked and matches a damage pattern
            if should_mask and (ref_base, read_base) in damage_patterns:
                fastaSeq[pos] = "-"

    return fastaSeq

def generate_fasta_error_profile(fastaSeq, sequence, readNum, error_profile_bases=None, prob_matrix=None, ref_msa=None, is_reverse_complement=False):
    """
    Generate error profile for FASTA files using phylogenetic conditional likelihood calculations.

    Args:
        fastaSeq: DNA sequence
        sequence: Alignment positions
        readNum: Read number for output
        error_profile_bases: If specified, only generate profiles for first and last N bases
        prob_matrix: Probability matrix from error correction (optional)
        ref_msa: Reference MSA sequence (optional)
        is_reverse_complement: Whether read is reverse complement (for output transformation)

    Returns:
        List of error profile entries to write
    """
    if args.no_error_profile:
        return []

    error_profile_entries = []
    error_profile_pos = 0
    last_align_pos = None

    # Find first non-insertion position
    first_pos = None
    for i, align_pos in enumerate(sequence):
        if align_pos != -1:
            first_pos = i
            last_align_pos = align_pos - 1
            break

    if first_pos is None:
        return error_profile_entries

    # Build mapping from gapped position to ungapped position
    gapped_to_ungapped = {}
    ungapped_pos = 0
    for gapped_pos, base in enumerate(fastaSeq):
        if base != '-':
            gapped_to_ungapped[gapped_pos] = ungapped_pos
            ungapped_pos += 1
    ungapped_read_length = ungapped_pos

    # Determine which ungapped positions to include in error profile
    ungapped_positions_to_include = set()
    if error_profile_bases is not None:
        if ungapped_read_length <= 2 * error_profile_bases:
            # If sequence is shorter than 2*error_profile_bases, include all positions
            ungapped_positions_to_include = set(range(ungapped_read_length))
        else:
            # Include first N and last N bases (ungapped)
            ungapped_positions_to_include = set(range(error_profile_bases)) | set(range(ungapped_read_length - error_profile_bases, ungapped_read_length))
    else:
        # Include all positions
        ungapped_positions_to_include = set(range(ungapped_read_length))

    for i, (base, align_pos) in enumerate(zip(fastaSeq, sequence)):
        if align_pos == -1:  # Skip insertions
            continue

        # Skip masked sites (bases changed to "-" due to deamination masking)
        if base == "-":
            continue

        # Get ungapped position for this base
        ungapped_i = gapped_to_ungapped.get(i)
        if ungapped_i is None:
            continue

        # Check if this position should be included in error profile (using ungapped position)
        if ungapped_i not in ungapped_positions_to_include:
            continue

        # Check for deletions (gaps in alignment)
        if last_align_pos is not None and align_pos > last_align_pos + 1:
            # There's a deletion gap, skip these positions in error profile
            error_profile_pos += (align_pos - last_align_pos - 1)

        # Use error rate from command line argument
        error_prob = args.fasta_error_rate
        log_quality_score = math.log(error_prob) if error_prob > 0 else float('-inf')

        # Initialize log likelihoods for all bases
        log_likelihoods = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}

        if prob_matrix is not None and ref_msa is not None and align_pos < len(ref_msa):
            # Use phylogenetic conditional likelihood calculations
            # Formula: log(L(base)) = log(sum_x P(base|x) * L(x)) + log(base_quality_score)
            # Where L(x) = 1 if x is reference base, 0 otherwise

            ref_base = ref_msa[align_pos].upper()
            bases = ['A', 'C', 'G', 'T']
            base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

            # Calculate log likelihoods for each possible base
            for observed_base in bases:
                observed_idx = base_to_idx[observed_base]

                # Calculate sum_x P(base|x) * L(x)
                # L(x) = 1 if x is reference base, 0 otherwise
                if ref_base in base_to_idx:
                    ref_idx = base_to_idx[ref_base]
                    prob_sum = prob_matrix[ref_idx, observed_idx]  # P(observed_base|ref_base) * 1
                else:
                    prob_sum = 0.0

                # Apply the formula: log(L(base)) = log(sum_x P(base|x) * L(x)) + log(base_quality_score)
                if prob_sum > 0:
                    log_likelihoods[observed_base] = math.log(prob_sum) + log_quality_score
                else:
                    log_likelihoods[observed_base] = float('-inf')
        else:
            # Fallback to simple error-rate-based calculation.
            # Use molecule-orientation base: complement for RC reads.
            complement_map_fasta = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            obs_base = complement_map_fasta.get(base.upper(), base.upper()) if is_reverse_complement else base.upper()
            correct_prob = 1.0 - error_prob
            error_share = error_prob / 3.0

            # Convert to log likelihoods (log of probability)
            log_likelihoods[obs_base] = math.log(correct_prob)
            for other_base in ['A', 'C', 'G', 'T']:
                if other_base != obs_base:
                    log_likelihoods[other_base] = math.log(error_share)

        # Add error profile entry (store as dict with position, will format at end)
        error_profile_entries.append({
            'pos': error_profile_pos,
            'A': log_likelihoods['A'],
            'C': log_likelihoods['C'],
            'G': log_likelihoods['G'],
            'T': log_likelihoods['T']
        })

        error_profile_pos += 1
        last_align_pos = align_pos

    # Format entries for output
    # For RC reads: keep entries in alignment order; only swap A↔T and C↔G columns
    # (converts molecule-orientation likelihoods back to alignment orientation for the placer)
    if is_reverse_complement and error_profile_entries:
        return [
            f"{readNum}\t{entry['pos']}\t{entry['T']:.6f}\t{entry['G']:.6f}\t{entry['C']:.6f}\t{entry['A']:.6f}"
            for entry in error_profile_entries
        ]
    else:
        return [
            f"{readNum}\t{entry['pos']}\t{entry['A']:.6f}\t{entry['C']:.6f}\t{entry['G']:.6f}\t{entry['T']:.6f}"
            for entry in error_profile_entries
        ]

def generate_error_profile(fastaSeq, quality, sequence, readNum, encoding, error_profile_bases=None, prob_matrix=None, ref_msa=None, check_fp=None, positional_matrices=None, is_reverse_complement=False, positional_bases=None):
    """
    Generate error profile using phylogenetic conditional likelihood calculations.

    Formula: *P(x | y) = P(calling x | y is true) * (1 - error) if x == y, P(calling x | y is true) + P(calling x | x is true) * error / 3.0
    Where x is the observed base, y is the latent true base, and error is the base quality score

    Args:
        fastaSeq: DNA sequence
        quality: Quality string
        sequence: Alignment positions
        readNum: Read number for output
        encoding: Quality score encoding
        error_profile_bases: If specified, only generate profiles for first and last N bases
        prob_matrix: Probability matrix from error correction (optional)
        ref_msa: Reference MSA sequence (optional)
        is_reverse_complement: Whether read is reverse complement (for positional matrix selection)

    Returns:
        List of error profile entries to write
    """
    if args.no_error_profile:
        return []

    # Debug: Check if phylogenetic parameters are provided
    if readNum < 3:
        phylo_available = prob_matrix is not None and ref_msa is not None
        print(f"DEBUG: generate_error_profile called for read {readNum}, phylogenetic params available: {phylo_available}")

    error_profile_entries = []
    error_profile_pos = 0
    last_align_pos = None

    # Find first non-insertion position
    first_pos = None
    for i, align_pos in enumerate(sequence):
        if align_pos != -1:
            first_pos = i
            last_align_pos = align_pos - 1
            break

    if first_pos is None:
        return error_profile_entries

    # Build mapping from gapped position to ungapped position
    # This is needed for correct positional matrix selection
    gapped_to_ungapped = {}
    ungapped_pos = 0
    for gapped_pos, base in enumerate(fastaSeq):
        if base != '-':
            gapped_to_ungapped[gapped_pos] = ungapped_pos
            ungapped_pos += 1
    ungapped_read_length = ungapped_pos  # Total non-gap bases

    #print(f"Total ungapped: {ungapped_pos}\tTotal: {len(fastaSeq)}")
    #print(f"{gapped_to_ungapped}")

    # Determine which ungapped positions to include in error profile
    ungapped_positions_to_include = set()
    if error_profile_bases is not None:
        if ungapped_read_length <= 2 * error_profile_bases:
            # If sequence is shorter than 2*error_profile_bases, include all positions
            ungapped_positions_to_include = set(range(ungapped_read_length))
        else:
            # Include first N and last N bases (ungapped)
            ungapped_positions_to_include = set(range(error_profile_bases)) | set(range(ungapped_read_length - error_profile_bases, ungapped_read_length))
    else:
        # Include all positions
        ungapped_positions_to_include = set(range(ungapped_read_length))

    for i, (base, qual_char, align_pos) in enumerate(zip(fastaSeq, quality, sequence)):
        #print(f"Base: {base}\tAlign_pos: {align_pos}\ti: {i}\tUngapped_i: {gapped_to_ungapped.get(i)}")
        if align_pos == -1:  # Skip insertions
            continue

        # Skip masked sites (bases changed to "-" due to deamination masking)
        if base == "-":
            continue

        # Get ungapped position for this base
        ungapped_i = gapped_to_ungapped.get(i)
        if ungapped_i is None:
            continue  # This shouldn't happen since we skip gaps above

        # Check if this position should be included in error profile (using ungapped position)
        if ungapped_i not in ungapped_positions_to_include:
            continue

        # Check for deletions (gaps in alignment)
        if last_align_pos is not None and align_pos > last_align_pos + 1:
            # There's a deletion gap, skip these positions in error profile
            error_profile_pos += (align_pos - last_align_pos - 1)

        # Calculate error probability from quality score
        error_prob = quality_to_error_prob(qual_char, encoding)

        # Initialize log likelihoods for all bases
        log_likelihoods = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}

        if not args.quality_score_only and prob_matrix is not None and ref_msa is not None and align_pos < len(ref_msa):
            # Use phylogenetic conditional likelihood calculations
            # New formula: P(x|y) = sum_j P(y->j) * P(x|j)

            bases = ['A', 'C', 'G', 'T']
            base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

            # Select the appropriate probability matrix based on molecule-orientation position.
            # For RC reads, position 0 of the molecule is the 3' end of the alignment, so
            # we flip the position index before looking up the positional matrix.
            complement_map_local = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            if is_reverse_complement:
                observed_base = complement_map_local.get(base.upper(), base.upper())
                pos_for_matrix = ungapped_read_length - ungapped_i - 1
            else:
                observed_base = base.upper()
                pos_for_matrix = ungapped_i

            current_prob_matrix = prob_matrix  # Default to fallback
            _positional_bases = positional_bases if positional_bases is not None else args.positional_matrix_bases
            if positional_matrices is not None:

                if pos_for_matrix < _positional_bases:
                    key = ('5prime', pos_for_matrix)
                    if key in positional_matrices and positional_matrices[key] is not None:
                        current_prob_matrix = positional_matrices[key]
                elif ungapped_read_length - pos_for_matrix <= _positional_bases:
                    pos_from_3prime = ungapped_read_length - pos_for_matrix - 1
                    key = ('3prime', pos_from_3prime)
                    if key in positional_matrices and positional_matrices[key] is not None:
                        current_prob_matrix = positional_matrices[key]
                else:
                    key = 'middle'
                    if key in positional_matrices and positional_matrices[key] is not None:
                        current_prob_matrix = positional_matrices[key]

            # Get observed base index
            if observed_base not in base_to_idx:
                # Skip non-standard bases
                continue
            observed_idx = base_to_idx[observed_base]

            # Get reference base at this position
            ref_base = ref_msa[align_pos].upper()

            # Always use phylogenetic error rate correction
            # Formula: P(X|Y) = sum_j P(Y->j) * P(X|j)
            # where:
            #   - Y is the latent/true base (what we're calculating the likelihood for)
            #   - j is an intermediate base after damage/non-evolutionary change of base
            #   - X is the observed base in the sequencing read
            #   - P(Y->j) is the transition probability from Y to j (from prob_matrix)
            #   - P(X|j) is the sequencing error model:
            #       P(X|j) = (1 - error_prob) when X == j  (correct sequencing)
            #       P(X|j) = error_prob / 3 when X != j    (sequencing error)

            error_profile_probs = {}
            calculation_details = {}  # For debug output

            for latent_base_Y in bases:
                latent_idx = base_to_idx[latent_base_Y]

                # Sum over all possible intermediate bases j
                p_x_given_y = 0.0
                terms = []  # For debug output

                for j_base in bases:
                    j_idx = base_to_idx[j_base]

                    # P(Y->j): transition probability from Y to j
                    p_y_to_j = current_prob_matrix[latent_idx, j_idx]

                    # P(X|j): sequencing error probability
                    if observed_base == j_base:
                        # Sequencing correctly reads j as X (high probability when quality is high)
                        p_x_given_j = 1.0 - error_prob
                    else:
                        # Sequencing error: true base is j but sequencer calls it X
                        p_x_given_j = error_prob / 3.0

                    term_value = p_y_to_j * p_x_given_j
                    p_x_given_y += term_value
                    terms.append((j_base, p_y_to_j, p_x_given_j, term_value))

                error_profile_probs[latent_base_Y] = p_x_given_y
                calculation_details[latent_base_Y] = terms

            # Convert to log likelihoods
            for latent_base_Y in bases:
                if error_profile_probs[latent_base_Y] > 0:
                    log_likelihoods[latent_base_Y] = math.log(error_profile_probs[latent_base_Y])
                else:
                    log_likelihoods[latent_base_Y] = float('-inf')

            # Print examples for the first few positions of first 10 reads to verify math
            if readNum < 10 and error_profile_pos < 3:
                # Determine which matrix was used (based on molecule-orientation position)
                matrix_used = "fallback (middle)"
                if positional_matrices is not None:
                    if pos_for_matrix < _positional_bases:
                        matrix_used = f"5prime position {pos_for_matrix}"
                    elif ungapped_read_length - pos_for_matrix <= _positional_bases:
                        pos_from_3prime = ungapped_read_length - pos_for_matrix - 1
                        matrix_used = f"3prime position {pos_from_3prime}"
                    else:
                        matrix_used = "middle"

                print(f"\n{'='*80}")
                print(f"Example calculation for read {readNum}, position {error_profile_pos}:")
                print(f"  Observed base X: {observed_base}")
                print(f"  Reference base: {ref_base}")
                print(f"  Method: Phylogenetic error rate correction (always)")
                print(f"  Quality character: '{qual_char}' (ASCII: {ord(qual_char)})")
                print(f"  Quality score error prob: {error_prob:.6f}")
                print(f"  Position in read: {i}, align_pos: {align_pos}")
                print(f"  Probability matrix used: {matrix_used}")

                # Show the probability matrix being used
                print(f"\n  Probability matrix P(from->to):")
                print(f"        A        C        G        T")
                for from_base in bases:
                    from_idx = base_to_idx[from_base]
                    print(f"  {from_base}: ", end="")
                    for to_base in bases:
                        to_idx = base_to_idx[to_base]
                        print(f"{current_prob_matrix[from_idx, to_idx]:8.6f} ", end="")
                    print()

                print(f"\n  Phylogenetic calculation details:")
                print(f"  Formula: P(X={observed_base}|Y) = sum_j P(Y->j) * P(X={observed_base}|j)")
                print(f"  Where P(X|j) = {1.0-error_prob:.6f} if j==X, else {error_prob/3.0:.6f}")

                for y_base in bases:
                    print(f"\n  For latent base Y={y_base}:")
                    terms = calculation_details[y_base]
                    for j_base, p_y_to_j, p_x_given_j, term_value in terms:
                        marker = " <-- X==j" if j_base == observed_base else ""
                        print(f"    j={j_base}: P({y_base}->{j_base})={p_y_to_j:.6f} × P(X={observed_base}|{j_base})={p_x_given_j:.6f} = {term_value:.6f}{marker}")
                    print(f"    Sum: P(X={observed_base}|Y={y_base}) = {error_profile_probs[y_base]:.6f}")
                    print(f"    Log likelihood: {log_likelihoods[y_base]:.6f}")

                print(f"\n  Final log likelihoods:")
                for y_base in bases:
                    print(f"    Y={y_base}: {log_likelihoods[y_base]:.6f}")
                print(f"{'='*80}")
        else:
            # Fallback to simple quality-based calculation
            correct_prob = 1.0 - error_prob
            error_share = error_prob / 3.0  # Divide error probability among 3 incorrect bases

            # Convert to log likelihoods (log of probability)
            log_likelihoods[base] = float('-inf') if correct_prob == 0 else math.log(correct_prob)
            for other_base in ['A', 'C', 'G', 'T']:
                if other_base != base:
                    log_likelihoods[other_base] = float('-inf') if error_share == 0 else math.log(error_share)

        # Add error profile entry (store as dict with position, will format at end)
        error_profile_entries.append({
            'pos': error_profile_pos,
            'A': log_likelihoods['A'],
            'C': log_likelihoods['C'],
            'G': log_likelihoods['G'],
            'T': log_likelihoods['T']
        })

        error_profile_pos += 1
        last_align_pos = align_pos

    # Format entries for output
    # For RC reads: keep entries in alignment order; only swap A↔T and C↔G columns
    # (converts molecule-orientation likelihoods back to alignment orientation for the placer)
    if is_reverse_complement and error_profile_entries:
        return [
            f"{readNum}\t{entry['pos']}\t{entry['T']:.6f}\t{entry['G']:.6f}\t{entry['C']:.6f}\t{entry['A']:.6f}"
            for entry in error_profile_entries
        ]
    else:
        return [
            f"{readNum}\t{entry['pos']}\t{entry['A']:.6f}\t{entry['C']:.6f}\t{entry['G']:.6f}\t{entry['T']:.6f}"
            for entry in error_profile_entries
        ]

def reprocess_error_profiles_with_phylogenetic_calculations(valid_reads_data, prob_matrix, ref_msa, positional_prob_matrices=None):
    """
    Reprocess error profiles using phylogenetic conditional likelihood calculations.

    Args:
        valid_reads_data: List of valid reads data
        prob_matrix: Probability matrix from substitution analysis (fallback)
        ref_msa: Reference MSA sequence
        positional_prob_matrices: Dictionary of position-specific probability matrices (2N + 1 matrices)
    """
    import numpy as np

    # Convert probability matrix dictionary to numpy array
    bases = ['A', 'C', 'G', 'T']
    prob_matrix_array = np.zeros((4, 4))

    for i, from_base in enumerate(bases):
        for j, to_base in enumerate(bases):
            prob_key = f"{from_base}→{to_base}"
            prob_matrix_array[i, j] = prob_matrix.get(prob_key, 0.0)

    print(f"Converted fallback probability matrix to {prob_matrix_array.shape} numpy array")

    # Convert positional probability matrices to numpy arrays
    positional_matrices_arrays = None
    if positional_prob_matrices is not None:
        positional_matrices_arrays = {}
        for key, matrix_dict in positional_prob_matrices.items():
            if matrix_dict is not None:
                matrix_array = np.zeros((4, 4))
                for i, from_base in enumerate(bases):
                    for j, to_base in enumerate(bases):
                        prob_key = f"{from_base}→{to_base}"
                        matrix_array[i, j] = matrix_dict.get(prob_key, 0.0)
                positional_matrices_arrays[key] = matrix_array
            else:
                positional_matrices_arrays[key] = None  # No data, will use fallback

        print(f"Converted {len([k for k, v in positional_matrices_arrays.items() if v is not None])} positional probability matrices")

    # Reprocess error profiles for each read - use the actual error profile output file
    error_profile_output = args.outError

    # Create check fails output file
    import os
    file_prefix = os.path.splitext(args.outError)[0]
    check_fails_file = f"{file_prefix}_check_fails.txt"

    with open(error_profile_output, 'w') as ep, open(check_fails_file, 'w') as check_fp:
        # print(f"DEBUG: Reprocessing {len(valid_reads_data)} reads for phylogenetic error profiles")

        # The validReadsData should be in the same order as processed originally
        for i, (read_name, fastaSeq, sequence, read_ref_msa, ref_id, is_reverse_complement, quality) in enumerate(valid_reads_data):
            original_readNum = i  # Use the index as the read number to match original processing 

            #print(f"{fastaSeq}")

            # if i < 3:  # Debug first few reads
            #     print(f"DEBUG: Processing read {i} ({read_name}), sequence length: {len(sequence)}, fastaSeq length: {len(fastaSeq)}")

            # generate_error_profile() now handles molecule orientation internally,
            # so we always pass the alignment-orientation sequence as-is.
            orientation_fastaSeq = list(fastaSeq)

            if len(quality) != len(orientation_fastaSeq):
                print(f"Warning: Length of quality and sequence are not the same length")

            ## Apply masking
            #strand_type = 'singlestrand' if args.singlestrand else 'doublestrand'

            ## Apply transition masking
            #masked_fastaSeq = apply_transition_masking(masked_fastaSeq, sequence, read_ref_msa)

            ## Apply damage masking (deamination and/or artifacts)
            #masked_fastaSeq = apply_damage_masking(masked_fastaSeq, sequence, read_ref_msa, args.damage_masking, args.damage_masking_bases, strand_type, is_reverse_complement, args.damage_types)

            # Generate error profile with phylogenetic calculations using masked sequence
            if quality:
                error_entries = generate_error_profile(
                    orientation_fastaSeq, quality, sequence, original_readNum, args.encoding,
                    args.error_profile_bases, prob_matrix_array, read_ref_msa, check_fp, positional_matrices_arrays,
                    is_reverse_complement
                )

                for entry in error_entries:
                    ep.write(entry + "\n")

    print(f"Reprocessed error profiles written to: {error_profile_output}")

def write_error_profiles_with_damage_matrix(validReadsData, positional_matrices_arrays, fallback_array, n_positions):
    """
    Write error profiles using a pre-loaded damage probability matrix.

    This function mirrors reprocess_error_profiles_with_phylogenetic_calculations but
    accepts numpy arrays directly (no dict→numpy conversion step) and uses n_positions
    (inferred from the .npy matrix shape) as the positional boundary instead of
    args.positional_matrix_bases.

    Args:
        validReadsData: List of (read_name, fastaSeq, sequence, ref_msa, ref_id,
                        is_reverse_complement, quality) tuples
        positional_matrices_arrays: Dict mapping ('5prime', i), 'middle', ('3prime', j)
                                    to (4,4) numpy arrays from the damage .npy
        fallback_array: (4,4) numpy array used when no positional matrix matches
                        (typically the middle slice of the damage matrix)
        n_positions: Number of end positions N (auto-detected from npy shape)
    """
    import os
    error_profile_output = args.outError
    file_prefix = os.path.splitext(args.outError)[0]
    check_fails_file = f"{file_prefix}_check_fails.txt"

    print(f"\nWriting error profiles using damage matrix (N={n_positions} end positions)...")

    with open(error_profile_output, 'w') as ep, open(check_fails_file, 'w') as check_fp:
        for i, (read_name, fastaSeq, sequence, read_ref_msa, ref_id,
                is_reverse_complement, quality) in enumerate(validReadsData):
            orientation_fastaSeq = list(fastaSeq)

            if len(quality) != len(orientation_fastaSeq):
                print(f"Warning: Length of quality and sequence differ for read {i}")

            if quality:
                error_entries = generate_error_profile(
                    orientation_fastaSeq, quality, sequence, i, args.encoding,
                    args.error_profile_bases, fallback_array, read_ref_msa,
                    check_fp, positional_matrices_arrays, is_reverse_complement,
                    positional_bases=n_positions
                )
                for entry in error_entries:
                    ep.write(entry + "\n")

    print(f"Error profiles written to: {error_profile_output}")

def reverse_complement(sequence):
    """
    Generate reverse complement of a DNA sequence.

    Args:
        sequence: DNA sequence string

    Returns:
        Reverse complement sequence
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def reverse_quality_scores(quality):
    """
    Reverse quality scores for reverse complement reads.

    Args:
        quality: Quality score string

    Returns:
        Reversed quality score string
    """
    return quality[::-1]

def count_mismatches_after_masking(fastaSeq, sequence, ref_msa):
    """
    Count mismatches between masked read sequence and reference.

    Args:
        fastaSeq: DNA sequence (after masking)
        sequence: Alignment positions
        ref_msa: Reference MSA sequence

    Returns:
        Number of mismatches (excluding masked positions marked as "-")
    """
    mismatches = 0

    for pos, ref_index in enumerate(sequence):
        if ref_index == -1:  # Skip insertions
            continue

        if ref_index < len(ref_msa):
            ref_base = ref_msa[ref_index].upper()
            read_base = fastaSeq[pos].upper()

            # Skip masked positions (marked as "-")
            if read_base == "-":
                continue

            # Count mismatch if bases are different
            if ref_base != read_base and ref_base != 'N' and read_base != 'N':
                mismatches += 1

    return mismatches

def write_mismatched_read_to_fasta(read_name, original_sequence, mismatch_count, fasta_file_handle):
    """
    Write a read with mismatches to a FASTA file.

    Args:
        read_name: Name of the read
        original_sequence: Original DNA sequence (before masking)
        mismatch_count: Number of mismatches found after masking
        fasta_file_handle: Open file handle for writing FASTA output
    """
    fasta_file_handle.write(f">{read_name} mismatches={mismatch_count}\n")
    fasta_file_handle.write(f"{original_sequence}\n")

def perform_substitution_analysis(valid_reads_data, output_prefix="substitution_analysis"):
    """
    Perform substitution analysis similar to analyze_substitutions.py

    Args:
        valid_reads_data: List of tuples (read, fastaSeq, sequence, ref_msa, ref_id, is_reverse_complement, quality)
        output_prefix: Prefix for output files

    Note on orientation
    -------------------
    Substitution counts and per-position damage profiles produced here are in
    **molecule orientation**.  For reverse-complement (RC) reads, both the read
    base and reference base are complemented, and positions are flipped so that
    position 0 is the 5' end of the original molecule for all reads.

    This means C→T deamination always appears at the 5' end (position 0) for
    both forward and RC reads, and the positional matrices are built consistently
    across both strand orientations.
    """
    from collections import defaultdict, Counter
    import math

    if args.debug:
        print("\nPerforming substitution analysis...")

    # Initialize data structures for all, forward, and reverse complement reads
    data_sets = {
        'all': {'substitutions': Counter(), 'base_counts': Counter(), 'ref_base_counts': Counter(), 'total_substitutions': 0, 'total_comparisons': 0},
        'forward': {'substitutions': Counter(), 'base_counts': Counter(), 'ref_base_counts': Counter(), 'total_substitutions': 0, 'total_comparisons': 0},
        'reverse': {'substitutions': Counter(), 'base_counts': Counter(), 'ref_base_counts': Counter(), 'total_substitutions': 0, 'total_comparisons': 0}
    }

    # Position-specific data structures
    position_data = {
        'all': {'position_subs': defaultdict(lambda: defaultdict(int)), 'reverse_position_subs': defaultdict(lambda: defaultdict(int))},
        'forward': {'position_subs': defaultdict(lambda: defaultdict(int)), 'reverse_position_subs': defaultdict(lambda: defaultdict(int))},
        'reverse': {'position_subs': defaultdict(lambda: defaultdict(int)), 'reverse_position_subs': defaultdict(lambda: defaultdict(int))}
    }

    # Position-specific data for rate matrix calculation (2N + 1 matrices)
    # N positions from 5' end, N positions from 3' end, 1 for middle
    positional_matrix_bases = args.positional_matrix_bases
    positional_substitutions = {
        '5prime': {},  # Keys: position 0 to N-1 from 5' end
        '3prime': {},  # Keys: position 0 to N-1 from 3' end
        'middle': {'substitutions': Counter(), 'ref_base_counts': Counter()}  # Middle positions
    }

    # Initialize dictionaries for each position
    for pos in range(positional_matrix_bases):
        positional_substitutions['5prime'][pos] = {'substitutions': Counter(), 'ref_base_counts': Counter()}
        positional_substitutions['3prime'][pos] = {'substitutions': Counter(), 'ref_base_counts': Counter()}

    # Count reads by orientation
    forward_count = 0
    reverse_count = 0

    # Get reference MSA from first read for error profile calculations
    first_ref_msa = None
    if valid_reads_data:
        first_ref_msa = valid_reads_data[0][3]  # ref_msa is index 3 in the tuple

    for read, fastaSeq, sequence, ref_msa, ref_id, is_reverse_complement, quality in valid_reads_data:
        #if forward_count > 10:
        #    continue
        #if reverse_count > 10:
        #    continue

        # Determine read type
        read_type = 'reverse' if is_reverse_complement else 'forward'
        if is_reverse_complement:
            reverse_count += 1
        else:
            forward_count += 1

        # Find substitutions for this read (returns ungapped positions and length)
        substitutions, position_subs, reverse_position_subs, ref_base_counts, ungapped_read_length = find_read_substitutions(
            fastaSeq, sequence, ref_msa, is_reverse_complement
        )

        #print(f"{substitutions}")
        #print(f"{position_subs}")
        #print(f"{reverse_position_subs}")

        # Update data for 'all' and specific orientation
        for dataset in ['all', read_type]:
            # Update reference base counts
            for base, count in ref_base_counts.items():
                data_sets[dataset]['ref_base_counts'][base] += count

            # Update substitution counts
            for ref_pos, ref_base, read_base, read_pos, reverse_pos in substitutions:
                sub_type = f"{ref_base}→{read_base}"
                data_sets[dataset]['substitutions'][sub_type] += 1
                data_sets[dataset]['total_substitutions'] += 1

            # Update position-specific counts
            for pos, sub_dict in position_subs.items():
                for sub_type, count in sub_dict.items():
                    position_data[dataset]['position_subs'][pos][sub_type] += count

            # Update reverse position-specific counts
            for pos, sub_dict in reverse_position_subs.items():
                for sub_type, count in sub_dict.items():
                    position_data[dataset]['reverse_position_subs'][pos][sub_type] += count

        # Update positional rate matrix data (for 2N + 1 matrices)
        # For RC reads, flip positions so 5prime/3prime refer to original molecule orientation
        # Use ungapped read length for correct positional calculations
        read_length = ungapped_read_length
        for ref_pos, ref_base, read_base, read_pos, reverse_pos in substitutions:
            sub_type = f"{ref_base}→{read_base}"

            if read_pos < positional_matrix_bases:
                positional_substitutions['5prime'][read_pos]['substitutions'][sub_type] += 1
            elif read_length - read_pos <= positional_matrix_bases:
                pos_from_3prime = read_length - read_pos - 1
                positional_substitutions['3prime'][pos_from_3prime]['substitutions'][sub_type] += 1
            else:
                positional_substitutions['middle']['substitutions'][sub_type] += 1

        # Update reference base counts for each positional zone.
        # Use molecule orientation: for RC reads, complement ref_base and flip position.
        complement_map_sub = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        gapped_to_ungapped_ref = {}
        ungapped_pos = 0
        for gapped_pos, base in enumerate(fastaSeq):
            if base != '-':
                gapped_to_ungapped_ref[gapped_pos] = ungapped_pos
                ungapped_pos += 1

        for i, base in enumerate(fastaSeq):
            if base == '-':  # Skip gap positions
                continue
            if i < len(sequence) and sequence[i] != -1:
                ref_pos = sequence[i]
                if ref_pos < len(ref_msa):
                    ref_base = ref_msa[ref_pos].upper()
                    if ref_base in ['A', 'C', 'G', 'T']:
                        ungapped_i = gapped_to_ungapped_ref.get(i, i)
                        if is_reverse_complement:
                            count_pos_ref = ungapped_read_length - ungapped_i - 1
                            count_ref_base = complement_map_sub.get(ref_base, ref_base)
                        else:
                            count_pos_ref = ungapped_i
                            count_ref_base = ref_base
                        if count_pos_ref < positional_matrix_bases:
                            positional_substitutions['5prime'][count_pos_ref]['ref_base_counts'][count_ref_base] += 1
                        elif read_length - count_pos_ref <= positional_matrix_bases:
                            pos_from_3prime = read_length - count_pos_ref - 1
                            positional_substitutions['3prime'][pos_from_3prime]['ref_base_counts'][count_ref_base] += 1
                        else:
                            positional_substitutions['middle']['ref_base_counts'][count_ref_base] += 1

        # Count valid base comparisons and observed base counts for all datasets.
        # Use molecule orientation: for RC reads, complement the observed base.
        for i, base in enumerate(fastaSeq):
            if i < len(sequence) and sequence[i] != -1:
                ref_pos = sequence[i]
                if ref_pos < len(ref_msa):
                    read_base = base.upper()
                    ref_base = ref_msa[ref_pos].upper()
                    if read_base in 'ATCG' and ref_base in 'ATCG' and read_base != "-":
                        data_sets['all']['total_comparisons'] += 1
                        data_sets[read_type]['total_comparisons'] += 1
                        # Count observed read bases in molecule orientation
                        count_obs_base = complement_map_sub.get(read_base, read_base) if is_reverse_complement else read_base
                        data_sets['all']['base_counts'][count_obs_base] += 1
                        data_sets[read_type]['base_counts'][count_obs_base] += 1

    # Write results for all datasets
    print(f"\nRead orientation breakdown:")
    print(f"  Forward reads: {forward_count}")
    print(f"  Reverse complement reads: {reverse_count}")
    print(f"  Total reads: {forward_count + reverse_count}")

    # Print detailed read examples for analysis (if requested)
    if args.show_read_examples:
        print_read_examples(valid_reads_data, num_examples=3)

    # Write comprehensive results and get probability matrix
    prob_matrix, _ = write_orientation_substitution_results(data_sets, position_data, output_prefix)

    # Calculate positional probability matrices (2N + 1 matrices)
    positional_prob_matrices = None
    if args.gtrParams and positional_matrix_bases > 0:
        print("\n" + "=" * 60)
        print("CALCULATING POSITIONAL RATE MATRICES")
        print("=" * 60)

        gtr_params = parse_gtr_parameters(args.gtrParams)
        positional_prob_matrices = {}
        bases = ['A', 'C', 'G', 'T']

        # Calculate matrices for 5' end positions
        print(f"\nCalculating rate matrices for first {positional_matrix_bases} positions (5' end):")
        for pos in range(positional_matrix_bases):
            data = positional_substitutions['5prime'][pos]
            if sum(data['substitutions'].values()) > 0:  # Only if we have data
                error_matrix = calculate_error_rate_matrix_from_data(data, gtr_params, args.error_fudge_factor)
                # rate_matrix = create_error_rate_matrix(error_matrix)
                prob_matrix_pos = create_error_prob_matrix(error_matrix)

                if prob_matrix_pos is None:
                    print(f"\n  Error: Could not calculate probability matrix for 5' Position {pos}")
                    print(f"  This usually means numpy/scipy are not available.")
                    positional_prob_matrices[('5prime', pos)] = None
                    continue

                positional_prob_matrices[('5prime', pos)] = prob_matrix_pos
                print(f"\n  5' Position {pos}: {sum(data['substitutions'].values())} substitutions")

                # Print detailed substitution statistics
                print(f"    " + "-" * 108)
                print(f"      Type    Count  % of Subs    Ref Norm     Obs/Exp Ratio")
                print(f"    " + "-" * 108)
                total_subs = sum(data['substitutions'].values())
                for sub_type in sorted(error_matrix['obs_exp_ratios'].keys(), key=lambda x: error_matrix['obs_exp_ratios'][x], reverse=True):
                    count = data['substitutions'][sub_type]
                    pct = (count / total_subs) * 100
                    ref_base = sub_type.split('→')[0]
                    ref_norm = count / data['ref_base_counts'][ref_base] if data['ref_base_counts'][ref_base] > 0 else 0
                    obs_exp_ratio = error_matrix['obs_exp_ratios'][sub_type]
                    print(f"       {sub_type}: {count:8d} ({pct:5.1f}%)  {ref_norm:8.6f}     {obs_exp_ratio:8.6f}")
                print()

                # Print rate matrix
                # print(f"    Rate Matrix:")
                # for from_base in bases:
                #     print(f"      {from_base}: ", end="")
                #     for to_base in bases:
                #         rate_key = f"{from_base}→{to_base}"
                #         rate_value = rate_matrix.get(rate_key, 0.0)
                #         print(f"{to_base}={rate_value:>8.6f} ", end="")
                #     print()

                # Print probability matrix
                print(f"    Probability Matrix:")
                for from_base in bases:
                    print(f"      {from_base}: ", end="")
                    for to_base in bases:
                        prob_key = f"{from_base}→{to_base}"
                        prob_value = prob_matrix_pos.get(prob_key, 0.0)
                        print(f"{to_base}={prob_value:>8.6f} ", end="")
                    print()
            else:
                print(f"  Position {pos}: No data, will use fallback")
                positional_prob_matrices[('5prime', pos)] = None

        # Calculate matrix for middle positions
        print(f"\nCalculating rate matrix for middle positions:")
        middle_data = positional_substitutions['middle']
        if sum(middle_data['substitutions'].values()) > 0:
            error_matrix = calculate_error_rate_matrix_from_data(middle_data, gtr_params, args.error_fudge_factor)
            # rate_matrix = create_error_rate_matrix(error_matrix)
            prob_matrix_middle = create_error_prob_matrix(error_matrix)

            if prob_matrix_middle is None:
                print(f"\n  Error: Could not calculate probability matrix for middle positions")
                print(f"  This usually means numpy/scipy are not available.")
                positional_prob_matrices['middle'] = None
            else:
                positional_prob_matrices['middle'] = prob_matrix_middle
                print(f"\n  Middle positions: {sum(middle_data['substitutions'].values())} substitutions")

                # Print detailed substitution statistics
                print(f"    " + "-" * 108)
                print(f"      Type    Count  % of Subs   Ref Norm   Ref/1000     Obs/Exp Ratio")
                print(f"    " + "-" * 108)
                total_subs = sum(middle_data['substitutions'].values())
                for sub_type in sorted(error_matrix['obs_exp_ratios'].keys(), key=lambda x: error_matrix['obs_exp_ratios'][x], reverse=True):
                    count = middle_data['substitutions'][sub_type]
                    pct = (count / total_subs) * 100
                    ref_base = sub_type.split('→')[0]
                    ref_norm = count / middle_data['ref_base_counts'][ref_base] if middle_data['ref_base_counts'][ref_base] > 0 else 0
                    obs_exp_ratio = error_matrix['obs_exp_ratios'][sub_type]
                    print(f"       {sub_type}: {count:8d} ({pct:5.1f}%)  {ref_norm:8.6f}  {ref_norm*1000:8.3f}     {obs_exp_ratio:8.6f}")
                print()

                # Print rate matrix
                # print(f"    Rate Matrix:")
                # for from_base in bases:
                #     print(f"      {from_base}: ", end="")
                #     for to_base in bases:
                #         rate_key = f"{from_base}→{to_base}"
                #         rate_value = rate_matrix.get(rate_key, 0.0)
                #         print(f"{to_base}={rate_value:>8.6f} ", end="")
                #     print()

                # Print probability matrix
                print(f"    Probability Matrix:")
                for from_base in bases:
                    print(f"      {from_base}: ", end="")
                    for to_base in bases:
                        prob_key = f"{from_base}→{to_base}"
                        prob_value = prob_matrix_middle.get(prob_key, 0.0)
                        print(f"{to_base}={prob_value:>8.6f} ", end="")
                    print()
        else:
            print(f"  Middle: No data, will use fallback")
            positional_prob_matrices['middle'] = None

        # Calculate matrices for 3' end positions
        print(f"\nCalculating rate matrices for last {positional_matrix_bases} positions (3' end):")
        for pos in range(positional_matrix_bases):
            data = positional_substitutions['3prime'][pos]
            if sum(data['substitutions'].values()) > 0:
                error_matrix = calculate_error_rate_matrix_from_data(data, gtr_params, args.error_fudge_factor)
                # rate_matrix = create_error_rate_matrix(error_matrix)
                prob_matrix_pos = create_error_prob_matrix(error_matrix)

                if prob_matrix_pos is None:
                    print(f"\n  Error: Could not calculate probability matrix for 3' Position {pos}")
                    print(f"  This usually means numpy/scipy are not available.")
                    positional_prob_matrices[('3prime', pos)] = None
                    continue

                positional_prob_matrices[('3prime', pos)] = prob_matrix_pos
                print(f"\n  3' Position {pos}: {sum(data['substitutions'].values())} substitutions")

                # Print detailed substitution statistics
                print(f"    " + "-" * 108)
                print(f"      Type    Count  % of Subs   Ref Norm   Ref/1000     Obs/Exp Ratio")
                print(f"    " + "-" * 108)
                total_subs = sum(data['substitutions'].values())
                for sub_type in sorted(error_matrix['obs_exp_ratios'].keys(), key=lambda x: error_matrix['obs_exp_ratios'][x], reverse=True):
                    count = data['substitutions'][sub_type]
                    pct = (count / total_subs) * 100
                    ref_base = sub_type.split('→')[0]
                    ref_norm = count / data['ref_base_counts'][ref_base] if data['ref_base_counts'][ref_base] > 0 else 0
                    obs_exp_ratio = error_matrix['obs_exp_ratios'][sub_type]
                    print(f"       {sub_type}: {count:8d} ({pct:5.1f}%)  {ref_norm:8.6f}  {ref_norm*1000:8.3f}     {obs_exp_ratio:8.6f}")
                print()

                # Print rate matrix
                # print(f"    Rate Matrix:")
                # for from_base in bases:
                #     print(f"      {from_base}: ", end="")
                #     for to_base in bases:
                #         rate_key = f"{from_base}→{to_base}"
                #         rate_value = rate_matrix.get(rate_key, 0.0)
                #         print(f"{to_base}={rate_value:>8.6f} ", end="")
                #     print()

                # Print probability matrix
                print(f"    Probability Matrix:")
                for from_base in bases:
                    print(f"      {from_base}: ", end="")
                    for to_base in bases:
                        prob_key = f"{from_base}→{to_base}"
                        prob_value = prob_matrix_pos.get(prob_key, 0.0)
                        print(f"{to_base}={prob_value:>8.6f} ", end="")
                    print()
            else:
                print(f"  Position {pos} from 3' end: No data, will use fallback")
                positional_prob_matrices[('3prime', pos)] = None

        print("\n" + "=" * 60)

    # Return probability matrices and reference MSA for error profile calculations
    return prob_matrix, first_ref_msa, positional_prob_matrices

def find_read_substitutions(fastaSeq, sequence, ref_msa, is_reverse_complement):
    """
    Find substitutions between read and reference sequences.

    IMPORTANT: Positions are tracked in ungapped coordinates (counting only actual bases,
    not gap characters) to ensure proper damage pattern analysis.

    Bases and positions are reported in molecule orientation for RC reads:
    - read_base and ref_base are complemented for RC reads
    - position 0 is the 5' end of the molecule for all reads
    - positions and bases are in molecule orientation for RC reads.

    Args:
        fastaSeq: List of bases in the read sequence (alignment orientation, may contain gaps)
        sequence: Alignment positions
        ref_msa: Reference MSA sequence
        is_reverse_complement: Whether the read was reverse complemented in alignment

    Returns:
        substitutions, position_subs, reverse_position_subs, ref_base_counts, ungapped_read_length
    """
    from collections import defaultdict

    substitutions = []
    position_subs = defaultdict(lambda: defaultdict(int))
    reverse_position_subs = defaultdict(lambda: defaultdict(int))
    ref_base_counts = defaultdict(int)
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-'}

    # Calculate ungapped read length (count actual bases, not gaps)
    ungapped_read_length = sum(1 for base in fastaSeq if base != '-')

    # All positions and bases are kept in alignment orientation for both forward and
    # RC reads.  calc_subs handles original-molecule orientation; gRPI only needs the
    # alignment-relative view for downstream error profile / placement.
    gapped_to_ungapped = {}
    ungapped_pos = 0
    for gapped_pos, base in enumerate(fastaSeq):
        if base != '-':
            gapped_to_ungapped[gapped_pos] = ungapped_pos
            ungapped_pos += 1

    gapped_to_ungapped_for_loop = gapped_to_ungapped

    ref_sequence = [ref_msa[ref_pos].upper() for ref_pos in sequence]
    newFasta = fastaSeq.copy()

    for gapped_read_pos, ref_pos in enumerate(sequence):
        if ref_pos == -1:  # Skip insertions
            continue

        if gapped_read_pos >= len(fastaSeq) or ref_pos >= len(ref_msa):
            continue

        # For reverse complement reads, the alignment positions are for the RC'd read
        # but we want to use the original read, so we need to map positions correctly
        read_base = newFasta[gapped_read_pos].upper()
        ref_base = ref_sequence[gapped_read_pos].upper()
        # if is_reverse_complement:
        #     print(f"rc: {read_base} vs {ref_base}")
        # else:
        #     print(f"oo: {read_base} vs {ref_base}")


        # Skip masked positions and ambiguous characters
        if read_base == "-":
            continue

        if read_base in 'ATCG' and ref_base in 'ATCG':
            # Convert to molecule orientation for RC reads
            if is_reverse_complement:
                count_read_base = complement_map.get(read_base, read_base)
                count_ref_base  = complement_map.get(ref_base, ref_base)
            else:
                count_read_base = read_base
                count_ref_base  = ref_base

            # Count reference bases (what was available to substitute from)
            ref_base_counts[count_ref_base] += 1

            if count_read_base != count_ref_base:
                # Get ungapped position for this base
                if gapped_read_pos in gapped_to_ungapped_for_loop:
                    ungapped_read_pos = gapped_to_ungapped_for_loop[gapped_read_pos]
                else:
                    continue  # Gap position, skip

                # Calculate reverse position (from end of read) using ungapped position
                reverse_pos = ungapped_read_length - ungapped_read_pos - 1

                # For RC reads: reverse_pos is distance from 5' of molecule (molecule orientation)
                if is_reverse_complement:
                    count_pos     = reverse_pos       # distance from 5' of molecule
                    count_rev_pos = ungapped_read_pos
                else:
                    count_pos     = ungapped_read_pos
                    count_rev_pos = reverse_pos

                # Store substitution with molecule-orientation positions and bases
                substitutions.append((ref_pos, count_ref_base, count_read_base, count_pos, count_rev_pos))
                position_subs[count_pos][f"{count_ref_base}→{count_read_base}"] += 1
                reverse_position_subs[count_rev_pos][f"{count_ref_base}→{count_read_base}"] += 1

    return substitutions, position_subs, reverse_position_subs, ref_base_counts, ungapped_read_length

def write_orientation_substitution_results(data_sets, position_data, output_prefix):
    """
    Write substitution analysis results broken down by orientation.
    """
    # Initialize probability matrix for return
    prob_matrix = None
    ref_msa = None

    # Print console results only for 'all' dataset
    dataset = data_sets['all']

    print("=" * 60)
    print("SUBSTITUTION ANALYSIS RESULTS - ALL READS")
    print("=" * 60)

    # Print GTR substitution matrix if parameters are available
    if gtr_params:
        print_gtr_substitution_matrix(gtr_params)

    print(f"Total substitutions found: {dataset['total_substitutions']}")
    print(f"Total base comparisons: {dataset['total_comparisons']}")
    if dataset['total_comparisons'] > 0:
        print(f"Overall substitution rate: {dataset['total_substitutions']/dataset['total_comparisons']*100:.3f}%")
    print()

    print("Observed read base counts:")
    print("-" * 30)
    for base in ['A', 'T', 'C', 'G']:
        count = dataset['base_counts'].get(base, 0)
        percentage = count / dataset['total_comparisons'] * 100 if dataset['total_comparisons'] > 0 else 0
        print(f"  {base}: {count:>8} ({percentage:5.1f}%)")
    print()

    print("Reference base counts:")
    print("-" * 30)
    for base in ['A', 'T', 'C', 'G']:
        count = dataset['ref_base_counts'].get(base, 0)
        percentage = count / dataset['total_comparisons'] * 100 if dataset['total_comparisons'] > 0 else 0
        print(f"  {base}: {count:>8} ({percentage:5.1f}%)")
    print()

    print("Substitution type counts and normalized rates:")
    print("(Ref Norm = substitutions / reference source base count)")
    if gtr_params:
        print("(Obs/Exp Ratio = observed rate / expected rate from GTR with data frequencies)")
        print("-" * 113)
        print(f"{'Type':>6} {'Count':>8} {'% of Subs':>10} {'Ref Norm':>10} {'Ref Count':>10} {'Expected Rate':>10} {'Expected Count':>10} {'Obs/Exp Ratio':>15} {'Obs/Exp Count Ratio':>15}")
        print("-" * 113)
    else:
        print("-" * 90)
        print(f"{'Type':>6} {'Count':>8} {'% of Subs':>10} {'Ref Norm':>10} {'Ref Count':>10}")
        print("-" * 90)

    for sub_type, count in dataset['substitutions'].most_common():
        percentage = count / dataset['total_substitutions'] * 100 if dataset['total_substitutions'] > 0 else 0

        # Extract bases for normalization
        ref_base, observed_base = sub_type.split('→')

        # # Observed base normalization (target)
        # observed_count = dataset['base_counts'].get(observed_base, 0)
        # obs_normalized_rate = count / observed_count if observed_count > 0 else 0
        # obs_per_thousand = obs_normalized_rate * 1000

        # Reference base normalization (source)
        ref_count = dataset['ref_base_counts'].get(ref_base, 0)
        ref_normalized_rate = count / ref_count if ref_count > 0 else 0
        ref_per_thousand = ref_normalized_rate * 1000

        # Obs/Exp ratio calculation if parameters available
        if gtr_params:
            expected_rate = get_expected_gtr_rate_from_data(sub_type, gtr_params, dataset['ref_base_counts'])
            obs_exp_ratio = ref_normalized_rate / expected_rate if expected_rate and expected_rate > 0 else float('inf')
            expected_count = expected_rate * ref_count
            obs_exp_count = count / expected_count
            print(f"{sub_type:>6}: {count:>8} ({percentage:5.1f}%) {ref_normalized_rate:9.6f} {ref_count:>8} {expected_rate:14.6} {expected_count:14.8} {obs_exp_ratio:14.6f} {obs_exp_count:14.6f}")
        else:
            print(f"{sub_type:>6}: {count:>8} ({percentage:5.1f}%) {ref_normalized_rate:9.6f} {ref_count:>8}")
    print()

    # Calculate and display error rate matrix if GTR parameters available
    if gtr_params:
        error_matrix = calculate_error_rate_matrix(dataset, gtr_params, args.error_fudge_factor)

        if error_matrix:
            print("=" * 80)
            print("ERROR RATE ANALYSIS")
            print("=" * 80)
            print(f"Baseline substitution (smallest Obs/Exp ratio): {error_matrix['smallest_substitution']}")
            print(f"Smallest Obs/Exp ratio: {error_matrix['smallest_ratio']:.6f}")
            print(f"Baseline ref_norm: {error_matrix['ref_norm_small']:.6f}")
            print(f"Baseline expected_rate: {error_matrix['expected_rate_small']:.6f}")
            print(f"Correction factor: {error_matrix['correction_factor']:.6f}")
            if error_matrix['fudge_factor'] != 0.0:
                print(f"Fudge factor: {error_matrix['fudge_factor']:.6f}")
            print()

            print("Error Rate Matrix:")
            formula = "Formula: error_rate = ref_norm - expected_rate × correction_factor"
            if error_matrix['fudge_factor'] != 0.0:
                formula += f" + {error_matrix['fudge_factor']:.6f}"
            print(formula)
            print("-" * 85)
            print(f"{'Type':>6} {'Ref Norm':>12} {'Expected Rate':>14} {'Error Rate':>12} {'Obs/Exp Ratio':>15}")
            print("-" * 85)

            # Sort by error rate (highest first) to highlight problematic substitutions
            sorted_subs = sorted(error_matrix['error_rates'].items(), key=lambda x: x[1], reverse=True)

            for sub_type, error_rate in sorted_subs:
                ref_norm = error_matrix['ref_norms'][sub_type]
                expected_rate = error_matrix['expected_rates'][sub_type]
                obs_exp_ratio = error_matrix['obs_exp_ratios'][sub_type]

                print(f"{sub_type:>6}: {ref_norm:>11.6f} {expected_rate:>13.6f} {error_rate:>11.6f} {obs_exp_ratio:>14.6f}")

            print("-" * 70)
            print("Positive error rates indicate excess substitutions (damage/sequencing error)")
            print("Negative error rates indicate depleted substitutions relative to GTR expectation")
            print()

            # Create and display proper rate matrix
            prob_matrix = create_error_prob_matrix(error_matrix)
            #if rate_matrix:
            #    print("Error Rate Matrix (rows sum to zero):")
            #    print("Diagonal elements = -1 × (sum of off-diagonal elements in row)")
            #    print("-" * 60)

            #    bases = ['A', 'C', 'G', 'T']
            #    print(f"{'From/To':>8}", end="")
            #    for to_base in bases:
            #        print(f"{to_base:>12}", end="")
            #    print(f"{'Row Sum':>12}")
            #    print("-" * 60)

            #    for from_base in bases:
            #        print(f"{from_base:>8}", end="")
            #        row_sum = 0.0
            #        for to_base in bases:
            #            rate_key = f"{from_base}→{to_base}"
            #            rate_value = rate_matrix.get(rate_key, 0.0)
            #            print(f"{rate_value:>12.6f}", end="")
            #            row_sum += rate_value
            #        print(f"{row_sum:>12.6f}")

            #    print("-" * 60)
            #    print("Note: All row sums should be approximately 0.000000")
            #    print()

            #    # Calculate and display probability matrix
            #    prob_matrix = matrix_exponential(rate_matrix, time=1.0)
            if prob_matrix:
                bases = ['A', 'C', 'G', 'T']
                print("Probability Matrix P = exp(Rate Matrix × t) with t = 1.0:")
                print("Each row represents P(i→j) probabilities and should sum to 1.0")
                print("-" * 60)

                print(f"{'From/To':>8}", end="")
                for to_base in bases:
                    print(f"{to_base:>12}", end="")
                print(f"{'Row Sum':>12}")
                print("-" * 60)

                for from_base in bases:
                    print(f"{from_base:>8}", end="")
                    row_sum = 0.0
                    for to_base in bases:
                        prob_key = f"{from_base}→{to_base}"
                        prob_value = prob_matrix.get(prob_key, 0.0)
                        print(f"{prob_value:>12.6f}", end="")
                        row_sum += prob_value
                    print(f"{row_sum:>12.6f}")

                print("-" * 60)
                print("Note: All row sums should be approximately 1.000000")
                print("Diagonal elements represent P(no change)")
                print("Off-diagonal elements represent P(substitution)")
                print()
            else:
                print("Could not calculate probability matrix (missing numpy/scipy)")
                print()

    # Write detailed results to files for each dataset
    for dataset_name in ['all', 'forward', 'reverse']:
        dataset = data_sets[dataset_name]
        pos_data = position_data[dataset_name]

        # Generate file names with orientation suffix
        if dataset_name == 'all':
            suffix = ""
        else:
            suffix = f"_{dataset_name}"

        position_output_file = f"{output_prefix}{suffix}_position_analysis.txt"
        reverse_output_file = f"{output_prefix}{suffix}_reverse_position_analysis.txt"
        summary_output_file = f"{output_prefix}{suffix}_summary.txt"

        # Write position analysis
        with open(position_output_file, 'w') as f:
            f.write("Read_Position\tTotal_Substitutions\tSubstitution_Details\n")
            for pos in sorted(pos_data['position_subs'].keys()):
                sub_dict = pos_data['position_subs'][pos]
                total = sum(sub_dict.values())
                details = "; ".join([f"{sub_type}:{count}" for sub_type, count in sub_dict.items()])
                f.write(f"{pos}\t{total}\t{details}\n")

        # Write reverse position analysis
        with open(reverse_output_file, 'w') as f:
            f.write("Position_From_End\tTotal_Substitutions\tSubstitution_Details\n")
            for pos in sorted(pos_data['reverse_position_subs'].keys()):
                sub_dict = pos_data['reverse_position_subs'][pos]
                total = sum(sub_dict.values())
                details = "; ".join([f"{sub_type}:{count}" for sub_type, count in sub_dict.items()])
                f.write(f"{pos}\t{total}\t{details}\n")

        # Write summary
        with open(summary_output_file, 'w') as f:
            orientation_label = dataset_name.upper() if dataset_name != 'all' else 'ALL READS'
            f.write(f"SUBSTITUTION ANALYSIS SUMMARY - {orientation_label}\n")
            f.write("=" * 50 + "\n")
            f.write(f"Total substitutions: {dataset['total_substitutions']}\n")
            f.write(f"Total base comparisons: {dataset['total_comparisons']}\n")
            if dataset['total_comparisons'] > 0:
                f.write(f"Overall substitution rate: {dataset['total_substitutions']/dataset['total_comparisons']*100:.3f}%\n")
            f.write("\nObserved read base counts:\n")
            for base in ['A', 'T', 'C', 'G']:
                count = dataset['base_counts'].get(base, 0)
                percentage = count / dataset['total_comparisons'] * 100 if dataset['total_comparisons'] > 0 else 0
                f.write(f"  {base}: {count} ({percentage:.1f}%)\n")
            f.write("\nReference base counts:\n")
            for base in ['A', 'T', 'C', 'G']:
                count = dataset['ref_base_counts'].get(base, 0)
                percentage = count / dataset['total_comparisons'] * 100 if dataset['total_comparisons'] > 0 else 0
                f.write(f"  {base}: {count} ({percentage:.1f}%)\n")
            f.write("\nSubstitution type counts:\n")
            for sub_type, count in dataset['substitutions'].most_common():
                percentage = count / dataset['total_substitutions'] * 100 if dataset['total_substitutions'] > 0 else 0

                # Extract bases for both normalization methods
                ref_base, observed_base = sub_type.split('→')

                # Observed base normalization (target)
                observed_count = dataset['base_counts'].get(observed_base, 0)
                obs_normalized_rate = count / observed_count if observed_count > 0 else 0
                obs_per_thousand = obs_normalized_rate * 1000

                # Reference base normalization (source)
                ref_count = dataset['ref_base_counts'].get(ref_base, 0)
                ref_normalized_rate = count / ref_count if ref_count > 0 else 0
                ref_per_thousand = ref_normalized_rate * 1000

                # Obs/Exp ratio calculation if parameters available
                if gtr_params:
                    expected_rate = get_expected_gtr_rate_from_data(sub_type, gtr_params, dataset['ref_base_counts'])
                    obs_exp_ratio = ref_normalized_rate / expected_rate if expected_rate and expected_rate > 0 else float('inf')
                    f.write(f"  {sub_type}: {count} ({percentage:.1f}%, obs_norm: {obs_normalized_rate:.6f}, obs_per_1000: {obs_per_thousand:.3f}, ref_norm: {ref_normalized_rate:.6f}, ref_per_1000: {ref_per_thousand:.3f}, obs_exp_ratio: {obs_exp_ratio:.6f})\n")
                else:
                    f.write(f"  {sub_type}: {count} ({percentage:.1f}%, obs_norm: {obs_normalized_rate:.6f}, obs_per_1000: {obs_per_thousand:.3f}, ref_norm: {ref_normalized_rate:.6f}, ref_per_1000: {ref_per_thousand:.3f})\n")

    # Summary of generated files
    print(f"\nDetailed results saved to:")
    for dataset_name in ['all', 'forward', 'reverse']:
        suffix = "" if dataset_name == 'all' else f"_{dataset_name}"
        orientation_label = dataset_name.title() if dataset_name != 'all' else 'All reads'
        print(f"\n{orientation_label}:")
        print(f"  - Position analysis: {output_prefix}{suffix}_position_analysis.txt")
        print(f"  - Reverse position analysis: {output_prefix}{suffix}_reverse_position_analysis.txt")
        print(f"  - Summary: {output_prefix}{suffix}_summary.txt")

    # Return probability matrix and reference MSA for error profile calculations
    return prob_matrix, ref_msa

def print_read_examples(valid_reads_data, num_examples=3):
    """
    Print detailed examples of forward and reverse complement reads for analysis.

    Args:
        valid_reads_data: List of tuples (read, fastaSeq, sequence, ref_msa, ref_id, is_reverse_complement, quality)
        num_examples: Number of examples to show for each orientation
    """
    print("\n" + "=" * 80)
    print("DETAILED READ EXAMPLES FOR ANALYSIS")
    print("=" * 80)

    forward_examples = []
    reverse_examples = []

    # Collect examples
    for read_data in valid_reads_data:
        read, fastaSeq, sequence, ref_msa, ref_id, is_reverse_complement = read_data

        if is_reverse_complement and len(reverse_examples) < num_examples:
            reverse_examples.append(read_data)
        elif not is_reverse_complement and len(forward_examples) < num_examples:
            forward_examples.append(read_data)

        if len(forward_examples) >= num_examples and len(reverse_examples) >= num_examples:
            break

    # Print forward examples
    print(f"\nFORWARD READS ({len(forward_examples)} examples):")
    print("-" * 60)
    for i, (read, fastaSeq, sequence, ref_msa, ref_id, is_reverse_complement) in enumerate(forward_examples):
        print_single_read_details(read, fastaSeq, sequence, ref_msa, ref_id, is_reverse_complement, i+1)

    # Print reverse complement examples
    print(f"\nREVERSE COMPLEMENT READS ({len(reverse_examples)} examples):")
    print("-" * 60)
    for i, (read, fastaSeq, sequence, ref_msa, ref_id, is_reverse_complement) in enumerate(reverse_examples):
        print_single_read_details(read, fastaSeq, sequence, ref_msa, ref_id, is_reverse_complement, i+1)

def print_single_read_details(read, fastaSeq, sequence, ref_msa, ref_id, is_reverse_complement, example_num):
    """
    Print detailed information about a single read.
    """
    print(f"\nExample {example_num}: {read}")
    print(f"Reference ID: {ref_id}")
    print(f"Reverse complement: {is_reverse_complement}")
    print(f"Read length: {len(fastaSeq)}")

    # Find substitutions for this read
    substitutions, position_subs, reverse_position_subs, ref_base_counts = find_read_substitutions(
        fastaSeq, sequence, ref_msa, is_reverse_complement
    )

    print(f"Total substitutions: {len(substitutions)}")

    # For display purposes, show sequences in the orientation used for alignment
    # This means RC reads are shown as reverse complement for consistency with analysis
    if is_reverse_complement:
        # Show RC reads in their reverse complement form (how they're analyzed)
        from collections import deque
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-'}
        read_display = [complement_map.get(base.upper(), base.upper()) for base in reversed(fastaSeq)]
        display_note = " (showing reverse complement orientation used in analysis)"
    else:
        # Show forward reads in original form
        read_display = list(fastaSeq)
        display_note = " (showing original forward orientation)"

    print(f"Sequence orientation{display_note}")

    ref_display = []

    # Build reference sequence aligned to read (this is already in the correct orientation for analysis)
    for i, ref_pos in enumerate(sequence):
        if ref_pos == -1:  # Insertion in read
            if i < len(read_display):
                ref_display.append("-")
            continue

        if ref_pos < len(ref_msa):
            # Reference is always shown in the orientation used for alignment
            ref_base = ref_msa[ref_pos].upper()
            ref_display.append(ref_base)
        else:
            ref_display.append("?")

    # Pad to same length
    max_len = max(len(read_display), len(ref_display))
    while len(read_display) < max_len:
        read_display.append(" ")
    while len(ref_display) < max_len:
        ref_display.append(" ")

    # Show sequence alignment in chunks of 60
    chunk_size = 60
    for start in range(0, max_len, chunk_size):
        end = min(start + chunk_size, max_len)

        # Position numbers
        pos_line = ""
        for i in range(start, end):
            if i % 10 == 0:
                pos_line += str(i // 10)
            else:
                pos_line += " "
        print(f"Pos: {pos_line}")

        pos_line = ""
        for i in range(start, end):
            pos_line += str(i % 10)
        print(f"     {pos_line}")

        # Reference and read sequences
        ref_chunk = "".join(ref_display[start:end])
        read_chunk = "".join(read_display[start:end])

        print(f"Ref: {ref_chunk}")
        print(f"Read:{read_chunk}")

        # Mark substitutions
        diff_line = ""
        for i in range(start, end):
            if i < len(read_display) and i < len(ref_display):
                if read_display[i].upper() != ref_display[i].upper() and read_display[i] != "-" and ref_display[i] != "-" and read_display[i] != " " and ref_display[i] != " ":
                    diff_line += "^"
                else:
                    diff_line += " "
            else:
                diff_line += " "
        print(f"     {diff_line}")
        print()

    # List all substitutions with positions matching the display orientation
    if substitutions:
        print("Substitutions found (positions match display above):")
        for ref_pos, ref_base, read_base, read_pos, reverse_pos in substitutions[:10]:  # Show first 10
            if is_reverse_complement:
                # For RC reads, show positions as they appear in the RC display
                display_pos = len(fastaSeq) - 1 - read_pos
                display_reverse_pos = len(fastaSeq) - 1 - reverse_pos
                print(f"  Position {display_pos}: {ref_base}→{read_base} (from end: {display_reverse_pos}) [original pos: {read_pos}]")
            else:
                print(f"  Position {read_pos}: {ref_base}→{read_base} (from end: {reverse_pos})")
        if len(substitutions) > 10:
            print(f"  ... and {len(substitutions) - 10} more")
    else:
        print("No substitutions found")

    print("-" * 40)

# Parse GTR model parameters if provided
gtr_params = parse_gtr_parameters(args.gtrParams)

# Initialize data structures
seqInfo = {}
treesUsed = {}

# Initialize read counting variables
totalReadsInAssignment = 0
readsWithAssignments = 0
readsWithSamFlags = 0
readsPassingTaxaFilter = 0
readsWithCompleteInfo = 0
readsPassingMismatchFilter = 0

# seqInfo[read][0] tree
# seqInfo[read][1] assignment
# seqInfo[read][2] fastq quality scores
# seqInfo[read][3] fasta seq (with gaps)
# seqInfo[read][4] alignment ID (ref name)
# seqInfo[read][5] alignment (semicolon-separated positions)
# seqInfo[read][6] is_reverse_complement (boolean)

def parse_msa_fastq_header(header):
    """Parse header from MSA-aligned FASTQ: @read_name start=X len=Y ref=Z rc=0|1"""
    parts = header.lstrip('@').split()
    read_name = parts[0]
    start = None
    length = None
    ref = None
    rc = None

    for part in parts[1:]:
        if part.startswith('start='):
            start = int(part.split('=')[1])
        elif part.startswith('len='):
            length = int(part.split('=')[1])
        elif part.startswith('ref='):
            ref = part.split('=')[1]
        elif part.startswith('rc='):
            rc = int(part.split('=')[1]) == 1

    return read_name, start, length, ref, rc

# Read SAM/BAM file for RC flags if provided
samFlags = {}
if args.samFile:
    if args.debug:
        print(f"\n=== READING SAM/BAM FILE FOR RC FLAGS ===")
    is_bam = args.samFile.endswith('.bam')
    if is_bam:
        try:
            import pysam
            with pysam.AlignmentFile(args.samFile, 'rb') as fp:
                for read in fp:
                    samFlags[read.query_name] = bool(read.flag & 16)
        except ImportError:
            print("Error: pysam is required to read BAM files. Install with: pip install pysam")
            raise
    else:
        with open(args.samFile, 'r') as fp:
            for line in fp:
                if line.startswith('@'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    read_name = fields[0]
                    flag = int(fields[1])
                    samFlags[read_name] = bool(flag & 16)

# All reads use tree 0
treesUsed['0'] = 0

# Read the MSA-aligned FASTQ file
if args.debug:
    print(f"\n=== READING MSA-ALIGNED FASTQ: {args.inputReads} ===")

with open(args.inputReads, 'r') as fp:
    while True:
        header = fp.readline().strip()
        if not header:
            break
        if not header.startswith('@'):
            continue

        seq = fp.readline().strip()
        fp.readline()  # + line
        qual = fp.readline().strip()

        totalReadsInAssignment += 1
        read_name, start_pos, length, ref_name, rc_from_header = parse_msa_fastq_header(header)

        if start_pos is None or ref_name is None:
            continue

        readsWithAssignments += 1

        # Build alignment string: positions from start to start+len-1
        alignment_str = ';'.join(str(start_pos + i) for i in range(length)) + ';'

        # Use RC from FASTQ header if available, otherwise fall back to SAM file
        if rc_from_header is not None:
            is_rc = rc_from_header
            readsWithSamFlags += 1  # Count reads with RC info from header
        elif read_name in samFlags:
            is_rc = samFlags[read_name]
            readsWithSamFlags += 1  # Count reads with RC info from SAM
        else:
            is_rc = False  # Default to forward if no RC info available

        seqInfo[read_name] = [
            0,                              # tree (always 0)
            args.node,                      # assignment (user-specified node)
            qual,                           # quality scores
            seq,                            # sequence (with gaps)
            ref_name,                       # alignment ID
            alignment_str,                  # alignment positions
            is_rc                           # is_reverse_complement
        ]

if args.debug:
    print(f"Reads processed: {totalReadsInAssignment}, with valid data: {readsWithAssignments}")

# Apply taxa filter if provided
allowedTaxa = set()
if args.taxaFilter:
    with open(args.taxaFilter, 'r') as fp:
        for line in fp:
            line = line.strip()
            if line and not line.startswith('#'):
                ref_id = line.split()[0] if line.split() else line
                allowedTaxa.add(ref_id)
    if args.debug:
        print(f"Loaded {len(allowedTaxa)} allowed taxa from filter file")

    # Filter out reads not in allowed taxa
    reads_to_filter = [r for r in seqInfo if seqInfo[r][4] not in allowedTaxa]
    for r in reads_to_filter:
        del seqInfo[r]
    readsPassingTaxaFilter = len(seqInfo)
else:
    readsPassingTaxaFilter = len(seqInfo)

# Get the number of references in each tree
# Count number of taxa from reference MSA
with open(args.refMSA, "r") as f:
    numTaxa = sum(1 for line in f if line.startswith('>'))
treesUsed['0'] = numTaxa

# Remove reads with incomplete information and count the remaining ones
reads_to_remove = []
filtering_reasons = {}
if args.debug:
    print("\n=== DEBUGGING READ FILTERING ===")
    print(f"Total reads in seqInfo before filtering: {len(seqInfo)}")
    print(f"SAM file provided: {bool(args.samFile)}")

for read in seqInfo:
    reasons = []

    # Check each condition individually
    if len(seqInfo[read]) < 7:
        reasons.append(f"insufficient fields ({len(seqInfo[read])}/7)")

    if seqInfo[read][2] == "":
        reasons.append("missing quality scores")

    if seqInfo[read][3] == "":
        reasons.append("missing sequence")

    if seqInfo[read][4] == "":
        reasons.append("missing alignment ID")

    if seqInfo[read][5] == "":
        reasons.append("missing alignment")

    # If SAM file provided, also check for SAM information
    if args.samFile and read not in samFlags:
        reasons.append("not found in SAM file")

    if reasons:
        reads_to_remove.append(read)
        filtering_reasons[read] = reasons

# Print first few filtering reasons for debugging
if args.debug:
    print(f"\nReads to be removed: {len(reads_to_remove)}")
    if reads_to_remove:
        print("First 5 examples of filtering reasons:")
        for i, read in enumerate(reads_to_remove[:5]):
            print(f"  {read}: {', '.join(filtering_reasons[read])}")
            # Also print the actual seqInfo structure for this read
            print(f"    seqInfo structure: {seqInfo[read]}")

# Write removed reads report to file
if reads_to_remove:
    removed_reads_report_file = args.outMap.replace('.txt', '_removed_reads_report.txt')
    print(f"\nWriting removed reads report to: {removed_reads_report_file}")

    with open(removed_reads_report_file, 'w') as report_fp:
        # Write header
        report_fp.write(f"REMOVED READS REPORT\n")
        report_fp.write(f"{'='*80}\n")
        report_fp.write(f"Total reads removed: {len(reads_to_remove)}\n")
        report_fp.write(f"Total reads retained: {len(seqInfo) - len(reads_to_remove)}\n")
        report_fp.write(f"\n")

        # Count reasons
        reason_counts = {}
        for read in reads_to_remove:
            for reason in filtering_reasons[read]:
                reason_counts[reason] = reason_counts.get(reason, 0) + 1

        report_fp.write(f"Summary of filtering reasons:\n")
        report_fp.write(f"{'-'*80}\n")
        for reason, count in sorted(reason_counts.items(), key=lambda x: x[1], reverse=True):
            report_fp.write(f"  {reason}: {count} reads\n")
        report_fp.write(f"\n")

        # Write detailed list
        report_fp.write(f"Detailed list of removed reads:\n")
        report_fp.write(f"{'-'*80}\n")
        report_fp.write(f"{'Read Name':<50} {'Missing Information'}\n")
        report_fp.write(f"{'-'*80}\n")

        for read in sorted(reads_to_remove):
            reasons_str = ', '.join(filtering_reasons[read])
            report_fp.write(f"{read:<50} {reasons_str}\n")

    print(f"Removed reads report written successfully.")

# Remove incomplete reads
for read in reads_to_remove:
    del seqInfo[read]

if args.debug:
    print(f"\nReads remaining after filtering: {len(seqInfo)}")

readsWithCompleteInfo = len(seqInfo)
fullInfoReads = len(seqInfo)

readNum = 0

# Initialize global probability matrix and reference MSA for error profile calculations
global_prob_matrix = None
global_ref_msa = None

# If GTR parameters are provided and substitution analysis is enabled,
# pre-calculate the probability matrix for error profile calculations
if args.gtrParams and args.substitution_analysis:
    print("\nPre-calculating probability matrix for error profile calculations...")

    # This is a placeholder - we'll calculate it properly after gathering substitution data
    # For now, initialize as None and calculate during substitution analysis
    pass

# MEMORY OPTIMIZATION: Load MSA file once into dictionary
print("\nLoading MSA reference sequences into memory...")
msa_sequences = {}

def load_msa_sequences():
    """Load all MSA reference sequences into a dictionary for fast lookup"""
    try:
        with gzip.open(args.refMSA, "rt") as mp:
            return load_msa_from_file(mp)
    except (gzip.BadGzipFile, OSError):
        # File is not gzipped, open as plain text
        with open(args.refMSA, "r") as mp:
            return load_msa_from_file(mp)

def load_msa_from_file(fp):
    """Load MSA sequences from file handle (handles multi-line FASTA)"""
    sequences = {}
    current_header = None
    current_seq = []

    for line in fp:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            # Save previous sequence if exists
            if current_header is not None:
                sequences[current_header] = ''.join(current_seq)
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)

    # Save last sequence
    if current_header is not None:
        sequences[current_header] = ''.join(current_seq)

    return sequences

msa_sequences = load_msa_sequences()
print(f"Loaded {len(msa_sequences)} MSA reference sequences")

# Look up MSA sequences by exact name: the header up to the first whitespace, which is how
# filter_msa_by_sam.py writes ref= in the FASTQ headers
msa_by_name = {header[1:].split()[0]: sequence for header, sequence in msa_sequences.items()}

if args.debug:
    print(f"First 5 MSA headers: {list(msa_sequences.keys())[:5]}")

# Process reads and collect valid ones for output
validReads = []
validReadsData = []  # For substitution analysis
mismatchedReads = []  # For FASTA output of reads with mismatches
for read in seqInfo.keys():
    # Look up MSA reference sequence from pre-loaded dictionary
    ref_id = seqInfo[read][4]
    ref_msa = None

    #print(msa_sequences.items())

    # Find matching sequence in MSA dictionary
    ref_msa = msa_by_name.get(ref_id)

    if not ref_msa:
        print(f"Error: Ref {ref_id} missing for {read}, skipping")
        continue

    # Process alignment for mismatch checking
    sequence = seqInfo[read][5].split(";")  # convert to list
    sequence = [i for i in sequence if i != '']     # remove empty strings
    sequence = [int(i) for i in sequence]   # convert to int list

    # Find first and last non-insertion positions for trimming
    firstPos = 0
    lastPos = len(sequence) - 1
    for i in range(len(sequence)):
        if sequence[i] != -1:
            firstPos = i
            break
    for i in range(len(sequence) - 1, -1, -1):
        if sequence[i] != -1:
            lastPos = i
            break

    # Get RC flag for damage masking and positional analysis
    # NOTE: The sequence in the FASTQ is already in alignment orientation (RC reads are already
    # reverse complemented by the aligner). Do NOT reverse complement again - just use as-is.
    # The alignment positions (start=X) were calculated for this orientation.
    is_reverse_complement = seqInfo[read][6]
    fastaSeq = list(seqInfo[read][3])

    # Apply trimming from both ends (counting actual bases, not gaps)
    trim_amount = args.trim_ends
    if trim_amount > 0:
        # Find new firstPos after trimming (skip trim_amount actual bases from start)
        trimmed_count = 0
        for i in range(firstPos, lastPos + 1):
            if i < len(fastaSeq) and fastaSeq[i] != '-':
                trimmed_count += 1
                if trimmed_count > trim_amount:
                    firstPos = i
                    break
        else:
            continue  # Read too short after trimming

        # Find new lastPos after trimming (skip trim_amount actual bases from end)
        trimmed_count = 0
        for i in range(lastPos, firstPos - 1, -1):
            if i < len(fastaSeq) and fastaSeq[i] != '-':
                trimmed_count += 1
                if trimmed_count > trim_amount:
                    lastPos = i
                    break
        else:
            continue  # Read too short after trimming

        if firstPos >= lastPos:
            continue  # Nothing left after trimming

    # Slice to trimmed range
    fastaSeq = fastaSeq[firstPos:lastPos + 1]
    sequence = sequence[firstPos:lastPos + 1]
    quality = seqInfo[read][2][firstPos:lastPos + 1]

    # Apply masking
    strand_type = 'singlestrand' if args.singlestrand else 'doublestrand'

    # Apply transition masking
    fastaSeq = apply_transition_masking(fastaSeq, sequence, ref_msa)

    # Apply damage masking (deamination and/or artifacts)
    fastaSeq = apply_damage_masking(fastaSeq, sequence, ref_msa, args.damage_masking, args.damage_masking_bases, strand_type, is_reverse_complement, args.damage_types)

    # Trim ends if they start/end with - and apply trimming to quality scores
    start = 0
    while start < len(fastaSeq) and fastaSeq[start] == "-":
        start += 1

    end = len(fastaSeq)
    while end > start and fastaSeq[end - 1] == "-":
        end -= 1

    fastaSeq = fastaSeq[start:end]
    quality = quality[start:end]
    sequence = sequence[start:end]

    if readNum == 1271:
        print(f"{fastaSeq}\n{quality}\n{sequence}")

    # Count mismatches after all masking is applied
    mismatch_count = count_mismatches_after_masking(fastaSeq, sequence, ref_msa)

    # Collect reads with mismatches (if any) for FASTA output
    if args.output_mismatched_fasta and mismatch_count > 0:
        mismatchedReads.append((read, seqInfo[read][3], mismatch_count))  # (read_name, original_sequence, mismatch_count)

    # Check if read passes mismatch threshold
    if mismatch_count <= args.max_mismatches:
        validReads.append(read)
        # Note: Masking has already been applied above, so we just need to copy the data
        # The fastaSeq already has transition and damage masking applied

        validReadsData.append((read, fastaSeq, sequence, ref_msa, ref_id, is_reverse_complement, quality))

readsPassingMismatchFilter = len(validReads)

# Debug: Check RC distribution in validReadsData after collection
if args.debug and args.substitution_analysis:
    debug_rc_count = sum(1 for item in validReadsData if item[5])
    debug_fwd_count = sum(1 for item in validReadsData if not item[5])
    print(f"\nDEBUG after validation loop: validReadsData size={len(validReadsData)}, forward={debug_fwd_count}, reverse={debug_rc_count}")
    # Check first 5 entries
    for i, item in enumerate(validReadsData[:5]):
        print(f"  Entry {i}: read={item[0][:30]}..., rc={item[5]}")

# Write output files
with open(args.outAlign, "w") as fp:     # alignment file
    with open(args.outAssign, "w") as gp:    # assignment file
        with open(args.outMap, "w") as tp:   # read mapping file
            with open(args.outError, "w") as ep:  # error profile file
                fp.write(str(readsPassingMismatchFilter) + "\n")  # num assignments
                gp.write(str(readsPassingMismatchFilter) + "\n")

                for read in validReadsData:
                    # Process only reads that passed mismatch filtering

                    # Write read mapping: TO DO: fix seqInfo.keys(read[0]) is my guess
                    tp.write(read[0] + "\t" + str(readNum) + "\n")

                    # Write assignment
                    gp.write(str(seqInfo[read[0]][0]) + " " + str(seqInfo[read[0]][1]) + "\n")

                    sequence = read[2]
                    is_reverse_complement = read[5]
                    quality = read[6]
                    fastaSeq = read[1]

                    # Apply masking
                    strand_type = 'singlestrand' if args.singlestrand else 'doublestrand'

                    # Generate error profile only if we don't have GTR params (phylogenetic calculations)
                    # or a pre-computed damage matrix — both cases reprocess later.
                    if not (args.gtrParams and args.substitution_analysis) and not args.damage_matrix:
                        error_entries = generate_error_profile(fastaSeq, quality, sequence, readNum, args.encoding, args.error_profile_bases,
                                                               is_reverse_complement=is_reverse_complement)

                        # Write error profile entries
                        for entry in error_entries:
                            ep.write(entry + "\n")
                    else:
                        # Skip error profile generation - will be done later with phylogenetic calculations
                        pass

                    # Write alignment using trimmed data
                    non_ins = [pos for pos in sequence if pos != -1]
                    align_len = (non_ins[-1] - non_ins[0] + 1) if non_ins else len(fastaSeq)
                    fp.write(f"{align_len}\t")    # length of read to reference including deletion gaps
                    fp.write(str(sequence[0]) + "\t")    # starting position
                    fp.write(fastaSeq[0])   # first base

                    lastBase = 0
                    for i in range(1, len(sequence)):
                        if sequence[i] == -1:
                            continue    # skip insertions
                        gap = sequence[i] - sequence[lastBase]
                        if gap > 1:
                            fp.write("-" * (gap - 1))  # add deletion gaps
                        fp.write(fastaSeq[i])
                        lastBase = i
                    fp.write("\n")

                    readNum = readNum + 1

print("\n=== READ FILTERING SUMMARY ===")
print(f"Total reads in assignment file: {totalReadsInAssignment}")
print(f"Reads with assignments (not 'unassigned'): {readsWithAssignments} (-{totalReadsInAssignment - readsWithAssignments} removed)")
if args.samFile:
    print(f"Reads with SAM flags: {readsWithSamFlags} (-{readsWithAssignments - readsWithSamFlags} removed)")
else:
    print(f"RC detection from headers: {readsWithSamFlags} (no SAM file provided)")
if args.taxaFilter:
    print(f"Reads passing taxa filter: {readsPassingTaxaFilter} (-{readsWithSamFlags - readsPassingTaxaFilter} removed)")
    print(f"  - Allowed taxa: {len(allowedTaxa)}")
else:
    readsPassingTaxaFilter = readsWithSamFlags  # No taxa filtering applied
    print(f"Reads passing taxa filter: {readsPassingTaxaFilter} (no taxa filtering)")
print(f"Reads with complete information: {readsWithCompleteInfo} (-{readsPassingTaxaFilter - readsWithCompleteInfo} removed)")
print(f"Reads passing mismatch filter (≤{args.max_mismatches} mismatches): {readsPassingMismatchFilter} (-{readsWithCompleteInfo - readsPassingMismatchFilter} removed)")
print(f"Final reads processed: {readsPassingMismatchFilter}")

print(f"\n=== PROCESSING SETTINGS ===")
print(f"Input file: {args.inputReads}")
print(f"Node assignment: {args.node} (tree 0)")
rc_method = "FASTQ header (rc= field)"
if args.samFile:
    rc_method += " with SAM file fallback"
print(f"RC detection method: {rc_method}")
print(f"Taxa filtering: {'Enabled' if args.taxaFilter else 'Disabled'}")
print(f"Mismatch filtering: Enabled (max {args.max_mismatches} mismatches after masking)")
print(f"Read trimming: {f'Enabled ({args.trim_ends} bases from each end)' if args.trim_ends > 0 else 'Disabled'}")
print(f"Transition masking: {'Enabled (A<->G, C<->T masked)' if args.mask_transitions else 'Disabled'}")
if args.damage_masking != 'none':
    if args.damage_types == 'deamination':
        patterns = "C→T, G→A"
    elif args.damage_types == 'artifacts':
        patterns = "G→T, C→A"
    else:  # both
        patterns = "C→T, G→A, G→T, C→A"
    if args.damage_masking == 'all':
        damage_status = f"Enabled ({patterns} masked across full sequence)"
    else:  # ends
        damage_status = f"Enabled ({patterns} masked at {args.damage_masking_bases} bases from each end)"
    print(f"Damage masking: {damage_status}")
else:
    print(f"Damage masking: Disabled")
error_profile_status = "Disabled" if args.no_error_profile else f"Enabled from quality scores ({args.error_profile_bases if args.error_profile_bases else 'all bases'})"
print(f"Error profile generation: {error_profile_status}")
print(f"Strand type: {'Single-stranded' if args.singlestrand else 'Double-stranded'}")
print(f"Substitution analysis: {'Enabled' if args.substitution_analysis else 'Disabled'}")

# Count reverse complement reads from valid reads only
rc_count = sum(1 for read in validReads if len(seqInfo[read]) >= 7 and seqInfo[read][6])
print(f"Reverse complement reads: {rc_count} out of {readNum} final reads")

# Write mismatched reads to FASTA file if requested
if args.output_mismatched_fasta:
    if mismatchedReads:
        with open(args.output_mismatched_fasta, 'w') as fasta_fp:
            for read_name, original_sequence, mismatch_count in mismatchedReads:
                write_mismatched_read_to_fasta(read_name, original_sequence, mismatch_count, fasta_fp)
        print(f"\n=== MISMATCHED READS OUTPUT ===")
        print(f"Wrote {len(mismatchedReads)} reads with mismatches to: {args.output_mismatched_fasta}")
    else:
        print(f"\nNo reads with mismatches found to write to FASTA file.")

# Use pre-computed damage matrix if provided; otherwise fall back to substitution analysis
if args.damage_matrix:
    if validReadsData:
        damage_npy, n_positions = load_damage_matrix_from_npy(args.damage_matrix)
        positional_matrices_arrays = build_positional_matrices_from_npy(damage_npy, n_positions)
        fallback_array = damage_npy[n_positions]  # middle slice as global fallback
        write_error_profiles_with_damage_matrix(
            validReadsData, positional_matrices_arrays, fallback_array, n_positions
        )
    else:
        print("\nNo valid reads to process with damage matrix.")

# Perform substitution analysis if requested (skipped when --damage-matrix is provided)
elif args.substitution_analysis:
    if validReadsData:
        # Debug: Check RC distribution in validReadsData
        rc_in_data = sum(1 for item in validReadsData if item[5])
        fwd_in_data = sum(1 for item in validReadsData if not item[5])
        print(f"\nDEBUG: validReadsData RC distribution: forward={fwd_in_data}, reverse={rc_in_data}")
        # Generate output prefix based on output files
        import os
        output_prefix = os.path.splitext(args.outAlign)[0] + "_substitutions"
        prob_matrix, ref_msa_from_analysis, positional_prob_matrices = perform_substitution_analysis(validReadsData, output_prefix)

        # If we have a probability matrix, reprocess error profiles with phylogenetic calculations
        if prob_matrix is not None and ref_msa_from_analysis is not None:
            print("\nReprocessing error profiles with phylogenetic conditional likelihood calculations...")
            print(f"DEBUG: prob_matrix type: {type(prob_matrix)}")
            print(f"DEBUG: ref_msa_from_analysis type: {type(ref_msa_from_analysis)}")
            print(f"DEBUG: positional_prob_matrices type: {type(positional_prob_matrices)}")
            reprocess_error_profiles_with_phylogenetic_calculations(validReadsData, prob_matrix, ref_msa_from_analysis, positional_prob_matrices)
        else:
            print(f"\nDEBUG: Not reprocessing - prob_matrix: {prob_matrix is not None}, ref_msa: {ref_msa_from_analysis is not None}")
    else:
        print("\nNo valid reads available for substitution analysis.")
