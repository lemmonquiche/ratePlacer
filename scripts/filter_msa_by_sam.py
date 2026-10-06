#!/usr/bin/env python3
"""
Filter MSA to positions where SAM/BAM alignments occur.

This script takes a SAM/BAM file with alignments to ungapped references and an MSA
file containing the same sequences with gaps, then outputs a filtered MSA containing
only the positions where reads aligned.
"""

import argparse
import gzip
import sys
import re
import array
from collections import defaultdict

try:
    import pysam
    _PYSAM_AVAILABLE = True
except ImportError:
    _PYSAM_AVAILABLE = False


def iter_alignments(path):
    """
    Yield alignment tuples from a SAM or BAM file.

    Each tuple: (read_name, flag, ref_name, ref_pos_0based, cigar, seq, qual_str)
    qual_str is None when qualities are absent.

    Requires pysam for BAM files. SAM files fall back to plain-text parsing
    when pysam is unavailable.
    """
    if path.endswith('.bam'):
        if not _PYSAM_AVAILABLE:
            print("Error: pysam is required to read BAM files. Install with: pip install pysam",
                  file=sys.stderr)
            sys.exit(1)
        with pysam.AlignmentFile(path, 'rb') as af:
            for read in af.fetch(until_eof=True):
                cigar = read.cigarstring if read.cigarstring else '*'
                seq = read.query_sequence if read.query_sequence else ''
                if read.query_qualities is not None:
                    qual = ''.join(chr(q + 33) for q in read.query_qualities)
                else:
                    qual = None
                yield (read.query_name, read.flag,
                       read.reference_name if read.reference_name else '*',
                       read.reference_start,  # already 0-based
                       cigar, seq, qual)
    else:
        # Plain SAM (or pysam-assisted SAM)
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue
                qual = fields[10] if fields[10] != '*' else None
                yield (fields[0], int(fields[1]), fields[2],
                       int(fields[3]) - 1,  # convert to 0-based
                       fields[5], fields[9], qual)


def parse_cigar(cigar_string):
    """
    Parse CIGAR string to get alignment length on reference.

    Args:
        cigar_string: CIGAR string (e.g., "36M", "12M1D24M")

    Returns:
        int: Length of alignment on reference
    """
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    ref_length = 0

    for length, op in ops:
        length = int(length)
        if op in 'MDN=X':  # Operations that consume reference
            ref_length += length

    return ref_length


def parse_cigar_alignment(cigar_string, ref_start, seq, qual=None):
    """
    Parse CIGAR string and return list of (ref_pos, base, qual) tuples.

    Only returns positions that consume reference (M, D, N, =, X).
    Insertions are skipped (not added to output).
    Deletions are represented as '-' with '!' quality (lowest phred33).

    Args:
        cigar_string: CIGAR string
        ref_start: 0-based reference start position
        seq: Read sequence
        qual: Quality string (optional, same length as seq)

    Returns:
        list: List of (ref_pos, base, qual) tuples
    """
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)

    result = []
    ref_pos = ref_start
    seq_pos = 0

    for length_str, op in ops:
        length = int(length_str)

        if op in 'M=X':  # Alignment match/mismatch
            for _ in range(length):
                if seq_pos < len(seq):
                    q = qual[seq_pos] if qual and seq_pos < len(qual) else '!'
                    result.append((ref_pos, seq[seq_pos], q))
                ref_pos += 1
                seq_pos += 1
        elif op == 'I':  # Insertion - skip (don't add to output)
            seq_pos += length
        elif op in 'DN':  # Deletion/skip - consume ref, output gap
            for _ in range(length):
                result.append((ref_pos, '-', '!'))  # Gap with lowest quality
                ref_pos += 1
        elif op == 'S':  # Soft clip
            seq_pos += length
        # H and P don't consume either

    return result


def build_position_map(msa_seq, verbose=False):
    """
    Build mapping from reference position (ungapped) to MSA position (with gaps).

    Args:
        msa_seq: MSA sequence string (may contain gaps)
        verbose: Print progress messages

    Returns:
        array: Array mapping reference position to MSA position
    """
    if verbose:
        print(f"    Building position map...", file=sys.stderr)

    positions = array.array('L')  # unsigned long

    for msa_pos, char in enumerate(msa_seq):
        if char != '-':  # Non-gap character
            positions.append(msa_pos)

        if verbose and msa_pos > 0 and msa_pos % 50000000 == 0:
            print(f"      ... {msa_pos:,} / {len(msa_seq):,} positions", file=sys.stderr)

    # Add sentinel value for end
    positions.append(len(msa_seq))

    if verbose:
        print(f"    Position map built: {len(positions):,} reference positions", file=sys.stderr)

    return positions


def bisect_left(arr, x):
    """Binary search to find leftmost position where x could be inserted."""
    lo, hi = 0, len(arr)
    while lo < hi:
        mid = (lo + hi) // 2
        if arr[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


def extract_reads_to_filtered_msa(sam_file, msa_sequences, final_msa_positions,
                                   output_file, verbose=False, trim_5=0, trim_3=0):
    """
    Extract reads from SAM and convert to filtered MSA coordinates.

    Outputs reads in FASTQ format with headers containing start position and
    length (including gaps from MSA). Insertions relative to reference are removed.
    Gaps are assigned '!' quality (lowest phred33 score).

    Memory-optimized: uses binary search instead of large dictionaries.

    Args:
        sam_file: Path to SAM file
        msa_sequences: Dict of original MSA sequences
        final_msa_positions: Sorted list of original MSA positions in the final output
        output_file: Path to output FASTQ file
        verbose: Print verbose output
        trim_5: Bases to trim from 5' end
        trim_3: Bases to trim from 3' end
    """
    # Convert to list for binary search (already sorted)
    final_positions_list = list(final_msa_positions)

    # Cache for position maps only (ref_pos -> msa_pos)
    position_maps = {}

    total_reads = 0
    extracted_reads = 0
    skipped_refs = set()

    if verbose:
        print(f"\nExtracting reads to filtered MSA coordinates (FASTQ)...", file=sys.stderr)

    with open(output_file, 'w') as out_f:
        for read_name, flag, ref_name, ref_pos, cigar, seq, qual in iter_alignments(sam_file):
            is_reverse_complement = 1 if (flag & 16) else 0

            total_reads += 1

            # Skip unmapped
            if flag & 4 or ref_name == '*' or cigar == '*':
                continue

            if ref_name not in msa_sequences:
                if ref_name not in skipped_refs:
                    skipped_refs.add(ref_name)
                continue

            # Build position map if not cached (ref_pos -> msa_pos)
            if ref_name not in position_maps:
                position_maps[ref_name] = build_position_map(msa_sequences[ref_name], False)

            pos_map = position_maps[ref_name]

            # Parse CIGAR to get (ref_pos, base, qual) tuples
            alignment = parse_cigar_alignment(cigar, ref_pos, seq, qual)

            if not alignment:
                continue

            # Apply trimming
            if trim_5 + trim_3 >= len(alignment):
                continue

            if trim_3 > 0:
                trimmed = alignment[trim_5:-trim_3]
            else:
                trimmed = alignment[trim_5:]

            if not trimmed:
                continue

            # Get MSA position range for this read
            first_ref_pos = trimmed[0][0]
            last_ref_pos = trimmed[-1][0]

            if first_ref_pos >= len(pos_map) - 1 or last_ref_pos >= len(pos_map) - 1:
                continue

            msa_start = pos_map[first_ref_pos]
            msa_end = pos_map[last_ref_pos]

            # Build ref_pos -> (base, qual) mapping for this read
            ref_to_base_qual = {rp: (base, q) for rp, base, q in trimmed}

            # Find the range of final positions that overlap with this read's MSA range
            # Using binary search to find start and end indices
            start_idx = bisect_left(final_positions_list, msa_start)
            end_idx = bisect_left(final_positions_list, msa_end + 1)

            if start_idx >= end_idx:
                continue

            # Build output sequence and quality
            output_bases = []
            output_quals = []
            first_final_pos = start_idx

            # Track current reference position as we iterate through MSA
            # We use the pos_map to check if an MSA position corresponds to a ref position
            for final_idx in range(start_idx, end_idx):
                orig_msa_pos = final_positions_list[final_idx]

                # Use binary search on pos_map to find if this msa_pos has a ref_pos
                # pos_map[ref_pos] = msa_pos, so we search for msa_pos in pos_map values
                # Since pos_map is sorted by msa_pos values, we can binary search
                ref_idx = bisect_left(pos_map, orig_msa_pos)

                if ref_idx < len(pos_map) and pos_map[ref_idx] == orig_msa_pos:
                    # This MSA position corresponds to ref_idx
                    if ref_idx in ref_to_base_qual:
                        base, q = ref_to_base_qual[ref_idx]
                        output_bases.append(base)
                        output_quals.append(q)
                    else:
                        # Reference position exists but not in read's range
                        output_bases.append('-')
                        output_quals.append('!')
                else:
                    # Gap in this reference at this MSA position
                    output_bases.append('-')
                    output_quals.append('!')

            if not output_bases:
                continue

            output_seq = ''.join(output_bases)
            output_qual = ''.join(output_quals)
            gapped_length = len(output_seq)

            # Write FASTQ entry
            # Header format: @read_name start=X len=Y ref=Z rc=0|1
            out_f.write(f'@{read_name} start={first_final_pos} len={gapped_length} ref={ref_name} rc={is_reverse_complement}\n')
            out_f.write(f'{output_seq}\n')
            out_f.write('+\n')
            out_f.write(f'{output_qual}\n')

            extracted_reads += 1

            if verbose and extracted_reads % 10000 == 0:
                print(f"  Extracted {extracted_reads:,} reads...", file=sys.stderr)

    if verbose:
        print(f"\nRead extraction complete:", file=sys.stderr)
        print(f"  Total reads in SAM: {total_reads:,}", file=sys.stderr)
        print(f"  Reads extracted: {extracted_reads:,}", file=sys.stderr)
        print(f"  Output written to: {output_file}", file=sys.stderr)


def parse_fasta(fasta_file):
    """
    Parse FASTA file and return sequences.
    Handles both single-line and multi-line (wrapped) sequences.
    Accepts plain or gzip-compressed files (.gz).

    Args:
        fasta_file: Path to FASTA file (plain or .gz)

    Returns:
        dict: Dictionary mapping sequence IDs to sequences
    """
    sequences = {}
    current_id = None
    current_seq = []

    opener = gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file, 'r')
    with opener as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            elif current_id:
                current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = ''.join(current_seq)

    return sequences


def collect_msa_positions_with_depth(sam_file, msa_sequences, verbose=False, trim_5=0, trim_3=0):
    """
    Parse SAM file and collect MSA positions with depth information.

    Args:
        sam_file: Path to SAM file
        msa_sequences: Dict of MSA sequences
        verbose: Print verbose output
        trim_5: Number of bases to trim from 5' end of each read
        trim_3: Number of bases to trim from 3' end of each read

    Returns:
        tuple: (msa_positions set, depth_per_position dict, stats dict)
    """
    msa_positions = set()
    depth_per_position = defaultdict(int)  # MSA position -> read count
    position_maps = {}

    total_reads = 0
    mapped_reads = 0
    trimmed_reads = 0
    total_bases = 0  # Total aligned bases (after trimming)
    skipped_refs = set()

    # Track per-reference stats
    ref_stats = defaultdict(lambda: {'reads': 0, 'bases': 0, 'positions': set()})

    for read_name, flag, ref_name, ref_pos, cigar, seq, qual in iter_alignments(sam_file):
        total_reads += 1

        # Skip unmapped reads
        if flag & 4 or ref_name == '*':
            continue

        # Check if reference exists in MSA
        if ref_name not in msa_sequences:
            if ref_name not in skipped_refs:
                if verbose:
                    print(f"Warning: {ref_name} not found in MSA, skipping alignments to it",
                          file=sys.stderr)
                skipped_refs.add(ref_name)
            continue

        # Build position map if not cached
        if ref_name not in position_maps:
            if verbose:
                print(f"  Processing reference: {ref_name}", file=sys.stderr)
                msa_len = len(msa_sequences[ref_name])
                ungapped = len(msa_sequences[ref_name].replace('-', ''))
                print(f"    MSA length: {msa_len:,} bp, ungapped: {ungapped:,} bp",
                      file=sys.stderr)

            position_maps[ref_name] = build_position_map(msa_sequences[ref_name], verbose)

        # Get alignment length from CIGAR
        ref_length = parse_cigar(cigar)

        # Apply trimming
        trimmed_ref_pos = ref_pos + trim_5
        trimmed_ref_length = ref_length - trim_5 - trim_3

        # Skip if trimming removes entire read
        if trimmed_ref_length <= 0:
            trimmed_reads += 1
            continue

        # Map to MSA positions
        pos_map = position_maps[ref_name]

        # Check bounds
        if trimmed_ref_pos >= len(pos_map):
            continue

        end_ref_pos = min(trimmed_ref_pos + trimmed_ref_length, len(pos_map) - 1)

        # Skip if end position is before start (shouldn't happen but safety check)
        if end_ref_pos <= trimmed_ref_pos:
            continue

        msa_start = pos_map[trimmed_ref_pos]
        msa_end = pos_map[end_ref_pos]

        # Add all positions in this range and track depth
        read_bases = 0
        for pos in range(msa_start, msa_end):
            msa_positions.add(pos)
            depth_per_position[pos] += 1
            ref_stats[ref_name]['positions'].add(pos)
            read_bases += 1

        total_bases += read_bases
        ref_stats[ref_name]['reads'] += 1
        ref_stats[ref_name]['bases'] += read_bases
        mapped_reads += 1

        # Progress indicator
        if verbose and mapped_reads % 10000 == 0:
            print(f"  Processed {mapped_reads:,} alignments, {len(msa_positions):,} MSA positions",
                  file=sys.stderr)

    # Compile statistics
    stats = {
        'total_reads': total_reads,
        'mapped_reads': mapped_reads,
        'trimmed_reads': trimmed_reads,
        'total_bases': total_bases,
        'unique_positions': len(msa_positions),
        'ref_stats': dict(ref_stats)
    }

    if verbose:
        print(f"\nTotal reads in SAM: {total_reads:,}", file=sys.stderr)
        print(f"Mapped reads processed: {mapped_reads:,}", file=sys.stderr)
        if trim_5 > 0 or trim_3 > 0:
            print(f"Reads trimmed away completely: {trimmed_reads:,}", file=sys.stderr)
            print(f"Trimming applied: {trim_5} bp from 5' end, {trim_3} bp from 3' end", file=sys.stderr)
        print(f"Unique MSA positions with alignments: {len(msa_positions):,}", file=sys.stderr)

    return msa_positions, depth_per_position, stats


def collect_msa_positions(sam_file, msa_sequences, verbose=False, trim_5=0, trim_3=0):
    """
    Parse SAM file and collect all MSA positions that have alignments.

    Args:
        sam_file: Path to SAM file
        msa_sequences: Dict of MSA sequences
        verbose: Print verbose output
        trim_5: Number of bases to trim from 5' end of each read
        trim_3: Number of bases to trim from 3' end of each read

    Returns:
        set: Set of MSA positions (0-based) that have alignments
    """
    msa_positions = set()
    position_maps = {}

    total_reads = 0
    mapped_reads = 0
    trimmed_reads = 0
    skipped_refs = set()

    for read_name, flag, ref_name, ref_pos, cigar, seq, qual in iter_alignments(sam_file):
        total_reads += 1

        # Skip unmapped reads
        if flag & 4 or ref_name == '*':
            continue

        # Check if reference exists in MSA
        if ref_name not in msa_sequences:
            if ref_name not in skipped_refs:
                if verbose:
                    print(f"Warning: {ref_name} not found in MSA, skipping alignments to it",
                          file=sys.stderr)
                skipped_refs.add(ref_name)
            continue

        # Build position map if not cached
        if ref_name not in position_maps:
            if verbose:
                print(f"  Processing reference: {ref_name}", file=sys.stderr)
                msa_len = len(msa_sequences[ref_name])
                ungapped = len(msa_sequences[ref_name].replace('-', ''))
                print(f"    MSA length: {msa_len:,} bp, ungapped: {ungapped:,} bp",
                      file=sys.stderr)

            position_maps[ref_name] = build_position_map(msa_sequences[ref_name], verbose)

        # Get alignment length from CIGAR
        ref_length = parse_cigar(cigar)

        # Apply trimming
        trimmed_ref_pos = ref_pos + trim_5
        trimmed_ref_length = ref_length - trim_5 - trim_3

        # Skip if trimming removes entire read
        if trimmed_ref_length <= 0:
            trimmed_reads += 1
            continue

        # Map to MSA positions
        pos_map = position_maps[ref_name]

        # Check bounds
        if trimmed_ref_pos >= len(pos_map):
            continue

        end_ref_pos = min(trimmed_ref_pos + trimmed_ref_length, len(pos_map) - 1)

        # Skip if end position is before start (shouldn't happen but safety check)
        if end_ref_pos <= trimmed_ref_pos:
            continue

        msa_start = pos_map[trimmed_ref_pos]
        msa_end = pos_map[end_ref_pos]

        # Add all positions in this range
        for pos in range(msa_start, msa_end):
            msa_positions.add(pos)

        mapped_reads += 1

        # Progress indicator
        if verbose and mapped_reads % 10000 == 0:
            print(f"  Processed {mapped_reads:,} alignments, {len(msa_positions):,} MSA positions",
                  file=sys.stderr)

    if verbose:
        print(f"\nTotal reads in SAM: {total_reads:,}", file=sys.stderr)
        print(f"Mapped reads processed: {mapped_reads:,}", file=sys.stderr)
        if trim_5 > 0 or trim_3 > 0:
            print(f"Reads trimmed away completely: {trimmed_reads:,}", file=sys.stderr)
            print(f"Trimming applied: {trim_5} bp from 5' end, {trim_3} bp from 3' end", file=sys.stderr)
        print(f"Unique MSA positions with alignments: {len(msa_positions):,}", file=sys.stderr)

    return msa_positions


def filter_msa_by_positions(sequences, positions, verbose=False):
    """
    Filter MSA sequences to only include specified positions.

    Args:
        sequences: Dictionary mapping sequence IDs to sequences
        positions: Set of positions (0-based) to include
        verbose: Print progress

    Returns:
        tuple: (filtered_sequences dict, sorted_positions list)
    """
    if verbose:
        print(f"\nFiltering MSA to {len(positions):,} positions...", file=sys.stderr)

    filtered_sequences = {}
    sorted_positions = sorted(positions)

    for seq_id, seq in sequences.items():
        filtered_seq = []
        for pos in sorted_positions:
            if pos < len(seq):
                filtered_seq.append(seq[pos])
            else:
                filtered_seq.append('-')

        filtered_sequences[seq_id] = ''.join(filtered_seq)

        if verbose:
            print(f"  {seq_id}: {len(seq):,} bp -> {len(filtered_seq):,} bp", file=sys.stderr)

    return filtered_sequences, sorted_positions


def remove_all_gap_columns(sequences, verbose=False):
    """
    Remove columns that contain only gaps or undetermined values across all sequences.

    Args:
        sequences: Dictionary mapping sequence IDs to sequences
        verbose: Print progress

    Returns:
        tuple: (cleaned_sequences dict, positions_to_keep list of indices kept)
    """
    if not sequences:
        return sequences, []

    # Characters treated as gaps or undetermined (missing data)
    # These should match what downstream phylogenetic software treats as missing
    undetermined_chars = set('-NnXx?.')

    # Get alignment length
    seq_length = len(next(iter(sequences.values())))

    if verbose:
        print(f"\nChecking for all-gap/undetermined columns in {seq_length:,} positions...", file=sys.stderr)

    # Find positions to keep (at least one determined character)
    positions_to_keep = []

    for pos in range(seq_length):
        has_determined = False
        for seq in sequences.values():
            if pos < len(seq) and seq[pos] not in undetermined_chars:
                has_determined = True
                break

        if has_determined:
            positions_to_keep.append(pos)

    # Build new sequences with only kept positions
    cleaned_sequences = {}
    for seq_id, seq in sequences.items():
        cleaned_seq = ''.join(seq[pos] for pos in positions_to_keep)
        cleaned_sequences[seq_id] = cleaned_seq

    if verbose:
        removed = seq_length - len(positions_to_keep)
        print(f"  Removed {removed:,} all-gap/undetermined columns", file=sys.stderr)
        print(f"  Final MSA length: {len(positions_to_keep):,} positions", file=sys.stderr)

    return cleaned_sequences, positions_to_keep


def calculate_coverage_stats(depth_per_position, msa_sequences, sample_name="Sample"):
    """
    Calculate coverage statistics from depth information.

    Args:
        depth_per_position: Dict mapping MSA position to read depth
        msa_sequences: Dict of MSA sequences (to get reference length)
        sample_name: Name of sample for output

    Returns:
        dict: Coverage statistics
    """
    if not depth_per_position:
        return {
            'sample': sample_name,
            'covered_positions': 0,
            'total_reference_length': 0,
            'coverage_percent': 0.0,
            'mean_depth': 0.0,
            'median_depth': 0.0,
            'max_depth': 0,
            'min_depth': 0,
            'depth_histogram': {}
        }

    # Get total ungapped reference length (union of all references)
    total_ref_length = 0
    for seq in msa_sequences.values():
        total_ref_length += len(seq.replace('-', ''))

    # Calculate depth statistics
    depths = list(depth_per_position.values())
    covered_positions = len(depths)
    mean_depth = sum(depths) / len(depths) if depths else 0
    sorted_depths = sorted(depths)
    median_depth = sorted_depths[len(sorted_depths) // 2] if sorted_depths else 0
    max_depth = max(depths) if depths else 0
    min_depth = min(depths) if depths else 0

    # Coverage as percentage of MSA positions (not ungapped ref length)
    # Use the number of unique MSA positions covered
    msa_length = len(next(iter(msa_sequences.values()))) if msa_sequences else 0
    coverage_percent = (covered_positions / msa_length * 100) if msa_length > 0 else 0

    # Depth histogram (binned)
    depth_histogram = defaultdict(int)
    for d in depths:
        depth_histogram[d] += 1

    return {
        'sample': sample_name,
        'covered_positions': covered_positions,
        'msa_length': msa_length,
        'total_reference_length': total_ref_length,
        'coverage_percent': coverage_percent,
        'mean_depth': mean_depth,
        'median_depth': median_depth,
        'max_depth': max_depth,
        'min_depth': min_depth,
        'depth_histogram': dict(depth_histogram)
    }


def print_coverage_stats(stats_list, combined_stats=None, file=sys.stderr):
    """
    Print coverage statistics in a formatted table.

    Args:
        stats_list: List of per-sample statistics dicts
        combined_stats: Combined statistics dict (optional)
        file: Output file handle
    """
    print(f"\n{'='*70}", file=file)
    print(f"COVERAGE AND DEPTH STATISTICS", file=file)
    print(f"{'='*70}", file=file)

    # Header
    print(f"\n{'Sample':<30} {'Reads':>10} {'Covered':>12} {'Coverage':>10} {'Mean':>8} {'Median':>8} {'Max':>8}",
          file=file)
    print(f"{'':<30} {'':>10} {'Positions':>12} {'(%)':>10} {'Depth':>8} {'Depth':>8} {'Depth':>8}",
          file=file)
    print(f"{'-'*30} {'-'*10} {'-'*12} {'-'*10} {'-'*8} {'-'*8} {'-'*8}", file=file)

    for stats in stats_list:
        sample = stats.get('sample', 'Unknown')[:30]
        reads = stats.get('mapped_reads', 0)
        covered = stats.get('covered_positions', 0)
        coverage = stats.get('coverage_percent', 0)
        mean_d = stats.get('mean_depth', 0)
        median_d = stats.get('median_depth', 0)
        max_d = stats.get('max_depth', 0)

        print(f"{sample:<30} {reads:>10,} {covered:>12,} {coverage:>10.2f} {mean_d:>8.1f} {median_d:>8.1f} {max_d:>8,}",
              file=file)

    if combined_stats:
        print(f"{'-'*30} {'-'*10} {'-'*12} {'-'*10} {'-'*8} {'-'*8} {'-'*8}", file=file)
        sample = "COMBINED"
        reads = combined_stats.get('mapped_reads', 0)
        covered = combined_stats.get('covered_positions', 0)
        coverage = combined_stats.get('coverage_percent', 0)
        mean_d = combined_stats.get('mean_depth', 0)
        median_d = combined_stats.get('median_depth', 0)
        max_d = combined_stats.get('max_depth', 0)

        print(f"{sample:<30} {reads:>10,} {covered:>12,} {coverage:>10.2f} {mean_d:>8.1f} {median_d:>8.1f} {max_d:>8,}",
              file=file)

    print(f"{'='*70}\n", file=file)


def write_depth_file(depth_per_position, output_file, msa_length=None):
    """
    Write per-position depth to a file.

    Args:
        depth_per_position: Dict mapping MSA position to read depth
        output_file: Path to output file
        msa_length: Total MSA length (to include zero-depth positions)
    """
    with open(output_file, 'w') as f:
        f.write("position\tdepth\n")

        if msa_length:
            # Output all positions including zero depth
            for pos in range(msa_length):
                depth = depth_per_position.get(pos, 0)
                f.write(f"{pos}\t{depth}\n")
        else:
            # Only output covered positions
            for pos in sorted(depth_per_position.keys()):
                f.write(f"{pos}\t{depth_per_position[pos]}\n")


def write_fasta(sequences, output_file, line_width=80):
    """
    Write sequences to FASTA file.

    Args:
        sequences: Dictionary mapping sequence IDs to sequences
        output_file: Path to output FASTA file
        line_width: Maximum line width for sequence lines (0 = no wrapping)
    """
    with open(output_file, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f'>{seq_id}\n')

            if line_width > 0:
                # Write sequence with line breaks
                for i in range(0, len(seq), line_width):
                    f.write(seq[i:i+line_width] + '\n')
            else:
                # Write as single line
                f.write(seq + '\n')


def iter_alignments_with_tags(path):
    """
    Like iter_alignments but also yields md_tag (str or None).

    Yields: (read_name, flag, ref_name, ref_pos_0based, cigar, seq, qual, md_tag)
    """
    if path.endswith('.bam'):
        if not _PYSAM_AVAILABLE:
            print("Error: pysam is required to read BAM files. Install with: pip install pysam",
                  file=sys.stderr)
            sys.exit(1)
        with pysam.AlignmentFile(path, 'rb') as af:
            for read in af.fetch(until_eof=True):
                cigar = read.cigarstring if read.cigarstring else '*'
                seq = read.query_sequence if read.query_sequence else ''
                if read.query_qualities is not None:
                    qual = ''.join(chr(q + 33) for q in read.query_qualities)
                else:
                    qual = None
                md_tag = read.get_tag('MD') if read.has_tag('MD') else None
                yield (read.query_name, read.flag,
                       read.reference_name if read.reference_name else '*',
                       read.reference_start,
                       cigar, seq, qual, md_tag)
    else:
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue
                qual = fields[10] if fields[10] != '*' else None
                md_tag = None
                for field in fields[11:]:
                    if field.startswith('MD:Z:'):
                        md_tag = field[5:]
                        break
                yield (fields[0], int(fields[1]), fields[2],
                       int(fields[3]) - 1,
                       fields[5], fields[9], qual, md_tag)


def count_md_mismatches(md_tag):
    """Count substitution mismatches from MD tag string (deletions not counted)."""
    md_no_del = re.sub(r'\^[ACGT]+', '', md_tag, flags=re.IGNORECASE)
    return len(re.findall(r'[ACGT]', md_no_del, re.IGNORECASE))


def build_sam_record_lookup(sam_file):
    """
    Pre-scan SAM and return {read_name: (flag, ref_name, ref_pos_0based, cigar, seq)}.
    Skips unmapped reads. On duplicate names, keeps first occurrence.
    """
    records = {}
    for read_name, flag, ref_name, ref_pos, cigar, seq, qual in iter_alignments(sam_file):
        if flag & 4:
            continue
        if read_name not in records:
            records[read_name] = (flag, ref_name, ref_pos, cigar, seq)
    return records


def count_fastq_msa_mismatches(fastq_seq, ref_span):
    """
    Count mismatches between a FASTQ sequence and its filtered-MSA reference span.
    Both strings are in filtered-MSA coordinates (same length). Gaps and Ns are skipped.
    """
    count = 0
    for rb, fb in zip(fastq_seq, ref_span):
        if rb == '-' or fb == '-':
            continue
        if rb.upper() == 'N' or fb.upper() == 'N':
            continue
        if rb.upper() != fb.upper():
            count += 1
    return count


def count_sam_path_mismatches(flag, ref_name, ref_pos, cigar, seq,
                               msa_sequences, pos_map, trim_5, trim_3):
    """
    Recompute mismatch count via the SAM path with trimming applied.
    Mirrors the downstream SAM analysis logic.
    """
    alignment = parse_cigar_alignment(cigar, ref_pos, seq, None)
    if not alignment:
        return 0

    if trim_5 + trim_3 >= len(alignment):
        return 0

    if trim_3 > 0:
        trimmed = alignment[trim_5:-trim_3]
    else:
        trimmed = alignment[trim_5:]

    if not trimmed:
        return 0

    msa_seq = msa_sequences[ref_name]
    count = 0
    for rp, base, _ in trimmed:
        if base == '-':
            continue
        if rp >= len(pos_map) - 1:
            continue
        msa_pos = pos_map[rp]
        if msa_pos >= len(msa_seq):
            continue
        msa_base = msa_seq[msa_pos].upper()
        if base.upper() != msa_base and base.upper() != 'N' and msa_base not in ('N', '-'):
            count += 1
    return count


def validate_fastq_against_msa(sam_file, fastq_file, filtered_msa_file, msa_sequences,
                                output_file, verbose, trim_5, trim_3, mode='w'):
    """
    Compare per-read mismatch counts from the SAM path vs the extracted FASTQ path.

    SAM path: parse CIGAR, apply trim_5/trim_3, compare trimmed bases vs original MSA
              via pos_map (same logic as downstream SAM analysis).
    FASTQ path: read extracted FASTQ, compare bases vs filtered MSA at [start:start+len)
                (same logic as downstream FASTQ analysis).

    Discordant reads (where the two counts differ) are written to output_file.

    Args:
        sam_file: Path to SAM/BAM file
        fastq_file: Path to already-extracted FASTQ file
        filtered_msa_file: Path to the pipeline's output filtered MSA FASTA
        msa_sequences: Dict of original MSA sequences
        output_file: Path to output file for discordant reads
        verbose: Print progress messages
        trim_5: Bases trimmed from 5' end during extraction
        trim_3: Bases trimmed from 3' end during extraction
        mode: File open mode ('w' or 'a')
    """
    sam_records = build_sam_record_lookup(sam_file)
    filtered_msa = parse_fasta(filtered_msa_file)

    position_maps = {}
    total_checked = 0
    total_skipped = 0
    total_discordant = 0
    BASES = 'ACGT'
    sub_tally = defaultdict(lambda: defaultdict(int))  # sub_tally[ref_base][obs_base]

    if verbose:
        print(f"\nValidating FASTQ against MSA: {fastq_file}", file=sys.stderr)

    with open(output_file, mode) as out_f, open(fastq_file, 'r') as fq_f:
        while True:
            header = fq_f.readline()
            if not header:
                break
            fastq_seq = fq_f.readline().rstrip('\n')
            fq_f.readline()  # '+'
            fq_f.readline()  # quality

            header = header.rstrip('\n')
            if not header.startswith('@'):
                continue

            # Parse header: @read_name start=X len=Y ref=Z rc=W
            parts = header[1:].split()
            read_name = parts[0]

            start = None
            gapped_len = None
            ref_name = None
            rc = None
            for part in parts[1:]:
                if part.startswith('start='):
                    start = int(part[6:])
                elif part.startswith('len='):
                    gapped_len = int(part[4:])
                elif part.startswith('ref='):
                    ref_name = part[4:]
                elif part.startswith('rc='):
                    rc = part[3:]

            if start is None or gapped_len is None or ref_name is None:
                total_skipped += 1
                continue

            if read_name not in sam_records:
                total_skipped += 1
                continue

            if ref_name not in filtered_msa:
                total_skipped += 1
                continue

            ref_span = filtered_msa[ref_name][start: start + gapped_len]

            # Accumulate substitution tally (FASTQ path, no masking)
            for rb, fb in zip(fastq_seq, ref_span):
                if rb == '-' or fb == '-':
                    continue
                rb_u, fb_u = rb.upper(), fb.upper()
                if rb_u not in BASES or fb_u not in BASES:
                    continue
                sub_tally[fb_u][rb_u] += 1

            fastq_mm = count_fastq_msa_mismatches(fastq_seq, ref_span)

            flag, ref_name_sam, ref_pos, cigar, seq = sam_records[read_name]

            if ref_name_sam not in msa_sequences:
                total_skipped += 1
                continue

            if ref_name_sam not in position_maps:
                position_maps[ref_name_sam] = build_position_map(msa_sequences[ref_name_sam], False)

            pos_map = position_maps[ref_name_sam]
            sam_mm = count_sam_path_mismatches(
                flag, ref_name_sam, ref_pos, cigar, seq,
                msa_sequences, pos_map, trim_5, trim_3
            )

            total_checked += 1

            if fastq_mm != sam_mm:
                total_discordant += 1

                # Collect SAM-path mismatch positions for display
                alignment = parse_cigar_alignment(cigar, ref_pos, seq, None)
                if trim_3 > 0:
                    trimmed = alignment[trim_5:-trim_3]
                else:
                    trimmed = alignment[trim_5:]

                msa_seq = msa_sequences[ref_name_sam]
                sam_mismatch_details = []
                for rp, base, _ in trimmed:
                    if base == '-':
                        continue
                    if rp >= len(pos_map) - 1:
                        continue
                    msa_pos = pos_map[rp]
                    if msa_pos >= len(msa_seq):
                        continue
                    msa_base = msa_seq[msa_pos].upper()
                    if base.upper() != msa_base and base.upper() != 'N' and msa_base not in ('N', '-'):
                        sam_mismatch_details.append((rp, base, msa_base))

                # Build FASTQ-path match string
                match_chars = []
                for rb, fb in zip(fastq_seq, ref_span):
                    if rb == '-' or fb == '-':
                        match_chars.append(' ')
                    elif rb.upper() == 'N' or fb.upper() == 'N':
                        match_chars.append(' ')
                    elif rb.upper() == fb.upper():
                        match_chars.append('|')
                    else:
                        match_chars.append('X')
                match_str = ''.join(match_chars)

                out_f.write(f"READ: {read_name}  ref={ref_name}  start={start}  len={gapped_len}"
                            f"  rc={rc}  SAM_mm={sam_mm}  FASTQ_mm={fastq_mm}\n")
                out_f.write(f"  SAM path (trimmed alignment vs original MSA):\n")
                for rp, base, msa_base in sam_mismatch_details:
                    out_f.write(f"    ref_pos={rp}  read={base}  msa={msa_base}  X\n")
                out_f.write(f"  FASTQ path (extracted read vs filtered MSA):\n")
                out_f.write(f"    read: {fastq_seq}\n")
                out_f.write(f"          {match_str}\n")
                out_f.write(f"    ref:  {ref_span}\n")
                out_f.write(f"\n")

    # Print substitution matrix
    total_bases = sum(sub_tally[rb][ob] for rb in BASES for ob in BASES)
    total_mm    = sum(sub_tally[rb][ob] for rb in BASES for ob in BASES if rb != ob)

    print(f"\nSubstitution matrix (FASTQ path vs filtered MSA, no masking):", file=sys.stderr)
    print(f"  {'':3} {'A':>10} {'C':>10} {'G':>10} {'T':>10}  {'Total ref':>10}",
          file=sys.stderr)
    for ref in BASES:
        row_total = sum(sub_tally[ref][ob] for ob in BASES)
        cells = ''.join(
            f"{'—':>10} " if ref == ob else f"{sub_tally[ref][ob]:>10,} "
            for ob in BASES
        )
        print(f"  {ref}: {cells} {row_total:>10,}", file=sys.stderr)
    print(f"  Total mismatches : {total_mm:,}", file=sys.stderr)
    print(f"  Total bases      : {total_bases:,}", file=sys.stderr)
    print(f"  {fastq_file}: {total_checked:,} reads checked, {total_skipped:,} skipped "
          f"(no SAM match), {total_discordant:,} discordant", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Filter MSA to positions where SAM/BAM alignments occur',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Filter MSA to positions with alignments from single SAM file
  %(prog)s -s alignments.sam -m msa.fasta -o filtered_msa.fa -v

  # Filter MSA using multiple SAM files (positions combined)
  %(prog)s -s sample1.sam sample2.sam -m msa.fasta -o filtered_msa.fa -v

  # Extract reads to separate FASTQ files for each SAM input
  %(prog)s -s sample1.sam sample2.sam -m msa.fasta -o filtered_msa.fa --extract-reads sample1.fq sample2.fq -v

  # Trim 5bp from both ends of reads (useful for aDNA with deamination damage)
  %(prog)s -s alignments.sam -m msa.fasta -o filtered_msa.fa --trim 5 -v

  # Calculate coverage and depth statistics
  %(prog)s -s sample1.sam sample2.sam -m msa.fasta -o filtered_msa.fa --coverage-stats -v

  # Output per-position depth to file
  %(prog)s -s alignments.sam -m msa.fasta -o filtered_msa.fa --depth-file depth.tsv

Notes:
  - Multiple SAM files can be provided; MSA positions are combined from all
  - When using --extract-reads with multiple SAM files, provide matching output files
  - The SAM files should contain alignments to the same ungapped reference sequences
  - The MSA file should contain the same sequences but with gaps (-)
  - Trimming removes bases from read ends before extracting MSA positions
  - Coverage stats show per-sample and combined read depth and coverage percentage
        """
    )

    parser.add_argument('-s', '--sam', required=True, nargs='+',
                        help='Input SAM file(s) - multiple files can be provided')
    parser.add_argument('-m', '--msa', required=True,
                        help='Input MSA FASTA file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output filtered MSA FASTA file')
    parser.add_argument('-w', '--line-width', type=int, default=80,
                        help='Line width for output FASTA (0 = no wrapping, default: 80)')
    parser.add_argument('--trim-5', type=int, default=0,
                        help='Trim N bases from 5\' end of each read (default: 0)')
    parser.add_argument('--trim-3', type=int, default=0,
                        help='Trim N bases from 3\' end of each read (default: 0)')
    parser.add_argument('--trim', type=int,
                        help='Trim N bases from both ends of each read (shortcut for --trim-5 N --trim-3 N)')
    parser.add_argument('--extract-reads', metavar='FILE', nargs='+',
                        help='Extract reads to FASTQ file(s) - provide one output file per input SAM file')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print verbose output')
    parser.add_argument('--coverage-stats', action='store_true',
                        help='Calculate and print coverage/depth statistics')
    parser.add_argument('--depth-file', metavar='FILE',
                        help='Output per-position depth to file (TSV format)')
    parser.add_argument('--depth-all-positions', action='store_true',
                        help='Include zero-depth positions in depth file')
    parser.add_argument('--validate-mismatches', metavar='FILE',
                        help='Validate MD-tag vs MSA mismatch counts; write discordant reads to FILE')

    args = parser.parse_args()

    # Handle --trim shortcut
    trim_5 = args.trim if args.trim is not None else args.trim_5
    trim_3 = args.trim if args.trim is not None else args.trim_3

    # Validate extract-reads count matches sam count
    if args.extract_reads:
        if len(args.extract_reads) != len(args.sam):
            print(f"Error: number of --extract-reads files ({len(args.extract_reads)}) must match "
                  f"number of --sam files ({len(args.sam)})", file=sys.stderr)
            sys.exit(1)

    if args.validate_mismatches and not args.extract_reads:
        print("Error: --validate-mismatches requires --extract-reads", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print(f"{'='*60}", file=sys.stderr)
        print(f"Filter MSA by SAM Alignments", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)
        print(f"SAM file(s): {len(args.sam)} input(s)", file=sys.stderr)
        for sam_file in args.sam:
            print(f"  - {sam_file}", file=sys.stderr)
        print(f"MSA file: {args.msa}", file=sys.stderr)
        print(f"Output file: {args.output}", file=sys.stderr)
        if args.extract_reads:
            print(f"Extract reads to: {len(args.extract_reads)} output(s)", file=sys.stderr)
            for i, fq_file in enumerate(args.extract_reads):
                print(f"  - {args.sam[i]} -> {fq_file}", file=sys.stderr)
        if trim_5 > 0 or trim_3 > 0:
            print(f"Read trimming: {trim_5} bp from 5' end, {trim_3} bp from 3' end", file=sys.stderr)
        print(f"", file=sys.stderr)

    # Load MSA sequences
    if args.verbose:
        print(f"Loading MSA sequences...", file=sys.stderr)

    msa_sequences = parse_fasta(args.msa)

    if args.verbose:
        print(f"Loaded {len(msa_sequences)} sequences from MSA", file=sys.stderr)
        for seq_id in list(msa_sequences.keys())[:5]:
            print(f"  - {seq_id}", file=sys.stderr)
        if len(msa_sequences) > 5:
            print(f"  ... and {len(msa_sequences) - 5} more", file=sys.stderr)
        print(f"", file=sys.stderr)

    # Collect MSA positions from all SAM alignments (with depth tracking if requested)
    if args.verbose:
        print(f"Collecting MSA positions from SAM alignments...", file=sys.stderr)

    msa_positions = set()
    combined_depth = defaultdict(int)
    per_sample_stats = []
    total_mapped_reads = 0

    for i, sam_file in enumerate(args.sam):
        if args.verbose and len(args.sam) > 1:
            print(f"\nProcessing SAM file {i+1}/{len(args.sam)}: {sam_file}", file=sys.stderr)

        # Use depth-tracking version if coverage stats requested
        if args.coverage_stats or args.depth_file:
            sam_positions, sample_depth, sample_stats = collect_msa_positions_with_depth(
                sam_file, msa_sequences, args.verbose, trim_5, trim_3
            )

            # Calculate per-sample coverage stats
            import os
            sample_name = os.path.basename(sam_file)
            coverage_stats = calculate_coverage_stats(sample_depth, msa_sequences, sample_name)
            coverage_stats['mapped_reads'] = sample_stats['mapped_reads']
            per_sample_stats.append(coverage_stats)

            # Accumulate combined depth
            for pos, depth in sample_depth.items():
                combined_depth[pos] += depth

            total_mapped_reads += sample_stats['mapped_reads']
        else:
            sam_positions = collect_msa_positions(sam_file, msa_sequences, args.verbose, trim_5, trim_3)

        msa_positions.update(sam_positions)

    if args.verbose and len(args.sam) > 1:
        print(f"\nCombined MSA positions from all SAM files: {len(msa_positions):,}", file=sys.stderr)

    # Calculate and print coverage statistics
    if args.coverage_stats or args.depth_file:
        combined_coverage_stats = calculate_coverage_stats(combined_depth, msa_sequences, "COMBINED")
        combined_coverage_stats['mapped_reads'] = total_mapped_reads

        if args.coverage_stats:
            print_coverage_stats(per_sample_stats, combined_coverage_stats if len(args.sam) > 1 else None)

        # Write depth file if requested
        if args.depth_file:
            msa_length = len(next(iter(msa_sequences.values()))) if msa_sequences else None
            if args.depth_all_positions:
                write_depth_file(combined_depth, args.depth_file, msa_length)
            else:
                write_depth_file(combined_depth, args.depth_file)
            if args.verbose:
                print(f"Depth file written to: {args.depth_file}", file=sys.stderr)

    # Filter MSA to selected positions
    filtered_msa, sorted_msa_positions = filter_msa_by_positions(msa_sequences, msa_positions, args.verbose)

    # Remove all-gap columns
    cleaned_msa, kept_indices = remove_all_gap_columns(filtered_msa, args.verbose)

    # Compute final MSA positions (original MSA positions that made it to output)
    final_msa_positions = [sorted_msa_positions[i] for i in kept_indices]

    # Write output
    if args.verbose:
        print(f"\nWriting filtered MSA to {args.output}...", file=sys.stderr)

    write_fasta(cleaned_msa, args.output, args.line_width)

    # Extract reads if requested (one FASTQ per SAM file)
    if args.extract_reads:
        for i, (sam_file, fq_file) in enumerate(zip(args.sam, args.extract_reads)):
            if args.verbose and len(args.sam) > 1:
                print(f"\nExtracting reads from SAM {i+1}/{len(args.sam)}: {sam_file}", file=sys.stderr)
            extract_reads_to_filtered_msa(
                sam_file, msa_sequences, final_msa_positions,
                fq_file, args.verbose, trim_5, trim_3
            )
            if args.validate_mismatches:
                validate_fastq_against_msa(
                    sam_file, fq_file, args.output,
                    msa_sequences, args.validate_mismatches,
                    args.verbose, trim_5, trim_3,
                    mode='w' if i == 0 else 'a'
                )

    if args.verbose:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"Done!", file=sys.stderr)
        print(f"Filtered MSA written to: {args.output}", file=sys.stderr)
        if args.extract_reads:
            print(f"Reads extracted to:", file=sys.stderr)
            for fq_file in args.extract_reads:
                print(f"  - {fq_file}", file=sys.stderr)
        if len(msa_positions) > 0:
            first_msa = next(iter(cleaned_msa.values()))
            print(f"Output MSA length: {len(first_msa):,} positions", file=sys.stderr)


if __name__ == '__main__':
    main()
