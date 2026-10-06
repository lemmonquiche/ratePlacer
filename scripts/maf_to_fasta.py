#!/usr/bin/env python3
"""
Convert MAF (Multiple Alignment Format) to FASTA format.
The reference sequence is assumed to be the first in each alignment block.
"""

import sys
import gzip
from collections import OrderedDict
from pathlib import Path


def get_genome_name(seq_name):
    """Extract genome name from sequence name (everything before the first dot)."""
    return seq_name.split('.')[0]


def is_gzipped(filename):
    """Check if file is gzipped by looking at magic bytes."""
    with open(filename, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def parse_maf(maf_file):
    """
    Parse a MAF file and extract sequences, preserving MSA format.

    All sequences will be the same length with gaps added for missing sequences.

    Args:
        maf_file: Path to MAF file (can be gzipped)

    Returns:
        OrderedDict mapping sequence names to concatenated aligned sequences
    """
    # First pass: identify all genome names and build seq→genome mapping
    seq_to_genome = {}  # e.g., "Pop_deltoides.CM125684.1" -> "Pop_deltoides"
    all_genome_names = OrderedDict()

    if is_gzipped(maf_file):
        f = gzip.open(maf_file, 'rt')
    else:
        f = open(maf_file, 'r')

    try:
        for line in f:
            if line.startswith('s'):
                parts = line.split()
                if len(parts) >= 7:
                    seq_name = parts[1]
                    genome_name = get_genome_name(seq_name)
                    seq_to_genome[seq_name] = genome_name
                    if genome_name not in all_genome_names:
                        all_genome_names[genome_name] = None
    finally:
        f.close()

    # Initialize sequences dict keyed by genome name
    sequences = OrderedDict((name, []) for name in all_genome_names)

    # Second pass: process alignment blocks
    if is_gzipped(maf_file):
        f = gzip.open(maf_file, 'rt')
    else:
        f = open(maf_file, 'r')

    try:
        current_block = {}
        block_length = 0

        for line in f:
            line = line.strip()

            # Skip comments and empty lines
            if line.startswith('#') or not line:
                continue

            # Start of new alignment block
            if line.startswith('a'):
                # Process previous block if it exists
                if current_block:
                    # Add sequences for this block (or gaps if missing)
                    for seq_name in sequences:
                        if seq_name in current_block:
                            sequences[seq_name].append(current_block[seq_name])
                        else:
                            # Add gaps for missing sequences
                            sequences[seq_name].append('-' * block_length)

                current_block = {}
                block_length = 0
                continue

            # Sequence line
            if line.startswith('s'):
                parts = line.split()
                if len(parts) >= 7:
                    # Format: s src start size strand srcSize sequence
                    seq_name = parts[1]
                    sequence = parts[6]
                    genome_name = seq_to_genome[seq_name]
                    block_length = len(sequence)  # All seqs in block should be same length
                    if genome_name in current_block:
                        print(f"  Warning: genome '{genome_name}' has multiple sequences in one block "
                              f"('{seq_name}' vs earlier entry). Keeping first.", file=sys.stderr)
                    else:
                        current_block[genome_name] = sequence

        # Don't forget the last block
        if current_block:
            for seq_name in sequences:
                if seq_name in current_block:
                    sequences[seq_name].append(current_block[seq_name])
                else:
                    sequences[seq_name].append('-' * block_length)

    finally:
        f.close()

    # Concatenate all sequences for each organism
    concatenated = OrderedDict()
    for seq_name, seq_list in sequences.items():
        concatenated[seq_name] = ''.join(seq_list)

    return concatenated


def is_valid_base(base):
    """
    Check if a base is a valid, determined nucleotide.
    Returns True for A, C, G, T (and lowercase), False for gaps, N, ?, X, etc.
    """
    return base.upper() in {'A', 'C', 'G', 'T'}


def filter_and_write_fasta(sequences, output_file):
    """
    Remove columns with only undetermined/missing data and write to FASTA.
    Filters out columns where all sequences have gaps, N's, or other ambiguous codes.
    Optimized for speed by processing chunks at a time.

    Args:
        sequences: OrderedDict mapping sequence names to sequences
        output_file: Path to output FASTA file
    """
    if not sequences:
        return

    seq_names = list(sequences.keys())
    seq_length = len(sequences[seq_names[0]])

    print("Filtering columns with only undetermined/missing data...")
    print("  (Removing columns where all sequences have gaps, N's, or ambiguous bases)")

    # Open output file handles for each sequence
    import tempfile
    temp_files = OrderedDict()
    for seq_name in seq_names:
        temp_files[seq_name] = tempfile.NamedTemporaryFile(mode='w', delete=False, buffering=8192*1024)

    try:
        # Process alignment in chunks - check which positions to keep
        chunk_size = 100000  # Process 100k positions at a time
        total_removed = 0
        total_kept = 0

        for chunk_start in range(0, seq_length, chunk_size):
            chunk_end = min(chunk_start + chunk_size, seq_length)
            chunk_len = chunk_end - chunk_start

            # Extract chunk for each sequence
            chunks = [sequences[seq_name][chunk_start:chunk_end] for seq_name in seq_names]

            # Zip sequences together to get columns, filter out uninformative columns
            # Keep only columns where at least one sequence has a determined nucleotide
            filtered_columns = []
            for column in zip(*chunks):
                # Keep column if any sequence has a valid determined base (A, C, G, T)
                if any(is_valid_base(base) for base in column):
                    filtered_columns.append(column)

            kept_count = len(filtered_columns)
            removed_count = chunk_len - kept_count

            # Transpose back to get sequences and write
            if kept_count > 0:
                # Transpose: convert list of columns back to list of sequences
                filtered_seqs = [''.join(bases) for bases in zip(*filtered_columns)]
                for seq_name, filtered_seq in zip(seq_names, filtered_seqs):
                    temp_files[seq_name].write(filtered_seq)

            total_kept += kept_count
            total_removed += removed_count

            # Print progress every 10M positions
            if chunk_end % 10000000 == 0 or chunk_end == seq_length:
                print(f"  Processed {chunk_end:,} / {seq_length:,} positions ({total_removed:,} uninformative columns removed, {total_kept:,} kept)")

        print(f"  Total uninformative columns removed: {total_removed:,}")
        print(f"  Final alignment length: {total_kept:,} bp")

        # Close temp files
        temp_file_paths = {}
        for seq_name, f in temp_files.items():
            temp_file_paths[seq_name] = f.name
            f.close()

        # Write final FASTA file
        print(f"Writing FASTA file: {output_file}")
        with open(output_file, 'w') as out:
            for seq_name in seq_names:
                out.write(f'>{seq_name}\n')

                # Read temp file and write in 60-char lines
                with open(temp_file_paths[seq_name], 'r') as temp:
                    sequence = temp.read()
                    for i in range(0, len(sequence), 60):
                        out.write(sequence[i:i+60] + '\n')

        # Clean up temp files
        import os
        for temp_path in temp_file_paths.values():
            os.unlink(temp_path)

    except Exception as e:
        # Clean up temp files on error
        import os
        for f in temp_files.values():
            try:
                f.close()
                os.unlink(f.name)
            except:
                pass
        raise e


def write_fasta(sequences, output_file):
    """
    Write sequences to FASTA format.

    Args:
        sequences: OrderedDict mapping sequence names to sequences
        output_file: Path to output FASTA file
    """
    with open(output_file, 'w') as f:
        for seq_name, sequence in sequences.items():
            f.write(f'>{seq_name}\n')
            # Write sequence in lines of 60 characters
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')


def main():
    if len(sys.argv) < 2:
        print("Usage: python maf_to_fasta.py <input.maf> [output.fasta]")
        print("  input.maf can be gzipped (.maf.gz)")
        print("  If output file not specified, uses input name with .fasta extension")
        sys.exit(1)

    input_file = sys.argv[1]

    # Determine output filename
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        # Remove .gz if present, then replace .maf with .fasta
        base = Path(input_file).stem
        if base.endswith('.maf'):
            base = base[:-4]
        output_file = base + '.fasta'

    print(f"Reading MAF file: {input_file}")
    sequences = parse_maf(input_file)

    print(f"Found {len(sequences)} genomes:")
    for genome_name, seq in sequences.items():
        print(f"  {genome_name}: {len(seq):,} bp (before filtering)")

    # Filter all-gap columns and write output
    filter_and_write_fasta(sequences, output_file)
    print("Done!")


if __name__ == '__main__':
    main()
