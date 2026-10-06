#!/usr/bin/env python3
"""
Filter reads from a multi-reference bowtie2 SAM/BAM into per-compartment SAM files.

Filtering logic (per QNAME group):
1. Discard unmapped records (FLAG 0x4). If all unmapped → discard read.
2. Outgroup filter: any mapped RNAME in outgroup → discard read.
3. Unknown reference filter: any mapped RNAME not in compartments (and not outgroup) → discard.
4. Within-reference multimap filter: any RNAME appears more than once → discard.
5. Cross-compartment filter: mapped records span more than one compartment → discard.
6. Primary alignment check: must have exactly one primary alignment → discard if not.
7. KEEP: write the single primary alignment record to <compartment>.sam.
"""

import argparse
import collections
import csv
import os
import tempfile
from pathlib import Path

import pysam


def parse_compartments(path: str) -> dict[str, str]:
    """Parse compartments CSV → dict[reference, compartment]."""
    ref_to_comp = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ref = row["reference"].strip()
            comp = row["compartment"].strip()
            if ref:
                ref_to_comp[ref] = comp
    return ref_to_comp


def parse_outgroup(path: str) -> set[str]:
    """Parse outgroup TXT → set of reference names."""
    outgroup = set()
    with open(path) as f:
        for line in f:
            ref = line.strip()
            if ref:
                outgroup.add(ref)
    return outgroup


def build_compartment_headers(in_header: pysam.AlignmentHeader, ref_to_comp: dict[str, str]) -> dict[str, dict]:
    """
    Build per-compartment header dicts with filtered @SQ lines.
    Returns dict[compartment, header_dict].
    """
    header_dict = in_header.to_dict()

    # Group @SQ entries by compartment
    comp_sq: dict[str, list] = collections.defaultdict(list)
    for sq in header_dict.get("SQ", []):
        sn = sq.get("SN", "")
        if sn in ref_to_comp:
            comp_sq[ref_to_comp[sn]].append(sq)

    comp_headers = {}
    for comp, sq_list in comp_sq.items():
        h = {}
        if "HD" in header_dict:
            h["HD"] = header_dict["HD"]
        h["SQ"] = sq_list
        if "PG" in header_dict:
            h["PG"] = header_dict["PG"]
        if "CO" in header_dict:
            h["CO"] = header_dict["CO"]
        comp_headers[comp] = h

    return comp_headers


def classify_read_group(
    records: list,
    ref_to_comp: dict[str, str],
    outgroup: set[str],
    min_score: int = 24,
) -> tuple[str | None, object | None]:
    """
    Apply all 6 filters to a QNAME group.
    Returns (compartment, primary_record) or (None, None) if discarded.
    """
    # Step 1: collect mapped records
    mapped = [r for r in records if not r.flag & 0x4]
    if not mapped:
        return None, None

    # Step 2: outgroup filter
    for r in mapped:
        if r.reference_name in outgroup:
            return None, None

    # Step 3: unknown reference filter
    for r in mapped:
        if r.reference_name not in ref_to_comp:
            return None, None

    # Step 4: within-reference multimap filter
    rname_counts = collections.Counter(r.reference_name for r in mapped)
    if any(c > 1 for c in rname_counts.values()):
        return None, None

    # Step 5: cross-compartment filter
    compartments = {ref_to_comp[r.reference_name] for r in mapped}
    if len(compartments) > 1:
        return None, None
    compartment = compartments.pop()

    # Step 6: primary alignment check
    # Primary = not secondary (0x100) and not supplementary (0x800)
    primaries = [r for r in mapped if not (r.flag & 0x100) and not (r.flag & 0x800)]
    if len(primaries) != 1:
        return None, None
    primary = primaries[0]

    # Step 7: alignment score filter on the primary alignment
    try:
        as_val = primary.get_tag("AS")
        if as_val < min_score:
            return None, None
    except KeyError:
        return None, None

    return compartment, primary


def main():
    parser = argparse.ArgumentParser(
        description="Filter SAM/BAM reads into per-compartment SAM files."
    )
    parser.add_argument("--input", required=True, help="Input SAM or BAM file")
    parser.add_argument("--compartments", required=True, help="Compartments CSV (reference,compartment)")
    parser.add_argument("--outgroup", required=True, help="Outgroup reference names (one per line)")
    parser.add_argument("--outdir", required=True, help="Output directory (created if absent)")
    parser.add_argument("--prefix", default="", help="Optional prefix for output file names")
    parser.add_argument("--min-score", type=int, default=-24, dest="min_score",
                        help="Minimum AS tag value for the primary alignment (default: -24)")
    args = parser.parse_args()

    ref_to_comp = parse_compartments(args.compartments)
    outgroup = parse_outgroup(args.outgroup)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    temp_path = None
    try:
        # Open input file
        in_bam = pysam.AlignmentFile(args.input, "r", check_sq=False)
        in_header = in_bam.header

        # Check if already sorted by queryname
        hd = in_header.to_dict().get("HD", {})
        so = hd.get("SO", "")
        need_sort = so != "queryname"

        if need_sort:
            print(f"Input is not queryname-sorted (SO={so!r}). Sorting by name to temp file...")
            tmp_fd, temp_path = tempfile.mkstemp(suffix=".bam")
            os.close(tmp_fd)
            in_bam.close()
            pysam.sort("-n", "-o", temp_path, args.input)
            in_bam = pysam.AlignmentFile(temp_path, "rb", check_sq=False)
            in_header = in_bam.header
            print("Sort complete.")

        # Build per-compartment headers
        comp_headers = build_compartment_headers(in_header, ref_to_comp)

        # Open output SAM files and build ref-name → ref-id maps for each compartment
        prefix = args.prefix
        comp_writers: dict[str, pysam.AlignmentFile] = {}
        comp_ref_index: dict[str, dict[str, int]] = {}  # compartment → {refname: idx}
        for comp, hdr in comp_headers.items():
            fname = outdir / f"{prefix}{comp}.sam"
            pysam_hdr = pysam.AlignmentHeader.from_dict(hdr)
            comp_writers[comp] = pysam.AlignmentFile(str(fname), "w", header=pysam_hdr)
            comp_ref_index[comp] = {sq["SN"]: i for i, sq in enumerate(hdr.get("SQ", []))}
            print(f"Opened output: {fname}")

        # AS score accumulators
        as_scores: dict[str, list[int]] = collections.defaultdict(list)

        # Streaming pass: accumulate per-QNAME group
        current_qname = None
        group: list = []

        stats = collections.Counter()

        def process_group(grp):
            comp, primary = classify_read_group(grp, ref_to_comp, outgroup, args.min_score)
            if comp is None:
                stats["discarded"] += 1
                return
            if comp not in comp_writers:
                # Compartment exists in data but not in header (shouldn't happen)
                stats["discarded"] += 1
                return
            # Remap reference_id to match the compartment's header
            ref_name = primary.reference_name
            primary.reference_id = comp_ref_index[comp][ref_name]
            comp_writers[comp].write(primary)
            stats["kept"] += 1
            stats[f"kept_{comp}"] += 1
            # AS tag is guaranteed present (checked in classify_read_group)
            as_scores[comp].append(int(primary.get_tag("AS")))

        for record in in_bam:
            qname = record.query_name
            if qname != current_qname:
                if group:
                    process_group(group)
                    group = []
                current_qname = qname
            group.append(record)

        if group:
            process_group(group)

        in_bam.close()

        # Close output writers
        for w in comp_writers.values():
            w.close()

        # Write AS distribution TSVs
        for comp, scores in as_scores.items():
            counter = collections.Counter(scores)
            tsv_path = outdir / f"{prefix}{comp}_AS_distribution.tsv"
            with open(tsv_path, "w") as f:
                f.write("AS_score\tcount\n")
                for score in sorted(counter):
                    f.write(f"{score}\t{counter[score]}\n")
            print(f"Wrote AS distribution: {tsv_path}")

        # Print summary
        print("\n=== Summary ===")
        print(f"Total groups processed: {stats['kept'] + stats['discarded']}")
        print(f"  Kept:      {stats['kept']}")
        print(f"  Discarded: {stats['discarded']}")
        for key, val in sorted(stats.items()):
            if key.startswith("kept_"):
                comp = key[len("kept_"):]
                print(f"    {comp}: {val}")

    finally:
        if temp_path and os.path.exists(temp_path):
            os.remove(temp_path)
            print(f"Removed temp file: {temp_path}")


if __name__ == "__main__":
    main()
