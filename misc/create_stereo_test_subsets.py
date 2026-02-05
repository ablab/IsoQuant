#!/usr/bin/env python3
"""
Create stereo-seq barcode subsets and corresponding filtered read files.

From the full 10M stereo barcode whitelist:
1. Randomly select 1M barcodes -> stereo.1M.tsv
2. Randomly select 80K barcodes (subset of the 1M) -> stereo.80K.tsv
3. Filter reads to keep only those matching each subset

Read ID format: >READ_<idx>_<transcript>_<BARCODE>_<UMI>_<extra>...
Barcode is field 3 (0-indexed, underscore-delimited).
"""

import os
import random
import sys

DATA_DIR = "/abga/work/andreyp/ci_isoquant/data/barcodes"

FULL_BARCODE_FILE = os.path.join(DATA_DIR, "D04620C2.10M.tsv")
FULL_READS_FILE = os.path.join(DATA_DIR, "Mouse.StereoSeq.D04620C2.10M.800K.fa")

SUBSET_1M_FILE = os.path.join(DATA_DIR, "stereo.1M.tsv")
SUBSET_80K_FILE = os.path.join(DATA_DIR, "stereo.80K.tsv")
READS_1M_FILE = os.path.join(DATA_DIR, "Mouse.StereoSeq.custom_sc.1M.fa")
READS_80K_FILE = os.path.join(DATA_DIR, "Mouse.StereoSeq.custom_sc.80K.fa")

SUBSET_1M_SIZE = 1000000
SUBSET_80K_SIZE = 80000


def load_barcodes_with_lines(filepath: str) -> tuple[list[str], list[str]]:
    """Load barcodes (first column) and full lines from TSV file."""
    barcodes = []
    lines = []
    with open(filepath) as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            bc = stripped.split()[0]
            barcodes.append(bc)
            lines.append(line)
    return barcodes, lines


def extract_barcode_from_read_id(read_id: str) -> str:
    """Extract barcode from FASTA read ID.

    Format: >READ_idx_transcript_BARCODE_UMI_...
    """
    parts = read_id.lstrip(">").split("_")
    if len(parts) >= 4:
        return parts[3]
    return ""


def filter_reads_fasta(input_fasta: str, output_fasta: str, barcode_set: set[str]) -> int:
    """Filter FASTA file, keeping reads with barcodes in the given set."""
    kept = 0
    total = 0
    with open(input_fasta) as inf, open(output_fasta, "w") as outf:
        keep_current = False
        for line in inf:
            if line.startswith(">"):
                total += 1
                bc = extract_barcode_from_read_id(line.strip())
                keep_current = bc in barcode_set
                if keep_current:
                    kept += 1
            if keep_current:
                outf.write(line)
    return kept, total


def main():
    random.seed(42)

    print("Loading full barcode set...")
    barcodes, lines = load_barcodes_with_lines(FULL_BARCODE_FILE)
    print(f"  Loaded {len(barcodes)} barcodes")

    # Select 1M subset
    print(f"Selecting {SUBSET_1M_SIZE} barcodes for 1M subset...")
    indices_1m = sorted(random.sample(range(len(barcodes)), SUBSET_1M_SIZE))
    barcodes_1m = set(barcodes[i] for i in indices_1m)

    with open(SUBSET_1M_FILE, "w") as f:
        for i in indices_1m:
            f.write(lines[i])
    print(f"  Written to {SUBSET_1M_FILE}")

    # Select 80K subset from the 1M subset
    print(f"Selecting {SUBSET_80K_SIZE} barcodes for 80K subset (from 1M)...")
    indices_80k_in_1m = sorted(random.sample(range(len(indices_1m)), SUBSET_80K_SIZE))
    barcodes_80k = set(barcodes[indices_1m[i]] for i in indices_80k_in_1m)

    with open(SUBSET_80K_FILE, "w") as f:
        for i in indices_80k_in_1m:
            orig_idx = indices_1m[i]
            f.write(lines[orig_idx])
    print(f"  Written to {SUBSET_80K_FILE}")

    # Filter reads for 1M subset
    print("Filtering reads for 1M subset...")
    kept_1m, total = filter_reads_fasta(FULL_READS_FILE, READS_1M_FILE, barcodes_1m)
    print(f"  Kept {kept_1m}/{total} reads -> {READS_1M_FILE}")

    # Filter reads for 80K subset
    print("Filtering reads for 80K subset...")
    kept_80k, total = filter_reads_fasta(FULL_READS_FILE, READS_80K_FILE, barcodes_80k)
    print(f"  Kept {kept_80k}/{total} reads -> {READS_80K_FILE}")

    print("Done.")


if __name__ == "__main__":
    main()
