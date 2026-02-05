#!/usr/bin/env python3
"""
Generate inflated barcode whitelist files for testing large k-mer indexer paths.

Appends random unique DNA sequences to existing barcode files to push
the total count above 100K, triggering Array2BitKmerIndexer / Dict2BitKmerIndexer
in the universal barcode calling algorithm.
"""

import os
import random
import sys

NUCLEOTIDES = "ACGT"
DATA_DIR = "/abga/work/andreyp/ci_isoquant/data/barcodes"

# Configuration: (input_file, output_file, dummy_count, barcode_length)
CONFIGS = [
    ("10xMultiome_5K.tsv", "10x.large_barcodes.tsv", 120000, 16),
    ("A0079_044_BeadBarcodes.tsv", "curio.large_barcodes.tsv", 50000, 14),
    ("visium_v1_hd.slide1.barcodes1.tsv", "visium.large_barcodes1.tsv", 120000, 15),
    ("visium_v1_hd.slide1.barcodes2.tsv", "visium.large_barcodes2.tsv", 120000, 14),
]


def load_existing_barcodes(filepath: str) -> set[str]:
    """Load first column from TSV barcode file."""
    barcodes = set()
    with open(filepath) as f:
        for line in f:
            bc = line.strip().split()[0]
            if bc:
                barcodes.add(bc)
    return barcodes


def generate_random_barcode(length: int) -> str:
    return "".join(random.choice(NUCLEOTIDES) for _ in range(length))


def generate_unique_barcodes(count: int, length: int, existing: set[str]) -> list[str]:
    """Generate random unique barcodes not in the existing set."""
    new_barcodes = []
    seen = set(existing)
    while len(new_barcodes) < count:
        bc = generate_random_barcode(length)
        if bc not in seen:
            seen.add(bc)
            new_barcodes.append(bc)
    return new_barcodes


def main():
    random.seed(42)  # Reproducible results

    for input_name, output_name, dummy_count, bc_length in CONFIGS:
        input_path = os.path.join(DATA_DIR, input_name)
        output_path = os.path.join(DATA_DIR, output_name)

        print(f"Processing {input_name} -> {output_name}")
        existing = load_existing_barcodes(input_path)
        print(f"  Loaded {len(existing)} existing barcodes")

        dummy = generate_unique_barcodes(dummy_count, bc_length, existing)
        print(f"  Generated {len(dummy)} dummy barcodes ({bc_length}bp)")

        # Write output: copy original file + append dummy barcodes
        with open(output_path, "w") as out:
            with open(input_path) as inp:
                for line in inp:
                    out.write(line)
            for bc in dummy:
                out.write(bc + "\n")

        total = len(existing) + len(dummy)
        print(f"  Total: {total} barcodes -> {output_path}")

    print("Done.")


if __name__ == "__main__":
    main()
