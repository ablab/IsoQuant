
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

#!/usr/bin/env python3

"""Generate barcode-to-spot-ID mapping for Visium HD composite barcodes.

Produces a TSV mapping composite barcodes (BC1|BC2) to spot IDs at multiple
resolutions by computing the cartesian product of two per-part coordinate files.

Usage:
    python prepare_visium_spot_ids.py part1_to_y.tsv part2_to_x.tsv -o barcode2spot.tsv
"""

import argparse
import csv
import sys
from itertools import product
from typing import Optional


def load_coordinate_file(filepath: str) -> list[tuple[str, list[str]]]:
    """Load TSV mapping barcode sequences to coordinate values.

    Args:
        filepath: Path to TSV file. Column 0 is barcode sequence,
                  columns 1..N are coordinate values at different resolutions.

    Returns:
        List of (barcode, [coord_val_1, coord_val_2, ...]) tuples.
    """
    entries: list[tuple[str, list[str]]] = []
    with open(filepath, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                sys.exit(f"Error: line has fewer than 2 columns in {filepath}: {row}")
            entries.append((row[0], row[1:]))
    return entries


def generate_spot_ids(
    file1: str,
    file2: str,
    output: str,
    prefixes: list[str],
    suffixes: list[str],
    delimiter: str,
    barcode_delimiter: str,
) -> None:
    """Generate cartesian product of barcode parts and write spot ID mapping.

    Args:
        file1: Path to TSV mapping barcode part1 to Y coordinates.
        file2: Path to TSV mapping barcode part2 to X coordinates.
        output: Path to output TSV file ("-" for stdout).
        prefixes: List of prefix strings, one per spot column.
        suffixes: List of suffix strings, one per spot column.
        delimiter: Delimiter between Y and X coordinate values within a spot ID.
        barcode_delimiter: Delimiter joining barcode part1 and part2.
    """
    entries1 = load_coordinate_file(file1)
    entries2 = load_coordinate_file(file2)

    if not entries1:
        sys.exit(f"Error: no entries loaded from {file1}")
    if not entries2:
        sys.exit(f"Error: no entries loaded from {file2}")

    n_cols = len(entries1[0][1])
    n_cols2 = len(entries2[0][1])
    if n_cols != n_cols2:
        sys.exit(
            f"Error: number of coordinate columns differs between files: "
            f"{n_cols} in {file1} vs {n_cols2} in {file2}"
        )

    if len(prefixes) == 1 and n_cols > 1:
        prefixes = prefixes * n_cols
    if len(suffixes) == 1 and n_cols > 1:
        suffixes = suffixes * n_cols

    if len(prefixes) != n_cols:
        sys.exit(
            f"Error: number of prefixes ({len(prefixes)}) does not match "
            f"number of coordinate columns ({n_cols})"
        )
    if len(suffixes) != n_cols:
        sys.exit(
            f"Error: number of suffixes ({len(suffixes)}) does not match "
            f"number of coordinate columns ({n_cols})"
        )

    out_handle = sys.stdout if output == "-" else open(output, "w", newline="")
    try:
        writer = csv.writer(out_handle, delimiter="\t")
        for (bc1, coords1), (bc2, coords2) in product(entries1, entries2):
            composite_barcode = f"{bc1}{barcode_delimiter}{bc2}"
            spot_ids = [
                f"{prefixes[i]}{coords1[i]}{delimiter}{coords2[i]}{suffixes[i]}"
                for i in range(n_cols)
            ]
            writer.writerow([composite_barcode] + spot_ids)
    finally:
        if out_handle is not sys.stdout:
            out_handle.close()


def main(args: Optional[list[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Generate barcode-to-spot-ID mapping for Visium HD composite barcodes. "
        "Computes the cartesian product of two per-part coordinate files."
    )
    parser.add_argument(
        "part1",
        help="TSV file mapping barcode part1 to coordinate values (Y). "
        "Column 0: barcode sequence, columns 1..N: coordinate at each resolution.",
    )
    parser.add_argument(
        "part2",
        help="TSV file mapping barcode part2 to coordinate values (X). "
        "Column 0: barcode sequence, columns 1..N: coordinate at each resolution.",
    )
    parser.add_argument(
        "-o", "--output",
        default="-",
        help="Output TSV file path (default: stdout).",
    )
    parser.add_argument(
        "--prefix",
        nargs="+",
        default=["s_002um_", "s_008um_", "s_016um_"],
        help="Prefix for spot IDs, one per coordinate column. "
        "If a single prefix is given, it is used for all columns. "
        "(default: s_002um_ s_008um_ s_016um_)",
    )
    parser.add_argument(
        "--suffix",
        nargs="+",
        default=[""],
        help="Suffix for spot IDs, one per coordinate column. "
        "If a single suffix is given, it is used for all columns. "
        "(default: empty)",
    )
    parser.add_argument(
        "--delimiter",
        default="_",
        help="Delimiter between Y and X values within a spot ID (default: _).",
    )
    parser.add_argument(
        "--barcode_delimiter",
        default="|",
        help="Delimiter joining barcode part1 and part2 (default: |).",
    )

    parsed = parser.parse_args(args)
    generate_spot_ids(
        parsed.part1,
        parsed.part2,
        parsed.output,
        parsed.prefix,
        parsed.suffix,
        parsed.delimiter,
        parsed.barcode_delimiter,
    )


if __name__ == "__main__":
    main()
