#!/usr/bin/env python3
import argparse
import csv
import os

def load_mapping(file):
    """Load simple 2-column mapping into dict."""
    mapping = {}
    with open(file, newline='') as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            mapping[row[0]] = row[1]
    return mapping

def load_spot_map(file):
    """Load mapping 2um -> (8um, 16um)."""
    mapping = {}
    with open(file, newline='') as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            spot2um, spot8um, spot16um = row
            mapping[spot2um[:-2]] = {"8um": spot8um[:-2], "16um": spot16um[:-2]}
    return mapping

def convert_barcode(barcode, part1_to_y, part2_to_x, spot_map):
    """Convert barcode (p1|p2) to 2um/8um/16um spot ids."""
    try:
        p1, p2 = barcode.split("|")
    except ValueError:
        return None, None, None  # malformed

    if p1 not in part1_to_y or p2 not in part2_to_x:
        return None, None, None

    x = part2_to_x[p2]
    y = part1_to_y[p1]

    spot2um = f"s_002um_{x}_{y}"
    if spot2um not in spot_map:
        return spot2um, None, None

    return spot2um, spot_map[spot2um]["8um"], spot_map[spot2um]["16um"]

def transform_reads(reads_file, barcode_col, part1_to_y, part2_to_x, spot_map, output_prefix):
    with open(reads_file, newline='') as infile, \
         open(f"{output_prefix}_2um.tsv", "w", newline='') as f2um, \
         open(f"{output_prefix}_8um.tsv", "w", newline='') as f8um, \
         open(f"{output_prefix}_16um.tsv", "w", newline='') as f16um:

        reader = csv.reader(infile, delimiter="\t")
        w2 = csv.writer(f2um, delimiter="\t")
        w8 = csv.writer(f8um, delimiter="\t")
        w16 = csv.writer(f16um, delimiter="\t")

        header = next(reader)
        w2.writerow(header)
        w8.writerow(header)
        w16.writerow(header)

        for row in reader:
            barcode = row[barcode_col]
            spot2um, spot8um, spot16um = convert_barcode(barcode, part1_to_y, part2_to_x, spot_map)

            row2, row8, row16 = row[:], row[:], row[:]
            row2[barcode_col] = spot2um if spot2um else "*"
            row8[barcode_col] = spot8um if spot8um else "*"
            row16[barcode_col] = spot16um if spot16um else "*"

            w2.writerow(row2)
            w8.writerow(row8)
            w16.writerow(row16)

def main():
    parser = argparse.ArgumentParser(
        description="Replace barcodes with 2um, 8um, 16um spot IDs in separate TSV files."
    )
    parser.add_argument("reads_tsv", help="Input TSV file with read_id, barcode, other columns")
    parser.add_argument("part1_to_y", help="TSV mapping barcode part1 -> Y coordinate")
    parser.add_argument("part2_to_x", help="TSV mapping barcode part2 -> X coordinate")
    parser.add_argument("spot_map", help="TSV mapping 2um -> 8um,16um")
    parser.add_argument("-c", "--barcode-col", type=int, default=2,
                        help="Column index (1-based) of barcode in reads file (default=2)")
    parser.add_argument("-o", "--output-prefix", default="reads_converted",
                        help="Prefix for output files (default=reads_converted)")
    args = parser.parse_args()

    # convert to 0-based index
    barcode_col = args.barcode_col - 1

    part1_to_y = load_mapping(args.part1_to_y)
    part2_to_x = load_mapping(args.part2_to_x)
    spot_map = load_spot_map(args.spot_map)

    transform_reads(args.reads_tsv, barcode_col, part1_to_y, part2_to_x, spot_map, args.output_prefix)

if __name__ == "__main__":
    main()
