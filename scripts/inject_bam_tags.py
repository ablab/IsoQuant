#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Inject barcode (CB) and UMI (UB) tags into a BAM file from barcoded_reads.tsv
# Useful for creating test BAMs that simulate cellranger-style tagged output

import argparse
import logging
import os
import sys

import pysam

logger = logging.getLogger(__name__)


def load_barcoded_reads(barcoded_reads_path: str) -> dict[str, tuple[str, str]]:
    """Load barcoded_reads.tsv into {read_id: (barcode, umi)} dict.

    Skips reads with '*' barcode (unassigned).
    Expected TSV columns: read_id, barcode, UMI, [additional columns ignored]
    """
    barcode_dict: dict[str, tuple[str, str]] = {}
    with open(barcoded_reads_path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("read_id\t"):
                continue
            tokens = line.strip().split("\t")
            if len(tokens) < 3:
                continue
            read_id, barcode, umi = tokens[0], tokens[1], tokens[2]
            if barcode == "*":
                continue
            barcode_dict[read_id] = (barcode, umi)
    logger.info("Loaded %d barcoded reads from %s", len(barcode_dict), barcoded_reads_path)
    return barcode_dict


def inject_tags(
    input_bam: str,
    output_bam: str,
    barcode_dict: dict[str, tuple[str, str]],
    barcode_tag: str = "CB",
    umi_tag: str = "UB",
    suffix: str = "",
) -> None:
    """Read input BAM, add barcode/UMI tags, write output BAM."""
    tagged_count = 0
    total_count = 0

    with pysam.AlignmentFile(input_bam, "rb") as infile:
        with pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
            for read in infile:
                total_count += 1
                read_id = read.query_name
                if read_id in barcode_dict:
                    barcode, umi = barcode_dict[read_id]
                    read.set_tag(barcode_tag, barcode + suffix)
                    if umi and umi != "*":
                        read.set_tag(umi_tag, umi)
                    tagged_count += 1
                outfile.write(read)

    logger.info("Tagged %d / %d reads", tagged_count, total_count)

    logger.info("Indexing output BAM...")
    pysam.index(output_bam)
    logger.info("Done: %s", output_bam)


def main():
    parser = argparse.ArgumentParser(
        description="Inject barcode/UMI tags from barcoded_reads.tsv into a BAM file"
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--barcoded_reads", required=True, help="barcoded_reads.tsv file")
    parser.add_argument("--output", required=True, help="Output BAM file")
    parser.add_argument("--barcode_tag", default="CB", help="BAM tag for barcode (default: CB)")
    parser.add_argument("--umi_tag", default="UB", help="BAM tag for UMI (default: UB)")
    parser.add_argument("--suffix", default="", help="Suffix to append to barcode (e.g. '-1' for cellranger)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    if not os.path.exists(args.bam):
        logger.error("Input BAM not found: %s", args.bam)
        sys.exit(1)
    if not os.path.exists(args.barcoded_reads):
        logger.error("Barcoded reads file not found: %s", args.barcoded_reads)
        sys.exit(1)

    barcode_dict = load_barcoded_reads(args.barcoded_reads)
    inject_tags(args.bam, args.output, barcode_dict, args.barcode_tag, args.umi_tag, args.suffix)


if __name__ == "__main__":
    main()
