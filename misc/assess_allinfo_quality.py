#!/usr/bin/env python3

############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Assess quality of UMI-filtered allinfo output by independently computing
# statistics from the allinfo file and cross-validating against stats.tsv.

import argparse
import gzip
import logging
import sys
from collections import Counter
from typing import Dict, Tuple

logger = logging.getLogger('AllinfoQA')


def setup_logger():
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Assess quality of UMI-filtered allinfo output")
    parser.add_argument("--allinfo", type=str, required=True,
                        help="Path to allinfo file (plain or .gz)")
    parser.add_argument("--stats", type=str, required=True,
                        help="Path to stats.tsv file from UMI filtering")
    parser.add_argument("--output", "-o", type=str, required=True,
                        help="Output quality report file (TSV)")
    return parser.parse_args()


def load_tsv_config(filepath: str) -> Dict[str, str]:
    config = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            tokens = line.strip().split('\t')
            if len(tokens) >= 2:
                config[tokens[0]] = tokens[1]
    return config


def compute_allinfo_stats(allinfo_path: str) -> Tuple[Dict[str, int], int]:
    """Parse allinfo file and compute statistics.

    Returns:
        Tuple of (stats_dict, duplicate_triplet_count)
    """
    open_fn = gzip.open if allinfo_path.endswith('.gz') else open
    mode = 'rt' if allinfo_path.endswith('.gz') else 'r'

    total_reads = 0
    spliced_reads = 0
    gene_barcode_pairs: set = set()
    triplet_counts: Counter = Counter()

    with open_fn(allinfo_path, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            if len(fields) < 13:
                logger.warning("Skipping malformed line with %d fields: %s",
                               len(fields), line[:100])
                continue

            total_reads += 1

            # Columns (0-indexed):
            # 0: read_id, 1: gene_id, 2: cell_type, 3: barcode, 4: umi,
            # 5: introns, 6: TSS, 7: polyA, 8: exons, 9: read_type,
            # 10: num_introns, 11: transcript_id, 12: transcript_type
            gene_id = fields[1]
            barcode = fields[3]
            umi = fields[4]
            num_introns = int(fields[10])

            if num_introns > 0:
                spliced_reads += 1

            gene_barcode_pairs.add((gene_id, barcode))
            triplet_counts[(gene_id, barcode, umi)] += 1

    duplicate_triplets = sum(1 for count in triplet_counts.values() if count > 1)

    stats = {
        "Total reads saved": total_reads,
        "Spliced reads saved": spliced_reads,
        "Unique gene-barcodes pairs": len(gene_barcode_pairs),
    }
    return stats, duplicate_triplets


def cross_validate(computed_stats: Dict[str, int],
                   stats_tsv: Dict[str, str]) -> bool:
    """Cross-validate computed allinfo stats against stats.tsv values.

    Returns True if all checks pass.
    """
    all_match = True

    checks = [
        ("Total reads saved", "Total reads saved"),
        ("Spliced reads saved", "Spliced reads saved"),
        ("Unique gene-barcodes pairs", "Unique gene-barcodes pairs"),
    ]

    for computed_key, stats_key in checks:
        if stats_key not in stats_tsv:
            logger.error("Key '%s' not found in stats.tsv", stats_key)
            all_match = False
            continue

        computed_val = computed_stats[computed_key]
        stats_val = int(stats_tsv[stats_key])

        if computed_val != stats_val:
            logger.error("Mismatch for '%s': allinfo=%d, stats.tsv=%d",
                         computed_key, computed_val, stats_val)
            all_match = False
        else:
            logger.info("OK '%s': %d", computed_key, computed_val)

    return all_match


def main():
    args = parse_args()
    setup_logger()

    logger.info("Loading allinfo from %s", args.allinfo)
    computed_stats, duplicate_triplets = compute_allinfo_stats(args.allinfo)
    logger.info("Computed stats from allinfo: %s", computed_stats)

    logger.info("Loading stats from %s", args.stats)
    stats_tsv = load_tsv_config(args.stats)
    logger.info("Stats.tsv contains %d entries", len(stats_tsv))

    # Cross-validate
    stats_match = cross_validate(computed_stats, stats_tsv)

    if duplicate_triplets > 0:
        logger.error("Found %d duplicate gene-barcode-UMI triplets in allinfo",
                      duplicate_triplets)
    else:
        logger.info("No duplicate gene-barcode-UMI triplets found")

    # Write quality report
    exit_code = 0
    with open(args.output, "w") as outf:
        # Pass through all stats.tsv entries
        for key, val in stats_tsv.items():
            outf.write("%s\t%s\n" % (key, val))

        # Add computed validation metrics
        outf.write("duplicate_triplets\t%d\n" % duplicate_triplets)
        allinfo_match_val = 1 if stats_match else 0
        outf.write("allinfo_stats_match\t%d\n" % allinfo_match_val)

    logger.info("Quality report written to %s", args.output)

    if not stats_match:
        logger.error("FAILED: allinfo stats do not match stats.tsv")
        exit_code = 1
    if duplicate_triplets > 0:
        logger.error("FAILED: duplicate gene-barcode-UMI triplets detected")
        exit_code = 1

    if exit_code == 0:
        logger.info("All checks passed")

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
