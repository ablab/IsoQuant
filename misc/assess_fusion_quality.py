#!/usr/bin/env python3

############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import argparse
import logging
import sys
from typing import Set, Tuple

logger = logging.getLogger('FusionQA')


def setup_logger():
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Assess fusion detection quality against a truth set")
    parser.add_argument("--fusion_tsv", type=str, required=True,
                        help="IsoQuant fusion output TSV file")
    parser.add_argument("--truth", type=str, required=True,
                        help="Truth set TSV (columns: fusion_id, left_transcript_id, "
                             "left_gene, ..., right_gene, ...)")
    parser.add_argument("--output", "-o", type=str, required=True,
                        help="Output quality report file (TSV)")
    parser.add_argument("--min_support", type=int, default=2,
                        help="Minimum supporting reads to consider a prediction (default: 2)")
    return parser.parse_args()


def load_truth_set(truth_path: str) -> Set[Tuple[str, str]]:
    truth_pairs: Set[Tuple[str, str]] = set()
    with open(truth_path) as f:
        header = f.readline().strip().split('\t')
        left_idx = header.index("left_gene")
        right_idx = header.index("right_gene")
        for line in f:
            tokens = line.strip().split('\t')
            if len(tokens) <= max(left_idx, right_idx):
                continue
            pair = tuple(sorted([tokens[left_idx], tokens[right_idx]]))
            truth_pairs.add(pair)
    logger.info("Loaded %d unique fusion pairs from truth set" % len(truth_pairs))
    return truth_pairs


def load_predictions(fusion_tsv: str, min_support: int = 2) -> Set[Tuple[str, str]]:
    predicted_pairs: Set[Tuple[str, str]] = set()
    total_lines = 0
    filtered_lines = 0
    with open(fusion_tsv) as f:
        header = f.readline().strip().split('\t')
        left_gene_idx = header.index("LeftGene")
        right_gene_idx = header.index("RightGene")
        support_idx = header.index("SupportingReads")
        for line in f:
            tokens = line.strip().split('\t')
            if len(tokens) <= max(left_gene_idx, right_gene_idx, support_idx):
                continue
            total_lines += 1
            support = int(tokens[support_idx])
            if support < min_support:
                filtered_lines += 1
                continue
            pair = tuple(sorted([tokens[left_gene_idx], tokens[right_gene_idx]]))
            predicted_pairs.add(pair)
    logger.info("Loaded %d unique predicted fusion pairs (%d total calls, %d filtered by support < %d)" %
                (len(predicted_pairs), total_lines, filtered_lines, min_support))
    return predicted_pairs


def compute_metrics(predicted: Set[Tuple[str, str]], truth: Set[Tuple[str, str]]):
    tp = len(predicted & truth)
    fp = len(predicted - truth)
    fn = len(truth - predicted)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = (2 * precision * sensitivity / (precision + sensitivity)
           if (precision + sensitivity) > 0 else 0.0)

    return {
        "true_positives": tp,
        "false_positives": fp,
        "false_negatives": fn,
        "sensitivity": sensitivity,
        "precision": precision,
        "f1_score": f1,
        "truth_set_size": len(truth),
        "predicted_count": len(predicted),
    }


def main():
    setup_logger()
    args = parse_args()

    truth_pairs = load_truth_set(args.truth)
    predicted_pairs = load_predictions(args.fusion_tsv, args.min_support)
    metrics = compute_metrics(predicted_pairs, truth_pairs)

    logger.info("Results: TP=%d, FP=%d, FN=%d" %
                (metrics["true_positives"], metrics["false_positives"], metrics["false_negatives"]))
    logger.info("Sensitivity=%.4f, Precision=%.4f, F1=%.4f" %
                (metrics["sensitivity"], metrics["precision"], metrics["f1_score"]))

    with open(args.output, "w") as f:
        for k, v in metrics.items():
            if isinstance(v, float):
                f.write("%s\t%.4f\n" % (k, v))
            else:
                f.write("%s\t%d\n" % (k, v))

    logger.info("Report written to %s" % args.output)

    # Print missed fusions for debugging
    missed = truth_pairs - predicted_pairs
    if missed:
        logger.info("Top missed fusions (FN): %s" %
                    ", ".join("%s-%s" % (a, b) for a, b in sorted(missed)[:10]))

    # Print false positives for debugging
    false_pos = predicted_pairs - truth_pairs
    if false_pos:
        logger.info("Top false positives (FP): %s" %
                    ", ".join("%s-%s" % (a, b) for a, b in sorted(false_pos)[:10]))

    return 0


if __name__ == "__main__":
    sys.exit(main())
