#!/usr/bin/env python3
############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Barcode calling quality assessment for simulated data.

Compares detected barcodes against ground truth embedded in read IDs.
Supports multiple barcode calling modes: tenX_v3, visium_hd, curio, stereo, stereo_split.

Usage:
    python assess_barcode_quality.py --mode stereo --input barcodes.tsv --output metrics.tsv
    python assess_barcode_quality.py --mode stereo_split --input barcodes.tsv --output metrics.tsv
"""

import argparse
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

try:
    import editdistance
    EDITDISTANCE_AVAILABLE = True
except ImportError:
    EDITDISTANCE_AVAILABLE = False


# Mode-specific barcode lengths
MODE_BARCODE_LENGTHS = {
    'tenX_v3': 16,
    'visium_hd': 31,  # 16 + 15 (two-part)
    'curio': 14,      # 8 + 6
    'stereo': 25,
    'stereo_split': 25,
}


def get_score_thresholds(barcode_length: int) -> Tuple[int, List[int]]:
    """
    Compute score thresholds based on barcode length.

    Returns:
        (min_score, list of score thresholds for precision/recall curve)
    """
    min_score = barcode_length - 3
    scores = list(range(min_score - 2, barcode_length + 1))
    return min_score, scores


def extract_ground_truth_single(read_id: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Extract ground truth barcode and UMI from simulated read ID.

    Read ID format: READ_<index>_<transcript>_<BARCODE>_<UMI>_<extra>...

    Returns:
        (barcode, umi) tuple
    """
    parts = read_id.split('_')
    if len(parts) < 4:
        return None, None
    if parts[0] == "PacBio":
        if len(parts) < 9:
            return None, None
        return parts[7], parts[8]
    barcode = parts[3]
    umi = parts[4] if len(parts) > 4 else None
    return barcode, umi


def extract_ground_truth_split(read_id: str, barcode_length: int = 25) -> Set[str]:
    """
    Extract multiple ground truth barcodes from simulated read ID.

    Read ID format: READ_<index>_<transcript>_<BC1>_<UMI1>_<extra>_<BC2>_<UMI2>_...

    Barcodes are at positions 3, 6, 9, ... and identified by length.

    Returns:
        Set of ground truth barcodes
    """
    parts = read_id.split('_')
    barcodes = set()
    for i in range(3, len(parts), 3):
        if i < len(parts) and len(parts[i]) == barcode_length:
            barcodes.add(parts[i])
    return barcodes


def assess_single_barcode_mode(
    input_file: str,
    barcode_length: int,
    barcode_col: int = 1,
    umi_col: int = 2,
    score_col: int = 3,
    min_score: Optional[int] = None
) -> Dict[str, float]:
    """
    Assess barcode calling quality for single-barcode modes.

    Modes: tenX_v3, visium_hd, curio, stereo

    Args:
        input_file: Path to barcode TSV file
        barcode_length: Expected barcode length
        barcode_col: Column index for detected barcode
        umi_col: Column index for detected UMI
        score_col: Column index for barcode score
        min_score: Minimum score threshold for filtering. If None or -1,
                   all barcodes are loaded (* treated as false negative).
    """
    _, score_thresholds = get_score_thresholds(barcode_length)

    # Use provided min_score or None for no filtering
    filter_by_score = min_score is not None and min_score >= 0

    # Statistics
    total_reads = 0
    barcoded_reads = 0
    correct_barcodes = 0

    score_barcoded = defaultdict(int)
    score_correct = defaultdict(int)

    umi_distances = [0] * 15
    true_barcodes = set()
    detected_barcodes = set()

    # Process input file
    with open(input_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            cols = line.strip().split('\t')
            if len(cols) <= max(barcode_col, score_col):
                continue

            total_reads += 1
            detected_bc = cols[barcode_col]
            detected_umi = cols[umi_col] if len(cols) > umi_col else None

            # Extract ground truth
            true_bc, true_umi = extract_ground_truth_single(cols[0])
            if true_bc:
                true_barcodes.add(true_bc)

            # Skip if no barcode detected
            if detected_bc == '*' or not detected_bc:
                continue

            # Parse score
            try:
                score = int(cols[score_col])
            except (ValueError, IndexError):
                score = 0

            # Skip if score is below threshold (when filtering enabled)
            if filter_by_score and score < min_score:
                continue

            barcoded_reads += 1
            detected_barcodes.add(detected_bc)

            score_barcoded[score] += 1

            # Check correctness
            if true_bc and detected_bc == true_bc:
                correct_barcodes += 1
                score_correct[score] += 1

            # UMI edit distance
            if EDITDISTANCE_AVAILABLE and true_umi and detected_umi:
                if len(detected_umi) > 12:
                    detected_umi = detected_umi[:13]
                ed = editdistance.eval(true_umi, detected_umi)
                if ed < 15:
                    umi_distances[ed] += 1

    # Compute metrics
    metrics = {
        'total_reads': total_reads,
        'barcoded_reads': barcoded_reads,
        'correct_barcodes': correct_barcodes,
        'unique_true_barcodes': len(true_barcodes),
        'unique_detected_barcodes': len(detected_barcodes),
    }

    # Precision and recall
    if barcoded_reads > 0:
        metrics['precision'] = 100.0 * correct_barcodes / barcoded_reads
    else:
        metrics['precision'] = 0.0

    if total_reads > 0:
        metrics['recall'] = 100.0 * correct_barcodes / total_reads
    else:
        metrics['recall'] = 0.0

    # Precision/recall at score thresholds
    for threshold in score_thresholds:
        total_at_threshold = sum(score_barcoded[s] for s in range(threshold, barcode_length + 1))
        correct_at_threshold = sum(score_correct[s] for s in range(threshold, barcode_length + 1))

        if total_at_threshold > 0:
            metrics[f'precision_at_{threshold}'] = 100.0 * correct_at_threshold / total_at_threshold
        else:
            metrics[f'precision_at_{threshold}'] = 0.0

        if total_reads > 0:
            metrics[f'recall_at_{threshold}'] = 100.0 * correct_at_threshold / total_reads
        else:
            metrics[f'recall_at_{threshold}'] = 0.0

    # UMI accuracy
    if EDITDISTANCE_AVAILABLE:
        metrics['umi_within_1'] = sum(umi_distances[:2])
        metrics['umi_within_2'] = sum(umi_distances[:3])
        metrics['umi_within_3'] = sum(umi_distances[:4])

    return metrics


def assess_split_barcode_mode(
    input_file: str,
    barcode_length: int = 25,
    barcode_col: int = 1
) -> Dict[str, float]:
    """
    Assess barcode calling quality for stereo_split mode (multiple barcodes per read).
    """
    # Collect all detections per read
    read_detections = defaultdict(list)

    with open(input_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            cols = line.strip().split('\t')
            if len(cols) <= barcode_col:
                continue

            # Extract base read ID (remove trailing position info)
            parts = cols[0].split('_')
            if len(parts) >= 3:
                read_id = '_'.join(parts[:-3]) if len(parts) > 6 else cols[0]
            else:
                read_id = cols[0]

            detected_bc = cols[barcode_col]
            read_detections[read_id].append(detected_bc)

    # Calculate statistics
    total_reads = len(read_detections)
    total_true_barcodes = 0
    total_assignments = 0
    correct_assignments = 0
    incorrect_assignments = 0
    excessive_assignments = 0

    no_barcodes_assigned = 0
    all_correct = 0
    some_correct = 0

    for read_id, detected_barcodes in read_detections.items():
        ground_truth = extract_ground_truth_split(read_id, barcode_length)
        total_true_barcodes += len(ground_truth)

        # Filter out non-detections
        valid_detections = [bc for bc in detected_barcodes if bc != '*' and bc]

        if not valid_detections:
            no_barcodes_assigned += 1
            continue

        # Count correct and incorrect
        correct_in_read = sum(1 for bc in ground_truth if bc in valid_detections)
        incorrect_in_read = sum(1 for bc in valid_detections if bc not in ground_truth)

        total_assignments += min(len(valid_detections), len(ground_truth))
        correct_assignments += correct_in_read
        incorrect_assignments += incorrect_in_read
        excessive_assignments += max(0, len(valid_detections) - len(ground_truth))

        if correct_in_read == len(ground_truth) and len(ground_truth) > 0:
            all_correct += 1
        elif correct_in_read > 0:
            some_correct += 1

    # Compute metrics
    metrics = {
        'total_reads': total_reads,
        'total_true_barcodes': total_true_barcodes,
        'no_barcodes_assigned': no_barcodes_assigned,
        'all_correct_reads': all_correct,
        'some_correct_reads': some_correct,
        'correct_assignments': correct_assignments,
        'incorrect_assignments': incorrect_assignments,
        'excessive_assignments': excessive_assignments,
    }

    # Percentages
    if total_reads > 0:
        metrics['pct_no_barcodes'] = 100.0 * no_barcodes_assigned / total_reads
        metrics['pct_all_correct'] = 100.0 * all_correct / total_reads
        metrics['pct_some_correct'] = 100.0 * some_correct / total_reads

    # Precision and recall
    if total_assignments > 0:
        metrics['assignment_precision'] = 100.0 * correct_assignments / total_assignments
    else:
        metrics['assignment_precision'] = 0.0

    if total_true_barcodes > 0:
        metrics['assignment_recall'] = 100.0 * correct_assignments / total_true_barcodes
    else:
        metrics['assignment_recall'] = 0.0

    return metrics


def write_metrics(metrics: Dict[str, float], output_file: Optional[str] = None):
    """Write metrics to TSV file or stdout."""
    lines = []
    for key, value in sorted(metrics.items()):
        if isinstance(value, float):
            lines.append(f'{key}\t{value:.2f}')
        else:
            lines.append(f'{key}\t{value}')

    output = '\n'.join(lines)

    if output_file:
        with open(output_file, 'w') as f:
            f.write(output + '\n')
    else:
        print(output)


def print_summary(metrics: Dict[str, float], mode: str):
    """Print human-readable summary to stderr."""
    print(f"\n=== Barcode Quality Assessment ({mode}) ===", file=sys.stderr)
    print(f"Total reads: {metrics.get('total_reads', 0)}", file=sys.stderr)

    if mode != 'stereo_split':
        print(f"Barcoded reads: {metrics.get('barcoded_reads', 0)}", file=sys.stderr)
        print(f"Correct barcodes: {metrics.get('correct_barcodes', 0)}", file=sys.stderr)
        print(f"Precision: {metrics.get('precision', 0):.2f}%", file=sys.stderr)
        print(f"Recall: {metrics.get('recall', 0):.2f}%", file=sys.stderr)

        if 'umi_within_1' in metrics:
            print(f"UMIs within 1 edit: {metrics.get('umi_within_1', 0)}", file=sys.stderr)
    else:
        print(f"Total true barcodes: {metrics.get('total_true_barcodes', 0)}", file=sys.stderr)
        print(f"Reads with all correct: {metrics.get('all_correct_reads', 0)} "
              f"({metrics.get('pct_all_correct', 0):.2f}%)", file=sys.stderr)
        print(f"Assignment precision: {metrics.get('assignment_precision', 0):.2f}%", file=sys.stderr)
        print(f"Assignment recall: {metrics.get('assignment_recall', 0):.2f}%", file=sys.stderr)

    print("", file=sys.stderr)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Assess barcode calling quality against simulated ground truth',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Supported modes:
  tenX_v3    - 10x Genomics single-cell v3 (16bp barcode)
  visium_hd    - Visium HD spatial (16+15bp two-part barcode)
  curio        - Curio spatial (14bp = 8+6 split barcode)
  stereo       - Stereo-seq spatial (25bp barcode)
  stereo_split - Stereo-seq split detection mode (multiple barcodes per read)

Ground truth is extracted from read IDs in format:
  READ_<index>_<transcript>_<BARCODE>_<UMI>_...
        """
    )
    parser.add_argument('--mode', '-m', required=True,
                        choices=['tenX_v3', 'visium_hd', 'curio', 'stereo', 'stereo_split'],
                        help='Barcode calling mode')
    parser.add_argument('--input', '-i', required=True,
                        help='Input barcode TSV file from detect_barcodes.py')
    parser.add_argument('--output', '-o',
                        help='Output metrics TSV file (default: stdout)')
    parser.add_argument('--barcode_col', type=int, default=1,
                        help='Column index for detected barcode (default: 1)')
    parser.add_argument('--umi_col', type=int, default=2,
                        help='Column index for detected UMI (default: 2)')
    parser.add_argument('--score_col', type=int, default=3,
                        help='Column index for barcode score (default: 3)')
    parser.add_argument('--barcode_length', type=int,
                        help='Override barcode length (default: auto from mode)')
    parser.add_argument('--min_score', type=int, default=-1,
                        help='Minimum score threshold for filtering barcodes. '
                             'If -1 (default), all barcodes are loaded for QA '
                             '(* treated as false negative)')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress summary output to stderr')

    return parser.parse_args()


def main():
    args = parse_args()

    # Determine barcode length
    if args.barcode_length:
        barcode_length = args.barcode_length
    else:
        barcode_length = MODE_BARCODE_LENGTHS.get(args.mode, 25)

    # Convert min_score: -1 means no filtering (None)
    min_score = args.min_score if args.min_score >= 0 else None

    # Run assessment
    if args.mode == 'stereo_split':
        metrics = assess_split_barcode_mode(
            args.input,
            barcode_length=barcode_length,
            barcode_col=args.barcode_col
        )
    else:
        metrics = assess_single_barcode_mode(
            args.input,
            barcode_length=barcode_length,
            barcode_col=args.barcode_col,
            umi_col=args.umi_col,
            score_col=args.score_col,
            min_score=min_score
        )

    # Output
    write_metrics(metrics, args.output)

    if not args.quiet:
        print_summary(metrics, args.mode)

    # Return exit code based on whether we got any results
    return 0 if metrics.get('total_reads', 0) > 0 else 1


if __name__ == '__main__':
    sys.exit(main())
