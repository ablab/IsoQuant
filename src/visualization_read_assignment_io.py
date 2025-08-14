import logging
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union
import numpy as np

from src.visualization_cache_utils import (
    build_read_assignment_cache_file,
    build_length_effects_cache_file,
    build_length_hist_cache_file,
    save_cache,
    load_cache,
    validate_read_assignment_data,
    validate_length_effects_data,
    validate_length_hist_data,
)


def _smart_open(path_str: str):
    import gzip
    try:
        with open(path_str, 'rb') as bf:
            if bf.read(2) == b'\x1f\x8b':
                return gzip.open(path_str, 'rt')
    except Exception:
        pass
    return open(path_str, 'rt')


def _calc_length_bp(exons_str: str) -> int:
    if not isinstance(exons_str, str) or not exons_str:
        return 0
    total = 0
    for part in exons_str.split(','):
        if '-' not in part:
            continue
        try:
            s, e = part.split('-')
            total += int(e) - int(s) + 1
        except Exception:
            continue
    return total


def get_read_assignment_counts(config, cache_dir: Path):
    """
    Returns read-assignment classification and assignment_type counts, using cache.
    Return format mirrors previous behavior in DictionaryBuilder:
      - If config.read_assignments is a list: ({sample: class_counts}, {sample: assign_type_counts})
      - Else: (class_counts, assign_type_counts)
    """
    logger = logging.getLogger('IsoQuant.visualization.read_assignment_io')
    if not config.read_assignments:
        raise FileNotFoundError("No read assignments file(s) found.")

    cache_file = build_read_assignment_cache_file(
        config.read_assignments, config.ref_only, cache_dir
    )

    if cache_file.exists():
        cached_data = load_cache(cache_file)
        if cached_data and validate_read_assignment_data(cached_data, config.read_assignments):
            logger.info("Using cached read assignment data.")
            if isinstance(config.read_assignments, list):
                return (
                    cached_data["classification_counts"],
                    cached_data["assignment_type_counts"],
                )
            return cached_data

    logger.info("Building read assignment data from scratch.")

    def process_file(file_path: str):
        classification_counts: Dict[str, int] = {}
        assignment_type_counts: Dict[str, int] = {}
        with _smart_open(file_path) as fh:
            # Skip header lines starting with '#'
            while True:
                pos = fh.tell()
                line = fh.readline()
                if not line:
                    break
                if not line.startswith('#'):
                    fh.seek(pos)
                    break
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                additional_info = parts[8]
                assignment_type = parts[5]
                classification = (
                    additional_info.split('Classification=')[-1].split(';')[0].strip()
                    if 'Classification=' in additional_info else 'Unknown'
                )
                classification_counts[classification] = classification_counts.get(classification, 0) + 1
                assignment_type_counts[assignment_type] = assignment_type_counts.get(assignment_type, 0) + 1
        return classification_counts, assignment_type_counts

    if isinstance(config.read_assignments, list):
        classification_counts_dict: Dict[str, Dict[str, int]] = {}
        assignment_type_counts_dict: Dict[str, Dict[str, int]] = {}
        for sample_name, file_path in config.read_assignments:
            c_counts, a_counts = process_file(file_path)
            classification_counts_dict[sample_name] = c_counts
            assignment_type_counts_dict[sample_name] = a_counts
        to_cache = {
            "classification_counts": classification_counts_dict,
            "assignment_type_counts": assignment_type_counts_dict,
        }
        save_cache(cache_file, to_cache)
        return classification_counts_dict, assignment_type_counts_dict
    else:
        counts = process_file(config.read_assignments)
        save_cache(cache_file, counts)
        return counts


def get_read_length_effects(config, cache_dir: Path) -> Dict[str, Any]:
    """
    Compute and cache read-length effects aggregates:
      - by length bin vs assignment_type and vs classification
      - dynamic keys for observed categories
    Returns dict with keys: bins, by_bin_assignment, by_bin_classification, assignment_keys, classification_keys, totals
    """
    logger = logging.getLogger('IsoQuant.visualization.read_assignment_io')
    if not config.read_assignments:
        raise FileNotFoundError("No read assignments file(s) found.")

    # Fixed bin order focused on 0-15 kb
    bin_order = ['<1kb','1-2kb','2-3kb','3-4kb','4-5kb','5-6kb','6-7kb','7-8kb','8-9kb','9-10kb','10-12kb','12-15kb','>15kb']
    cache_file = build_length_effects_cache_file(
        config.read_assignments, config.ref_only, cache_dir, bin_order
    )

    if cache_file.exists():
        cached = load_cache(cache_file)
        if cached and validate_length_effects_data(cached, expected_bins=bin_order):
            logger.info("Using cached read length effects.")
            return cached

    from collections import defaultdict
    by_bin_assignment: Dict[str, Dict[str, int]] = {b: defaultdict(int) for b in bin_order}
    by_bin_classification: Dict[str, Dict[str, int]] = {b: defaultdict(int) for b in bin_order}
    assignment_keys = set()
    classification_keys = set()
    totals: Dict[str, int] = {b: 0 for b in bin_order}

    def assign_bin(length_bp: int) -> str:
        if length_bp < 1000: return '<1kb'
        if length_bp < 2000: return '1-2kb'
        if length_bp < 3000: return '2-3kb'
        if length_bp < 4000: return '3-4kb'
        if length_bp < 5000: return '4-5kb'
        if length_bp < 6000: return '5-6kb'
        if length_bp < 7000: return '6-7kb'
        if length_bp < 8000: return '7-8kb'
        if length_bp < 9000: return '8-9kb'
        if length_bp < 10000: return '9-10kb'
        if length_bp < 12000: return '10-12kb'
        if length_bp < 15000: return '12-15kb'
        return '>15kb'

    def process_file(file_path: str):
        with _smart_open(file_path) as fh:
            # Skip header lines
            while True:
                pos = fh.tell()
                line = fh.readline()
                if not line:
                    break
                if not line.startswith('#'):
                    fh.seek(pos)
                    break
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                assignment_type = parts[5]
                exons_str = parts[7]
                addi = parts[8]
                classification = (
                    addi.split('Classification=')[-1].split(';')[0].strip()
                    if 'Classification=' in addi else 'unknown'
                )
                length_bp = _calc_length_bp(exons_str)
                b = assign_bin(length_bp)
                totals[b] += 1
                by_bin_assignment[b][assignment_type] += 1
                by_bin_classification[b][classification] += 1
                assignment_keys.add(assignment_type)
                classification_keys.add(classification)

    if isinstance(config.read_assignments, list):
        for _sample, file_path in config.read_assignments:
            process_file(file_path)
    else:
        process_file(config.read_assignments)

    # Convert defaultdicts to dicts for safer pickling/validation
    by_bin_assignment = {b: dict(d) for b, d in by_bin_assignment.items()}
    by_bin_classification = {b: dict(d) for b, d in by_bin_classification.items()}

    result = {
        'bins': bin_order,
        'by_bin_assignment': by_bin_assignment,
        'by_bin_classification': by_bin_classification,
        'assignment_keys': sorted(list(assignment_keys)),
        'classification_keys': sorted(list(classification_keys)),
        'totals': totals,
    }
    save_cache(cache_file, result)
    return result


def get_read_length_histogram(config, cache_dir: Path, bin_edges: List[int] = None) -> Dict[str, Any]:
    """
    Compute and cache a histogram of read lengths derived from the exons column.
    Returns dict with keys: edges, counts, total
    """
    logger = logging.getLogger('IsoQuant.visualization.read_assignment_io')
    if not config.read_assignments:
        raise FileNotFoundError("No read assignments file(s) found.")

    if bin_edges is None:
        bin_edges = [
            0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000,
            5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000,
            12000, 15000,
        ]

    cache_file = build_length_hist_cache_file(
        config.read_assignments, config.ref_only, cache_dir, bin_edges
    )

    if cache_file.exists():
        cached = load_cache(cache_file)
        if cached and validate_length_hist_data(cached, expected_edges=bin_edges):
            logger.info("Using cached read length histogram.")
            return cached

    lengths: List[int] = []

    def process_file(file_path: str):
        with _smart_open(file_path) as fh:
            for line in fh:
                if not line or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 8:
                    continue
                exons = parts[7]
                lengths.append(_calc_length_bp(exons))

    if isinstance(config.read_assignments, list):
        for _sample, file_path in config.read_assignments:
            process_file(file_path)
    else:
        process_file(config.read_assignments)

    counts, edges = np.histogram(np.array(lengths, dtype=np.int64), bins=np.array(bin_edges))
    result = {
        'edges': edges.tolist(),
        'counts': counts.tolist(),
        'total': int(len(lengths)),
    }
    save_cache(cache_file, result)
    return result


