############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Shared polyA / TSS peak detection.

Single source of truth for the peak-detection tunables and a reusable
detector that turns a per-feature histogram of terminal read positions into a
set of accepted peaks. Both the side-output counters (`terminal_counter.py`)
and transcript discovery (`intron_graph.py`,
`graph_based_model_construction.py`) consume this module so the feature
computation cannot drift between training-time and the two inference sites.

The detection mirrors `TerminalCounter`'s pipeline for a single histogram:
build a zero-padded base-1 histogram, compute summary features, run
`scipy.signal.find_peaks`, fall back to the distribution mode when no peak is
found, rank peaks by height, and keep the ones the trained XGBoost classifier
accepts. Coordinates returned are genomic.
"""

import logging
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional

import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import find_peaks, peak_widths
from xgboost import XGBClassifier

logger = logging.getLogger('IsoQuant')

_DATA_DIR = Path(__file__).parent / "data"
POLYA_MODEL_PATH = _DATA_DIR / "model_polya.json"
TSS_MODEL_PATH = _DATA_DIR / "model_tss.json"

# Peak finding + classification parameters. Match values used to train the
# shipped XGBoost models -- changing these without retraining will break the
# feature alignment between fit time and inference time.
PEAK_DISTANCE = 10
PEAK_REL_HEIGHT = 0.98
HISTOGRAM_PAD = 10
ANNOTATION_TOLERANCE = 10

# XGBoost feature columns (must match training).
FEATURE_COLUMNS = ['var', 'skew', 'peak_count', 'peak_width', 'entropy',
                   'mean_height', 'peak_heights', 'relative_height']

# Process-level model cache. The models are loaded lazily on first use so the
# XGBoost / OpenMP runtime is only initialized inside the fork() workers that
# actually call detect_peaks -- never in the parent, where the inherited
# OpenMP state would deadlock the pool. One instance per worker process.
_MODEL_CACHE: Dict[str, XGBClassifier] = {}


def _get_model(model_path: Path) -> XGBClassifier:
    key = str(model_path)
    model = _MODEL_CACHE.get(key)
    if model is None:
        model = XGBClassifier()
        model.load_model(key)
        _MODEL_CACHE[key] = model
    return model


def get_polya_model() -> XGBClassifier:
    return _get_model(POLYA_MODEL_PATH)


def get_tss_model() -> XGBClassifier:
    return _get_model(TSS_MODEL_PATH)


class Peak(NamedTuple):
    """An accepted terminal-position peak in genomic coordinates."""
    position: int   # predicted site
    count: int      # reads supporting the peak window [left, right]
    left: int       # genomic left bound of the peak window
    right: int      # genomic right bound of the peak window


def detect_peaks(position_counts: Dict[int, int],
                 model: XGBClassifier) -> List[Peak]:
    """Detect accepted peaks in a single ``{genomic_position: count}`` histogram."""
    return detect_peaks_batch([position_counts], model)[0]


def detect_peaks_batch(histograms: List[Optional[Dict[int, int]]],
                       model: XGBClassifier) -> List[List[Peak]]:
    """Detect peaks for several histograms with a single batched ``predict``.

    Returns a list parallel to ``histograms``; each element is the list of
    accepted :class:`Peak` for that histogram (empty when the histogram is
    empty or every peak is rejected).
    """
    results: List[List[Peak]] = [[] for _ in histograms]
    records: Dict[int, dict] = {}
    feature_rows: List[list] = []
    row_owner: List[tuple] = []  # (histogram_index, peak_dict)

    for i, counts in enumerate(histograms):
        if not counts:
            continue
        record = _histogram_record(counts)
        records[i] = record

        if record['peak_count'] == 0:
            # No detected peak: fall back to the distribution mode (always
            # accepted, never scored by the classifier).
            mode = record['mode']
            position = record['start'] + mode
            count = _window_count(record['histogram'], mode, mode)
            results[i].append(Peak(position, count, position, position))
            continue

        for peak in _rank_peaks(record):
            feature_rows.append([record['var'], record['skew'],
                                 record['peak_count'], peak['peak_width'],
                                 record['entropy'], record['mean_height'],
                                 peak['peak_heights'], peak['relative_height']])
            row_owner.append((i, peak))

    if feature_rows:
        features = pd.DataFrame(feature_rows, columns=FEATURE_COLUMNS).astype(float)
        predicted = model.predict(features)
        for (i, peak), keep in zip(row_owner, predicted):
            if int(keep) != 1:
                continue
            record = records[i]
            start = record['start']
            position = start + int(peak['peak_location'])
            count = _window_count(record['histogram'],
                                  peak['peak_left'], peak['peak_right'])
            results[i].append(Peak(position, count,
                                   start + int(peak['peak_left']),
                                   start + int(peak['peak_right'])))
    return results


def _histogram_record(position_counts: Dict[int, int]) -> dict:
    """Build the per-histogram features and raw peaks (mirrors
    ``TerminalCounter._build_peak_dataframe`` for one transcript)."""
    data: List[int] = []
    for position, count in position_counts.items():
        if count > 0:
            data.extend([int(position)] * int(count))

    data_min = min(data)
    data_max = max(data)
    hist_counts = np.histogram(
        data, bins=(data_max - data_min) + 1, range=(data_min, data_max))[0]
    padded = [0] * HISTOGRAM_PAD + list(hist_counts) + [0] * HISTOGRAM_PAD

    var = float(np.var(data))
    record = {
        'start': data_min,
        'mode': int(stats.mode(data, keepdims=False).mode) - data_min,
        'var': var,
        'skew': (float(stats.skew(data))
                 if len(data) > 3 and var >= 1e-12 else float('nan')),
        'entropy': float(stats.entropy(padded)),
        'mean_height': float(np.mean(padded)),
        'histogram': padded,
    }

    # find_peaks returns padded-array indices; subtract HISTOGRAM_PAD to get
    # start-relative coordinates.
    peaks_idx = find_peaks(padded, distance=PEAK_DISTANCE)[0]
    record['peak_count'] = len(peaks_idx)
    if record['peak_count'] == 0:
        return record

    widths = peak_widths(padded, peaks_idx, rel_height=PEAK_REL_HEIGHT)
    peak_location = [int(j - HISTOGRAM_PAD) for j in peaks_idx]
    record['raw_peaks'] = {
        'peak_location': peak_location,
        'peak_width': list(widths[0]),
        'peak_left': [int(j - HISTOGRAM_PAD) for j in widths[2]],
        'peak_right': [int(j - HISTOGRAM_PAD) for j in widths[3]],
        'peak_heights': [padded[int(p) + HISTOGRAM_PAD] for p in peak_location],
    }
    return record


def _rank_peaks(record: dict) -> List[dict]:
    """Order a transcript's peaks by descending height and attach
    ``relative_height`` (mirrors ``TerminalCounter._rank_peaks``)."""
    raw = record['raw_peaks']
    heights = np.asarray(raw['peak_heights'], dtype=float)
    if heights.size > 1:
        order = np.argsort(-heights)
        top = float(heights[order][0])
    else:
        order = [0]
        top = None

    ranked = []
    for idx in order:
        height = raw['peak_heights'][idx]
        ranked.append({
            'peak_location': raw['peak_location'][idx],
            'peak_width': raw['peak_width'][idx],
            'peak_left': raw['peak_left'][idx],
            'peak_right': raw['peak_right'][idx],
            'peak_heights': height,
            'relative_height': (1.0 if top is None
                                else (height / top if top else 0.0)),
        })
    return ranked


def _window_count(histogram: List[int], left: int, right: int) -> int:
    """Reads inside the padded-histogram window ``[left, right]`` (mirrors
    ``TerminalCounter._counts_for_peak``); left/right are start-relative."""
    low = max(0, int(left) + HISTOGRAM_PAD)
    high = min(len(histogram), int(right) + HISTOGRAM_PAD + 1)
    if high <= low:
        return 0
    return int(np.asarray(histogram[low:high]).sum())
