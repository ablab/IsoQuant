############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Unit tests for the shared polyA / TSS peak detector (terminal_peaks).

The detector is the single source of truth used both by the side-output
prediction counters (terminal_counter) and by transcript discovery
(intron_graph / graph_based_model_construction). These tests pin its output
against terminal_counter's own pandas pipeline so the feature computation can
not silently drift between the two call sites.
"""

from collections import Counter

import numpy as np
import pytest

from isoquant_lib import terminal_counter as tc
from isoquant_lib import terminal_peaks as tp
from isoquant_lib.isoform_assignment import ReadAssignmentType

# Reuse the counter test harness (stub model + read-assignment builders).
from isoquant_tests.test_polya_prediction import (
    _StubModel,
    _make_args,
    _make_read_assignment,
    _read_tsv,
    stub_model,  # noqa: F401  (pytest fixture)
)


def _counter_peaks(positions, tmp_path, name):
    """Run the real PolyACounter pipeline over reads at ``positions`` and
    return the set of (prediction, counts) it emits."""
    out = tmp_path / name
    counter = tc.PolyACounter(_make_args(), str(out))
    exons = [(100, 200), (300, 400)]
    for pos in positions:
        counter.add_read_info(_make_read_assignment(polya_pos=pos, exons=exons))
    counter.dump()
    df = _read_tsv(out)
    if df.empty:
        return set()
    return set(zip(df["prediction"].astype(int), df["counts"].astype(int)))


def _detect_peaks(positions, accept=True):
    histogram = dict(Counter(positions))
    peaks = tp.detect_peaks(histogram, _StubModel(accept=accept))
    return {(p.position, p.count) for p in peaks}


@pytest.mark.parametrize("positions", [
    [420],                                                   # single read
    [420 + i % 3 for i in range(30)],                       # one tight peak
    [420 + i % 3 for i in range(30)] +
    [500 + i % 3 for i in range(20)],                        # two peaks (ranking)
    [400, 401],                                              # adjacent-bin plateau
])
def test_detect_peaks_parity_with_counter(stub_model, tmp_path, positions):
    counter_peaks = _counter_peaks(positions, tmp_path, "parity.tsv")
    assert _detect_peaks(positions) == counter_peaks
    # And the union of supporting counts is conserved either way.
    assert sum(c for _, c in _detect_peaks(positions)) == \
           sum(c for _, c in counter_peaks)


def test_detect_peaks_empty_returns_empty():
    assert tp.detect_peaks({}, _StubModel(accept=True)) == []


def test_detect_peaks_batch_matches_single():
    histograms = [
        {420 + i % 3: 1 for i in range(1)},
        {420 + i % 3: 10 for i in range(3)},
        {},
    ]
    model = _StubModel(accept=True)
    batched = tp.detect_peaks_batch(histograms, model)
    singles = [tp.detect_peaks(h, model) for h in histograms]
    assert batched == singles
    assert batched[2] == []


def test_detect_peaks_rejected_peak_dropped():
    # A clear peak is dropped when the model rejects it (no zero-peak fallback).
    positions = [420 + i % 3 for i in range(30)]
    assert _detect_peaks(positions, accept=False) == set()
