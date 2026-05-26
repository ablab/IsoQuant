############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Unit tests for the polyA / TSS prediction counters.

The XGBoost classifier is patched with a stub that returns ``1`` for every
peak. This keeps the tests independent of the shipped model file and lets us
exercise the counter's bookkeeping, feature extraction, and output layout.
"""

import os
import types
from argparse import Namespace

import numpy as np
import pandas as pd
import pytest

from isoquant_lib import terminal_counter as tc
from isoquant_lib.isoform_assignment import ReadAssignmentType


class _StubModel:
    """Minimal stand-in for ``XGBClassifier`` used at inference time."""

    def __init__(self, accept=True):
        self._accept = accept

    def load_model(self, _path):
        return None

    def predict(self, X):
        value = 1 if self._accept else 0
        return np.full(len(X), value, dtype=int)


@pytest.fixture
def stub_model(monkeypatch):
    model = _StubModel(accept=True)
    monkeypatch.setattr(tc, "XGBClassifier", lambda *a, **kw: model)
    return model


def _make_args():
    return Namespace(genedb=True, fl_data=True, read_group=None)


def _make_gene_info(transcript_id="T1", exons=((100, 200), (300, 400))):
    return types.SimpleNamespace(all_isoforms_exons={transcript_id: list(exons)})


def _make_isoform_match(transcript_id="T1", gene_id="G1"):
    return types.SimpleNamespace(assigned_transcript=transcript_id,
                                 assigned_gene=gene_id)


def _make_polya_info(external_polya_pos=-1, external_polyt_pos=-1):
    return types.SimpleNamespace(
        external_polya_pos=external_polya_pos,
        external_polyt_pos=external_polyt_pos,
        internal_polya_pos=-1,
        internal_polyt_pos=-1,
    )


def _make_read_assignment(
    *, strand="+", chr_id="chr1", polya_pos=420, exons=None,
    assignment_type=ReadAssignmentType.unique, polyA_found=True,
    transcript_id="T1", gene_id="G1", read_group_ids=None,
    gene_info=None,
):
    exons = exons or [(100, 200), (300, polya_pos)]
    polya_info = _make_polya_info(
        external_polya_pos=polya_pos if strand == "+" else -1,
        external_polyt_pos=polya_pos if strand == "-" else -1,
    )
    return types.SimpleNamespace(
        assignment_type=assignment_type,
        strand=strand,
        chr_id=chr_id,
        polyA_found=polyA_found,
        polya_info=polya_info,
        corrected_exons=exons,
        isoform_matches=[_make_isoform_match(transcript_id, gene_id)],
        read_group_ids=list(read_group_ids) if read_group_ids else [],
        gene_info=gene_info or _make_gene_info(transcript_id, exons),
    )


def _read_tsv(path):
    return pd.read_csv(path, sep="\t")


def test_polya_counter_filters_invalid_reads(stub_model, tmp_path):
    counter = tc.PolyACounter(_make_args(), str(tmp_path / "p.tsv"))

    # polyA_found=False is rejected
    counter.add_read_info(_make_read_assignment(polyA_found=False))
    # Wrong assignment type is rejected
    counter.add_read_info(_make_read_assignment(
        assignment_type=ReadAssignmentType.noninformative))
    # Strand without matching polya_info position is rejected
    bad = _make_read_assignment(strand="+", polya_pos=420)
    bad.polya_info = _make_polya_info(external_polya_pos=-1)
    counter.add_read_info(bad)

    assert counter.transcripts == {}


def test_polya_counter_accumulates_positions_per_group(stub_model, tmp_path):
    counter = tc.PolyACounter(_make_args(), str(tmp_path / "p.tsv"),
                              string_pools=object(), group_index=0)
    # Single transcript, three reads in two groups
    for pos, gid in [(410, 0), (411, 0), (455, 1)]:
        counter.add_read_info(_make_read_assignment(
            polya_pos=pos, read_group_ids=[gid]))

    entry = counter.transcripts["T1"]
    assert entry["data"] == [410, 411, 455]
    assert entry[0] == [410, 411]
    assert entry[1] == [455]


def test_tss_counter_extracts_strand_aware_start(stub_model, tmp_path):
    counter = tc.TSSCounter(_make_args(), str(tmp_path / "t.tsv"))
    counter.add_read_info(_make_read_assignment(strand="+",
                                                exons=[(105, 200), (300, 400)]))
    counter.add_read_info(_make_read_assignment(strand="-",
                                                exons=[(100, 200), (300, 405)]))
    assert counter.transcripts["T1"]["data"] == [105, 405]


def test_dump_empty_writes_header_only(stub_model, tmp_path):
    out = tmp_path / "empty.tsv"
    counter = tc.PolyACounter(_make_args(), str(out))
    counter.dump()

    df = _read_tsv(out)
    assert list(df.columns) == tc.EMPTY_COLUMNS
    assert df.empty


def test_dump_ungrouped_emits_one_row_per_peak(stub_model, tmp_path):
    out = tmp_path / "ungrouped.tsv"
    counter = tc.PolyACounter(_make_args(), str(out))
    # Cluster 30 reads tightly around 420 to form a clear peak that the stub
    # model accepts; annotated end is exon[-1][1]+1 = 401, so the peak is Novel.
    exons = [(100, 200), (300, 400)]
    for offset in range(30):
        counter.add_read_info(_make_read_assignment(
            polya_pos=420 + offset % 3, exons=exons))
    counter.dump()

    df = _read_tsv(out)
    assert list(df.columns) == tc.EMPTY_COLUMNS
    assert len(df) >= 1
    assert (df["transcript_id"] == "T1").all()
    assert df["counts"].sum() == 30
    assert df["flag"].iloc[0] == "Novel"


def test_dump_flags_peak_as_known_within_tolerance(stub_model, tmp_path):
    out = tmp_path / "known.tsv"
    counter = tc.PolyACounter(_make_args(), str(out))
    # Annotated polyA for strand '+' is exons[-1][1]+1 = 401, so reads at 402
    # land within ANNOTATION_TOLERANCE (10).
    exons = [(100, 200), (300, 400)]
    for _ in range(30):
        counter.add_read_info(_make_read_assignment(polya_pos=402, exons=exons))
    counter.dump()

    df = _read_tsv(out)
    assert (df["flag"] == "Known").all()


def test_dump_grouped_emits_group_id_column(stub_model, tmp_path):
    # Minimal stub for string_pools: an object that returns a pool with a
    # ``get_str`` method translating integer ids to names.
    pool = types.SimpleNamespace(get_str=lambda i: f"g{i}")
    string_pools = types.SimpleNamespace(get_read_group_pool=lambda _idx: pool)

    out = tmp_path / "grouped.tsv"
    counter = tc.PolyACounter(_make_args(), str(out),
                              string_pools=string_pools, group_index=0)
    exons = [(100, 200), (300, 400)]
    for offset in range(30):
        counter.add_read_info(_make_read_assignment(
            polya_pos=420, exons=exons, read_group_ids=[offset % 2]))
    counter.dump()

    df = _read_tsv(out)
    assert list(df.columns) == tc.EMPTY_COLUMNS_GROUPED
    assert set(df["group_id"]) == {"g0", "g1"}
    # Total counts_byGroup should match the 30 contributing reads.
    assert df["counts_byGroup"].sum() == 30


def test_model_rejection_yields_empty_output(monkeypatch, tmp_path):
    # When the model predicts 0 for every peak, only zero-peak fallback rows
    # could survive. A 30-read cluster produces a peak, so output is empty.
    monkeypatch.setattr(tc, "XGBClassifier", lambda *a, **kw: _StubModel(accept=False))
    out = tmp_path / "rejected.tsv"
    counter = tc.PolyACounter(_make_args(), str(out))
    exons = [(100, 200), (300, 400)]
    for offset in range(30):
        counter.add_read_info(_make_read_assignment(
            polya_pos=420 + offset % 3, exons=exons))
    counter.dump()

    df = _read_tsv(out)
    assert df.empty


def test_rank_peaks_orders_by_height_descending():
    counter_cls = tc.TerminalCounter
    row = pd.Series({
        "peak_heights": np.array([3.0, 8.0, 5.0]),
        "peak_location": [10, 20, 30],
        "peak_prominence": [0.3, 0.8, 0.5],
        "peak_width": [1.0, 2.0, 3.0],
        "peak_left": [9, 19, 29],
        "peak_right": [11, 21, 31],
    })
    ranked = counter_cls._rank_peaks(None, row)
    assert ranked["peak_location"] == [20, 30, 10]
    assert ranked["rank"] == [1, 2, 3]
    assert ranked["relative_height"][0] == 1.0
    assert ranked["relative_height"][1] == pytest.approx(5.0 / 8.0)


def test_rank_peaks_handles_single_peak():
    counter_cls = tc.TerminalCounter
    row = pd.Series({
        "peak_heights": np.array([7.0]),
        "peak_location": [42],
        "peak_prominence": [0.7],
        "peak_width": [1.5],
        "peak_left": [40],
        "peak_right": [44],
    })
    ranked = counter_cls._rank_peaks(None, row)
    assert ranked["rank"] == [0]
    assert ranked["relative_height"] == [1.0]


def test_counts_for_peak_clamps_to_histogram_bounds():
    counter_cls = tc.TerminalCounter
    histogram = ([0] * tc.HISTOGRAM_PAD + [1, 2, 3, 4, 5] +
                 [0] * tc.HISTOGRAM_PAD)
    row = pd.Series({
        "histogram": histogram,
        "peak_left": -2,   # before histogram start, must clamp to 0
        "peak_right": 10,  # past histogram end, must clamp
    })
    assert counter_cls._counts_for_peak(None, row) == sum([1, 2, 3, 4, 5])


@pytest.mark.parametrize("strand,expected", [("+", 401), ("-", 99)])
def test_polya_counter_annotated_position(stub_model, tmp_path, strand, expected):
    counter = tc.PolyACounter(_make_args(), str(tmp_path / "p.tsv"))
    read = _make_read_assignment(strand=strand,
                                 exons=[(100, 200), (300, 400)])
    assert counter._annotated_position(read, "T1") == expected


@pytest.mark.parametrize("strand,expected", [("+", 99), ("-", 401)])
def test_tss_counter_annotated_position(stub_model, tmp_path, strand, expected):
    counter = tc.TSSCounter(_make_args(), str(tmp_path / "t.tsv"))
    read = _make_read_assignment(strand=strand,
                                 exons=[(100, 200), (300, 400)])
    assert counter._annotated_position(read, "T1") == expected


def test_tss_counter_rejects_dot_strand(stub_model, tmp_path):
    counter = tc.TSSCounter(_make_args(), str(tmp_path / "t.tsv"))
    counter.add_read_info(_make_read_assignment(strand="."))
    assert counter.transcripts == {}


def test_tss_counter_rejects_empty_exons(stub_model, tmp_path):
    counter = tc.TSSCounter(_make_args(), str(tmp_path / "t.tsv"))
    read = _make_read_assignment()
    read.corrected_exons = []
    counter.add_read_info(read)
    assert counter.transcripts == {}


def test_polya_counter_rejects_negative_strand_without_polyt(stub_model, tmp_path):
    counter = tc.PolyACounter(_make_args(), str(tmp_path / "p.tsv"))
    bad = _make_read_assignment(strand="-", polya_pos=420)
    bad.polya_info = _make_polya_info(external_polyt_pos=-1)
    counter.add_read_info(bad)
    assert counter.transcripts == {}


def test_dump_emits_rows_for_multiple_transcripts(stub_model, tmp_path):
    out = tmp_path / "multi.tsv"
    counter = tc.PolyACounter(_make_args(), str(out))
    exons_a = [(100, 200), (300, 400)]
    exons_b = [(1000, 1100), (1300, 1500)]
    gene_info = types.SimpleNamespace(
        all_isoforms_exons={"T1": exons_a, "T2": exons_b})
    for _ in range(30):
        counter.add_read_info(_make_read_assignment(
            transcript_id="T1", gene_id="G1",
            exons=exons_a, polya_pos=420, gene_info=gene_info))
        counter.add_read_info(_make_read_assignment(
            transcript_id="T2", gene_id="G2",
            exons=exons_b, polya_pos=1520, gene_info=gene_info))
    counter.dump()

    df = _read_tsv(out)
    assert set(df["transcript_id"]) == {"T1", "T2"}
    assert set(df["gene_id"]) == {"G1", "G2"}


def test_dump_zero_peak_fallback_uses_histogram_mode(stub_model, tmp_path):
    # A single read produces a one-bin histogram; find_peaks returns nothing,
    # so the zero-peak branch fires and emits a row using the mode position.
    out = tmp_path / "mode.tsv"
    counter = tc.PolyACounter(_make_args(), str(out))
    counter.add_read_info(_make_read_assignment(
        polya_pos=420, exons=[(100, 200), (300, 400)]))
    counter.dump()

    df = _read_tsv(out)
    assert len(df) == 1
    assert df.iloc[0]["prediction"] == 420


def test_no_op_methods_do_not_raise(stub_model, tmp_path):
    counter = tc.PolyACounter(_make_args(), str(tmp_path / "p.tsv"))
    counter.add_read_info(None)
    counter.add_read_info_raw("rid", ["F1"], [0])
    counter.add_confirmed_features({"F1"})
    counter.add_unassigned(_make_read_assignment())
    counter.add_unaligned(3)
    counter.finalize()
    counter.finalize(_make_args())
    assert counter.transcripts == {}


def test_init_truncates_existing_output_file(stub_model, tmp_path):
    path = tmp_path / "stale.tsv"
    path.write_text("garbage from a previous run\n")
    tc.PolyACounter(_make_args(), str(path))
    assert path.read_text() == ""
