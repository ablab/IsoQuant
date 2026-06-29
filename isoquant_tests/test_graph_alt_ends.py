############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Unit tests for alternative-end NIC discovery in
``graph_based_model_construction`` (branch ``transcript_model_ends``).

The constructor is heavy to build (it wires up an assigner + profile
constructor), so the helpers under test are exercised on a bare instance
created with ``__new__`` and only the attributes each method touches.
``detect_peaks`` / ``get_polya_model`` / ``get_tss_model`` are stubbed so the
tests pin the gating/wiring logic, not scipy's peak finder or the shipped
XGBoost model.
"""

import types
from argparse import Namespace

import pytest

from isoquant_lib import graph_based_model_construction as gbmc
from isoquant_lib.gene_info import TranscriptModel, TranscriptModelType
from isoquant_lib.terminal_peaks import Peak

# Reuse the polyA-counter harness (stub model + read-assignment builder).
from isoquant_tests.test_polya_prediction import _StubModel, _make_read_assignment


class _IdDist:
    def __init__(self):
        self.n = 0

    def increment(self):
        self.n += 1
        return self.n


def _make_constructor(use_tss_model=True, apa_delta=10, min_novel_count=2,
                      terminal_position_rel=0.1, chr_id="chr1", all_isoforms_exons=None):
    c = gbmc.GraphBasedModelConstructor.__new__(gbmc.GraphBasedModelConstructor)
    c.args = Namespace(apa_delta=apa_delta, min_novel_count=min_novel_count,
                       terminal_position_rel=terminal_position_rel)
    c.use_tss_model = use_tss_model
    c.create_nics = True
    c.novel_apa = False
    c.id_distributor = _IdDist()
    c.gene_info = types.SimpleNamespace(chr_id=chr_id,
                                        all_isoforms_exons=all_isoforms_exons or {})
    c.internal_counter = {}
    return c


def _model(exon_blocks, strand="+", ttype=TranscriptModelType.known, tid="T1",
           gene_id="G1", chr_id="chr1"):
    return TranscriptModel(chr_id, strand, tid, gene_id,
                           [tuple(e) for e in exon_blocks], ttype)


@pytest.fixture
def stub_models(monkeypatch):
    model = _StubModel(accept=True)
    monkeypatch.setattr(gbmc, "get_polya_model", lambda: model)
    monkeypatch.setattr(gbmc, "get_tss_model", lambda: model)
    return model


# -- static helpers -----------------------------------------------------------

def test_intron_chain_key_multi_exon():
    m = _model([(100, 200), (300, 400), (500, 600)])
    assert gbmc.GraphBasedModelConstructor._intron_chain_key(m) == ((200, 300), (400, 500))


def test_intron_chain_key_monoexon_is_empty():
    assert gbmc.GraphBasedModelConstructor._intron_chain_key(_model([(100, 600)])) == ()


def test_closest_inward_greater_and_less():
    f = gbmc.GraphBasedModelConstructor._closest_inward
    assert f([10, 20, 30], 15, True) == 20      # first strictly greater (ascending)
    assert f([10, 20, 30], 30, True) is None
    assert f([30, 20, 10], 25, False) == 20     # first strictly less (descending)
    assert f([30, 20, 10], 5, False) is None


def test_confirmed_polya_pos_plus_and_minus():
    f = gbmc.GraphBasedModelConstructor._confirmed_polya_pos
    assert f(_make_read_assignment(strand="+", polya_pos=500), "+") == 500
    assert f(_make_read_assignment(strand="-", polya_pos=120), "-") == 120


def test_confirmed_polya_pos_requires_polya_found():
    a = _make_read_assignment(strand="+", polya_pos=500, polyA_found=False)
    assert gbmc.GraphBasedModelConstructor._confirmed_polya_pos(a, "+") is None


def test_confirmed_polya_pos_missing_strand_field_returns_none():
    # '+' read carries external_polya_pos but external_polyt_pos == -1
    a = _make_read_assignment(strand="+", polya_pos=500)
    assert gbmc.GraphBasedModelConstructor._confirmed_polya_pos(a, "-") is None


# -- _terminal_model ----------------------------------------------------------

def test_terminal_model_polya_side_always_available(stub_models):
    c = _make_constructor(use_tss_model=True)
    assert c._terminal_model("+", left=False) is stub_models   # 3' polyA of '+'
    assert c._terminal_model("-", left=True) is stub_models    # 3' polyA of '-'


def test_terminal_model_tss_side_needs_fl_data(stub_models):
    assert _make_constructor(use_tss_model=True)._terminal_model("+", left=True) is stub_models
    assert _make_constructor(use_tss_model=False)._terminal_model("+", left=True) is None


def test_terminal_model_dot_strand_is_none(stub_models):
    c = _make_constructor()
    assert c._terminal_model(".", left=True) is None
    assert c._terminal_model(".", left=False) is None


# -- _alternative_end_positions (relative-support + apa_delta gate) ------------

def test_alternative_end_positions_gate(monkeypatch):
    c = _make_constructor(apa_delta=10, min_novel_count=2, terminal_position_rel=0.1)
    peaks = [
        Peak(position=400, count=50, left=395, right=405),   # dominant, == annotated -> drop
        Peak(position=500, count=8, left=495, right=505),    # strong alt -> keep (>= cutoff 5)
        Peak(position=600, count=3, left=595, right=605),    # minor -> below cutoff -> drop
        Peak(position=405, count=40, left=400, right=410),   # within apa_delta of annotated -> drop
    ]
    monkeypatch.setattr(gbmc, "detect_peaks", lambda hist, model: peaks)
    out = c._alternative_end_positions({1: 1}, _StubModel(), annotated_pos=400, clamp=lambda p: True)
    assert out == [500]


def test_alternative_end_positions_respects_clamp(monkeypatch):
    c = _make_constructor(apa_delta=10, min_novel_count=2, terminal_position_rel=0.1)
    peaks = [Peak(400, 50, 395, 405), Peak(500, 8, 495, 505)]
    monkeypatch.setattr(gbmc, "detect_peaks", lambda hist, model: peaks)
    # clamp rejects the 500 peak -> nothing left
    out = c._alternative_end_positions({1: 1}, _StubModel(), annotated_pos=400, clamp=lambda p: p < 450)
    assert out == []


# -- derive_alternative_end_models + _nic_model_with_boundary ------------------

def _fake_detect(hist, model):
    # end histogram carries the alternative polyA at 500; start histogram carries
    # only the annotated start at 100 (filtered out by the apa_delta gate).
    if 500 in hist:
        return [Peak(500, 30, 495, 505)]
    if 100 in hist:
        return [Peak(100, 30, 95, 105)]
    return []


def test_derive_alternative_end_models_known_source(monkeypatch, stub_models):
    c = _make_constructor(use_tss_model=True, apa_delta=10, min_novel_count=2,
                          terminal_position_rel=0.1)
    monkeypatch.setattr(gbmc, "detect_peaks", _fake_detect)
    source = _model([(100, 200), (300, 400)], strand="+",
                    ttype=TranscriptModelType.known, tid="ENST1")
    reads = [_make_read_assignment(strand="+", polya_pos=500, exons=[(100, 200), (300, 500)])
             for _ in range(20)]

    out = c.derive_alternative_end_models(source, reads)

    assert len(out) == 1
    nic = out[0]
    assert nic.exon_blocks[0] == (100, 200)               # start unchanged
    assert nic.exon_blocks[-1] == (300, 500)              # 3' end moved to alt polyA
    assert nic.transcript_type == TranscriptModelType.novel_in_catalog
    assert nic.transcript_id.endswith(".nic")
    assert nic.gene_id == source.gene_id and nic.strand == "+"
    # source model is left untouched
    assert source.exon_blocks[-1] == (300, 400)


def test_derive_alternative_end_models_empty_reads():
    assert _make_constructor().derive_alternative_end_models(
        _model([(100, 200), (300, 400)]), []) == []


def test_nic_model_with_boundary_replaces_single_end():
    c = _make_constructor()
    source = _model([(100, 200), (300, 400), (500, 600)], strand="-", tid="ENST2", gene_id="G2")
    nic = c._nic_model_with_boundary(source, TranscriptModelType.novel_in_catalog, end=650)
    assert nic.exon_blocks == [(100, 200), (300, 400), (500, 650)]
    assert nic.transcript_id.endswith(".nic") and nic.strand == "-" and nic.gene_id == "G2"
    assert source.exon_blocks[-1] == (500, 600)            # source untouched


def test_nic_model_with_boundary_nnic_suffix_for_novel_source():
    c = _make_constructor()
    source = _model([(100, 200), (300, 400)], ttype=TranscriptModelType.novel_not_in_catalog)
    m = c._nic_model_with_boundary(source, TranscriptModelType.novel_not_in_catalog, start=50)
    assert m.exon_blocks[0] == (50, 200)
    assert m.transcript_id.endswith(".nnic")


# -- _drop_duplicate_alt_end_models -------------------------------------------

def test_drop_duplicate_alt_end_models():
    ref_exons = {"REF1": [(100, 200), (300, 400)]}
    c = _make_constructor(apa_delta=10, all_isoforms_exons=ref_exons)
    known = _model([(100, 200), (300, 400)], ttype=TranscriptModelType.known, tid="REF1")
    # same intron chain + both ends within apa_delta of REF1 -> it IS the known -> drop
    dup = _model([(105, 200), (300, 398)], ttype=TranscriptModelType.novel_in_catalog, tid="DUP")
    # genuinely different 3' end -> a real alt-end NIC -> keep
    keep = _model([(100, 200), (300, 800)], ttype=TranscriptModelType.novel_in_catalog, tid="KEEP")
    c.transcript_model_storage = [known, dup, keep]

    c._drop_duplicate_alt_end_models()

    ids = {m.transcript_id for m in c.transcript_model_storage}
    assert ids == {"REF1", "KEEP"}
