############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Unit tests for intron-graph terminal-vertex refinement
(``IntronGraph._refine_positions`` / ``_attach_side``) on branch
``transcript_model_ends``.

The graph is heavy to build, so the methods under test run on a bare instance
created with ``__new__`` plus the few attributes/stubs they touch.
"""

import types
from collections import defaultdict

from isoquant_lib.intron_graph import IntronGraph, VERTEX_read_end


def _bare_graph(apa_delta=10, polya_predictions=None, tss_predictions=None):
    g = IntronGraph.__new__(IntronGraph)
    g.params = types.SimpleNamespace(
        apa_delta=apa_delta,
        terminal_position_abs=1,
        terminal_position_rel=0.0,
        terminal_internal_position_rel=0.0,
    )
    g.polya_predictions = polya_predictions
    g.tss_predictions = tss_predictions
    g.outgoing_edges = defaultdict(set)
    g.incoming_edges = defaultdict(set)
    return g


# -- _refine_positions --------------------------------------------------------

def test_refine_positions_snaps_within_delta():
    g = _bare_graph(apa_delta=10)
    assert g._refine_positions({105: 3}, [100, 500]) == {100: 3}


def test_refine_positions_leaves_when_far():
    g = _bare_graph(apa_delta=10)
    assert g._refine_positions({130: 3}, [100, 500]) == {130: 3}


def test_refine_positions_merges_counts_on_collision():
    g = _bare_graph(apa_delta=10)
    assert g._refine_positions({98: 2, 103: 5}, [100]) == {100: 7}


def test_refine_positions_identity_without_predictions():
    g = _bare_graph(apa_delta=10)
    assert g._refine_positions({105: 3}, None) == {105: 3}
    assert g._refine_positions({105: 3}, []) == {105: 3}


# -- _attach_side side-selection (regression for review fix #1) ----------------

def test_attach_side_3prime_readend_refines_with_polya_not_tss():
    # read_end=True is the genomic 3' side: a VERTEX_read_end position must be
    # refined toward the polyA predictions, never the (5') TSS predictions.
    # Place the read-end cluster at 1005, equidistant from polya (1000) and tss
    # (1010); the chosen target tells us which set was used.
    intron = (100, 200)
    g = _bare_graph(apa_delta=10, polya_predictions=[1000], tss_predictions=[1010])
    g.intron_collector = types.SimpleNamespace(clustered_introns={intron: 10})
    # No polyA-confirmed vertices; one read-end cluster at 1005.
    g.cluster_polya_positions = lambda positions, i, read_end: {}
    g.cluster_terminal_positions = lambda extra, read_end, cutoff: {1005: 4}

    polya_confirmed = {intron: {}}
    read_terminal = {intron: {1005: 4}}
    g._attach_side([intron], polya_confirmed, read_terminal, read_end=True)

    vertices = g.outgoing_edges[intron]
    assert (VERTEX_read_end, 1000) in vertices    # snapped to polyA prediction
    assert (VERTEX_read_end, 1010) not in vertices  # NOT the TSS prediction
    assert (VERTEX_read_end, 1005) not in vertices  # and it was refined, not left raw
