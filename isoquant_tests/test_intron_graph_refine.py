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

from isoquant_lib.model_construction.intron_graph import IntronGraph, TerminalVertex


def _bare_graph(apa_delta=10, polya_predictions=None, tss_predictions=None,
                transcript_introns=None, gene_id_map=None):
    g = IntronGraph.__new__(IntronGraph)
    g.params = types.SimpleNamespace(
        apa_delta=apa_delta,
        terminal_position_abs=1,
        terminal_position_rel=0.0,
        terminal_internal_position_rel=0.0,
    )
    # Predictions are (transcript_id, position); admissibility resolves the
    # transcript to its gene and that gene's terminal introns.
    g.polya_predictions = polya_predictions
    g.tss_predictions = tss_predictions
    g.gene_info = types.SimpleNamespace(
        all_isoforms_introns=transcript_introns or {},
        gene_id_map=gene_id_map or {},
    )
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

def _one_gene_graph(intron, apa_delta=10, polya_predictions=None, tss_predictions=None):
    # Single gene G whose only isoform ends (both sides) at `intron`, so every
    # prediction of G is admissible on it.
    g = _bare_graph(apa_delta=apa_delta,
                    polya_predictions=polya_predictions,
                    tss_predictions=tss_predictions,
                    transcript_introns={"T": [intron]},
                    gene_id_map={"T": "G"})
    g.intron_collector = types.SimpleNamespace(clustered_introns={intron: 10},
                                               substitute=lambda v: v)
    return g


def test_attach_side_3prime_readend_refines_with_polya_not_tss():
    # read_end=True is the genomic 3' side: a TerminalVertex.read_end position must be
    # refined toward the polyA predictions, never the (5') TSS predictions.
    # Place the read-end cluster at 1005, equidistant from polya (1000) and tss
    # (1010); the chosen target tells us which set was used.
    intron = (100, 200)
    g = _one_gene_graph(intron, apa_delta=10,
                        polya_predictions=[("T", 1000)], tss_predictions=[("T", 1010)])
    # No polyA-confirmed vertices; one read-end cluster at 1005.
    g.cluster_polya_positions = lambda positions, i, read_end: {}
    g.cluster_terminal_positions = lambda extra, read_end, cutoff: {1005: 4}

    polya_confirmed = {intron: {}}
    read_terminal = {intron: {1005: 4}}
    g._attach_side([intron], polya_confirmed, read_terminal, read_end=True)

    vertices = g.outgoing_edges[intron]
    assert (TerminalVertex.read_end, 1000) in vertices    # snapped to polyA prediction
    assert (TerminalVertex.read_end, 1010) not in vertices  # NOT the TSS prediction
    assert (TerminalVertex.read_end, 1005) not in vertices  # and it was refined, not left raw


# -- prediction admissibility (regression for issue #415) ----------------------

def test_attach_side_uses_own_gene_tss_prediction():
    # A TSS prediction of this intron's own gene creates a tss_right vertex for
    # the read termini that sit at it.
    intron = (100, 200)
    g = _one_gene_graph(intron, apa_delta=10, tss_predictions=[("T", 1000)])
    g.cluster_polya_positions = lambda positions, i, read_end: {}
    g.cluster_terminal_positions = lambda extra, read_end, cutoff: {}

    g._attach_side([intron], {intron: {}}, {intron: {1005: 2}}, read_end=True)
    assert (TerminalVertex.tss_right, 1000) in g.outgoing_edges[intron]


def test_attach_side_ignores_neighbouring_gene_tss_prediction():
    # Issue #415: a neighbouring gene's TSS (here RAG1's, next to TRAF6) must not
    # spawn a terminal vertex on this gene's intron -- it would sit past the
    # dominant read-end cluster and disqualify all of its reads in thread_ends.
    intron = (100, 200)
    g = _bare_graph(apa_delta=10, tss_predictions=[("OTHER_T", 1000)],
                    transcript_introns={"T": [intron], "OTHER_T": [(5000, 6000)]},
                    gene_id_map={"T": "G", "OTHER_T": "OTHER_G"})
    g.intron_collector = types.SimpleNamespace(clustered_introns={intron: 10},
                                               substitute=lambda v: v)
    g.cluster_polya_positions = lambda positions, i, read_end: {}
    g.cluster_terminal_positions = lambda extra, read_end, cutoff: {1005: 2}

    g._attach_side([intron], {intron: {}}, {intron: {1005: 2}}, read_end=True)
    vertices = g.outgoing_edges[intron]
    assert not any(v[0] == TerminalVertex.tss_right for v in vertices)
    # the reads stay in the ordinary clustering path instead
    assert (TerminalVertex.read_end, 1005) in vertices


def test_attach_side_accepts_sibling_isoform_tss_prediction():
    # Predictions are computed per reference isoform but a TSS is usually shared
    # by the gene's other isoforms, whose terminal introns differ -- a sibling's
    # terminal intron must still be able to use it.
    intron, sibling_intron = (100, 200), (100, 250)
    g = _bare_graph(apa_delta=10, tss_predictions=[("T2", 1000)],
                    transcript_introns={"T1": [intron], "T2": [sibling_intron]},
                    gene_id_map={"T1": "G", "T2": "G"})
    g.intron_collector = types.SimpleNamespace(clustered_introns={intron: 10},
                                               substitute=lambda v: v)
    g.cluster_polya_positions = lambda positions, i, read_end: {}
    g.cluster_terminal_positions = lambda extra, read_end, cutoff: {}

    g._attach_side([intron], {intron: {}}, {intron: {1005: 2}}, read_end=True)
    assert (TerminalVertex.tss_right, 1000) in g.outgoing_edges[intron]


def test_attach_side_ignores_prediction_inside_the_intron():
    # A predicted position must lie beyond the intron, inside the terminal exon.
    intron = (100, 2000)
    g = _one_gene_graph(intron, apa_delta=10, tss_predictions=[("T", 1000)])
    g.cluster_polya_positions = lambda positions, i, read_end: {}
    g.cluster_terminal_positions = lambda extra, read_end, cutoff: {2005: 2}

    g._attach_side([intron], {intron: {}}, {intron: {2005: 2}}, read_end=True)
    assert not any(v[0] == TerminalVertex.tss_right for v in g.outgoing_edges[intron])
