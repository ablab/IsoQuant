############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Unit tests for ``IntronPathProcessor.thread_ends`` / ``thread_starts``
terminal-vertex selection.

Regression cover for issue #415: a confirmed-TSS vertex is attached without a
coverage cutoff, so it must claim only the reads that end at it and must not
enter the extreme / second-last comparison that decides whether the ordinary
read-end population is full-length.
"""

import types
from collections import defaultdict

from isoquant_lib.model_construction.intron_graph import TerminalVertex
from isoquant_lib.model_construction.intron_path import IntronPathProcessor


def _processor(intron, outgoing=(), incoming=(), apa_delta=50, delta=5):
    graph = types.SimpleNamespace(
        intron_collector=types.SimpleNamespace(clustered_introns={intron: 10}),
        outgoing_edges=defaultdict(set),
        incoming_edges=defaultdict(set),
    )
    graph.outgoing_edges[intron] = set(outgoing)
    graph.incoming_edges[intron] = set(incoming)
    graph.get_outgoing = lambda i, vertex_type=None: {
        v for v in graph.outgoing_edges[i]
        if (vertex_type is None and v[0] >= 0) or v[0] == vertex_type}
    graph.get_incoming = lambda i, vertex_type=None: {
        v for v in graph.incoming_edges[i]
        if (vertex_type is None and v[0] >= 0) or v[0] == vertex_type}
    params = types.SimpleNamespace(apa_delta=apa_delta, delta=delta)
    return IntronPathProcessor(params, graph)


# -- thread_ends --------------------------------------------------------------

def test_dominant_read_end_survives_a_more_extreme_tss_vertex():
    # Issue #415: the 835-read TRAF6 terminus at 36510313 with a 2-read TSS
    # vertex 57 bp further right. Reads of the dominant cluster must still get
    # their vertex instead of being rejected for not passing the TSS one.
    intron = (36501538, 36510047)
    p = _processor(intron, outgoing=[(TerminalVertex.read_end, 36510313),
                                     (TerminalVertex.tss_right, 36510370)])
    assert p.thread_ends(intron, 36510300, trusted=False) == (TerminalVertex.read_end, 36510313)
    assert p.thread_ends(intron, 36510313, trusted=False) == (TerminalVertex.read_end, 36510313)


def test_tss_vertex_claims_reads_that_end_at_it():
    intron = (36501538, 36510047)
    p = _processor(intron, outgoing=[(TerminalVertex.read_end, 36510313),
                                     (TerminalVertex.tss_right, 36510370)])
    assert p.thread_ends(intron, 36510385, trusted=False) == (TerminalVertex.tss_right, 36510370)


def test_read_ending_before_an_inner_read_end_cluster_is_still_rejected():
    # The second-last rule itself is preserved for genuine read-end / polyA
    # candidates: a read stopping at the inner terminus is ambiguous.
    intron = (100, 200)
    p = _processor(intron, outgoing=[(TerminalVertex.read_end, 1000),
                                     (TerminalVertex.polya, 900)])
    assert p.thread_ends(intron, 850, trusted=False) is None


# -- thread_starts ------------------------------------------------------------

def test_dominant_read_start_survives_a_more_extreme_tss_vertex():
    intron = (1000, 2000)
    p = _processor(intron, incoming=[(TerminalVertex.read_start, 800),
                                     (TerminalVertex.tss_left, 740)])
    assert p.thread_starts(intron, 810, trusted=False) == (TerminalVertex.read_start, 800)


def test_tss_left_vertex_claims_reads_that_start_at_it():
    intron = (1000, 2000)
    p = _processor(intron, incoming=[(TerminalVertex.read_start, 800),
                                     (TerminalVertex.tss_left, 740)])
    assert p.thread_starts(intron, 730, trusted=False) == (TerminalVertex.tss_left, 740)
