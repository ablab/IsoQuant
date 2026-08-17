############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Unit tests for ``IntronPathProcessor.thread_ends`` / ``thread_starts``
terminal-vertex selection.

Regression cover for issue #415. A confirmed-TSS vertex claims the reads that
terminate at it (proximity, mirroring the trusted-polyA lookup) before the
extreme / second-last comparison runs. It also takes part in that comparison,
which is safe because ``_attach_side`` only keeps bare read termini that lie
beyond every confirmed site on the same side -- so a bare vertex is always the
outermost candidate, never an inner truncation cluster.
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

def test_tss_vertex_claims_reads_that_end_at_it():
    # The proximity claim runs first, so a read terminating at the TSS reaches it
    # even though the bare read-end vertex is the outermost candidate. Without it
    # the extreme rule would hand this read the 1300 vertex instead.
    intron = (100, 200)
    p = _processor(intron, outgoing=[(TerminalVertex.tss_right, 1000),
                                     (TerminalVertex.read_end, 1300)])
    assert p.thread_ends(intron, 1010, trusted=False) == (TerminalVertex.tss_right, 1000)


def test_bare_read_end_beyond_a_tss_vertex_is_reachable():
    # The arrangement _attach_side now guarantees: the bare cluster lies beyond
    # every confirmed site, so its own reads reach it through the extreme rule.
    intron = (100, 200)
    p = _processor(intron, outgoing=[(TerminalVertex.tss_right, 1000),
                                     (TerminalVertex.read_end, 1300)])
    assert p.thread_ends(intron, 1290, trusted=False) == (TerminalVertex.read_end, 1300)


def test_read_ending_before_an_inner_read_end_cluster_is_still_rejected():
    # The second-last rule itself is preserved for genuine read-end / polyA
    # candidates: a read stopping at the inner terminus is ambiguous.
    intron = (100, 200)
    p = _processor(intron, outgoing=[(TerminalVertex.read_end, 1000),
                                     (TerminalVertex.polya, 900)])
    assert p.thread_ends(intron, 850, trusted=False) is None


# -- thread_starts ------------------------------------------------------------

def test_tss_left_vertex_claims_reads_that_start_at_it():
    intron = (1000, 2000)
    p = _processor(intron, incoming=[(TerminalVertex.tss_left, 800),
                                     (TerminalVertex.read_start, 500)])
    assert p.thread_starts(intron, 790, trusted=False) == (TerminalVertex.tss_left, 800)


def test_bare_read_start_beyond_a_tss_vertex_is_reachable():
    intron = (1000, 2000)
    p = _processor(intron, incoming=[(TerminalVertex.tss_left, 800),
                                     (TerminalVertex.read_start, 500)])
    assert p.thread_starts(intron, 510, trusted=False) == (TerminalVertex.read_start, 500)
