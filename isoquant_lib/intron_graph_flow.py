############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Convert IsoQuant's :class:`IntronGraph` into an "algorithmic" flow network:
integer-labelled vertices with a super-source (id 0) and a super-target
(id N), plus an edge list and an edge -> flow dict. Intended for external
flow-decomposition / ILP tooling that expects consecutive int vertices.

Ported from the ``src/encode_ilp_gurobi.py::Intron2Graph`` prototype on the
``ilp_models`` branch.
"""

import logging
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .intron_graph import (
    VERTEX_polya,
    VERTEX_polyt,
    VERTEX_read_end,
    VERTEX_read_start,
    is_starting_vertex,
    is_terminal_vertex,
)

logger = logging.getLogger('IsoQuant')


_VERTEX_TYPE = {
    VERTEX_polya: "polya",
    VERTEX_read_end: "read_end",
    VERTEX_polyt: "polyt",
    VERTEX_read_start: "read_start",
}


def _vertex_info(v: Tuple[int, int]) -> Tuple[str, int, int]:
    """Return ``(type, start, end)`` for an IntronGraph vertex."""
    if v[0] in _VERTEX_TYPE:
        pos = v[1]
        return _VERTEX_TYPE[v[0]], pos, pos
    return "intron", v[0], v[1]


class Intron2Graph:
    """Integer-labelled flow network built from an :class:`IntronGraph`.

    Vertex 0 is the super-source, vertex ``self.target`` is the super-target;
    all intron, starting and terminal vertices carry ids in ``[1, target)``.
    """

    def __init__(self, intron_graph) -> None:
        self.intron2vertex: Dict[Tuple[int, int], int] = {}
        self.vertex2intron: Dict[int, Tuple[int, int]] = {}
        self.source: int = 0
        self.edge_list: List[Tuple[int, int]] = []
        self.flow_dict: Dict[Tuple[int, int], int] = defaultdict(int)

        next_id = 1

        # 1. real intron vertices
        for intron in intron_graph.intron_collector.clustered_introns.keys():
            self.intron2vertex[intron] = next_id
            self.vertex2intron[next_id] = intron
            next_id += 1

        # 2. starting vertices (polyT / read_start) referenced from the graph
        for intron, preds in intron_graph.incoming_edges.items():
            for preceding in preds:
                if is_starting_vertex(preceding) and preceding not in self.intron2vertex:
                    self.intron2vertex[preceding] = next_id
                    self.vertex2intron[next_id] = preceding
                    next_id += 1

        # 3. terminal vertices (polyA / read_end) referenced from the graph
        for intron, succs in intron_graph.outgoing_edges.items():
            for subsequent in succs:
                if is_terminal_vertex(subsequent) and subsequent not in self.intron2vertex:
                    self.intron2vertex[subsequent] = next_id
                    self.vertex2intron[next_id] = subsequent
                    next_id += 1

        self.target: int = next_id

        # 4. intron -> intron / starting -> intron / intron -> terminal edges
        edge_set = set()
        starting_totals: Dict[Tuple[int, int], int] = defaultdict(int)
        terminal_totals: Dict[Tuple[int, int], int] = defaultdict(int)

        for intron, preds in intron_graph.incoming_edges.items():
            for preceding in preds:
                u = self.intron2vertex[preceding]
                v = self.intron2vertex[intron]
                if (u, v) not in edge_set:
                    self.edge_list.append((u, v))
                    edge_set.add((u, v))
                self.flow_dict[(u, v)] = intron_graph.edge_weights[(preceding, intron)]
                if is_starting_vertex(preceding):
                    starting_totals[preceding] += intron_graph.edge_weights[(preceding, intron)]

        for intron, succs in intron_graph.outgoing_edges.items():
            for subsequent in succs:
                if not is_terminal_vertex(subsequent):
                    continue
                u = self.intron2vertex[intron]
                v = self.intron2vertex[subsequent]
                if (u, v) not in edge_set:
                    self.edge_list.append((u, v))
                    edge_set.add((u, v))
                self.flow_dict[(u, v)] = intron_graph.edge_weights[(intron, subsequent)]
                terminal_totals[subsequent] += intron_graph.edge_weights[(intron, subsequent)]

        # 5. super-source -> each starting vertex, weight = sum of outgoing flow
        for start_v, total in starting_totals.items():
            sv = self.intron2vertex[start_v]
            self.edge_list.append((self.source, sv))
            self.flow_dict[(self.source, sv)] = total

        # 6. each terminal vertex -> super-target
        for end_v, total in terminal_totals.items():
            tv = self.intron2vertex[end_v]
            self.edge_list.append((tv, self.target))
            self.flow_dict[(tv, self.target)] = total

    def transcript_to_path(self, transcript) -> List[Tuple[int, int]]:
        vertex_list = [self.intron2vertex[t] for t in transcript]
        return list(zip(vertex_list[:-1], vertex_list[1:]))

    def path_to_transcript(self, path) -> List[Tuple[int, int]]:
        return [self.vertex2intron[v] for v in path]


def _thread_transcript(intron_collector, introns) -> Tuple[str, List[Tuple[int, int]]]:
    """Project an annotated transcript's intron tuples onto the simplified graph.

    Mirrors ``IntronPathProcessor.thread_introns``: discarded introns abort,
    surviving introns are substituted via ``intron_correction_map``. Returns
    ``(status, path)`` where ``path`` is empty unless ``status == "ok"``.
    """
    if not introns:
        return "monoexonic", []
    path: List[Tuple[int, int]] = []
    for intron in introns:
        if intron in intron_collector.discarded_introns:
            return "discarded", []
        substituted = intron_collector.substitute(intron)
        if substituted not in intron_collector.clustered_introns:
            return "unmapped", []
        path.append(substituted)
    return "ok", path


def _dump_paths(
    flow: "Intron2Graph",
    intron_graph,
    gene_info,
    counts_map: Dict[str, float],
    chr_id: str,
    paths_path: str,
) -> None:
    """Write per-gene ground-truth paths TSV.

    Columns: ``transcript_id  count  status  path`` (path is comma-separated
    int vertex ids). Emits one row for every transcript in ``gene_info`` that
    appears in ``counts_map``.
    """
    rows = []
    for t_id, introns in gene_info.all_isoforms_introns.items():
        if t_id not in counts_map:
            continue
        count = counts_map[t_id]
        status, tuple_path = _thread_transcript(intron_graph.intron_collector, introns)
        if status == "ok":
            try:
                vertex_path = [flow.intron2vertex[t] for t in tuple_path]
                path_str = ",".join(str(v) for v in vertex_path)
            except KeyError:
                # Substituted intron slipped through; treat as unmapped.
                status = "unmapped"
                path_str = ""
        else:
            path_str = ""
        rows.append((t_id, count, status, path_str))

    if not rows:
        return

    with open(paths_path, "w") as pf:
        pf.write("transcript_id\tcount\tstatus\tpath\n")
        for t_id, count, status, path_str in rows:
            pf.write("%s\t%g\t%s\t%s\n" % (t_id, count, status, path_str))


def dump_flow_graph(
    intron_graph,
    chr_id: str,
    gene_id: str,
    out_dir: str,
    gene_info=None,
    ground_truth_counts: Optional[Dict[str, float]] = None,
) -> None:
    """Write the gene's flow network to ``out_dir`` as TSVs.

    - ``<chr>.<gene>.vertices.tsv``: ``vertex_id  type  chr  start  end``
      (source / target rows have empty start/end columns)
    - ``<chr>.<gene>.edges.tsv``: ``u  v  weight``
    - ``<chr>.<gene>.paths.tsv`` (only when ``ground_truth_counts`` and
      ``gene_info`` are given and at least one transcript in the gene
      appears in the counts map): ``transcript_id  count  status  path``
    """
    if not intron_graph.intron_collector.clustered_introns \
            and not intron_graph.outgoing_edges \
            and not intron_graph.incoming_edges:
        return

    os.makedirs(out_dir, exist_ok=True)
    flow = Intron2Graph(intron_graph)

    safe_gene = gene_id.replace("/", "_")
    vertices_path = os.path.join(out_dir, "%s.%s.vertices.tsv" % (chr_id, safe_gene))
    edges_path = os.path.join(out_dir, "%s.%s.edges.tsv" % (chr_id, safe_gene))

    with open(vertices_path, "w") as vf:
        vf.write("vertex_id\ttype\tchr\tstart\tend\n")
        vf.write("%d\tsource\t%s\t\t\n" % (flow.source, chr_id))
        for vid in sorted(flow.vertex2intron.keys()):
            vtype, start, end = _vertex_info(flow.vertex2intron[vid])
            vf.write("%d\t%s\t%s\t%d\t%d\n" % (vid, vtype, chr_id, start, end))
        vf.write("%d\ttarget\t%s\t\t\n" % (flow.target, chr_id))

    with open(edges_path, "w") as ef:
        ef.write("u\tv\tweight\n")
        for u, v in flow.edge_list:
            ef.write("%d\t%d\t%d\n" % (u, v, flow.flow_dict[(u, v)]))

    if ground_truth_counts and gene_info is not None:
        paths_path = os.path.join(out_dir, "%s.%s.paths.tsv" % (chr_id, safe_gene))
        _dump_paths(flow, intron_graph, gene_info, ground_truth_counts, chr_id, paths_path)

    logger.debug("Dumped flow graph for %s / %s (%d vertices, %d edges) to %s",
                 chr_id, gene_id, flow.target + 1, len(flow.edge_list), out_dir)
