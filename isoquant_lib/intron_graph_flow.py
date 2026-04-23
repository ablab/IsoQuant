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


GAP_MARKER = "*"


def _thread_transcript(
    flow: "Intron2Graph",
    intron_graph,
    introns,
) -> Tuple[str, str, int, int]:
    """Project an annotated transcript's intron tuples onto the simplified graph.

    Each annotated intron is run through ``intron_collector.substitute`` and
    mapped to a graph vertex; introns in ``discarded_introns`` or whose
    substituted form is not in ``clustered_introns`` map to ``None``. Every
    consecutive pair is then checked against ``intron_graph.outgoing_edges``.

    Returns ``(status, path_str, missing_vertices, missing_edges)``:
    - ``status`` is one of ``monoexonic`` / ``ok`` / ``disconnected`` /
      ``partial`` (partial = at least one vertex missing).
    - ``path_str`` is a comma-separated rendering where surviving vertex ids
      are joined normally and any break (missing vertex or disconnected edge
      between two surviving vertices) is marked by a single ``*`` token.
    """
    if not introns:
        return "monoexonic", "", 0, 0

    intron_collector = intron_graph.intron_collector
    mapped: List[Optional[int]] = []
    for intron in introns:
        if intron in intron_collector.discarded_introns:
            mapped.append(None)
            continue
        substituted = intron_collector.substitute(intron)
        if substituted not in intron_collector.clustered_introns:
            mapped.append(None)
            continue
        mapped.append(flow.intron2vertex.get(substituted))

    missing_vertices = sum(1 for v in mapped if v is None)

    edge_ok: List[bool] = []
    for i in range(len(mapped) - 1):
        u, v = mapped[i], mapped[i + 1]
        if u is None or v is None:
            edge_ok.append(False)
            continue
        u_tuple = flow.vertex2intron[u]
        v_tuple = flow.vertex2intron[v]
        edge_ok.append(v_tuple in intron_graph.outgoing_edges.get(u_tuple, ()))
    missing_edges = sum(1 for ok in edge_ok if not ok)

    tokens: List[str] = []

    def push_gap() -> None:
        if not tokens or tokens[-1] != GAP_MARKER:
            tokens.append(GAP_MARKER)

    for i, v in enumerate(mapped):
        if v is None:
            push_gap()
        else:
            if i > 0 and mapped[i - 1] is not None and not edge_ok[i - 1]:
                push_gap()
            tokens.append(str(v))

    path_str = ",".join(tokens)

    if missing_vertices > 0:
        status = "partial"
    elif missing_edges > 0:
        status = "disconnected"
    else:
        status = "ok"
    return status, path_str, missing_vertices, missing_edges


def _dump_paths(
    flow: "Intron2Graph",
    intron_graph,
    gene_info,
    counts_map: Dict[str, float],
    chr_id: str,
    paths_path: str,
) -> None:
    """Write per-gene ground-truth paths TSV.

    Columns: ``transcript_id  count  status  path  missing_vertices
    missing_edges``. ``path`` is comma-separated; surviving vertex ids are
    rendered as ints and any break (missing vertex or disconnected edge)
    collapses into a single ``*`` token.
    """
    rows = []
    for t_id, introns in gene_info.all_isoforms_introns.items():
        if t_id not in counts_map:
            continue
        count = counts_map[t_id]
        status, path_str, missing_vertices, missing_edges = _thread_transcript(
            flow, intron_graph, introns
        )
        rows.append((t_id, count, status, path_str, missing_vertices, missing_edges))

    if not rows:
        return

    with open(paths_path, "w") as pf:
        pf.write("transcript_id\tcount\tstatus\tpath\tmissing_vertices\tmissing_edges\n")
        for t_id, count, status, path_str, mv, me in rows:
            pf.write("%s\t%g\t%s\t%s\t%d\t%d\n" % (t_id, count, status, path_str, mv, me))


def _dump_ref_data(
    flow: "Intron2Graph",
    intron_graph,
    gene_info,
    ref_vertices_path: str,
    ref_edges_path: str,
) -> None:
    """Emit per-gene reference-vs-graph diff TSVs.

    ``<chr>.<gene>.ref_vertices.tsv`` — one row per unique annotated intron
    across all transcripts in this gene:
    ``ref_start  ref_end  status  vertex_id  graph_start  graph_end``
    - ``status`` ∈ {``in_graph``, ``discarded``, ``unmapped``}
    - ``vertex_id`` is the integer graph id (empty for non-``in_graph``)
    - ``graph_start/graph_end`` are the post-substitution coords, repeating
      the ref coords when no substitution occurred (empty for ``discarded``)

    ``<chr>.<gene>.ref_edges.tsv`` — one row per unique consecutive intron
    pair appearing in at least one transcript:
    ``u_start  u_end  v_start  v_end  status  u_id  v_id``
    - ``status`` ∈ {``in_graph``, ``missing_edge``, ``missing_vertex``}
    """
    intron_collector = intron_graph.intron_collector

    vertex_info: Dict[Tuple[int, int], Tuple[Optional[int], str, Optional[Tuple[int, int]]]] = {}
    edge_info: Dict[Tuple[Tuple[int, int], Tuple[int, int]], Tuple[Optional[int], Optional[int], str]] = {}

    def resolve(intron: Tuple[int, int]) -> Optional[int]:
        if intron in vertex_info:
            return vertex_info[intron][0]
        if intron in intron_collector.discarded_introns:
            vertex_info[intron] = (None, "discarded", None)
            return None
        substituted = intron_collector.substitute(intron)
        if substituted not in intron_collector.clustered_introns:
            vertex_info[intron] = (None, "unmapped", substituted)
            return None
        vid = flow.intron2vertex.get(substituted)
        if vid is None:
            vertex_info[intron] = (None, "unmapped", substituted)
            return None
        vertex_info[intron] = (vid, "in_graph", substituted)
        return vid

    for _t_id, introns in gene_info.all_isoforms_introns.items():
        if not introns:
            continue
        mapped = [resolve(intron) for intron in introns]
        for i in range(len(introns) - 1):
            key = (introns[i], introns[i + 1])
            if key in edge_info:
                continue
            u, v = mapped[i], mapped[i + 1]
            if u is None or v is None:
                edge_info[key] = (u, v, "missing_vertex")
                continue
            u_tuple = flow.vertex2intron[u]
            v_tuple = flow.vertex2intron[v]
            if v_tuple in intron_graph.outgoing_edges.get(u_tuple, ()):
                edge_info[key] = (u, v, "in_graph")
            else:
                edge_info[key] = (u, v, "missing_edge")

    with open(ref_vertices_path, "w") as vf:
        vf.write("ref_start\tref_end\tstatus\tvertex_id\tgraph_start\tgraph_end\n")
        for (rs, re_), (vid, status, graph_tuple) in sorted(vertex_info.items()):
            vid_str = "" if vid is None else str(vid)
            if graph_tuple is None:
                gs_str = ge_str = ""
            else:
                gs_str, ge_str = str(graph_tuple[0]), str(graph_tuple[1])
            vf.write("%d\t%d\t%s\t%s\t%s\t%s\n" % (rs, re_, status, vid_str, gs_str, ge_str))

    with open(ref_edges_path, "w") as ef:
        ef.write("u_start\tu_end\tv_start\tv_end\tstatus\tu_id\tv_id\n")
        for ((us, ue), (vs, ve_)), (u_id, v_id, status) in sorted(edge_info.items()):
            u_str = "" if u_id is None else str(u_id)
            v_str = "" if v_id is None else str(v_id)
            ef.write("%d\t%d\t%d\t%d\t%s\t%s\t%s\n" % (us, ue, vs, ve_, status, u_str, v_str))


def dump_flow_graph(
    intron_graph,
    chr_id: str,
    gene_id: str,
    out_dir: str,
    gene_info=None,
    ground_truth_counts: Optional[Dict[str, float]] = None,
    dump_ref_data: bool = False,
) -> None:
    """Write the gene's flow network to ``out_dir`` as TSVs.

    - ``<chr>.<gene>.vertices.tsv``: ``vertex_id  type  chr  start  end``
      (source / target rows have empty start/end columns)
    - ``<chr>.<gene>.edges.tsv``: ``u  v  weight``
    - ``<chr>.<gene>.paths.tsv`` (only when ``ground_truth_counts`` and
      ``gene_info`` are given and at least one transcript in the gene
      appears in the counts map): ``transcript_id  count  status  path
      missing_vertices  missing_edges``
    - ``<chr>.<gene>.ref_vertices.tsv`` and ``<chr>.<gene>.ref_edges.tsv``
      (only when ``dump_ref_data`` is true and ``gene_info`` carries
      annotated transcripts): every unique annotated intron / consecutive
      intron pair with its graph status (see ``_dump_ref_data`` for
      column layouts).
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

    if dump_ref_data and gene_info is not None and gene_info.all_isoforms_introns:
        ref_vertices_path = os.path.join(out_dir, "%s.%s.ref_vertices.tsv" % (chr_id, safe_gene))
        ref_edges_path = os.path.join(out_dir, "%s.%s.ref_edges.tsv" % (chr_id, safe_gene))
        _dump_ref_data(flow, intron_graph, gene_info, ref_vertices_path, ref_edges_path)

    logger.debug("Dumped flow graph for %s / %s (%d vertices, %d edges) to %s",
                 chr_id, gene_id, flow.target + 1, len(flow.edge_list), out_dir)
