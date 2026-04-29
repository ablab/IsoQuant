############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Convert IsoQuant's :class:`IntronGraph` into an "algorithmic" flow network:
integer-labelled consecutive vertices plus an edge list and an edge -> flow
dict. By default the network does not include a super-source or
super-target — opt in via ``Intron2Graph(..., add_super_source_target=True)``
when the downstream tool needs a single-source / single-sink wrap.

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

    Vertex ids start at ``0`` and are consecutive. By default no super-source
    or super-target is added; set ``add_super_source_target=True`` to wrap the
    graph with a single super-source (id 0) and super-target (last id).
    """

    def __init__(self, intron_graph, add_super_source_target: bool = False) -> None:
        self.intron2vertex: Dict[Tuple[int, int], int] = {}
        self.vertex2intron: Dict[int, Tuple[int, int]] = {}
        self.add_super_source_target: bool = add_super_source_target
        self.source: Optional[int] = None
        self.target: Optional[int] = None
        self.edge_list: List[Tuple[int, int]] = []
        self.flow_dict: Dict[Tuple[int, int], int] = defaultdict(int)

        next_id = 0
        if add_super_source_target:
            self.source = next_id
            next_id += 1

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

        if add_super_source_target:
            self.target = next_id
            next_id += 1

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

        if add_super_source_target:
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
EDGE_PRESENT = "-"
EDGE_MISSING = "|"


def _edge_in_graph(intron_graph, u_tuple, v_tuple) -> bool:
    """True iff ``u_tuple -> v_tuple`` exists in the IntronGraph.

    Starting-vertex edges only live in ``incoming_edges[intron]`` (they are
    never recorded in ``outgoing_edges[starting_vertex]``), so consulting
    ``outgoing_edges`` alone misses them.
    """
    if v_tuple in intron_graph.outgoing_edges.get(u_tuple, ()):
        return True
    if u_tuple in intron_graph.incoming_edges.get(v_tuple, ()):
        return True
    return False


def _terminal_candidates(
    flow: "Intron2Graph",
) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """Return ``(low_candidates, high_candidates)`` as sorted ``(pos, vid)`` lists.

    ``low_candidates`` covers polyT / read_start vertices (incoming side of an
    intron); ``high_candidates`` covers polyA / read_end (outgoing side).
    """
    low: List[Tuple[int, int]] = []
    high: List[Tuple[int, int]] = []
    for v_tuple, vid in flow.intron2vertex.items():
        if is_starting_vertex(v_tuple):
            low.append((v_tuple[1], vid))
        elif is_terminal_vertex(v_tuple):
            high.append((v_tuple[1], vid))
    low.sort()
    high.sort()
    return low, high


def _match_terminal(
    candidates: List[Tuple[int, int]],
    position: int,
    apa_delta: int,
) -> Optional[int]:
    """Return the closest candidate vertex id within ``apa_delta`` of ``position``.

    Linear scan; graph-local terminal lists are short. Ties broken by first
    occurrence in the sorted list (lowest position wins).
    """
    best_vid: Optional[int] = None
    best_diff = apa_delta + 1
    for cand_pos, vid in candidates:
        diff = abs(cand_pos - position)
        if diff <= apa_delta and diff < best_diff:
            best_diff = diff
            best_vid = vid
    return best_vid


def _render_path(
    flow: "Intron2Graph",
    intron_graph,
    mapped: List[Optional[int]],
) -> Tuple[str, str, str, int, int]:
    """Render a vertex-id chain into ``(status, path, path_simple, mv, me)``.

    ``mapped`` is a list of flow vertex ids with ``None`` for missing slots.
    ``status`` is ``ok`` (all vertices present, all edges present),
    ``disconnected`` (vertices all present, at least one edge missing) or
    ``partial`` (at least one vertex missing). ``path`` interleaves vertex
    tokens with ``-`` (edge present), ``|`` (both real but no edge), or
    ``-`` (default whenever either side is ``*``); ``path_simple`` is a
    comma-separated list of resolved vertex ids only.
    """
    missing_vertices = sum(1 for v in mapped if v is None)

    edge_ok: List[bool] = []
    for i in range(len(mapped) - 1):
        u, v = mapped[i], mapped[i + 1]
        if u is None or v is None:
            edge_ok.append(False)
            continue
        u_tuple = flow.vertex2intron[u]
        v_tuple = flow.vertex2intron[v]
        edge_ok.append(_edge_in_graph(intron_graph, u_tuple, v_tuple))
    missing_edges = sum(1 for ok in edge_ok if not ok)

    vertex_tokens = [GAP_MARKER if v is None else str(v) for v in mapped]
    simple_path_str = ",".join(str(v) for v in mapped if v is not None)

    parts: List[str] = []
    for i, token in enumerate(vertex_tokens):
        if i > 0:
            prev_real = mapped[i - 1] is not None
            cur_real = mapped[i] is not None
            if prev_real and cur_real:
                parts.append(EDGE_PRESENT if edge_ok[i - 1] else EDGE_MISSING)
            else:
                parts.append(EDGE_PRESENT)
        parts.append(token)
    path_str = "".join(parts)

    if missing_vertices > 0:
        status = "partial"
    elif missing_edges > 0:
        status = "disconnected"
    else:
        status = "ok"
    return status, path_str, simple_path_str, missing_vertices, missing_edges


def _thread_transcript(
    flow: "Intron2Graph",
    intron_graph,
    gene_info,
    t_id: str,
    introns,
    apa_delta: int,
    low_candidates: List[Tuple[int, int]],
    high_candidates: List[Tuple[int, int]],
) -> Tuple[str, str, str, int, int]:
    """Project an annotated transcript onto the simplified graph.

    Each annotated intron runs through ``intron_collector.substitute`` and
    maps to a graph vertex; the transcript's 5'/3' exon boundaries (low/high
    genomic ends) are matched to the closest starting / terminal vertex in
    the graph within ``apa_delta``. Anything that doesn't resolve becomes a
    ``*`` slot — no synthetic vertices are introduced.

    Returns ``(status, path_str, simple_path_str, missing_vertices,
    missing_edges)`` — see :func:`_render_path`.
    """
    if not introns:
        return "monoexonic", "", "", 0, 0

    intron_collector = intron_graph.intron_collector
    intron_mapped: List[Optional[int]] = []
    for intron in introns:
        if intron in intron_collector.discarded_introns:
            intron_mapped.append(None)
            continue
        substituted = intron_collector.substitute(intron)
        if substituted not in intron_collector.clustered_introns:
            intron_mapped.append(None)
            continue
        intron_mapped.append(flow.intron2vertex.get(substituted))

    exons = gene_info.all_isoforms_exons[t_id]
    low_vid = _match_terminal(low_candidates, exons[0][0], apa_delta)
    high_vid = _match_terminal(high_candidates, exons[-1][1], apa_delta)

    mapped: List[Optional[int]] = [low_vid] + intron_mapped + [high_vid]
    return _render_path(flow, intron_graph, mapped)


def _dump_paths(
    flow: "Intron2Graph",
    intron_graph,
    gene_info,
    counts_map: Dict[str, float],
    chr_id: str,
    paths_path: str,
    coverage_scale_factor: int = 1,
) -> None:
    """Write per-gene ground-truth paths TSV.

    Columns: ``transcript_id  count  count_scaled  status  path
    path_simple  missing_vertices  missing_edges``. ``count`` is the
    original (pre-downsampling) abundance from
    ``--ground_truth_counts``; ``count_scaled = round(count /
    coverage_scale_factor)`` is the expected post-downsampling load
    that matches the observed weights in ``vertices.tsv`` /
    ``read_subpaths.tsv`` when IsoQuant is processing 1 read out of
    every ``coverage_scale_factor`` (driven by
    ``--max_coverage_small_chr`` / ``--max_coverage_normal_chr``).
    With no downsampling (``scale == 1``) the two columns are equal.
    """
    apa_delta = intron_graph.params.apa_delta
    low_candidates, high_candidates = _terminal_candidates(flow)
    scale = max(1, int(coverage_scale_factor))

    rows = []
    for t_id, introns in gene_info.all_isoforms_introns.items():
        if t_id not in counts_map:
            continue
        count = counts_map[t_id]
        scaled_count = int(round(count / scale)) if scale > 1 else count
        status, path_str, simple_str, missing_vertices, missing_edges = _thread_transcript(
            flow, intron_graph, gene_info, t_id, introns,
            apa_delta, low_candidates, high_candidates,
        )
        rows.append((t_id, count, scaled_count, status, path_str, simple_str, missing_vertices, missing_edges))

    if not rows:
        return

    with open(paths_path, "w") as pf:
        pf.write("transcript_id\tcount\tcount_scaled\tstatus\tpath\tpath_simple\tmissing_vertices\tmissing_edges\n")
        for t_id, count, scaled_count, status, path_str, simple_str, mv, me in rows:
            pf.write("%s\t%g\t%g\t%s\t%s\t%s\t%d\t%d\n" %
                     (t_id, count, scaled_count, status, path_str, simple_str, mv, me))


def _dump_read_subpaths(
    flow: "Intron2Graph",
    intron_graph,
    paths: Dict[Tuple, int],
    fl_paths,
    out_path: str,
) -> None:
    """Write ``read_subpaths.tsv`` — one row per distinct threaded read path.

    Columns: ``read_count  is_fl  status  path  path_simple
    missing_vertices  missing_edges``. ``read_count`` is the number of
    reads that produced the path (from ``IntronPathStorage.paths``);
    ``is_fl`` is 1 when both a starting and a terminal vertex were
    threaded (the FL subset); otherwise 0. Path / path_simple use the
    same connectivity encoding as ``paths.tsv``.
    """
    if not paths:
        return

    rows = []
    for path_tuple, count in paths.items():
        mapped: List[Optional[int]] = [flow.intron2vertex.get(v) for v in path_tuple]
        status, path_str, simple_str, mv, me = _render_path(flow, intron_graph, mapped)
        is_fl = 1 if path_tuple in fl_paths else 0
        rows.append((count, is_fl, status, path_str, simple_str, mv, me))

    rows.sort(key=lambda r: (-r[0], -r[1]))

    with open(out_path, "w") as f:
        f.write("read_count\tis_fl\tstatus\tpath\tpath_simple\tmissing_vertices\tmissing_edges\n")
        for count, is_fl, status, path_str, simple_str, mv, me in rows:
            f.write("%d\t%d\t%s\t%s\t%s\t%d\t%d\n" %
                    (count, is_fl, status, path_str, simple_str, mv, me))


def _dump_ref_data(
    flow: "Intron2Graph",
    intron_graph,
    gene_info,
    ref_vertices_path: str,
    ref_edges_path: str,
) -> None:
    """Emit per-gene reference-vs-graph diff TSVs.

    ``<chr>.<gene>.ref_vertices.tsv`` — one row per unique annotated slot
    across all transcripts in this gene:
    ``kind  ref_start  ref_end  status  vertex_id  graph_start  graph_end``
    - ``kind`` ∈ {``intron``, ``starting``, ``terminal``} — ``starting`` is
      the low-coord exon boundary (matched against polyT / read_start),
      ``terminal`` is the high-coord boundary (polyA / read_end)
    - ``status`` for introns ∈ {``in_graph``, ``discarded``, ``unmapped``};
      for terminal slots ∈ {``in_graph``, ``unmapped``}
    - ``vertex_id`` is the integer graph id (``*`` when unresolved)
    - ``graph_start/graph_end`` are the post-substitution intron coords or
      the matched terminal position; ``*`` when no graph counterpart

    Terminal slots are matched by position within
    ``intron_graph.params.apa_delta``; if no graph terminal sits inside the
    radius the slot is recorded as ``unmapped`` rather than synthesized.

    ``<chr>.<gene>.ref_edges.tsv`` — one row per unique consecutive pair
    appearing in at least one transcript (intron→intron, plus
    starting→first_intron and last_intron→terminal):
    ``kind  u_start  u_end  v_start  v_end  status  u_id  v_id``
    - ``status`` ∈ {``in_graph``, ``missing_edge``, ``missing_vertex``}
    """
    intron_collector = intron_graph.intron_collector
    apa_delta = intron_graph.params.apa_delta
    low_candidates, high_candidates = _terminal_candidates(flow)

    vertex_info: Dict[
        Tuple[str, int, int],
        Tuple[Optional[int], str, Optional[Tuple[int, int]]],
    ] = {}
    edge_info: Dict[
        Tuple, Tuple[Optional[int], Optional[int], str]
    ] = {}

    def resolve_intron(intron: Tuple[int, int]) -> Optional[int]:
        key = ("intron", intron[0], intron[1])
        if key in vertex_info:
            return vertex_info[key][0]
        if intron in intron_collector.discarded_introns:
            vertex_info[key] = (None, "discarded", None)
            return None
        substituted = intron_collector.substitute(intron)
        if substituted not in intron_collector.clustered_introns:
            vertex_info[key] = (None, "unmapped", substituted)
            return None
        vid = flow.intron2vertex.get(substituted)
        if vid is None:
            vertex_info[key] = (None, "unmapped", substituted)
            return None
        vertex_info[key] = (vid, "in_graph", substituted)
        return vid

    def resolve_terminal(position: int, side: str) -> Optional[int]:
        kind = "starting" if side == "low" else "terminal"
        key = (kind, position, position)
        if key in vertex_info:
            return vertex_info[key][0]
        candidates = low_candidates if side == "low" else high_candidates
        vid = _match_terminal(candidates, position, apa_delta)
        if vid is None:
            vertex_info[key] = (None, "unmapped", None)
            return None
        graph_tuple = flow.vertex2intron[vid]
        vertex_info[key] = (vid, "in_graph", (graph_tuple[1], graph_tuple[1]))
        return vid

    def classify_edge(
        u_vid: Optional[int], v_vid: Optional[int]
    ) -> Tuple[Optional[int], Optional[int], str]:
        if u_vid is None or v_vid is None:
            return u_vid, v_vid, "missing_vertex"
        u_tuple = flow.vertex2intron[u_vid]
        v_tuple = flow.vertex2intron[v_vid]
        if _edge_in_graph(intron_graph, u_tuple, v_tuple):
            return u_vid, v_vid, "in_graph"
        return u_vid, v_vid, "missing_edge"

    for t_id, introns in gene_info.all_isoforms_introns.items():
        if not introns:
            continue
        mapped = [resolve_intron(intron) for intron in introns]
        for i in range(len(introns) - 1):
            key = ("intron",
                   introns[i][0], introns[i][1],
                   introns[i + 1][0], introns[i + 1][1])
            if key in edge_info:
                continue
            edge_info[key] = classify_edge(mapped[i], mapped[i + 1])

        exons = gene_info.all_isoforms_exons[t_id]
        low_pos = exons[0][0]
        high_pos = exons[-1][1]
        low_vid = resolve_terminal(low_pos, "low")
        high_vid = resolve_terminal(high_pos, "high")

        first_intron = introns[0]
        last_intron = introns[-1]
        start_key = ("starting", low_pos, low_pos,
                     first_intron[0], first_intron[1])
        if start_key not in edge_info:
            edge_info[start_key] = classify_edge(low_vid, mapped[0])
        end_key = ("terminal", last_intron[0], last_intron[1],
                   high_pos, high_pos)
        if end_key not in edge_info:
            edge_info[end_key] = classify_edge(mapped[-1], high_vid)

    with open(ref_vertices_path, "w") as vf:
        vf.write("kind\tref_start\tref_end\tstatus\tvertex_id\tgraph_start\tgraph_end\n")
        for (kind, rs, re_), (vid, status, graph_tuple) in sorted(vertex_info.items()):
            vid_str = GAP_MARKER if vid is None else str(vid)
            if graph_tuple is None:
                gs_str = ge_str = GAP_MARKER
            else:
                gs_str, ge_str = str(graph_tuple[0]), str(graph_tuple[1])
            vf.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\n"
                     % (kind, rs, re_, status, vid_str, gs_str, ge_str))

    with open(ref_edges_path, "w") as ef:
        ef.write("kind\tu_start\tu_end\tv_start\tv_end\tstatus\tu_id\tv_id\n")
        for key, (u_id, v_id, status) in sorted(edge_info.items()):
            kind = key[0]
            us, ue, vs, ve_ = key[1], key[2], key[3], key[4]
            u_str = GAP_MARKER if u_id is None else str(u_id)
            v_str = GAP_MARKER if v_id is None else str(v_id)
            ef.write("%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n"
                     % (kind, us, ue, vs, ve_, status, u_str, v_str))


def dump_flow_graph(
    intron_graph,
    chr_id: str,
    gene_id: str,
    out_dir: str,
    gene_info=None,
    ground_truth_counts: Optional[Dict[str, float]] = None,
    dump_ref_data: bool = False,
    path_storage=None,
    add_super_source_target: bool = False,
    include_edge_weights: bool = False,
) -> None:
    """Write the gene's flow network to ``<out_dir>/<chr>.<gene>/`` as TSVs.

    - ``vertices.tsv``: ``vertex_id  type  chr  start  end  weight``
      (intron rows carry the clustered-intron count; source / target /
      terminal / starting rows mark missing fields with ``*``)
    - ``edges.tsv``: ``u  v`` (default) or ``u  v  weight`` when
      ``include_edge_weights=True``
    - ``paths.tsv`` (only when ``ground_truth_counts`` and ``gene_info`` are
      given and at least one transcript in the gene appears in the counts
      map): ``transcript_id  count  status  path  path_simple
      missing_vertices  missing_edges``
    - ``ref_vertices.tsv`` and ``ref_edges.tsv`` (only when
      ``dump_ref_data`` is true and ``gene_info`` carries annotated
      transcripts): every unique annotated intron / consecutive intron
      pair with its graph status (see ``_dump_ref_data`` for column
      layouts).
    - ``read_subpaths.tsv`` (only when ``path_storage`` is provided):
      every distinct read-threaded path with its supporting read count
      and FL flag (see ``_dump_read_subpaths``).

    ``add_super_source_target`` and ``include_edge_weights`` are internal
    knobs (no CLI exposure) for downstream tooling that wants the wrapped
    flow network or per-edge weights.
    """
    if not intron_graph.intron_collector.clustered_introns \
            and not intron_graph.outgoing_edges \
            and not intron_graph.incoming_edges:
        return

    flow = Intron2Graph(intron_graph, add_super_source_target=add_super_source_target)
    clustered_introns = intron_graph.intron_collector.clustered_introns

    safe_gene = gene_id.replace("/", "_")
    gene_dir = os.path.join(out_dir, "%s.%s" % (chr_id, safe_gene))
    os.makedirs(gene_dir, exist_ok=True)

    vertices_path = os.path.join(gene_dir, "vertices.tsv")
    edges_path = os.path.join(gene_dir, "edges.tsv")

    incoming_totals: Dict[int, int] = defaultdict(int)
    outgoing_totals: Dict[int, int] = defaultdict(int)
    for (u, v), w in flow.flow_dict.items():
        incoming_totals[v] += w
        outgoing_totals[u] += w

    def _vertex_weight(intron_tuple: Tuple[int, int]) -> str:
        if intron_tuple in clustered_introns:
            return str(clustered_introns[intron_tuple])
        vid = flow.intron2vertex.get(intron_tuple)
        if vid is not None:
            if is_terminal_vertex(intron_tuple):
                return str(incoming_totals[vid])
            if is_starting_vertex(intron_tuple):
                return str(outgoing_totals[vid])
        return GAP_MARKER

    with open(vertices_path, "w") as vf:
        vf.write("vertex_id\ttype\tchr\tstart\tend\tweight\n")
        if flow.source is not None:
            vf.write("%d\tsource\t%s\t%s\t%s\t%s\n" % (
                flow.source, chr_id, GAP_MARKER, GAP_MARKER, GAP_MARKER))
        for vid in sorted(flow.vertex2intron.keys()):
            intron_tuple = flow.vertex2intron[vid]
            vtype, start, end = _vertex_info(intron_tuple)
            vf.write("%d\t%s\t%s\t%d\t%d\t%s\n" % (
                vid, vtype, chr_id, start, end, _vertex_weight(intron_tuple)))
        if flow.target is not None:
            vf.write("%d\ttarget\t%s\t%s\t%s\t%s\n" % (
                flow.target, chr_id, GAP_MARKER, GAP_MARKER, GAP_MARKER))

    with open(edges_path, "w") as ef:
        if include_edge_weights:
            ef.write("u\tv\tweight\n")
            for u, v in flow.edge_list:
                ef.write("%d\t%d\t%d\n" % (u, v, flow.flow_dict[(u, v)]))
        else:
            ef.write("u\tv\n")
            for u, v in flow.edge_list:
                ef.write("%d\t%d\n" % (u, v))

    if ground_truth_counts and gene_info is not None:
        paths_path = os.path.join(gene_dir, "paths.tsv")
        scale = getattr(gene_info, "coverage_scale_factor", 1)
        _dump_paths(flow, intron_graph, gene_info, ground_truth_counts, chr_id, paths_path,
                    coverage_scale_factor=scale)

    if dump_ref_data and gene_info is not None and gene_info.all_isoforms_introns:
        ref_vertices_path = os.path.join(gene_dir, "ref_vertices.tsv")
        ref_edges_path = os.path.join(gene_dir, "ref_edges.tsv")
        _dump_ref_data(flow, intron_graph, gene_info, ref_vertices_path, ref_edges_path)

    if path_storage is not None and getattr(path_storage, "paths", None):
        read_subpaths_path = os.path.join(gene_dir, "read_subpaths.tsv")
        _dump_read_subpaths(flow, intron_graph,
                            path_storage.paths,
                            getattr(path_storage, "fl_paths", set()),
                            read_subpaths_path)

    logger.debug("Dumped flow graph for %s / %s (%d vertices, %d edges) to %s",
                 chr_id, gene_id, len(flow.vertex2intron), len(flow.edge_list), gene_dir)
