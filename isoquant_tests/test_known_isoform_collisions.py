############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Unit tests for annotated isoforms that share (or nearly share) an intron chain.

Two annotated isoforms whose chains differ by a few bases at one splice site
(NAGNAG acceptors) used to be merged in the intron graph, after which both
threaded to one graph path, and the path -> isoform map kept only the isoform
listed last. The other one was never built, even with >100k supporting reads.
The same map collision hit isoforms with identical chains and different ends.

The graph and context are heavy to build, so the methods under test run on bare
instances created with ``__new__`` plus the few attributes they touch.
"""

import types
from collections import defaultdict

from isoquant_lib.assignment.isoform_assignment import MatchEventSubtype, ReadAssignmentType
from isoquant_lib.gene_info import TranscriptModel, TranscriptModelType
from isoquant_lib.model_construction import model_filter
from isoquant_lib.model_construction.context import ModelConstructionContext
from isoquant_lib.model_construction.fl_graph_constructor import GraphBasedConstructor
from isoquant_lib.model_construction.intron_graph import IntronGraph, TerminalVertex
from isoquant_lib.model_construction.model_filter import ModelFilter

# Fcgr3: ENSMUST00000164044.8 and ENSMUST00000169017.8 differ only by a 3 bp acceptor shift
MAJOR_INTRON = (170885518, 170886154)
MINOR_INTRON = (170885515, 170886154)


def _bare_graph(counts, known_introns):
    g = IntronGraph.__new__(IntronGraph)
    g.params = types.SimpleNamespace(graph_clustering_distance=10, graph_clustering_ratio=0.5)
    g.intron_collector = types.SimpleNamespace(clustered_introns=dict(counts), known_introns=set(known_introns))
    return g


def _bare_context(isoform_introns=None, isoform_exons=None, substitutions=None):
    ctx = ModelConstructionContext.__new__(ModelConstructionContext)
    ctx.gene_info = types.SimpleNamespace(all_isoforms_introns=isoform_introns or {},
                                          all_isoforms_exons=isoform_exons or {})
    substitutions = substitutions or {}
    ctx.path_processor = types.SimpleNamespace(
        thread_introns=lambda introns: [substitutions.get(i, i) for i in introns])
    return ctx


def _read(isoform_id, assignment_type=ReadAssignmentType.unique):
    return types.SimpleNamespace(assignment_type=assignment_type,
                                 isoform_matches=[types.SimpleNamespace(assigned_transcript=isoform_id)])


# -- IntronGraph.collapse_vertex_set ------------------------------------------

def test_annotated_introns_a_few_bases_apart_are_not_collapsed():
    graph = _bare_graph({MAJOR_INTRON: 420, MINOR_INTRON: 26}, known_introns=[MAJOR_INTRON, MINOR_INTRON])
    assert graph.collapse_vertex_set({MAJOR_INTRON, MINOR_INTRON}) == {}


def test_novel_minor_intron_is_still_collapsed_into_annotated_one():
    graph = _bare_graph({MAJOR_INTRON: 420, MINOR_INTRON: 26}, known_introns=[MAJOR_INTRON])
    assert graph.collapse_vertex_set({MAJOR_INTRON, MINOR_INTRON}) == {MINOR_INTRON: MAJOR_INTRON}


def test_novel_minor_intron_is_still_collapsed_into_novel_one():
    graph = _bare_graph({MAJOR_INTRON: 420, MINOR_INTRON: 26}, known_introns=[])
    assert graph.collapse_vertex_set({MAJOR_INTRON, MINOR_INTRON}) == {MINOR_INTRON: MAJOR_INTRON}


# -- ModelConstructionContext.get_known_spliced_isoforms ----------------------

def test_isoforms_sharing_a_path_are_all_kept():
    chain = [(100, 200), (300, 400)]
    ctx = _bare_context(isoform_introns={"T1": chain, "T2": list(chain), "T3": [(100, 200)]})
    known = ctx.get_known_spliced_isoforms(ctx.gene_info)
    assert known == {tuple(chain): ["T1", "T2"], ((100, 200),): ["T3"]}


def test_isoforms_merged_by_substitution_are_all_kept():
    ctx = _bare_context(isoform_introns={"major": [(10, 20), MAJOR_INTRON], "minor": [(10, 20), MINOR_INTRON]},
                        substitutions={MINOR_INTRON: MAJOR_INTRON})
    known = ctx.get_known_spliced_isoforms(ctx.gene_info)
    assert known == {((10, 20), MAJOR_INTRON): ["major", "minor"]}


# -- ModelConstructionContext.select_known_isoform ----------------------------

EXONS = {"T1": [(100, 199), (401, 500)], "T2": [(80, 199), (401, 520)]}


def test_single_candidate_is_returned():
    ctx = _bare_context(isoform_exons=EXONS)
    assert ctx.select_known_isoform(["T2"], [], (100, 500)) == "T2"


def test_path_assignment_decides_when_it_names_one_candidate():
    ctx = _bare_context(isoform_exons=EXONS)
    reads = [_read("T2"), _read("T2")]
    assert ctx.select_known_isoform(["T1", "T2"], reads, (100, 500), path_assignment=_read("T1")) == "T1"


def test_read_majority_decides_for_an_ambiguous_path():
    ctx = _bare_context(isoform_exons=EXONS)
    ambiguous_path = types.SimpleNamespace(isoform_matches=[types.SimpleNamespace(assigned_transcript="T1"),
                                                            types.SimpleNamespace(assigned_transcript="T2")])
    reads = [_read("T2"), _read("T2"), _read("T1"), _read("T1", ReadAssignmentType.ambiguous)]
    assert ctx.select_known_isoform(["T1", "T2"], reads, (100, 500), path_assignment=ambiguous_path) == "T2"


def test_closest_ends_break_a_tie():
    ctx = _bare_context(isoform_exons=EXONS)
    assert ctx.select_known_isoform(["T1", "T2"], [], (82, 518)) == "T2"
    assert ctx.select_known_isoform(["T1", "T2"], [], (100, 500)) == "T1"


# -- GraphBasedConstructor._reference_isoform_for_path -------------------------

def test_exact_intron_chain_beats_a_shifted_sibling_with_better_ends():
    # ENSMUST00000047179.7 vs ENSMUST00000234956.2: one donor 9 bp apart, the sibling starts 60 bp
    # upstream. A path with the first chain but the sibling's start was assigned to the sibling as
    # an fsm with an intron_shift, so the exactly matching isoform lost all its reads.
    exact = [(200, 300), (409, 500)]
    shifted = [(200, 300), (400, 500)]
    exons = {"exact": [(160, 199), (301, 408), (501, 600)], "shifted": [(100, 199), (301, 399), (501, 600)]}
    path = ((TerminalVertex.read_start, 100),) + tuple(exact) + ((TerminalVertex.read_end, 600),)

    ctx = _bare_context(isoform_exons=exons)
    ctx.known_isoforms_in_graph = {tuple(exact): ["exact"], tuple(shifted): ["shifted"]}
    ctx.path_storage = types.SimpleNamespace(paths_to_reads={path: [_read("exact")]})
    constructor = GraphBasedConstructor.__new__(GraphBasedConstructor)
    constructor.ctx = ctx
    constructor.gene_info = types.SimpleNamespace(all_isoforms_exons=exons, isoform_strands={})
    constructor.args = types.SimpleNamespace(apa_delta=50)
    constructor.use_tss_model = False

    events = [types.SimpleNamespace(event_type=e) for e in
              (MatchEventSubtype.fsm, MatchEventSubtype.intron_shift)]
    path_assignment = types.SimpleNamespace(
        assignment_type=ReadAssignmentType.unique_minor_difference,
        isoform_matches=[types.SimpleNamespace(assigned_transcript="shifted", match_subclassifications=events)])
    assert constructor._reference_isoform_for_path(path_assignment, path, tuple(exact), (100, 600)) == "exact"


# -- ModelFilter.detect_similar_isoforms ---------------------------------------
# Two annotated introns a few bases apart stay separate graph vertices, so reads split between them
# make two novel models that differ only at that splice site. The weaker one must not remove the
# stronger one just because its path came first.

def _novel_model(transcript_id, exons, transcript_type=TranscriptModelType.novel_in_catalog):
    model = TranscriptModel("chr18", "-", transcript_id, "gene", exons, transcript_type)
    model.intron_path = tuple((exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1))
    return model


def _model_filter(monkeypatch, read_counts, delta=6):
    # a single-model "gene" matches another model when exon counts agree and every boundary is within delta
    class Assigner:
        def __init__(self, model, *args, **kwargs):
            self.model = model

        def assign_to_isoform(self, read_id, exons):
            matching = len(exons) == len(self.model.exon_blocks) and \
                all(abs(x - y) <= delta for e1, e2 in zip(exons, self.model.exon_blocks) for x, y in zip(e1, e2))
            return types.SimpleNamespace(matching=matching, assignment_type=ReadAssignmentType.inconsistent)

    class ProfileConstructor:
        def __init__(self, *args):
            pass

        def construct_profiles(self, exons, polya_info, polya_hints):
            return exons

    monkeypatch.setattr(model_filter, "GeneInfo", types.SimpleNamespace(from_models=lambda models, d: models[0]))
    monkeypatch.setattr(model_filter, "LongReadAssigner", Assigner)
    monkeypatch.setattr(model_filter, "CombinedProfileConstructor", ProfileConstructor)
    monkeypatch.setattr(model_filter, "is_matching_assignment", lambda assignment: assignment.matching)
    f = ModelFilter.__new__(ModelFilter)
    f.args = types.SimpleNamespace(delta=delta)
    f.string_pools = None
    f.use_tss_model = False
    f.store = types.SimpleNamespace(internal_counter=defaultdict(int, read_counts))
    return f


# Mouse.ONT_simulated.R9.trunc, ENSMUST00000035804.9: 40 reads align at acceptor 46853509, 6 at the
# annotated alternative 46853506; the 6-read model sorted first and removed the 37-read one
WEAK_EXONS = [(46846258, 46847028), (46848664, 46848833), (46853354, 46853505), (46856385, 46856462),
              (46861040, 46861465)]
STRONG_EXONS = [(46846258, 46847028), (46848664, 46848833), (46853354, 46853508), (46856385, 46856462),
                (46861040, 46861465)]


def test_weaker_near_identical_novel_model_is_removed(monkeypatch):
    f = _model_filter(monkeypatch, {"weak": 6, "strong": 37})
    models = [_novel_model("weak", WEAK_EXONS), _novel_model("strong", STRONG_EXONS)]
    assert f.detect_similar_isoforms(models) == {"weak": "strong"}


def test_known_model_absorbs_a_better_supported_novel_one(monkeypatch):
    f = _model_filter(monkeypatch, {"novel": 37, "known": 6})
    models = [_novel_model("novel", STRONG_EXONS),
              _novel_model("known", WEAK_EXONS, transcript_type=TranscriptModelType.known)]
    assert f.detect_similar_isoforms(models) == {"novel": "known"}
