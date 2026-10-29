############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Unit tests for known spliced isoforms built from read assignments (non-FL).

A known isoform without a full-length graph path (e.g. no trusted read end because
polyA tails were trimmed) is still built from its uniquely assigned reads when they
reach both annotated ends. The fl_pacbio and reliable presets used to skip this step,
so such isoforms vanished from the models however many reads they had.
"""

import types

import pytest

import isoquant
from isoquant_lib.assignment.isoform_assignment import MatchEventSubtype, ReadAssignmentType
from isoquant_lib.model_construction.assignment_based_constructor import AssignmentBasedConstructor
from isoquant_lib.model_construction.context import ModelStore

PRESETS = ["reliable", "default_pacbio", "sensitive_pacbio", "pacbio_all", "default_ont", "sensitive_ont",
           "fl_pacbio", "all", "assembly"]
EXONS = [(100, 199), (301, 400), (501, 600)]
INTRONS = ((200, 300), (401, 500))


def _preset_args(strategy: str) -> types.SimpleNamespace:
    args = types.SimpleNamespace(model_construction_strategy=strategy, graph_clustering_distance=None,
                                 report_novel_unspliced=False, no_model_construction=False,
                                 polya_requirement="auto", polya_trimmed="none", report_canonical="auto")
    isoquant.set_model_construction_options(args)
    args.simple_alignments_mapq_cutoff = 1
    args.requires_polya_for_construction = False
    args.apa_delta = 50
    return args


def _unique_read(read_id: str, exons: list) -> types.SimpleNamespace:
    match = types.SimpleNamespace(assigned_transcript="T1",
                                  match_subclassifications=[types.SimpleNamespace(event_type=MatchEventSubtype.fsm)])
    return types.SimpleNamespace(read_id=read_id, corrected_exons=exons, multimapper=False, mapping_quality=60,
                                 polyA_found=False, assignment_type=ReadAssignmentType.unique,
                                 isoform_matches=[match])


def _constructor(args) -> AssignmentBasedConstructor:
    gene_info = types.SimpleNamespace(chr_id="chr1", all_isoforms_exons={"T1": EXONS},
                                      all_isoforms_introns={"T1": list(INTRONS)}, isoform_strands={"T1": "+"},
                                      gene_id_map={"T1": "G1"}, sources={"T1": "test"}, other_features={"T1": []})
    ctx = types.SimpleNamespace(store=ModelStore(gene_info), args=args, gene_info=gene_info,
                                known_isoforms_in_graph_ids={"T1": INTRONS})
    return AssignmentBasedConstructor(ctx)


@pytest.fixture(autouse=True)
def fresh_detected_known_isoforms(monkeypatch):
    # process-wide set shared by all stores
    monkeypatch.setattr(ModelStore, "detected_known_isoforms", set())


@pytest.mark.parametrize("strategy", PRESETS)
def test_known_isoform_without_fl_path_is_built_from_assignments(strategy):
    constructor = _constructor(_preset_args(strategy))
    # 5'- and 3'-shortened reads together cover both annotated ends
    reads = [_unique_read("r1", [(110, 199), (301, 400), (501, 560)]),
             _unique_read("r2", [(140, 199), (301, 400), (501, 595)]),
             _unique_read("r3", [(105, 199), (301, 400), (501, 520)])]
    constructor.construct_assignment_based_isoforms(reads)

    assert [m.transcript_id for m in constructor.store.transcript_model_storage] == ["T1"]
    assert constructor.store.internal_counter["T1"] == 3


@pytest.mark.parametrize("strategy", PRESETS)
def test_known_isoform_is_not_built_when_no_read_reaches_an_annotated_end(strategy):
    constructor = _constructor(_preset_args(strategy))
    reads = [_unique_read("r%d" % i, [(110, 199), (301, 400), (501, 520)]) for i in range(5)]
    constructor.construct_assignment_based_isoforms(reads)

    assert constructor.store.transcript_model_storage == []
