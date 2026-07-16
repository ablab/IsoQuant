############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict

from .common import junctions_from_blocks
from .assignment_io import ReadAssignmentType
from .gene_info import TranscriptModelType
from .isoform_assignment import (
    match_subtype_to_str_with_additional_info,
    ReadAssignment,
    MatchClassification,
    IsoformMatch,
)
from .polya_finder import PolyAInfo
from .transcript_to_gene_joiner import TranscriptToGeneJoiner
from .model_construction.context import ModelStore, ModelConstructionContext
# Re-exported for isoquant.py's `from ...graph_based_model_construction import StrandnessReportingLevel`.
from .model_construction.context import StrandnessReportingLevel  # noqa: F401
from .model_construction.fl_graph_constructor import FLGraphConstructor
from .model_construction.assignment_based_constructor import AssignmentBasedConstructor
from .model_construction.end_processor import TranscriptEndProcessor
from .model_construction.model_filter import ModelFilter
from .model_construction.read_assigner import ConstructionReadAssigner
from .model_construction.model_read_counter import ModelReadCounter

logger = logging.getLogger('IsoQuant')


class GraphBasedModelConstructor:

    def __init__(self, gene_info, chr_record, args, transcript_counter, gene_counter, id_distributor,
                 grouping_strategy_names=None,
                 use_technical_replicas=False,
                 string_pools=None,
                 polya_predictions=None,
                 tss_predictions=None):
        # Read-only inputs + intron graph + the shared ModelStore live on the ctx;
        # the extracted construction stages each receive this ctx.
        store = ModelStore(gene_info)
        self.ctx = ModelConstructionContext(
            gene_info, chr_record, args, id_distributor, string_pools,
            polya_predictions, tss_predictions, store,
            grouping_strategy_names=grouping_strategy_names,
            use_technical_replicas=use_technical_replicas)
        ctx = self.ctx
        self.store = store
        # Inputs used by this orchestrator itself (gene_info is also read by
        # parallel_workers; the rest drive process()/compare_models_with_known).
        self.gene_info = ctx.gene_info
        self.args = ctx.args
        self.string_pools = ctx.string_pools
        self.id_distributor = ctx.id_distributor
        self.assigner = ctx.assigner
        self.profile_constructor = ctx.profile_constructor
        self.transcript2transcript = []

        # Construction stages, each operating on the shared ctx / ctx.store.
        # Stage 1: full-length isoforms from intron-graph paths.
        self.fl_constructor = FLGraphConstructor(ctx)
        # Stage 2: models from per-read reference assignments.
        self.assignment_constructor = AssignmentBasedConstructor(ctx)
        # Stage 3: polyA/TSS end refinement + alternative-end NIC discovery.
        self.end_processor = TranscriptEndProcessor(ctx)
        # Stage 4a: structural filtering + splice-site correction (uses stage 3).
        self.model_filter = ModelFilter(ctx, self.end_processor)
        # Stage 4b: construction-phase read (re)assignment + ownership resolution.
        self.read_assigner = ConstructionReadAssigner(ctx)
        # Stage 4c: final honest per-read assignment + counting + read_info.
        self.model_read_counter = ModelReadCounter(ctx, transcript_counter, gene_counter)

    @property
    def transcript_model_storage(self):
        # Public accessor for the final model list (read by parallel_workers /
        # printers); the list itself lives on the shared ModelStore.
        return self.store.transcript_model_storage

    @property
    def model_read_assignments(self):
        # Honest per-read ReadAssignments against the final model set (read by
        # parallel_workers for the model read_info output); produced by stage 4c.
        return self.model_read_counter.model_read_assignments

    def get_transcript_id(self):
        return self.id_distributor.increment()

    def process(self, read_assignment_storage):
        self.ctx.build_graph(read_assignment_storage)

        self.fl_constructor.construct_fl_isoforms()
        self.assignment_constructor.construct_assignment_based_isoforms(read_assignment_storage)

        self.model_filter.pre_filter_transcripts()
        self.read_assigner.assign_reads_to_models(read_assignment_storage)
        self.model_filter.filter_transcripts()
        self.model_filter.correct_transcript_splice_sites()
        # reassign reads
        self.read_assigner.assign_reads_to_models(read_assignment_storage)
        # hand reads shared between an FL-path model and an assignment-based known back
        # to the FL-path model (novel-vs-known and known-vs-known), dropping empty knowns
        self.read_assigner.resolve_read_ownership_conflicts()

        transcript_joiner = TranscriptToGeneJoiner(self.store.transcript_model_storage, self.gene_info)
        self.store.transcript_model_storage = transcript_joiner.join_transcripts()
        self.model_read_counter.build_model_read_assignments(read_assignment_storage)

        if self.args.sqanti_output:
            self.compare_models_with_known()

    def compare_models_with_known(self):
        if not self.gene_info.all_isoforms_exons:
            for model in self.store.transcript_model_storage:
                # create intergenic
                assignment = ReadAssignment(model.transcript_id,
                                            ReadAssignmentType.intergenic,
                                            self.string_pools,
                                            match=IsoformMatch(MatchClassification.intergenic, string_pools=self.string_pools))
                if model.strand == "-":
                    polya_info = PolyAInfo(-1, model.exon_blocks[0][0], -1, -1)
                else:
                    polya_info = PolyAInfo(model.exon_blocks[-1][1], -1, -1, -1)

                assignment.polya_info = polya_info
                assignment.cage_found = False
                assignment.exons = model.exon_blocks
                assignment.strand = model.strand
                assignment.chr_id = model.chr_id
                assignment.set_additional_info("indel_count", "NA")
                assignment.set_additional_info("junctions_with_indels", "NA")
                assignment.introns_match = False
                assignment.gene_info = self.gene_info

                FSM_class = "C"
                assignment.set_additional_info("FSM_class", FSM_class)
                self.transcript2transcript.append(assignment)
            return

        gene_to_model_dict = defaultdict(list)
        for model in self.store.transcript_model_storage:
            gene_to_model_dict[model.gene_id].append(model.transcript_id)

        self.transcript2transcript = []
        for model in self.store.transcript_model_storage:
            if model.transcript_type == TranscriptModelType.known:
                continue
            if model.strand == "-":
                polya_info = PolyAInfo(-1, model.exon_blocks[0][0], -1, -1)
            else:
                polya_info = PolyAInfo(model.exon_blocks[-1][1], -1, -1, -1)

            combined_profile = self.profile_constructor.construct_profiles(model.exon_blocks, polya_info, [])
            assignment = self.assigner.assign_to_isoform(model.transcript_id, combined_profile)
            if assignment is None:
                continue

            assignment.polya_info = polya_info
            assignment.cage_found = False
            assignment.exons = model.exon_blocks
            assignment.strand = model.strand
            assignment.chr_id = model.chr_id
            assignment.set_additional_info("indel_count", "NA")
            assignment.set_additional_info("junctions_with_indels", "NA")
            assignment.introns_match = all(e == 1 for e in combined_profile.read_intron_profile.read_profile)
            assignment.gene_info = self.gene_info

            if assignment.assignment_type.is_unassigned() or \
                    not assignment.isoform_matches:
                # create intergenic
                assignment.assignment_type = ReadAssignmentType.intergenic
                FSM_class = "C"
                assignment.set_additional_info("FSM_class", FSM_class)
                self.transcript2transcript.append(assignment)
                continue

            if len(gene_to_model_dict[assignment.isoform_matches[0].assigned_gene]) == 1:
                FSM_class = "A"
            else:
                FSM_class = "C"
            assignment.set_additional_info("FSM_class", FSM_class)

            assigned_transcript_id = assignment.isoform_matches[0].assigned_transcript
            if not assigned_transcript_id or assigned_transcript_id not in self.gene_info.all_isoforms_introns:
                continue

            reference_introns = self.gene_info.all_isoforms_introns[assigned_transcript_id]
            isoform_introns = junctions_from_blocks(model.exon_blocks)
            event_string = ",".join([match_subtype_to_str_with_additional_info(x, model.strand,
                                                                               isoform_introns, reference_introns)
                                     for x in assignment.isoform_matches[0].match_subclassifications])
            model.add_additional_attribute("similar_reference_id", assigned_transcript_id)
            model.add_additional_attribute("alternatives", event_string)
            self.transcript2transcript.append(assignment)
