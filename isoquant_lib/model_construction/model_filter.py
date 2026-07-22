############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Stage 4a: structural filtering of constructed models. Prefilters low-coverage
# simple models, drops similar / unconfirmed-truncation isoforms, corrects novel
# splice sites, and interleaves the terminal end-refinement stage (delegated to
# TranscriptEndProcessor) at the same points as before.

import logging

from ..gene_info import GeneInfo, TranscriptModelType
from isoquant_lib.model_construction.intron_graph import TerminalVertex
from ..assignment.isoform_assignment import is_matching_assignment, ReadAssignmentType
from ..assignment.long_read_assigner import LongReadAssigner
from ..assignment.long_read_profiles import CombinedProfileConstructor
from ..terminal_prediction.polya_finder import PolyAInfo
from isoquant_lib.model_construction.transcript_splice_site_corrector import (
    count_deletions_for_splice_site_locations,
    correct_splice_site_errors,
    generate_updated_exon_list,
    SUPPORTED_STRANDS,
)

logger = logging.getLogger('IsoQuant')


class ModelFilter:
    """Structural filtering + splice-site correction of constructed models.

    Delegates polyA/TSS end refinement to the TranscriptEndProcessor at the same
    points inside filter_transcripts as the original monolith did.
    """

    def __init__(self, ctx, end_processor):
        self.ctx = ctx
        self.store = ctx.store
        self.args = ctx.args
        self.chr_record = ctx.chr_record
        self.string_pools = ctx.string_pools
        self.create_nics = ctx.create_alt_ts_nics
        self.novel_apa = ctx.novel_apa
        self.use_tss_model = ctx.use_tss_model
        self.end_processor = end_processor

    def pre_filter_transcripts(self):
        internal_count_values = []
        for model in self.store.transcript_model_storage:
            if len(model.exon_blocks) <= 2:
                internal_count_values.append(self.store.internal_counter[model.transcript_id])

        internal_count_values = sorted(internal_count_values, reverse=True)
        coverage_cutoff = self.args.min_novel_count
        # dirty hack to avoid slow assignment in chrM
        if len(internal_count_values) > 50:
            coverage_cutoff = internal_count_values[50]

        filtered_storage = []
        for model in self.store.transcript_model_storage:
            if len(model.exon_blocks) > 2:
                filtered_storage.append(model)
                continue

            if (model.transcript_type != TranscriptModelType.known and
                    self.store.internal_counter[model.transcript_id] < coverage_cutoff):
                self.store.delete_from_storage(model.transcript_id)
                continue

            if (model.transcript_type != TranscriptModelType.known and
                    self.store.mapping_quality(model.transcript_id) < self.args.simple_models_mapq_cutoff):
                self.store.delete_from_storage(model.transcript_id)
                continue
            filtered_storage.append(model)

        self.store.transcript_model_storage = filtered_storage

    def filter_transcripts(self):
        pre_filtered_storage = []
        to_substitute = self.detect_similar_isoforms(self.store.transcript_model_storage)

        for model in self.store.transcript_model_storage:
            if model.transcript_type == TranscriptModelType.known:
                pre_filtered_storage.append(model)
                continue
            # check coverage
            component_coverage = self.ctx.intron_graph.get_max_component_coverage(model.intron_path)
            if component_coverage == 0 or len(model.intron_path) == 0 or \
                    (len(model.intron_path) == 1 and self.ctx.intron_graph.is_monointron(model.intron_path[0])):
                component_coverage = self.ctx.intron_graph.get_overlapping_component_max_coverage((model.get_start(), model.get_end()))
                novel_isoform_cutoff = max(self.args.min_novel_count,
                                           self.args.min_mono_count_rel * component_coverage)
            else:
                novel_isoform_cutoff = max(self.args.min_novel_count,
                                           self.args.min_novel_count_rel * component_coverage)

            if model.transcript_id in to_substitute:
                # logger.debug("Novel model %s has a similar isoform %s" % (model.transcript_id, to_substitute[model.transcript_id]))
                self.store.delete_from_storage(model.transcript_id)
                continue

            if self.store.internal_counter[model.transcript_id] < novel_isoform_cutoff:
                # logger.debug("Novel model %s has coverage %d < %.2f, component cov = %d" % (model.transcript_id,
                #                                                        self.store.internal_counter[model.transcript_id],
                #                                                        novel_isoform_cutoff, component_coverage))
                self.store.delete_from_storage(model.transcript_id)
                continue

            if len(model.exon_blocks) <= 2:
                mapq = self.store.mapping_quality(model.transcript_id)
                # logger.debug("Novel model %s has quality %.2f" % (model.transcript_id, mapq))
                if mapq < self.args.simple_models_mapq_cutoff:
                    # logger.debug("Novel model %s has poor quality" % model.transcript_id)
                    self.store.delete_from_storage(model.transcript_id)
                    continue

            self.end_processor.correct_novel_transcript_ends(model, self.store.transcript_read_ids[model.transcript_id])
            self.end_processor._mark_terminal_confirmation(model)
            pre_filtered_storage.append(model)

        filtered_storage = []
        to_substitute = self.detect_similar_isoforms(pre_filtered_storage, remove_unconfirmed=True)

        for model in pre_filtered_storage:
            if model.transcript_type == TranscriptModelType.known:
                filtered_storage.append(model)
                continue

            if model.transcript_id in to_substitute:
                # logger.debug("Novel model %s has a similar isoform %s" % (model.transcript_id, to_substitute[model.transcript_id]))
                self.store.delete_from_storage(model.transcript_id)
                continue

            filtered_storage.append(model)

        self.store.transcript_model_storage = filtered_storage
        if self.create_nics or self.novel_apa:
            self.end_processor._add_known_alternative_end_models()
            self.end_processor._drop_duplicate_alt_end_models()

    def correct_transcript_splice_sites(self) -> None:
        """Shift novel transcript-model splice sites onto canonical motifs when a
        consistent cluster of read deletions next to a junction indicates the
        aligner mis-placed it. Known reference junctions come from the annotation
        and are left untouched."""
        if getattr(self.args, "no_splice_site_correction", False):
            return
        for model in self.store.transcript_model_storage:
            if model.transcript_type == TranscriptModelType.known:
                continue
            corrected_exons = self._correct_one_model_splice_sites(
                model, self.store.transcript_read_ids[model.transcript_id])
            if corrected_exons and corrected_exons != model.exon_blocks and \
                    self._valid_exon_chain(corrected_exons):
                logger.debug("Corrected splice sites for %s: %s -> %s" %
                             (model.transcript_id, model.exon_blocks, corrected_exons))
                model.exon_blocks = corrected_exons

    def _correct_one_model_splice_sites(self, model, reads):
        """Return a corrected exon-block list for a single model, or None."""
        exons = model.exon_blocks
        if len(exons) < 2 or not reads or model.strand not in SUPPORTED_STRANDS:
            return None

        splice_site_cases = {}
        for read in reads:
            if not read.exons:
                continue
            count_deletions_for_splice_site_locations(
                read.exons[0][0], read.exons[-1][1], read.boundary_deletions or [],
                exons, splice_site_cases)

        if not splice_site_cases:
            return None
        locations_with_errors = correct_splice_site_errors(
            splice_site_cases, self.chr_record, model.strand)
        if not locations_with_errors:
            return None
        return generate_updated_exon_list(splice_site_cases, locations_with_errors, exons)

    @staticmethod
    def _valid_exon_chain(exons) -> bool:
        """Guard against a correction that inverts an exon or collides with a
        neighbour (introns must keep at least one base)."""
        for i, (start, end) in enumerate(exons):
            if start > end:
                return False
            if i > 0 and start <= exons[i - 1][1]:
                return False
        return True

    def detect_similar_isoforms(self, model_storage, remove_unconfirmed=False):
        to_substitute = {}
        for model in model_storage:
            if len(model.exon_blocks) <= 2 or model.transcript_id in to_substitute:
                continue
            transcript_model_gene_info = GeneInfo.from_models([model], self.args.delta)
            assigner = LongReadAssigner(transcript_model_gene_info, self.args, self.string_pools)
            profile_constructor = CombinedProfileConstructor(transcript_model_gene_info, self.args)

            for m in model_storage:
                if m.transcript_type == TranscriptModelType.known or m.transcript_id == model.transcript_id or \
                        m.transcript_id in to_substitute or len(m.exon_blocks) == 1 or not m.intron_path or \
                        len(m.exon_blocks) > len(model.exon_blocks):
                    continue

                if m.intron_path[0][0] == TerminalVertex.polyt:
                    polya_info = PolyAInfo(-1, m.intron_path[0][1], -1, -1)
                elif m.intron_path[-1][0] == TerminalVertex.polya:
                    polya_info = PolyAInfo(m.intron_path[-1][1], -1, -1, -1)
                else:
                    polya_info = PolyAInfo(-1, -1, -1, -1)
                combined_profile = profile_constructor.construct_profiles(m.exon_blocks, polya_info, [])
                assignment = assigner.assign_to_isoform(m.transcript_id, combined_profile)

                if is_matching_assignment(assignment):
                    to_substitute[m.transcript_id] = model.transcript_id
                elif (remove_unconfirmed and
                      assignment.assignment_type in (ReadAssignmentType.unique_minor_difference,
                                                     ReadAssignmentType.inconsistent_non_intronic) and
                      not (m.polya_confirmed or (self.use_tss_model and m.tss_confirmed))):
                    # m is a terminal/ISM sub-structure of the longer model: its intron chain is a
                    # subset and only the terminal exons differ (unique_minor_difference with an
                    # ISM/elongation event, or inconsistent_non_intronic), and its polyA/TSS
                    # terminus is not confirmed by a read-terminus peak -> truncation artifact, not
                    # a real alternative-end transcript. A genuinely different intron chain gives
                    # plain `inconsistent` and is left alone. Remove it (reads reassigned to model).
                    to_substitute[m.transcript_id] = model.transcript_id

        return to_substitute
