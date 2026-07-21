############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Stage 1 of model construction: full-length isoforms from intron-graph paths.
# Each FL path becomes a known reference isoform, an alternative-end NIC, or a
# novel (NIC/NNIC) model; reads are bound to the resulting model in the store.

import logging
from functools import cmp_to_key
from typing import Optional, Tuple

from ..common import TranscriptNaming, cmp, get_exons
from ..gene_info import TranscriptModel, TranscriptModelType
from isoquant_lib.model_construction.intron_graph import TerminalVertex
from ..assignment.isoform_assignment import is_matching_assignment, ReadAssignment
from ..terminal_prediction.polya_finder import PolyAInfo
from .context import StrandnessReportingLevel

logger = logging.getLogger('IsoQuant')


class GraphBasedConstructor:
    """Builds full-length transcript models from intron-graph FL paths."""

    def __init__(self, ctx):
        self.ctx = ctx
        self.store = ctx.store
        self.args = ctx.args
        self.gene_info = ctx.gene_info
        self.use_tss_model = ctx.use_tss_model
        self.create_alt_ts_nics = ctx.create_alt_ts_nics
        self.assigner = ctx.assigner
        self.profile_constructor = ctx.profile_constructor
        self.strand_detector = ctx.strand_detector
        self.use_technical_replicas = ctx.use_technical_replicas
        self.file_name_group_idx = ctx.file_name_group_idx

    def _reference_isoform_for_path(self, assignment: ReadAssignment, path: tuple, intron_path: tuple,
                                    transcript_range: Tuple[int, int]) -> Optional[str]:
        # Reference isoform whose intron chain matches this path. It is reported
        # as the annotated known UNLESS a *detected polyA* terminal vertex
        # (TerminalVertex.polya / TerminalVertex.polyt) disagrees with the annotated end by more
        # than apa_delta -> then the caller emits a novel-in-catalog isoform with
        # the refined ends. A bare read_end / read_start (no polyA evidence, e.g.
        # a degraded ONT terminus) never triggers reclassification, so known
        # transcripts are not lost to noise.
        if is_matching_assignment(assignment):
            matched_reference_id = assignment.isoform_matches[0].assigned_transcript
        elif intron_path in self.ctx.known_isoforms_in_graph:
            matched_reference_id = self.ctx.known_isoforms_in_graph[intron_path]
        else:
            return None

        if matched_reference_id is None or \
                matched_reference_id not in self.gene_info.all_isoforms_exons:
            return None

        ref_exons = self.gene_info.all_isoforms_exons[matched_reference_id]
        left_diff = abs(transcript_range[0] - ref_exons[0][0]) > self.args.apa_delta
        right_diff = abs(transcript_range[1] - ref_exons[-1][1]) > self.args.apa_delta

        # PolyA side: a detected polyA vertex (TerminalVertex.polyt on the genomic-left
        # 3' end of a '-' transcript, TerminalVertex.polya on the genomic-right 3' end of
        # a '+' transcript) that disagrees with the annotation -> alternative
        # polyA NIC. Bare read termini never trigger here (degraded-end safety).
        if path[0][0] == TerminalVertex.polyt and left_diff:
            return None
        if path[-1][0] == TerminalVertex.polya and right_diff:
            return None

        # TSS side (only with --fl_data): a *confirmed* TSS vertex (snapped to a TSS
        # prediction, tss_left on the genomic-left 5' end of a '+' transcript /
        # tss_right on the genomic-right 5' end of a '-' transcript) that disagrees
        # with the annotation -> alternative-TSS NIC. Symmetric to the polyA block:
        # a bare read_start / read_end never triggers, so ordinary 5' read-start
        # scatter/extension does not spawn a NIC that steals the known's reads.
        if self.use_tss_model:
            strand = self.gene_info.isoform_strands.get(matched_reference_id, '.')
            if strand == '+' and path[0][0] == TerminalVertex.tss_left and left_diff:
                return None
            if strand == '-' and path[-1][0] == TerminalVertex.tss_right and right_diff:
                return None

        return matched_reference_id

    def construct_fl_isoforms(self):
        # a minor trick to compare tuples of pairs, whose starting and terminating elements have different type
        logger.debug("Total FL paths %d" % len(self.ctx.path_storage.fl_paths))
        for path in sorted(self.ctx.path_storage.fl_paths,
                           key=cmp_to_key(lambda x,y: cmp(x,y) if len(x)==len(y) else cmp(len(y), len(x)))):
            # do not include terminal vertices
            # logger.debug(">>> Considering path " + str(path))
            intron_path = path[1:-1]
            if not intron_path: continue
            transcript_range = (path[0][1], path[-1][1])
            novel_exons = get_exons(transcript_range, list(intron_path))
            count = self.ctx.path_storage.paths[path]
            new_transcript_id = TranscriptNaming.transcript_prefix + str(self.ctx.get_transcript_id())
            # logger.debug("uuu %s: %s" % (new_transcript_id, str(novel_exons)))

            reference_isoform = None
            # check if new transcript matches a reference one
            if intron_path[0][0] == TerminalVertex.polyt:
                polya_info = PolyAInfo(-1, intron_path[0][1], -1, -1)
            elif intron_path[-1][0] == TerminalVertex.polya:
                polya_info = PolyAInfo(intron_path[-1][1], -1, -1, -1)
            else:
                polya_info = PolyAInfo(-1, -1, -1, -1)
            combined_profile = self.profile_constructor.construct_profiles(novel_exons, polya_info, [])
            assignment = self.assigner.assign_to_isoform(new_transcript_id, combined_profile)
            # check that no serious contradiction occurs
            # logger.debug("uuu Checking novel transcript %s: %s; assignment type %s" %
            #             (new_transcript_id, str(novel_exons), str(assignment.assignment_type)))

            if self.create_alt_ts_nics:
                # use the path's goal-1-refined terminal vertices: keep the
                # reference only when both ends agree with the annotation within
                # apa_delta, otherwise fall through to the novel branch below to
                # emit an alternative-end NIC built from the refined exons
                reference_isoform = self._reference_isoform_for_path(assignment, path, intron_path, transcript_range)
            elif is_matching_assignment(assignment):
                reference_isoform = assignment.isoform_matches[0].assigned_transcript
                # logger.debug("uuu Substituting with known isoform %s" % reference_isoform)
            elif intron_path in self.ctx.known_isoforms_in_graph:
                # path was not assigned to any known isoform but intron chain still matches
                continue

            new_model = None
            if reference_isoform:
                # adding FL reference isoform
                if reference_isoform in self.store.detected_known_isoforms:
                    pass
                elif count < self.args.min_known_count:
                    pass # logger.debug("uuu Isoform %s has low coverage %d" % (reference_isoform, count))
                else:
                    new_model = self.store.transcript_from_reference(reference_isoform)
                    self.store.detected_known_isoforms.add(reference_isoform)
                    #logger.debug("Adding known spliced isoform %s" % reference_isoform)
                    #logger.debug("Annotated positions: %d, %d, %s" % (new_model.exon_blocks[0][0], new_model.exon_blocks[-1][1], new_model.strand))
                    #logger.debug("Graph positions: %s, %s" % (str(path[0]), str(path[-1])))
            else:
                # adding FL novel isoform
                # component_coverage = self.intron_graph.get_max_component_coverage(intron_path)
                novel_isoform_cutoff = self.args.min_novel_count

                has_polyt = path[0][0] == TerminalVertex.polyt
                has_polya = path[-1][0] == TerminalVertex.polya
                polya_site = has_polya or has_polyt
                transcript_strand = self.strand_detector.get_strand(intron_path, has_polya, has_polyt)
                transcript_clean_strand = self.strand_detector.get_clean_strand(intron_path)

                #logger.debug("uuu Novel isoform %s has coverage: %d cutoff = %d, component cov = %d, max_coverage = %d"
                #             % (new_transcript_id, count, novel_isoform_cutoff, component_coverage, self.intron_graph.max_coverage))
                if count < novel_isoform_cutoff:
                    # logger.debug("uuu Novel isoform %s has low coverage: %d\t%d" % (new_transcript_id, count, novel_isoform_cutoff))
                    pass
                elif (len(novel_exons) == 2 and
                      ((self.args.require_monointronic_polya and not polya_site) or transcript_clean_strand == '.')):
                    # logger.debug("uuu Avoiding single intron %s isoform: %d\t%s" % (new_transcript_id, count, str(path)))
                    pass
                elif ((self.args.report_canonical_strategy == StrandnessReportingLevel.only_canonical
                       and transcript_clean_strand == '.')
                      or (self.args.report_canonical_strategy == StrandnessReportingLevel.only_stranded
                          and transcript_strand == '.')):
                    logger.debug("Avoiding unreliable transcript with %d exons (strand cannot be detected)" % len(novel_exons))
                    pass
                else:
                    if self.use_technical_replicas and self.file_name_group_idx >= 0:
                        # Check if reads come from same file (technical replicates)
                        read_assignments = self.ctx.path_storage.paths_to_reads[path]
                        if read_assignments:
                            file_groups = set([a.read_group[self.file_name_group_idx]
                                             for a in read_assignments
                                             if self.file_name_group_idx < len(a.read_group)])
                            if len(file_groups) <= 1:
                                #logger.debug("%s was suspended due to technical replicas check" % new_transcript_id)
                                continue

                    transcript_gene = self.ctx.select_reference_gene(intron_path, transcript_range, transcript_strand)
                    if transcript_gene is None:
                        transcript_gene = (TranscriptNaming.novel_gene_prefix + self.gene_info.chr_id +
                                           "_" + str(self.ctx.get_transcript_id()))
                    elif transcript_strand == '.':
                        transcript_strand = self.gene_info.gene_strands[transcript_gene]

                    if all(intron in self.ctx.known_introns for intron in intron_path):
                        transcript_type = TranscriptModelType.novel_in_catalog
                        id_suffix = TranscriptNaming.nic_transcript_suffix
                    else:
                        transcript_type = TranscriptModelType.novel_not_in_catalog
                        id_suffix = TranscriptNaming.nnic_transcript_suffix

                    new_model = TranscriptModel(self.gene_info.chr_id, transcript_strand,
                                                new_transcript_id + ".%s" % self.gene_info.chr_id + id_suffix,
                                                transcript_gene, novel_exons, transcript_type)
                    new_model.intron_path = intron_path
                    logger.debug("uuu Adding novel spliced isoform %s : %d\t%d" % (new_transcript_id, count, novel_isoform_cutoff))

            if new_model:
                self.store.transcript_model_storage.append(new_model)
                for read_assignment in self.ctx.path_storage.paths_to_reads[path]:
                    self.store.save_assigned_read(read_assignment, new_model.transcript_id)
