############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Stage 2 of model construction: models from per-read reference assignments.
# Partitions unique reads and emits known spliced (non-FL) isoforms, known
# mono-exon isoforms, and novel mono-exon isoforms, binding reads in the store.

import logging
from collections import defaultdict

from ..common import TranscriptNaming, interval_len, intersection_len
from ..gene_info import TranscriptModel, TranscriptModelType
from isoquant_lib.model_construction.intron_graph import TerminalVertex
from ..isoform_assignment import MatchEventSubtype

logger = logging.getLogger('IsoQuant')


class AssignmentBasedConstructor:
    """Builds transcript models from per-read reference assignments."""

    def __init__(self, ctx):
        self.ctx = ctx
        self.store = ctx.store
        self.args = ctx.args
        self.gene_info = ctx.gene_info

    def construct_assignment_based_isoforms(self, read_assignment_storage):
        spliced_isoform_reads = defaultdict(list)
        isoform_left_support = defaultdict(int)
        isoform_right_support = defaultdict(int)
        polya_sites = defaultdict(int)
        mono_exon_isoform_reads = defaultdict(list)
        mono_exon_isoform_coverage = {}
        novel_mono_exon_reads = []

        for read_assignment in read_assignment_storage:
            if len(read_assignment.corrected_exons) <= 2 and \
                    (read_assignment.multimapper or read_assignment.mapping_quality < self.args.simple_alignments_mapq_cutoff):
                continue

            if not read_assignment:
                continue
            if (len(read_assignment.corrected_exons) == 1 and
                    read_assignment.polyA_found and not read_assignment.multimapper and
                    (read_assignment.assignment_type.is_unassigned() or
                     read_assignment.assignment_type.is_inconsistent())):
                novel_mono_exon_reads.append(read_assignment)

            if not read_assignment.assignment_type.is_unique():
                continue

            refrenence_isoform_id = read_assignment.isoform_matches[0].assigned_transcript
            if refrenence_isoform_id in self.store.detected_known_isoforms:
                continue

            events = read_assignment.isoform_matches[0].match_subclassifications
            if any(e.event_type == MatchEventSubtype.mono_exon_match for e in events):
                mono_exon_isoform_reads[refrenence_isoform_id].append(read_assignment)
                assert len(self.gene_info.all_isoforms_introns[refrenence_isoform_id]) == 0
                transcript_start = self.gene_info.all_isoforms_exons[refrenence_isoform_id][0][0]
                transcript_end = self.gene_info.all_isoforms_exons[refrenence_isoform_id][-1][1]
                t_len = transcript_end - transcript_start + 1

                if refrenence_isoform_id not in mono_exon_isoform_coverage:
                    mono_exon_isoform_coverage[refrenence_isoform_id] = [0 for _ in range(t_len)]
                start = max(0, read_assignment.corrected_exons[0][0] - transcript_start)
                end = min(t_len, read_assignment.corrected_exons[-1][1] - transcript_start + 1)
                for i in range(start, end):
                    mono_exon_isoform_coverage[refrenence_isoform_id][i] = 1

                if self.gene_info.isoform_strands[refrenence_isoform_id] == '+':
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_right for x in events):
                        polya_sites[refrenence_isoform_id] += 1
                else:
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_left for x in events):
                        polya_sites[refrenence_isoform_id] += 1
            elif len(self.gene_info.all_isoforms_exons[refrenence_isoform_id]) > 1:
                if self.store.read_assignment_counts[read_assignment.read_id] > 0:
                    pass
                    # logger.debug("Spliced read %s was previously used for construction, assigned id %s" %
                    #             (read_assignment.read_id, refrenence_isoform_id))
                spliced_isoform_reads[refrenence_isoform_id].append(read_assignment)

                if self.args.requires_polya_for_construction and self.gene_info.isoform_strands[refrenence_isoform_id] == '-':
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_left for x in events):
                        isoform_left_support[refrenence_isoform_id] += 1
                elif abs(self.gene_info.all_isoforms_exons[refrenence_isoform_id][0][0] - read_assignment.corrected_exons[0][0]) <= self.args.apa_delta:
                    isoform_left_support[refrenence_isoform_id] += 1

                if self.args.requires_polya_for_construction and self.gene_info.isoform_strands[refrenence_isoform_id] == '+':
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_right for x in events):
                        isoform_right_support[refrenence_isoform_id] += 1
                elif abs(self.gene_info.all_isoforms_exons[refrenence_isoform_id][-1][1] - read_assignment.corrected_exons[-1][1]) <= self.args.apa_delta:
                    isoform_right_support[refrenence_isoform_id] += 1

        self.construct_monoexon_isoforms(mono_exon_isoform_reads, mono_exon_isoform_coverage, polya_sites)
        if not self.args.fl_only:
            logger.debug("Constructing nonFL isoforms")
            self.construct_nonfl_isoforms(spliced_isoform_reads, isoform_left_support, isoform_right_support)
        if self.args.report_novel_unspliced:
            self.construct_monoexon_novel(novel_mono_exon_reads)

    def collect_terminal_exons_from_graph(self):
        polya_exons = []
        polyt_exons = []
        for intron in self.ctx.intron_graph.outgoing_edges.keys():
            for v in self.ctx.intron_graph.outgoing_edges[intron]:
                if v[0] == TerminalVertex.polya:
                    polya_exons.append((intron[1], v[1]))
        for intron in self.ctx.intron_graph.incoming_edges.keys():
            for v in self.ctx.intron_graph.incoming_edges[intron]:
                if v[0] == TerminalVertex.polyt:
                    polyt_exons.append((v[1], intron[0]))
        # logger.debug("PolyA terminal exons: " + str(polya_exons))
        # logger.debug("PolyT terminal exons: " + str(polyt_exons))
        return polya_exons, polyt_exons

    def is_internal_monoexonic_read(self, alignment, terminal_exons, forward=True):
        read_coordinates = alignment.corrected_exons[0]
        if forward:
            for e in terminal_exons:
                if abs(e[1] - alignment.corrected_exons[-1][1]) <= self.args.apa_delta and \
                        read_coordinates[0] >= e[0] - self.args.delta:
                    return True
        else:
            for e in terminal_exons:
                if abs(e[0] - alignment.corrected_exons[0][0]) <= self.args.apa_delta and \
                        read_coordinates[1] <= e[1] + self.args.delta:
                    return True
        return False

    def construct_monoexon_novel(self, novel_mono_exon_reads):
        logger.debug("Constructing novel monoexon")
        polya_exons, polyt_exons = self.collect_terminal_exons_from_graph()
        polya_reads = defaultdict(list)
        polyt_reads = defaultdict(list)
        for a in novel_mono_exon_reads:
            if a.polya_info.external_polya_pos != -1:
                if not self.is_internal_monoexonic_read(a, polya_exons, forward=True):
                    polya_reads[a.polya_info.external_polya_pos].append(a)
            if a.polya_info.external_polyt_pos != -1:
                if not self.is_internal_monoexonic_read(a, polyt_exons, forward=False):
                    polyt_reads[a.polya_info.external_polyt_pos].append(a)

        novel_monoexon = set()
        clustered_polya_reads = self.cluster_monoexons(polya_reads)
        novel_monoexon.update(self.generate_monoexon_from_clustered(clustered_polya_reads, True))
        clustered_polyt_reads = self.cluster_monoexons(polyt_reads)
        novel_monoexon.update(self.generate_monoexon_from_clustered(clustered_polyt_reads, False))

    def generate_monoexon_from_clustered(self, clustered_reads, forward=True):
        cutoff = self.args.min_novel_count
        result = set()
        for three_prime_pos in clustered_reads.keys():
            count = len(clustered_reads[three_prime_pos])
            if count < cutoff:
                continue
            # TODO improve
            if forward:
                five_prime_pos = min([a.corrected_exons[0][0] for a in clustered_reads[three_prime_pos]])
            else:
                five_prime_pos = max([a.corrected_exons[-1][1] for a in clustered_reads[three_prime_pos]])

            strand = '+' if forward else '-'
            coordinates = (five_prime_pos, three_prime_pos) if forward else (three_prime_pos, five_prime_pos)
            new_transcript_id = TranscriptNaming.transcript_prefix + str(self.ctx.get_transcript_id())
            transcript_gene = (TranscriptNaming.novel_gene_prefix + self.gene_info.chr_id +
                               "_" + str(self.ctx.get_transcript_id()))
            transcript_type = TranscriptModelType.novel_not_in_catalog
            id_suffix = TranscriptNaming.nnic_transcript_suffix

            is_valid = True
            half_len = interval_len(coordinates) / 2
            for existing_model in self.store.transcript_model_storage:
                if any(intersection_len(exon, coordinates) > half_len for exon in existing_model.exon_blocks):
                    is_valid = False
                    break
            if not is_valid:
                continue

            new_model = TranscriptModel(self.gene_info.chr_id, strand,
                                        new_transcript_id + ".%s" % self.gene_info.chr_id + id_suffix,
                                        transcript_gene, [coordinates], transcript_type)
            # logger.debug("uuu Adding novel MONOEXON isoform %s : %s, %d\t%d" % (new_transcript_id, str(coordinates), count, cutoff))
            result.add(coordinates)

            self.store.transcript_model_storage.append(new_model)
            for read_assignment in clustered_reads[three_prime_pos]:
                self.store.save_assigned_read(read_assignment, new_model.transcript_id)
        return result

    def cluster_monoexons(self, grouped_reads):
        clustered_counts = defaultdict(list)
        while grouped_reads:
            best_pair = max(grouped_reads.items(), key=lambda x:len(x[1]))
            top_position = best_pair[0]
            for pos in range(top_position - self.args.apa_delta, top_position + self.args.apa_delta + 1):
                if pos in grouped_reads:
                    clustered_counts[top_position] += grouped_reads[pos]
                    del grouped_reads[pos]

        return clustered_counts

    def construct_monoexon_isoforms(self, mono_exon_isoform_reads, mono_exon_isoform_coverage, polya_sites):
        for isoform_id in mono_exon_isoform_reads.keys():
            count = len(mono_exon_isoform_reads[isoform_id])
            coverage = float(mono_exon_isoform_coverage[isoform_id].count(1)) / \
                       float(len(mono_exon_isoform_coverage[isoform_id]))
            polya_support = polya_sites[isoform_id]

            # logger.debug(">> Monoexon transcript %s: %d\t%d\t%.4f\t%d" % (isoform_id, self.ctx.intron_graph.max_coverage, count, coverage, polya_support))
            if (count < self.args.min_known_count or coverage < self.args.min_mono_exon_coverage or
                    (self.args.require_monoexonic_polya and polya_support == 0)):
                pass # logger.debug(">> Will NOT be added, abs cutoff=%d" % (self.params.min_known_count))
            elif isoform_id not in self.store.detected_known_isoforms:
                new_model = self.store.transcript_from_reference(isoform_id)
                self.store.transcript_model_storage.append(new_model)
                self.store.detected_known_isoforms.add(isoform_id)
                for read_assignment in mono_exon_isoform_reads[isoform_id]:
                    self.store.save_assigned_read(read_assignment, new_model.transcript_id)
                # logger.debug(">> Adding known monoexon isoform %s, %s, count = %d: %s" %
                #             (self.store.transcript_model_storage[-1].transcript_id, isoform_id,
                #              count, str(self.gene_info.all_isoforms_exons[isoform_id])))

    def construct_nonfl_isoforms(self, spliced_isoform_reads, spliced_isoform_left_support, spliced_isoform_right_support):
        logger.debug("Constructing nonFL isoforms")
        for isoform_id in spliced_isoform_reads.keys():
            if isoform_id in self.store.detected_known_isoforms:
                continue
            count = len(spliced_isoform_reads[isoform_id])
            if isoform_id not in self.ctx.known_isoforms_in_graph_ids:
                #logger.debug("<< Isoform %s has %d assignments but is not in the graph" % (isoform_id, count))
                continue

            intron_path = self.ctx.known_isoforms_in_graph_ids[isoform_id]
            #logger.debug("Known non-FL spliced isoform %s" % isoform_id)
            if count < self.args.min_known_count or \
                    spliced_isoform_left_support[isoform_id] < 1 or \
                    spliced_isoform_right_support[isoform_id] < 1:
                pass
                # logger.debug("Will not be added")
            else:
                logger.debug("<< Adding known non-FL spliced isoform %s" % isoform_id)
                new_model = self.store.transcript_from_reference(isoform_id)
                self.store.transcript_model_storage.append(new_model)
                self.store.detected_known_isoforms.add(isoform_id)
                # Remember which reads built this known so that, after filtering,
                # resolve_read_ownership_conflicts can hand back any read still
                # shared with a surviving novel FL model.
                self.store.nonfl_known_reads[new_model.transcript_id] = list(spliced_isoform_reads[isoform_id])
                for read_assignment in spliced_isoform_reads[isoform_id]:
                    self.store.save_assigned_read(read_assignment, new_model.transcript_id)
