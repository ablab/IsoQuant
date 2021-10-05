############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
from collections import namedtuple
from functools import reduce
from functools import cmp_to_key
import copy

from src.common import *
from src.assignment_io import *
from src.isoform_assignment import *
from src.intron_graph import *
from src.gene_info import *
from src.alignment_info import *

logger = logging.getLogger('IsoQuant')


class GraphBasedModelConstructor:
    transcript_id_counter = AtomicCounter()
    transcript_prefix = "transcript_"
    known_transcript_suffix = ".known"
    nic_transcript_suffix = ".nic"
    nnic_transcript_suffix = ".nnic"

    def __init__(self, gene_info, chr_record, params, transcript_counter):
        self.gene_info = gene_info
        self.chr_record = chr_record
        self.params = params

        self.strand_detector = StrandDetector(self.chr_record)
        self.intron_genes = defaultdict(set)
        self.set_gene_properties()

        self.intron_graph = None
        self.path_processor = None
        self.path_storage = None
        self.detected_known_isoforms = set()

        self.known_isoforms_in_graph = {}
        self.known_introns = set()
        self.known_isoforms_in_graph_ids = {}
        self.assigner = LongReadAssigner(self.gene_info, self.params)
        self.profile_constructor = CombinedProfileConstructor(self.gene_info, self.params)

        self.transcript_model_storage = []
        self.transcript_read_ids = defaultdict(list)
        self.transcript_counter = transcript_counter
        self.reads_used_in_construction = set()
        self.unused_reads = []

    def get_transcript_id(self):
        return GraphBasedModelConstructor.transcript_id_counter.increment()

    def set_gene_properties(self):
        intron_strands_dicts = defaultdict(lambda: defaultdict(int))
        self.intron_genes = defaultdict(set)
        for t_id, introns in self.gene_info.all_isoforms_introns.items():
            strand = self.gene_info.isoform_strands[t_id]
            gene_id = self.gene_info.gene_id_map[t_id]
            for intron in introns:
                intron_strands_dicts[intron][strand] += 1
                self.intron_genes[intron].add(gene_id)

        for intron in intron_strands_dicts.keys():
            if len(intron_strands_dicts[intron].keys()) == 1:
                # intron has a single strand
                self.strand_detector.set_strand(intron, list(intron_strands_dicts[intron].keys())[0])
            else:
                self.strand_detector.set_strand(intron)

    def select_reference_gene(self, transcript_introns, transcript_range, transcript_strand):
        gene_counts = defaultdict(int)
        for intron in transcript_introns:
            if intron not in self.intron_genes:
                continue
            for g_id in self.intron_genes[intron]:
                gene_counts[g_id] += 1

        ordered_genes = sorted(gene_counts.items(), key=lambda x: (x[1], x[0]), reverse=True)
        for g in ordered_genes:
            gene_id = g[0]
            if transcript_strand == '.' or self.gene_info.gene_strands[gene_id] == transcript_strand:
                return gene_id

        overlap_dict = {}
        gene_regions = self.gene_info.get_gene_regions()
        for gene_id in gene_regions.keys():
            gene_coverage = read_coverage_fraction([transcript_range], [gene_regions[gene_id]])
            if gene_coverage > 0.0 and \
                    (transcript_strand == '.' or self.gene_info.gene_strands[gene_id] == transcript_strand):
                overlap_dict[gene_id] = gene_coverage

        if overlap_dict:
            return get_top_count(overlap_dict)
        return None

    def process(self, read_assignment_storage):
        self.intron_graph = IntronGraph(self.params, self.gene_info, read_assignment_storage)
        self.path_processor = IntronPathProcessor(self.params, self.intron_graph)
        self.path_storage = IntronPathStorage(self.params, self.path_processor)
        self.path_storage.fill(read_assignment_storage)
        self.known_isoforms_in_graph = self.get_known_spliced_isoforms(self.gene_info)
        self.known_introns = set(self.gene_info.intron_profiles.features)

        for intron_path, isoform_id in self.known_isoforms_in_graph.items():
            self.known_isoforms_in_graph_ids[isoform_id] = intron_path

        self.construct_fl_isoforms()
        self.construnct_assignment_based_isoforms(read_assignment_storage)
        self.assign_reads_to_models(read_assignment_storage)
        for model in self.transcript_model_storage:
            if model.transcript_id.endswith(self.known_transcript_suffix):
                continue
            self.correct_novel_transcrip_ends(model, self.transcript_read_ids[model.transcript_id])

    def get_known_spliced_isoforms(self, gene_info, s="known"):
        known_isoforms = {}
        for isoform_id in gene_info.all_isoforms_introns:
            isoform_introns = gene_info.all_isoforms_introns[isoform_id]
            intron_path = self.path_processor.thread_introns(isoform_introns)
            if not intron_path:
                # logger.debug("== No path found for %s isoform %s: %s" % (s, isoform_id, gene_info.all_isoforms_introns[isoform_id]))
                continue
            # logger.debug("== Path found for %s isoform %s: %s" % (s, isoform_id, gene_info.all_isoforms_introns[isoform_id]))
            known_isoforms[tuple(intron_path)] = isoform_id
        return known_isoforms

    def is_reference_isoform(self, isoform_assignment):
        if isoform_assignment.assignment_type == ReadAssignmentType.unique:
            return True
        elif isoform_assignment.assignment_type == ReadAssignmentType.unique_minor_difference:
            allowed_set = {MatchEventSubtype.none, MatchEventSubtype.exon_elongation_left, MatchEventSubtype.exon_elongation_right}
            return all(m.event_type in allowed_set for m in isoform_assignment.isoform_matches[0].match_subclassifications)
        return False

    def save_assigned_read(self, read_assignment, transcript_id):
        read_id = read_assignment.read_id
        self.transcript_read_ids[transcript_id].append(read_assignment)
        self.transcript_counter.add_read_info_raw(read_id, [transcript_id], read_assignment.read_group)

    def construct_fl_isoforms(self):
        novel_isoform_cutoff = max(self.params.min_novel_count, self.params.min_novel_count_rel * self.intron_graph.max_coverage)
        self.detected_known_isoforms = set()

        # a minor trick to compare tuples of pairs, whose starting and terminating elements have different type
        logger.debug("Total FL paths %d" % len(self.path_storage.fl_paths))
        for path in sorted(self.path_storage.fl_paths,
                           key=cmp_to_key(lambda x,y: cmp(x,y) if len(x)==len(y) else cmp(len(y), len(x)))):
            # do not include terminal vertices
            intron_path = path[1:-1]
            transcript_range = (path[0][1], path[-1][1])
            novel_exons = get_exons(transcript_range, list(intron_path))
            count = self.path_storage.paths[path]
            new_transcript_id = self.transcript_prefix + str(self.get_transcript_id())

            reference_isoform = None
            if intron_path not in self.known_isoforms_in_graph:
                # check if new transcript matches a reference one
                combined_profile = self.profile_constructor.construct_profiles(novel_exons, PolyAInfo(-1, -1, -1, -1), [])
                assignment = self.assigner.assign_to_isoform(new_transcript_id, combined_profile)
                # check that no serious contradiction occurs

                if self.is_reference_isoform(assignment):
                    reference_isoform = assignment.isoform_matches[0].assigned_transcript
                    logger.debug("Substituting with known isoform %s" % reference_isoform)
            else:
                # path matches reference exactly
                reference_isoform = self.known_isoforms_in_graph[intron_path]
                logger.debug("Matches with known isoform %s" % reference_isoform)

            new_model = None
            if reference_isoform:
                # adding FL reference isoform
                if reference_isoform not in self.detected_known_isoforms and count >= self.params.min_known_count:
                    new_model = self.transcript_from_reference(reference_isoform, new_transcript_id)
                    self.detected_known_isoforms.add(reference_isoform)
                    logger.debug("Adding known spliced isoform %s" % reference_isoform)
                    logger.debug("Annotated positions: %d, %d, %s" % (new_model.exon_blocks[0][0], new_model.exon_blocks[-1][1], new_model.strand))
                    logger.debug("Graph positions: %s, %s" % (str(path[0]), str(path[-1])))
            else:
                # adding FL novel isoform
                polya_site = (path[0][0] == VERTEX_polyt or path[-1][0] == VERTEX_polya)
                transcript_strand = self.strand_detector.get_strand(intron_path)
                if count < novel_isoform_cutoff:
                    logger.debug("Novel isoform %s has low coverage: %d\t%d" %
                                 (new_transcript_id, count, novel_isoform_cutoff))
                elif len(novel_exons) == 2 and (not polya_site or transcript_strand == '.'):
                    logger.debug("Avoiding single intron %s isoform: %d\t%s" %
                                 (new_transcript_id, count, str(path)))
                else:
                    if self.params.use_technical_replicas and \
                            len(set([a.read_group for a in self.path_storage.paths_to_reads[path]])) <= 1:
                        logger.debug("%s was suspended due to technical replicas check" % new_transcript_id)
                        continue

                    transcript_gene = self.select_reference_gene(intron_path, transcript_range, transcript_strand)
                    if transcript_gene is None:
                        transcript_gene = "novel_gene_" + str(self.get_transcript_id())
                    elif transcript_strand == '.':
                        transcript_strand = self.gene_info.gene_strands[transcript_gene]

                    if all(intron in self.known_introns for intron in intron_path):
                        transcript_type = TranscriptModelType.novel_in_catalog
                        id_suffix = self.nic_transcript_suffix
                    else:
                        transcript_type = TranscriptModelType.novel_not_in_catalog
                        id_suffix = self.nnic_transcript_suffix

                    new_model = TranscriptModel(self.gene_info.chr_id, transcript_strand,
                                                new_transcript_id + id_suffix, "novel", transcript_gene,
                                                novel_exons, transcript_type)
                    logger.debug("Adding novel spliced isoform %s : %d\t%d" %
                                 (new_transcript_id, count, novel_isoform_cutoff))

            if new_model:
                self.transcript_model_storage.append(new_model)
                for read_assignment in self.path_storage.paths_to_reads[path]:
                    self.save_assigned_read(read_assignment, new_model.transcript_id)
                    self.reads_used_in_construction.add(read_assignment.read_id)

    def construnct_assignment_based_isoforms(self, read_assignment_storage):
        spliced_isoform_reads = defaultdict(list)
        isoform_left_support = defaultdict(int)
        isoform_right_support = defaultdict(int)
        polya_sites = defaultdict(int)
        mono_exon_isoform_reads = defaultdict(list)
        mono_exon_isoform_coverage = {}

        for read_assignment in read_assignment_storage:
            if not read_assignment:
                continue
            if read_assignment.assignment_type not in {ReadAssignmentType.unique,
                                                       ReadAssignmentType.unique_minor_difference}:
                continue
            refrenence_isoform_id = read_assignment.isoform_matches[0].assigned_transcript
            if refrenence_isoform_id in self.detected_known_isoforms:
                continue

            events = read_assignment.isoform_matches[0].match_subclassifications
            if any(e.event_type == MatchEventSubtype.mono_exon_match for e in events):
                mono_exon_isoform_reads[refrenence_isoform_id].append(read_assignment)
                assert len(self.gene_info.all_isoforms_exons[refrenence_isoform_id]) == 1
                transcript_exon = self.gene_info.all_isoforms_exons[refrenence_isoform_id][0]
                t_len = transcript_exon[1] - transcript_exon[0] + 1

                if refrenence_isoform_id not in mono_exon_isoform_coverage:
                    mono_exon_isoform_coverage[refrenence_isoform_id] = [0 for _ in range(t_len)]
                start = max(0, read_assignment.corrected_exons[0][0] - transcript_exon[0])
                end = min(t_len, read_assignment.corrected_exons[-1][1] - transcript_exon[0] + 1)
                for i in range(start, end):
                    mono_exon_isoform_coverage[refrenence_isoform_id][i] = 1

                if self.gene_info.isoform_strands[refrenence_isoform_id] == '+':
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_right for x in events):
                        polya_sites[refrenence_isoform_id] += 1
                else:
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_left for x in events):
                        polya_sites[refrenence_isoform_id] += 1
            elif len(self.gene_info.all_isoforms_exons[refrenence_isoform_id]) > 1:
                spliced_isoform_reads[refrenence_isoform_id].append(read_assignment)

                if self.params.needs_polya_for_construction and self.gene_info.isoform_strands[refrenence_isoform_id] == '-':
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_left for x in events):
                        isoform_left_support[refrenence_isoform_id] += 1
                elif abs(self.gene_info.all_isoforms_exons[refrenence_isoform_id][0][0] - read_assignment.corrected_exons[0][0]) <= self.params.apa_delta:
                    isoform_left_support[refrenence_isoform_id] += 1

                if self.params.needs_polya_for_construction and self.gene_info.isoform_strands[refrenence_isoform_id] == '+':
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_right for x in events):
                        isoform_right_support[refrenence_isoform_id] += 1
                elif abs(self.gene_info.all_isoforms_exons[refrenence_isoform_id][-1][1] - read_assignment.corrected_exons[-1][1]) <= self.params.apa_delta:
                    isoform_right_support[refrenence_isoform_id] += 1

        self.construct_monoexon_isoforms(mono_exon_isoform_reads, mono_exon_isoform_coverage, polya_sites)
        if not self.params.fl_only:
            logger.debug("Constructing nonFL isoforms")
            self.construct_nonfl_isoforms(spliced_isoform_reads, isoform_left_support, isoform_right_support)

    def construct_monoexon_isoforms(self, mono_exon_isoform_reads, mono_exon_isoform_coverage, polya_sites):
        novel_isoform_cutoff = max(self.params.min_novel_count, self.params.min_novel_count_rel * self.intron_graph.max_coverage)

        for isoform_id in mono_exon_isoform_reads.keys():
            count = len(mono_exon_isoform_reads[isoform_id])
            coverage = float(mono_exon_isoform_coverage[isoform_id].count(1)) / \
                       float(len(mono_exon_isoform_coverage[isoform_id]))
            polya_support = polya_sites[isoform_id]

            # logger.debug(">> Monoexon transcript %s: %d\t%d\t%.4f\t%d" % (isoform_id, self.intron_graph.max_coverage, count, coverage, polya_support))
            if count < self.params.min_known_count or coverage < self.params.min_mono_exon_coverage or polya_support == 0:
                logger.debug("Will not be added %s" % isoform_id)
            else:
                new_model = self.transcript_from_reference(isoform_id)
                self.transcript_model_storage.append(new_model)
                self.detected_known_isoforms.add(isoform_id)
                for read_assignment in mono_exon_isoform_reads[isoform_id]:
                    self.save_assigned_read(read_assignment, new_model.transcript_id)
                    self.reads_used_in_construction.add(read_assignment.read_id)
                logger.debug(">> Adding known monoexon isoform %s, %s" %
                             (self.transcript_model_storage[-1].transcript_id, isoform_id))

    def construct_nonfl_isoforms(self, spliced_isoform_reads, spliced_isoform_left_support, spliced_isoform_right_support):
        logger.debug("Constructing nonFL isoforms")
        for isoform_id in spliced_isoform_reads.keys():
            if isoform_id in self.detected_known_isoforms:
                continue
            count = len(spliced_isoform_reads[isoform_id])
            if isoform_id not in self.known_isoforms_in_graph_ids:
                logger.debug("<< Isoform %s is not in the graph" % isoform_id)
                continue

            intron_path = self.known_isoforms_in_graph_ids[isoform_id]
            logger.debug("Known non-FL spliced isoform %s" % isoform_id)
            if count < self.params.min_known_count or \
                    spliced_isoform_left_support[isoform_id] < 1 or \
                    spliced_isoform_right_support[isoform_id] < 1:
                logger.debug("Will not be added")
            else:
                logger.debug("Adding known non-FL spliced isoform %s" % isoform_id)
                new_model = self.transcript_from_reference(isoform_id)
                self.transcript_model_storage.append(new_model)
                self.detected_known_isoforms.add(isoform_id)
                for read_assignment in spliced_isoform_reads[isoform_id]:
                    self.save_assigned_read(read_assignment, new_model.transcript_id)
                    self.reads_used_in_construction.add(read_assignment.read_id)

        for isoform_id in self.detected_known_isoforms:
            if isoform_id not in self.known_isoforms_in_graph_ids:
                continue
            path = self.known_isoforms_in_graph_ids[isoform_id]

    # create transcript model object from reference isoforms
    def transcript_from_reference(self, isoform_id, transcript_id=None):
        if not transcript_id:
            transcript_id = self.transcript_prefix + str(self.get_transcript_id())
        new_transcript_id = transcript_id + self.known_transcript_suffix
        return TranscriptModel(self.gene_info.chr_id, self.gene_info.isoform_strands[isoform_id],
                               new_transcript_id, isoform_id, self.gene_info.gene_id_map[isoform_id],
                               self.gene_info.all_isoforms_exons[isoform_id], TranscriptModelType.known)

    # assign reads back to constructed isoforms
    def assign_reads_to_models(self, read_assignments):
        if not self.transcript_model_storage:
            logger.debug("No transcripts were assigned")
            self.unused_reads = [a.read_id for a in read_assignments]
            return
        
        logger.debug("Creating artificial GeneInfo from %d transcript models" % len(self.transcript_model_storage))
        transcript_model_gene_info = GeneInfo.from_models(self.transcript_model_storage, self.params.delta)
        assigner = LongReadAssigner(transcript_model_gene_info, self.params, quick_mode=True)
        profile_constructor = CombinedProfileConstructor(transcript_model_gene_info, self.params)

        for assignment in read_assignments:
            read_id = assignment.read_id
            if read_id in self.reads_used_in_construction:
                continue

            read_exons = assignment.corrected_exons
            model_combined_profile = profile_constructor.construct_profiles(read_exons, assignment.polya_info, [])
            model_assignment = assigner.assign_to_isoform(assignment.read_id, model_combined_profile)
            model_assignment.read_group = assignment.read_group
            # check that no serious contradiction occurs
            if model_assignment.assignment_type in [ReadAssignmentType.unique,
                                                    ReadAssignmentType.unique_minor_difference,
                                                    ReadAssignmentType.ambiguous]:
                self.transcript_counter.add_read_info_raw(read_id,
                                                          [m.assigned_transcript for m in model_assignment.isoform_matches],
                                                          model_assignment.read_group)
                for m in model_assignment.isoform_matches:
                    self.transcript_read_ids[m.assigned_transcript].append(assignment)
            else:
                self.unused_reads.append(read_id)

    def correct_novel_transcrip_ends(self, transcript_model, assigned_reads):
        logger.debug("Verifying ends for transcript %s" % transcript_model.transcript_id)
        transcript_end = transcript_model.exon_blocks[-1][1]
        transcript_start = transcript_model.exon_blocks[0][0]
        start_supported = False
        read_starts = set()
        end_supported = False
        read_ends = defaultdict(int)

        for assignment in assigned_reads:
            read_exons = assignment.corrected_exons
            if abs(read_exons[0][0] - transcript_start) <= self.params.apa_delta:
                start_supported = True
            if not start_supported and read_exons[0][0] < transcript_model.exon_blocks[0][1]:
                read_starts.add(read_exons[0][0])
            if abs(read_exons[-1][1] - transcript_end) <= self.params.apa_delta:
                end_supported = True
            if not end_supported and read_exons[-1][1] > transcript_model.exon_blocks[-1][0]:
                read_ends[read_exons[-1][1]] += 1

        new_transcript_start = None
        if not start_supported:
            read_starts = sorted(read_starts)
            for read_start in read_starts:
                if read_start > transcript_start:
                    new_transcript_start = read_start
        if new_transcript_start and new_transcript_start < transcript_model.exon_blocks[0][1]:
            logger.debug("Changed start for transcript %s: from %d to %d" %
                         (transcript_model.transcript_id, transcript_model.exon_blocks[0][0], new_transcript_start))
            transcript_model.exon_blocks[0] = (new_transcript_start, transcript_model.exon_blocks[0][1])

        new_transcript_end = None
        if not end_supported:
            read_ends = sorted(read_ends)
            for read_end in read_ends:
                if read_end < transcript_end:
                    new_transcript_end = read_end
        if new_transcript_end and new_transcript_end > transcript_model.exon_blocks[-1][0]:
            logger.debug("Changed end for transcript %s: from %d to %d" %
                         (transcript_model.transcript_id, transcript_model.exon_blocks[-1][1], new_transcript_end))
            transcript_model.exon_blocks[-1] = (transcript_model.exon_blocks[-1][0], new_transcript_end)


class IntronPathStorage:
    def __init__(self, params, path_processor):
        self.params = params
        self.path_processor = path_processor
        self.intron_graph = path_processor.intron_graph
        self.paths = defaultdict(int)
        self.fl_paths = set()
        self.paths_to_reads = defaultdict(list)

    def fill(self, read_assignments):
        for a in read_assignments:
            if a.multimapper:
                continue
            intron_path = self.path_processor.thread_introns(a.corrected_introns)
            if not intron_path:
                continue
            read_end = a.corrected_exons[-1][1]
            is_end_trusted = a.strand == '+' and \
                             (a.polya_info.external_polya_pos != -1 or
                              a.polya_info.internal_polya_pos != -1)
            terminal_vertex = self.path_processor.thread_ends(intron_path[-1], read_end, is_end_trusted)
            if terminal_vertex:
                intron_path.append(terminal_vertex)
                
            read_start = a.corrected_exons[0][0]
            is_start_trusted = a.strand == '-' and \
                               (a.polya_info.external_polyt_pos != -1 or
                                a.polya_info.internal_polyt_pos != -1)
            starting_vertex = self.path_processor.thread_starts(intron_path[0], read_start, is_start_trusted)
            if starting_vertex:
                intron_path = [starting_vertex] + intron_path

            path_tuple = tuple(intron_path)
            self.paths[path_tuple] += 1
            if terminal_vertex and starting_vertex:
                if not self.params.needs_polya_for_construction or\
                        (terminal_vertex[0] == VERTEX_polya or starting_vertex[0] == VERTEX_polyt):
                    self.fl_paths.add(path_tuple)
            self.paths_to_reads[path_tuple].append(a)


class IntronPathProcessor:
    def __init__(self, params, intron_graph):
        self.params = params
        self.intron_graph = intron_graph
        self.all_vertices = set()
        self.all_vertices.update(self.intron_graph.intron_collector.clustered_introns.keys())
        for edge_set in self.intron_graph.outgoing_edges.values():
            self.all_vertices.update(edge_set)
        for edge_set in self.intron_graph.incoming_edges.values():
            self.all_vertices.update(edge_set)
        self.visited = set()

    def visit_vertex(self, v):
        if v in self.all_vertices:
            self.visited.add(v)

    def visit_path(self, p):
        for v in p:
            self.visit_vertex(v)

    def thread_introns(self, introns):
        path = []
        for intron in introns:
            if intron in self.intron_graph.intron_collector.discarded_introns:
                return None
            path.append(self.intron_graph.intron_collector.substitute(intron))
        return path

    def thread_ends(self, intron, end, trusted=False):
        possible_polyas = self.intron_graph.get_outgoing(intron, VERTEX_polya)
        if trusted:
            # find closes polyA
            for v in possible_polyas:
                if abs(v[1] - end) <= self.params.apa_delta:
                    return v

        outgoing_introns = self.intron_graph.get_outgoing(intron)
        if len(outgoing_introns) > 0:
            # intron has outgoing edges
            rightmost_exon_end = max([intron[0] for intron in outgoing_introns]) - 1
            if not trusted and end <= rightmost_exon_end + self.params.delta:
                # read end lies within next exon and has no polyA
                return None

        # consider all terminal position available for intron
        all_possible_ends = sorted(list(self.intron_graph.get_outgoing(intron, VERTEX_read_end)) +
                                   list(possible_polyas), key=lambda x:x[1])
        if len(all_possible_ends) == 0:
            return None

        rightmost_end = all_possible_ends[-1]
        if trusted and end >= rightmost_end[1] and rightmost_end[0] == VERTEX_read_end:
            # if we have trusted read, in cannot stop earlier that rightmost end (otherwise it should match polyA)
            return rightmost_end
        elif not trusted and end <= rightmost_end[1] + self.params.apa_delta and \
                (len(all_possible_ends) <= 1 or end > all_possible_ends[-2][1]):
            # non trusted should end before rightmost position + apa_delta but not earlier than second last
            return rightmost_end
        return None

    def thread_starts(self, intron, start, trusted=False):
        possible_polyas = self.intron_graph.get_incoming(intron, VERTEX_polyt)
        if trusted:
            # find closes polyT
            for v in possible_polyas:
                if abs(v[1] - start) <= self.params.apa_delta:
                    return v

        incoming_introns = self.intron_graph.get_incoming(intron)
        if len(incoming_introns) > 0:
            # intron has outgoing edges
            leftmost_exon_start = min([intron[1] for intron in incoming_introns]) + 1
            if not trusted and start >= leftmost_exon_start - self.params.delta:
                # read start lies within previous exon and has no polyA
                return None

        all_possible_starts = sorted(list(self.intron_graph.get_incoming(intron, VERTEX_read_start)) +
                                     list(possible_polyas), key=lambda x: x[1])
        if len(all_possible_starts) == 0:
            return None

        leftmost_start = all_possible_starts[0]
        if trusted and start <= leftmost_start[1] and leftmost_start[0] == VERTEX_read_start:
            return leftmost_start
        elif not trusted and start >= leftmost_start[1] and \
                (len(all_possible_starts) <= 1 or start < all_possible_starts[1][1]):
            return leftmost_start
        return None
