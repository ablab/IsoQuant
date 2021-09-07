############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
from collections import namedtuple
from functools import reduce
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

    def __init__(self, gene_info, params, expressed_gene_info=None):
        self.gene_info = gene_info
        self.params = params
        self.expressed_gene_info = expressed_gene_info

        self.intron_graph = None
        self.path_processor = None
        self.path_storage = None
        self.known_isoforms = {}
        self.known_introns = {}

        self.transcript_model_storage = []
        self.transcript_read_ids = defaultdict(set)
        self.unused_reads = []
        self.transcript_counts = defaultdict(float)


    def get_transcript_id(self):
        return GraphBasedModelConstructor.transcript_id_counter.increment()

    def process(self, read_assignment_storage):
        self.intron_graph = IntronGraph(self.params, self.gene_info, read_assignment_storage)
        self.path_processor = IntronPathProcessor(self.params, self.intron_graph)
        self.path_storage = IntronPathStorage(self.params, self.path_processor)
        self.path_storage.fill(read_assignment_storage)
        self.known_isoforms, self.known_introns = self.get_known_spliced_isoforms(self.gene_info)
        if not self.expressed_gene_info:
            self.expressed_isoforms = {}
            self.expressed_introns = {}
        else:
            self.expressed_isoforms, self.expressed_introns = self.get_known_spliced_isoforms(self.expressed_gene_info, "expressed")
        self.expressed_detected_set = set()
        self.visited_introns = set()

        self.construct_fl_isoforms()
        self.construct_monoexon_isoforms(read_assignment_storage)
        self.assign_reads_to_models(read_assignment_storage)

        for intron_chain, isoform_id in self.expressed_isoforms.items():
            if isoform_id in self.expressed_detected_set:
                continue
            else:
                logger.debug("++ Isoform %s was not detected: %s" % (isoform_id, str(intron_chain)))

        logger.debug("<< UNVISITED VERTICES <<")
        for v, count in self.intron_graph.intron_collector.clustered_introns.items():
            if v not in self.visited_introns:
                logger.debug("<< Vertex %s: %d was NOT visited" % (v, count))
        # split reads into clusters
        # self.construct_isoform_groups(read_assignment_storage)

    def get_known_spliced_isoforms(self, gene_info, s="known"):
        known_isoforms = {}
        for isoform_id in gene_info.all_isoforms_introns:
            isoform_introns = gene_info.all_isoforms_introns[isoform_id]
            intron_path = self.path_processor.thread_introns(isoform_introns)
            if not intron_path:
                logger.debug("== No path found for %s isoform %s: %s" % (s, isoform_id, gene_info.all_isoforms_introns[isoform_id]))
                continue
            else:
                logger.debug("== Path found for %s isoform %s: %s" % (s, isoform_id, gene_info.all_isoforms_introns[isoform_id]))
            known_isoforms[tuple(intron_path)] = isoform_id
        known_introns = set(gene_info.intron_profiles.features)
        return known_isoforms, known_introns

    def construct_fl_isoforms(self):
        # TODO refactor
        novel_cutoff = max(self.params.min_novel_count, self.params.min_novel_count_rel * self.intron_graph.max_coverage)
        added_known_isoforms = set()

        for path in self.path_storage.fl_paths:
            # do not include terminal vertices
            intron_path = path[1:-1]
            transcript_range = (path[0][1], path[-1][1])
            known_path = intron_path in self.known_isoforms
            if known_path:
                if self.path_storage.paths[path] < self.params.min_known_count:
                    continue
                transcript_type = TranscriptModelType.known
                id_suffix = self.known_transcript_suffix
                isoform_id = self.known_isoforms[intron_path]
                transcript_strand = self.gene_info.isoform_strands[isoform_id]
                transcript_gene = self.gene_info.gene_id_map[isoform_id]
                logger.debug("uuu Adding known spliced isoform %s" % isoform_id)
            else:
                if self.path_storage.paths[path] < novel_cutoff:
                    continue
                isoform_id = "novel"
                transcript_strand = self.path_storage.get_path_strand(path)
                transcript_gene = self.path_storage.get_path_gene(path)
                if transcript_gene is None:
                    transcript_gene = self.select_reference_gene(transcript_range)
                    transcript_strand = self.gene_info.gene_strands[transcript_gene]
                if all(intron in self.known_introns for intron in intron_path):
                    transcript_type = TranscriptModelType.novel_in_catalog
                    id_suffix = self.nic_transcript_suffix
                else:
                    transcript_type = TranscriptModelType.novel_not_in_catalog
                    id_suffix = self.nnic_transcript_suffix
                logger.debug("uuu Adding novel spliced isoform")

            transcript_num = self.get_transcript_id()
            new_transcript_id =  self.transcript_prefix + str(transcript_num) + id_suffix
            novel_exons = get_exons(transcript_range, list(intron_path))
            count = self.path_storage.paths[path]

            new_model = None
            if not known_path:
                assigner = LongReadAssigner(self.gene_info, self.params)
                profile_constructor = CombinedProfileConstructor(self.gene_info, self.params)

                combined_profile = profile_constructor.construct_profiles(novel_exons, PolyAInfo(-1, -1, -1, -1), [])
                assignment = assigner.assign_to_isoform(new_transcript_id, combined_profile)
                # check that no serious contradiction occurs
                logger.debug("uuu Checking novel transcript %s: %s; assignment type %s" %
                             (new_transcript_id, str(novel_exons), str(assignment.assignment_type)))

                if assignment.assignment_type == ReadAssignmentType.unique:
                    ref_tid = assignment.isoform_matches[0].assigned_transcript
                    if ref_tid not in added_known_isoforms:
                        logger.debug("uuu Substituting with known isoform %s" % ref_tid)
                        new_model = self.transcript_from_reference(ref_tid, count, transcript_num)
                else:
                    logger.debug("uuu Adding new model %s" % new_transcript_id)
                    new_model = TranscriptModel(self.gene_info.chr_id, transcript_strand,
                                                new_transcript_id, isoform_id, transcript_gene,
                                                novel_exons, transcript_type,
                                                additional_info="count %d" % count)
            elif isoform_id not in added_known_isoforms:
                added_known_isoforms.add(isoform_id)
                logger.debug("uuu Adding with known isoform %s -> %s" % (isoform_id, new_transcript_id))
                new_model = TranscriptModel(self.gene_info.chr_id, transcript_strand,
                                            new_transcript_id, isoform_id, transcript_gene,
                                            novel_exons, transcript_type,
                                            additional_info="count %d" % count)
            logger.debug("uuu %s: %s" % (new_transcript_id, str(novel_exons)))
            if new_model:
                self.transcript_model_storage.append(new_model)

            if intron_path in self.expressed_isoforms:
                ref_id = self.expressed_isoforms[intron_path]
                self.expressed_detected_set.add(ref_id)
                logger.debug("## Isoform %s matches reference chain %s, count = %d" % (new_transcript_id, ref_id, count))
            elif self.expressed_gene_info and isoform_id in self.expressed_gene_info.all_isoforms_exons.keys():
                self.expressed_detected_set.add(isoform_id)
                logger.debug("## Isoform %s is the reference %s, count = %d" % (new_transcript_id, isoform_id, count))
            else:
                logger.debug("## Isoform %s does match a reference chain, count = %d " % (new_transcript_id, count))

            if not known_path and len(novel_exons) == 2:
                # novel single intron transcrtipt
                logger.debug("uuu Added single intron isoform")
                
            for v in path:
                self.visited_introns.add(v)

    def select_reference_gene(self, transcript_rage):
        overlap_dict = {}
        gene_regions = self.gene_info.get_gene_regions()
        for gene_id in gene_regions.keys():
            overlap_dict[gene_id] = read_coverage_fraction([transcript_rage], [gene_regions[gene_id]])
        return get_top_count(overlap_dict)

    # group reads by the isoforms and modification events
    def construct_monoexon_isoforms(self, read_assignment_storage):
        logger.debug("Constructing isoform groups")
        mono_exon_isoform_counts = defaultdict(int)
        mono_exon_isoform_coverage = {}

        for read_assignment in read_assignment_storage:
            if not read_assignment:
                continue
            if read_assignment.assignment_type in {ReadAssignmentType.unique,
                                                   ReadAssignmentType.unique_minor_difference} and \
                    any(e.event_type == MatchEventSubtype.mono_exon_match for e in read_assignment.isoform_matches[0].match_subclassifications):
                t_id = read_assignment.isoform_matches[0].assigned_transcript
                mono_exon_isoform_counts[t_id] += 1
                assert len(self.gene_info.all_isoforms_exons[t_id]) == 1
                t_range = self.gene_info.all_isoforms_exons[t_id][0]
                t_len = t_range[1] - t_range[0] + 1
                if t_id not in mono_exon_isoform_coverage:
                    mono_exon_isoform_coverage[t_id] = [0 for _ in range(t_len)]
                start = max(0, read_assignment.corrected_exons[0][0] - t_range[0])
                end = min(t_len, read_assignment.corrected_exons[-1][1] - t_range[0] + 1)
                for i in range(start, end):
                    mono_exon_isoform_coverage[t_id][i] = 1

        for isoform_id in mono_exon_isoform_counts:
            count = mono_exon_isoform_counts[isoform_id]
            coverage = float(mono_exon_isoform_coverage[isoform_id].count(1)) / float(len(mono_exon_isoform_coverage[isoform_id]))
            logger.debug(">> Transcript %s, count %d, coverage %.4f" % (isoform_id, count, coverage))
            if count < self.params.min_known_count or coverage < self.params.min_mono_exon_coverage:
                logger.debug(">> Will not add")
                if self.expressed_gene_info and isoform_id in self.expressed_gene_info.all_isoforms_exons.keys():
                    logger.debug("## But found in the reference set!")
                    self.expressed_detected_set.add(isoform_id)
                continue
            self.transcript_model_storage.append(self.transcript_from_reference(isoform_id, count))
            logger.debug(">> Adding known monoexon isoform %s, %s, count = %d: %s" %
                         (self.transcript_model_storage[-1].transcript_id, isoform_id,
                          count, str(self.gene_info.all_isoforms_exons[isoform_id])))

            if self.expressed_gene_info and isoform_id in self.expressed_gene_info.all_isoforms_exons.keys():
                logger.debug("## Found in the reference set")
                self.expressed_detected_set.add(isoform_id)
            else:
                logger.debug("## Cannot be found in the reference set")

    # create transcript model object from reference isoforms
    def transcript_from_reference(self, isoform_id, count=0, transcript_id=None):
        if not transcript_id:
            transcript_id = self.get_transcript_id()
        new_transcript_id = self.transcript_prefix + str(transcript_id) + self.known_transcript_suffix
        return TranscriptModel(self.gene_info.chr_id, self.gene_info.isoform_strands[isoform_id],
                               new_transcript_id, isoform_id, self.gene_info.gene_id_map[isoform_id],
                               self.gene_info.all_isoforms_exons[isoform_id], TranscriptModelType.known,
                               additional_info="count %d" % count)

    # assign reads back to constructed isoforms
    def assign_reads_to_models(self, read_assignments):
        if not self.transcript_model_storage:
            logger.debug("No transcripts were assigned")
            self.unused_reads = [a.read_id for a in read_assignments]
            return
        
        logger.debug("Verifying transcript models")
        self.transcript_read_ids = defaultdict(set)
        self.unused_reads = []
        self.transcript_counts = defaultdict(float)

        logger.debug("Creating aritificial GeneInfo from %d transcript models" % len(self.transcript_model_storage))
        transcript_model_gene_info = GeneInfo.from_models(self.transcript_model_storage, self.params.delta)
        assigner = LongReadAssigner(transcript_model_gene_info, self.params)
        profile_constructor = CombinedProfileConstructor(transcript_model_gene_info, self.params)

        for assignment in read_assignments:
            read_exons = assignment.corrected_exons
            read_id = assignment.read_id
            logger.debug("Checking read %s: %s" % (assignment.read_id, str(read_exons)))
            model_combined_profile = profile_constructor.construct_profiles(read_exons, assignment.polya_info, [])
            model_assignment = assigner.assign_to_isoform(assignment.read_id, model_combined_profile)
            # check that no serious contradiction occurs
            if model_assignment.assignment_type in [ReadAssignmentType.unique,
                                                    ReadAssignmentType.unique_minor_difference]:
                t_id = model_assignment.isoform_matches[0].assigned_transcript
                self.transcript_read_ids[t_id].add(read_id)
                self.transcript_counts[t_id] += 1.0
            elif model_assignment.assignment_type == ReadAssignmentType.ambiguous:
                # FIXME: add qunatification options
                total_matches = len(model_assignment.isoform_matches)
                for m in model_assignment.isoform_matches:
                    self.transcript_counts[m.assigned_transcript] += 1.0 / total_matches
            else:
                self.unused_reads.append(read_id)


class IntronPathStorage:
    def __init__(self, params, path_processor):
        self.params = params
        self.path_processor = path_processor
        self.intron_graph = path_processor.intron_graph
        self.paths = defaultdict(int)
        self.path_strands = defaultdict(str)
        self.path_genes = defaultdict(str)
        self.fl_paths = set()

    def fill(self, read_assignments):
        path_strands = defaultdict(lambda: defaultdict(int))
        path_genes = defaultdict(lambda: defaultdict(int))
        for a in read_assignments:
            intron_path = self.path_processor.thread_introns(a.corrected_introns)
            if not intron_path:
                continue
            read_end = a.corrected_exons[-1][1]
            is_end_trusted = a.strand == '+' and a.polyA_found
            terminal_vertex = self.path_processor.thread_ends(intron_path[-1], read_end, is_end_trusted)
            if terminal_vertex:
                intron_path.append(terminal_vertex)
                
            read_start = a.corrected_exons[0][0]
            is_start_trusted = a.strand == '-' and a.polyA_found
            starting_vertex = self.path_processor.thread_starts(intron_path[0], read_start, is_start_trusted)
            if starting_vertex:
                intron_path = [starting_vertex] + intron_path

            path_tuple = tuple(intron_path)
            self.paths[path_tuple] += 1
            if terminal_vertex and starting_vertex:
                self.fl_paths.add(path_tuple)

            path_strands[path_tuple][a.strand] += 1
            assigned_genes = set([m.assigned_gene for m in a.isoform_matches])
            for g in assigned_genes:
                if g:
                    path_genes[path_tuple][g] += 1

        self.set_path_info(path_strands, path_genes)

        for p in self.paths.keys():
            is_fl = "   FL" if p in self.fl_paths else "nonFL"
            logger.debug("%s path: %s: %d, %s, %s" %
                         (is_fl, str(p), self.paths[p], self.get_path_strand(p), self.get_path_gene(p)))

    def set_path_info(self, path_strands, path_genes):
        for path in self.paths.keys():
            logger.debug(str(path) + " = " + str(self.paths[path]) + " = " + str(path_strands[path]))
            associated_strands = sorted(path_strands[path].items(), key=lambda x: x[1], reverse=True)
            associated_genes = sorted(path_genes[path].items(), key=lambda x: x[1], reverse=True)

            selected_strand = associated_strands[0][0]
            if selected_strand == '.' or \
                    (len(associated_strands) > 1 and associated_strands[0][1] / associated_strands[1][1] < 2):
                # two similar strands appear
                if associated_genes:
                    selected_strand = self.intron_graph.gene_info.gene_strands[associated_genes[0][0]]
            self.path_strands[path] = selected_strand

            for g in associated_genes:
                logger.debug(str(g))
                gene_id = g[0]
                if self.intron_graph.gene_info.gene_strands[gene_id] == self.path_strands[path]:
                    self.path_genes[path] = gene_id
                    break
            if path not in self.path_genes:
                logger.debug(
                    "Did not find suitable gene for path " + str(path) + ", strand: " + self.path_strands[path])
                self.path_genes[path] = 'novel_gene'

    def get_path_strand(self, intron_path):
        return self.path_strands[intron_path]

    def get_path_gene(self, intron_path):
        return self.path_genes[intron_path]


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
        possible_polyas = self.intron_graph.get_outgoing(intron, VertexType.polya)
        if trusted:
            # find closes polyA
            for v in possible_polyas:
                if abs(v[1] - end) <= self.params.apa_delta:
                    return v

        # consider all terminal position available for intron
        all_possible_ends = sorted(list(self.intron_graph.get_outgoing(intron, VertexType.read_end)) +
                                   list(possible_polyas), key=lambda x:x[1])
        if len(all_possible_ends) == 0:
            return None
        elif len(all_possible_ends) == 1:
            return all_possible_ends[0]

        rightmost_end = all_possible_ends[-1]
        if trusted and end >= rightmost_end[1] and rightmost_end[0] == VertexType.read_end:
            # if we have trusted read, in cannot stop earlier that rightmost end (otherwise it should match polyA)
            return rightmost_end
        elif not trusted and end <= rightmost_end[1] and end > all_possible_ends[-2][1]:
            # non trusted should end before rightmost position but not earlier than second last
            return rightmost_end
        return None

    def thread_starts(self, intron, start, trusted=False):
        possible_polyas = self.intron_graph.get_incoming(intron, VertexType.polyt)
        if trusted:
            # find closes polyT
            for v in possible_polyas:
                if abs(v[1] - start) <= self.params.apa_delta:
                    return v

        all_possible_starts = sorted(
            list(self.intron_graph.get_incoming(intron, VertexType.read_start)) + list(possible_polyas),
            key=lambda x: x[1])
        if len(all_possible_starts) == 0:
            return None
        elif len(all_possible_starts) == 1:
            return all_possible_starts[0]

        leftmost_start = all_possible_starts[0]
        if trusted and start <= leftmost_start[1] and leftmost_start[0] == VertexType.read_start:
            return leftmost_start
        elif not trusted and start >= leftmost_start[1] and start < all_possible_starts[1][1]:
            return leftmost_start
        return None









