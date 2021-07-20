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

logger = logging.getLogger('IsoQuant')


class GraphBasedModelConstructor:
    transcript_id_counter = AtomicCounter()
    transcript_prefix = "transcript_"
    known_transcript_suffix = ".known"
    nic_transcript_suffix = ".nic"
    nnic_transcript_suffix = ".nnic"

    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params

        self.intron_graph = None
        self.path_processor = None
        self.path_storage = None
        self.known_isoforms = {}
        self.known_introns = {}
        self.mono_exon_isoforms = defaultdict(int)

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
        self.get_known_spliced_isoforms()
        self.construct_fl_isoforms()
        self.construct_monoexon_isoforms(read_assignment_storage)

        # split reads into clusters
        # self.construct_isoform_groups(read_assignment_storage)

    def get_known_spliced_isoforms(self):
        self.known_isoforms = {}
        for isoform_id in self.gene_info.all_isoforms_introns:
            isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
            intron_path = self.path_processor.thread_introns(isoform_introns)
            if not intron_path:
                logger.debug("No path founds for isoform: %s" % isoform_id)
                continue
            self.known_isoforms[tuple(intron_path)] = isoform_id
        self.known_introns = set(self.gene_info.intron_profiles.features)

    def construct_fl_isoforms(self):
        novel_cutoff = max(self.params.min_novel_count, self.params.min_novel_count_rel * self.intron_graph.max_coverage)
        for path in self.path_storage.fl_paths:
            # do not include terminal vertices
            intron_path = path[1:-1]
            transcript_range = (path[0][1], path[-1][1])
            if intron_path in self.known_isoforms:
                if self.path_storage.paths[path] < self.params.min_known_count:
                    continue
                transcript_type = TranscriptModelType.known
                id_suffix = self.known_transcript_suffix
                isoform_id = self.known_isoforms[intron_path]
                transcript_strand = self.gene_info.isoform_strands[isoform_id]
                transcript_gene = self.gene_info.gene_id_map[isoform_id]
            else:
                if self.path_storage.paths[path] < novel_cutoff:
                    continue
                isoform_id = "novel"
                transcript_strand = self.path_storage.get_path_strand(intron_path)
                transcript_gene = self.path_storage.get_path_gene(intron_path)
                if transcript_gene is None:
                    transcript_gene = self.select_reference_gene(transcript_range)
                    transcript_strand = self.gene_info.gene_strands[transcript_gene]
                if all(intron in self.known_introns for intron in intron_path):
                    transcript_type = TranscriptModelType.novel_in_catalog
                    id_suffix = self.nic_transcript_suffix
                else:
                    transcript_type = TranscriptModelType.novel_not_in_catalog
                    id_suffix = self.nnic_transcript_suffix
            new_transcript_id =  self.transcript_prefix + str(self.get_transcript_id()) + id_suffix

            novel_exons = get_exons(transcript_range, list(intron_path))
            new_model = TranscriptModel(self.gene_info.chr_id, transcript_strand,
                                        new_transcript_id, isoform_id, transcript_gene,
                                        novel_exons, transcript_type,
                                        additional_info="count %d" % self.path_storage.paths[path])
            self.transcript_model_storage.append(new_model)

    def select_reference_gene(self, transcript_rage):
        overlap_dict = {}
        gene_regions = self.gene_info.get_gene_regions()
        for gene_id in gene_regions.keys():
            overlap_dict[gene_id] = read_coverage_fraction([transcript_rage], [gene_regions[gene_id]])
        return get_top_count(overlap_dict)

    # group reads by the isoforms and modification events
    def construct_monoexon_isoforms(self, read_assignment_storage):
        logger.debug("Constructing isoform groups")
        self.mono_exon_isoforms = defaultdict(int)

        for read_assignment in read_assignment_storage:
            if not read_assignment:
                continue
            if read_assignment.assignment_type in {ReadAssignmentType.unique,
                                                   ReadAssignmentType.unique_minor_difference} and \
                    any(e.event_type == MatchEventSubtype.mono_exon_match for e in read_assignment.isoform_matches[0].match_subclassifications):
                self.mono_exon_isoforms[read_assignment.isoform_matches[0].assigned_transcript] += 1

        for isoform_id in self.mono_exon_isoforms:
            count = self.mono_exon_isoforms[isoform_id]
            if count < self.params.min_known_count:
                continue
            logger.debug("Adding known monoexon isoform %s" % isoform_id)
            self.transcript_model_storage.append(self.transcript_from_reference(isoform_id, count))

    # create transcript model object from reference isoforms
    def transcript_from_reference(self, isoform_id, count=0):
        new_transcript_id = self.transcript_prefix + str(self.get_transcript_id()) + self.known_transcript_suffix
        return TranscriptModel(self.gene_info.chr_id, self.gene_info.isoform_strands[isoform_id],
                               new_transcript_id, isoform_id, self.gene_info.gene_id_map[isoform_id],
                               self.gene_info.all_isoforms_exons[isoform_id], TranscriptModelType.known,
                               additional_info="count %d" % count)


class IntronPathStorage:
    def __init__(self, params, path_processor):
        self.params = params
        self.path_processor = path_processor
        self.intron_graph = path_processor.intron_graph
        self.paths = defaultdict(int)
        self.intron_strands = defaultdict(str)
        self.intron_genes = defaultdict(str)
        self.fl_paths = set()

    def fill(self, read_assignments):
        intron_strands = defaultdict(lambda: defaultdict(int))
        intron_genes = defaultdict(lambda: defaultdict(int))
        for a in read_assignments:
            intron_path = self.path_processor.thread_introns(a.corrected_introns)
            if not intron_path:
                continue
            for intron in intron_path:
                intron_strands[intron][a.strand] += 1
                if a.isoform_matches:
                    intron_genes[intron][a.isoform_matches[0].assigned_gene] += 1

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

        self.intron_strands = get_best_from_count_dicts(intron_strands)
        self.intron_genes = get_best_from_count_dicts(intron_genes)
        for p in self.paths.keys():
            is_fl = "FL" if p in self.fl_paths else "NO"
            logger.debug("%s path: %s: %d, %s, %s" %
                         (is_fl, str(p), self.paths[p], self.get_path_strand(p), self.get_path_gene(p)))

    def get_path_strand(self, intron_path):
        return get_collective_property(intron_path, self.intron_strands)

    def get_path_gene(self, intron_path):
        return get_collective_property(intron_path, self.intron_genes)


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









