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
        self.transcript_model_storage = []
        self.transcript_read_ids = defaultdict(set)
        self.unused_reads = []
        self.transcript_counts = defaultdict(float)
        self.fsm = defaultdict(list)
        self.consistent = []
        self.inconsistent = []
        self.used_reads = set()
        self.assigned_reads = set()

    def get_transcript_id(self):
        return GraphBasedModelConstructor.transcript_id_counter.increment()

    def process(self, read_assignment_storage):
        self.intron_graph = IntronGraph(self.params, self.gene_info, read_assignment_storage)
        self.path_processor = IntronPathProcessor(self.intron_graph)
        # split reads into clusters
        self.construct_isoform_groups(read_assignment_storage)

    # group reads by the isoforms and modification events
    def construct_isoform_groups(self, read_assignment_storage):
        logger.debug("Constructing isoform groups")
        self.fsm = defaultdict(list)
        self.consistent = []
        self.inconsistent = []

        for read_assignment in read_assignment_storage:
            if not read_assignment:
                continue
            if read_assignment.assignment_type in {ReadAssignmentType.unique,
                                                   ReadAssignmentType.unique_minor_difference} and \
                    MatchEventSubtype.fsm in read_assignment.isoform_matches[0].match_subclassifications:
                self.fsm[read_assignment.isoform_matches[0].assigned_transcript].append(read_assignment)
            elif read_assignment.assignment_type in {ReadAssignmentType.unique,
                                                     ReadAssignmentType.unique_minor_difference,
                                                     ReadAssignmentType.ambiguous}:
                self.consistent.append(read_assignment)
            else:
                self.inconsistent.append(read_assignment)

    def construct_known_isoforms(self):
        for isoform_id in self.fsm.keys():
            isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
            intron_path = self.path_processor.thread_introns(isoform_introns)
            if not intron_path:
                logger.debug("No path founds for isoform: %s" % isoform_id)
                continue



class IntronPathStorage:
    def __init__(self, params, path_processor):
        self.params = params
        self.path_processor = path_processor
        self.intron_graph = path_processor.intron_graph
        self.paths = defaultdict(int)
        self.fl_paths = set()

    def fill(self, read_assignments):
        for a in read_assignments:
            intron_path = self.path_processor.thread_introns(a.corrected_introns)
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


class IntronPathProcessor:
    def __init__(self, params, intron_graph):
        self.params = params
        self.intron_graph = intron_graph
        self.all_vertices = set()
        self.all_vertices.update(self.intron_graph.intron_collector.clustered_introns.keys())
        self.all_vertices.update(self.intron_graph.outgoing_edges.values())
        self.all_vertices.update(self.intron_graph.incoming_edges.values())
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









