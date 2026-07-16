############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Intron-path threading over the intron graph: turns a read's corrected introns
# into a graph path (with terminal polyA/TSS/read-end vertices) and buckets the
# reads by path. Moved verbatim out of graph_based_model_construction.

from collections import defaultdict

from isoquant_lib.model_construction.intron_graph import TerminalVertex


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
                if not self.params.requires_polya_for_construction or\
                        (terminal_vertex[0] == TerminalVertex.polya or starting_vertex[0] == TerminalVertex.polyt):
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
        possible_polyas = self.intron_graph.get_outgoing(intron, TerminalVertex.polya)
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

        # consider all terminal position available for intron (bare read ends,
        # confirmed 5' TSS ends for '-' transcripts, and confirmed polyA ends)
        all_possible_ends = sorted(list(self.intron_graph.get_outgoing(intron, TerminalVertex.read_end)) +
                                   list(self.intron_graph.get_outgoing(intron, TerminalVertex.tss_right)) +
                                   list(possible_polyas), key=lambda x:x[1])
        if len(all_possible_ends) == 0:
            return None

        rightmost_end = all_possible_ends[-1]
        if trusted and end >= rightmost_end[1] and rightmost_end[0] == TerminalVertex.read_end:
            # if we have trusted read, in cannot stop earlier that rightmost end (otherwise it should match polyA)
            return rightmost_end
        elif not trusted and end <= rightmost_end[1] + self.params.apa_delta and \
                (len(all_possible_ends) <= 1 or end > all_possible_ends[-2][1]):
            # non trusted should end before rightmost position + apa_delta but not earlier than second last
            return rightmost_end
        return None

    def thread_starts(self, intron, start, trusted=False):
        possible_polyas = self.intron_graph.get_incoming(intron, TerminalVertex.polyt)
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

        all_possible_starts = sorted(list(self.intron_graph.get_incoming(intron, TerminalVertex.read_start)) +
                                     list(self.intron_graph.get_incoming(intron, TerminalVertex.tss_left)) +
                                     list(possible_polyas), key=lambda x: x[1])
        if len(all_possible_starts) == 0:
            return None

        leftmost_start = all_possible_starts[0]
        if trusted and start <= leftmost_start[1] and leftmost_start[0] == TerminalVertex.read_start:
            return leftmost_start
        elif not trusted and start >= leftmost_start[1] and \
                (len(all_possible_starts) <= 1 or start < all_possible_starts[1][1]):
            return leftmost_start
        return None
