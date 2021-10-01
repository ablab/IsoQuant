############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
from enum import Enum, unique

from src.common import *

logger = logging.getLogger('IsoQuant')

VERTEX_polya = -10
VERTEX_read_end = -11
VERTEX_polyt = -20
VERTEX_read_start = -21


class IntronCollector:
    def __init__(self, gene_info, delta=0):
        self.gene_info = gene_info
        self.known_introns = set(self.gene_info.intron_profiles.features)
        self.delta = delta
        # all intron counts
        # clustered intron counts
        self.clustered_introns = defaultdict(int)
        # how introns were corrected after clustering
        self.intron_correction_map = {}
        self.discarded_introns = set()

    def collect_introns(self, read_assignments):
        all_introns = defaultdict(int)
        for assignment in read_assignments:
            if not assignment.corrected_introns or assignment.multimapper:
                continue
            for intron in assignment.corrected_introns:
                all_introns[intron] += 1
        return all_introns

    def construct_similar_intron_map(self, all_introns):
        ordered_introns = sorted(all_introns.keys())
        similar_intron_map = defaultdict(list)
        # create a dict of similar introns
        for i, intron in enumerate(ordered_introns):
            j = i + 1
            while j < len(ordered_introns) and abs(ordered_introns[j][0] - intron[0]) <= self.delta:
                if abs(ordered_introns[j][1] - intron[1]) <= self.delta:
                    similar_intron_map[intron].append(ordered_introns[j])
                    similar_intron_map[ordered_introns[j]].append(intron)
                j += 1
        return similar_intron_map

    def cluster_introns(self, all_introns, min_count):
        similar_intron_map = self.construct_similar_intron_map(all_introns)
        introns_sorted_by_counts = sorted([(v, k) for k, v in all_introns.items()], reverse=True)
        for count, intron in introns_sorted_by_counts:
            if intron in self.known_introns:
                # known intron is always added as trustworthy
                self.clustered_introns[intron] = count
            elif intron in similar_intron_map:
                # intron has a similar intron
                similar_introns = []
                for similar_intron in similar_intron_map[intron]:
                    if similar_intron in self.clustered_introns:
                        # if similar intron was already added to the cluster with the higher count
                        similar_introns.append((count, similar_intron))

                if similar_introns:
                    # take the best one as a substitute
                    substitute_intron = sorted(similar_introns, reverse=True)[0]
                    self.clustered_introns[substitute_intron[1]] += substitute_intron[0]
                    self.intron_correction_map[intron] = substitute_intron[1]
                else:
                    # no introns were found
                    self.clustered_introns[intron] = count
            elif count < min_count:
                self.discarded_introns.add(intron)
                continue
            else:
                self.clustered_introns[intron] = count

    def process(self, read_assignments, min_count):
        logger.debug("Processing introns")
        all_introns = self.collect_introns(read_assignments)
        self.cluster_introns(all_introns, min_count)

    def add_substitute(self, original_intron, substitute_intron):
        self.clustered_introns[substitute_intron] += self.clustered_introns[original_intron]
        del self.clustered_introns[original_intron]
        self.intron_correction_map[original_intron] = substitute_intron

    def discard(self, intron):
        self.discarded_introns.add(intron)
        if intron in self.clustered_introns:
            del self.clustered_introns[intron]

    def substitute(self, v):
        if v in self.intron_correction_map:
            v = self.intron_correction_map[v]
        return v

    def simplify_correction_map(self):
        all_introns = sorted(self.intron_correction_map.keys())
        to_remove = set()
        for intron in all_introns:
            subs = self.intron_correction_map[intron]
            if subs in self.discarded_introns:
                to_remove.add(intron)
                continue
            if subs not in self.intron_correction_map:
                continue
            while subs in self.intron_correction_map:
                subs = self.intron_correction_map[subs]
            if subs in self.discarded_introns:
                to_remove.add(intron)
                continue
            self.intron_correction_map[intron] = subs

        for intron in to_remove:
            self.discard(intron)
            del self.intron_correction_map[intron]


class IntronGraph:
    def __init__(self, params, gene_info, read_assignments):
        self.params = params
        self.gene_info = gene_info
        self.read_assignments = read_assignments

        self.incoming_edges = defaultdict(set)
        self.outgoing_edges = defaultdict(set)
        self.intron_collector = IntronCollector(gene_info, params.delta)
        self.max_coverage = 0

        self.starting_known_positions = defaultdict(list)
        self.terminal_known_positions = defaultdict(list)
        for t, introns in self.gene_info.all_isoforms_introns.items():
            if not introns:
                continue
            self.starting_known_positions[introns[0]].append(self.gene_info.all_isoforms_exons[t][0][0])
            self.terminal_known_positions[introns[-1]].append(self.gene_info.all_isoforms_exons[t][-1][1])

        logger.debug("Collecting introns for %s" % self.gene_info.gene_db_list[0].id)
        self.intron_collector.process(read_assignments, self.params.min_novel_intron_count)
        self.construct()
        self.simplify()
        self.attach_terminal_positions()

    def add_edge(self, v1, v2):
        if v1 in self.intron_collector.intron_correction_map:
            v1 = self.intron_collector.intron_correction_map[v1]
        if v2 in self.intron_collector.intron_correction_map:
            v2 = self.intron_collector.intron_correction_map[v2]
        self.outgoing_edges[v1].add(v2)
        self.incoming_edges[v2].add(v1)

    def is_isolated(self, v):
        return not self.outgoing_edges[v] and not self.incoming_edges[v]

    def get_outgoing(self, intron, v_type):
        res = []
        for v in self.outgoing_edges[intron]:
            if v[0] == v_type:
                res.append(v)
        return sorted(res)

    def get_incoming(self, intron, v_type):
        res = []
        for v in self.incoming_edges[intron]:
            if v[0] == v_type:
                res.append(v)
        return sorted(res)

    def is_dead_end(self, intron):
        return not any(x[0] > 0 for x in self.outgoing_edges[intron])

    def is_dead_start(self, intron):
        return not any(x[0] > 0 for x in self.incoming_edges[intron])

    # merge vertex to its substitute, remove if isolated
    def collapse_vertex(self, to_collapse, substitute_vertex):
        self.outgoing_edges[substitute_vertex].update(self.outgoing_edges[to_collapse])
        for i in self.outgoing_edges[to_collapse]:
            self.incoming_edges[i].remove(to_collapse)
            self.incoming_edges[i].add(substitute_vertex)

        self.incoming_edges[substitute_vertex].update(self.incoming_edges[to_collapse])
        for i in self.incoming_edges[to_collapse]:
            self.outgoing_edges[i].remove(to_collapse)
            self.outgoing_edges[i].add(substitute_vertex)

        self.intron_collector.add_substitute(to_collapse, substitute_vertex)

    def construct(self):
        logger.debug("Constructing graph for %s" % self.gene_info.gene_db_list[0].id)
        for assignment in self.read_assignments:
            if assignment.multimapper or any(intron in self.intron_collector.discarded_introns for intron in assignment.corrected_introns):
                continue

            for i in range(len(assignment.corrected_introns) - 1):
                intron1 = assignment.corrected_introns[i]
                intron2 = assignment.corrected_introns[i + 1]
                self.add_edge(intron1, intron2)

        self.max_coverage =  max(self.intron_collector.clustered_introns.values()) if self.intron_collector.clustered_introns else 0

    def simplify(self):
        logger.debug("Simplifying graph")
        # check all outgoing edges
        to_remove = set()
        logger.debug("Removing outgoing tips and bulges")
        for current_intron in sorted(self.outgoing_edges.keys()):
            out_introns = self.outgoing_edges[current_intron]
            substitute_dict = self.collapse_vertex_set(out_introns)
            for i in sorted(substitute_dict.keys()):
                if i in to_remove:
                    continue
                to_remove.add(i)
                self.collapse_vertex(i, substitute_dict[i])

        for i in to_remove:
            del self.outgoing_edges[i]
            del self.incoming_edges[i]
        to_remove.clear()

        # check all incoming edges
        logger.debug("Removing incoming tips and bulges")
        for current_intron in sorted(self.incoming_edges.keys()):
            inc_introns = self.incoming_edges[current_intron]
            substitute_dict = self.collapse_vertex_set(inc_introns)
            for i in sorted(substitute_dict.keys()):
                if i in to_remove:
                    continue
                to_remove.add(i)
                self.collapse_vertex(i, substitute_dict[i])

        for i in to_remove:
            del self.outgoing_edges[i]
            del self.incoming_edges[i]
        to_remove.clear()

        # check all isolated vertices
        isolated = set()
        logger.debug("Collapsing isolated introns")
        for intron in self.intron_collector.clustered_introns.keys():
            if self.is_isolated(intron):
                isolated.add(intron)
        substitute_dict = self.collapse_vertex_set(isolated)
        for i in sorted(substitute_dict.keys()):
            to_remove.add(i)
            self.collapse_vertex(i, substitute_dict[i])

        for i in to_remove:
            del self.outgoing_edges[i]
            del self.incoming_edges[i]
        to_remove.clear()
        # self.print_graph()

        # remove low covered isolated vertices
        logger.debug("Removing isolated introns")
        self.max_coverage = max(
            self.intron_collector.clustered_introns.values()) if self.intron_collector.clustered_introns else 0
        count_cutoff = max(self.params.min_novel_isolated_intron_abs,
                           self.max_coverage * self.params.min_novel_isolated_intron_rel)
        for intron in isolated:
            if intron not in self.intron_collector.clustered_introns or intron in self.intron_collector.known_introns:
                # already removed or known
                continue
            if self.intron_collector.clustered_introns[intron] < count_cutoff:
                self.intron_collector.discard(intron)
                to_remove.add(intron)

        for i in to_remove:
            del self.outgoing_edges[i]
            del self.incoming_edges[i]

        self.intron_collector.simplify_correction_map()

    def collapse_vertex_set(self, vertex_set):
        if len(vertex_set) <= 1:
            return {}

        approved_set = set()
        substitute_dict = {}
        # get counts
        vertex_counts = sorted([(self.intron_collector.clustered_introns[i], i) for i in vertex_set], reverse=True)
        for count, vertex in vertex_counts:
            similar_vertices = []
            for i in sorted(approved_set):
                start_dist = abs(i[0] - vertex[0])
                end_dist = abs(i[1] - vertex[1])
                if start_dist < self.params.graph_clustering_distance and \
                        end_dist < self.params.graph_clustering_distance and \
                        count < self.intron_collector.clustered_introns[i] * self.params.graph_clustering_ratio:
                    # found similar intron among approved
                    similar_vertices.append((start_dist + end_dist, i))

            if not similar_vertices:
                # no similar introns
                approved_set.add(vertex)
            else:
                # selecting the most similar intron among approved
                substitute_vertex = sorted(similar_vertices)[0][1]
                substitute_dict[vertex] = substitute_vertex

        return substitute_dict

    def attach_terminal_positions(self):
        logger.debug("Setting terminal positions paths for %s" % self.gene_info.gene_db_list[0].id)
        polya_ends, read_ends, polyt_starts, read_starts = self.collect_terminal_positions()

        for intron in sorted(self.intron_collector.clustered_introns):
            self.attach_transcpt_ends(intron, polya_ends, read_ends, read_end=True)
            self.attach_transcpt_ends(intron, polyt_starts, read_starts, read_end=False)

    def attach_transcpt_ends(self, intron, polya_confirmed_positions, read_terminal_positions, read_end=True):
        read_ends_cutoff = self.params.terminal_position_abs
        clustered_polyas = self.cluster_polya_positions(polya_confirmed_positions[intron])
        if clustered_polyas:
            read_ends_cutoff = max(read_ends_cutoff, max(clustered_polyas.values()) * self.params.terminal_position_rel)
            extra_end_positions = {}
            furthers_confirmed_position = max(clustered_polyas.keys()) if read_end else min(clustered_polyas.keys())
            for position, count in read_terminal_positions[intron].items():
                if read_end and position >= furthers_confirmed_position + self.params.apa_delta:
                    extra_end_positions[position] = count
                elif not read_end and position <= furthers_confirmed_position - self.params.apa_delta :
                    extra_end_positions[position] = count
        else:
            extra_end_positions = read_terminal_positions[intron]

        if read_end and intron in self.outgoing_edges and len(self.outgoing_edges[intron]) > 0:
            # intron has outgoing edges, hard cut off
            read_ends_cutoff = max(read_ends_cutoff, self.intron_collector.clustered_introns[intron] *
                                   self.params.terminal_internal_position_rel)
        elif not read_end and intron in self.incoming_edges and len(self.incoming_edges[intron]) > 0:
            # intron has incoming edges, hard cut off
            read_ends_cutoff = max(read_ends_cutoff, self.intron_collector.clustered_introns[intron] *
                                   self.params.terminal_internal_position_rel)

        terminal_positions = self.cluster_terminal_positions(extra_end_positions,
                                                             read_end=read_end,
                                                             cutoff=read_ends_cutoff)
        if read_end:
            if intron in self.terminal_known_positions:
                logger.debug("Annotated terminal positions: " + str(sorted(self.terminal_known_positions[intron])))
            logger.debug("PolyA terminal positions: " + str(sorted(clustered_polyas.keys())))
            logger.debug("Simple terminal positions: " + str(sorted(terminal_positions.keys())))
            for pos in clustered_polyas.keys():
                self.outgoing_edges[intron].add((VERTEX_polya, pos))
            for pos in terminal_positions.keys():
                self.outgoing_edges[intron].add((VERTEX_read_end, pos))
        else:
            if intron in self.starting_known_positions:
                logger.debug("Annotated terminal positions: " + str(sorted(self.starting_known_positions[intron])))
            else:
                logger.debug("No annotated terminal positions")
            logger.debug("PolyA terminal positions: " + str(sorted(clustered_polyas.keys())))
            logger.debug("Simple terminal positions: " + str(sorted(terminal_positions.keys())))
            for pos in clustered_polyas.keys():
                self.incoming_edges[intron].add((VERTEX_polyt, pos))
            for pos in terminal_positions.keys():
                self.incoming_edges[intron].add((VERTEX_read_start, pos))

    def collect_terminal_positions(self):
        polya_ends = defaultdict(lambda: defaultdict(int))
        read_ends = defaultdict(lambda: defaultdict(int))
        polyt_starts = defaultdict(lambda: defaultdict(int))
        read_starts = defaultdict(lambda: defaultdict(int))

        for assignment in self.read_assignments:
            if assignment.multimapper or not assignment.corrected_introns:
                continue
            if any(intron in self.intron_collector.discarded_introns for intron in
                   assignment.corrected_introns):
                continue

            starting_intron = self.intron_collector.substitute(assignment.corrected_introns[0])
            read_start = assignment.corrected_exons[0][0]
            polyt_detected = assignment.strand == '-' and \
                             (assignment.polya_info.external_polyt_pos != -1 or
                              assignment.polya_info.internal_polyt_pos != -1)
            if polyt_detected:
                polyt_starts[starting_intron][read_start] += 1
            elif not self.is_start_internal(starting_intron, read_start):
                read_starts[starting_intron][read_start] += 1

            terminating_intron = self.intron_collector.substitute(assignment.corrected_introns[-1])
            read_end = assignment.corrected_exons[-1][1]
            polya_detected = assignment.strand == '+' and \
                             (assignment.polya_info.external_polya_pos != -1 or
                              assignment.polya_info.internal_polya_pos != -1)
            if polya_detected:
                polya_ends[terminating_intron][read_end] += 1
            elif not self.is_end_internal(terminating_intron, read_end):
                read_ends[terminating_intron][read_end] += 1

        self.print_terminal(read_starts, "Read starts")
        self.print_terminal(polyt_starts, "Read polyT")
        self.print_terminal(read_ends, "Read ends")
        self.print_terminal(polya_ends, "Read polyA")

        return polya_ends, read_ends, polyt_starts, read_starts

    def cluster_polya_positions(self, position_dict):
        clustered_counts = {}
        if not position_dict:
            return clustered_counts

        while position_dict:
            best_pair = max(position_dict.items(), key=lambda x:x[1])
            top_position = best_pair[0]
            total_count = 0
            for pos in range(top_position - self.params.apa_delta, top_position + self.params.apa_delta + 1):
                if pos in position_dict:
                    total_count += position_dict[pos]
                    del position_dict[pos]
            clustered_counts[top_position] = total_count

        max_count = max(clustered_counts.values())
        if max_count == self.params.terminal_position_abs:
            return clustered_counts
        cutoff = max(max_count * self.params.terminal_position_rel, self.params.terminal_position_abs)
        result = {}
        for k,v in clustered_counts.items():
            if v >= cutoff:
                result[k] = v
        return result

    def cluster_terminal_positions(self, position_dict, read_end, cutoff=0):
        # simple solution for now
        if not position_dict:
            return {}
        total_count = sum(position_dict.values())
        pos = max(position_dict.keys()) if read_end else min(position_dict.keys())
        if total_count < cutoff:
            return {}
        return {pos: total_count}

    def is_start_internal(self, intron, read_start):
        is_internal = False
        # checking previous introns
        for inc in self.incoming_edges[intron]:
            if inc[1] - self.params.delta <= read_start:
                is_internal = True
                break
        return is_internal

    def is_end_internal(self, intron, read_end):
        is_internal = False
        # checking previous introns
        for out in self.outgoing_edges[intron]:
            if out[0] + self.params.delta >= read_end:
                is_internal = True
                break
        return is_internal

    def print_graph(self):
        logger.debug("Printing graph")
        logger.debug("Vertices: %d, substituted: %d, discarded: %d" % (len(self.intron_collector.clustered_introns),
                                                                       len(self.intron_collector.intron_correction_map),
                                                                       len(self.intron_collector.discarded_introns)))
        logger.debug("Outgoing: %d, incoming %d" % (len(self.incoming_edges), len(self.outgoing_edges)))
        logger.debug("Discarded")
        logger.debug(self.intron_collector.discarded_introns)
        logger.debug("Collected")
        logger.debug(self.intron_collector.clustered_introns)
        logger.debug("Substituted")
        logger.debug(self.intron_collector.intron_correction_map)
        for intron in sorted(self.intron_collector.clustered_introns.keys()):
            logger.debug("Intron %s, count %d -> %s" % (str(intron), self.intron_collector.clustered_introns[intron],
                                                        ",".join([str(x) for x in self.outgoing_edges[intron]])))

        for intron in sorted(self.intron_collector.clustered_introns.keys()):
            logger.debug("Intron %s, count %d <- %s" % (str(intron), self.intron_collector.clustered_introns[intron],
                                                        ",".join([str(x) for x in self.incoming_edges[intron]])))

    def print_terminal(self, pos_dict, s=""):
        logger.debug(s)
        for intron in pos_dict:
            logger.debug("Intron %s -> %s" % (str(intron),
                                              ",".join([str(x) + ":" + str(pos_dict[intron][x])
                                                        for x in sorted(pos_dict[intron].keys())])))
