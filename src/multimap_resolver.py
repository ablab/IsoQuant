############################################################################
# Copyright (c) 2022-2025 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import math
from enum import Enum

from .isoform_assignment import ReadAssignmentType
from .common import intersection_len


logger = logging.getLogger('IsoQuant')


class MultimapResolvingStrategy(Enum):
    ignore_multimapper = 0
    merge = 1
    take_best = 2


class MultimapResolver:
    duplicate_counter = 0

    def __init__(self, strategy):
        self.strategy = strategy

    # list of ShortReadAssignments
    def resolve(self, assignment_list):
        if assignment_list is None or len(assignment_list) <= 1:
            return assignment_list

        if self.strategy == MultimapResolvingStrategy.ignore_multimapper:
            for a in assignment_list:
                a.assignment_type = ReadAssignmentType.suspended
            return assignment_list

        elif self.strategy == MultimapResolvingStrategy.merge:
            return self.merge_assignments(assignment_list)

        elif self.strategy == MultimapResolvingStrategy.take_best:
            return self.select_best_assignment(assignment_list)

        else:
            raise ValueError("Unsupported multimap strategy")

    def select_best_assignment(self, assignment_list):
        logger.debug("Resolving read %s" % assignment_list[0].read_id)
        primary_unique = []
        consistent_assignments = []
        inconsistent_assignments = []
        primary_inconsistent = []
        noninformative = []

        for i, a in enumerate(assignment_list):
            if a.assignment_type.is_inconsistent():
                if not a.multimapper:
                    primary_inconsistent.append(i)
                inconsistent_assignments.append(i)
            elif a.assignment_type.is_consistent():
                consistent_assignments.append(i)
                if not a.multimapper and not a.assignment_type == ReadAssignmentType.ambiguous:
                    primary_unique.append(i)
            else:
                noninformative.append(i)

        if primary_unique:
            return self.filter_assignments(assignment_list, primary_unique)

        if consistent_assignments:
            return self.filter_assignments(assignment_list, consistent_assignments)

        if primary_inconsistent:
            return self.select_best_inconsistent(assignment_list, primary_inconsistent)

        if inconsistent_assignments:
            return self.select_best_inconsistent(assignment_list, inconsistent_assignments)

        if noninformative:
            return self.select_noninformative(assignment_list, noninformative)

        logger.warning("Unexpected assignments in multimap resolution %s" % assignment_list[0].read_id)
        return assignment_list[0:1]

    @staticmethod
    def select_best_inconsistent(assignment_list, inconsistent_assignments):
        if len(inconsistent_assignments) <= 1:
            return MultimapResolver.filter_assignments(assignment_list, inconsistent_assignments)

        assignment_scores = []
        for i in inconsistent_assignments:
            assignment_scores.append((assignment_list[i].penalty_score, i))
        best_score = min(assignment_scores, key=lambda x: x[0])[0]
        best_assignments = [x[1] for x in filter(lambda x: x[0] == best_score, assignment_scores)]
        return MultimapResolver.filter_assignments(assignment_list, best_assignments)

    @staticmethod
    def merge_assignments(assignment_list):
        informative_assignments = set()
        for i in range(len(assignment_list)):
            if assignment_list[i].assignment_type != ReadAssignmentType.noninformative:
                informative_assignments.add(i)

        return MultimapResolver.filter_assignments(assignment_list, informative_assignments)

    # select all given assignments from the list, makr rest as suspended
    # if selected assignments feature more than 1 gene/isoform - mark as ambiguous accordingly
    @staticmethod
    def filter_assignments(assignment_list, assignments_to_keep):
        assignments_to_keep = MultimapResolver.find_duplicates(assignment_list, assignments_to_keep)

        all_genes = set()
        all_isoforms = set()
        for i in assignments_to_keep:
            assignment = assignment_list[i]
            all_genes.update(assignment.genes)
            all_isoforms.update(assignment.isoforms)

        change_transcript_assignment_type = len(all_isoforms) > 1
        change_gene_assignment_type = len(all_genes) > 1

        assignments_to_keep_set = set(assignments_to_keep)
        for i in range(len(assignment_list)):
            assignment = assignment_list[i]
            if i in assignments_to_keep_set:
                if assignment.assignment_type.is_inconsistent():
                    ambiguity_assignment_type = ReadAssignmentType.inconsistent_ambiguous
                else:
                    ambiguity_assignment_type = ReadAssignmentType.ambiguous

                if change_transcript_assignment_type:
                    assignment.assignment_type = ambiguity_assignment_type
                    assignment.multimapper = True
                if change_gene_assignment_type:
                    assignment.gene_assignment_type = ambiguity_assignment_type
                    assignment.multimapper = True
            else:
                assignment.assignment_type = ReadAssignmentType.suspended
                assignment.gene_assignment_type = ReadAssignmentType.suspended

        return assignment_list

    # finds groups of duplicated assignments
    @staticmethod
    def find_duplicates(assignment_list, assignment_indices):
        if len(assignment_indices) <= 1:
            return assignment_indices

        selected_assignments = []
        discarded_duplicates = set()
        for i in range(len(assignment_indices)):
            index1 = assignment_indices[i]
            if index1 in discarded_duplicates:
                # this assignment was already discarded - no need to compare it
                continue
            # if it was not discarded earlier - keep anyway; next assignments will either duplicate it or won't
            selected_assignments.append(index1)

            for j in range(i + 1, len(assignment_indices)):
                # find duplicating assignments among the next onex
                index2 = assignment_indices[j]
                if index2 in discarded_duplicates:
                    # this assignment was already discarded - no need to compare it
                    continue

                if assignment_list[index1] == assignment_list[index2]:
                    # discard duplicating assignment with larger index
                    discarded_duplicates.add(index2)

                    # check if duplicates come from the same region,
                    # if they are - it means we have identical records in the BAM files, print a warning
                    if assignment_list[index1].genomic_region == assignment_list[index2].genomic_region:
                        MultimapResolver.duplicate_counter += 1
                        example_assignment = assignment_list[index1]
                        if MultimapResolver.duplicate_counter == 5:
                            logger.warning(
                                "More possible duplicates were detected but will not be printed. "
                                "Use --debug to see all of them in the isoquant.log."
                                "All duplicated entries will be ignored. "
                                "This maybe caused either by the same read having more than one alignment with "
                                "the same reference id, start and end coordinates, "
                                "or by duplicated read names in the original files. ")
                        if MultimapResolver.duplicate_counter >= 5:
                            logger.debug(
                                "Read %s seems to have duplicated alignment entries at %s: %d and will be ignored. "
                                % (example_assignment.read_id, example_assignment.chr_id, example_assignment.start))
                        else:
                            logger.info(
                                "Read %s seems to have duplicated alignment entries at %s: %d and will be ignored. "
                                % (example_assignment.read_id, example_assignment.chr_id, example_assignment.start))

        return selected_assignments

    @staticmethod
    def select_noninformative(assignment_list, assignment_indices):
        # triplets (overlap_length, genomic_region_start, index)
        overlap_index_list = []
        max_overlap_len = 0

        for i in assignment_indices:
            assignment = assignment_list[i]
            overlap_len = intersection_len(assignment.genomic_region, (assignment.start, assignment.end))
            max_overlap_len = max(overlap_len, max_overlap_len)
            overlap_index_list.append((overlap_len, assignment.genomic_region[0], i))

        # select assignment with the best overlap with genic region and lowest region start (for reproducibility)
        min_region_start = math.inf
        best_assignment = -1
        for info in overlap_index_list:
            if info[0] == max_overlap_len and info[1] < min_region_start:
                min_region_start = info[1]
                best_assignment = info[2]

        assert best_assignment != -1
        return MultimapResolver.filter_assignments(assignment_list, {best_assignment})

