############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from enum import Enum

from .isoform_assignment import ReadAssignmentType

logger = logging.getLogger('IsoQuant')


class MultimapResolvingStrategy(Enum):
    ignore_multimapper = 0
    merge = 1
    take_best = 2


class MultimapResolver:

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
        primary_unique = set()
        consistent_assignments = set()
        inconsistent_assignments = set()
        primary_inconsistent = set()

        for i, a in enumerate(assignment_list):
            if a.assignment_type.is_inconsistent():
                if not a.multimapper:
                    primary_inconsistent.add(i)
                inconsistent_assignments.add(i)
            elif a.assignment_type.is_consistent():
                consistent_assignments.add(i)
                if not a.multimapper and not a.assignment_type == ReadAssignmentType.ambiguous:
                    primary_unique.add(i)

        if primary_unique:
            return self.filter_assignments(assignment_list, primary_unique)

        if consistent_assignments:
            return self.filter_assignments(assignment_list, consistent_assignments)

        if primary_inconsistent:
            return self.select_best_inconsistent(assignment_list, primary_inconsistent)

        if inconsistent_assignments:
            return self.select_best_inconsistent(assignment_list, inconsistent_assignments)

        return assignment_list

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
    def filter_assignments(assignment_list, assignments_to_keep, mark_as_ambiguous=True):
        all_genes = set()
        all_isoforms = set()
        for i in assignments_to_keep:
            assignment = assignment_list[i]
            all_genes.update(assignment.genes)
            all_isoforms.update(assignment.isoforms)

        change_transcript_assignment_type = mark_as_ambiguous and len(all_isoforms) > 1
        change_gene_assignment_type = mark_as_ambiguous and len(all_genes) > 1

        for i in range(len(assignment_list)):
            assignment = assignment_list[i]
            if i in assignments_to_keep:
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
