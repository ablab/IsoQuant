############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from enum import Enum
from collections import defaultdict

from src.common import *
from src.isoform_assignment import *

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
                a.assignment_type = ReadAssignmentType.noninformative
            return assignment_list

        elif self.strategy == MultimapResolvingStrategy.merge:
            informative_assignments = set()
            for i in range(len(assignment_list)):
                if assignment_list[i].assignment_type != ReadAssignmentType.noninformative:
                    informative_assignments.add(i)
            if len(informative_assignments) == 1:
                return self.suspend_assignments(assignment_list, informative_assignments)
            return self.suspend_assignments(assignment_list, informative_assignments, ReadAssignmentType.ambiguous)

        elif self.strategy == MultimapResolvingStrategy.take_best:
            return self.select_best_assignment(assignment_list)
        else:
            raise ValueError("Unsupported multimap strategy")

    def select_best_assignment(self, assignment_list):
        logger.debug("Resolving read %s" % assignment_list[0].read_id)
        logger.debug("Read assignment types %s" % " ".join(a.assignment_type.name for a in assignment_list))
        primary_unique = set()
        consistent_assignments = set()
        inconsistent_assignments = set()
        primary_inconsistent = set()
        for i, a in enumerate(assignment_list):
            if a.assignment_type in [ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference] and \
                    not a.multimapper:
                primary_unique.add(i)
            elif a.assignment_type == ReadAssignmentType.inconsistent and not a.multimapper:
                primary_inconsistent.add(i)
            if a.assignment_type in [ReadAssignmentType.unique,
                                     ReadAssignmentType.unique_minor_difference,
                                     ReadAssignmentType.ambiguous]:
                consistent_assignments.add(i)
            elif a.assignment_type == ReadAssignmentType.inconsistent:
                inconsistent_assignments.add(i)

        if primary_unique:
            if len(primary_unique) > 1:
                # TODO silence warn
                logger.warning("Multiple primary unique " +
                               ",".join([assignment_list[i].gene_id for i in primary_unique]))
                return self.suspend_assignments(assignment_list, primary_unique, ReadAssignmentType.ambiguous)
            # primary unique is found, rest is noninformative
            logger.debug("Primary unique assignment selected: %s" %
                         " ".join([assignment_list[i].gene_id for i in primary_unique]))
            return self.suspend_assignments(assignment_list, primary_unique)

        if consistent_assignments:
            logger.debug("Merging %d consistent assignments " % len(consistent_assignments))
            return self.suspend_assignments(assignment_list, consistent_assignments, ReadAssignmentType.ambiguous)

        if primary_inconsistent:
            return self.select_best_inconsistent(assignment_list, primary_inconsistent)

        logger.debug("Merging inconsistent from %d assignments" % len(inconsistent_assignments))
        return self.select_best_inconsistent(assignment_list, inconsistent_assignments)

    def select_best_inconsistent(self, assignment_list, inconsistent_assignments):
        if len(inconsistent_assignments) > 1:
            logger.debug("= Multiple primary inconsistent " +
                         ",".join([assignment_list[i].gene_id for i in inconsistent_assignments]))
            assignment_scores = []
            for i in inconsistent_assignments:
                assignment_scores.append((assignment_list[i].score, i))
            best_score = min(assignment_scores, key=lambda x: x[0])[0]
            logger.debug("= Best score " + str(best_score))
            best_isoforms = [x[1] for x in filter(lambda x: x[0] == best_score, assignment_scores)]
            return self.suspend_assignments(assignment_list, best_isoforms,
                                            ReadAssignmentType.inconsistent if len(best_isoforms) > 1 else None)

        logger.debug("Primary inconsistent assignment selected")
        return self.suspend_assignments(assignment_list, inconsistent_assignments)

    def suspend_assignments(self, assignment_list, assignments_to_keep, new_type=None):
        for i in range(len(assignment_list)):
            if i in assignments_to_keep:
                if new_type is not None:
                    assignment_list[i].assignment_type = new_type
                    assignment_list[i].multimapper = True
                continue
            assignment_list[i].assignment_type = ReadAssignmentType.noninformative
        return assignment_list
