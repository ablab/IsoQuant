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
            primary_unique = None
            consistent_assignments = set()
            inconsistent_assignments = set()
            primary_inconsistent = None
            for  i, a in enumerate(assignment_list):
                if a.assignment_type in [ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference] and \
                        not a.multimapper:
                    primary_unique = i
                elif a.assignment_type == ReadAssignmentType.inconsistent and not a.multimapper:
                    primary_inconsistent = i
                if a.assignment_type not in [ReadAssignmentType.inconsistent, ReadAssignmentType.noninformative]:
                    consistent_assignments.add(i)
                elif a.assignment_type == ReadAssignmentType.inconsistent:
                    inconsistent_assignments.add(i)

            if primary_unique is not None:
                # primary unique is found, rest is noninformative
                logger.debug("Primary unique assignment selected: %d" % assignment_list[primary_unique].start)
                return self.suspend_assignments(assignment_list, [primary_unique])

            if consistent_assignments:
                logger.debug("Merging %d consistent assignments " % len(consistent_assignments))
                return self.suspend_assignments(assignment_list, consistent_assignments, ReadAssignmentType.ambiguous)

            if primary_inconsistent is not None:
                logger.debug("Primary inconsistent assignment selected: %d" % assignment_list[primary_inconsistent].start)
                return self.suspend_assignments(assignment_list, [primary_inconsistent])

            logger.debug("Merging inconsistent from %d assignments" % len(inconsistent_assignments))
            return self.suspend_assignments(assignment_list, inconsistent_assignments, ReadAssignmentType.inconsistent)
        else:
            raise ValueError("Unsupported multimap strategy")

    def suspend_assignments(self, assignment_list, assignments_to_keep, new_type=None):
        for i in range(len(assignment_list)):
            if i in assignments_to_keep:
                if new_type is not None:
                    assignment_list[i].assignment_type = new_type
                    assignment_list[i].multimapper = True
                continue
            assignment_list[i].assignment_type = ReadAssignmentType.noninformative
        return assignment_list
