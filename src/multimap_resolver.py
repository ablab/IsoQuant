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

    def resolve(self, assignment_list):
        if assignment_list is None or len(assignment_list) <= 1:
            return assignment_list

        if self.strategy == MultimapResolvingStrategy.ignore_multimapper:
            return []
        elif self.strategy == MultimapResolvingStrategy.merge:
            informative_assignments = list(filter(lambda x: x.assignment_type != ReadAssignmentType.noninformative, assignment_list))
            if len(informative_assignments) == 1:
                return informative_assignments
            return self.merge(informative_assignments)
        elif self.strategy == MultimapResolvingStrategy.take_best:
            classified_assignments = defaultdict(list)
            for a in assignment_list:
                if a.assignment_type == ReadAssignmentType.noninformative:
                    continue
                classified_assignments[a.assignment_type].append(a)

            consistent_assignemnts = classified_assignments[ReadAssignmentType.unique] +\
                                     classified_assignments[ReadAssignmentType.unique_minor_difference]
            primary_unique_assignment = self.select_primary(consistent_assignemnts)
            if primary_unique_assignment is not None:
                # logger.debug("Primary unique assignment selected: %s" % primary_unique_assignment.isoform_matches[0].assigned_transcript)
                return [primary_unique_assignment]

            consistent_assignemnts += classified_assignments[ReadAssignmentType.ambiguous]
            if consistent_assignemnts:
                logger.debug("Merging %d consistent assignments " % (len(consistent_assignemnts)))
                return self.merge(consistent_assignemnts)

            inconsistent_assignments = classified_assignments[ReadAssignmentType.inconsistent]
            primary_inconsistent_assignment = self.select_primary(inconsistent_assignments)
            if primary_inconsistent_assignment is not None:
                primary_inconsistent_assignment.multimapper = True
                # logger.debug("Selected primary inconsistent %s" %                            primary_inconsistent_assignment.isoform_matches[0].assigned_transcript)
                return [primary_inconsistent_assignment]
            logger.debug("Merging inconsistent from %d assignments" % len(inconsistent_assignments))
            return self.merge(inconsistent_assignments, is_inconsistent=True)
        else:
            raise ValueError("Unsupported multimap strategy")

    def select_primary(self, assignment_list):
        primary_assignment = None
        for a in assignment_list:
            if not a.multimapper:
                primary_assignment = a
        return primary_assignment

    def merge(self, assignment_list, is_inconsistent=False):
        if not assignment_list:
            return assignment_list

        for a in assignment_list:
            # set all as multimappers
            a.multimapper = True
            a.assignment_type = ReadAssignmentType.inconsistent if is_inconsistent else ReadAssignmentType.ambiguous

        return assignment_list
