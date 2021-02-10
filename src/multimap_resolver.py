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
        read_id = assignment_list[0].read_id
        if self.strategy == MultimapResolvingStrategy.ignore_multimapper:
            return None
        elif self.strategy == MultimapResolvingStrategy.merge:
            return self.merge(list(filter(lambda x: x.assignment_type != ReadAssignmentType.noninformative, assignment_list)))
        elif self.strategy == MultimapResolvingStrategy.take_best:
            # TODO: improve take_best strategy
            classified_assignments = defaultdict(list)
            for ra in assignment_list:
                if ra.assignment_type == ReadAssignmentType.noninformative:
                    continue
                classified_assignments[ra.assignment_type].append(ra)

            primary_unique_assignment = \
                self.select_primary(classified_assignments[ReadAssignmentType.unique] +
                                    classified_assignments[ReadAssignmentType.unique_minor_difference])
            if primary_unique_assignment is not None:
                logger.debug("Primary unique assignment selected: %s" % primary_unique_assignment.assigned_transcript)
                return primary_unique_assignment

            for assignment_type in [ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference,
                                    ReadAssignmentType.ambiguous]:
                if len(classified_assignments[assignment_type]) > 0:
                    assignments = classified_assignments[ReadAssignmentType.unique]
                    logger.debug("Merging %d assignments of type %s" % (len(assignments), str(assignment_type)))
                    return self.merge(assignments)
            return self.merge(classified_assignments[ReadAssignmentType.inconsistent], is_inconsistent=True)
        else:
            raise ValueError("Unsupported multimap strategy")

    def select_primary(self, assignment_list):
        primary_assignment = None
        for a in assignment_list:
            if not a.is_secondary:
                primary_assignment = a
        return primary_assignment

    def merge(self, assignment_list, is_inconsistent=False):
        if len(assignment_list) == 0:
            return None
        read_id = assignment_list[0].read_id
        if len(assignment_list) == 1:
            return assignment_list[0]

        event_type = ReadAssignmentType.inconsistent if is_inconsistent else ReadAssignmentType.ambiguous
        read_assignment = ReadAssignment(read_id, event_type)
        for ra in assignment_list:
            for match in ra.isoform_matches:
                read_assignment.add_match(match)

        return read_assignment
