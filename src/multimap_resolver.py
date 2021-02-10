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

            if len(classified_assignments[ReadAssignmentType.unique]) > 0:
                return self.merge(classified_assignments[ReadAssignmentType.unique])
            elif len(classified_assignments[ReadAssignmentType.unique_minor_difference]) > 0:
                return self.merge(classified_assignments[ReadAssignmentType.unique_minor_difference])
            elif len(classified_assignments[ReadAssignmentType.ambiguous]) > 0:
                return self.merge(classified_assignments[ReadAssignmentType.ambiguous])
            else:
                return self.merge(classified_assignments[ReadAssignmentType.inconsistent], is_inconsistent=True)
        else:
            raise ValueError("Unsupported multimap strategy")

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
