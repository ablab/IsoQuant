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
            return ReadAssignment(read_id, ReadAssignmentType.ambiguous)
        elif self.strategy == MultimapResolvingStrategy.merge:
            return self.merge(list(filter(lambda x: x.assignment_type != ReadAssignmentType.empty, assignment_list)))
        elif self.strategy == MultimapResolvingStrategy.take_best:
            # TODO: improve take_best strategy
            classified_assignments = defaultdict(list)
            for ra in assignment_list:
                if ra.assignment_type == ReadAssignmentType.empty:
                    continue
                classified_assignments[ra.assignment_type.name].append(ra)

            if len(classified_assignments["unique"]) > 0:
                return self.merge(classified_assignments["unique"])
            elif len(classified_assignments["unique_minor_difference"]) > 0:
                return self.merge(classified_assignments["unique_minor_difference"])
            elif len(classified_assignments["contradictory"]) > 0:
                return self.merge(classified_assignments["contradictory"])
            else:
                return self.merge(classified_assignments["ambiguous"])
        else:
            raise ValueError("Unsupported multimap strategy")

    def merge(self, assignment_list):
        if len(assignment_list) == 0:
            return None
        read_id = assignment_list[0].read_id
        if len(assignment_list) == 1:
            return assignment_list[0]

        read_assignment = ReadAssignment(read_id, ReadAssignmentType.ambiguous)
        for ra in assignment_list:
            for match in ra.isoform_matches:
                read_assignment.add_match(match)

        return read_assignment
