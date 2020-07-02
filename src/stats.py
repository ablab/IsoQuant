############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict


logger = logging.getLogger('IsoQuant')


class EnumStats:
    def __init__(self):
        self.stats_dict = defaultdict(int)

    # element must be Enum
    def add(self, element):
        self.stats_dict[element] += 1

    def print_start(self, header_string = ""):
        if header_string:
            logger.info(header_string)
        for e in sorted(self.stats_dict.keys(), key=lambda x:x.name):
            logger.info("%s: %d" % (e.name, self.stats_dict[e]))

