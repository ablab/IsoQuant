############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from common import concat_gapless_blocks

logger = logging.getLogger('10XProfile')

# TODO discuss how to correctly process 10x
class ReadMapping:
    def __init__(self):
        self.exon_counts = {}
        self.intron_counts = {}

    def add_read(self, alignment):
        blocks = concat_gapless_blocks(alignment.get_blocks(), alignment.cigartuples)
        for b in blocks:
            if b not in self.exon_counts:
                self.exon_counts[b] = 0
            self.exon_counts[b] += 1

    def get_blocks(self):
        return sorted(self.exon_counts.keys())

