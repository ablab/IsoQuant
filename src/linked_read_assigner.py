############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import namedtuple

from src.isoform_assignment import *
from src.gene_info import *
from src.junction_comparator import *

logger = logging.getLogger('IsoQuant')

IsoformDiff = namedtuple("IsoformDiff", ("id", "diff"))


class LinkedReadProfileConstructor:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params
        # TODO Chi-Lam: cloud read profile (e.g. count/span of introns exons)

        
    def construct_intron_profile(self, intron_profile, read_blocks, boundaries):
        '''read_block is exon block; boundaries are intron boundaries'''

        blocks = list(map(lambda x: (x[0], x[1]+1), read_blocks))
        blocks = set(sum(blocks, ()))
        intersect = blocks.intersection(boundaries)
        if len(intersect) < 2:
            # within exon or new splice site
            blocks = sorted(list(blocks))
            tmp_boudaries = sorted(boundaries + blocks)
            block_idx = [tmp_boudaries.index(b) for b in blocks]
            if len(block_idx) > 2:
                print(block_idx)
                for idx in block_idx[::2]:
                    if int(idx / 2) <= len(intron_profile) - 1: 
                        intron_profile[int(idx / 2)] -= 1 # spanning introns
        elif len(intersect) >= 2:
            block_idx = sorted([boundaries.index(i) for i in intersect])
            print(block_idx)
            for i, idx in enumerate(block_idx):
                if i == len(block_idx):
                    break
                if idx % 2 == 0:
                    if i < len(block_idx) - 1 and block_idx[i+1] == idx + 1:
                        intron_profile[int(idx / 2)] += 1
                    else:
                        intron_profile[int(idx / 2)] -= 1
        return intron_profile




class LinkedReadAssigner:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params
        # TODO Chi-Lam: assigns clouds to known isoforms