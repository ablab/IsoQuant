############################################################################
# Copyright (c) 2020 Saint Petersburg State University
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


class LinkedProfileConstructor:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params
        # TODO Chi-Lam: cloud read profile (e.g. count/span of introns exons)


class LinkedReadAssigner:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params
        # TODO Chi-Lam: assigns clouds to known isoforms