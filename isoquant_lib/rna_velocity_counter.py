import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import loompy
import scipy.sparse as sparse
from scipy import stats
from xgboost import XGBClassifier

from .isoform_assignment import (
    IsoformMatch,
    ReadAssignment,
    ReadAssignmentType,
)
from .long_read_counter import AbstractCounter
from .read_groups import AbstractReadGrouper
from .gene_info import FeatureInfo, GeneInfo

logger = logging.getLogger('IsoQuant')

from .long_read_counter import AbstractCounter

SPLICED_ASSIGNMENT_TYPES = [ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference, ReadAssignmentType.ambiguous]
UNSPLICED_ASSIGNMENT_TYPES = [ReadAssignmentType.inconsistent_ambiguous, ReadAssignmentType.inconsistent]

ACCEPTED_ASSIGNMENT_TYPES = [ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference, ReadAssignmentType.ambiguous, ReadAssignmentType.inconsistent_ambiguous, ReadAssignmentType.inconsistent]

class RNAVelocityCounter(AbstractCounter):
    def __init__(self, args, output_prefix: str,
                 string_pools=None, group_index: int = 0) -> None:


        self.ignore_read_groups = string_pools is None
        self.output_prefix = output_prefix
        self.output_file = output_prefix
        self.output_counts_file_name = output_prefix
        self.output_tpm_file_name = None
        self.output_stats_file_name = None
        self.usable_file_name = None
        

        self.args = args
        self.string_pools = string_pools
        self.group_index = group_index
        self.spliced_col = []
        self.spliced_row = []
        self.spliced_val = []

        self.unspliced_col = []
        self.unspliced_row = []
        self.unspliced_val = []

        self.gene_ids = {}
        self.cell_ids = {}

        
    
    def add_read_info(self, read_assignment: ReadAssignment):
        if read_assignment is None:
            return
        if read_assignment.assignment_type not in ACCEPTED_ASSIGNMENT_TYPES:
            return

        isoform_match: IsoformMatch = read_assignment.isoform_matches[0]
        
        gene_id = isoform_match.assigned_gene
        group_id = read_assignment.read_group_ids[self.group_index]

        
        if gene_id not in self.gene_ids:
            self.gene_ids[gene_id] = len(self.gene_ids)

        if group_id not in self.cell_ids:
            self.cell_ids[group_id] = len(self.cell_ids)

        

        if read_assignment.assignment_type in SPLICED_ASSIGNMENT_TYPES:
            self.spliced_col.append(self.cell_ids[group_id])
            self.spliced_row.append(self.gene_ids[gene_id])
            self.spliced_val.append(1)
        else:
            self.unspliced_col.append(self.cell_ids[group_id])
            self.unspliced_row.append(self.gene_ids[gene_id])
            self.unspliced_val.append(1)

    
    def create_loom(self):
        self.total_genes = len(self.gene_ids)
        
        self.total_cells = len(self.cell_ids)


        sorted_group_ids = list(self.cell_ids.keys())
        pool = self.string_pools.get_read_group_pool(self.group_index)
        barcode_strings = [pool.get_str(gid) for gid in sorted_group_ids]

        spliced_matrix = sparse.coo_matrix(
                                            (self.spliced_val, (self.spliced_row, self.spliced_col)), 
                                            shape=(self.total_genes, self.total_cells)
                                        )
        spliced_matrix.sum_duplicates()
        
        unspliced_matrix = sparse.coo_matrix(
                                            (self.unspliced_val, (self.unspliced_row, self.unspliced_col)), 
                                            shape=(self.total_genes, self.total_cells)
                                        )
        unspliced_matrix.sum_duplicates()

        row_attrs = {'gene_id': np.array(list(self.gene_ids.keys()))}

        col_attrs = {'cell_id': np.array(barcode_strings)}

        loompy.create(self.output_file, spliced_matrix, row_attrs, col_attrs)
        with loompy.connect(self.output_file) as ds:
            ds.layers['spliced'] = spliced_matrix
            ds.layers['unspliced'] = unspliced_matrix



    def dump(self):
        return

    def finalize(self, args=None):
        self.create_loom()
        return
