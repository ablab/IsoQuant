import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import loompy
from scipy import sparse


from .isoform_assignment import (
    IsoformMatch,
    ReadAssignment,
    ReadAssignmentType,
)
from .long_read_counter import AbstractCounter
from .read_groups import AbstractReadGrouper
from .gene_info import FeatureInfo, GeneInfo
import csv

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
  
  


        self.spliced = {}
        self.unspliced = {}

        
    
    def add_read_info(self, read_assignment: ReadAssignment):
        if read_assignment is None:
            return
        if read_assignment.assignment_type not in ACCEPTED_ASSIGNMENT_TYPES:
            return

        isoform_match: IsoformMatch = read_assignment.isoform_matches[0]
        
        gene_id = isoform_match.assigned_gene

        group_id = read_assignment.read_group_ids[self.group_index]
        pool = self.string_pools.get_read_group_pool(self.group_index)
        cell_id = pool.get_str(group_id)

        if read_assignment.assignment_type in SPLICED_ASSIGNMENT_TYPES:
            if (cell_id, gene_id) not in self.spliced:
                self.spliced[(cell_id, gene_id)] = 1
            else:
                self.spliced[(cell_id, gene_id)] += 1
        else:

            if (cell_id, gene_id) not in self.unspliced:
                self.unspliced[(cell_id, gene_id)] = 1
            else:
                self.unspliced[(cell_id, gene_id)] += 1

    
    
    def dump(self):

        all_keys = set(self.spliced.keys()).union(set(self.unspliced.keys()))
    
        with open(self.output_file, "a") as f:
            writer = csv.writer(f, delimiter='\t')
            
            
            for cell_id, gene_id in all_keys:
                s_count = self.spliced.get((cell_id, gene_id), 0)
                u_count = self.unspliced.get((cell_id, gene_id), 0)
                
                writer.writerow([cell_id, gene_id, s_count, u_count])
                
        return       
    

    def create_loom(self):
        self.total_genes = self.velocity_counts_df['gene_id'].nunique()
        self.total_cells = self.velocity_counts_df['cell_id'].nunique()
        self.velocity_counts_df['cell_id_encoded'], self.unique_cells = pd.factorize(self.velocity_counts_df['cell_id'])
        self.velocity_counts_df['gene_id_encoded'], self.unique_genes = pd.factorize(self.velocity_counts_df['gene_id'])
        self.spliced_val = self.velocity_counts_df['spliced']
        self.spliced_row = self.velocity_counts_df['gene_id_encoded']
        self.spliced_col = self.velocity_counts_df['cell_id_encoded']

        self.unspliced_val = self.velocity_counts_df['unspliced']
        self.unspliced_row = self.velocity_counts_df['gene_id_encoded']
        self.unspliced_col = self.velocity_counts_df['cell_id_encoded']

        spliced_matrix = sparse.coo_matrix(
                                            (self.spliced_val, (self.spliced_row, self.spliced_col)), 
                                            shape=(self.total_genes, self.total_cells)
                                        )
        
        unspliced_matrix = sparse.coo_matrix(
                                            (self.unspliced_val, (self.unspliced_row, self.unspliced_col)), 
                                            shape=(self.total_genes, self.total_cells)
                                        )

        row_attrs = {'gene_id': self.unique_genes.to_numpy()}

        col_attrs = {'cell_id': self.unique_cells.to_numpy()}

        loompy.create(self.output_file + ".loom", spliced_matrix, row_attrs, col_attrs)
        with loompy.connect(self.output_file + ".loom") as ds:
            ds.layers['spliced'] = spliced_matrix
            ds.layers['unspliced'] = unspliced_matrix
        return
    

    def finalize(self, args=None):
        self.velocity_counts_df = pd.read_csv(self.output_file, sep='\t', header = None, names = ['cell_id', 'gene_id', 'spliced', 'unspliced'])
        self.create_loom()
        return
        