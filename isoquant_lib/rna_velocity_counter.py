import csv
import logging
import os

import pandas as pd
from scipy import sparse

# NOTE: loompy is imported lazily inside create_loom(). Importing it at module
# load time pulls in h5py/HDF5 (and its OpenMP runtime) into the main process,
# which breaks the fork-based ProcessPoolExecutor used for per-chromosome
# processing ("fork() called from a process already using GNU OpenMP").

from .isoform_assignment import (
    IsoformMatch,
    ReadAssignment,
    ReadAssignmentType,
)
from .long_read_counter import AbstractCounter

logger = logging.getLogger('IsoQuant')

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
        # Truncate any stale output from a previous run, like AbstractCounter.
        # Per-chr counters are only built when a chromosome is (re)processed, so
        # this does not clobber finished chromosomes on --resume.
        open(self.output_file, "w").close()


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

        # Write a one-line header for each per-chromosome fragment. merge_counts()
        # keeps the header of the first fragment and strips one line from every
        # other fragment (header_lines=1); without a header the first data row of
        # each subsequent fragment would be dropped during the merge.
        write_header = not os.path.exists(self.output_file) or os.path.getsize(self.output_file) == 0
        with open(self.output_file, "a") as f:
            writer = csv.writer(f, delimiter='\t')
            if write_header:
                writer.writerow(["#cell_id", "gene_id", "spliced", "unspliced"])
            for cell_id, gene_id in all_keys:
                s_count = self.spliced.get((cell_id, gene_id), 0)
                u_count = self.unspliced.get((cell_id, gene_id), 0)
                writer.writerow([cell_id, gene_id, s_count, u_count])
    

    def create_loom(self):
        import loompy  # lazy import, see note at top of module
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

        loom_file = self.output_file + ".loom"
        # loompy.create() refuses to overwrite; remove any stale file from a
        # previous run (e.g. when re-running with --force).
        if os.path.exists(loom_file):
            os.remove(loom_file)
        loompy.create(loom_file, spliced_matrix, row_attrs, col_attrs)
        with loompy.connect(loom_file) as ds:
            ds.layers['spliced'] = spliced_matrix
            ds.layers['unspliced'] = unspliced_matrix
        return
    

    def finalize(self, args=None):
        self.velocity_counts_df = pd.read_csv(self.output_file, sep='\t', comment='#', header=None,
                                              names=['cell_id', 'gene_id', 'spliced', 'unspliced'])
        if self.velocity_counts_df.empty:
            logger.info("No RNA velocity counts collected, skipping loom creation for %s", self.output_file)
            return
        self.create_loom()
        