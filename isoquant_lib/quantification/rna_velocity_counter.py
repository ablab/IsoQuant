############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""RNA velocity (spliced / unspliced) per-cell gene counts.

One counter instance handles a single chromosome (when constructed inside a
worker) or the merged sample (when constructed at finalization). Reads are
bucketed per (cell barcode, gene) into a spliced and an unspliced tally, the
per-chromosome tallies are written as TSV fragments in :meth:`dump`, and
:meth:`finalize` turns the merged sample-level TSV into a velocyto-style
``.loom`` with ``spliced`` / ``unspliced`` layers.

Splicing status is currently derived from the read's assignment type: a read
consistent with some annotated isoform is spliced, an inconsistent one is
unspliced. See ``.claude/RNA_VELOCITY.md`` for the limitations of that proxy.
"""

import csv
import logging
import os
from typing import Dict, Optional, Tuple

import pandas as pd
from scipy import sparse

# NOTE: loompy is imported lazily inside create_loom(). Importing it at module
# load time pulls in h5py/HDF5 (and its OpenMP runtime) into the main process,
# which breaks the fork-based ProcessPoolExecutor used for per-chromosome
# processing ("fork() called from a process already using GNU OpenMP").

from isoquant_lib.assignment.isoform_assignment import (
    ReadAssignment,
    ReadAssignmentType,
)
from isoquant_lib.quantification.long_read_counter import AbstractCounter

logger = logging.getLogger('IsoQuant')

SPLICED_ASSIGNMENT_TYPES = frozenset((
    ReadAssignmentType.unique,
    ReadAssignmentType.unique_minor_difference,
    ReadAssignmentType.ambiguous,
))

# Note: inconsistent_non_intronic is deliberately absent -- it marks reads whose
# disagreement with the annotation is not intronic, so it carries no unspliced
# signal. inconsistent_genic / inconsistent_multigenic never appear here either:
# they are only ever set as gene_assignment_type, while assignment_type stays
# noninformative (see .claude/READ_ASSIGNMENT_LIFECYCLE.md).
UNSPLICED_ASSIGNMENT_TYPES = frozenset((
    ReadAssignmentType.inconsistent,
    ReadAssignmentType.inconsistent_ambiguous,
))

ACCEPTED_ASSIGNMENT_TYPES = SPLICED_ASSIGNMENT_TYPES | UNSPLICED_ASSIGNMENT_TYPES

VELOCITY_COLUMNS = ["cell_id", "gene_id", "spliced", "unspliced"]


class RNAVelocityCounter(AbstractCounter):
    """Per-(cell, gene) spliced / unspliced read tallies with loom export."""

    def __init__(self, args, output_prefix: str,
                 string_pools=None, group_index: int = 0) -> None:
        # Skip AbstractCounter.__init__ -- we don't want counts_file_name's
        # suffix machinery; the velocity TSV path is already the full name.
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

        # (int group id, int gene id) -> count. Interned ids are kept instead of
        # strings: a chromosome's worth of (cell, gene) pairs is millions of keys
        # in a single-cell run, and two ints are far cheaper than two strings.
        # Resolved back to strings once, in dump().
        self.spliced: Dict[Tuple[int, int], int] = {}
        self.unspliced: Dict[Tuple[int, int], int] = {}

    # -- AbstractCounter interface --------------------------------------------
    # No-ops: the velocity counter only consumes finished read assignments and
    # emits its own TSV; it takes no part in the confirmed-feature / unassigned
    # / unaligned flows.

    def add_confirmed_features(self, features) -> None:
        return

    def add_unassigned(self, read_assignment=None) -> None:
        return

    def add_unaligned(self, n_reads: int = 1) -> None:
        return

    def _get_group_id(self, read_assignment: ReadAssignment) -> Optional[int]:
        # A read with no barcode has no cell to be counted against; unlike the
        # feature counters (which fold such reads into group 0) velocity output
        # is per cell, so drop it rather than invent a cell.
        if self.ignore_read_groups or not read_assignment.read_group_ids:
            return None
        if self.group_index >= len(read_assignment.read_group_ids):
            return None
        group_id = read_assignment.read_group_ids[self.group_index]
        # read_group_to_ids stores -1 when a strategy yielded no value for the
        # read. StringPool.get_str indexes a list, so -1 would silently resolve
        # to whatever barcode happens to be last in the pool.
        if group_id < 0:
            return None
        return group_id

    @staticmethod
    def _get_gene_id(read_assignment: ReadAssignment) -> Optional[int]:
        # Inconsistent reads may carry a match with no gene at all (e.g. the
        # quick_mode exit of match_inconsistent), so take the first match that
        # actually names a gene, as the feature counters do.
        for m in read_assignment.isoform_matches:
            if m.assigned_gene_id is not None:
                return m.assigned_gene_id
        return None

    def add_read_info(self, read_assignment: Optional[ReadAssignment] = None) -> None:
        if read_assignment is None:
            return
        if read_assignment.assignment_type not in ACCEPTED_ASSIGNMENT_TYPES:
            return

        gene_id = self._get_gene_id(read_assignment)
        if gene_id is None:
            return
        group_id = self._get_group_id(read_assignment)
        if group_id is None:
            return

        key = (group_id, gene_id)
        if read_assignment.assignment_type in SPLICED_ASSIGNMENT_TYPES:
            self.spliced[key] = self.spliced.get(key, 0) + 1
        else:
            self.unspliced[key] = self.unspliced.get(key, 0) + 1

    def dump(self) -> None:
        all_keys = set(self.spliced.keys()) | set(self.unspliced.keys())

        # Write a one-line header for each per-chromosome fragment. merge_counts()
        # keeps the header of the first fragment and strips one line from every
        # other fragment (header_lines=1); without a header the first data row of
        # each subsequent fragment would be dropped during the merge.
        write_header = not os.path.exists(self.output_file) or os.path.getsize(self.output_file) == 0
        group_pool = self.string_pools.get_read_group_pool(self.group_index)
        gene_pool = self.string_pools.gene_pool
        with open(self.output_file, "a") as f:
            writer = csv.writer(f, delimiter='\t')
            if write_header:
                writer.writerow(["#" + VELOCITY_COLUMNS[0]] + VELOCITY_COLUMNS[1:])
            for key in all_keys:
                group_id, gene_id = key
                writer.writerow([group_pool.get_str(group_id), gene_pool.get_str(gene_id),
                                 self.spliced.get(key, 0), self.unspliced.get(key, 0)])
        self.spliced.clear()
        self.unspliced.clear()

    def create_loom(self, counts_df: pd.DataFrame) -> None:
        import loompy  # lazy import, see note at top of module

        cell_codes, unique_cells = pd.factorize(counts_df['cell_id'])
        gene_codes, unique_genes = pd.factorize(counts_df['gene_id'])
        shape = (len(unique_genes), len(unique_cells))

        # Duplicate (gene, cell) entries are summed by coo_matrix, so a gene seen
        # in more than one merged fragment still ends up with its total.
        spliced_matrix = sparse.coo_matrix(
            (counts_df['spliced'], (gene_codes, cell_codes)), shape=shape)
        unspliced_matrix = sparse.coo_matrix(
            (counts_df['unspliced'], (gene_codes, cell_codes)), shape=shape)

        row_attrs = {'Gene': unique_genes.to_numpy()}
        col_attrs = {'CellID': unique_cells.to_numpy()}

        loom_file = self.output_file + ".loom"
        # loompy.create() refuses to overwrite; remove any stale file from a
        # previous run (e.g. when re-running with --force).
        if os.path.exists(loom_file):
            os.remove(loom_file)
        # The main (unnamed) layer holds spliced counts, matching velocyto's
        # convention, and is repeated as an explicit 'spliced' layer so scVelo
        # finds it by name.
        loompy.create(loom_file, spliced_matrix, row_attrs, col_attrs)
        with loompy.connect(loom_file) as ds:
            ds.layers['spliced'] = spliced_matrix
            ds.layers['unspliced'] = unspliced_matrix
        logger.info("RNA velocity counts are stored in %s and %s", self.output_file, loom_file)

    def finalize(self, args=None) -> None:
        counts_df = pd.read_csv(self.output_file, sep='\t', comment='#', header=None,
                                names=VELOCITY_COLUMNS)
        if counts_df.empty:
            logger.info("No RNA velocity counts collected, skipping loom creation for %s",
                        self.output_file)
            return
        try:
            self.create_loom(counts_df)
        except ImportError:
            # Do not fail a finished run over an optional export: the TSV, which
            # holds the same counts, has already been written.
            logger.warning("loompy is not installed, skipping the RNA velocity loom export. "
                           "Counts remain available in %s", self.output_file)
