############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""RNA velocity (spliced / unspliced / ambiguous) per-cell gene counts.

One counter instance handles a single chromosome (when constructed inside a
worker) or the merged sample (when constructed at finalization). Reads are
bucketed per (cell barcode, gene) into spliced / unspliced / ambiguous tallies,
the per-chromosome tallies are written as TSV fragments in :meth:`dump`, and
:meth:`finalize` turns the merged sample-level TSV into a velocyto-style
``.loom`` with ``spliced`` / ``unspliced`` / ``ambiguous`` layers.

Splicing status is derived from **intronic evidence**, not from whether the read
agrees with the annotation: a read is unspliced when it lies inside an intron
(``genic_intron``) or retains one (an intron-retention match event), and spliced
otherwise. A read that looks retained against one isoform and mature against
another is ambiguous. See ``.claude/RNA_VELOCITY.md``.
"""

import csv
import logging
import os
from enum import Enum
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

# NOTE: loompy is imported lazily inside create_loom(). Importing it at module
# load time pulls in h5py/HDF5 (and its OpenMP runtime) into the main process,
# which breaks the fork-based ProcessPoolExecutor used for per-chromosome
# processing ("fork() called from a process already using GNU OpenMP").

from isoquant_lib.assignment.isoform_assignment import (
    IsoformMatch,
    MatchClassification,
    ReadAssignment,
    ReadAssignmentType,
)
# INTRON_RETENTION_EVENTS (shared with IntronRetentionCounter) is defined next to
# the counter classes; see the comment there for which events count and why.
from isoquant_lib.quantification.long_read_counter import AbstractCounter, INTRON_RETENTION_EVENTS
from isoquant_lib.utils.string_pools import UNASSIGNED_GROUP_ID, UNASSIGNED_GROUP_NAME

logger = logging.getLogger('IsoQuant')

# Reads that were matched against isoforms; their splicing status comes from the
# match events. inconsistent_non_intronic is included: its reads disagree with
# the annotation for non-intronic reasons, so they are spliced molecules.
MATCHED_ASSIGNMENT_TYPES = frozenset((
    ReadAssignmentType.unique,
    ReadAssignmentType.unique_minor_difference,
    ReadAssignmentType.ambiguous,
    ReadAssignmentType.inconsistent,
    ReadAssignmentType.inconsistent_ambiguous,
    ReadAssignmentType.inconsistent_non_intronic,
))

# A read lying inside an intron never reaches a transcript: assign_to_isoform's
# "EMPTY - intronic" branch routes it through assign_to_overlapping_genes, which
# leaves assignment_type == noninformative and records the gene(s) it overlaps in
# matches classed genic_intron. These are the pre-mRNA reads velocity is built
# on, so noninformative is accepted when -- and only when -- such a match is
# present. See .claude/READ_ASSIGNMENT_LIFECYCLE.md.
UNMATCHED_ASSIGNMENT_TYPES = frozenset((
    ReadAssignmentType.noninformative,
))

ACCEPTED_ASSIGNMENT_TYPES = MATCHED_ASSIGNMENT_TYPES | UNMATCHED_ASSIGNMENT_TYPES

# Matches that name a gene but say nothing about splicing: the read overlaps the
# gene body yet resembles no isoform, so there is no intron/exon structure to read
# a verdict off. Counting them either way would be a guess.
#
# MatchClassification.undefined is deliberately *not* here: it means the SQANTI-style
# classifier had nothing to say, not that the read resembles no isoform, and such a
# read still carries match events the verdict can be computed from. Treating it as
# undecidable would silently drop reads.
UNDECIDABLE_CLASSIFICATIONS = frozenset((
    MatchClassification.genic,
    MatchClassification.intergenic,
))


class SplicingStatus(Enum):
    spliced = 0
    unspliced = 1
    ambiguous = 2


VELOCITY_COLUMNS = ["cell_id", "gene_id", "spliced", "unspliced", "ambiguous"]

# Peak dense memory one loom column block may cost, per layer. The loom write
# holds four of these at a time, so ~256 MB, independent of the cell count.
LOOM_WINDOW_BYTES = 64 * 1024 ** 2


class RNAVelocityCounter(AbstractCounter):
    """Per-(cell, gene) spliced / unspliced / ambiguous read tallies with loom export."""

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
        # nothing is written here: the owner empties output_paths() before counting

        self.args = args
        self.string_pools = string_pools
        self.group_index = group_index

        # (int group id, int gene id) -> count. Interned ids are kept instead of
        # strings: a chromosome's worth of (cell, gene) pairs is millions of keys
        # in a single-cell run, and two ints are far cheaper than two strings.
        # Resolved back to strings once, in dump().
        self.spliced: Dict[Tuple[int, int], int] = {}
        self.unspliced: Dict[Tuple[int, int], int] = {}
        self.ambiguous: Dict[Tuple[int, int], int] = {}
        self._tallies = {SplicingStatus.spliced: self.spliced,
                         SplicingStatus.unspliced: self.unspliced,
                         SplicingStatus.ambiguous: self.ambiguous}

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
        # A read with no barcode has no cell to be counted against. The feature
        # counters give it the UNASSIGNED_GROUP_NAME bucket like any other group;
        # velocity output is per cell and 'NA' is not a cell, so it is dropped
        # instead. The id-shaped half of that policy lives here, the name-shaped
        # half in dump() -- see _is_unassigned_group.
        if self.ignore_read_groups or not read_assignment.read_group_ids:
            return None
        if self.group_index >= len(read_assignment.read_group_ids):
            return None
        group_id = read_assignment.read_group_ids[self.group_index]
        # Legacy sentinel, only reachable from assignments serialized before
        # read_group_to_ids started interning UNASSIGNED_GROUP_NAME. A negative id
        # indexes the pool list from the end, naming a real barcode.
        if group_id == UNASSIGNED_GROUP_ID:
            return None
        return group_id

    @staticmethod
    def _is_unassigned_group(cell_id: str) -> bool:
        """True for the bucket every grouper falls back to when it cannot place a read.

        BarcodeGrouper returns UNASSIGNED_GROUP_NAME whenever a read carries no
        barcode (e.g. it is absent from --barcoded_reads), and since the group-id
        unification that name is also what an unset value interns as. Either way it
        is a pseudo-cell pooling every unplaceable read, not a real one, so it is
        kept out of the counts and out of the loom.
        """
        return cell_id == UNASSIGNED_GROUP_NAME

    @staticmethod
    def _get_gene_id(read_assignment: ReadAssignment) -> Optional[int]:
        # Inconsistent reads may carry a match with no gene at all (e.g. the
        # quick_mode exit of match_inconsistent), so take the first match that
        # actually names a gene, as the feature counters do.
        for m in read_assignment.isoform_matches:
            if m.assigned_gene_id is not None:
                return m.assigned_gene_id
        return None

    @staticmethod
    def _match_status(match: IsoformMatch) -> Optional[SplicingStatus]:
        """Splicing verdict for one read-to-isoform match, or None if undecidable."""
        # Read lies inside an intron -- pre-mRNA, whatever else it looks like.
        if match.match_classification == MatchClassification.genic_intron:
            return SplicingStatus.unspliced
        if any(e.event_type in INTRON_RETENTION_EVENTS for e in match.match_subclassifications):
            return SplicingStatus.unspliced
        if match.match_classification in UNDECIDABLE_CLASSIFICATIONS:
            return None
        return SplicingStatus.spliced

    @classmethod
    def _read_status(cls, read_assignment: ReadAssignment) -> Optional[SplicingStatus]:
        """Aggregate per-match verdicts into one status for the read.

        A read that retains an intron against every isoform it matched is
        unspliced; one that matches only mature isoforms is spliced; one that
        does both is ambiguous -- it is genuinely compatible with either model,
        which is what velocyto's third category is for.
        """
        statuses: List[SplicingStatus] = []
        for match in read_assignment.isoform_matches:
            if match.assigned_gene_id is None:
                continue
            status = cls._match_status(match)
            if status is not None:
                statuses.append(status)
        if not statuses:
            return None
        if all(st == SplicingStatus.unspliced for st in statuses):
            return SplicingStatus.unspliced
        if all(st == SplicingStatus.spliced for st in statuses):
            return SplicingStatus.spliced
        return SplicingStatus.ambiguous

    def add_read_info(self, read_assignment: Optional[ReadAssignment] = None) -> None:
        if read_assignment is None:
            return
        if read_assignment.assignment_type not in ACCEPTED_ASSIGNMENT_TYPES:
            return

        status = self._read_status(read_assignment)
        if status is None:
            return
        gene_id = self._get_gene_id(read_assignment)
        if gene_id is None:
            return
        group_id = self._get_group_id(read_assignment)
        if group_id is None:
            return

        key = (group_id, gene_id)
        tally = self._tallies[status]
        tally[key] = tally.get(key, 0) + 1

    def dump(self) -> None:
        all_keys = set(self.spliced.keys()) | set(self.unspliced.keys()) | set(self.ambiguous.keys())

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
                cell_id = group_pool.get_str(group_id)
                # Filtering by name happens here rather than per read in
                # _get_group_id: this loop already resolves the name, so the
                # hot path stays an integer compare.
                if self._is_unassigned_group(cell_id):
                    continue
                writer.writerow([cell_id, gene_pool.get_str(gene_id),
                                 self.spliced.get(key, 0), self.unspliced.get(key, 0),
                                 self.ambiguous.get(key, 0)])
        for tally in self._tallies.values():
            tally.clear()

    def create_loom(self, counts_df: pd.DataFrame) -> None:
        import loompy  # lazy import, see note at top of module

        cell_codes, unique_cells = pd.factorize(counts_df['cell_id'])
        gene_codes, unique_genes = pd.factorize(counts_df['gene_id'])
        n_genes, n_cells = len(unique_genes), len(unique_cells)

        # Duplicate (gene, cell) entries are summed by coo_matrix, so a gene seen
        # in more than one merged fragment still ends up with its total. csc so
        # the column slicing below is cheap. int32 halves the dense block a
        # window costs; per-(cell, gene) read counts cannot approach 2^31.
        def to_csc(column: str) -> sparse.csc_matrix:
            return sparse.coo_matrix(
                (counts_df[column].to_numpy(dtype=np.int32), (gene_codes, cell_codes)),
                shape=(n_genes, n_cells), dtype=np.int32).tocsc()

        unique_layers = {'spliced': to_csc('spliced'), 'unspliced': to_csc('unspliced'),
                         'ambiguous': to_csc('ambiguous')}
        # The main (unnamed) layer repeats the spliced matrix, matching
        # velocyto's convention, and 'spliced' is written explicitly so scVelo
        # finds it by name.
        layer_sources = {'': 'spliced', 'spliced': 'spliced',
                         'unspliced': 'unspliced', 'ambiguous': 'ambiguous'}

        loom_file = self.output_file + ".loom"
        # loompy refuses to overwrite; remove any stale file from a previous run
        # (e.g. when re-running with --force).
        if os.path.exists(loom_file):
            os.remove(loom_file)

        # Written one column block at a time rather than via loompy.create().
        # Loom is a dense HDF5 format, and assigning a sparse matrix to a layer
        # densifies the whole thing in memory: loompy's own window calculation
        # (layer_manager.py, `1024**3 // 8 * m.shape[0]`) multiplies where it
        # means to divide, so the window always clamps to every remaining column.
        # Measured at 5k genes x 20k cells that is +1.4 GB peak RSS; a spatial
        # run at 30k genes x 500k cells would ask for ~60 GB. add_columns() is
        # loompy's streaming API and caps the cost at one window instead.
        window = max(1, LOOM_WINDOW_BYTES // (n_genes * np.dtype(np.int32).itemsize))
        row_attrs = {'Gene': unique_genes.to_numpy()}
        cell_ids = unique_cells.to_numpy()
        with loompy.new(loom_file) as ds:
            for start in range(0, n_cells, window):
                stop = min(start + window, n_cells)
                # '' and 'spliced' are the same matrix, so densify each block once.
                block = {name: m[:, start:stop].toarray() for name, m in unique_layers.items()}
                ds.add_columns({name: block[key] for name, key in layer_sources.items()},
                               {'CellID': cell_ids[start:stop]}, row_attrs=row_attrs)
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
