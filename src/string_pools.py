############################################################################
# Copyright (c) 2025 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
String interning for memory optimization.

Replaces duplicated strings with integer indices referencing shared string pools.
Expected memory reduction: ~78% for billion-read datasets.

Architecture:
- StringPool: Bidirectional mapping between strings and integers
- StringPoolManager: Manages all pools for a worker
  - Global pools: gene, transcript (from annotation), chromosome (from chr_ids list)
  - Per-chromosome pools: barcode, UMI (from split files)
  - Read group pools: Multiple types depending on strategy

Usage:
- Workers build chromosome pool from chr_ids list
- Workers build gene/transcript pools from annotation database
- All IsoformMatch/ReadAssignment objects store integer IDs
- Properties provide transparent string access via pools
- Serialization writes integer format for disk savings

Implementation Details:
- All data classes (IsoformMatch, ReadAssignment, BasicReadAssignment)
  require string_pools parameter - no fallback to strings
- Pools are built deterministically (sorted order) for consistency
- Per-worker pools ensure no sharing/locking between parallel workers
"""

import os
import glob
import re
import logging
from typing import List, Dict

from .read_groups import AbstractReadGrouper, get_grouping_pool_types
from .assignment_loader import load_genedb
from .serialization import write_int, write_string, read_int, read_string

logger = logging.getLogger('IsoQuant')


class StringPool:
    """
    Bidirectional mapping between strings and integers for memory optimization.

    Each unique string is stored once and referenced by a unique integer ID.
    Provides O(1) lookup in both directions.
    """

    def __init__(self):
        """Initialize empty string pool."""
        self.str_to_int: Dict[str, int] = {}
        self.int_to_str: List[str] = []
        self.next_id: int = 0

    def add(self, string: str) -> int:
        """
        Add string to pool and return its integer ID.
        If string already exists, returns existing ID.

        Args:
            string: String to add to pool

        Returns:
            Integer ID for the string
        """
        if string in self.str_to_int:
            return self.str_to_int[string]

        int_id = self.next_id
        self.str_to_int[string] = int_id
        self.int_to_str.append(string)
        self.next_id += 1
        return int_id

    def get_int(self, string: str) -> int:
        """
        Get integer ID for string. Adds string if not present.

        Args:
            string: String to look up

        Returns:
            Integer ID for the string
        """
        return self.add(string)

    def get_str(self, int_id: int) -> str:
        """
        Get string for integer ID.

        Args:
            int_id: Integer ID to look up

        Returns:
            String corresponding to the ID

        Raises:
            IndexError: If int_id is out of range
        """
        return self.int_to_str[int_id]

    def __len__(self) -> int:
        """Return number of unique strings in pool."""
        return len(self.int_to_str)

    def __contains__(self, string: str) -> bool:
        """Check if string is in pool."""
        return string in self.str_to_int


class StringPoolManager:
    """
    Manages all string pools for a chromosome processing worker.

    Pools are categorized by data type:
    - Global pools: gene, transcript (from annotation), chromosome (from chr_ids list)
    - Per-chromosome pools: barcode, UMI (from split files)
    - Read group pools: Multiple types depending on grouping strategy
    """

    def __init__(self):
        """Initialize all string pools."""
        # Global pools
        self.gene_pool = StringPool()  # from annotation
        self.transcript_pool = StringPool()  # from annotation
        self.chromosome_pool = StringPool()  # from chr_ids list

        # Per-chromosome pools (from split files)
        self.barcode_pool = StringPool()
        self.umi_pool = StringPool()

        # Read group pools by type
        self.read_group_tsv_pools: Dict[str, StringPool] = {}  # pool_key -> pool (static, from TSV files)
        self.read_group_dynamic_pools: Dict[int, StringPool] = {}  # spec_index -> pool (dynamic, for tag/read_id)
        self.file_name_pool = StringPool()  # For file_name grouping (static)
        self.barcode_spot_pools: Dict[int, StringPool] = {}  # For barcode_spot grouping (multi-column support)

        # Mapping from group spec index to pool type ('file_name', 'barcode_spot', 'tsv:N', 'dynamic')
        self.group_spec_pool_types: Dict[int, str] = {}

    def build_from_gffutils(self, gffutils_db):
        """
        Build global pools from gffutils database.

        Pools are built in deterministic sorted order to ensure consistency
        across workers.

        Args:
            gffutils_db: gffutils.FeatureDB object
        """
        logger.debug("Building string pools from gffutils database")

        # Gene pool
        gene_ids = set()
        for gene in gffutils_db.features_of_type('gene'):
            gene_ids.add(gene.id)
        for gene_id in sorted(gene_ids):
            self.gene_pool.add(gene_id)
        logger.debug(f"Gene pool: {len(self.gene_pool)} unique genes")

        # Transcript pool
        transcript_ids = set()
        # Query different transcript feature types
        for feature_type in ['transcript', 'RNA', 'mRNA', 'miRNA', 'ncRNA', 'tRNA', 'rRNA',
                             'snRNA', 'lnc_RNA', 'snoRNA', 'telomerase_RNA', 'Y_RNA', 'scRNA']:
            try:
                for transcript in gffutils_db.features_of_type(feature_type):
                    transcript_ids.add(transcript.id)
            except:
                pass  # Feature type might not exist
        for transcript_id in sorted(transcript_ids):
            self.transcript_pool.add(transcript_id)
        logger.debug(f"Transcript pool: {len(self.transcript_pool)} unique transcripts")

    def build_from_gene_db(self, gene_db):
        """
        Build global pools from gene annotation database.

        Pools are built in deterministic sorted order to ensure consistency
        across workers.

        Args:
            gene_db: GeneInfo database from gene_info.py
        """
        logger.debug("Building string pools from gene database")

        # Gene pool
        gene_ids = sorted(gene_db.genes.keys())
        for gene_id in gene_ids:
            self.gene_pool.add(gene_id)
        logger.debug(f"Gene pool: {len(self.gene_pool)} unique genes")

        # Transcript pool
        transcript_ids = sorted(gene_db.transcripts.keys())
        for transcript_id in transcript_ids:
            self.transcript_pool.add(transcript_id)
        logger.debug(f"Transcript pool: {len(self.transcript_pool)} unique transcripts")

    def build_chromosome_pool(self, chr_ids):
        """
        Build chromosome pool from list of chromosome IDs.

        Args:
            chr_ids: List of chromosome IDs to process
        """
        # Sort alphabetically for deterministic order across workers
        for chr_id in sorted(chr_ids):
            self.chromosome_pool.add(chr_id)
        logger.debug(f"Chromosome pool: {len(self.chromosome_pool)} chromosomes")

    def build_file_name_pool(self, sample):
        """
        Build file name pool from sample file list.

        Uses same logic as FileNameGrouper: uses readable_names_dict (labels)
        if available, otherwise uses basename without extension.
        Pool is built in sorted order for deterministic IDs across workers.
        Also includes 'NA' for reads without filename info.

        Args:
            sample: Sample object with file_list and readable_names_dict attributes
        """
        # Collect all file names first
        file_names = set()

        # If readable_names_dict is available (e.g., when -l labels are used),
        # use its values - matching FileNameGrouper behavior
        if sample.readable_names_dict:
            for readable_name in sample.readable_names_dict.values():
                file_names.add(readable_name)
        else:
            # file_list is a list of lists (libraries), each library has one or more files
            for lib in sample.file_list:
                # Get basename without extension, matching FileNameGrouper logic
                readable_name = os.path.splitext(os.path.basename(lib[0]))[0]
                file_names.add(readable_name)

        # Add in sorted order for deterministic IDs
        for name in sorted(file_names):
            self.file_name_pool.add(name)
        logger.debug(f"File name pool: {len(self.file_name_pool)} unique readable names")

    def build_barcode_spot_pool(self, barcode2spot_spec: str, column_index: int = 0,
                                 file_column: int = 1):
        """
        Build barcode-to-spot pool for a specific column.

        Pool is built in sorted order for deterministic IDs across workers.
        Also includes 'NA' for reads without barcode mapping.

        Args:
            barcode2spot_spec: barcode2spot spec (may include :col:cols suffix)
            column_index: 0-based pool index (for storage in barcode_spot_pools dict)
            file_column: Actual file column to read (from parsed spec)
        """
        # Extract filename from spec
        bc2spot_file = barcode2spot_spec.split(':')[0]

        # Collect all spots for this column
        spots = set()
        with open(bc2spot_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) > file_column:
                    spots.add(parts[file_column])

        # Add default group ID for reads without barcode mapping
        spots.add(AbstractReadGrouper.default_group_id)

        # Create and populate pool
        pool = StringPool()
        for spot in sorted(spots):
            pool.add(spot)

        self.barcode_spot_pools[column_index] = pool
        logger.debug(f"Barcode spot pool (col {column_index}): {len(pool)} unique spots")

    def load_barcode_pool(self, barcode_file: str):
        """
        Load barcode pool from split barcode file.

        Pool is built in sorted order for deterministic IDs across workers.
        Also includes 'NA' for reads without barcode.

        Args:
            barcode_file: Path to split barcode file for chromosome
                         Format: read_id\tbarcode\tumi
        """
        if not os.path.exists(barcode_file):
            logger.debug(f"Barcode file not found: {barcode_file}")
            return

        # Collect all barcodes and UMIs first
        # Use split('\t') to preserve empty fields (empty UMI = consecutive tabs)
        barcodes = set()
        umis = set()
        with open(barcode_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2 and parts[1]:
                    barcodes.add(parts[1])
                if len(parts) >= 3 and parts[2]:
                    umis.add(parts[2])

        # Add default group ID for reads without barcode
        barcodes.add(AbstractReadGrouper.default_group_id)

        # Add in sorted order for deterministic IDs
        for barcode in sorted(barcodes):
            self.barcode_pool.add(barcode)
        for umi in sorted(umis):
            self.umi_pool.add(umi)

        logger.debug(f"Barcode pool: {len(self.barcode_pool)} unique barcodes")
        logger.debug(f"UMI pool: {len(self.umi_pool)} unique UMIs")

    def load_read_group_tsv_pool(self, read_group_file: str, pool_key: str, col_index: int = 1, delimiter: str = '\t'):
        """
        Load read group pool from split TSV/CSV file.

        Pool is built in sorted order for deterministic IDs across workers.

        Args:
            read_group_file: Path to split read_group file for chromosome
            pool_key: Unique key for this pool (e.g., 'spec2_col1' for multi-column files)
            col_index: Column index to read group values from (1-based, default 1)
            delimiter: Field delimiter (from file spec, defaults to tab)
        """
        if not os.path.exists(read_group_file):
            logger.debug(f"Read group file not found: {read_group_file}")
            return

        # Collect all group values first
        # Always include 'NA' (default group ID) in case some reads are not in the TSV file
        group_values = {AbstractReadGrouper.default_group_id}
        with open(read_group_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip()
                if not line:
                    continue
                parts = line.split(delimiter)
                if len(parts) > col_index:
                    group_values.add(parts[col_index])

        # Add in sorted order for deterministic IDs
        pool = StringPool()
        for value in sorted(group_values):
            pool.add(value)

        self.read_group_tsv_pools[pool_key] = pool
        logger.debug(f"Read group TSV pool ({pool_key}): {len(pool)} unique values")

    def set_group_spec_pool_type(self, spec_index: int, pool_type: str):
        """
        Set the pool type for a grouping specification index.

        Args:
            spec_index: Index of the grouping spec (position in read_group list)
            pool_type: One of 'file_name', 'barcode_spot', 'tsv:N', 'dynamic'
        """
        self.group_spec_pool_types[spec_index] = pool_type

    def get_read_group_pool(self, spec_index: int) -> StringPool:
        """
        Get the appropriate pool for a grouping spec index.

        Args:
            spec_index: Index of the grouping spec

        Returns:
            The StringPool to use for this grouping spec
        """
        if spec_index not in self.group_spec_pool_types:
            # Default to dynamic pool if not specified
            return self._get_or_create_dynamic_pool(spec_index)

        pool_type = self.group_spec_pool_types[spec_index]

        if pool_type == 'file_name':
            return self.file_name_pool
        elif pool_type.startswith('barcode_spot:'):
            # Multi-column barcode_spot: 'barcode_spot:COL_INDEX'
            col_idx = int(pool_type.split(':')[1])
            if col_idx in self.barcode_spot_pools:
                return self.barcode_spot_pools[col_idx]
            return self._get_or_create_dynamic_pool(spec_index)
        elif pool_type == 'barcode_spot':
            # Legacy single-column (use column 0, first spot column)
            if 0 in self.barcode_spot_pools:
                return self.barcode_spot_pools[0]
            return self._get_or_create_dynamic_pool(spec_index)
        elif pool_type == 'barcode':
            # Use barcode pool for direct barcode grouping
            return self.barcode_pool
        elif pool_type.startswith('tsv:'):
            # Extract pool_key from 'tsv:SPEC_INDEX:COL_INDEX:DELIMITER'
            # pool_key format is 'SPEC_INDEX:COL_INDEX'
            parts = pool_type.split(':')
            if len(parts) >= 3:
                pool_key = f"{parts[1]}:{parts[2]}"
            else:
                pool_key = parts[1]  # Fallback for old format 'tsv:N'
            if pool_key in self.read_group_tsv_pools:
                return self.read_group_tsv_pools[pool_key]
            else:
                # Pool not loaded yet, use dynamic pool for this spec
                return self._get_or_create_dynamic_pool(spec_index)
        else:  # 'dynamic' or unknown
            return self._get_or_create_dynamic_pool(spec_index)

    def _get_or_create_dynamic_pool(self, spec_index: int) -> StringPool:
        """Get or create a dynamic pool for a spec index."""
        if spec_index not in self.read_group_dynamic_pools:
            self.read_group_dynamic_pools[spec_index] = StringPool()
        return self.read_group_dynamic_pools[spec_index]

    def read_group_to_ids(self, group_strings: List[str]) -> List[int]:
        """
        Convert list of group strings to list of integer IDs.

        Args:
            group_strings: List of group ID strings

        Returns:
            List of integer IDs
        """
        if not group_strings:
            return []

        ids = []
        for spec_index, group_str in enumerate(group_strings):
            if group_str is not None:
                pool = self.get_read_group_pool(spec_index)
                ids.append(pool.get_int(group_str))
            else:
                ids.append(-1)  # Use -1 for None values
        return ids

    def read_group_from_ids(self, group_ids: List[int]) -> List[str]:
        """
        Convert list of integer IDs to list of group strings.

        Args:
            group_ids: List of integer IDs

        Returns:
            List of group ID strings
        """
        if not group_ids:
            return []

        strings = []
        for spec_index, group_id in enumerate(group_ids):
            if group_id < 0:
                strings.append(None)  # -1 means None
            else:
                pool = self.get_read_group_pool(spec_index)
                try:
                    strings.append(pool.get_str(group_id))
                except IndexError:
                    # Pool not populated yet or ID out of range
                    # This can happen during deserialization if pools weren't loaded
                    # (e.g., default grouper with 'NA' or dynamic pools not yet built)
                    # Fall back to returning the ID as a string
                    pool_type = self.group_spec_pool_types.get(spec_index, 'unknown')
                    logger.warning(f"Read group ID {group_id} not found in pool for spec {spec_index} (type={pool_type}, pool_size={len(pool)}), using ID as string")
                    strings.append(str(group_id))
        return strings

    def has_dynamic_pools(self) -> bool:
        """Check if any dynamic pools have data."""
        return any(len(pool) > 0 for pool in self.read_group_dynamic_pools.values())

    def serialize_dynamic_pools(self, outfile):
        """
        Serialize dynamic read group pools to binary file.

        Format:
            num_specs (int)
            for each spec:
                spec_index (int)
                num_strings (int)
                for each string:
                    string (length-prefixed)
        """

        # Only serialize non-empty dynamic pools
        non_empty_pools = {idx: pool for idx, pool in self.read_group_dynamic_pools.items()
                          if len(pool) > 0}

        write_int(len(non_empty_pools), outfile)
        for spec_index in sorted(non_empty_pools.keys()):
            pool = non_empty_pools[spec_index]
            write_int(spec_index, outfile)
            write_int(len(pool), outfile)
            # Write strings in ID order (0, 1, 2, ...) to preserve mapping
            for i in range(len(pool)):
                write_string(pool.get_str(i), outfile)

    def deserialize_dynamic_pools(self, infile):
        """
        Deserialize dynamic read group pools from binary file.

        Rebuilds pools with same ID mappings as when serialized.
        """
        num_specs = read_int(infile)
        for _ in range(num_specs):
            spec_index = read_int(infile)
            num_strings = read_int(infile)
            pool = StringPool()
            for _ in range(num_strings):
                string = read_string(infile)
                pool.add(string)  # IDs assigned in order: 0, 1, 2, ...
            self.read_group_dynamic_pools[spec_index] = pool

    def get_stats(self) -> str:
        """
        Get statistics about all pools.

        Returns:
            Formatted string with pool sizes
        """
        lines = ["String Pool Statistics:"]
        lines.append(f"  Genes: {len(self.gene_pool)}")
        lines.append(f"  Transcripts: {len(self.transcript_pool)}")
        lines.append(f"  Chromosomes: {len(self.chromosome_pool)}")
        lines.append(f"  Barcodes: {len(self.barcode_pool)}")
        lines.append(f"  UMIs: {len(self.umi_pool)}")
        lines.append(f"  File names: {len(self.file_name_pool)}")
        total_spots = sum(len(p) for p in self.barcode_spot_pools.values())
        lines.append(f"  Barcode spots: {total_spots} (across {len(self.barcode_spot_pools)} pools)")
        for spec_idx, pool in self.read_group_dynamic_pools.items():
            lines.append(f"  Read group dynamic (spec {spec_idx}): {len(pool)}")
        for spec_idx, pool in self.read_group_tsv_pools.items():
            lines.append(f"  Read group TSV (spec {spec_idx}): {len(pool)}")
        return "\n".join(lines)


def setup_string_pools(args, sample, chr_ids, chr_id=None, gffutils_db=None,
                       load_barcode_pool=False, load_tsv_pools=False, read_group_file_prefix=None):
    """
    Set up string pools for memory optimization during parallel processing.

    Args:
        args: Command-line arguments
        sample: SampleData object
        chr_ids: List of all chromosome IDs to process (for chromosome pool)
        chr_id: Current chromosome ID (required if loading per-chromosome pools)
        gffutils_db: Pre-loaded gffutils database (if None, will load from args.genedb)
        load_barcode_pool: Whether to load per-chromosome barcode pool
        load_tsv_pools: Whether to load per-chromosome TSV read group pools
        read_group_file_prefix: Override for read_group_file path (used when sample has modified prefix)

    Returns:
        StringPoolManager with pools configured
    """

    string_pools = StringPoolManager()

    # Build chromosome pool from the definitive list of chromosomes to process
    string_pools.build_chromosome_pool(chr_ids)

    # Build gene/transcript pools from annotation (if available)
    if gffutils_db is None and args.genedb:
        gffutils_db = load_genedb(args.genedb)
    if gffutils_db:
        string_pools.build_from_gffutils(gffutils_db)

    # Set up read group pool type mapping
    pool_types = get_grouping_pool_types(args)
    for spec_idx, pool_type in pool_types.items():
        string_pools.set_group_spec_pool_type(spec_idx, pool_type)

    # Build file_name pool if needed
    if any(pt == 'file_name' for pt in pool_types.values()):
        string_pools.build_file_name_pool(sample)

    # Build barcode_spot pools if needed (supports multi-column)
    if any(pt == 'barcode_spot' or pt.startswith('barcode_spot:') for pt in pool_types.values()):
        if args.barcode2spot:
            from .read_groups import parse_barcode2spot_spec
            filename, barcode_col, spot_cols = parse_barcode2spot_spec(args.barcode2spot)
            for col_idx, file_col in enumerate(spot_cols):
                string_pools.build_barcode_spot_pool(args.barcode2spot, column_index=col_idx,
                                                      file_column=file_col)

    # Load per-chromosome pools if requested
    if chr_id:
        if load_barcode_pool and sample.barcodes_split_reads:
            barcode_file = sample.get_barcodes_split_file(chr_id)
            string_pools.load_barcode_pool(barcode_file)

        # Use override if provided, otherwise use sample's read_group_file
        read_group_file = read_group_file_prefix if read_group_file_prefix else sample.read_group_file
        if load_tsv_pools and read_group_file:
            # Build mapping from tsv spec_index to (col_index, delimiter) from pool_types
            # pool_type format: 'tsv:SPEC_INDEX:COL_INDEX:DELIMITER'
            tsv_pool_info = {}  # spec_index -> list of (col_index, delimiter, pool_key)
            for pool_type in pool_types.values():
                if pool_type.startswith('tsv:'):
                    parts = pool_type.split(':')
                    if len(parts) >= 4:
                        tsv_spec_idx = int(parts[1])
                        col_idx = int(parts[2])
                        delimiter = parts[3]
                        # pool_key is spec_idx:col_idx to uniquely identify multi-column pools
                        pool_key = f"{tsv_spec_idx}:{col_idx}"
                        if tsv_spec_idx not in tsv_pool_info:
                            tsv_pool_info[tsv_spec_idx] = []
                        tsv_pool_info[tsv_spec_idx].append((col_idx, delimiter, pool_key))

            base_pattern = read_group_file + "_spec*_" + chr_id
            spec_files = sorted(glob.glob(base_pattern))
            for spec_file in spec_files:
                match = re.search(r'_spec(\d+)_', spec_file)
                if match:
                    spec_index = int(match.group(1))
                    # Load pool for each column in this spec file
                    for col_idx, delimiter, pool_key in tsv_pool_info.get(spec_index, [(1, '\t', f"{spec_index}:1")]):
                        string_pools.load_read_group_tsv_pool(spec_file, pool_key, col_idx, delimiter)

    return string_pools
