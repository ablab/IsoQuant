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

import logging
from typing import Optional, List, Dict

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
        self.barcode_spot_pool = StringPool()  # For barcode_spot grouping (static)

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

        Uses same logic as FileNameGrouper: basename without extension.
        Pool is built in sorted order for deterministic IDs across workers.

        Args:
            sample: Sample object with file_list attribute
        """
        import os
        # Collect all file names first
        file_names = set()
        # file_list is a list of lists (libraries), each library has one or more files
        for lib in sample.file_list:
            # Get basename without extension, matching FileNameGrouper logic
            readable_name = os.path.splitext(os.path.basename(lib[0]))[0]
            file_names.add(readable_name)
        # Add in sorted order for deterministic IDs
        for name in sorted(file_names):
            self.file_name_pool.add(name)
        logger.debug(f"File name pool: {len(self.file_name_pool)} unique readable names")

    def build_barcode_spot_pool(self, barcode2spot_files: List[str]):
        """
        Build barcode-to-spot pool from mapping files.

        Pool is built in sorted order for deterministic IDs across workers.

        Args:
            barcode2spot_files: List of TSV files mapping barcode -> spot/cell_type
        """
        # Collect all spots first
        spots = set()
        for bc2spot_file in barcode2spot_files:
            with open(bc2spot_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        spots.add(parts[1])
        # Add in sorted order for deterministic IDs
        for spot in sorted(spots):
            self.barcode_spot_pool.add(spot)
        logger.debug(f"Barcode spot pool: {len(self.barcode_spot_pool)} unique spots")

    def load_barcode_pool(self, barcode_file: str):
        """
        Load barcode pool from split barcode file.

        Pool is built in sorted order for deterministic IDs across workers.

        Args:
            barcode_file: Path to split barcode file for chromosome
                         Format: read_id\tbarcode\tumi
        """
        import os
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
        import os
        if not os.path.exists(read_group_file):
            logger.debug(f"Read group file not found: {read_group_file}")
            return

        # Collect all group values first
        # Always include 'NA' (default group ID) in case some reads are not in the TSV file
        from .read_groups import AbstractReadGrouper
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
        elif pool_type == 'barcode_spot':
            return self.barcode_spot_pool
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
        from .serialization import write_int, write_string

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
        from .serialization import read_int, read_string

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
        lines.append(f"  Barcode spots: {len(self.barcode_spot_pool)}")
        for spec_idx, pool in self.read_group_dynamic_pools.items():
            lines.append(f"  Read group dynamic (spec {spec_idx}): {len(pool)}")
        for spec_idx, pool in self.read_group_tsv_pools.items():
            lines.append(f"  Read group TSV (spec {spec_idx}): {len(pool)}")
        return "\n".join(lines)
