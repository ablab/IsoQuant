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
  - Global pools: gene, transcript, chromosome (from annotation)
  - Per-chromosome pools: barcode, UMI (from split files)
  - Read group pools: Multiple types depending on strategy

Usage:
- Workers build pools from gene database at start
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
    - Global pools: gene, transcript, chromosome (from annotation)
    - Per-chromosome pools: barcode, UMI (from split files)
    - Read group pools: Multiple types depending on grouping strategy
    """

    def __init__(self):
        """Initialize all string pools."""
        # Global pools (from annotation)
        self.gene_pool = StringPool()
        self.transcript_pool = StringPool()
        self.chromosome_pool = StringPool()

        # Per-chromosome pools (from split files)
        self.barcode_pool = StringPool()
        self.umi_pool = StringPool()

        # Read group pools by type
        self.read_group_tsv_pools: Dict[int, StringPool] = {}  # spec_index -> pool
        self.read_group_dynamic_pool = StringPool()  # For tag/read_id grouping
        self.file_name_pool = StringPool()  # For file_name grouping
        self.barcode_spot_pool = StringPool()  # For barcode_spot grouping

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

        # Chromosome pool
        chromosomes = set()
        for gene in gffutils_db.features_of_type('gene'):
            chromosomes.add(gene.seqid)
        for chr_id in sorted(chromosomes):
            self.chromosome_pool.add(chr_id)
        logger.debug(f"Chromosome pool: {len(self.chromosome_pool)} unique chromosomes")

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

        # Chromosome pool
        chromosomes = set()
        for gene in gene_db.genes.values():
            chromosomes.add(gene.chr_id)
        for chr_id in sorted(chromosomes):
            self.chromosome_pool.add(chr_id)
        logger.debug(f"Chromosome pool: {len(self.chromosome_pool)} unique chromosomes")

    def build_file_name_pool(self, sample):
        """
        Build file name pool from sample file list.

        Args:
            sample: Sample object with file_list attribute
        """
        import os
        for file_entry in sample.file_list:
            # file_entry can be either a tuple (file_path, type) or just file_path
            if isinstance(file_entry, (list, tuple)):
                file_path = file_entry[0]
            else:
                file_path = file_entry
            filename = os.path.basename(file_path)
            self.file_name_pool.add(filename)
        logger.debug(f"File name pool: {len(self.file_name_pool)} unique filenames")

    def build_barcode_spot_pool(self, barcode2spot_files: List[str]):
        """
        Build barcode-to-spot pool from mapping files.

        Args:
            barcode2spot_files: List of TSV files mapping barcode -> spot/cell_type
        """
        for bc2spot_file in barcode2spot_files:
            with open(bc2spot_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        spot = parts[1]
                        self.barcode_spot_pool.add(spot)
        logger.debug(f"Barcode spot pool: {len(self.barcode_spot_pool)} unique spots")

    def load_barcode_pool(self, barcode_file: str):
        """
        Load barcode pool from split barcode file.

        Args:
            barcode_file: Path to split barcode file for chromosome
                         Format: read_id\tbarcode\tumi
        """
        import os
        if not os.path.exists(barcode_file):
            logger.debug(f"Barcode file not found: {barcode_file}")
            return

        with open(barcode_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    barcode = parts[1]
                    self.barcode_pool.add(barcode)
                if len(parts) >= 3:
                    umi = parts[2]
                    self.umi_pool.add(umi)

        logger.debug(f"Barcode pool: {len(self.barcode_pool)} unique barcodes")
        logger.debug(f"UMI pool: {len(self.umi_pool)} unique UMIs")

    def load_read_group_tsv_pool(self, read_group_file: str, spec_index: int):
        """
        Load read group pool from split TSV file.

        Args:
            read_group_file: Path to split read_group file for chromosome
            spec_index: Index of the file spec (for multiple file specs)
        """
        import os
        if not os.path.exists(read_group_file):
            logger.debug(f"Read group file not found: {read_group_file}")
            return

        pool = StringPool()
        with open(read_group_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    group_value = parts[1]
                    pool.add(group_value)

        self.read_group_tsv_pools[spec_index] = pool
        logger.debug(f"Read group TSV pool (spec {spec_index}): {len(pool)} unique values")

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
            return self.read_group_dynamic_pool

        pool_type = self.group_spec_pool_types[spec_index]

        if pool_type == 'file_name':
            return self.file_name_pool
        elif pool_type == 'barcode_spot':
            return self.barcode_spot_pool
        elif pool_type.startswith('tsv:'):
            # Extract TSV spec index from 'tsv:N'
            tsv_spec_idx = int(pool_type.split(':')[1])
            if tsv_spec_idx in self.read_group_tsv_pools:
                return self.read_group_tsv_pools[tsv_spec_idx]
            else:
                # Pool not loaded yet, return empty dynamic pool
                return self.read_group_dynamic_pool
        else:  # 'dynamic' or unknown
            return self.read_group_dynamic_pool

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
                    logger.debug(f"Read group ID {group_id} not found in pool for spec {spec_index}, using ID as string")
                    strings.append(str(group_id))
        return strings

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
        lines.append(f"  Read group dynamic: {len(self.read_group_dynamic_pool)}")
        for spec_idx, pool in self.read_group_tsv_pools.items():
            lines.append(f"  Read group TSV (spec {spec_idx}): {len(pool)}")
        return "\n".join(lines)
