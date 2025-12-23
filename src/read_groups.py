############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import gzip
import os
import pysam
from collections import defaultdict

from .table_splitter import split_read_table_parallel


logger = logging.getLogger('IsoQuant')


class AbstractReadGrouper:
    default_group_id = 'NA'

    def __init__(self):
        self.read_groups = set()

    def get_group_id(self, alignment, filename=None):
        raise NotImplementedError()


class DefaultReadGrouper(AbstractReadGrouper):
    def __init__(self):
        AbstractReadGrouper.__init__(self)
        self.read_groups = {self.default_group_id}

    def get_group_id(self, alignment, filename=None):
        return self.default_group_id


class AlignmentTagReadGrouper(AbstractReadGrouper):
    def __init__(self, tag='RG'):
        AbstractReadGrouper.__init__(self)
        self.tag = tag

    def get_group_id(self, alignment, filename=None):
        try:
            tag_value = alignment.get_tag(self.tag)
        except KeyError:
            logger.warning("Tag %s is not present for read %s, skipping" % (self.tag, alignment.query_name))
            self.read_groups.add(self.default_group_id)
            return self.default_group_id
        self.read_groups.add(tag_value)
        return tag_value


class ReadIdSplitReadGrouper(AbstractReadGrouper):
    def __init__(self, delim):
        AbstractReadGrouper.__init__(self)
        self.delim = delim

    def get_group_id(self, alignment, filename=None):
        read_id = alignment.query_name
        values = read_id.split(self.delim)
        if len(values) == 1:
            logger.warning("Delimiter %s is not present in read id %s, skipping" % (self.delim, read_id))
            return ""

        self.read_groups.add(values[-1])
        return values[-1]


class FileNameGrouper(AbstractReadGrouper):
    def __init__(self, args, sample):
        AbstractReadGrouper.__init__(self)
        if sample.readable_names_dict:
            self.readable_names_dict = sample.readable_names_dict
            return

        self.readable_names_dict = {}
        for sample in args.input_data.samples:
            for lib in sample.file_list:
                readable_name = os.path.splitext(os.path.basename(lib[0]))[0]
                for f in lib:
                    self.readable_names_dict[f] = readable_name

    def get_group_id(self, alignment, filename=None):
        if filename in self.readable_names_dict:
            self.read_groups.add(self.readable_names_dict[filename])
            return self.readable_names_dict[filename]
        if not filename:
            self.read_groups.add(AlignmentTagReadGrouper.default_group_id)
            return AlignmentTagReadGrouper.default_group_id
        self.read_groups.add(filename)
        return filename


# TODO: remove barocde table, use assignment.barcode
class BarcodeSpotGrouper(AbstractReadGrouper):
    """Grouper that maps reads to spots/cell types via barcodes"""
    def __init__(self, barcode_file, barcode2spot_files):
        """
        Initialize barcode-to-spot grouper.

        Args:
            barcode_file: Path to split barcode file for chromosome (read_id -> barcode, umi)
            barcode2spot_files: List of TSV files mapping barcode -> spot/cell type
        """
        AbstractReadGrouper.__init__(self)
        logger.debug(f"Reading barcodes from {barcode_file}")

        # Load barcode dict: read_id -> (barcode, umi)
        self.read_to_barcode = {}
        if os.path.exists(barcode_file):
            for line in open(barcode_file):
                if line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    # Store just the barcode (second column)
                    self.read_to_barcode[parts[0]] = parts[1]

        # Load barcode2spot mapping: barcode -> spot/cell type
        self.barcode_to_spot = {}
        for barcode2spot_file in barcode2spot_files:
            logger.debug(f"Reading barcode-to-spot mapping from {barcode2spot_file}")
            self.barcode_to_spot.update(load_table(barcode2spot_file, 0, 1, '\t'))

        logger.info(f"Loaded {len(self.read_to_barcode)} read-barcode mappings and "
                   f"{len(self.barcode_to_spot)} barcode-spot mappings")

    def get_group_id(self, alignment, filename=None):
        """Map read to spot via barcode"""
        read_id = alignment.query_name

        # Look up barcode for this read
        if read_id not in self.read_to_barcode:
            self.read_groups.add(self.default_group_id)
            return self.default_group_id

        barcode = self.read_to_barcode[read_id]

        # Look up spot for this barcode
        if barcode not in self.barcode_to_spot:
            self.read_groups.add(self.default_group_id)
            return self.default_group_id

        spot = self.barcode_to_spot[barcode]
        self.read_groups.add(spot)
        return spot


class SharedTableData:
    """
    Shared storage for table data to avoid loading the same file multiple times.

    Allows multiple ReadTableGrouper instances to reference the same
    loaded table data, reducing memory usage and I/O when using multiple columns
    from the same file.
    """
    def __init__(self, table_tsv_file, read_id_column_index=0, group_id_column_indices=None, delim='\t'):
        """
        Load table data once for sharing across multiple groupers.

        Args:
            table_tsv_file: Path to TSV file
            read_id_column_index: Column index for read IDs
            group_id_column_indices: List of column indices to load
            delim: Column delimiter
        """
        if group_id_column_indices is None:
            group_id_column_indices = [1]
        logger.debug("Loading shared table data from " + table_tsv_file)
        # Load all specified columns into a dict: read_id -> [col1_val, col2_val, ...]
        self.read_map = load_multicolumn_table(table_tsv_file, read_id_column_index,
                                               group_id_column_indices, delim)
        self.num_columns = len(group_id_column_indices)


class ReadTableGrouper(AbstractReadGrouper):
    """
    Grouper that uses a single column from shared table data.

    Multiple instances can reference the same SharedTableData to avoid
    redundant file loading and memory usage.
    """
    def __init__(self, shared_data, column_index):
        """
        Initialize grouper for a specific column.

        Args:
            shared_data: SharedTableData instance
            column_index: Index within the shared data's column list (0-based)
        """
        AbstractReadGrouper.__init__(self)
        self.shared_data = shared_data
        self.column_index = column_index

    def get_group_id(self, alignment, filename=None):
        """Returns the group ID for this column"""
        if alignment.query_name not in self.shared_data.read_map:
            self.read_groups.add(self.default_group_id)
            return self.default_group_id

        group_ids = self.shared_data.read_map[alignment.query_name]
        group_id = group_ids[self.column_index]
        self.read_groups.add(group_id)
        return group_id


class MultiReadGrouper:
    """Manages multiple read groupers and returns list of group IDs"""
    def __init__(self, groupers):
        """
        Initialize with a list of groupers
        Args:
            groupers: list of AbstractReadGrouper objects
        """
        if not isinstance(groupers, list):
            groupers = [groupers]
        self.groupers = groupers
        self.read_groups = [set() for _ in groupers]

    def get_group_id(self, alignment, filename=None):
        """Returns a list of group IDs from all groupers"""
        group_ids = []
        for i, grouper in enumerate(self.groupers):
            gid = grouper.get_group_id(alignment, filename)
            group_ids.append(gid)
            self.read_groups[i].add(gid)
        return group_ids

    def get_all_groups(self):
        """Returns all unique groups across all groupers"""
        all_groups = []
        for groups in self.read_groups:
            all_groups.append(list(groups))
        return all_groups


def get_file_grouping_properties(values):
    assert len(values) >= 2
    if len(values) > 4:
        return values[1], int(values[2]), int(values[3]), values[4]
    elif len(values) > 3:
        return values[1], int(values[2]), int(values[3]), "\t"
    else:
        return values[1], 0, 1, "\t"


def prepare_read_groups(args, sample):
    """
    Prepare read group files by splitting them by chromosome for memory efficiency.
    Handles both single and multiple grouping specifications.
    Uses improved parallel algorithm for better performance.

    args.read_group should be a list of grouping specifications (nargs='+').
    """
    if not hasattr(args, "read_group") or args.read_group is None:
        return

    # Handle both list (nargs='+') and string (backward compatibility)
    if isinstance(args.read_group, str):
        specs = args.read_group.split(';')
    else:
        specs = args.read_group

    # Collect chromosome names from BAM files to create split files
    bam_files = list(map(lambda x: x[0], sample.file_list))
    chromosomes = set()
    for bam_file in bam_files:
        try:
            bam = pysam.AlignmentFile(bam_file, "rb")
            chromosomes.update(bam.references)
            bam.close()
        except:
            pass
    chromosomes = list(chromosomes)

    for spec in specs:
        spec = spec.strip()
        if not spec:
            continue

        values = spec.split(':')
        if values[0] != 'file':
            continue

        # Parse specification
        table_filename = values[1]
        read_id_column_index = int(values[2]) if len(values) > 2 else 0

        if len(values) >= 4 and ',' in values[3]:
            # Multi-column TSV: file:filename:read_col:group_cols:delim
            group_id_column_indices = [int(x) for x in values[3].split(',')]
            delim = values[4] if len(values) > 4 else '\t'
            logger.info("Splitting multi-column read group file %s for better memory consumption" % table_filename)
        else:
            # Single column TSV
            group_id_column_index = int(values[3]) if len(values) > 3 else 1
            group_id_column_indices = [group_id_column_index]
            delim = values[4] if len(values) > 4 else '\t'
            logger.info("Splitting read group file %s for better memory consumption" % table_filename)

        # Build output file names for each chromosome
        split_reads_file_names = {chr_id: sample.read_group_file + "_" + chr_id for chr_id in chromosomes}

        # Use improved parallel splitting with line-by-line streaming
        num_threads = args.threads if hasattr(args, 'threads') else 4

        split_read_table_parallel(sample, table_filename, split_reads_file_names,
                                  num_threads,
                                  read_column=read_id_column_index,
                                  group_columns=tuple(group_id_column_indices),
                                  delim=delim)


def parse_grouping_spec(spec_string, args, sample, chr_id):
    """
    Parse a single grouping specification and return the appropriate grouper(s).

    Returns:
        - Single grouper for most cases
        - List of groupers for multi-column file specifications (to create separate grouped counts)
    """
    values = spec_string.split(':')

    if values[0] == "file_name":
        return FileNameGrouper(args, sample)
    elif values[0] == 'barcode_spot':
        # Format: barcode_spot:file1.tsv or just use --barcode2spot files
        if len(values) >= 2:
            # Explicit file(s) specified
            barcode2spot_files = values[1:]
        elif hasattr(args, 'barcode2spot') and args.barcode2spot:
            # Use --barcode2spot files
            barcode2spot_files = args.barcode2spot
        else:
            logger.critical("barcode_spot grouping requires --barcode2spot or explicit file specification")
            return None

        # Get split barcode file for this chromosome
        if not hasattr(sample, 'barcodes_split_reads') or not sample.barcodes_split_reads:
            logger.critical("barcode_spot grouping requires barcoded reads (use --barcoded_reads)")
            return None
        # FIXME
        barcode_file = sample.barcodes_split_reads + "_" + chr_id
        return BarcodeSpotGrouper(barcode_file, barcode2spot_files)
    elif values[0] == 'tag':
        if len(values) < 2:
            return AlignmentTagReadGrouper(tag="RG")
        return AlignmentTagReadGrouper(tag=values[1])
    elif values[0] == 'read_id':
        return ReadIdSplitReadGrouper(delim=values[1])
    elif values[0] == 'file':
        # Format: file:filename:read_col:group_cols:delim
        # group_cols can be comma-separated like "1,2,3"
        if len(values) < 2:
            logger.critical("group specification %s is too short, specifiy at least a file name" % spec_string)
            return None

        read_id_column_index = int(values[2]) if len(values) > 2 else 0
        delim = values[4] if len(values) > 4 else '\t'
        group_col_spec = values[3] if len(values) > 3 else "1"

        if ',' in group_col_spec:
            # Multiple columns - create separate groupers sharing the same table data
            # This makes file:table.tsv:0:1,2,3 equivalent to three separate --read_group arguments
            group_id_column_indices = [int(x) for x in group_col_spec.split(',')]
            read_group_chr_filename = sample.read_group_file + "_" + chr_id

            # Create shared table data once
            shared_data = SharedTableData(read_group_chr_filename, read_id_column_index,
                                         group_id_column_indices, delim)

            # Create a separate grouper for each column
            groupers = []
            for i in range(len(group_id_column_indices)):
                groupers.append(ReadTableGrouper(shared_data, i))

            return groupers  # Return list of groupers
        else:
            # Single column - use ReadTableGrouper
            group_id_column_index = int(values[3])
            read_group_chr_filename = sample.read_group_file + "_" + chr_id

            # Create shared data with single column for consistency
            shared_data = SharedTableData(read_group_chr_filename, read_id_column_index,
                                         [group_id_column_index], delim)
            return ReadTableGrouper(shared_data, 0)
    else:
        logger.critical("Unsupported read grouping option: %s" % values[0])
        return None


def create_read_grouper(args, sample, chr_id):
    """
    Create read grouper(s) based on args.read_group specification.

    args.read_group should be a list of grouping specifications (nargs='+'):
    - Single grouper: ["tag:RG"] or ["file_name"]
    - Multiple groupers: ["tag:RG", "file_name", "read_id:_"]
    - Multi-column TSV: ["file:table.tsv:0:1,2,3"] (columns 1,2,3 as separate groups)

    Returns:
        MultiReadGrouper if multiple specifications, otherwise single grouper
    """
    if not hasattr(args, "read_group") or args.read_group is None:
        return DefaultReadGrouper()

    specs = args.read_group
    groupers = []
    for spec in specs:
        spec = spec.strip()
        if spec:
            grouper = parse_grouping_spec(spec, args, sample, chr_id)
            if grouper:
                # parse_grouping_spec can return either a single grouper or a list of groupers
                # (multi-column file specs return lists)
                if isinstance(grouper, list):
                    groupers.extend(grouper)
                else:
                    groupers.append(grouper)

    if not groupers:
        logger.warning("No valid groupers specified, using default")
        return DefaultReadGrouper()
    elif len(groupers) == 1:
        return groupers[0]
    else:
        return MultiReadGrouper(groupers)


def get_grouping_strategy_names(args) -> list:
    """
    Extract descriptive names for each grouping strategy from args.read_group.

    Returns a list of strategy names like: ["tag_CB", "tag_UB", "file_name", "read_id"]
    If no read_group is specified, returns ["default"].

    For multi-column TSV files with N columns, returns N separate names like:
    ["file_col1", "file_col2", "file_col3"]
    """
    if not hasattr(args, "read_group") or args.read_group is None:
        return ["default"]

    # Handle both list (nargs='+') and string (backward compatibility)
    if isinstance(args.read_group, str):
        specs = args.read_group.split(';')
    else:
        specs = args.read_group

    strategy_names = []
    for spec in specs:
        spec = spec.strip()
        if not spec:
            continue

        values = spec.split(':')
        spec_type = values[0]

        if spec_type == "file_name":
            strategy_names.append("file_name")
        elif spec_type == 'barcode_spot':
            strategy_names.append("barcode_spot")
        elif spec_type == 'tag':
            tag_name = values[1] if len(values) > 1 else "RG"
            strategy_names.append(f"tag_{tag_name}")
        elif spec_type == 'read_id':
            delim = values[1] if len(values) > 1 else "_"
            # Sanitize delimiter for filename
            safe_delim = delim.replace('/', '_').replace('\\', '_')
            strategy_names.append(f"read_id_{safe_delim}")
        elif spec_type == 'file':
            # Check if multi-column
            if len(values) >= 4 and ',' in values[3]:
                group_col_indices = values[3].split(',')
                for col_idx in group_col_indices:
                    strategy_names.append(f"file_col{col_idx}")
            else:
                col_idx = values[3] if len(values) > 3 else "1"
                strategy_names.append(f"file_col{col_idx}")

    return strategy_names if strategy_names else ["default"]


def load_table(table_tsv_file, read_id_column_index, group_id_column_index, delim):
    min_columns = max(read_id_column_index, group_id_column_index)
    _, outer_ext = os.path.splitext(table_tsv_file)
    if outer_ext.lower() in ['.gz', '.gzip']:
        handle = gzip.open(table_tsv_file, "rt")
    else:
        handle = open(table_tsv_file, 'r')

    read_map = {}
    for line in handle:
        line = line.strip()
        if line.startswith('#') or not line:
            continue

        column_values = line.split(delim)
        if len(column_values) <= min_columns:
            logger.warning("Malformed input read information table, minimum, of %d columns expected, "
                           "file %s, line: %s" % (min_columns, table_tsv_file, line))
            continue

        read_id = column_values[read_id_column_index]
        if read_id in read_map:
            logger.warning("Duplicate information for read %s" % read_id)

        group_id = column_values[group_id_column_index]
        read_map[read_id] = group_id
    return read_map


def load_multicolumn_table(table_tsv_file, read_id_column_index, group_id_column_indices, delim):
    min_columns = max(read_id_column_index, max(group_id_column_indices))
    _, outer_ext = os.path.splitext(table_tsv_file)
    if outer_ext.lower() in ['.gz', '.gzip']:
        handle = gzip.open(table_tsv_file, "rt")
    else:
        handle = open(table_tsv_file, 'r')

    read_map = {}
    for line in handle:
        line = line.strip()
        if line.startswith('#') or not line:
            continue

        column_values = line.split(delim)
        if len(column_values) <= min_columns:
            logger.warning("Malformed input read information table, minimum, of %d columns expected, "
                           "file %s, line: %s" % (min_columns, table_tsv_file, line))
            continue

        read_id = column_values[read_id_column_index]
        if read_id in read_map:
            logger.warning("Duplicate information for read %s" % read_id)

        column_data = [column_values[i] for i in group_id_column_indices]
        read_map[read_id] = column_data
    return read_map


def split_table(table_file, sample, out_prefix, read_id_column_index, group_id_column_indices, delim):
    read_groups = load_multicolumn_table(table_file, read_id_column_index, group_id_column_indices, delim)
    read_group_files = {}
    processed_reads = defaultdict(set)
    bam_files = list(map(lambda x: x[0], sample.file_list))

    for bam_file in bam_files:
        bam = pysam.AlignmentFile(bam_file, "rb")
        for chr_id in bam.references:
            if chr_id not in read_group_files:
                read_group_files[chr_id] = open(out_prefix + "_" + chr_id, "w")
        for read_alignment in bam:
            chr_id = read_alignment.reference_name
            if not chr_id:
                continue

            read_id = read_alignment.query_name
            if read_id in read_groups and read_id not in processed_reads[chr_id]:
                # read_groups[read_id] is a list, join with delimiter
                group_values = delim.join(read_groups[read_id])
                read_group_files[chr_id].write("%s\t%s\n" % (read_id, group_values))
                processed_reads[chr_id].add(read_id)

    for f in read_group_files.values():
        f.close()


def prepare_barcoded_reads(args, sample):
    logger.info("Splitting barcoded reads %s for better memory consumption" % sample.out_barcodes_tsv)
    split_table(sample.out_barcodes_tsv, sample, sample.barcodes_split_reads, 0, [1, 2, 3, 4], '\t')
