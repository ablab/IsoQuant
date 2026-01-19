############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import gzip
import os
import sys

import pysam
from collections import defaultdict

from .error_codes import IsoQuantExitCode
from .table_splitter import split_read_table_parallel


logger = logging.getLogger('IsoQuant')


class AbstractReadGrouper:
    default_group_id = 'NA'

    def __init__(self):
        self.read_groups = set()

    def get_group_id(self, alignment, read_assignment=None, filename=None):
        raise NotImplementedError()


class DefaultReadGrouper(AbstractReadGrouper):
    def __init__(self):
        AbstractReadGrouper.__init__(self)
        self.read_groups = {self.default_group_id}

    def get_group_id(self, alignment, read_assignment=None, filename=None):
        return self.default_group_id


class AlignmentTagReadGrouper(AbstractReadGrouper):
    def __init__(self, tag='RG'):
        AbstractReadGrouper.__init__(self)
        self.tag = tag

    def get_group_id(self, alignment, read_assignment=None, filename=None):
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

    def get_group_id(self, alignment, read_assignment=None, filename=None):
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

    def get_group_id(self, alignment, read_assignment=None, filename=None):
        if filename in self.readable_names_dict:
            self.read_groups.add(self.readable_names_dict[filename])
            return self.readable_names_dict[filename]
        if not filename:
            self.read_groups.add(AlignmentTagReadGrouper.default_group_id)
            return AlignmentTagReadGrouper.default_group_id
        self.read_groups.add(filename)
        return filename


class BarcodeSpotGrouper(AbstractReadGrouper):
    """Grouper that maps reads to spots/cell types via barcodes from read_assignment"""
    def __init__(self, barcode2spot_files):
        """
        Initialize barcode-to-spot grouper.

        Args:
            barcode2spot_files: List of TSV files mapping barcode -> spot/cell type
        """
        AbstractReadGrouper.__init__(self)

        # Load barcode2spot mapping only (barcode comes from read_assignment)
        self.barcode_to_spot = {}
        for barcode2spot_file in barcode2spot_files:
            logger.debug(f"Reading barcode-to-spot mapping from {barcode2spot_file}")
            self.barcode_to_spot.update(load_table(barcode2spot_file, 0, 1, '\t'))

        logger.debug(f"Loaded {len(self.barcode_to_spot)} barcode-spot mappings")

    def get_group_id(self, alignment, read_assignment=None, filename=None):
        """Map read to spot via barcode from read_assignment"""
        if read_assignment is None or read_assignment.barcode is None:
            self.read_groups.add(self.default_group_id)
            return self.default_group_id

        barcode = read_assignment.barcode

        # Look up spot for this barcode
        if barcode not in self.barcode_to_spot:
            self.read_groups.add(self.default_group_id)
            return self.default_group_id

        spot = self.barcode_to_spot[barcode]
        self.read_groups.add(spot)
        return spot


class BarcodeGrouper(AbstractReadGrouper):
    """Grouper that uses barcode directly from read_assignment"""
    def __init__(self):
        AbstractReadGrouper.__init__(self)

    def get_group_id(self, alignment, read_assignment=None, filename=None):
        """Return barcode as group ID"""
        if read_assignment is None or read_assignment.barcode is None:
            self.read_groups.add(self.default_group_id)
            return self.default_group_id

        barcode = read_assignment.barcode
        self.read_groups.add(barcode)
        return barcode


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

    def get_group_id(self, alignment, read_assignment=None, filename=None):
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

    def get_group_id(self, alignment, read_assignment=None, filename=None):
        """Returns a list of group IDs from all groupers"""
        group_ids = []
        for i, grouper in enumerate(self.groupers):
            gid = grouper.get_group_id(alignment, read_assignment, filename)
            group_ids.append(gid)
            self.read_groups[i].add(gid)
        return group_ids

    def get_all_groups(self):
        """Returns all unique groups across all groupers"""
        all_groups = []
        for groups in self.read_groups:
            all_groups.append(list(groups))
        return all_groups


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

    for spec_index, spec in enumerate(specs):
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
            logger.info("Splitting multi-column read group file %s (spec %d) for better memory consumption" %
                       (table_filename, spec_index))
        else:
            # Single column TSV
            group_id_column_index = int(values[3]) if len(values) > 3 else 1
            group_id_column_indices = [group_id_column_index]
            delim = values[4] if len(values) > 4 else '\t'
            logger.info("Splitting read group file %s (spec %d) for better memory consumption" %
                       (table_filename, spec_index))

        # Build output file names for each chromosome with spec_index to avoid overwrites
        split_reads_file_names = {chr_id: sample.read_group_file + "_spec" + str(spec_index) + "_" + chr_id
                                 for chr_id in chromosomes}

        # Use improved parallel splitting with line-by-line streaming
        num_threads = args.threads if hasattr(args, 'threads') else 4

        split_read_table_parallel(sample, table_filename, split_reads_file_names,
                                  num_threads,
                                  read_column=read_id_column_index,
                                  group_columns=tuple(group_id_column_indices),
                                  delim=delim)


def parse_grouping_spec(spec_string, args, sample, chr_id, spec_index=0):
    """
    Parse a single grouping specification and return the appropriate grouper(s).

    Args:
        spec_string: Grouping specification string
        args: Command line arguments
        sample: Sample object
        chr_id: Chromosome ID
        spec_index: Index of this spec in the list (for file-based grouping)

    Returns:
        - Single grouper for most cases
        - List of groupers for multi-column file specifications (to create separate grouped counts)
    """
    values = spec_string.split(':')

    if values[0] == "file_name":
        return FileNameGrouper(args, sample)
    elif values[0] == 'barcode_spot':
        # Always uses --barcode2spot files (no inline file specification)
        if not hasattr(args, 'barcode2spot') or not args.barcode2spot:
            logger.critical("barcode_spot grouping requires --barcode2spot")
            sys.exit(IsoQuantExitCode.MISSING_REQUIRED_OPTION)

        if not hasattr(sample, 'barcodes_split_reads') or not sample.barcodes_split_reads:
            logger.critical("barcode_spot grouping requires barcoded reads (use --barcoded_reads)")
            sys.exit(IsoQuantExitCode.MISSING_REQUIRED_OPTION)

        return BarcodeSpotGrouper(args.barcode2spot)
    elif values[0] == 'barcode':
        # Direct barcode grouping - uses barcode from read_assignment
        if not hasattr(sample, 'barcodes_split_reads') or not sample.barcodes_split_reads:
            logger.critical("barcode grouping requires barcoded reads (use --barcoded_reads)")
            sys.exit(IsoQuantExitCode.MISSING_REQUIRED_OPTION)

        return BarcodeGrouper()
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
            sys.exit(IsoQuantExitCode.INVALID_PARAMETER)

        read_id_column_index = int(values[2]) if len(values) > 2 else 0
        delim = values[4] if len(values) > 4 else '\t'
        group_col_spec = values[3] if len(values) > 3 else "1"

        if ',' in group_col_spec:
            # Multiple columns - create separate groupers sharing the same table data
            # This makes file:table.tsv:0:1,2,3 equivalent to three separate --read_group arguments
            group_id_column_indices = [int(x) for x in group_col_spec.split(',')]
            read_group_chr_filename = sample.read_group_file + "_spec" + str(spec_index) + "_" + chr_id

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
            group_id_column_index = int(group_col_spec)
            read_group_chr_filename = sample.read_group_file + "_spec" + str(spec_index) + "_" + chr_id

            # Create shared data with single column for consistency
            shared_data = SharedTableData(read_group_chr_filename, read_id_column_index,
                                         [group_id_column_index], delim)
            return ReadTableGrouper(shared_data, 0)
    else:
        logger.critical("Unsupported read grouping option: %s" % values[0])
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)


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
    for spec_index, spec in enumerate(specs):
        spec = spec.strip()
        if spec:
            grouper = parse_grouping_spec(spec, args, sample, chr_id, spec_index)
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


def get_grouping_pool_types(args) -> dict:
    """
    Build mapping from grouping spec index to pool type for string interning.

    Returns dict mapping spec_index -> pool_type, where pool_type is one of:
    - 'file_name': Uses file_name_pool (known in advance from file list)
    - 'barcode_spot': Uses barcode_spot_pool (known in advance from barcode2spot files)
    - 'tsv:SPEC_INDEX:COL_INDEX:DELIMITER': Uses read_group_tsv_pools[SPEC_INDEX:COL_INDEX]
      (loaded per-chromosome from split TSV files)
    - 'dynamic': Uses read_group_dynamic_pool (collected during processing)

    Args:
        args: Command-line arguments with read_group specification

    Returns:
        Dict[int, str]: Mapping from group spec index to pool type
    """
    if not hasattr(args, "read_group") or args.read_group is None:
        return {}

    # Handle both list (nargs='+') and string (backward compatibility)
    if isinstance(args.read_group, str):
        specs = args.read_group.split(';')
    else:
        specs = args.read_group

    pool_types = {}
    grouper_index = 0  # Track actual grouper index (expands for multi-column files)

    for spec_index, spec in enumerate(specs):
        spec = spec.strip()
        if not spec:
            continue

        values = spec.split(':')
        spec_type = values[0]

        if spec_type == "file_name":
            pool_types[grouper_index] = 'file_name'
            grouper_index += 1
        elif spec_type == 'barcode_spot':
            pool_types[grouper_index] = 'barcode_spot'
            grouper_index += 1
        elif spec_type in ['tag', 'read_id']:
            # BAM tags and read_id suffixes are discovered dynamically
            pool_types[grouper_index] = 'dynamic'
            grouper_index += 1
        elif spec_type == 'file':
            # TSV/CSV files: extract delimiter (default to tab)
            delimiter = values[4] if len(values) >= 5 else '\t'
            # Check if multi-column
            if len(values) >= 4 and ',' in values[3]:
                # Multi-column: each column is a separate grouper with its own column index
                # Format: tsv:SPEC_INDEX:COL_INDEX:DELIMITER
                group_col_indices = values[3].split(',')
                for col_idx in group_col_indices:
                    pool_types[grouper_index] = f'tsv:{spec_index}:{col_idx}:{delimiter}'
                    grouper_index += 1
            else:
                # Single column
                col_idx = values[3] if len(values) > 3 else '1'
                pool_types[grouper_index] = f'tsv:{spec_index}:{col_idx}:{delimiter}'
                grouper_index += 1

    return pool_types


def get_grouping_strategy_names(args) -> list:
    """
    Extract descriptive names for each grouping strategy from args.read_group.

    Returns a list of strategy names like: ["tag_CB", "tag_UB", "file_name", "read_id"]
    If no read_group is specified, returns ["default"].

    For multi-column TSV files with N columns, returns N separate names like:
    ["file0_col1", "file0_col2", "file0_col3"]

    For multiple file specs, includes spec index to avoid name collisions:
    ["file0_col1", "file1_col1", "file2_col1"]
    """
    if not hasattr(args, "read_group") or args.read_group is None:
        return ["default"]

    # Handle both list (nargs='+') and string (backward compatibility)
    if isinstance(args.read_group, str):
        specs = args.read_group.split(';')
    else:
        specs = args.read_group

    strategy_names = []
    for spec_index, spec in enumerate(specs):
        spec = spec.strip()
        if not spec:
            continue

        values = spec.split(':')
        spec_type = values[0]

        if spec_type == "file_name":
            strategy_names.append("file_name")
        elif spec_type == 'barcode_spot':
            strategy_names.append("barcode_spot")
        elif spec_type == 'barcode':
            strategy_names.append("barcode")
        elif spec_type == 'tag':
            tag_name = values[1] if len(values) > 1 else "RG"
            strategy_names.append(f"tag_{tag_name}")
        elif spec_type == 'read_id':
            delim = values[1] if len(values) > 1 else "_"
            # Sanitize delimiter for filename
            safe_delim = delim.replace('/', '_').replace('\\', '_')
            strategy_names.append(f"read_id_{safe_delim}")
        elif spec_type == 'file':
            # Include spec_index to distinguish between different file specs
            # Check if multi-column
            if len(values) >= 4 and ',' in values[3]:
                group_col_indices = values[3].split(',')
                for col_idx in group_col_indices:
                    strategy_names.append(f"file{spec_index}_col{col_idx}")
            else:
                col_idx = values[3] if len(values) > 3 else "1"
                strategy_names.append(f"file{spec_index}_col{col_idx}")

    return strategy_names if strategy_names else ["default"]


def load_table(table_tsv_file, read_id_column_index, group_id_column_index, delim):
    min_columns = max(read_id_column_index, group_id_column_index)
    _, outer_ext = os.path.splitext(table_tsv_file)
    if outer_ext.lower() in ['.gz', '.gzip']:
        handle = gzip.open(table_tsv_file, "rt")
    else:
        handle = open(table_tsv_file, 'r')

    read_map = {}
    warn_count = 0
    for line in handle:
        line = line.strip()
        if line.startswith('#') or not line:
            continue

        column_values = line.split(delim)
        if len(column_values) <= min_columns:
            if warn_count == 0:
                logger.warning("Malformed input read information table, minimum, of %d columns expected, "
                               "file %s, line: %s" % (min_columns, table_tsv_file, line))
            warn_count += 1
            continue

        read_id = column_values[read_id_column_index]
        if read_id in read_map:
            logger.warning("Duplicate information for read %s" % read_id)

        group_id = column_values[group_id_column_index]
        read_map[read_id] = group_id

    if warn_count > 0:
        logger.warning("Total number of malformed lines in %s: %d" % (table_tsv_file, warn_count))
    return read_map


def load_multicolumn_table(table_tsv_file, read_id_column_index, group_id_column_indices, delim):
    min_columns = max(read_id_column_index, max(group_id_column_indices))
    _, outer_ext = os.path.splitext(table_tsv_file)
    if outer_ext.lower() in ['.gz', '.gzip']:
        handle = gzip.open(table_tsv_file, "rt")
    else:
        handle = open(table_tsv_file, 'r')

    read_map = {}
    warn_count = 0
    for line in handle:
        line = line.strip()
        if line.startswith('#') or not line:
            continue

        column_values = line.split(delim)
        if len(column_values) <= min_columns:
            if warn_count == 0:
                logger.warning("Malformed input read information table, minimum, of %d columns expected, "
                               "file %s, line: %s" % (min_columns, table_tsv_file, line))
            warn_count += 1
            continue

        read_id = column_values[read_id_column_index]
        if read_id in read_map:
            logger.warning("Duplicate information for read %s" % read_id)

        column_data = [column_values[i] for i in group_id_column_indices]
        read_map[read_id] = column_data

    if warn_count > 0:
        logger.warning("Total number of malformed lines in %s: %d" % (table_tsv_file, warn_count))
    return read_map
