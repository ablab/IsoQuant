############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import gzip
import os
import pysam
from collections import defaultdict


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
            return

        self.read_groups.add(values[-1])
        return values[-1]


class ReadTableGrouper(AbstractReadGrouper):
    def __init__(self, table_tsv_file, read_id_column_index=0, group_id_column_index=1, delim='\t'):
        AbstractReadGrouper.__init__(self)
        logger.debug("Reading read groups from " + table_tsv_file)
        self.read_map = load_table(table_tsv_file, read_id_column_index, group_id_column_index, delim)

    def get_group_id(self, alignment, filename=None):
        if alignment.query_name not in self.read_map:
            self.read_groups.add(self.default_group_id)
            return self.default_group_id
        self.read_groups.add(self.read_map[alignment.query_name])
        return self.read_map[alignment.query_name]


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


def get_file_grouping_properties(values):
    assert len(values) >= 2
    if len(values) > 4:
        return values[1], int(values[2]), int(values[3]), values[4]
    elif len(values) > 3:
        return values[1], int(values[2]), int(values[3]), "\t"
    else:
        return values[1], 0, 1, "\t"


def prepare_read_groups(args, sample):
    if not hasattr(args, "read_group") or args.read_group is None:
        return
    option = args.read_group
    values = option.split(':')
    if values[0] != 'file':
        return
    table_filename, read_id_column_index, group_id_column_index, delim = get_file_grouping_properties(values)
    logger.info("Splitting read group file %s for better memory consumption" % table_filename)
    split_read_group_table(table_filename, sample, read_id_column_index, group_id_column_index, delim)


def create_read_grouper(args, sample, chr_id):
    if not hasattr(args, "read_group") or args.read_group is None:
        return DefaultReadGrouper()

    option = args.read_group
    values = option.split(':')
    if values[0] == "file_name":
        return FileNameGrouper(args, sample)
    elif values[0] == 'tag':
        if len(values) < 2:
            return AlignmentTagReadGrouper(tag="RG")
        return AlignmentTagReadGrouper(tag=values[1])
    elif values[0] == 'read_id':
        return ReadIdSplitReadGrouper(delim=values[1])
    elif values[0] == 'file':
        read_group_chr_filename = sample.read_group_file + "_" + chr_id
        return ReadTableGrouper(read_group_chr_filename, 0, 1, '\t')
    else:
        logger.critical("Unsupported read grouping option")
        return DefaultReadGrouper()


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


def split_read_group_table(table_file, sample, read_id_column_index, group_id_column_index, delim):
    read_groups = load_table(table_file, read_id_column_index, group_id_column_index, delim)
    read_group_files = {}
    processed_reads = defaultdict(set)
    bam_files = list(map(lambda x: x[0], sample.file_list))

    for bam_file in bam_files:
        bam = pysam.AlignmentFile(bam_file, "rb")
        for chr_id in bam.references:
            read_group_files[chr_id] = open(sample.read_group_file + "_" + chr_id, "w")
        for read_alignment in bam:
            chr_id = read_alignment.reference_name
            if not chr_id:
                continue

            read_id = read_alignment.query_name
            if read_id in read_groups and read_id not in processed_reads[chr_id]:
                read_group_files[chr_id].write("%s\t%s\n" % (read_id, read_groups[read_id]))
                processed_reads[chr_id].add(read_id)

    for f in read_group_files.values():
        f.close()
