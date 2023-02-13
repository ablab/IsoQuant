############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import gzip
import os

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
        self.read_map = {}
        min_columns = max(read_id_column_index, group_id_column_index)
        logger.info("Reading read groups from " + table_tsv_file)
        _, outer_ext = os.path.splitext(table_tsv_file)
        if outer_ext.lower() in ['.gz', '.gzip']:
            handle = gzip.open(table_tsv_file, "rt")
        else:
            handle = open(table_tsv_file, 'r')

        for l in handle:
            l = l.strip()
            if l.startswith('#') or not l:
                continue

            column_values = l.split(delim)
            if len(column_values) <= min_columns:
                logger.warning("Malformed input read information table, minimum, of %d columns expected, "
                               "file %s, line: %s" % (min_columns, table_tsv_file, l))
                continue

            read_id = column_values[read_id_column_index]
            if read_id in self.read_map:
                logger.warning("Duplicate information for read %s" % read_id)

            group_id = column_values[group_id_column_index]
            self.read_map[read_id] = group_id

    def get_group_id(self, alignment, filename=None):
        if alignment.query_name not in self.read_map:
            self.read_groups.add(self.default_group_id)
            return self.default_group_id
        self.read_groups.add(self.read_map[alignment.query_name])
        return self.read_map[alignment.query_name]


class FileNameGroupper(AbstractReadGrouper):
    def __init__(self, args):
        AbstractReadGrouper.__init__(self)
        if args.input_data.readable_names_dict:
            self.readable_names_dict = args.input_data.readable_names_dict
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


def create_read_grouper(args):
    if not hasattr(args, "read_group") or args.read_group is None:
        return DefaultReadGrouper()

    option = args.read_group
    values = option.split(':')
    if values[0] == "file_name":
        return FileNameGroupper(args)
    elif values[0] == 'tag':
        return AlignmentTagReadGrouper(tag=values[1])
    elif values[0] == 'read_id':
        return ReadIdSplitReadGrouper(delim=values[1])
    elif values[0] == 'file':
        if len(values) > 4:
            return ReadTableGrouper(values[1], int(values[2]), int(values[3]), values[4])
        elif len(values) > 3:
            return ReadTableGrouper(values[1], int(values[2]), int(values[3]))
        else:
            return ReadTableGrouper(values[1])
    else:
        logger.critical("Unsupported read groupping option")
        return DefaultReadGrouper()
