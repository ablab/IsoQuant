
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import os
import re
import shutil
from collections import defaultdict

from .common import  rreplace

logger = logging.getLogger('IsoQuant')


def merge_file_list(fname, label, chr_ids):
    return [rreplace(fname, label, f"{label}_{chr_id}") for chr_id in chr_ids]


def merge_files(file_name, label, chr_ids, merged_file_handler, copy_header=True):
    file_names = merge_file_list(file_name, label, chr_ids)
    file_names.sort(key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)])
    for i, file_name in enumerate(file_names):
        if not os.path.exists(file_name): continue
        header_count = 0
        with open(file_name, 'r') as f:
            while f.readline().startswith("#"):
                header_count += 1
        with open(file_name, 'rt') as f:
            if not (copy_header and i == 0):
                for j in range(header_count):
                    f.readline()
            shutil.copyfileobj(f, merged_file_handler)
    for file_name in file_names:
        os.remove(file_name)


def merge_counts(counter, label, chr_ids, unaligned_reads=0):
    file_name = counter.output_counts_file_name
    merged_file_handler = counter.get_output_file_handler()
    merge_files(file_name, label, chr_ids, merged_file_handler)

    if counter.usable_file_name:
        for f in merge_file_list(counter.usable_file_name, label, chr_ids):
            counter.load_usable(f)
            os.remove(f)

    stat_dict = defaultdict(int)
    if counter.output_stats_file_name and counter.ignore_read_groups:
        stats_file_names = merge_file_list(counter.output_stats_file_name, label, chr_ids)
        for file_name in stats_file_names:
            for l in open(file_name):
                v = l.strip().split()
                stat_dict[v[0]] += int(v[1])
            os.remove(file_name)

        if unaligned_reads > 0:
            stat_dict["__not_aligned"] = unaligned_reads
        for v in stat_dict.keys():
            merged_file_handler.write("%s\t%d\n" % (v, stat_dict[v]))


def normalize_path(config_path, file_path):
    if os.path.isabs(file_path):
        return os.path.normpath(file_path)
    else:
        return os.path.normpath(os.path.join(os.path.dirname(config_path), file_path))

