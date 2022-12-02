import logging
import os
import re
import shutil
from collections import defaultdict

from src.read_groups import AbstractReadGrouper

logger = logging.getLogger('IsoQuant')


def merge_files(file_names, merged_file_name, stats_file_names=None, ignore_read_groups=True, copy_header=True):
    file_names.sort(key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)])
    with open(merged_file_name, 'ab') as outf:
        for i, file_name in enumerate(file_names):
            if not os.path.exists(file_name): continue
            header_count = 0
            with open(file_name, 'r') as f:
                while f.readline().startswith("#"):
                    header_count += 1
            with open(file_name, 'rb') as f:
                if not (copy_header and i == 0):
                    for j in range(header_count):
                        f.readline()
                shutil.copyfileobj(f, outf)
        for file_name in file_names:
            os.remove(file_name)

    ambiguous_reads, not_assigned_reads, not_aligned_reads = 0, 0, 0
    if stats_file_names and ignore_read_groups:
        for file_name in stats_file_names:
            with open(file_name) as f:
                line = f.readline()
                ambiguous_reads += int(line.split()[1])
                line = f.readline()
                not_assigned_reads += int(line.split()[1])
                line = f.readline()
                not_aligned_reads += int(line.split()[1])
            os.remove(file_name)
        with open(merged_file_name, 'a') as outf:
            outf.write("__ambiguous\t%d\n" % ambiguous_reads)
            outf.write("__no_feature\t%d\n" % not_assigned_reads)
            outf.write("__not_aligned\t%d\n" % not_aligned_reads)

