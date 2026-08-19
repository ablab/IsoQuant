
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import glob
import gzip
import logging
import os
import re
import shutil
from collections import defaultdict

from isoquant_lib.common import rreplace

GZIP_SUFFIX = ".gz"


def open_text_write(file_name):
    """Open for text writing, compressing when the name says so."""
    if file_name.endswith(GZIP_SUFFIX):
        return gzip.open(file_name, "wt")
    return open(file_name, "w")


def open_text_read(file_name):
    """Open for text reading, decompressing when the name says so."""
    if file_name.endswith(GZIP_SUFFIX):
        return gzip.open(file_name, "rt")
    return open(file_name, "r")


def resolve_optionally_gzipped(file_name):
    """Return the existing path among <file_name> and <file_name>.gz.

    Outputs that are compressed once the run finishes are still referred to by their plain
    name (in resumed runs, for instance), so readers have to accept either.
    """
    if os.path.exists(file_name):
        return file_name
    gzipped = file_name + GZIP_SUFFIX
    if os.path.exists(gzipped):
        return gzipped
    return file_name


def gzip_file_in_place(file_name, keep_original=False):
    """Compress a finished output to <file_name>.gz. Returns the resulting path."""
    if file_name.endswith(GZIP_SUFFIX):
        return file_name
    if not os.path.exists(file_name):
        return file_name
    gzipped = file_name + GZIP_SUFFIX
    with open(file_name, "rb") as inf, gzip.open(gzipped, "wb") as outf:
        shutil.copyfileobj(inf, outf)
    if not keep_original:
        os.remove(file_name)
    return gzipped


def strip_compression_suffix(file_name):
    """Drop a trailing compression suffix so extension logic sees the real one."""
    for suffix in (GZIP_SUFFIX, ".gzip", ".bgz"):
        if file_name.endswith(suffix):
            return file_name[:-len(suffix)]
    return file_name

logger = logging.getLogger('IsoQuant')


def merge_file_list(fname, label, chr_ids):
    # Per-chromosome files are written as "<label>_<chr_id><rest>" (see
    # SampleData.get_chr_prefix). The label prefixes the *basename*, so insert
    # "_<chr_id>" right after it there. Using rreplace on the whole path would
    # match the label anywhere (e.g. a short prefix like "p" inside
    # "...transcript...") and reconstruct the wrong per-chr name, silently
    # dropping counts and crashing on the missing files.
    directory, base = os.path.split(fname)
    if base.startswith(label):
        return [os.path.join(directory, f"{label}_{chr_id}{base[len(label):]}") for chr_id in chr_ids]
    # Defensive fallback for unexpected callers where the label is not a
    # basename prefix.
    return [rreplace(fname, label, f"{label}_{chr_id}") for chr_id in chr_ids]


def merge_files(file_name, label, chr_ids, merged_file_handler, copy_header=True, header_lines=None):
    file_names = merge_file_list(file_name, label, chr_ids)
    file_names.sort(key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)])
    for i, file_name in enumerate(file_names):
        if not os.path.exists(file_name): continue
        if header_lines is not None:
            header_count = header_lines
        else:
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
        if os.path.exists(file_name):
            os.remove(file_name)


def merge_counts(counter, label, chr_ids, unaligned_reads=0):
    file_name = counter.output_counts_file_name
    merged_file_handler = counter.get_output_file_handler()
    merge_files(file_name, label, chr_ids, merged_file_handler, header_lines=1)

    if counter.usable_file_name:
        for f in merge_file_list(counter.usable_file_name, label, chr_ids):
            counter.load_usable(f)
            os.remove(f)

    stat_dict = defaultdict(int)
    if counter.output_stats_file_name and counter.ignore_read_groups:
        stats_file_names = merge_file_list(counter.output_stats_file_name, label, chr_ids)
        for file_name in stats_file_names:
            for line in open(file_name):
                v = line.strip().split()
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
