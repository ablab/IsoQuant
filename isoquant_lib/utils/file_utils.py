
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import gzip
import logging
import os
import re
import shutil
import sys
from collections import defaultdict
from typing import Dict, List, Optional

from isoquant_lib.common import rreplace
from isoquant_lib.utils.error_codes import IsoQuantExitCode
from isoquant_lib.utils.file_naming import convert_chr_id_to_file_name_str

logger = logging.getLogger('IsoQuant')

GZIP_SUFFIX = ".gz"

# Levels for the two kinds of output IsoQuant writes, measured on real ONT data (see below).
# Python's gzip defaults to 9, which is a bad trade for both.
#
# Tables (TSV, BED, MTX, allinfo): 6. On a real barcode table level 9 runs at 16 MB/s against
# 39 MB/s at 6 for 3% less output, and going below 6 saves little -- level 5 is 11% faster for
# 0.4% more, level 1 is 2.4x faster for 8.5% more.
GZIP_LEVEL = 6
# Sequences (FASTA/FASTQ): 4. Nucleotide data sits near gzip's entropy floor, so the high
# levels grind: 12.5 MB/s at 6 against 65.3 MB/s at 4, for 5.7% more output. These are also
# the largest files IsoQuant writes.
GZIP_LEVEL_SEQUENCES = 4

SEQUENCE_SUFFIXES = (".fa", ".fasta", ".fq", ".fastq")


def gzip_level_for(file_name):
    """Compression level for an output, chosen by the kind of data its name implies."""
    if strip_compression_suffix(file_name).endswith(SEQUENCE_SUFFIXES):
        return GZIP_LEVEL_SEQUENCES
    return GZIP_LEVEL


def open_text_write(file_name, compresslevel=None):
    """Open for text writing, compressing when the name says so.

    The level defaults to what the name implies; pass one to override.
    """
    if file_name.endswith(GZIP_SUFFIX):
        if compresslevel is None:
            compresslevel = gzip_level_for(file_name)
        return gzip.open(file_name, "wt", compresslevel=compresslevel)
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
    with open(file_name, "rb") as inf, \
            gzip.open(gzipped, "wb", compresslevel=gzip_level_for(file_name)) as outf:
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


def check_file_exists(file_path: str, description: str):
    """Check that a file exists, exit with error if not."""
    if not os.path.isfile(file_path):
        logger.critical(f"{description} {file_path} does not exist")
        sys.exit(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)


def read_stats_tsv(file_name: Optional[str], stats: Optional[Dict[str, int]] = None) -> Dict[str, int]:
    """Read one of the two-column "name<TAB>count" stat files the stages write.

    Counts are added to stats when one is given, so the per-chunk files of a stage can
    be summed by calling this once per file. A missing or empty path contributes
    nothing, and lines that are not a name and a number are skipped.
    """
    if stats is None:
        stats = {}
    if not file_name or not os.path.exists(file_name):
        return stats
    with open(file_name) as f:
        for line in f:
            values = line.rstrip("\n").split("\t")
            if len(values) != 2:
                continue
            try:
                count = int(values[1])
            except ValueError:
                continue
            stats[values[0]] = stats.get(values[0], 0) + count
    return stats



def merge_file_list(fname: str, label: str, chr_ids, fragment_dir: Optional[str] = None) -> List[str]:
    # Per-chromosome files are written as "<label>_<chr_id><rest>" (see
    # SampleData.get_chr_prefix). The label prefixes the *basename*, so insert
    # "_<chr_id>" right after it there. Using rreplace on the whole path would
    # match the label anywhere (e.g. a short prefix like "p" inside
    # "...transcript...") and reconstruct the wrong per-chr name, silently
    # dropping counts and crashing on the missing files.
    # The chromosome id is made file-name safe exactly like the writers do, and the
    # fragments live in fragment_dir (aux/per_chr) when it is given.
    directory, base = os.path.split(fname)
    if fragment_dir is not None:
        directory = fragment_dir
    if base.startswith(label):
        return [os.path.join(directory, f"{label}_{convert_chr_id_to_file_name_str(chr_id)}{base[len(label):]}")
                for chr_id in chr_ids]
    if fragment_dir is not None:
        raise ValueError("Cannot derive per-chromosome file names: %s does not start with %s" % (base, label))
    # Defensive fallback for unexpected callers where the label is not a
    # basename prefix.
    return [rreplace(fname, label, f"{label}_{convert_chr_id_to_file_name_str(chr_id)}") for chr_id in chr_ids]


def merge_files(file_name, label, chr_ids, merged_file_handler, copy_header=True, header_lines=None,
                remove_inputs: bool = True, fragment_dir: Optional[str] = None) -> List[str]:
    """Concatenate the per-chromosome fragments of file_name into merged_file_handler.

    Returns the fragment names; with remove_inputs=False they are kept, so the caller can
    delete them only after its checkpoint marker is written.
    """
    file_names = merge_file_list(file_name, label, chr_ids, fragment_dir)
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
    if remove_inputs:
        remove_files(file_names)
    return file_names


def merge_counts(counter, label, chr_ids, unaligned_reads=0, remove_inputs: bool = True,
                 fragment_dir: Optional[str] = None) -> List[str]:
    """Merge the per-chromosome counts and stats rows of counter into its output file.

    The .usable fragments are not touched here: load_usable_fragments reads them right
    before finalize, which is the only step that needs them. Returns the consumed
    fragment names (counts and stats).
    """
    file_name = counter.output_counts_file_name
    with counter.get_output_file_handler() as merged_file_handler:
        consumed = merge_files(file_name, label, chr_ids, merged_file_handler, header_lines=1,
                               remove_inputs=remove_inputs, fragment_dir=fragment_dir)

        stat_dict = defaultdict(int)
        if counter.output_stats_file_name and counter.ignore_read_groups:
            stats_file_names = merge_file_list(counter.output_stats_file_name, label, chr_ids, fragment_dir)
            for file_name in stats_file_names:
                for line in open(file_name):
                    v = line.strip().split()
                    stat_dict[v[0]] += int(v[1])
            consumed += stats_file_names
            if remove_inputs:
                remove_files(stats_file_names)

            if unaligned_reads > 0:
                stat_dict["__not_aligned"] = unaligned_reads
            for v in stat_dict.keys():
                merged_file_handler.write("%s\t%d\n" % (v, stat_dict[v]))
    return consumed


def usable_fragment_list(counter, label, chr_ids, fragment_dir: Optional[str] = None) -> List[str]:
    if not counter.usable_file_name:
        return []
    return merge_file_list(counter.usable_file_name, label, chr_ids, fragment_dir)


def load_usable_fragments(counter, label, chr_ids, fragment_dir: Optional[str] = None) -> List[str]:
    """Load the per-chromosome usable-read counts (TPM normalisation) into counter.

    Returns the fragment names, which the caller removes once finalize is done.
    """
    fragments = usable_fragment_list(counter, label, chr_ids, fragment_dir)
    for f in fragments:
        counter.load_usable(f)
    return fragments


def remove_files(file_names) -> None:
    for f in file_names:
        if os.path.exists(f):
            os.remove(f)


def normalize_path(config_path, file_path):
    if os.path.isabs(file_path):
        return os.path.normpath(file_path)
    else:
        return os.path.normpath(os.path.join(os.path.dirname(config_path), file_path))
