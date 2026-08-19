############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Writing BAM outputs: merging per-chromosome fragments, and copying reads with tags.

The pipeline processes chromosomes in parallel, so any BAM it emits is built as one fragment
per chromosome and merged at the end. Records are copied verbatim apart from the tags being
added, because downstream consumers depend on fields IsoQuant does not itself use -- fusion
detection in particular reads breakpoints from the SA tag and from terminal soft clips.
"""

import logging
import os
import shutil
from typing import Dict, Iterable, Optional, Tuple

import pysam

logger = logging.getLogger('IsoQuant')

# Tags for the gene and transcript a read was assigned to, following the 10x convention.
GENE_TAG = "GX"
TRANSCRIPT_TAG = "TX"

# read_id -> (barcode, umi, gene, transcript); any element may be None
ReadTags = Dict[str, Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]]

# Stand-ins for "nothing here" in the tables these tags are read from: "*" for a barcode or
# UMI that was not called, "." for a feature this code left unset, and the literal "None" that
# IsoQuant writes for a read with no assigned transcript (novel and ambiguous reads). A tag is
# omitted rather than carrying any of them.
PLACEHOLDERS = frozenset(("*", ".", "None"))


def index_bam(bam_file: str, threads: int = 1) -> None:
    """Index a BAM, falling back to CSI for references too long for BAI."""
    try:
        pysam.index('-@', str(threads), bam_file)
    except pysam.SamtoolsError as err:
        logger.info("Samtools failed to generate a .bai index: %s" % err)
        logger.info("Trying to build a CSI index instead")
        try:
            pysam.index('-@', str(threads), '-c', bam_file)
        except pysam.SamtoolsError as csi_err:
            logger.error("Failed to create CSI index: %s" % csi_err)


def merge_bam_files(output_bam: str, input_bams: Iterable[str], threads: int = 1,
                    remove_fragments: bool = True) -> Optional[str]:
    """Merge per-chromosome BAM fragments into one indexed BAM.

    The fragments share their header (each is written from the same template) and cover
    disjoint references, so a plain coordinate merge is correct.
    """
    fragments = [bam for bam in input_bams if bam and os.path.exists(bam)]
    if not fragments:
        logger.warning("No BAM fragments to merge into %s" % output_bam)
        return None

    if len(fragments) == 1:
        shutil.move(fragments[0], output_bam)
    else:
        pysam.merge("-f", "-@", str(threads), output_bam, *fragments)
        if remove_fragments:
            for fragment in fragments:
                os.remove(fragment)

    index_bam(output_bam, threads)
    logger.info("Merged %d BAM fragments into %s" % (len(fragments), output_bam))
    return output_bam


def load_survivor_tags(survivors_file: str) -> ReadTags:
    """Read a UMI-filtering survivors table into read_id -> (barcode, umi, gene, transcript).

    The table carries the tag columns only when a tagged BAM was requested; without them the
    read ids are still the definition of the surviving subset, so they map to empty tags.
    """
    read_tags = {}
    if not os.path.exists(survivors_file):
        return read_tags
    with open(survivors_file) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                read_tags[fields[0]] = (None, None, None, None)
            else:
                read_tags[fields[0]] = (fields[1], fields[2], fields[3], fields[4])
    return read_tags


def load_barcode_umi_tags(barcode_table: str) -> ReadTags:
    """Read a split barcode table into read_id -> (barcode, umi, None, None)."""
    read_tags = {}
    if not os.path.exists(barcode_table):
        return read_tags
    with open(barcode_table) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            read_tags[fields[0]] = (fields[1], fields[2], None, None)
    return read_tags


def _apply_tags(read, tags, barcode_tag: str, umi_tag: str) -> None:
    for tag_name, value in zip((barcode_tag, umi_tag, GENE_TAG, TRANSCRIPT_TAG), tags):
        if value and value not in PLACEHOLDERS:
            read.set_tag(tag_name, value)


def write_tagged_chromosome_bam(input_bams: Iterable[str], output_bam: str, chr_id: str,
                                read_tags: ReadTags, barcode_tag: str = "CB", umi_tag: str = "UB",
                                primary_only: bool = False, keep_untagged: bool = True) -> int:
    """Copy one chromosome's alignments into output_bam, adding tags from read_tags.

    primary_only drops secondary and supplementary records. keep_untagged decides whether a
    read absent from read_tags is copied without tags (a tagged copy of the input) or skipped
    (a subset restricted to the listed reads).
    """
    bam_list = list(input_bams)
    written = 0
    with pysam.AlignmentFile(bam_list[0], "rb") as template:
        with pysam.AlignmentFile(output_bam, "wb", template=template) as outf:
            for bam_file in bam_list:
                with pysam.AlignmentFile(bam_file, "rb") as inf:
                    for read in inf.fetch(chr_id):
                        if primary_only and (read.is_secondary or read.is_supplementary):
                            continue
                        tags = read_tags.get(read.query_name)
                        if tags is None and not keep_untagged:
                            continue
                        if tags is not None:
                            _apply_tags(read, tags, barcode_tag, umi_tag)
                        outf.write(read)
                        written += 1
    return written


def write_unmapped_bam(input_bams: Iterable[str], output_bam: str) -> int:
    """Copy unmapped reads, which fetch() by chromosome never returns.

    Only truly unplaced records (reference_id == -1) are taken: an unmapped read parked at a
    reference position does come back from fetch(), and copying it here too would duplicate it.
    """
    bam_list = list(input_bams)
    written = 0
    with pysam.AlignmentFile(bam_list[0], "rb") as template:
        with pysam.AlignmentFile(output_bam, "wb", template=template) as outf:
            for bam_file in bam_list:
                with pysam.AlignmentFile(bam_file, "rb") as inf:
                    for read in inf.fetch(until_eof=True):
                        if read.is_unmapped and read.reference_id == -1:
                            outf.write(read)
                            written += 1
    return written
