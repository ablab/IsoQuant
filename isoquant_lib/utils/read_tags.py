############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Read tags (e.g. Dorado `pt:i`) that are carried from unaligned reads into IsoQuant's own alignments.

Tags come from an unaligned BAM, or from FASTA/FASTQ headers whose comment consists of SAM tags
(as written by `dorado basecaller --emit-fastq` or `samtools fastq -T`). Tags minimap2 writes itself
are dropped to avoid duplicates (Dorado's `ts:i` would clash with minimap2's `ts:A`), as are tags
that refer to header lines (`RG`, `PG`), since the header of the unaligned file is not carried over.
"""

import gzip
import re
from typing import IO, Iterable, List, Optional, Set, Tuple

import pysam

TAG_SAMPLE_SIZE = 1000

MINIMAP2_TAGS = frozenset(["NM", "MD", "AS", "SA", "ms", "nn", "ts", "tp", "cm", "s1", "s2",
                           "de", "dv", "rl", "cg", "cs", "ds", "zd"])
HEADER_BOUND_TAGS = frozenset(["RG", "PG"])

_TAG_VALUE_PATTERNS = {
    "A": r"[!-~]",
    "i": r"[-+]?\d+",
    "f": r"[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?|[-+]?(inf|nan)",
    "Z": r"[ !-~]*",
    "H": r"([0-9A-Fa-f][0-9A-Fa-f])*",
    "B": r"[cCsSiIf](,[-+0-9.eE]+)*",
}
_SAM_TAG_RE = re.compile(r"[A-Za-z][A-Za-z0-9]:(?:%s)" %
                         "|".join("%s:(?:%s)" % (t, p) for t, p in _TAG_VALUE_PATTERNS.items()))


def tags_from_comment(comment: str) -> Optional[Set[str]]:
    """Tag names in a read header comment; None when the comment is not entirely tab-separated SAM tags."""
    tags = set()
    if not comment:
        return tags
    for field in comment.split("\t"):
        if not _SAM_TAG_RE.fullmatch(field):
            return None
        tags.add(field[:2])
    return tags


def split_header(header_line: str) -> Tuple[str, str]:
    """Read name and comment of a FASTA/FASTQ header line."""
    parts = re.split(r"[ \t]", header_line[1:].rstrip("\r\n"), maxsplit=1)
    return parts[0], parts[1].lstrip(" \t") if len(parts) > 1 else ""


def tags_to_keep(tags: Iterable[str]) -> List[str]:
    return sorted(t for t in tags if t not in MINIMAP2_TAGS and t not in HEADER_BOUND_TAGS)


def _open_fastx(path: str) -> IO[str]:
    with open(path, "rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "r")


def _fastx_headers(handle: IO[str], sample_size: int) -> Iterable[str]:
    """Header lines of the first reads; FASTQ is expected in the usual 4-line layout."""
    first = handle.readline()
    if not first:
        return
    if first.startswith(">"):
        yield first
        count = 1
        for line in handle:
            if count >= sample_size:
                return
            if line.startswith(">"):
                yield line
                count += 1
    elif first.startswith("@"):
        yield first
        for count, line in enumerate(handle, start=1):
            if count // 4 >= sample_size:
                return
            if count % 4 == 0:
                yield line


def is_fasta(path: str) -> bool:
    with _open_fastx(path) as handle:
        return handle.read(1) == ">"


def fastx_tags_to_keep(path: str, sample_size: int = TAG_SAMPLE_SIZE) -> List[str]:
    """Tags to carry over from a FASTA/FASTQ file; empty when the first reads have no SAM-tag comments,
    or when any of them has a comment that is not made of SAM tags."""
    found: Set[str] = set()
    with _open_fastx(path) as handle:
        for header in _fastx_headers(handle, sample_size):
            tags = tags_from_comment(split_header(header)[1])
            if tags is None:
                return []
            found.update(tags)
    return tags_to_keep(found)


def bam_tags_to_keep(path: str, sample_size: int = TAG_SAMPLE_SIZE) -> List[str]:
    """Tags to carry over from an unaligned BAM, collected from its first reads."""
    found: Set[str] = set()
    # unaligned BAMs have no index, which htslib would report as an error on opening
    verbosity = pysam.set_verbosity(0)
    try:
        with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
            for i, read in enumerate(bam):
                if i >= sample_size:
                    break
                found.update(tag for tag, _ in read.get_tags())
    finally:
        pysam.set_verbosity(verbosity)
    return tags_to_keep(found)


def count_tagged_alignments(path: str, tag: str, sample_size: int = TAG_SAMPLE_SIZE) -> Tuple[int, int]:
    """Number of primary mapped alignments among the first sample_size that carry the tag, and how many were checked."""
    checked = 0
    tagged = 0
    with pysam.AlignmentFile(path, "rb") as bam:
        for alignment in bam.fetch(until_eof=True):
            if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                continue
            checked += 1
            if alignment.has_tag(tag):
                tagged += 1
            if checked >= sample_size:
                break
    return tagged, checked
