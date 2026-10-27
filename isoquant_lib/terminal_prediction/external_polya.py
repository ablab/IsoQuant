############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""PolyA status from an external source (`--polya_trimmed tag:/list:/flnc:`).

Dorado (`pt:i` tag) and `isoseq refine` (`flnc.report.csv`) find polyA tails by anchoring on
adapters/primers, which is more reliable than the sequence scan in `PolyAFinder`. The sequence
scan still runs; where the two disagree, the external source wins (see `reconcile_polya`).
"""

import logging
from enum import Enum, unique
from typing import Dict, Iterator, List, Optional, Tuple

from isoquant_lib.terminal_prediction.polya_finder import PolyAInfo
from isoquant_lib.utils.file_utils import open_text_read

logger = logging.getLogger('IsoQuant')

LIST_SOURCE = "list"
FLNC_SOURCE = "flnc"

# second column of the normalized table
NO_HINT = "."
NO_TAIL = "0"
STRAND_HINTS = ("+", "-")

FLNC_ID_COLUMN = "id"
FLNC_POLYA_COLUMN = "polyAlen"
FLNC_STRAND_COLUMN = "strand"


@unique
class ExternalPolyAStatus(Enum):
    unknown = 0
    present = 1
    no_tail = 2


class ExternalPolyA:
    """What the external source says about one read; hint is the transcript strand relative to the read."""
    __slots__ = ("status", "hint")

    def __init__(self, status: ExternalPolyAStatus, hint: Optional[str] = None):
        self.status: ExternalPolyAStatus = status
        self.hint: Optional[str] = hint

    def __eq__(self, other) -> bool:
        return isinstance(other, ExternalPolyA) and self.status == other.status and self.hint == other.hint

    def __repr__(self) -> str:
        return "ExternalPolyA(%s, %s)" % (self.status.name, self.hint)


# shared instances, so that per-read dictionaries do not hold one object per read
UNKNOWN = ExternalPolyA(ExternalPolyAStatus.unknown)
PRESENT_NO_HINT = ExternalPolyA(ExternalPolyAStatus.present)
PRESENT_FORWARD = ExternalPolyA(ExternalPolyAStatus.present, "+")
PRESENT_REVERSE = ExternalPolyA(ExternalPolyAStatus.present, "-")
NO_TAIL_FOUND = ExternalPolyA(ExternalPolyAStatus.no_tail)

_FROM_TABLE_VALUE: Dict[str, ExternalPolyA] = {
    NO_HINT: PRESENT_NO_HINT,
    "+": PRESENT_FORWARD,
    "-": PRESENT_REVERSE,
    NO_TAIL: NO_TAIL_FOUND,
}


_PRESENT_WITH_HINT: Dict[Optional[str], ExternalPolyA] = {
    None: PRESENT_NO_HINT,
    "+": PRESENT_FORWARD,
    "-": PRESENT_REVERSE,
}


def hint_from_strandedness(stranded: Optional[str]) -> Optional[str]:
    """--stranded forward/reverse gives the transcript strand relative to the read for every read."""
    return {"forward": "+", "reverse": "-"}.get(stranded)


def external_polya_from_tag(alignment, tag: str, default_hint: Optional[str] = None,
                            strand_tag: str = "TS") -> ExternalPolyA:
    """Dorado style: tag >= 0 means a tail was found (0 = anchor found, length not estimated).
    -1 means the primer anchor was not found, which says nothing about the tail.
    TS:A (written only when Dorado primer trimming is on, normally absent for dRNA) gives the side;
    without it default_hint (from --stranded) is used."""
    if not alignment.has_tag(tag):
        return UNKNOWN
    value = alignment.get_tag(tag)
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        return UNKNOWN
    if alignment.has_tag(strand_tag):
        hint = alignment.get_tag(strand_tag)
        if hint == "+":
            return PRESENT_FORWARD
        if hint == "-":
            return PRESENT_REVERSE
    return _PRESENT_WITH_HINT[default_hint]


def resolve_external_strand(external: ExternalPolyA, mapped_strand: str) -> Optional[str]:
    """Genomic strand of the transcript, when the source gives the side of the tail."""
    if external.status != ExternalPolyAStatus.present or external.hint is None:
        return None
    if external.hint == "+":
        return mapped_strand
    return "-" if mapped_strand == "+" else "+"


def has_any_polya(polya_info: PolyAInfo) -> bool:
    return (polya_info.external_polya_pos != -1 or polya_info.external_polyt_pos != -1 or
            polya_info.internal_polya_pos != -1 or polya_info.internal_polyt_pos != -1)


def reconcile_polya(polya_info: PolyAInfo, external_strand: Optional[str], external_no_tail: bool) -> None:
    """Drop sequence-detected tails the external source contradicts.

    Runs before PolyAFixer, so a contradicted tail does not trim exons either.
    A tail on the side the source confirms is kept, since its position is exact.
    """
    if external_no_tail:
        polya_info.external_polya_pos = -1
        polya_info.internal_polya_pos = -1
        polya_info.external_polyt_pos = -1
        polya_info.internal_polyt_pos = -1
    elif external_strand == "+":
        polya_info.external_polyt_pos = -1
        polya_info.internal_polyt_pos = -1
    elif external_strand == "-":
        polya_info.external_polya_pos = -1
        polya_info.internal_polya_pos = -1


def inject_external_polya(polya_info: PolyAInfo, read_exons: List[Tuple[int, int]],
                          strand: Optional[str]) -> None:
    """Mark a tail at the read end of the given genomic strand, unless the sequence already found one there."""
    if strand == "+":
        if polya_info.external_polya_pos == -1 and polya_info.internal_polya_pos == -1:
            polya_info.external_polya_pos = read_exons[-1][1] + 1
    elif strand == "-":
        if polya_info.external_polyt_pos == -1 and polya_info.internal_polyt_pos == -1:
            polya_info.external_polyt_pos = read_exons[0][0] - 1


def _normalize_list_line(line: str) -> Optional[Tuple[str, str]]:
    if not line.strip() or line.startswith("#"):
        return None
    columns = line.split()
    hint = columns[1] if len(columns) > 1 and columns[1] in STRAND_HINTS else NO_HINT
    return columns[0], hint


def read_flnc_header(line: str) -> Tuple[str, int, int, Optional[int]]:
    """Returns delimiter, read id, polyAlen and strand (if present) columns; raises ValueError for a non-flnc header."""
    delim = "," if "," in line else "\t"
    columns = [c.strip() for c in line.rstrip("\n").split(delim)]
    if FLNC_ID_COLUMN not in columns or FLNC_POLYA_COLUMN not in columns:
        raise ValueError("header must contain '%s' and '%s' columns, found: %s" %
                         (FLNC_ID_COLUMN, FLNC_POLYA_COLUMN, line.strip()))
    strand_col = columns.index(FLNC_STRAND_COLUMN) if FLNC_STRAND_COLUMN in columns else None
    return delim, columns.index(FLNC_ID_COLUMN), columns.index(FLNC_POLYA_COLUMN), strand_col


def check_flnc_report(file_name: str) -> None:
    """Fail early on a file that is not an isoseq refine report; raises ValueError."""
    with open_text_read(file_name) as inf:
        read_flnc_header(inf.readline())


def _iterate_flnc(inf) -> Iterator[Tuple[str, str]]:
    delim, id_col, polya_col, strand_col = read_flnc_header(inf.readline())
    min_columns = max(id_col, polya_col, strand_col if strand_col is not None else 0) + 1
    for line in inf:
        columns = line.rstrip("\n").split(delim)
        if len(columns) < min_columns:
            continue
        try:
            polya_len = int(columns[polya_col])
        except ValueError:
            continue
        if polya_len <= 0:
            yield columns[id_col], NO_TAIL
            continue
        # FLNC reads are oriented, so the strand is normally '+'
        strand = columns[strand_col].strip() if strand_col is not None else "+"
        yield columns[id_col], strand if strand in STRAND_HINTS else NO_HINT


def normalize_polya_reads(src: str, source: str, dst: str) -> int:
    """Convert a read list or flnc report into `read_id<TAB>{+,-,.,0}`; returns the number of reads written."""
    count = 0
    with open_text_read(src) as inf, open(dst, "w") as outf:
        if source == FLNC_SOURCE:
            records = _iterate_flnc(inf)
        elif source == LIST_SOURCE:
            records = filter(None, map(_normalize_list_line, inf))
        else:
            raise ValueError("unknown external polyA source " + source)
        for read_id, value in records:
            outf.write("%s\t%s\n" % (read_id, value))
            count += 1

    if count == 0:
        logger.warning("No reads were read from %s, check the file format for --polya_trimmed %s:" % (src, source))
    else:
        logger.info("Loaded external polyA information for %d reads from %s" % (count, src))
    return count


def load_polya_read_dict(file_name: str) -> Dict[str, ExternalPolyA]:
    polya_dict: Dict[str, ExternalPolyA] = {}
    with open(file_name) as inf:
        for line in inf:
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 2:
                continue
            polya_dict[columns[0]] = _FROM_TABLE_VALUE.get(columns[1], PRESENT_NO_HINT)
    return polya_dict
