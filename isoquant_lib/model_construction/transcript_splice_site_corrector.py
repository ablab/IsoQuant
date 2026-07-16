############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Splice-site error correction for constructed transcript models.
#
# Noisy long reads (especially ONT) carry clustered deletions right next to
# splice junctions, so the aligner frequently places a junction a few bases off
# from the true canonical site. For each constructed (novel) transcript model we
# look at the reads assigned to it and, at every exon boundary that falls inside a
# read's aligned span, count deleted reference bases in a small window next to the
# junction. If a distinct deletion shift (of an accepted size, supported by enough
# reads) would move the boundary onto a canonical splice motif, the boundary is
# shifted.
#
# Ported from PR #108 ("Keep cigartuples") and reworked so that we do NOT keep
# full per-read CIGARs in memory: a model splice site is always at an exon
# boundary and the reads assigned to it share that junction structure, so only
# deletions near a read's own exon boundaries can ever land in a junction window.
# We therefore store, per read, just those boundary-proximal deletions as
# (ref_pos, length) pairs (usually an empty list) and count deleted bases in
# reference space. Canonical motifs reuse the verified convention from
# common.get_intron_strand (CANONICAL_FWD_SITES / CANONICAL_REV_SITES).

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from isoquant_lib.common import CANONICAL_FWD_SITES, CANONICAL_REV_SITES

logger = logging.getLogger('IsoQuant')


@dataclass  # not slots=True: CI runs Python 3.8
class SpliceSiteCase:
    """Per-splice-site accumulator: whether the boundary is an exon end, the
    histogram of per-read deleted-base counts in the adjacent window, and the
    resolved most common deletion shift (filled by correct_splice_site_errors)."""
    location_is_end: bool
    deletions: Dict[int, int] = field(default_factory=dict)
    most_common_del: int = -1

# --- Tunable constants ---------------------------------------------------------
# Accepted absolute deletion shifts (in bp) considered real splice-site errors.
ACCEPTED_DEL_CASES: Tuple[int, ...] = (3, 4, 5, 6)
# Number of reference bases inspected next to each splice site.
WINDOW_SIZE: int = 8
# Minimum number of reads spanning a splice site before it is considered.
MIN_N_OF_ALIGNED_READS: int = 5

SUPPORTED_STRANDS: Tuple[str, ...] = ('+', '-')

# A deletion can only affect a splice-site window if it lies within roughly a
# window plus the largest accepted shift of a read exon boundary; deletions
# further inside an exon are irrelevant and are dropped at collection time.
BOUNDARY_MARGIN: int = WINDOW_SIZE + max(ACCEPTED_DEL_CASES)

# CIGAR operations that consume the reference (M, D, N, =, X); D is singled out.
_BAM_CDEL: int = 2
_REF_CONSUMING_OPS: frozenset = frozenset((0, 2, 3, 7, 8))

# Valid splice dinucleotides derived from the verified canonical sets in common.
# A boundary that is an exon *end* determines the downstream intron's left
# (donor-side) dinucleotide; an exon *start* determines the upstream intron's
# right (acceptor-side) dinucleotide. The strand is already baked into the sets.
_FWD_LEFT_SITES: frozenset = frozenset(donor for donor, _ in CANONICAL_FWD_SITES)
_FWD_RIGHT_SITES: frozenset = frozenset(acceptor for _, acceptor in CANONICAL_FWD_SITES)
_REV_LEFT_SITES: frozenset = frozenset(donor for donor, _ in CANONICAL_REV_SITES)
_REV_RIGHT_SITES: frozenset = frozenset(acceptor for _, acceptor in CANONICAL_REV_SITES)


def extract_relevant_deletions(cigartuples: List[Tuple[int, int]],
                               ref_start: int,
                               exons: List[Tuple[int, int]],
                               margin: int = BOUNDARY_MARGIN) -> List[Tuple[int, int]]:
    """Return the boundary-proximal deletions of an alignment as (ref_pos, length)
    pairs (1-based reference coordinates of the first deleted base).

    ``ref_start`` is the 1-based reference position of the first aligned base
    (i.e. ``read_exons[0][0]`` / ``alignment.reference_start + 1``). Only
    deletions within ``margin`` of some exon boundary are kept; interior
    deletions can never fall in a splice-site window."""
    if not cigartuples or not exons:
        return []
    boundaries = [b for exon in exons for b in exon]
    deletions: List[Tuple[int, int]] = []
    ref_pos = ref_start
    for op, length in cigartuples:
        if op == _BAM_CDEL:
            del_end = ref_pos + length - 1
            if any(ref_pos - margin <= b <= del_end + margin for b in boundaries):
                deletions.append((ref_pos, length))
            ref_pos += length
        elif op in _REF_CONSUMING_OPS:
            ref_pos += length
        # insertions / soft/hard clips do not consume the reference
    return deletions


def extract_splice_site_locations_within_aligned_read(read_start: int,
                                                      read_end: int,
                                                      exons: List[Tuple[int, int]]
                                                      ) -> List[Tuple[int, bool]]:
    """Return (location, location_is_end) for every exon boundary that falls
    within [read_start, read_end]. A read may start/end in the middle of an exon,
    so both boundaries of each exon are tested independently."""
    matching_locations: List[Tuple[int, bool]] = []
    for exon_start, exon_end in exons:
        if read_start <= exon_start <= read_end:
            matching_locations.append((exon_start, False))
        if read_start <= exon_end <= read_end:
            matching_locations.append((exon_end, True))
        if read_end <= exon_end:
            break
    return matching_locations


def _deleted_bases_in_window(deletions: List[Tuple[int, int]], low: int, high: int) -> int:
    """Number of deleted reference bases within the inclusive window [low, high]."""
    total = 0
    for del_start, del_len in deletions:
        del_end = del_start + del_len - 1
        overlap = min(high, del_end) - max(low, del_start) + 1
        if overlap > 0:
            total += overlap
    return total


def count_deletions_for_splice_site_locations(read_start: int,
                                              read_end: int,
                                              deletions: List[Tuple[int, int]],
                                              exons: List[Tuple[int, int]],
                                              splice_site_cases: Dict[int, "SpliceSiteCase"],
                                              window_size: int = WINDOW_SIZE) -> None:
    """Accumulate, per splice site, a histogram of the number of deleted bases a
    single read carries in the window next to that boundary."""
    for splice_site_location, location_is_end in \
            extract_splice_site_locations_within_aligned_read(read_start, read_end, exons):
        if location_is_end:
            low, high = splice_site_location - window_size + 1, splice_site_location
        else:
            low, high = splice_site_location, splice_site_location + window_size - 1
        count_of_deletions = _deleted_bases_in_window(deletions, low, high)

        splice_site_data = splice_site_cases.get(splice_site_location)
        if splice_site_data is None:
            splice_site_data = SpliceSiteCase(location_is_end=location_is_end)
            splice_site_cases[splice_site_location] = splice_site_data
        deletion_counts = splice_site_data.deletions
        deletion_counts[count_of_deletions] = deletion_counts.get(count_of_deletions, 0) + 1


def compute_most_common_case_of_deletions(deletions: Dict[int, int], location_is_end: bool) -> int:
    """Return the single most common deletion count, signed by direction (negative
    for exon ends, positive for exon starts). Returns -1 when there is no unique
    mode."""
    max_count = max(deletions.values())
    modes = [k for k, v in deletions.items() if v == max_count]
    if len(modes) != 1:
        return -1
    return -modes[0] if location_is_end else modes[0]


def _corrected_site_is_canonical(chr_record,
                                 boundary: int,
                                 most_common_del: int,
                                 location_is_end: bool,
                                 strand: str) -> bool:
    """Check whether shifting a 1-based exon boundary by ``most_common_del`` lands
    the adjacent intron dinucleotide on a canonical splice motif for ``strand``.

    Uses the same indexing convention as common.get_intron_strand: for an intron
    starting at 1-based position S the left dinucleotide is chr_record[S-1:S+1];
    for an intron ending at 1-based position E the right dinucleotide is
    chr_record[E-2:E].
    """
    corrected = boundary + most_common_del
    if location_is_end:
        # Downstream intron starts at corrected + 1 -> left (donor-side) site.
        start = corrected            # 0-based index of 1-based (corrected + 1)
        valid = _FWD_LEFT_SITES if strand == '+' else _REV_LEFT_SITES
    else:
        # Upstream intron ends at corrected - 1 -> right (acceptor-side) site.
        start = corrected - 3        # 0-based index of 1-based (corrected - 2)
        valid = _FWD_RIGHT_SITES if strand == '+' else _REV_RIGHT_SITES
    if start < 0:
        return False
    try:
        site = str(chr_record[start:start + 2]).upper()
    except (KeyError, IndexError):
        return False
    return site in valid


def correct_splice_site_errors(splice_site_cases: Dict[int, "SpliceSiteCase"],
                               chr_record,
                               strand: str,
                               min_n_of_aligned_reads: int = MIN_N_OF_ALIGNED_READS,
                               accepted_del_cases: Tuple[int, ...] = ACCEPTED_DEL_CASES
                               ) -> List[int]:
    """Return the splice-site locations whose most common deletion shift is an
    accepted size and lands on a canonical motif. As a side effect fills in
    ``most_common_del`` for each qualifying case (used when rebuilding exons)."""
    locations_with_errors: List[int] = []
    for splice_site_location, splice_site_data in splice_site_cases.items():
        if sum(splice_site_data.deletions.values()) < min_n_of_aligned_reads:
            continue

        most_common_del = compute_most_common_case_of_deletions(
            splice_site_data.deletions, splice_site_data.location_is_end)
        splice_site_data.most_common_del = most_common_del
        if abs(most_common_del) not in accepted_del_cases:
            continue

        if _corrected_site_is_canonical(chr_record, splice_site_location, most_common_del,
                                        splice_site_data.location_is_end, strand):
            locations_with_errors.append(splice_site_location)

    return locations_with_errors


def generate_updated_exon_list(splice_site_cases: Dict[int, "SpliceSiteCase"],
                               locations_with_errors: List[int],
                               exons: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Rebuild the exon list applying the corrections in ``locations_with_errors``."""
    corrections = set(locations_with_errors)
    updated_exons: List[Tuple[int, int]] = []
    for exon_start, exon_end in exons:
        if exon_start in corrections:
            exon_start = exon_start + splice_site_cases[exon_start].most_common_del
        if exon_end in corrections:
            exon_end = exon_end + splice_site_cases[exon_end].most_common_del
        updated_exons.append((exon_start, exon_end))
    return updated_exons
