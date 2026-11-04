############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""CAGE peaks indexed for nearest-peak distance queries."""

from functools import lru_cache
from typing import Dict, List, Optional, Tuple

import numpy as np

from isoquant_lib.utils.file_utils import open_text_read

NO_PEAK_DISTANCE = np.iinfo(np.int64).max
BED_STRAND_COLUMN = 5
BED_HEADER_PREFIXES = ('#', 'track', 'browser')


class CageFormatError(ValueError):
    """The CAGE file is not a BED file with a strand column."""


def _read_bed(bed_file: str) -> Dict[Tuple[str, str], List[Tuple[int, int]]]:
    peaks: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    with open_text_read(bed_file) as bed:
        for line_number, line in enumerate(bed, 1):
            if not line.strip() or line.startswith(BED_HEADER_PREFIXES):
                continue
            fields = line.rstrip('\r\n').split('\t')
            where = "%s, line %d: " % (bed_file, line_number)
            if len(fields) <= BED_STRAND_COLUMN:
                raise CageFormatError(where + "expected a tab-separated BED with at least 6 columns "
                                              "(strand is required), found %d" % len(fields))
            strand = fields[BED_STRAND_COLUMN]
            if strand not in ('+', '-', '.'):
                raise CageFormatError(where + "unexpected strand '%s'" % strand)
            try:
                # BED starts are 0-based, IsoQuant positions and BED ends are 1-based inclusive
                interval = (int(fields[1]) + 1, int(fields[2]))
            except ValueError:
                raise CageFormatError(where + "peak coordinates are not integers")
            peaks.setdefault((fields[0], strand), []).append(interval)
    if not peaks:
        raise CageFormatError(bed_file + " contains no peaks")
    return peaks


class CagePeaks:
    """Peaks of a BED file, kept per (chromosome, strand) and sorted by start.

    A running maximum of peak ends makes overlapping peaks need no special
    handling: of all peaks starting at or before a position, only the one
    reaching furthest can be the closest."""

    def __init__(self, bed_file: str) -> None:
        self._peaks: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]] = {}
        for key, intervals in _read_bed(bed_file).items():
            peaks = np.array(sorted(intervals), dtype=np.int64)
            self._peaks[key] = (peaks[:, 0], np.maximum.accumulate(peaks[:, 1]))

    def has_chromosome(self, chrom: str) -> bool:
        return any(chrom == peak_chrom for peak_chrom, _ in self._peaks)

    def distances(self, chrom: str, strand: str, positions: np.ndarray) -> Optional[np.ndarray]:
        """Distance from each position to the closest peak on the same strand,
        0 inside a peak. None if this chromosome strand has no peaks."""
        peaks = self._peaks.get((chrom, strand))
        if peaks is None:
            return None
        starts, max_ends = peaks
        positions = np.asarray(positions, dtype=np.int64)
        next_peak = np.searchsorted(starts, positions, side='right')
        to_previous = np.where(next_peak > 0,
                               positions - max_ends[np.maximum(next_peak - 1, 0)], NO_PEAK_DISTANCE)
        to_next = np.where(next_peak < len(starts),
                           starts[np.minimum(next_peak, len(starts) - 1)] - positions, NO_PEAK_DISTANCE)
        return np.maximum(np.minimum(to_previous, to_next), 0)


@lru_cache(maxsize=1)
def load_cage_peaks(bed_file: str) -> CagePeaks:
    """Parse the file once per process, it is shared by all its counters and chromosomes."""
    return CagePeaks(bed_file)
