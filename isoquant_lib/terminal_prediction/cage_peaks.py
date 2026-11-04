############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""CAGE peaks indexed for nearest-peak distance queries."""

from functools import lru_cache
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

NO_PEAK_DISTANCE = np.iinfo(np.int64).max


class CagePeaks:
    """Peaks of a BED file, kept per (chromosome, strand) and sorted by start.

    A running maximum of peak ends makes overlapping peaks need no special
    handling: of all peaks starting at or before a position, only the one
    reaching furthest can be the closest."""

    def __init__(self, bed_file: str) -> None:
        self._peaks: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]] = {}
        table = pd.read_csv(bed_file, sep='\t', header=None, usecols=[0, 1, 2, 5],
                            names=['chrom', 'start', 'end', 'strand'], comment='#')
        for key, peaks in table.groupby(['chrom', 'strand']):
            peaks = peaks.sort_values('start')
            self._peaks[key] = (peaks['start'].values.astype(np.int64),
                                np.maximum.accumulate(peaks['end'].values.astype(np.int64)))

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
