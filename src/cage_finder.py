############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import pybedtools as pbt


logger = logging.getLogger('IsoQuant')


class CagePeakFinder:
    def __init__(self, cage_file, shift_size=50, window_size=5):
        self.cage_peaks = self._load_cage_peaks(cage_file)
        self.shift_size = shift_size
        self.window_size = window_size

    def _load_cage_peaks(self, cage_file):
        return pbt.BedTool(cage_file)

    def _get_search_region(self, alignment, extended=False):
        contig = alignment.reference_name
        search_size = self.shift_size if extended else self.window_size
        if alignment.is_reverse:
            strand = '-'
            start = max(alignment.reference_end - self.window_size, 0)
            end = alignment.reference_end + search_size
        else:
            strand = '.'
            start = max(alignment.reference_start - search_size, 0)
            end = alignment.reference_start + self.window_size
        return contig, start, end, strand

    def find_cage_peak(self, alignment):
        logger.debug("Searching for cage peak for %s " % alignment.query_name)

        contig, start, end, strand = self._get_search_region(alignment, extended=True)
        alignment_interval = pbt.Interval(chrom=contig, start=start, end=end, strand=strand)
        cage_intersections = self.cage_peaks.all_hits(alignment_interval)

        if len(cage_intersections) > 0:
            logger.debug('CAGE peaks found: {}'.format(cage_intersections))
        else:
            logger.debug('No CAGE peaks found')

        return cage_intersections
