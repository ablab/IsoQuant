############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Detecting cell barcodes from extracted barcode counts.

A stock 10x whitelist has millions of entries but a run has a few thousand cells, and
matching every read against millions of candidates degenerates into exact matching --
`min_score` has to equal the barcode length for the result to mean anything, so any read
carrying a sequencing error in the barcode is lost.

Counting first turns the whitelist into a filter over a few thousand candidate cell
barcodes rather than a per-read search space. This module sits between the two barcode
calling passes: the first extracts barcode windows verbatim (TenXBarcodeDetector with
whitelist_matching=False), this decides which of them are cells, and the second is
ordinary barcode calling against that much shorter list.

The counting and selection are what remain of a port of Badger
(https://github.com/algbio/Badger); its edit-distance graph correction was dropped after
measuring it against the existing SSW matcher. See .claude/CELL_BARCODE_SELECTION.md.
"""

import logging
import math
from collections import defaultdict
from statistics import mean
from typing import Dict, Iterable, List, Optional, Sequence, Set, Union

import numpy

from .common import load_barcodes

logger = logging.getLogger('IsoQuant')

# Cell barcodes must be observed more often than mean(top n_cells counts) / this
CENTER_COUNT_FRACTION = 5.0
# ... and never less often than this, however sparse the data
MIN_CENTER_COUNT = 5
# Barcodes rarer than this are not considered when estimating the cell number
MIN_COUNT_FOR_ESTIMATION = 2
# How far below its chord the log-log count curve must bend before we call it a knee
KNEE_MIN_DEVIATION = 1e-9

NUCLEOTIDES = frozenset("ACGT")

NOSEQ = "*"
AUTO = "auto"


def estimate_cell_number(sorted_counts: Sequence[int]) -> int:
    """Estimate the number of cell-associated barcodes from the count distribution.

    Knee of the log-log count-versus-rank curve, located as the point furthest below the
    chord joining its ends. Used when --n_cells is "auto" or absent.
    """
    counts = [c for c in sorted_counts if c >= MIN_COUNT_FOR_ESTIMATION]
    if len(counts) < 10:
        return len(counts)

    x = numpy.log10(numpy.arange(1, len(counts) + 1, dtype=numpy.float64))
    y = numpy.log10(numpy.asarray(counts, dtype=numpy.float64))
    # distance from each point to the chord, signed so that points below it are positive
    dx, dy = x[-1] - x[0], y[-1] - y[0]
    norm = math.hypot(dx, dy)
    if norm == 0.0:
        return len(counts)
    distance = (dx * (y[0] - y) - (x[0] - x) * dy) / norm

    # A curve that never bends away from its chord has no knee to find -- every barcode
    # above the noise floor looks equally cell-like, so take all of them. Without this the
    # argmax of an all-zero array would silently report a single cell.
    if distance.max() <= KNEE_MIN_DEVIATION:
        return len(counts)
    return int(numpy.argmax(distance)) + 1


class CellBarcodeSelector:
    """Counts extracted barcodes and picks the ones belonging to real cells."""

    def __init__(self, barcode_length: int = 16):
        self.barcode_length: int = barcode_length
        self.counts: Dict[str, int] = defaultdict(int)
        self.centers: List[str] = []
        self.malformed_count: int = 0

    def add_barcode(self, barcode: str) -> bool:
        """Tally one extracted barcode. Returns False if it was unusable."""
        if len(barcode) == self.barcode_length + 1:
            barcode = barcode[:-1]
        if len(barcode) != self.barcode_length or not set(barcode) <= NUCLEOTIDES:
            self.malformed_count += 1
            return False
        self.counts[barcode] += 1
        return True

    def add_barcodes(self, barcodes: Iterable[str]) -> None:
        for barcode in barcodes:
            self.add_barcode(barcode)

    def merge_counts(self, counts: Dict[str, int], malformed: int = 0) -> None:
        """Merge a partial count table, e.g. produced by another process."""
        for barcode, count in counts.items():
            self.counts[barcode] += count
        self.malformed_count += malformed

    def sorted_barcodes(self) -> List[str]:
        """Distinct extracted barcodes, most frequent first."""
        return sorted(self.counts, key=lambda bc: (self.counts[bc], bc), reverse=True)

    def select(self, n_cells: Optional[int] = None,
               whitelist: Optional[Set[str]] = None,
               interval: int = 25) -> List[str]:
        """Pick the barcodes that represent actually sequenced cells."""
        by_count = self.sorted_barcodes()
        # Restrict to the candidate pool before anything else: an abundant barcode outside
        # the whitelist (ambient RNA, a chimera, a mis-anchored window) would otherwise
        # inflate the count cutoff below and push genuine cells under it.
        candidates = [bc for bc in by_count if whitelist is None or bc in whitelist]
        if not candidates:
            self.centers = []
            return []

        if n_cells is None:
            n_cells = estimate_cell_number([self.counts[bc] for bc in candidates])
            logger.info("Estimated number of cell-associated barcodes: %d" % n_cells)
        n_cells = max(1, min(n_cells, len(candidates)))

        cutoff = max(mean(self.counts[bc] for bc in candidates[:n_cells]) / CENTER_COUNT_FRACTION,
                     MIN_CENTER_COUNT)
        upper = n_cells + n_cells * interval / 100.0
        lower = n_cells - n_cells * interval / 100.0
        logger.info("Selecting up to %d cell barcodes out of %d candidates, minimal read count %.1f" %
                    (upper, len(candidates), cutoff))

        centers: List[str] = []
        i = 0
        while i < len(candidates) and self.counts[candidates[i]] > cutoff and len(centers) <= upper:
            centers.append(candidates[i])
            i += 1
        # If the cutoff was too aggressive, keep going until the lower bound is reached -- but
        # never down into the singletons, which is what injects spurious cells when n_cells is
        # far too high.
        while i < len(candidates) and len(centers) < lower and self.counts[candidates[i]] >= MIN_CENTER_COUNT:
            centers.append(candidates[i])
            i += 1

        logger.info("Detected %d cell barcodes, rarest one observed %d times" %
                    (len(centers), self.counts[centers[-1]] if centers else 0))
        if centers and len(centers) < lower:
            logger.warning("Detected %d cell barcodes, fewer than the %d expected. Check --n_cells, "
                           "or whether the whitelist matches the protocol." % (len(centers), n_cells))
        self.centers = centers
        return centers

    def stats(self) -> Dict[str, int]:
        reads_in_cells = sum(self.counts[bc] for bc in self.centers)
        return {
            "Distinct barcodes extracted": len(self.counts),
            "Malformed barcodes skipped": self.malformed_count,
            "Cell barcodes detected": len(self.centers),
            "Reads exactly matching a cell barcode": reads_in_cells,
        }


def load_whitelist(barcode_whitelist: Optional[List[str]]) -> Optional[Set[str]]:
    """Load the candidate pool cell barcodes are selected from, if one was supplied."""
    if not barcode_whitelist or barcode_whitelist == [AUTO]:
        return None
    whitelist: Set[str] = set()
    for file_name in barcode_whitelist:
        whitelist.update(load_barcodes(file_name))
    logger.info("Loaded %d whitelisted barcodes to select cell barcodes from" % len(whitelist))
    return whitelist


def select_cell_barcodes(selector: "CellBarcodeSelector", output_file: str,
                         barcode_whitelist: Optional[List[str]],
                         n_cells: Union[int, str, None] = AUTO,
                         n_cells_interval: int = 25,
                         stats_file: Optional[str] = None) -> str:
    """Pick the cell barcodes out of a filled count table and write them to output_file."""
    logger.info("Detecting cell barcodes from extracted barcode counts")
    if not selector.counts:
        logger.warning("No barcodes were extracted, cannot detect cell barcodes")

    whitelist = load_whitelist(barcode_whitelist)
    if whitelist is None:
        # Counts alone cannot tell a cell from a recurring extraction artifact: a
        # mis-anchored barcode window repeats across reads and looks exactly like an
        # abundant cell. A whitelist rejects those, because an artifact is not a valid
        # protocol barcode. Measured on ONT cDNA R10.4: 8 such artifacts among 5008
        # detected barcodes cost 5 points of precision.
        logger.warning("No barcode whitelist given, cell barcodes are selected on read counts alone. "
                       "Supplying the protocol whitelist as a candidate pool is more accurate.")
    requested_cells = None if n_cells == AUTO else n_cells
    centers = selector.select(n_cells=requested_cells, whitelist=whitelist,
                              interval=n_cells_interval)

    with open(output_file, "w") as out:
        for barcode in centers:
            out.write("%s\n" % barcode)
    logger.info("Detected %d cell barcodes, stored in %s" % (len(centers), output_file))

    stats = selector.stats()
    for key, value in stats.items():
        logger.info("  %s: %d" % (key, value))
    if stats_file:
        with open(stats_file, "w") as out_stats:
            for key, value in stats.items():
                out_stats.write("%s\t%d\n" % (key, value))
    return output_file
