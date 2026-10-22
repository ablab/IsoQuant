############################################################################
# Copyright (c) 2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Barcode-to-spot mapping used during Visium HD barcode calling.

A Visium HD barcode is a pair of parts, PART1|PART2, identifying a 2 um spot. With --barcode2barcode
the barcodes are reported as spot IDs of a coarser resolution, so a read whose part barcodes tie
between neighbours can still be assigned when all tied combinations fall into one spot.
"""

import gzip
import logging
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy

from ..assignment.read_groups import parse_barcode2spot_spec
from ..utils.error_codes import IsoQuantExitCode

logger = logging.getLogger('IsoQuant')

JOINT_BARCODE_SEPARATOR = "|"
UNMAPPED = -1


def parse_spot_mapping_spec(spec: str, option_name: str = "--barcode2barcode") -> Tuple[str, int, int]:
    """Parse FILE[:barcode_col:spot_col]; barcodes can be called at one resolution only."""
    filename, barcode_col, spot_cols = parse_barcode2spot_spec(spec)
    if len(spot_cols) != 1:
        logger.critical("%s selects %d spot columns (%s), but barcodes can be called at one resolution only; "
                        "select a single column, e.g. %s:%d:%d"
                        % (option_name, len(spot_cols), ",".join(map(str, spot_cols)),
                           filename, barcode_col, spot_cols[0]))
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
    return filename, barcode_col, spot_cols[0]


def _index_barcodes(barcodes: Sequence[str]) -> Dict[str, int]:
    index: Dict[str, int] = {}
    for barcode in barcodes:
        index.setdefault(barcode, len(index))
    return index


class SpotMapping:
    """PART1|PART2 -> spot ID, stored as a part-index matrix of spot indices.

    A dict keyed by composite barcode would take gigabytes for a full slide (11.2M barcodes) and be
    copied into every barcode calling worker; the matrix takes 4 bytes per barcode pair.
    """

    def __init__(self, part1_barcodes: Sequence[str], part2_barcodes: Sequence[str],
                 mapping_file: str, barcode_col: int = 0, spot_col: int = 1):
        self.part1_index = _index_barcodes(part1_barcodes)
        self.part2_index = _index_barcodes(part2_barcodes)
        self.part2_count = len(self.part2_index)
        self.spot_index = numpy.full(len(self.part1_index) * self.part2_count, UNMAPPED, dtype=numpy.int32)
        self.spot_names: List[str] = []
        self._load(mapping_file, barcode_col, spot_col)

    def _load(self, mapping_file: str, barcode_col: int, spot_col: int) -> None:
        name_to_index: Dict[str, int] = {}
        loaded = skipped = 0
        handle = gzip.open(mapping_file, "rt") if mapping_file.endswith(".gz") else open(mapping_file)
        with handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                columns = line.rstrip("\n").split("\t")
                if len(columns) <= max(barcode_col, spot_col) or not columns[spot_col]:
                    skipped += 1
                    continue
                parts = columns[barcode_col].split(JOINT_BARCODE_SEPARATOR)
                if len(parts) != 2 or parts[0] not in self.part1_index or parts[1] not in self.part2_index:
                    skipped += 1
                    continue
                spot = columns[spot_col]
                spot_index = name_to_index.get(spot)
                if spot_index is None:
                    spot_index = len(self.spot_names)
                    name_to_index[spot] = spot_index
                    self.spot_names.append(spot)
                self.spot_index[self.part1_index[parts[0]] * self.part2_count + self.part2_index[parts[1]]] = spot_index
                loaded += 1

        logger.info("Loaded %d barcode-to-spot mappings into %d spots from %s"
                    % (loaded, len(self.spot_names), mapping_file))
        if skipped:
            logger.warning("Skipped %d lines of %s whose barcode is not PART1%sPART2 of the whitelists "
                           "or has no spot ID" % (skipped, mapping_file, JOINT_BARCODE_SEPARATOR))
        if not loaded:
            logger.critical("No barcode in %s matches the barcode whitelists" % mapping_file)
            sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)

    def _spot_index_of(self, barcode1: str, barcode2: str) -> int:
        index1 = self.part1_index.get(barcode1)
        index2 = self.part2_index.get(barcode2)
        if index1 is None or index2 is None:
            return UNMAPPED
        return int(self.spot_index[index1 * self.part2_count + index2])

    def spot(self, barcode1: str, barcode2: str) -> Optional[str]:
        """Spot of a barcode pair, None when the pair is not in the mapping."""
        spot_index = self._spot_index_of(barcode1, barcode2)
        return None if spot_index == UNMAPPED else self.spot_names[spot_index]

    def unique_spot(self, barcodes1: Iterable[str], barcodes2: Iterable[str]) -> Optional[str]:
        """The spot shared by every mapped pair of candidates, None if there is none or more than one."""
        barcodes2 = list(barcodes2)
        found = UNMAPPED
        for barcode1 in barcodes1:
            for barcode2 in barcodes2:
                spot_index = self._spot_index_of(barcode1, barcode2)
                if spot_index == UNMAPPED:
                    continue
                if found == UNMAPPED:
                    found = spot_index
                elif spot_index != found:
                    return None
        return None if found == UNMAPPED else self.spot_names[found]
