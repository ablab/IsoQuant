############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Cell barcode detection.

Runs between the two barcode calling passes. The first pass extracts barcode windows
verbatim (TenXBarcodeDetector with whitelist_matching=False); this stage counts them and
decides which barcodes belong to real cells; the second pass is ordinary barcode calling
against that much smaller list.

The point is that a stock 10x whitelist has millions of entries but a run has a few
thousand cells, and matching every read against millions of candidates degenerates into
exact matching. Counting first turns the whitelist into a filter over a few thousand
candidate cell barcodes instead of a per-read search space.

See .claude/CELL_BARCODE_SELECTION.md.
"""

import logging
from typing import Dict, List, Optional, Set, Union

from .barcode_graph import BarcodeGraph
from .common import load_barcodes

logger = logging.getLogger('IsoQuant')

NOSEQ = "*"
READ_ID_COLUMN = 0
BARCODE_COLUMN = 1
AUTO = "auto"


def _is_header(fields: List[str]) -> bool:
    """Raw tables carry a plain (non-commented) header; a barcode is never "barcode"."""
    return fields[BARCODE_COLUMN] == "barcode" or fields[READ_ID_COLUMN].startswith("#")


def count_barcodes(raw_barcode_files: List[str], graph: BarcodeGraph) -> int:
    """Tally every extracted barcode across all raw tables. Returns the number of rows read."""
    total_rows = 0
    for file_name in raw_barcode_files:
        rows = 0
        with open(file_name) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= BARCODE_COLUMN or _is_header(fields):
                    continue
                rows += 1
                if fields[BARCODE_COLUMN] != NOSEQ:
                    graph.add_barcode(fields[BARCODE_COLUMN])
        logger.info("Read %d extracted barcodes from %s" % (rows, file_name))
        total_rows += rows
    return total_rows


def load_whitelist(barcode_whitelist: Optional[List[str]]) -> Optional[Set[str]]:
    """Load the candidate pool cell barcodes are selected from, if one was supplied."""
    if not barcode_whitelist or barcode_whitelist == [AUTO]:
        return None
    whitelist: Set[str] = set()
    for file_name in barcode_whitelist:
        whitelist.update(load_barcodes(file_name))
    logger.info("Loaded %d whitelisted barcodes to select cell barcodes from" % len(whitelist))
    return whitelist


def select_cell_barcodes(raw_barcode_files: List[str], output_file: str,
                         barcode_whitelist: Optional[List[str]],
                         barcode_length: int = 16,
                         n_cells: Union[int, str, None] = AUTO,
                         n_cells_interval: int = 25,
                         stats_file: Optional[str] = None) -> str:
    """Detect cell barcodes from extracted barcode counts and write them to output_file."""
    logger.info("Detecting cell barcodes from extracted barcode counts")
    graph = BarcodeGraph(barcode_length=barcode_length)

    count_barcodes(raw_barcode_files, graph)
    if not graph.counts:
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
    centers = graph.select_cluster_centers(n_cells=requested_cells, whitelist=whitelist,
                                           interval=n_cells_interval)

    with open(output_file, "w") as out:
        for barcode in centers:
            out.write("%s\n" % barcode)
    logger.info("Detected %d cell barcodes, stored in %s" % (len(centers), output_file))

    stats = selection_stats(graph)
    for key, value in stats.items():
        logger.info("  %s: %d" % (key, value))
    if stats_file:
        with open(stats_file, "w") as out_stats:
            for key, value in stats.items():
                out_stats.write("%s\t%d\n" % (key, value))
    return output_file


def selection_stats(graph: BarcodeGraph) -> Dict[str, int]:
    reads_in_cells = sum(graph.counts[bc] for bc in graph.centers)
    return {
        "Distinct barcodes extracted": len(graph.counts),
        "Malformed barcodes skipped": graph.malformed_count,
        "Cell barcodes detected": len(graph.centers),
        "Reads exactly matching a cell barcode": reads_in_cells,
    }
