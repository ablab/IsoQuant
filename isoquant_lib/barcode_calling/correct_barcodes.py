############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Graph-based barcode correction stage.

Runs after raw barcode extraction (TenXBarcodeDetector with whitelist_matching=False).
Reads the raw barcode tables, decides which barcodes correspond to real cells, corrects
the rest onto them, and writes the usual barcoded read tables so that everything
downstream (dataset_processor.split_read_barcode_table and beyond) is unchanged.

See .claude/BARCODE_GRAPH_CORRECTION.md.
"""

import logging
import os
from typing import Dict, List, Optional, Set

from .barcode_graph import BarcodeGraph
from .common import load_barcodes

logger = logging.getLogger('IsoQuant')

NOSEQ = "*"
READ_ID_COLUMN = 0
BARCODE_COLUMN = 1
RAW_BARCODE_HEADER = "raw_barcode"


def _is_header(fields: List[str]) -> bool:
    """Raw tables carry a plain (non-commented) header; a barcode is never "barcode"."""
    return fields[BARCODE_COLUMN] == "barcode" or fields[READ_ID_COLUMN].startswith("#")


def _read_barcode_column(file_name: str):
    """Yield (fields, barcode) for every data row of a raw barcode table."""
    with open(file_name) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= BARCODE_COLUMN or _is_header(fields):
                continue
            yield fields, fields[BARCODE_COLUMN]


def count_barcodes(raw_barcode_files: List[str], graph: BarcodeGraph) -> int:
    """Tally every observed barcode across all raw tables. Returns the number of rows read."""
    total_rows = 0
    for file_name in raw_barcode_files:
        rows = 0
        for _, barcode in _read_barcode_column(file_name):
            rows += 1
            if barcode != NOSEQ:
                graph.add_barcode(barcode)
        logger.info("Read %d barcoded reads from %s" % (rows, os.path.basename(file_name)))
        total_rows += rows
    return total_rows


def write_corrected_table(raw_barcode_file: str, output_file: str,
                          assignments: Dict[str, str], barcode_length: int) -> Dict[str, int]:
    """Rewrite a raw table with corrected barcodes, keeping the raw one in a new column."""
    stats = {"corrected": 0, "unchanged": 0, "dropped": 0, "no_barcode": 0}
    with open(raw_barcode_file) as inf, open(output_file, "w") as outf:
        for line in inf:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= BARCODE_COLUMN:
                continue
            if _is_header(fields):
                outf.write("%s\t%s\n" % ("\t".join(fields), RAW_BARCODE_HEADER))
                continue
            raw_barcode = fields[BARCODE_COLUMN]
            if raw_barcode == NOSEQ:
                stats["no_barcode"] += 1
                corrected = NOSEQ
            else:
                lookup = raw_barcode[:-1] if len(raw_barcode) == barcode_length + 1 else raw_barcode
                corrected = assignments.get(lookup, NOSEQ)
                if corrected == NOSEQ:
                    stats["dropped"] += 1
                elif corrected == lookup:
                    stats["unchanged"] += 1
                else:
                    stats["corrected"] += 1
            fields[BARCODE_COLUMN] = corrected
            outf.write("%s\t%s\n" % ("\t".join(fields), raw_barcode))
    return stats


def load_whitelist(barcode_whitelist: Optional[List[str]]) -> Optional[Set[str]]:
    if not barcode_whitelist:
        return None
    whitelist: Set[str] = set()
    for file_name in barcode_whitelist:
        whitelist.update(load_barcodes(file_name))
    logger.info("Loaded %d whitelisted barcodes" % len(whitelist))
    return whitelist


def correct_barcodes(raw_barcode_files: List[str], output_files: List[str],
                     barcode_whitelist: Optional[List[str]],
                     barcode_length: int = 16,
                     n_cells: Optional[int] = None,
                     n_cells_interval: int = 25,
                     threshold: int = 1,
                     rounds: int = 2,
                     kmer_size: Optional[int] = None,
                     threads: int = 1,
                     stats_file: Optional[str] = None,
                     implementation: str = "centers") -> Dict[str, int]:
    """Run the full correction stage over a sample's raw barcode tables."""
    logger.info("Correcting barcodes using the edit-distance graph")
    graph = BarcodeGraph(threshold=threshold, barcode_length=barcode_length, kmer_size=kmer_size)

    count_barcodes(raw_barcode_files, graph)
    if not graph.counts:
        logger.warning("No barcodes were extracted, nothing to correct")

    whitelist = load_whitelist(barcode_whitelist)
    if implementation == "full":
        # materialise the whole graph; same result, kept for cross-checking
        graph.construct_graph(threads=threads)
        graph.select_cluster_centers(n_cells=n_cells, whitelist=whitelist, interval=n_cells_interval)
        graph.cluster(rounds=rounds)
    else:
        graph.select_cluster_centers(n_cells=n_cells, whitelist=whitelist, interval=n_cells_interval)
        graph.cluster_from_centers(rounds=rounds)

    assignments = graph.get_assignments()
    stats = graph.assignment_stats()
    for raw_file, output_file in zip(raw_barcode_files, output_files):
        file_stats = write_corrected_table(raw_file, output_file, assignments, barcode_length)
        for key, value in file_stats.items():
            stats["Reads %s" % key] = stats.get("Reads %s" % key, 0) + value
        logger.info("Wrote corrected barcodes to %s" % output_file)

    for key, value in stats.items():
        logger.info("  %s: %d" % (key, value))
    if stats_file:
        with open(stats_file, "w") as out_stats:
            for key, value in stats.items():
                out_stats.write("%s\t%d\n" % (key, value))
    return stats
