#!/usr/bin/env python3

############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Convert read_info.tsv to legacy read_assignments.tsv or allinfo formats.

Usage:
    python -m isoquant_lib.convert_read_info \
        --read_info SAMPLE.read_info.tsv.gz \
        --format read_assignments \
        --output SAMPLE.read_assignments.tsv

    python -m isoquant_lib.convert_read_info \
        --read_info SAMPLE.read_info.tsv.gz \
        --format allinfo \
        --output SAMPLE.allinfo.tsv
"""

import argparse
import gzip
import sys
from typing import TextIO

from .common import junctions_from_blocks


# read_info column indices
RI_READ_ID = 0
RI_CHR = 1
RI_STRAND = 2
RI_GENE_ID = 3
RI_GENE_ASSIGNMENT_TYPE = 4
RI_ISOFORM_ID = 5
RI_ISOFORM_ASSIGNMENT_TYPE = 6
RI_ASSIGNMENT_EVENTS = 7
RI_CLASSIFICATION = 8
RI_EXONS = 9
RI_POLYA = 10
RI_CAGE = 11
RI_CANONICAL = 12
RI_BARCODE = 13
RI_UMI = 14
RI_CELL_TYPE = 15
RI_GROUPS = 16
RI_ADDITIONAL = 17


def _open_file(path: str, mode: str = "r") -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def _parse_exons(exons_str: str) -> list:
    """Parse '100-200,300-400' into [(100,200), (300,400)]."""
    if exons_str == ".":
        return []
    pairs = []
    for part in exons_str.split(","):
        s, e = part.split("-")
        pairs.append((int(s), int(e)))
    return pairs


def _exons_to_range_str(exons: list) -> str:
    """Convert [(100,200), (300,400)] to '100-200,300-400'."""
    if not exons:
        return "."
    return ",".join(f"{s}-{e}" for s, e in exons)


def _assignment_type_to_read_type(assignment_type: str) -> str:
    """Map isoform_assignment_type to allinfo read_type."""
    if assignment_type in ("unique", "unique_minor_difference"):
        return "known"
    elif assignment_type in ("ambiguous",):
        return "known_ambiguous"
    elif assignment_type.startswith("inconsistent"):
        return "novel"
    else:
        return "none"


def _format_coord(chr_id: str, start: int, end: int, strand: str) -> str:
    return f"{chr_id}_{start}_{end}_{strand}"


def convert_to_read_assignments(read_info_path: str, output_path: str):
    """Convert read_info.tsv to legacy read_assignments.tsv format."""
    with _open_file(read_info_path) as inf, _open_file(output_path, "w") as outf:
        # Copy comment lines from read_info header
        for line in inf:
            if line.startswith("#"):
                outf.write(line)
            else:
                break  # first non-comment line is the header, skip it

        # Write read_assignments header
        outf.write("read_id\tchr\tstrand\tisoform_id\tgene_id"
                    "\tassignment_type\tassignment_events\texons\tadditional_info\tgroups\n")

        # Process the first data line (already read above)
        if not line.startswith("#"):
            _convert_ra_line(line, outf)

        for line in inf:
            if line.startswith("#"):
                continue
            _convert_ra_line(line, outf)


def _convert_ra_line(line: str, outf: TextIO):
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 18:
        return

    # Reconstruct additional_info blob
    additional_parts = []
    additional_parts.append(f"gene_assignment={cols[RI_GENE_ASSIGNMENT_TYPE]};")
    additional_parts.append(f"PolyA={cols[RI_POLYA]};")
    if cols[RI_CAGE] != ".":
        additional_parts.append(f"CAGE={cols[RI_CAGE]};")
    if cols[RI_CANONICAL] != ".":
        additional_parts.append(f"Canonical={cols[RI_CANONICAL]};")
    additional_parts.append(f"Classification={cols[RI_CLASSIFICATION]};")
    # Append remaining additional attributes
    if cols[RI_ADDITIONAL] != "*":
        additional_parts.append(cols[RI_ADDITIONAL])
    additional_info = " ".join(additional_parts) if additional_parts else "*"

    groups = cols[RI_GROUPS] if cols[RI_GROUPS] != "." else ""

    out_cols = [
        cols[RI_READ_ID],
        cols[RI_CHR],
        cols[RI_STRAND],
        cols[RI_ISOFORM_ID],
        cols[RI_GENE_ID],
        cols[RI_ISOFORM_ASSIGNMENT_TYPE],
        cols[RI_ASSIGNMENT_EVENTS],
        cols[RI_EXONS],
        additional_info,
        groups,
    ]
    outf.write("\t".join(out_cols) + "\n")


def convert_to_allinfo(read_info_path: str, output_path: str):
    """Convert read_info.tsv to allinfo format (SC/spatial)."""
    with _open_file(read_info_path) as inf, _open_file(output_path, "w") as outf:
        for line in inf:
            if line.startswith("#") or line.startswith("read_id"):
                continue
            _convert_allinfo_line(line, outf)


def _convert_allinfo_line(line: str, outf: TextIO):
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 18:
        return

    chr_id = cols[RI_CHR]
    strand = cols[RI_STRAND]
    exons = _parse_exons(cols[RI_EXONS])
    introns = junctions_from_blocks(exons) if exons else []

    exons_str = ";%;" + ";%;".join(
        _format_coord(chr_id, e[0], e[1], strand) for e in exons) if exons else ""
    introns_str = ";%;" + ";%;".join(
        _format_coord(chr_id, i[0], i[1], strand) for i in introns) if introns else ""

    # Reconstruct TSS/polyA from events
    events = cols[RI_ASSIGNMENT_EVENTS]
    TSS = "NoTSS"
    polyA = "NoPolyA"
    if exons:
        if "terminal_site_match_left" in events or "tss_match" in events:
            tss_pos = exons[0][0] if strand == "+" else exons[-1][1]
            TSS = _format_coord(chr_id, tss_pos, tss_pos, strand)
        if "correct_polya_site" in events:
            polya_pos = exons[-1][1] if strand == "+" else exons[0][0]
            polyA = _format_coord(chr_id, polya_pos, polya_pos, strand)

    read_type = _assignment_type_to_read_type(cols[RI_ISOFORM_ASSIGNMENT_TYPE])

    # Get transcript_type from additional
    transcript_type = "unknown"
    if cols[RI_ADDITIONAL] != "*":
        for part in cols[RI_ADDITIONAL].split():
            if part.startswith("transcript_type="):
                transcript_type = part.split("=", 1)[1].rstrip(";")
                break

    out_cols = [
        cols[RI_READ_ID],
        cols[RI_GENE_ID],
        cols[RI_CELL_TYPE],
        cols[RI_BARCODE],
        cols[RI_UMI],
        introns_str,
        TSS,
        polyA,
        exons_str,
        read_type,
        str(len(introns)),
        cols[RI_ISOFORM_ID],
        transcript_type,
    ]
    outf.write("\t".join(out_cols) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Convert read_info.tsv to legacy formats")
    parser.add_argument("--read_info", required=True, help="Input read_info.tsv[.gz] file")
    parser.add_argument("--format", required=True, choices=["read_assignments", "allinfo"],
                        help="Output format")
    parser.add_argument("--output", required=True, help="Output file path")
    args = parser.parse_args()

    if args.format == "read_assignments":
        convert_to_read_assignments(args.read_info, args.output)
    elif args.format == "allinfo":
        convert_to_allinfo(args.read_info, args.output)


if __name__ == "__main__":
    main()
