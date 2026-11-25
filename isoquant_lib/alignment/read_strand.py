############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Transcript strand of a read before it is assigned.

A polyA site that `--fl_data` or `--polya_trimmed` adds sits at the 3' end of the read, so its side needs
the read strand. Predicting the strand before assignment lets the site go in before the read profiles are
built, as for a tail found in the read sequence. Reads left undecided get the site after assignment
(see `AlignmentCollector.process_genic`).
"""

from bisect import bisect_left, bisect_right
from enum import Enum, unique
from typing import Dict, List, Optional, Set, Tuple

from isoquant_lib.common import get_intron_strand, junctions_from_blocks
from isoquant_lib.terminal_prediction.external_polya import has_polya_on_strand, hint_from_strandedness
from isoquant_lib.terminal_prediction.polya_finder import PolyAInfo


@unique
class StrandSource(Enum):
    undecided = 0
    polya_tail = 1
    splice_sites = 2
    transcript_overlap = 3
    stranded_library = 4


class ReadStrandPredictor:
    """Signals, in order of priority:
    1. a tail found at one read end;
    2. spliced reads: strands of introns matching annotated introns exactly, and canonical splice sites of the
       other introns; all of them must agree;
    3. unspliced reads: strands of the annotated transcripts overlapping the read; all of them must agree;
    4. --stranded forward/reverse: the mapped strand or the opposite one (library strandness is rarely complete,
       so it is used only when the read itself says nothing)."""

    def __init__(self, gene_info, chr_record, stranded: Optional[str] = None):
        self.chr_record = chr_record
        # per region: reads share their introns, and chr_record is usually read from disk (pyfaidx)
        self.splice_site_strands: Dict[Tuple[int, int], str] = {}
        self.library_hint: Optional[str] = hint_from_strandedness(stranded)
        self.intron_strands: Dict[Tuple[int, int], Set[str]] = {}
        for isoform_id, introns in gene_info.all_isoforms_introns.items():
            strand = gene_info.isoform_strands[isoform_id]
            if strand not in ("+", "-"):
                continue
            for intron in introns:
                self.intron_strands.setdefault(intron, set()).add(strand)

        # transcript spans rather than gene records: a gene may hold transcripts on both strands
        transcripts = sorted((exons[0][0], exons[-1][1], gene_info.isoform_strands[isoform_id])
                             for isoform_id, exons in gene_info.all_isoforms_exons.items()
                             if exons and gene_info.isoform_strands[isoform_id] in ("+", "-"))
        self.transcript_starts: List[int] = [t[0] for t in transcripts]
        self.transcripts: List[Tuple[int, int, str]] = transcripts
        self.max_transcript_len: int = max((t[1] - t[0] for t in transcripts), default=0)

    def predict(self, read_exons: List[Tuple[int, int]], polya_info: PolyAInfo,
                mapped_strand: str) -> Tuple[Optional[str], StrandSource]:
        has_polya = has_polya_on_strand(polya_info, "+")
        has_polyt = has_polya_on_strand(polya_info, "-")
        if has_polya != has_polyt:
            return ("+" if has_polya else "-"), StrandSource.polya_tail

        introns = junctions_from_blocks(read_exons)
        if introns:
            strand = self.splice_strand(introns)
            if strand is not None:
                return strand, StrandSource.splice_sites
        else:
            strand = self.overlap_strand(read_exons[0][0], read_exons[-1][1])
            if strand is not None:
                return strand, StrandSource.transcript_overlap

        if self.library_hint is not None:
            if self.library_hint == "+":
                return mapped_strand, StrandSource.stranded_library
            return ("-" if mapped_strand == "+" else "+"), StrandSource.stranded_library
        return None, StrandSource.undecided

    def splice_strand(self, introns: List[Tuple[int, int]]) -> Optional[str]:
        # exact match only: a shifted annotated intron is checked by its splice sites
        strands = set()
        for intron in introns:
            annotated = self.intron_strands.get(intron)
            if annotated:
                strands.update(annotated)
            elif self.chr_record is not None:
                strand = self.splice_site_strands.get(intron)
                if strand is None:
                    strand = get_intron_strand(intron, self.chr_record)
                    self.splice_site_strands[intron] = strand
                if strand != ".":
                    strands.add(strand)
            if len(strands) > 1:
                return None
        return next(iter(strands)) if strands else None

    def overlap_strand(self, start: int, end: int) -> Optional[str]:
        strands = set()
        first = bisect_left(self.transcript_starts, start - self.max_transcript_len)
        last = bisect_right(self.transcript_starts, end)
        for _, transcript_end, strand in self.transcripts[first:last]:
            if transcript_end >= start:
                strands.add(strand)
                if len(strands) > 1:
                    return None
        return next(iter(strands)) if strands else None
