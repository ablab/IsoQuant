############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from functools import partial

import pybedtools as pbt
from Bio import Seq

logger = logging.getLogger('IsoQuant')


# given a cigar string of an alignment, move shift bases along the alignment from the start
# if shift is negative, moves from the alignment's end
# return value is always positive
def move_ref_coord_alogn_alignment(alignment, shift):
    if shift == 0:
        return 0

    cigar_tuples = alignment.cigartuples
    assert cigar_tuples

    direction = 1 if shift > 0 else -1 # 1 for forward -1 for backward
    shift = abs(shift)

    if direction == 1:
        current_pos = 0
        if len(cigar_tuples) > 1 and cigar_tuples[0][0] == 5 and cigar_tuples[1][0] == 4:
            # hard clipped
            current_pos = 2
        elif cigar_tuples[0][0] in {4, 5}:
            # soft clipped
            current_pos = 1
    else:
        current_pos = -1
        if len(cigar_tuples) > 1 and cigar_tuples[-1][0] == 5 and cigar_tuples[-2][0] == 4:
            # hard clipped
            current_pos -= 2
        elif cigar_tuples[-1][0] in {4, 5}:
            # soft clipped
            current_pos -= 1

    read_length_consumed = 0
    reference_length_consumed = 0
    while current_pos < len(cigar_tuples) and current_pos >= -len(cigar_tuples) and read_length_consumed < shift:
        cigar_event = cigar_tuples[current_pos][0]
        event_len = cigar_tuples[current_pos][1]
        if cigar_event == 1:
            # insertion
            read_length_consumed += event_len
        elif cigar_event in {2, 3}:
            # deletion or intron
            reference_length_consumed += event_len
        elif cigar_event == 0:
            # match
            remaining_bases = shift - read_length_consumed # how many nucleotides in read to complete shift
            #logger.debug("%d" % remaining_bases)
            if event_len < remaining_bases:
                reference_length_consumed += event_len
                read_length_consumed += event_len
            else:
                read_length_consumed += remaining_bases
                reference_length_consumed += remaining_bases
        elif cigar_event in {4, 5}:
            # met clipping on the other side
            break
        else:
            # unexpected event
            logger.warning("Unexpected event: " + cigar_event)
        #logger.debug("%d, %d, %d, %d" % (cigar_event, event_len, read_length_consumed, reference_length_consumed))

        current_pos += direction

    return reference_length_consumed


class PolyAInfo:
    def __init__(self, external_polya_pos, external_polyt_pos, internal_polya_pos, internal_polyt_pos):
        self.external_polya_pos = external_polya_pos
        self.external_polyt_pos = external_polyt_pos
        self.internal_polya_pos = internal_polya_pos
        self.internal_polyt_pos = internal_polyt_pos


class PolyAFinder:
    def __init__(self, window_size=16, min_polya_fraction=0.75):
        self.window_size = window_size
        self.min_polya_fraction = min_polya_fraction
        self.polyA_count = int(self.window_size * self.min_polya_fraction)

    def detect_polya(self, alignement):
        return PolyAInfo(self.find_polya_external(alignement), self.find_polyt_external(alignement),
                         self.find_polya_internal(alignement), self.find_polyt_internal(alignement))

    def find_polya_external(self, alignment):
        return self.find_polya_tail(alignment, 2, 2 * self.window_size)

    def find_polya_internal(self, alignment):
        return self.find_polya_tail(alignment, 3 * self.window_size, 2)

    def find_polyt_external(self, alignment):
        return self.find_polyt_head(alignment, 2, 2 * self.window_size)

    def find_polyt_internal(self, alignment):
        return self.find_polyt_head(alignment, 3 * self.window_size, 2)

    # == polyA stuff ==
    def find_polya_tail(self, alignment, from_pos, to_pos):
        logger.debug("Detecting polyA tail for %s " % alignment.query_name)
        cigar_tuples = alignment.cigartuples
        soft_clipped_tail_len = 0

        if len(cigar_tuples) > 1 and cigar_tuples[-1][0] == 5 and cigar_tuples[-2][0] == 4:
            # hard clipped
            soft_clipped_tail_len = cigar_tuples[-2][1]
        elif cigar_tuples[-1][0] == 4:
            # soft clipped
            soft_clipped_tail_len = cigar_tuples[-1][1]

        seq = alignment.seq
        if not seq:
            return -1
        assert soft_clipped_tail_len < len(seq)

        read_mapped_region_end = len(seq) - soft_clipped_tail_len
        to_check_start = max(0, read_mapped_region_end - from_pos)
        to_check_end = min(len(seq), read_mapped_region_end + to_pos + 1)
        sequence_to_check = alignment.seq[to_check_start:to_check_end].upper()
        pos = self.find_polya(sequence_to_check.upper())

        logger.debug("read start: %d, ckeck start: %d, check end: %d, pos: %d" % (read_mapped_region_end, to_check_start, to_check_end, pos))
        logger.debug(sequence_to_check)

        if pos == -1:
            logger.debug("No polyA found")
            return -1

        # add position of the region we check
        pos = to_check_start + pos
        if pos >= read_mapped_region_end: # poly A starts after last mapped base
            shift = pos - read_mapped_region_end
            reference_polya_start = alignment.reference_end + shift
        else:
            shift = pos - read_mapped_region_end
            ref_shift = move_ref_coord_alogn_alignment(alignment, shift)
            reference_polya_start = alignment.reference_end - ref_shift
            logger.debug("shift: %d, ref shift: %d, reference: %d" % (shift, ref_shift, reference_polya_start))

        logger.debug("PolyA found at position %d" % reference_polya_start)
        return reference_polya_start

    def find_polyt_head(self, alignment, from_pos, to_pos):
        logger.debug("Detecting polyT head for %s " % alignment.query_name)
        cigar_tuples = alignment.cigartuples
        soft_clipped_head_len = 0

        if len(cigar_tuples) > 1 and cigar_tuples[0][0] == 5 and cigar_tuples[1][0] == 4:
            # hard clipped
            soft_clipped_head_len = cigar_tuples[1][1]
        elif cigar_tuples[0][0] == 4:
            # soft clipped
            soft_clipped_head_len = cigar_tuples[0][1]

        seq = alignment.seq
        if not seq:
            return -1
        assert soft_clipped_head_len < len(seq)

        read_mapped_region_start = soft_clipped_head_len
        to_check_start = max(0, read_mapped_region_start - to_pos)
        to_check_end = min(len(seq), read_mapped_region_start + from_pos + 1)
        sequence_to_check = str(Seq.Seq(alignment.seq[to_check_start:to_check_end]).reverse_complement()).upper()

        pos = self.find_polya(sequence_to_check.upper())

        logger.debug("read start: %d, ckeck start: %d, check end: %d, pos: %d" % (read_mapped_region_start, to_check_start, to_check_end, pos))
        logger.debug(sequence_to_check)

        if pos == -1:
            logger.debug("No polyT found")
            return -1

        # reverse position
        pos = to_check_end - pos - 1
        if pos <= read_mapped_region_start:  # poly T starts to the left of first mapped base
            shift = read_mapped_region_start - pos
            reference_polyt_end = alignment.reference_start - shift
        else:
            shift = pos - read_mapped_region_start
            ref_shift = move_ref_coord_alogn_alignment(alignment, shift)
            reference_polyt_end = alignment.reference_start + ref_shift
            logger.debug("shift: %d, ref shift: %d, reference: %d" % (shift, ref_shift, reference_polyt_end))

        logger.debug("PolyT found at position %d" % reference_polyt_end)
        return reference_polyt_end

    # poly A tail detection
    def find_polya(self, seq):
        if len(seq) < self.window_size:
            return -1
        i = 0
        a_count = seq[0:self.window_size].count('A')
        while i < len(seq) - self.window_size:
            if a_count >= self.polyA_count:
                break
            first_base_a = seq[i] == 'A'
            new_base_a = i + self.window_size < len(seq) and seq[i + self.window_size] == 'A'
            if first_base_a and not new_base_a:
                a_count -= 1
            elif not first_base_a and new_base_a:
                a_count += 1
            i += 1

        if i >= len(seq) - self.window_size:
            return -1

        return i + max(0, seq[i:].find('AA'))


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
