############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging

import pybedtools as pbt
from Bio import Seq
from common import move_ref_coord_along_alignment, CigarEvent


logger = logging.getLogger('IsoQuant')


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
        return self.find_polya_tail(alignment, 6, 2 * self.window_size)

    def find_polya_internal(self, alignment):
        return self.find_polya_tail(alignment, 4 * self.window_size, 2, check_entire_tail=True)

    def find_polyt_external(self, alignment):
        return self.find_polyt_head(alignment, 6, 2 * self.window_size)

    def find_polyt_internal(self, alignment):
        return self.find_polyt_head(alignment, 4 * self.window_size, 2, check_entire_head=True)

    # == polyA stuff ==
    def find_polya_tail(self, alignment, from_pos, to_pos, check_entire_tail=False):
        # logger.debug("Detecting polyA tail for %s " % alignment.query_name)
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
        pos = self.find_polya(sequence_to_check)

        #logger.debug("read start: %d, ckeck start: %d, check end: %d, pos: %d" % (read_mapped_region_end, to_check_start, to_check_end, pos))
        #logger.debug(sequence_to_check)
        if check_entire_tail and pos != -1:
            entire_tail = sequence_to_check[pos:]
            if entire_tail.count('A') < len(entire_tail) * self.min_polya_fraction:
                logger.debug("Internal polyA seems unreliable")
                pos = -1

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
            ref_shift = move_ref_coord_along_alignment(alignment, shift)
            reference_polya_start = alignment.reference_end - ref_shift
            #logger.debug("shift: %d, ref shift: %d, reference: %d" % (shift, ref_shift, reference_polya_start))

        logger.debug("PolyA found at position %d" % reference_polya_start)
        return reference_polya_start

    def find_polyt_head(self, alignment, from_pos, to_pos, check_entire_head=False):
        # logger.debug("Detecting polyT head for %s " % alignment.query_name)
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

        pos = self.find_polya(sequence_to_check)
        #logger.debug("read start: %d, ckeck start: %d, check end: %d, pos: %d" % (read_mapped_region_start, to_check_start, to_check_end, pos))
        #logger.debug(sequence_to_check)

        if check_entire_head and pos != -1:
            entire_tail = sequence_to_check[pos:]
            if entire_tail.count('A') < len(entire_tail) * self.min_polya_fraction:
                #logger.debug("Internal polyT seems unreliable")
                pos = -1

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
            ref_shift = move_ref_coord_along_alignment(alignment, shift)
            reference_polyt_end = alignment.reference_start + ref_shift
            #logger.debug("shift: %d, ref shift: %d, reference: %d" % (shift, ref_shift, reference_polyt_end))

        logger.debug("PolyT found at position %d" % reference_polyt_end)
        return max(1, reference_polyt_end)

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

    def check_aligned_bases_near_polya(self, alignment, chr_str, reference_polya_start, bases_to_check=10):
        total_errors = 0
        remaining_bases = bases_to_check
        cigar_tuples = alignment.cigartuples
        cigar_pos = len(cigar_tuples) - 1
        read_pos = len(alignment.seq) - 1
        ref_pos = alignment.reference_end

        if alignment.reference_end < reference_polya_start:
            dist = reference_polya_start - alignment.reference_end
            total_errors += dist
            remaining_bases -= dist

        skip_bases = 0
        if reference_polya_start < alignment.reference_end:
            skip_bases = alignment.reference_end - reference_polya_start



        if CigarEvent(cigar_tuples[cigar_pos][0]) == CigarEvent.hard_clipping:
            cigar_pos -= 1
        if CigarEvent(cigar_tuples[cigar_pos][0]) == CigarEvent.soft_clipping:
            cigar_pos -= 1
            read_pos -= cigar_tuples[cigar_pos][1]

        remaining_cigar_event_len = cigar_tuples[cigar_pos][1]
        current_event = CigarEvent(cigar_tuples[cigar_pos][0])
        while remaining_bases > 0:
            if remaining_cigar_event_len == 0:
                if cigar_pos == 0:
                    break
                cigar_pos -= 1
                remaining_cigar_event_len = cigar_tuples[cigar_pos][1]
                current_event = CigarEvent(cigar_tuples[cigar_pos][0])

            if current_event in CigarEvent.get_match_events():
                if skip_bases == 0:
                    if current_event == CigarEvent.seq_mismatch:
                        total_errors += 1
                    elif current_event == CigarEvent.match and alignment.seq[read_pos] != chr_str[ref_pos]:
                        total_errors += 1
                ref_pos -= 1
                read_pos -= 1
            elif current_event == CigarEvent.deletion:
                if skip_bases == 0:
                    total_errors += 1
                ref_pos -= 1
            elif current_event == CigarEvent.insertion:
                if skip_bases == 0:
                    total_errors += 1
                read_pos -= 1
            else:
                break

            if skip_bases == 0:
                remaining_bases -= 1
            else:
                skip_bases -= 1
            remaining_cigar_event_len -= 1

        if remaining_bases > 0:
            total_errors += remaining_bases

    def concat_gapless_blocks(blocks, cigar_tuples):
        cigar_index = 0
        block_index = 0
        resulting_blocks = []

        current_block = None
        deletions_before_block = 0

        while cigar_index < len(cigar_tuples) and block_index < len(blocks):
            # init new block
            cigar_event = CigarEvent(cigar_tuples[cigar_index][0])
            if current_block is None:
                # init new block from match
                if cigar_event in CigarEvent.get_match_events():
                    current_block = (blocks[block_index][0] - deletions_before_block, blocks[block_index][1])
                    deletions_before_block = 0
                    block_index += 1
                # keep track of deletions before matched block
                elif cigar_event == CigarEvent.deletion:
                    deletions_before_block = cigar_tuples[cigar_index][1]
            # found intron, add current block
            elif cigar_event == CigarEvent.skipped:
                resulting_blocks.append(current_block)
                current_block = None
            # add deletion to block
            elif cigar_event == CigarEvent.deletion:
                current_block = (current_block[0], current_block[1] + cigar_tuples[cigar_index][1])
            # found match - merge blocks
            elif cigar_event in CigarEvent.get_match_events():
                # if abs(current_block[1] - blocks[block_index][0]) > 1:
                #    logger.debug("Distant blocks")
                #    logger.debug(current_block, blocks[block_index])
                current_block = (current_block[0], blocks[block_index][1])

                block_index += 1
            cigar_index += 1

        if current_block is not None:
            resulting_blocks.append(current_block)

        return resulting_blocks
