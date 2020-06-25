############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from functools import partial
from Bio import Seq

logger = logging.getLogger('IsoQuant')


class PolyAFinder:
    def __init__(self, window_size=20, min_polya_fraction = 0.85):
        self.window_size = window_size
        self.min_polya_fraction = min_polya_fraction

    # == polyA stuff ==
    def find_polya_tail(self, alignment):
        logger.debug("Detecting polyA tail for %s " % alignment.query_name)
        cigar_tuples = alignment.cigartuples
        clipped_size = 0
        # hard clipped

        if len(cigar_tuples) > 1 and cigar_tuples[-1][0] == 5 and cigar_tuples[-2][0] == 4:
            clipped_size = cigar_tuples[-2][1]
        elif cigar_tuples[-1][0] == 4:
            clipped_size = cigar_tuples[-1][1]

        seq = alignment.seq
        if not seq:
            return -1
        tail_len = min(clipped_size+self.window_size, len(seq))
        tail_start = len(seq) - tail_len
        pos = self.find_polya(alignment.seq[tail_start:].upper())
        if pos == -1:
            logger.debug("No polyA found")
            return -1
        # FIXME this does not include indels
        ref_tail_start = alignment.reference_end + clipped_size - tail_len
        ref_polya_start = ref_tail_start + pos + 1
        logger.debug("PolyA found at position %d" % ref_polya_start)
        return ref_polya_start

    def find_polyt_head(self, alignment):
        logger.debug("Detecting polyT head for %s " % alignment.query_name)
        cigar_tuples = alignment.cigartuples
        clipped_size = 0
        # hard clipped
        if len(cigar_tuples) > 1 and cigar_tuples[0][0] == 5 and cigar_tuples[1][0] == 4:
            clipped_size = cigar_tuples[1][1]
        elif cigar_tuples[0][0] == 4:
            clipped_size = cigar_tuples[0][1]

        seq = alignment.seq
        if not seq:
            return -1
        head_len = min(clipped_size+self.window_size, len(seq))
        rc_head = str(Seq.Seq(alignment.seq[:head_len]).reverse_complement()).upper()
        pos = self.find_polya(rc_head)
        if pos == -1:
            logger.debug("No polyT found")
            return -1
        # FIXME this does not include indels
        ref_head_end = alignment.reference_start - clipped_size + head_len
        ref_polyt_end = ref_head_end - pos
        logger.debug("PolyA found at position %d" % ref_polyt_end)
        return ref_polyt_end

    # poly A tail detection
    def find_polya(self, seq):
        polyA_count = int(self.window_size * self.min_polya_fraction)
        i = 0
        while i < len(seq) - self.window_size:
            if seq[i:i + self.window_size].count('A') >= polyA_count:
                break
            i += 1
        if i == len(seq) - self.window_size:
            return -1
        return i
