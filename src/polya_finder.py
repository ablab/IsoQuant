############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from functools import partial

logger = logging.getLogger('IsoQuant')


class PolyAFinder:

    # == polyA stuff ==
    def find_polya_tail(self, alignment):
        cigar_tuples = alignment.cigartuples
        clipped_size = 0
        # hard clipped
        if len(cigar_tuples) > 1 and cigar_tuples[-1][0] == 5 and cigar_tuples[-2][0] == 4:
            clipped_size = cigar_tuples[-2][1]
        elif cigar_tuples[-1][0] == 4:
            clipped_size = cigar_tuples[-1][1]

        pos = self.find_polya(alignment.seq[-clipped_size:].upper())
        if pos == -1:
            return -1
        return alignment.get_blocks()[-1][1] + pos + 1

    def find_polyt_head(self, alignment):
        cigar_tuples = alignment.cigartuples
        clipped_size = 0
        # hard clipped
        if len(cigar_tuples) > 1 and cigar_tuples[0][0] == 5 and cigar_tuples[1][0] == 4:
            clipped_size = cigar_tuples[1][1]
        elif cigar_tuples[0][0] == 4:
            clipped_size = cigar_tuples[0][1]

        pos = self.find_polya(str(Seq(alignment.seq[:clipped_size]).reverse_complement()).upper())
        if pos == -1:
            return -1
        return alignment.reference_start - pos - 1

    # poly A tail detection
    def find_polya(self, seq):
        window_size = 20
        polyA_fraction = 0.9
        polyA_count = int(window_size * polyA_fraction)
        i = 0
        while i < len(seq) - window_size - 1:
            if seq[i:i + window_size].count('A') >= polyA_count:
                break
            i += 1
        if i == len(seq) - window_size - 1:
            return -1
        return i

    def resolve_ambiguity_with_polya(self, read_id, combined_read_profile, matched_isoforms):
        strands = set()
        for isoform_id in matched_isoforms:
            strands.add(self.gene_info.isoform_strands[isoform_id])
        if len(strands) > 1:
            return matched_isoforms

        read_intron_profile = combined_read_profile.read_intron_profile
        read_split_exon_profile = combined_read_profile.read_split_exon_profile
        if list(strands)[0] == '+':
            pos = self.find_polya_tail(combined_read_profile.alignment)
            if pos == -1:
                return matched_isoforms
            elif pos <= combined_read_profile.alignment.get_blocks()[-1][1]:
                logger.warning("Wrong polyA tail position " + str(pos) + ", "
                               "while read ends at " + str(combined_read_profile.alignment.get_blocks()[-1][1]))

            logger.debug("PolyA site found " + str(pos))
            new_split_exon_profile = read_split_exon_profile.gene_profile
            for i in range(len(self.gene_info.split_exon_profiles.features)):
                if self.gene_info.split_exon_profiles.features[i][1] > pos and new_split_exon_profile[i] == 0:
                    new_split_exon_profile[i] = -1
            new_intron_profile = read_intron_profile.gene_profile
            for i in range(len(self.gene_info.intron_profiles.features)):
                if self.gene_info.intron_profiles.features[i][1] > pos and new_intron_profile[i] == 0:
                    new_intron_profile[i] = -1
        else:
            pos = self.find_polyt_head(combined_read_profile.alignment)
            if pos == -1:
                return matched_isoforms
            elif pos >= combined_read_profile.alignment.get_blocks()[0][0]:
                logger.warning("Wrong polyT head position " + str(pos) + ", "
                               "while read starts at " + str(combined_read_profile.alignment.get_blocks()[0][0]))

            logger.debug("PolyT site found " + str(pos))
            new_split_exon_profile = read_split_exon_profile.gene_profile
            for i in range(len(self.gene_info.split_exon_profiles.features)):
                if self.gene_info.split_exon_profiles.features[i][0] < pos and new_split_exon_profile[i] == 0:
                    new_split_exon_profile[i] = -1
            new_intron_profile = read_intron_profile.gene_profile
            for i in range(len(self.gene_info.intron_profiles.features)):
                if self.gene_info.intron_profiles.features[i][0] < pos and new_intron_profile[i] == 0:
                    new_intron_profile[i] = -1

        logger.debug("New exon profile " + str(new_split_exon_profile))
        logger.debug("New intron profile " + str(new_intron_profile))
        return self.match_non_contradictory(read_id, combined_read_profile, use_polya=False).assigned_features
    