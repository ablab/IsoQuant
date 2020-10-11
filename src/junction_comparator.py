############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import namedtuple

import pysam
from Bio import pairwise2

from src.isoform_assignment import *
from src.gene_info import *
from src.long_read_profiles import *


logger = logging.getLogger('IsoQuant')
pairwise2.MAX_ALIGNMENTS = 1


class JunctionComparator:
    absent = -10

    def __init__(self, params, intron_profile_constructor):
        self.params = params
        self.intron_profile_constructor = intron_profile_constructor
        self.misalignment_checker = MisalignmentChecker(params.reference, params.delta, params.max_intron_shift,
                                                        params.max_missed_exon_len)
        # self.exon_misalignment_checker = ExonMisalignmentChecker(params.reference)

    def compare_junctions(self, read_junctions, read_region, isoform_junctions, isoform_region, alignment=None):
        """ compare read splice junctions against similar isoform

        Parameters
        ----------
        read_junctions: list of tuple of int
        read_region: tuple of int
        isoform_junctions: list of tuple of int
        isoform_region: tuple of int

        Returns
        -------
        list of detected contradiction events
        """
        if not read_junctions:
            return [self.get_mono_exon_subtype(read_region, isoform_junctions)]

        read_pos = 0
        isoform_pos = 0
        read_features_present = [0 for i in range(0, len(read_junctions))]
        isoform_features_present = [0 for i in range(0, len(isoform_junctions))]
        contradictory_region_pairs = []
        current_contradictory_region = (self.absent, self.absent)

        while read_pos < len(read_junctions) and isoform_pos < len(isoform_junctions):
            if equal_ranges(isoform_junctions[isoform_pos], read_junctions[read_pos], self.params.delta):
                # junctions are equal
                read_features_present[read_pos] = 1
                isoform_features_present[isoform_pos] = 1
                if current_contradictory_region != (self.absent, self.absent):
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = (self.absent, self.absent)
                read_pos += 1
                isoform_pos += 1

            elif overlaps(isoform_junctions[isoform_pos], read_junctions[read_pos]):
                # junctions overlap, but are unequal
                read_features_present[read_pos] = -1
                isoform_features_present[isoform_pos] = -1
                if current_contradictory_region == (self.absent, self.absent):
                    current_contradictory_region = ((read_pos, read_pos), (isoform_pos, isoform_pos))
                else:
                    current_contradictory_region = (
                    (current_contradictory_region[0][0], read_pos), (current_contradictory_region[1][0], isoform_pos))
                if read_junctions[read_pos][1] < isoform_junctions[isoform_pos][1]:
                    read_pos += 1
                else:
                    isoform_pos += 1

            elif left_of(isoform_junctions[isoform_pos], read_junctions[read_pos]):
                # isoform junction is behind, move on
                if current_contradictory_region != (self.absent, self.absent):
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = (self.absent, self.absent)
                if read_pos > 0 or overlaps(read_region, isoform_junctions[isoform_pos]):
                    if isoform_features_present[isoform_pos] != -1:
                        contradictory_region_pairs.append(((self.absent, self.absent), (isoform_pos, isoform_pos)))
                    isoform_features_present[isoform_pos] = -1
                isoform_pos += 1

            else:
                # read junction is behind, move on
                if current_contradictory_region != (self.absent, self.absent):
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = (self.absent, self.absent)
                if isoform_pos > 0 or overlaps(isoform_region, read_junctions[read_pos]):
                    if read_features_present[read_pos] != -1:
                        contradictory_region_pairs.append(((read_pos, read_pos), (self.absent, isoform_pos)))
                    read_features_present[read_pos] = -1
                read_pos += 1

        if current_contradictory_region != (self.absent, self.absent):
            contradictory_region_pairs.append(current_contradictory_region)

        # check terminating regions
        while read_pos < len(read_junctions):
            if overlaps(isoform_region, read_junctions[read_pos]):
                if read_features_present[read_pos] != -1:
                    contradictory_region_pairs.append(((read_pos, read_pos), (self.absent, isoform_pos)))
                    read_features_present[read_pos] = -1
            else:
                break
            read_pos += 1

        while isoform_pos < len(isoform_junctions):
            if overlaps(read_region, isoform_junctions[isoform_pos]):
                if isoform_features_present[isoform_pos] != -1:
                    contradictory_region_pairs.append(((self.absent, self.absent), (isoform_pos, isoform_pos)))
                    isoform_features_present[isoform_pos] = -1
            else:
                break
            isoform_pos += 1

        logger.debug("+ + Inspected contradictory read")
        # logger.debug("+ + Read profile " + str(read_features_present))
        # logger.debug("+ + Read introns " + str(read_junctions))
        # logger.debug("+ + Read region " + str(read_region))

        # logger.debug("+ + Isoform profile " + str(isoform_features_present))
        # logger.debug("+ + Isoform introns " + str(isoform_junctions))
        # logger.debug("+ + Isoform region " + str(isoform_region))

        matching_events = []
        if any(el == -1 for el in read_features_present) or any(el == -1 for el in isoform_features_present):
            # classify contradictions
            logger.debug("+ + Classifying contradictions")
            matching_events = self.detect_contradiction_type(read_junctions, isoform_junctions,
                                                             contradictory_region_pairs, alignment)

        if read_features_present[0] == 0 or read_features_present[-1] == 0:
            if all(x == 0 for x in read_features_present):
                return [make_event(MatchEventSubtype.undefined)]
            logger.debug("+ + Found only extra terminal introns ")
            self.add_extra_out_exon_events(matching_events, read_features_present, read_junctions, isoform_region[0])

        if len(matching_events) == 0:
            logger.debug("No contradiction detected")
            return [make_event(MatchEventSubtype.none)]
        return matching_events

    def detect_contradiction_type(self, read_junctions, isoform_junctions, contradictory_region_pairs, alignment):
        """

        Parameters
        ----------
        read_junctions: list of tuples of int
        isoform_junctions: list of tuples of int
        contradictory_region_pairs: list of tuples of int

        Returns
        -------
        list of contradiction events
        """
        contradiction_events = []
        for read_cregion, isoform_cregion in contradictory_region_pairs:
            # classify each contradictory area separately
            event = self.compare_overlapping_contradictional_regions(read_junctions, isoform_junctions,
                                                                     read_cregion, isoform_cregion, alignment)
            contradiction_events.append(event)

        return contradiction_events

    def compare_overlapping_contradictional_regions(self, read_junctions, isoform_junctions, read_cregion,
                                                    isoform_cregion, alignment=None):
        if read_cregion[0] == self.absent:
            return make_event(MatchEventSubtype.intron_retention, isoform_cregion[0], read_cregion)
        elif isoform_cregion[0] == self.absent:
            # intron_start = read_junctions[read_cregion[0]]
            if self.are_known_introns(read_junctions, read_cregion):
                return make_event(MatchEventSubtype.extra_intron_known, isoform_cregion[1], read_cregion)
            return make_event(MatchEventSubtype.extra_intron, isoform_cregion[1], read_cregion)

        read_intron_total_len = sum(
            [read_junctions[i][1] - read_junctions[i][0] + 1 for i in range(read_cregion[0], read_cregion[1] + 1)])
        isoform_intron_total_len = sum(
            [isoform_junctions[i][1] - isoform_junctions[i][0] + 1 for i in range(isoform_cregion[0], isoform_cregion[1] + 1)])
        total_intron_len_diff = abs(read_intron_total_len - isoform_intron_total_len)

        read_introns_known = self.are_known_introns(read_junctions, read_cregion)

        if read_cregion[1] == read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            event = self.classify_single_intron_alternation(
                read_junctions, isoform_junctions, read_cregion[0], isoform_cregion[0],
                total_intron_len_diff, read_introns_known, alignment)

        elif read_cregion[1] - read_cregion[0] == isoform_cregion[1] - isoform_cregion[0] and \
                total_intron_len_diff < self.params.delta:
            if read_introns_known:
                event = MatchEventSubtype.mutually_exclusive_exons_known
            else:
                event = MatchEventSubtype.mutually_exclusive_exons_novel

        elif read_cregion[1] == read_cregion[0] and isoform_cregion[1] > isoform_cregion[0]:
            event = self.classify_skipped_exons(isoform_junctions, isoform_cregion, read_junctions, read_cregion,
                                                total_intron_len_diff, read_introns_known, alignment=alignment)

        elif read_cregion[1] > read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            if read_introns_known:
                event = MatchEventSubtype.exon_gain_known
            else:
                event = MatchEventSubtype.exon_gain_novel

        else:
            if read_introns_known:
                event = MatchEventSubtype.alternative_structure_known
            else:
                event = MatchEventSubtype.alternative_structure_novel
        return make_event(event, isoform_cregion[0], read_cregion)

    def classify_skipped_exons(self, isoform_junctions, isoform_cregion, read_junctions, read_cregion,
                               total_intron_len_diff, read_introns_known, alignment=None):
        total_exon_len = sum([isoform_junctions[i + 1][0] - isoform_junctions[i][1] + 1
                              for i in range(isoform_cregion[0], isoform_cregion[1])])

        return self.misalignment_checker.classify_skipped_exon(
            isoform_junctions, isoform_cregion, read_junctions, read_cregion,
            total_intron_len_diff, read_introns_known, total_exon_len, alignment)

    def classify_single_intron_alternation(self, read_junctions, isoform_junctions, read_cpos, isoform_cpos,
                                           total_intron_len_diff, read_introns_known, alignment=None):
        if total_intron_len_diff <= 2 * self.params.delta:
            if read_introns_known:
                return MatchEventSubtype.intron_migration
            else:
                return self.misalignment_checker.classify_intron_misalignment(
                    read_junctions, isoform_junctions, read_cpos, isoform_cpos, alignment)
        else:
            # TODO correct when strand is negative
            if abs(isoform_junctions[isoform_cpos][0] - read_junctions[read_cpos][0]) < self.params.delta:
                return MatchEventSubtype.alt_acceptor_site_known if read_introns_known \
                    else MatchEventSubtype.alt_acceptor_site_novel
            elif abs(isoform_junctions[isoform_cpos][1] - read_junctions[read_cpos][1]) < self.params.delta:
                return MatchEventSubtype.alt_donor_site_known if read_introns_known \
                    else MatchEventSubtype.alt_donor_site_novel
            else:
                return MatchEventSubtype.intron_alternation_known if read_introns_known \
                    else MatchEventSubtype.intron_alternation_novel

    def get_mono_exon_subtype(self, read_region, isoform_junctions):
        if len(isoform_junctions) == 0:
            event = MatchEventSubtype.mono_exon_match
        elif not any(overlaps(read_region, rj) for rj in isoform_junctions):
            event = MatchEventSubtype.mono_exonic
        elif any(contains(read_region, rj) for rj in isoform_junctions):
            # TODO save intron retention position
            event = MatchEventSubtype.unspliced_intron_retention
        else:
            event = MatchEventSubtype.unspliced_genic
        return make_event(event)

    def add_extra_out_exon_events(self, match_events, read_intron_read_profile, read_introns, isoform_start):
        extra_left = read_intron_read_profile[0] == 0
        extra_right = read_intron_read_profile[-1] == 0

        if all(x == 0 for x in read_intron_read_profile):
            if read_introns[0][0] < isoform_start:
                extra_right = False
            else:
                extra_left = False

        if extra_left:
            read_pos = 0
            while read_pos < len(read_intron_read_profile) and read_intron_read_profile[read_pos] == 0:
                read_pos += 1
            match_events.append(make_event(MatchEventSubtype.extra_intron_flanking_left,
                                           SupplementaryMatchConstansts.extra_left_mod_position, (0, read_pos - 1)))
        if extra_right:
            max_right = len(read_intron_read_profile) - 1
            read_pos = max_right
            while read_pos >= 0 and read_intron_read_profile[read_pos] == 0:
                read_pos -= 1
            match_events.append(make_event(MatchEventSubtype.extra_intron_flanking_right,
                                           SupplementaryMatchConstansts.extra_right_mod_position, (read_pos + 1, max_right)))

    def profile_for_junctions_introns(self, junctions, region):
        selected_junctions = []
        logger.debug("Checking for known introns " + str(region))
        for i in range(region[0], region[1] + 1):
            selected_junctions.append(junctions[i])

        selected_junctions_profile = self.intron_profile_constructor.construct_profile_for_features(selected_junctions)
        return selected_junctions_profile

    def are_known_introns(self, junctions, region):
        selected_junctions_profile = self.profile_for_junctions_introns(junctions, region)
        return all(el == 1 for el in selected_junctions_profile.read_profile)


class MisalignmentChecker:
    def __init__(self, reference, delta, max_intron_shift, max_missed_exon_len):
        if reference is not None:
            if not os.path.exists(reference + '.fai'):
                logger.info('Building fasta index with samtools')
                pysam.faidx(reference)
            reference = pysam.FastaFile(reference)
        self.reference = reference
        self.delta = delta
        self.max_intron_shift = max_intron_shift
        self.max_missed_exon_len = max_missed_exon_len
        self.scores = (2, -1, -1, -.2)

    def classify_skipped_exon(self, isoform_junctions, isoform_cregion, read_junctions, read_cregion,
                              total_intron_len_diff, read_introns_known, total_exon_len, alignment):
        assert len(set(read_cregion)) == 1

        read_cpos = read_cregion[0]
        read_left, read_right = read_junctions[read_cpos]
        l, r = isoform_cregion
        iso_left, iso_right = isoform_junctions[l][0], isoform_junctions[r][1]
        exon_left, exon_right = isoform_junctions[l][1], isoform_junctions[r][0]
        if total_intron_len_diff < 2 * self.delta and total_exon_len <= self.max_missed_exon_len:
            if alignment is None or self.reference is None:
                return MatchEventSubtype.exon_misallignment

            if iso_left - self.delta <= read_left and iso_right + self.delta >= read_right:
                seq = self.get_read_sequence(alignment, (iso_left, read_left)) \
                          + self.get_read_sequence(alignment, (read_right, iso_right))
                ref_seq = self.reference.fetch(alignment.reference_name, exon_left, exon_right)
                if self.sequences_match(seq, ref_seq):
                    return MatchEventSubtype.exon_misallignment
                else:
                    return MatchEventSubtype.exon_skipping_novel_intron

        if read_introns_known:
            return MatchEventSubtype.exon_skipping_known_intron
        return MatchEventSubtype.exon_skipping_novel_intron

    def classify_intron_misalignment(self, read_junctions, isoform_junctions, read_cpos, isoform_cpos, alignment=None):

        read_left, read_right = read_junctions[read_cpos]
        iso_left, iso_right = isoform_junctions[isoform_cpos]

        if abs(read_left - iso_left) + abs(read_right - iso_right) <= 2 * self.delta:
            return MatchEventSubtype.intron_shift

        if alignment is None or self.reference is None:
            if abs(read_left - iso_left) <= self.max_intron_shift:
                return MatchEventSubtype.intron_shift
            return MatchEventSubtype.intron_alternation_novel

        read_region, ref_region = self.get_aligned_regions_intron(read_left, read_right, iso_left, iso_right)
        seq = self.get_read_sequence(alignment, read_region)
        ref_seq = self.reference.fetch(alignment.reference_name, *ref_region)

        if self.sequences_match(seq, ref_seq):
            return MatchEventSubtype.intron_shift
        return MatchEventSubtype.intron_alternation_novel

    def sequences_match(self, seq, ref_seq):
        score, size = 0, 1
        for a in pairwise2.align.globalms(seq, ref_seq, *self.scores):
            score, size = a[2], a[4]
        if score > 0.7 * size:
            return True
        return False

    @staticmethod
    def get_aligned_regions_intron(read_left, read_right, iso_left, iso_right):
        if iso_left <= read_left:
            return (iso_left, read_left), (iso_right, read_right)
        return (read_right, iso_right), (read_left, iso_left)

    @staticmethod
    def get_read_sequence(alignment, ref_region):
        ref_start, ref_end = ref_region
        seq = ''
        last_ref_pos = 0
        for read_pos, ref_pos in alignment.get_aligned_pairs():
            if ref_start <= last_ref_pos <= ref_end:
                if read_pos is not None:
                    seq += alignment.seq[read_pos]
            last_ref_pos = ref_pos or last_ref_pos
        return seq

