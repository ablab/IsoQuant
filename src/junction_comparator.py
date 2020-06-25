############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import namedtuple

from src.isoform_assignment import *
from src.gene_info import *
from src.long_read_profiles import *


logger = logging.getLogger('IsoQuant')


class JunctionComparator():
    absent = -10

    def __init__(self, params, intron_profile_constructor):
        self.params = params
        self.intron_profile_constructor = intron_profile_constructor

    def compare_junctions(self, read_junctions, read_region, isoform_junctions, isoform_region):
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
                if (current_contradictory_region != (self.absent, self.absent)):
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
                if (read_junctions[read_pos][1] < isoform_junctions[isoform_pos][1]):
                    read_pos += 1
                else:
                    isoform_pos += 1

            elif left_of(isoform_junctions[isoform_pos], read_junctions[read_pos]):
                # isoform junction is behind, move on
                if (current_contradictory_region != (self.absent, self.absent)):
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = (self.absent, self.absent)
                if read_pos > 0 or overlaps(read_region, isoform_junctions[isoform_pos]):
                    if (isoform_features_present[isoform_pos] != -1):
                        contradictory_region_pairs.append(((self.absent, self.absent), (isoform_pos, isoform_pos)))
                    isoform_features_present[isoform_pos] = -1
                isoform_pos += 1

            else:
                # read junction is behind, move on
                if (current_contradictory_region != (self.absent, self.absent)):
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = (self.absent, self.absent)
                if isoform_pos > 0 or overlaps(isoform_region, read_junctions[read_pos]):
                    if (read_features_present[read_pos] != -1):
                        contradictory_region_pairs.append(((read_pos, read_pos), (self.absent, isoform_pos)))
                    read_features_present[read_pos] = -1
                read_pos += 1

        if (current_contradictory_region != (self.absent, self.absent)):
            contradictory_region_pairs.append(current_contradictory_region)

        # check terminating regions
        while read_pos < len(read_junctions):
            if overlaps(isoform_region, read_junctions[read_pos]):
                if (read_features_present[read_pos] != -1):
                    contradictory_region_pairs.append(((read_pos, read_pos), (self.absent, isoform_pos)))
                    read_features_present[read_pos] = -1
            else:
                break
            read_pos += 1

        while isoform_pos < len(isoform_junctions):
            if overlaps(read_region, isoform_junctions[isoform_pos]):
                if (isoform_features_present[isoform_pos] != -1):
                    contradictory_region_pairs.append(((self.absent, self.absent), (isoform_pos, isoform_pos)))
                    isoform_features_present[isoform_pos] = -1
            else:
                break
            isoform_pos += 1

        logger.debug("+ + Inspected contradictory read")
        #logger.debug("+ + Read profile " + str(read_features_present))
        #logger.debug("+ + Read introns " + str(read_junctions))
        #logger.debug("+ + Read region " + str(read_region))

        #logger.debug("+ + Isoform profile " + str(isoform_features_present))
        #logger.debug("+ + Isoform introns " + str(isoform_junctions))
        #logger.debug("+ + Isoform region " + str(isoform_region))

        matching_events = []
        if any(el == -1 for el in read_features_present) or any(el == -1 for el in isoform_features_present):
            # classify contradictions
            logger.debug("+ + Classifying contradictions")
            matching_events = self.detect_contradiction_type(read_junctions, isoform_junctions, contradictory_region_pairs)

        if read_features_present[0] == 0 or read_features_present[-1] == 0:
            if all(x == 0 for x in read_features_present):
                return [make_event(MatchEventSubtype.undefined)]
            logger.debug("+ + Found only extra terminal introns ")
            self.add_extra_out_exon_events(matching_events, read_features_present, read_junctions, isoform_region[0])

        if len(matching_events) == 0:
            logger.debug("No contradition detected")
            return [make_event(MatchEventSubtype.none)]
        return matching_events

    def detect_contradiction_type(self, read_junctions, isoform_junctions, contradictory_region_pairs):
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
        for pair in contradictory_region_pairs:
            # classify each contradictory area separately
            event = self.compare_overlapping_contradictional_regions(read_junctions, isoform_junctions, pair[0], pair[1])
            contradiction_events.append(event)

        return list(set(contradiction_events))

    def compare_overlapping_contradictional_regions(self, read_junctions, isoform_junctions, read_cregion, isoform_cregion):
        if read_cregion[0] == self.absent:
            return make_event(MatchEventSubtype.intron_retention, isoform_cregion[0], read_cregion)
        elif isoform_cregion[0] == self.absent:
            #intron_start = read_junctions[read_cregion[0]]
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
            event = self.classify_single_intron_alternation(read_junctions, isoform_junctions, read_cregion[0],
                                                           isoform_cregion[0], total_intron_len_diff, read_introns_known)

        elif read_cregion[1] - read_cregion[0] == isoform_cregion[1] - isoform_cregion[0] and \
                total_intron_len_diff < self.params.delta:
            if read_introns_known:
                event = MatchEventSubtype.mutually_exclusive_exons_known
            else:
                event = MatchEventSubtype.mutually_exclusive_exons_novel

        elif read_cregion[1] == read_cregion[0] and isoform_cregion[1] > isoform_cregion[0]:
            event = self.classify_skipped_exons(isoform_junctions, isoform_cregion, total_intron_len_diff, read_introns_known)

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

    def classify_skipped_exons(self, isoform_junctions, isoform_cregion,
                               total_intron_len_diff, read_introns_known):
        total_exon_len = sum([isoform_junctions[i + 1][0] - isoform_junctions[i][1] + 1
                              for i in range(isoform_cregion[0], isoform_cregion[1])])

        if total_intron_len_diff < 2 * self.params.delta and total_exon_len <= self.params.max_missed_exon_len:
            event = MatchEventSubtype.exon_misallignment
        else:
            if read_introns_known:
                event = MatchEventSubtype.exon_skipping_known_intron
            else:
                event = MatchEventSubtype.exon_skipping_novel_intron
        return event

    def classify_single_intron_alternation(self, read_junctions, isoform_junctions, read_cpos, isoform_cpos,
                                           total_intron_len_diff, read_introns_known):
        if total_intron_len_diff <= 2 * self.params.delta:
            if read_introns_known:
                event = MatchEventSubtype.intron_migration
            else:
                if abs(isoform_junctions[isoform_cpos][0] - read_junctions[read_cpos][0]) <= self.params.max_intron_shift:
                    event = MatchEventSubtype.intron_shift
                else:
                    event = MatchEventSubtype.intron_alternation_novel
        else:
            # TODO correct when strand is negative
            if abs(isoform_junctions[isoform_cpos][0] - read_junctions[read_cpos][0]) < self.params.delta:
                event = MatchEventSubtype.alt_acceptor_site_known if read_introns_known \
                    else MatchEventSubtype.alt_acceptor_site_novel
            elif abs(isoform_junctions[isoform_cpos][1] - read_junctions[read_cpos][1]) < self.params.delta:
                event = MatchEventSubtype.alt_donor_site_known if read_introns_known \
                    else MatchEventSubtype.alt_donor_site_novel
            else:
                event = MatchEventSubtype.intron_alternation_known if read_introns_known \
                    else MatchEventSubtype.intron_alternation_novel
        return event

    def get_mono_exon_subtype(self, read_region, isoform_junctions):
        if not any(overlaps(read_region, rj) for rj in isoform_junctions):
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
            match_events.append(make_event(MatchEventSubtype.extra_intron_out_left,
                                           SupplementaryMatchConstansts.extra_left_mod_position, (0, read_pos - 1)))
        if extra_right:
            max_right = len(read_intron_read_profile) - 1
            read_pos = max_right
            while read_pos >= 0 and read_intron_read_profile[read_pos] == 0:
                read_pos -= 1
            match_events.append(make_event(MatchEventSubtype.extra_intron_out_right,
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