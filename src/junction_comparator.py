############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging

from src.isoform_assignment import *
from src.gene_info import *
from src.long_read_profiles import *


logger = logging.getLogger('IsoQuant')


class JunctionComparator:
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
            return self.get_mono_exon_subtype(read_region, isoform_junctions)

        read_pos = 0
        isoform_pos = 0
        read_features_present = [0 for i in range(0, len(read_junctions))]
        isoform_features_present = [0 for i in range(0, len(isoform_junctions))]
        contradictory_region_pairs = []
        absent_region = (SupplementaryMatchConstansts.absent_position, SupplementaryMatchConstansts.absent_position)
        current_contradictory_region = absent_region

        while read_pos < len(read_junctions) and isoform_pos < len(isoform_junctions):
            if equal_ranges(isoform_junctions[isoform_pos], read_junctions[read_pos], self.params.delta):
                # junctions are equal
                read_features_present[read_pos] = 1
                isoform_features_present[isoform_pos] = 1
                if current_contradictory_region != absent_region:
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = absent_region
                read_pos += 1
                isoform_pos += 1

            elif overlaps(isoform_junctions[isoform_pos], read_junctions[read_pos]):
                # junctions overlap, but are unequal
                read_features_present[read_pos] = -1
                isoform_features_present[isoform_pos] = -1
                if current_contradictory_region == absent_region:
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
                if current_contradictory_region != absent_region:
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = absent_region
                if read_pos > 0 or overlaps(read_region, isoform_junctions[isoform_pos]):
                    if isoform_features_present[isoform_pos] != -1:
                        contradictory_region_pairs.append(((SupplementaryMatchConstansts.absent_position, read_pos), (isoform_pos, isoform_pos)))
                    isoform_features_present[isoform_pos] = -1
                isoform_pos += 1

            else:
                # read junction is behind, move on
                if current_contradictory_region != absent_region:
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = absent_region
                if isoform_pos > 0 or overlaps(isoform_region, read_junctions[read_pos]):
                    if read_features_present[read_pos] != -1:
                        contradictory_region_pairs.append(((read_pos, read_pos), (SupplementaryMatchConstansts.absent_position, isoform_pos)))
                    read_features_present[read_pos] = -1
                read_pos += 1

        if current_contradictory_region != absent_region:
            contradictory_region_pairs.append(current_contradictory_region)

        # check terminating regions
        while read_pos < len(read_junctions):
            if overlaps(isoform_region, read_junctions[read_pos]):
                if read_features_present[read_pos] != -1:
                    contradictory_region_pairs.append(((read_pos, read_pos), (SupplementaryMatchConstansts.absent_position, isoform_pos)))
                    read_features_present[read_pos] = -1
            else:
                break
            read_pos += 1

        while isoform_pos < len(isoform_junctions):
            if overlaps(read_region, isoform_junctions[isoform_pos]):
                if isoform_features_present[isoform_pos] != -1:
                    contradictory_region_pairs.append(((SupplementaryMatchConstansts.absent_position, read_pos), (isoform_pos, isoform_pos)))
                    isoform_features_present[isoform_pos] = -1
            else:
                break
            isoform_pos += 1

        # logger.debug("+ + Inspected contradictory read")
        # logger.debug("+ + Read profile " + str(read_features_present))
        # logger.debug("+ + Read introns " + str(read_junctions))
        # logger.debug("+ + Read region " + str(read_region))

        # logger.debug("+ + Isoform profile " + str(isoform_features_present))
        # logger.debug("+ + Isoform introns " + str(isoform_junctions))
        # logger.debug("+ + Isoform region " + str(isoform_region))

        matching_events = []
        # logger.debug(read_features_present)
        # logger.debug(isoform_features_present)
        if any(el == -1 for el in read_features_present) or any(el == -1 for el in isoform_features_present):
            # classify contradictions
            # logger.debug("+ + Classifying contradictions")
            matching_events = self.detect_contradiction_type(read_region, read_junctions, isoform_region, isoform_junctions, contradictory_region_pairs)

        if read_features_present[0] == 0 or read_features_present[-1] == 0:
            # if all(x == 0 for x in read_features_present):
            #    return [MatchEvent(MatchEventSubtype.undefined)]
            # logger.debug("+ + Found only extra terminal introns ")
            self.add_extra_out_exon_events(matching_events, read_features_present, read_region, read_junctions, isoform_region[0])

        if len(matching_events) == 0:
            # logger.debug("No contradiction detected")
            return [MatchEvent(MatchEventSubtype.none)]
        return matching_events

    def detect_contradiction_type(self, read_region, read_junctions, isoform_region, isoform_junctions, contradictory_region_pairs):
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
            event = self.compare_overlapping_contradictional_regions(read_region, read_junctions, isoform_region, isoform_junctions, pair[0], pair[1])
            if event is not None:
                contradiction_events.append(event)

        return contradiction_events

    def compare_overlapping_contradictional_regions(self, read_region, read_junctions,
                                                    isoform_region, isoform_junctions,
                                                    read_cregion, isoform_cregion):
        if read_cregion[0] == SupplementaryMatchConstansts.absent_position:
            if isoform_cregion[0] != isoform_cregion[1]:
                logger.warning("Multiple intron retentions in a single event:" + str(isoform_cregion))
            overlapped_isoform_intron = isoform_junctions[isoform_cregion[0]]
            if contains(read_region, overlapped_isoform_intron):
                if interval_len(overlapped_isoform_intron) <= self.params.micro_intron_length and \
                        contains_well_inside(
                            get_preceding_exon_from_junctions(read_region, read_junctions, read_cregion[1]),
                            overlapped_isoform_intron, self.params.minimal_exon_overlap):
                    return MatchEvent(MatchEventSubtype.fake_micro_intron_retention, isoform_cregion, read_cregion)
                return MatchEvent(MatchEventSubtype.intron_retention, isoform_cregion, read_cregion)
            elif overlaps_at_least(read_region, overlapped_isoform_intron, self.params.minor_exon_extension):
                if overlapped_isoform_intron[0] <= read_region[0]:
                    return MatchEvent(MatchEventSubtype.incomplete_intron_retention_left, isoform_cregion, read_cregion)
                else:
                    return MatchEvent(MatchEventSubtype.incomplete_intron_retention_right, isoform_cregion, read_cregion)
            else:
                return None

        elif isoform_cregion[0] == SupplementaryMatchConstansts.absent_position:
            if self.are_known_introns(read_junctions, read_cregion):
                return MatchEvent(MatchEventSubtype.extra_intron_known, isoform_cregion, read_cregion)
            elif self.are_suspicious_introns(read_region, read_junctions, read_cregion):
                return MatchEvent(MatchEventSubtype.none, isoform_cregion, read_cregion)
            elif read_cregion[0] == 0 and \
                    interval_len(get_exon(read_region, read_junctions, 0)) <= self.params.max_fake_terminal_exon_len:
                # we have extra intron to the left and first exons is short
                return MatchEvent(MatchEventSubtype.fake_terminal_exon_left,
                                  SupplementaryMatchConstansts.extra_left_region, read_cregion)
            elif read_cregion[1] == len(read_junctions) - 1 and \
                    interval_len(get_exon(read_region, read_junctions, -1)) <= self.params.max_fake_terminal_exon_len:
                # we have extra intron to the right and last exons is short
                return MatchEvent(MatchEventSubtype.fake_terminal_exon_right,
                                  SupplementaryMatchConstansts.extra_right_region, read_cregion)
            return MatchEvent(MatchEventSubtype.extra_intron_novel, isoform_cregion, read_cregion)

        read_intron_total_len = sum(interval_len(read_junctions[i])
                                    for i in range(read_cregion[0], read_cregion[1] + 1))
        isoform_intron_total_len = sum(interval_len(isoform_junctions[i])
                                       for i in range(isoform_cregion[0], isoform_cregion[1] + 1))
        total_intron_len_diff = abs(read_intron_total_len - isoform_intron_total_len)
        intron_length_is_similar = \
            total_intron_len_diff <= min(self.params.max_intron_abs_diff,
                                         self.params.max_intron_rel_diff * max(read_intron_total_len, isoform_intron_total_len))

        read_introns_known = self.are_known_introns(read_junctions, read_cregion)

        left_exons_overlap = overlaps(get_exon(read_region, read_junctions, read_cregion[0]),
                                      get_exon(isoform_region, isoform_junctions, isoform_cregion[0]))
        right_exons_overlap = overlaps(get_exon(read_region, read_junctions, read_cregion[1] + 1),
                                       get_exon(isoform_region, isoform_junctions, isoform_cregion[1] + 1))
        surrounded_by_exons = left_exons_overlap and right_exons_overlap
        read_left_site = read_junctions[read_cregion[0]][0]
        read_right_site = read_junctions[read_cregion[1]][1]
        isoform_left_site = isoform_junctions[isoform_cregion[0]][0]
        isoform_right_site = isoform_junctions[isoform_cregion[1]][1]
        similar_letf_site = abs(read_left_site - isoform_left_site) <= 2 * self.params.delta
        similar_right_site = abs(read_right_site - isoform_right_site) <= 2 * self.params.delta
        similar_bounds = similar_letf_site and similar_right_site
        read_introns_inside = contains_approx((isoform_left_site, isoform_right_site),
                                              (read_left_site, read_right_site), self.params.delta)
        isofrom_introns_inside = contains_approx((read_left_site, read_right_site),
                                                 (isoform_left_site, isoform_right_site), self.params.delta)

        event = None
        if surrounded_by_exons and read_cregion[1] == read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            # single exons is altered
            event = self.classify_single_intron_alternation(read_region, read_junctions, isoform_region,
                                                            isoform_junctions, read_cregion[0], isoform_cregion[0],
                                                            intron_length_is_similar, read_introns_known)

        elif len(read_junctions) > 1 and \
                read_cregion[1] == read_cregion[0] and isoform_cregion[1] == isoform_cregion[0] and \
                ((read_cregion[0] == 0 and isoform_cregion[0] == 0 and similar_right_site) or
                 (read_cregion[0] == len(read_junctions)-1 and isoform_cregion[0] == len(isoform_junctions) - 1
                  and similar_letf_site)):
            # terminal exon alternation
            if read_cregion[0] == 0 and isoform_cregion[0] == 0:
                read_exon = get_preceding_exon_from_junctions(read_region, read_junctions, 0)
                isoform_exon = get_preceding_exon_from_junctions(isoform_region, isoform_junctions, 0)
            else:
                read_exon = get_following_exon_from_junctions(read_region, read_junctions, -1)
                isoform_exon = get_following_exon_from_junctions(isoform_region, isoform_junctions, -1)
            if abs(interval_len(read_exon) - interval_len(isoform_exon)) < 2 * self.params.delta:
                # TODO: verify with alignment
                if read_cregion[0] == 0:
                    event = MatchEventSubtype.terminal_exon_misalignment_left
                else:
                    event = MatchEventSubtype.terminal_exon_misalignment_right
            elif read_introns_known:
                event = MatchEventSubtype.terminal_exon_shift_known
            else:
                event = MatchEventSubtype.terminal_exon_shift_novel

        elif surrounded_by_exons and similar_bounds and \
                read_cregion[1] - read_cregion[0] == isoform_cregion[1] - isoform_cregion[0] >= 1 and \
                total_intron_len_diff <= 2 * self.params.delta:
                # FIXME: increase threshold when re-alignment is implemented
            # several introns of the same total length as reference onces
            if read_introns_known:
                event = MatchEventSubtype.mutually_exclusive_exons_known
            else:
                event = MatchEventSubtype.mutually_exclusive_exons_novel

        elif surrounded_by_exons and read_introns_inside and \
                read_cregion[1] == read_cregion[0] and isoform_cregion[1] > isoform_cregion[0]:
            # skipped exon(s)
            event = self.classify_skipped_exons(isoform_junctions, isoform_cregion,
                                                intron_length_is_similar, read_introns_known, similar_bounds)

        elif surrounded_by_exons and similar_bounds and \
                read_cregion[1] > read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            # gain exon
            if read_introns_known:
                event = MatchEventSubtype.exon_gain_known
            elif self.are_suspicious_introns(read_region, read_junctions, read_cregion):
                event = MatchEventSubtype.intron_retention
            else:
                event = MatchEventSubtype.exon_gain_novel

        elif surrounded_by_exons and intron_length_is_similar and isofrom_introns_inside and \
                read_cregion[1] > read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            # exon detach
            if read_introns_known:
                event = MatchEventSubtype.exon_detach_known
            else:
                event = MatchEventSubtype.exon_detach_novel

        if event is None:
            # none of above, complex alternative structure
            if read_introns_known:
                event = MatchEventSubtype.alternative_structure_known
            elif surrounded_by_exons and self.are_suspicious_introns(read_region, read_junctions, read_cregion):
                event = MatchEventSubtype.intron_retention
            else:
                event = MatchEventSubtype.alternative_structure_novel
        return MatchEvent(event, isoform_cregion, read_cregion)

    def classify_skipped_exons(self, isoform_junctions, isoform_cregion,
                               intron_length_is_similar, read_introns_known, similar_bounds):
        total_exon_len = sum([isoform_junctions[i + 1][0] - isoform_junctions[i][1] + 1
                              for i in range(isoform_cregion[0], isoform_cregion[1])])
        event = None
        if intron_length_is_similar:
            if total_exon_len <= self.params.max_missed_exon_len:
                event = MatchEventSubtype.exon_misalignment
            elif read_introns_known:
                event = MatchEventSubtype.exon_merge_known
            else:
                event = MatchEventSubtype.exon_merge_novel
        elif similar_bounds:
            if read_introns_known:
                event = MatchEventSubtype.exon_skipping_known
            else:
                event = MatchEventSubtype.exon_skipping_novel
        return event

    def classify_single_intron_alternation(self, read_region, read_junctions, isoform_region, isoform_junctions,
                                           read_cpos, isoform_cpos, intron_length_is_similar, read_introns_known):
        if intron_length_is_similar:
            if abs(isoform_junctions[isoform_cpos][0] - read_junctions[read_cpos][0]) <= self.params.max_intron_shift:
                event = MatchEventSubtype.intron_shift
            elif read_introns_known:
                event = MatchEventSubtype.intron_migration
            else:
                event = MatchEventSubtype.intron_alternation_novel
        else:
            event = MatchEventSubtype.intron_alternation_known if read_introns_known \
                else MatchEventSubtype.intron_alternation_novel

            if abs(isoform_junctions[isoform_cpos][0] - read_junctions[read_cpos][0]) <= self.params.delta:
                # left splice sites are similar, check for exon overlap
                following_read_exon = get_following_exon_from_junctions(read_region, read_junctions, read_cpos)
                following_isoform_exon = get_following_exon_from_junctions(isoform_region, isoform_junctions, isoform_cpos)
                min_overlap = max(1, min(self.params.min_abs_exon_overlap,
                                         round(self.params.min_rel_exon_overlap * interval_len(following_isoform_exon))))
                if overlaps_at_least(following_read_exon, following_isoform_exon, min_overlap):
                    event = alternative_sites[("right", read_introns_known)]
            elif abs(isoform_junctions[isoform_cpos][1] - read_junctions[read_cpos][1]) <= self.params.delta:
                # right splice sites are similar, check for exon overlap
                preceding_read_exon = get_preceding_exon_from_junctions(read_region, read_junctions, read_cpos)
                preceding_isoform_exon = get_preceding_exon_from_junctions(isoform_region, isoform_junctions, isoform_cpos)
                min_overlap = max(1, min(self.params.min_abs_exon_overlap,
                                      round(self.params.min_rel_exon_overlap * interval_len(preceding_isoform_exon))))
                if overlaps_at_least(preceding_read_exon, preceding_isoform_exon, min_overlap):
                    event = alternative_sites[("left", read_introns_known)]

            if event in {MatchEventSubtype.intron_alternation_novel,
                         MatchEventSubtype.alt_left_site_novel,
                         MatchEventSubtype.alt_right_site_novel} and \
                self.are_suspicious_introns(read_region, read_junctions, (read_cpos, read_cpos)):
                event = MatchEventSubtype.intron_retention

        return event

    def are_suspicious_introns(self, read_region, read_junctions, read_cregion):
        total_intron_len = 0
        for cpos in range(read_cregion[0], read_cregion[1] + 1):
            if interval_len(read_junctions[cpos]) > self.params.max_suspicious_intron_abs_len:
                return False
            total_intron_len += interval_len(read_junctions[cpos])

        total_exon_len = 0
        for cpos in range(read_cregion[0], read_cregion[1] + 1):
            total_exon_len += interval_len(get_preceding_exon_from_junctions(read_region, read_junctions, cpos))
        total_exon_len += interval_len(get_following_exon_from_junctions(read_region, read_junctions, read_cregion[-1]))

        return float(total_intron_len) <= float(total_exon_len) * self.params.max_suspicious_intron_rel_len

    def get_mono_exon_subtype(self, read_region, isoform_junctions):
        if len(isoform_junctions) == 0:
            # both isoform and read are monoexon
            events = [MatchEvent(MatchEventSubtype.mono_exon_match)]
        else:
            # check if read captures some introns
            events = []
            for i, intron in enumerate(isoform_junctions):
                if contains(read_region, intron):
                    # full intron retention
                    if interval_len(intron) <= self.params.micro_intron_length and \
                            contains_well_inside(read_region, intron, self.params.minimal_exon_overlap):
                        events.append(MatchEvent(MatchEventSubtype.fake_micro_intron_retention, isoform_region=(i, i),
                                                 read_region=(SupplementaryMatchConstansts.absent_position, 0)))
                    else:
                        events.append(MatchEvent(MatchEventSubtype.unspliced_intron_retention, isoform_region=(i, i),
                                                 read_region=(SupplementaryMatchConstansts.absent_position, 0)))
                elif overlaps_at_least(read_region, intron, self.params.minor_exon_extension) and \
                        not contains(intron, read_region):
                    # partial IR
                    if intron[0] <= read_region[0]:
                        events.append(MatchEvent(MatchEventSubtype.incomplete_intron_retention_left, isoform_region=(i, i),
                                                 read_region=(SupplementaryMatchConstansts.absent_position, 0)))
                    else:
                        events.append(MatchEvent(MatchEventSubtype.incomplete_intron_retention_right, isoform_region=(i, i),
                                                 read_region=(SupplementaryMatchConstansts.absent_position, 0)))

        if not events:
            # monoexonic read without significant intron overlap
            events = [MatchEvent(MatchEventSubtype.mono_exonic)]
        return events

    def add_extra_out_exon_events(self, match_events, read_intron_read_profile, read_region, read_introns, isoform_start):
        extra_left = read_intron_read_profile[0] == 0
        extra_right = read_intron_read_profile[-1] == 0

        if all(x == 0 for x in read_intron_read_profile):
            if read_introns[0][0] < isoform_start:
                extra_right = False
            else:
                extra_left = False

        if extra_left:
            read_pos = 0
            if interval_len(get_exon(read_region, read_introns, read_pos)) <= self.params.max_fake_terminal_exon_len:
                match_events.append(MatchEvent(MatchEventSubtype.fake_terminal_exon_left,
                                               SupplementaryMatchConstansts.extra_left_region, (read_pos, read_pos)))
                read_pos += 1

            while read_pos < len(read_intron_read_profile) and read_intron_read_profile[read_pos] == 0:
                match_events.append(MatchEvent(MatchEventSubtype.extra_intron_flanking_left,
                                               SupplementaryMatchConstansts.extra_left_region, (read_pos, read_pos)))
                read_pos += 1

        if extra_right:
            read_pos = len(read_intron_read_profile) - 1
            if interval_len(get_exon(read_region, read_introns, read_pos + 1)) <= self.params.max_fake_terminal_exon_len:
                match_events.append(MatchEvent(MatchEventSubtype.fake_terminal_exon_right,
                                               SupplementaryMatchConstansts.extra_right_region, (read_pos, read_pos)))
                read_pos -= 1

            while read_pos >= 0 and read_intron_read_profile[read_pos] == 0:
                match_events.append(MatchEvent(MatchEventSubtype.extra_intron_flanking_right,
                                               SupplementaryMatchConstansts.extra_right_region, (read_pos, read_pos)))
                read_pos -= 1

    def profile_for_junctions_introns(self, junctions, region):
        selected_junctions = []
        # logger.debug("Checking for known introns " + str(region))
        for i in range(region[0], region[1] + 1):
            selected_junctions.append(junctions[i])

        selected_junctions_profile = self.intron_profile_constructor.construct_profile_for_features(selected_junctions)
        return selected_junctions_profile

    def are_known_introns(self, junctions, region):
        selected_junctions_profile = self.profile_for_junctions_introns(junctions, region)
        return all(el == 1 for el in selected_junctions_profile.read_profile)
