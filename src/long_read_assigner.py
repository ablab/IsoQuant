############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import copy
from collections import namedtuple

from src.common import *
from src.gene_info import *
from src.long_read_profiles import *


logger = logging.getLogger('IsoQuant')
IsoformDiff = namedtuple("IsoformDiff", ("id", "diff"))


class AssignmentType:
    unique = "unique"
    empty = "empty"
    ambiguous = "ambiguous"
    contradictory = "contradictory"


class MatchingEvent:
    none = "none"
    undefined = "undefined"
    unspliced = "unspliced"
    no_contradiction = "no_contradiction"

    fsm = "full_splice_match"
    ambiguous_isoform = "assigned_to_ambiguous_isoform"

    extra_intron = "additional_novel_intron"
    extra_intron_known = "additional_known_intron"
    extra_intron_out = "additional_terminal_intron"

    intron_retention = "intron_retention"
    intron_migration = "intron_migration"
    intron_alternation_novel = "intron_change_to_novel"
    intron_alternation_known = "intron_change_to_known"
    donor_site = "donor_site"
    acceptor_site = "acceptor_site"
    mutually_exclusive_exons_novel = "mutualy_exclusive_novel_exons"
    mutually_exclusive_exons_known = "mutualy_exclusive_known_exons"
    exon_skipping_known_intron = "exon_skipping_known_intron"
    exon_skipping_novel_intron = "exon_skipping_novel_intron"
    exon_gain_known = "gains_known_exon"
    exon_gain_novel = "gains_novel_exon"
    alternative_structure_novel = "alternative_structure_novel_introns"
    alternative_structure_known = "alternative_structure_known_introns"

    exon_elongation = "exon_elongation"
    intron_shift = "intron_shift"
    exon_misallignment = "small_exon_misallignment"

    @staticmethod
    def is_elongated_exon(event):
        return event.find(MatchingEvent.exon_elongation) != 1

    @staticmethod
    def is_extra(event):
        return any(event.find(e) != -1 for e in [MatchingEvent.extra_intron_out, MatchingEvent.extra_intron_known, MatchingEvent.extra_intron])

    @staticmethod
    def is_minor_error(event):
        return event == MatchingEvent.intron_shift or event == MatchingEvent.exon_misallignment

    @staticmethod
    def concat_events(event1, event2):
        if event1 == MatchingEvent.none:
            return event2
        elif event2 == MatchingEvent.none:
            return event1
        else:
            return event1 + '_' + event2


class ReadAssignment:
    def __init__(self, read_id, assignment_type, assigned_features = [], match_events = []):
        self.read_id = read_id
        self.assignment_type = copy.deepcopy(assignment_type)
        if len(match_events) == 0 and len(assigned_features) > 0:
            self.match_events = copy.deepcopy([copy.deepcopy(MatchingEvent.none)] * len(assigned_features))
        else:
            self.match_events = copy.deepcopy(match_events)
        self.assigned_features = copy.deepcopy(assigned_features)

    def add_feature(self, feature, event = MatchingEvent.none):
        self.assigned_features.append(feature)
        self.match_events.append(copy.deepcopy(event))

    def add_event(self, index, event):
        if self.match_events[index] == MatchingEvent.none:
            self.match_events[index] = copy.deepcopy(event)
        else:
            self.match_events[index] += "+" + copy.deepcopy(event)

    def add_common_event(self, event):
        for i in range(len(self.match_events)):
            self.add_event(i, event)

    def set_assignment_type(self, assignment_type):
        self.assignment_type = copy.deepcopy(assignment_type)


class LongReadAssigner:

    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params

    def match_profile(self, read_gene_profile, isoform_profiles, hint=None):
        """ match read profiles to a known isoform junction profile

        Parameters
        ----------
        read_gene_profile: list of int
        isoform_profiles: dict of str
        hint: set, optional
            potential candidates, all other profiles will be skipped

        Returns
        -------
        result: list of IsoformDiff
            list of tuples (id and difference)
        """
        isoforms = []
        for isoform_id, isoform_profile in isoform_profiles.items():
            if hint and isoform_id not in hint:
                continue
            diff = difference_in_present_features(isoform_profile, read_gene_profile)
            isoforms.append(IsoformDiff(isoform_id, diff))
        return sorted(isoforms, key=lambda x: x[1])

    def find_matching_isoforms(self, read_gene_profile, isoform_profiles, hint=None):
        """ match read profiles to a known isoform junction profile

        Parameters
        ----------
        read_gene_profile: list of int
        isoform_profiles: dict of str
        hint: set, optional
            potential candidates, all other profiles will be skipped

        Returns
        -------
        result: set of str
            matching isoforms

        """
        isoforms = self.match_profile(read_gene_profile, isoform_profiles, hint)
        return set(map(lambda x: x[0], filter(lambda x: x[1] == 0, isoforms)))

    # === Isoforom matching function ===
    def assign_to_isoform(self, read_id, combined_read_profile):
        """ assign read to isoform according to it

        Parameters
        ----------
        read_id: str
        combined_read_profile: CombinedReadProfiles

        Returns
        -------

        """
        read_intron_profile = combined_read_profile.read_intron_profile
        read_split_exon_profile = combined_read_profile.read_split_exon_profile

        logger.debug("> Read blocks" + str(read_split_exon_profile.read_features))
        logger.debug("Gene split exons" + str(self.gene_info.split_exon_profiles.features))
        logger.debug("Read exons profile" + str(read_split_exon_profile.read_profile))
        logger.debug("Gene exons profile" + str(read_split_exon_profile.gene_profile))

        logger.debug("< Read introns" + str(read_intron_profile.read_features))
        logger.debug("Gene introns" + str(self.gene_info.intron_profiles.features))
        logger.debug("Read intron profile" + str(read_intron_profile.read_profile))
        logger.debug("Gene intron profile" + str(read_intron_profile.gene_profile))

        if all(el == 0 for el in read_split_exon_profile.read_profile) or all(el == 0 for el in read_split_exon_profile.gene_profile):
            # none of the blocks matched
            logger.debug("EMPTY")
            return ReadAssignment(read_id, AssignmentType.empty)
        elif any(el == -1 for el in read_intron_profile.read_profile) or any(el == -1 for el in read_split_exon_profile.read_profile):
            # Read has missing exons / introns
            logger.debug("+ Has contradictory features")
            return self.resolve_contradictory(read_id, combined_read_profile)
        elif any(el == 0 for el in read_intron_profile.read_profile) or any(el == 0 for el in read_split_exon_profile.read_profile):
            # Read has extra flanking exons / introns
            logger.debug("+ Has extra flanking features")
            return self.resolve_extra_flanking(read_id, combined_read_profile)
        else:
            logger.debug("+ No contradictory features")
            assignment = self.match_non_contradictory(read_id, combined_read_profile)

            #check for extra flanking sequences
            if assignment.assignment_type == AssignmentType.unique and not MatchingEvent.is_extra(assignment.match_events[0]) \
                    and not MatchingEvent.is_elongated_exon(assignment.match_events[0]):
                self.resolve_exon_elongation(read_split_exon_profile, assignment)

            return assignment

    def is_fsm(self, read_intron_profile, isoform_id):
        read_profile = read_intron_profile.gene_profile
        isoform_profile = self.gene_info.intron_profiles.profiles[isoform_id]
        return all(el == 1 for el in read_intron_profile.read_profile) \
               and 1 in read_profile and 1 in isoform_profile and \
               read_profile.index(1) == isoform_profile.index(1) and \
               rindex(read_profile, 1) == rindex(isoform_profile, 1)

    def resolve_exon_elongation(self, read_split_exon_profile, assignment, compare_to_isoform=False, assignment_index=0):
        split_exons = self.gene_info.split_exon_profiles.features

        if compare_to_isoform:
            isofrom_id = assignment.assigned_features[assignment_index]
            isofrom_profile = self.gene_info.split_exon_profiles.profiles[isofrom_id]
            start_exon = split_exons[isofrom_profile.index(1)]
            end_exon = split_exons[rindex(isofrom_profile, 1)]
        else:
            start_exon = split_exons[read_split_exon_profile.gene_profile.index(1)]
            end_exon = split_exons[rindex(read_split_exon_profile.gene_profile, 1)]
        read_start = read_split_exon_profile.read_features[0][0]
        read_end = read_split_exon_profile.read_features[-1][1]

        extra5 = start_exon[0] - read_start
        extra3 = read_end - end_exon[1]

        logger.debug("+ + Checking exon elongation")
        logger.debug("Read: " + str(read_start) + "-" + str(read_end))
        logger.debug("Isoform: " + str(read_start) + "-" + str(read_end))
        if extra3 < self.params.delta and extra5 < self.params.delta:
            logger.debug("+ + None")
            return assignment
        else:
            logger.debug("+ + Minor")
            assignment.add_event(assignment_index, MatchingEvent.exon_elongation)
            if extra3 > self.params.max_exon_extension or extra5 > self.params.max_exon_extension:
                logger.debug("+ + No, serious")
                assignment.set_assignment_type(AssignmentType.contradictory)

    def match_non_contradictory(self, read_id, combined_read_profile, has_zeros=False):
        """ match profile when all read features are assigned

        Parameters
        ----------
        read_id: str
        combined_read_profile: CombinedReadProfiles
        has_zeros: bool, optional

        Returns
        -------

        """
        read_intron_profile = combined_read_profile.read_intron_profile
        read_split_exon_profile = combined_read_profile.read_split_exon_profile

        intron_matched_isoforms = self.find_matching_isoforms(read_intron_profile.gene_profile,
                                                              self.gene_info.intron_profiles.profiles)
        exon_matched_isoforms = self.find_matching_isoforms(read_split_exon_profile.gene_profile,
                                                            self.gene_info.split_exon_profiles.profiles)
        both_mathched_isoforms = intron_matched_isoforms.intersection(exon_matched_isoforms)

        logger.debug("Intron matched " + str(intron_matched_isoforms))
        logger.debug("Exon matched " + str(exon_matched_isoforms))
        logger.debug("Both matched " + str(both_mathched_isoforms))

        read_assignment = None
        if len(both_mathched_isoforms) == 1:
            isoform_id = list(both_mathched_isoforms)[0]

            logger.debug("+ + UNIQUE match found " + isoform_id)
            logger.debug("Exon profile: " + str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
            logger.debug("Intron profile: " + str(self.gene_info.intron_profiles.profiles[isoform_id]))

            if not has_zeros and self.is_fsm(read_intron_profile, isoform_id):
                logger.debug("+ + Full splice match " + isoform_id)
                read_assignment = ReadAssignment(read_id, AssignmentType.unique, list(both_mathched_isoforms), [MatchingEvent.fsm])
            else:
                read_assignment = ReadAssignment(read_id, AssignmentType.unique, list(both_mathched_isoforms))

        elif len(both_mathched_isoforms) > 1:

            counter = 0
            for isoform_id in both_mathched_isoforms:
                logger.debug("+ + AMBIGUOUS match " + str(counter) + ": " + isoform_id)
                counter += 1
                logger.debug("Exon profile: " + str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
                logger.debug("Intron profile: " + str(self.gene_info.intron_profiles.profiles[isoform_id]))

            if not self.params.resolve_ambiguous or len(self.gene_info.ambiguous_isoforms.intersection(both_mathched_isoforms)) > 0:
                read_assignment = ReadAssignment(read_id, AssignmentType.ambiguous, list(both_mathched_isoforms))
            else:
                logger.debug("+ + Resolving ambiguity ")
                possible_isoforms = self.gene_info.ambiguous_isoforms.intersection(both_mathched_isoforms)
                intron_matched_isoforms = \
                    self.resolve_ambiguous(read_intron_profile.gene_profile, self.gene_info.intron_profiles.profiles, possible_isoforms)
                exon_matched_isoforms = \
                    self.resolve_ambiguous(read_split_exon_profile.gene_profile, self.gene_info.split_exon_profiles.profile, possible_isoforms)

                both_mathched_isoforms = intron_matched_isoforms.intersection(exon_matched_isoforms)

                if len(both_mathched_isoforms) == 1:
                    read_assignment = ReadAssignment(read_id, AssignmentType.unique, list(both_mathched_isoforms), [MatchingEvent.ambiguous_isoform])
                else:
                    read_assignment = ReadAssignment(read_id, AssignmentType.ambiguous, list(both_mathched_isoforms))

        elif len(both_mathched_isoforms) == 0:
            logger.debug("+ + Inersection is empty ")

            if len(intron_matched_isoforms) > 0:
                #intron profile match, but extra exons
                if len(intron_matched_isoforms) == 1:
                    isoform_id = list(intron_matched_isoforms)[0]

                    logger.debug("+ + But there is unique intron profile match "  + isoform_id)
                    logger.debug(str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
                    logger.debug(str(self.gene_info.intron_profiles.profiles[isoform_id]))

                    read_assignment = ReadAssignment(read_id, AssignmentType.unique, list(intron_matched_isoforms))
                    if not has_zeros:
                        if self.is_fsm(read_intron_profile, isoform_id):
                            read_assignment.add_event(0, MatchingEvent.fsm)
                        self.resolve_exon_elongation(read_split_exon_profile, read_assignment, compare_to_isoform=True)
                else:
                    logger.debug("+ + But there is amb intron profile match")
                    counter = 0
                    for isoform_id in intron_matched_isoforms:
                        logger.debug("+ + AMBIGUOUS match " + str(counter) + ": " + isoform_id)
                        counter += 1
                        logger.debug("Exon profile: " + str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
                        logger.debug("Intron profile: " + str(self.gene_info.intron_profiles.profiles[isoform_id]))


                    read_assignment = ReadAssignment(read_id, AssignmentType.ambiguous, list(intron_matched_isoforms))
                    if not has_zeros:
                        for i in range(len(intron_matched_isoforms)):
                            self.resolve_exon_elongation(read_split_exon_profile, read_assignment,
                                                         compare_to_isoform=True, assignment_index=i)

            else:
                # alternative isoforms made of known introns/exons or intron retention
                logger.debug("+ + Resolving unmatched ")
                read_assignment = self.resolve_contradictory(read_id, combined_read_profile)
        return read_assignment

    # resolve assignment ambiguities caused by identical profiles
    def resolve_ambiguous(self, read_gene_profile, isoform_profiles, matched_isoforms):
        for t in matched_isoforms:
            matched_positions = find_matching_positions(isoform_profiles[t], read_gene_profile)

            all_positions_detected = True
            for i in range(len(matched_positions)):
                if matched_positions[i] == 0 and isoform_profiles[t][i] == 1:
                    all_positions_detected = False
                    break

            if all_positions_detected:
                return set([t])

        return matched_isoforms

    #resolve when there are 0s  at the ends of read profile
    def resolve_extra_flanking(self, read_id, combined_read_profile):
        read_intron_profile = combined_read_profile.read_intron_profile
        logger.debug("+ + " + str(read_intron_profile.read_profile))
        if read_intron_profile.read_profile[0] == 1 and read_intron_profile.read_profile[-1] == 1:
            logger.warning("+ + Both terminal introns present, odd case")

        assignment = self.match_non_contradictory(read_id, combined_read_profile, has_zeros=True)
        if not self.params.allow_extra_terminal_introns:
            assignment.set_assignment_type(AssignmentType.contradictory)
        assignment.add_common_event(MatchingEvent.extra_intron_out)

        return assignment

    #resolve when there are -1s in read profile or when there are no exactly matching isoforms, but no -1s in read profiles
    def resolve_contradictory(self, read_id, combined_read_profile):
        read_intron_profile = combined_read_profile.read_intron_profile
        read_split_exon_profile = combined_read_profile.read_split_exon_profile

        intron_matching_isoforms = self.match_profile(read_intron_profile.gene_profile, self.gene_info.intron_profiles.profiles)
        exon_matching_isoforms = self.match_profile(read_split_exon_profile.gene_profile, self.gene_info.split_exon_profiles.profiles)
        return self.detect_differences(read_id, read_intron_profile, read_split_exon_profile, intron_matching_isoforms, exon_matching_isoforms)

    def detect_differences(self, read_id, read_intron_profile, read_split_exon_profile,
                           intron_matching_isoforms, exon_matching_isoforms):
        """ compare read to closest matching isoforms

        Parameters
        ----------
        read_id: str
        read_intron_profile: MappedReadProfile
        read_split_exon_profile: MappedReadProfile
        intron_matching_isoforms: list of IsoformDiff
        exon_matching_isoforms: list of IsoformDiff

        Returns
        -------
        result: ReadAssignment
        """
        # get isoforms that have closes intron and exon profiles
        best_intron_isoform_ids = get_first_best_from_sorted(intron_matching_isoforms)
        best_exon_isoforms = list(filter(lambda x: x[0] in best_intron_isoform_ids, exon_matching_isoforms))
        best_isoform_ids = get_first_best_from_sorted(best_exon_isoforms)

        assignment = ReadAssignment(read_id, AssignmentType.contradictory)
        logger.debug("+ + Closest matching isoforms " + str(best_isoform_ids))
        isoform_id = list(best_isoform_ids)[0]
        logger.debug(str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
        logger.debug(str(self.gene_info.intron_profiles.profiles[isoform_id]))

        for isoform_id in best_isoform_ids:
            # get intron coordinates
            isoform_introns = get_blocks_from_profile(self.gene_info.intron_profiles.features,
                                                      self.gene_info.intron_profiles.profiles[isoform_id])
            # read start-end coordinates
            read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[-1][1])
            # isoform start-end
            isoform_exon_profile = self.gene_info.split_exon_profiles.profiles[isoform_id]
            isoform_start = self.gene_info.split_exon_profiles.features[isoform_exon_profile.index(1)][0]
            isoform_end = self.gene_info.split_exon_profiles.features[rindex(isoform_exon_profile, 1)][1]
            isoform_region = (isoform_start, isoform_end)
            matching_event = self.compare_junctions(read_intron_profile.read_features, read_region,
                                                    isoform_introns, isoform_region)
            logger.debug("+ + Existing features " + str(assignment.assigned_features))
            assignment.add_feature(isoform_id, matching_event)
            logger.debug("+ + Found contradiction for " + isoform_id + ": " + matching_event)

        if all(e == MatchingEvent.no_contradiction for e in assignment.match_events):
            # No contradiction
            new_assignment_type = AssignmentType.unique if len(best_isoform_ids) == 1 else AssignmentType.ambiguous
            return ReadAssignment(read_id, new_assignment_type, assignment.assigned_features, assignment.match_events)

        if self.params.correct_minor_errors:
            #correct assignment type if contradiction is minor
            if all(MatchingEvent.is_minor_error(e) for e in assignment.match_events):
                new_assignment_type = AssignmentType.unique if len(best_isoform_ids) == 1 else AssignmentType.ambiguous
                assignment.set_assignment_type(new_assignment_type)

        return assignment

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

        """
        if len(read_junctions) == 0:
            return MatchingEvent.unspliced
        read_pos = 0
        isoform_pos = 0
        read_features_present = [0 for i in range(0, len(read_junctions))]
        isoform_features_present = [0 for i in range(0, len(isoform_junctions))]
        contradictory_region_pairs = []
        current_contradictory_region = (None, None)

        while read_pos < len(read_junctions) and isoform_pos < len(isoform_junctions):
            if equal_ranges(isoform_junctions[isoform_pos], read_junctions[read_pos], self.params.delta):
                # junctions are equal
                read_features_present[read_pos] = 1
                isoform_features_present[isoform_pos] = 1
                if (current_contradictory_region != (None, None)):
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = (None, None)
                read_pos += 1
                isoform_pos += 1

            elif overlaps(isoform_junctions[isoform_pos], read_junctions[read_pos]):
                # junctions overlap, but are unequal
                read_features_present[read_pos] = -1
                isoform_features_present[isoform_pos] = -1
                if current_contradictory_region == (None, None):
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
                if (current_contradictory_region != (None, None)):
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = (None, None)
                if read_pos > 0 or contains(read_region, isoform_junctions[isoform_pos]):
                    if (isoform_features_present[isoform_pos] != -1):
                        contradictory_region_pairs.append((None, (isoform_pos, isoform_pos)))
                    isoform_features_present[isoform_pos] = -1
                isoform_pos += 1

            else:
                # read junction is behind, move on
                if (current_contradictory_region != (None, None)):
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = (None, None)
                if isoform_pos > 0 or contains(isoform_region, read_junctions[read_pos]):
                    if (read_features_present[read_pos] != -1):
                        contradictory_region_pairs.append(((read_pos, read_pos), None))
                    read_features_present[read_pos] = -1
                read_pos += 1

        if (current_contradictory_region != (None, None)):
            contradictory_region_pairs.append(current_contradictory_region)

        # check terminating regions
        while read_pos < len(read_junctions):
            if contains(isoform_region, read_junctions[read_pos]):
                if (read_features_present[read_pos] != -1):
                    contradictory_region_pairs.append(((read_pos, read_pos), None))
                    read_features_present[read_pos] = -1
            else:
                break
            read_pos += 1

        while isoform_pos < len(isoform_junctions):
            if contains(read_region, isoform_junctions[isoform_pos]):
                if (isoform_features_present[isoform_pos] != -1):
                    contradictory_region_pairs.append((None, (isoform_pos, isoform_pos)))
                    isoform_features_present[isoform_pos] = -1
            else:
                break
            isoform_pos += 1

        logger.debug("+ + Inspected contradictory read")
        logger.debug("+ + Read profile " + str(read_features_present))
        logger.debug("+ + Read introns " + str(read_junctions))
        logger.debug("+ + Read region " + str(read_region))

        logger.debug("+ + Isoform profile " + str(isoform_features_present))
        logger.debug("+ + Isoform introns " + str(isoform_junctions))
        logger.debug("+ + Isoform region " + str(isoform_region))

        if any(el == -1 for el in read_features_present) or any(el == -1 for el in isoform_features_present):
            # classify contradictions
            logger.debug("+ + Classifying contradictions")
            return self.detect_contradiction_type(read_junctions, isoform_junctions, contradictory_region_pairs)
        elif read_features_present[0] == 0 or read_features_present[-1] == 0:
            logger.debug("+ + Found only extra terminal introns")
            return MatchingEvent.extra_intron_out
        else:
            logger.debug("No contradition detected, odd case")
            return MatchingEvent.no_contradiction

    def detect_contradiction_type(self, read_junctions, isoform_junctions, contradictory_region_pairs):
        """

        Parameters
        ----------
        read_junctions: list of tuples of int
        isoform_junctions: list of tuples of int
        contradictory_region_pairs: list of tuples of int

        Returns
        -------
        result: str
        """
        contradiction_events = []
        for pair in contradictory_region_pairs:
            # classify each contradictory area separately
            event = self.compare_overlapping_contradictional_regions(read_junctions, isoform_junctions, pair[0], pair[1])
            contradiction_events.append(event)

        contradiction_events_set = set(contradiction_events)
        if len(contradiction_events_set) == 1:
            return contradiction_events[0]
        else:
            return "+".join(list(contradiction_events_set))

    def compare_overlapping_contradictional_regions(self, read_junctions, isoform_junctions, read_cregion, isoform_cregion):
        if read_cregion is None:
            return MatchingEvent.intron_retention
        elif isoform_cregion is None:
            if self.are_known_introns(read_junctions, read_cregion):
                return MatchingEvent.extra_intron_known
            return MatchingEvent.extra_intron

        read_intron_total_len = sum(
            [read_junctions[i][1] - read_junctions[i][0] for i in range(read_cregion[0], read_cregion[1] + 1)])
        isoform_intron_total_len = sum(
            [isoform_junctions[i][1] - isoform_junctions[i][0] for i in range(isoform_cregion[0], isoform_cregion[1] + 1)])
        total_intron_len_diff = abs(read_intron_total_len - isoform_intron_total_len)

        read_introns_known = self.are_known_introns(read_junctions, read_cregion)

        if read_cregion[1] == read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            return self.classify_single_intron_alteration(read_junctions, isoform_junctions, read_cregion[0],
                                                          isoform_cregion[0], total_intron_len_diff, read_introns_known)

        elif read_cregion[1] - read_cregion[0] == isoform_cregion[1] - isoform_cregion[0] and \
                total_intron_len_diff < self.params.delta:
            if read_introns_known:
                return MatchingEvent.mutually_exclusive_exons_known
            else:
                return MatchingEvent.mutually_exclusive_exons_novel

        elif read_cregion[1] == read_cregion[0] and isoform_cregion[1] > isoform_cregion[0]:
            return self.classify_skipped_exons(isoform_junctions, isoform_cregion, total_intron_len_diff, read_introns_known)

        elif read_cregion[1] > read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            if read_introns_known:
                return MatchingEvent.exon_gain_known
            else:
                return MatchingEvent.exon_gain_novel

        else:
            if read_introns_known:
                return MatchingEvent.alternative_structure_known
            else:
                return MatchingEvent.alternative_structure_novel

    def classify_skipped_exons(self, isoform_junctions, isoform_cregion,
                               total_intron_len_diff, read_introns_known):
        total_exon_len = sum([isoform_junctions[i + 1][0] - isoform_junctions[i][1] + 1
                              for i in range(isoform_cregion[0], isoform_cregion[1])])

        if total_intron_len_diff < 2 * self.params.delta and total_exon_len <= self.params.max_missed_exon_len:
            return MatchingEvent.exon_misallignment
        else:
            if read_introns_known:
                return MatchingEvent.exon_skipping_known_intron
            else:
                return MatchingEvent.exon_skipping_novel_intron

    def classify_single_intron_alteration(self, read_junctions, isoform_junctions, read_cpos, isoform_cpos,
                                          total_intron_len_diff, read_introns_known):
        if total_intron_len_diff <= 2 * self.params.delta:
            if read_introns_known:
                return MatchingEvent.intron_migration
            else:
                if abs(isoform_junctions[isoform_cpos][0] - read_junctions[read_cpos][0]) <= self.params.max_intron_shift:
                    return MatchingEvent.intron_shift
                else:
                    return MatchingEvent.intron_alternation_novel
        else:
            site = MatchingEvent.none
            if abs(isoform_junctions[isoform_cpos][0] - read_junctions[read_cpos][0]) < self.params.delta:
                site = MatchingEvent.acceptor_site
            elif abs(isoform_junctions[isoform_cpos][1] - read_junctions[read_cpos][1]) < self.params.delta:
                site = MatchingEvent.donor_site

            if read_introns_known:
                return MatchingEvent.concat_events(MatchingEvent.intron_alternation_known, site)
            else:
                return MatchingEvent.concat_events(MatchingEvent.intron_alternation_novel, site)

    def are_known_introns(self, junctions, region):
        selected_junctions = []
        logger.debug("Checking for known introns " + str(region))
        for i in range(region[0], region[1] + 1):
            selected_junctions.append(junctions[i])

        logger.debug(str(selected_junctions))
        intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.intron_profiles.features,
                                                  (self.gene_info.start, self.gene_info.end),
                                                  comparator = partial(equal_ranges, delta = self.params.delta))
        selected_junctions_profile = intron_profile_constructor.construct_profile_for_introns(selected_junctions)
        logger.debug(str(selected_junctions_profile.read_profile))
        return all(el == 1 for el in selected_junctions_profile.read_profile)


    # === Exon matching ==
    def assign_exons(self, combined_read_profile):
        pass