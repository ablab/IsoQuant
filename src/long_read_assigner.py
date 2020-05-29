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
IsoformDiff = namedtuple("IsoformDiff", ("id", "diff"))


class LongReadAssigner:

    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params

    ### ======== SUPPORT FUNCTIONS =======
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

    def get_gene_id(self, transcript_id):
        return self.gene_info.gene_id_map[transcript_id]

    # check for extra sequences and modify assignment accordingly
    def check_for_extra_terminal_seqs(self, read_split_exon_profile, assignment):
        for match in assignment.isoform_matches:
            exon_elongation_type = self.categorize_exon_elongation_subtype(read_split_exon_profile, match.assigned_transcript)
            if exon_elongation_type == MatchEventSubtype.extra_exon_out:
                # serious exon elongation
                match.add_subclassification(exon_elongation_type)
                if assignment.assignment_type == ReadAssignmentType.unique or \
                    assignment.assignment_type == ReadAssignmentType.minor:
                    assignment.set_assignment_type(ReadAssignmentType.contradictory)
            elif exon_elongation_type != MatchEventSubtype.none:
                # minor exon elongation
                match.add_subclassification(exon_elongation_type)
                if assignment.assignment_type == ReadAssignmentType.unique:
                    assignment.set_assignment_type(ReadAssignmentType.minor)

    # detect exon elongation subtyp
    def categorize_exon_elongation_subtype(self, read_split_exon_profile, isoform_id):
        split_exons = self.gene_info.split_exon_profiles.features
        isofrom_profile = self.gene_info.split_exon_profiles.profiles[isoform_id]

        # find first and last common exons
        common_first_exon = -1
        for i in range(len(split_exons)):
            if isofrom_profile[i] == 1 and read_split_exon_profile.gene_profile[i] == 1:
                common_first_exon = i
                break
        common_last_exon = -1
        for i in range(len(split_exons)):
            index = len(split_exons) - i - 1
            if isofrom_profile[index] == 1 and read_split_exon_profile.gene_profile[index] == 1:
                common_last_exon = index
                break
        if common_first_exon == -1 or common_last_exon == -1:
            logger.debug(" + Werid case for exon elongation, no matching exons")

        isoform_start = split_exons[common_first_exon][0]
        isoform_end = split_exons[common_last_exon][-1]
        read_start = read_split_exon_profile.read_features[0][0]
        read_end = read_split_exon_profile.read_features[-1][1]
        extra_left = isoform_start - read_start
        extra_right = read_end - isoform_end

        logger.debug("+ + Checking exon elongation")
        logger.debug("Read: " + str(read_start) + "-" + str(read_end))
        logger.debug("Isoform: " + str(read_start) + "-" + str(read_end))

        if extra_right <= self.params.delta and extra_left <= self.params.delta:
            logger.debug("+ + None")
            return MatchEventSubtype.none
        else:
            logger.debug("+ + Minor")
            if extra_right < self.params.max_exon_extension and extra_left < self.params.max_exon_extension:
                if extra_right > self.params.delta and extra_left > self.params.delta:
                    logger.debug(" + Exctra sequence on both ends")
                    return MatchEventSubtype.exon_elongation_both
                elif extra_right > self.params.delta:
                    logger.debug(" + Exctra sequence on left end")
                    if self.gene_info.isoform_strands[isoform_id] == "+":
                        return MatchEventSubtype.exon_elongation5
                    else:
                        return MatchEventSubtype.exon_elongation5
                else:
                    logger.debug(" + Exctra sequence on right end")
                    if self.gene_info.isoform_strands[isoform_id] == "-":
                        return MatchEventSubtype.exon_elongation5
                    else:
                        return MatchEventSubtype.exon_elongation5
            else:
                return MatchEventSubtype.extra_exon_out

    # get incompleteness type
    def detect_ism_subtype(self, read_intron_profile, isoform_id):
        if len(read_intron_profile.read_profile) == 0:
            logger.debug(" + Mono exon")
            return MatchEventSubtype.unspliced

        read_profile = read_intron_profile.gene_profile
        isoform_profile = self.gene_info.intron_profiles.profiles[isoform_id]
        is_left_truncated = left_truncated(read_profile, isoform_profile)
        is_right_truncated = right_truncated(read_profile, isoform_profile)

        if is_left_truncated and is_right_truncated:
            logger.debug(" + Internal")
            return MatchEventSubtype.ism_internal
        elif is_left_truncated:
            logger.debug(" + Truncated on the left")
            if self.gene_info.isoform_strands[isoform_id] == "+":
                return MatchEventSubtype.ism_5
            else:
                return MatchEventSubtype.ism_3
        elif is_right_truncated:
            logger.debug(" + Truncated on the right")
            if self.gene_info.isoform_strands[isoform_id] == "-":
                return MatchEventSubtype.ism_5
            else:
                return MatchEventSubtype.ism_3
        else:
            logger.debug(" + No ISM truncation, extra splice sites ")
            return MatchEventSubtype.unspliced

    # check where it is full splice match
    def is_fsm(self, read_intron_profile, isoform_id):
        read_profile = read_intron_profile.gene_profile
        isoform_profile = self.gene_info.intron_profiles.profiles[isoform_id]
        return all(el == 1 for el in read_intron_profile.read_profile) \
               and not left_truncated(read_profile, isoform_profile) \
               and not right_truncated(read_profile, isoform_profile)

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

    # make proper match subtype
    def categorize_correct_splice_match(self, read_intron_profile, isoform_id):
        if self.is_fsm(read_intron_profile, isoform_id):
            logger.debug("+ + Full splice match " + isoform_id)
            isoform_match = IsoformMatch(MatchClassification.fsm, self.get_gene_id(isoform_id), isoform_id)
        else:
            logger.debug("+ + Incomplete splice match " + isoform_id)
            isoform_match = IsoformMatch(MatchClassification.ism, self.get_gene_id(isoform_id), isoform_id,
                                         self.detect_ism_subtype(read_intron_profile, isoform_id))
        return isoform_match

    # make proper match subtype
    def categorize_multiple_splice_matches(self, read_intron_profile, isoform_ids):
        isoform_matches = []
        for isoform_id in isoform_ids:
            isoform_matches.append(self.categorize_correct_splice_match(read_intron_profile, isoform_id))
        return isoform_matches

    # =========== END SUPPORT ============


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
            if all(el == 0 for el in read_intron_profile.gene_profile):
                logger.debug("EMPTY - intergenic")
                return ReadAssignment(read_id, ReadAssignmentType.empty, IsoformMatch(MatchClassification.intergenic))
            else:
                logger.debug("EMPTY - intronic")
                return ReadAssignment(read_id, ReadAssignmentType.empty, IsoformMatch(MatchClassification.genic_intron))

        elif any(el == -1 for el in read_intron_profile.read_profile) or any(el == -1 for el in read_split_exon_profile.read_profile):
            # Read has missing exons / introns
            logger.debug("+ Has contradictory features")
            return self.match_contradictory(read_id, combined_read_profile)

        elif any(el == 0 for el in read_intron_profile.read_profile) or any(el == 0 for el in read_split_exon_profile.read_profile):
            # Read has extra flanking exons / introns
            logger.debug("+ Has extra flanking features")
            return self.match_with_extra_flanking(read_id, combined_read_profile)

        else:
            logger.debug("+ No contradictory features")
            assignment = self.match_non_contradictory(read_id, combined_read_profile)

            if assignment is None:
                # alternative isoforms made of known introns/exons or intron retention
                logger.debug("+ + Resolving unmatched ")
                assignment = self.match_contradictory(read_id, combined_read_profile)
            else:
                # check for extra flanking sequences
                self.check_for_extra_terminal_seqs(read_split_exon_profile, assignment)

            return assignment

    def match_non_contradictory(self, read_id, combined_read_profile):
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
            # unique intersection
            isoform_id = list(both_mathched_isoforms)[0]

            logger.debug("+ + UNIQUE match found " + isoform_id)
            logger.debug("Exon profile: " + str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
            logger.debug("Intron profile: " + str(self.gene_info.intron_profiles.profiles[isoform_id]))

            isoform_match = self.categorize_correct_splice_match(read_intron_profile, isoform_id)
            read_assignment = ReadAssignment(read_id, ReadAssignmentType.unique, isoform_match)

        elif len(both_mathched_isoforms) > 1:
            read_assignment = self.resolve_multiple_assignments(read_id, combined_read_profile, both_mathched_isoforms)

        elif len(both_mathched_isoforms) == 0:
            logger.debug("+ + Intersection is empty ")
            read_assignment = self.resolve_incompatible_assignments(read_id, combined_read_profile, intron_matched_isoforms)

        return read_assignment

    # resolve when intersection of intron and exon matched isoforms is empty
    def resolve_incompatible_assignments(self, read_id, combined_read_profile, intron_matched_isoforms):
        read_intron_profile = combined_read_profile.read_intron_profile
        read_split_exon_profile = combined_read_profile.read_split_exon_profile

        # return None for contradictory
        read_assignment = None
        if len(intron_matched_isoforms) > 0:
            # intron profile match, but extra exons
            counter = 0
            for isoform_id in intron_matched_isoforms:
                logger.debug("+ + AMBIGUOUS match " + str(counter) + ": " + isoform_id)
                counter += 1
                logger.debug("Exon profile: " + str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
                logger.debug("Intron profile: " + str(self.gene_info.intron_profiles.profiles[isoform_id]))

            if len(intron_matched_isoforms) == 1:
                isoform_id = list(intron_matched_isoforms)[0]
                logger.debug("+ + There is unique intron profile match " + isoform_id)
                isoform_match = self.categorize_correct_splice_match(read_intron_profile, list(intron_matched_isoforms)[0])
                read_assignment = ReadAssignment(read_id, ReadAssignmentType.minor, isoform_match)
            else:
                logger.debug("+ + There is ambiguous intron profile match")
                isoform_matches = self.categorize_multiple_splice_matches(read_intron_profile, intron_matched_isoforms)
                read_assignment = ReadAssignment(read_id, ReadAssignmentType.ambiguous, isoform_matches)

        return read_assignment

    # resolve case when multiple isoforms are matched
    def resolve_multiple_assignments(self, read_id, combined_read_profile, both_mathched_isoforms):
        read_intron_profile = combined_read_profile.read_intron_profile
        read_split_exon_profile = combined_read_profile.read_split_exon_profile

        counter = 0
        for isoform_id in both_mathched_isoforms:
            logger.debug("+ + AMBIGUOUS match " + str(counter) + ": " + isoform_id)
            counter += 1
            logger.debug("Exon profile: " + str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
            logger.debug("Intron profile: " + str(self.gene_info.intron_profiles.profiles[isoform_id]))

        if not self.params.resolve_ambiguous:
            isoform_matches = self.categorize_multiple_splice_matches(read_intron_profile, both_mathched_isoforms)
            read_assignment = ReadAssignment(read_id, ReadAssignmentType.ambiguous, isoform_matches)
        else:
            logger.debug("+ + Resolving ambiguity ")
            possible_isoforms = self.gene_info.ambiguous_isoforms.intersection(both_mathched_isoforms)
            intron_matched_isoforms = \
                self.pick_best_from_ambiguous(read_intron_profile.gene_profile, self.gene_info.intron_profiles.profiles,
                                              possible_isoforms)
            exon_matched_isoforms = \
                self.pick_best_from_ambiguous(read_split_exon_profile.gene_profile,
                                              self.gene_info.split_exon_profiles.profile, possible_isoforms)

            both_mathched_isoforms = intron_matched_isoforms.intersection(exon_matched_isoforms)

            if len(both_mathched_isoforms) == 1:
                isoform_match = self.categorize_correct_splice_match(read_intron_profile, list(both_mathched_isoforms)[0])
                read_assignment = ReadAssignment(read_id, ReadAssignmentType.unique, isoform_match)
            else:
                isoform_matches = self.categorize_multiple_splice_matches(read_intron_profile, both_mathched_isoforms)
                read_assignment = ReadAssignment(read_id, ReadAssignmentType.ambiguous, isoform_matches)

        return read_assignment

    # resolve assignment ambiguities caused by identical profiles
    def pick_best_from_ambiguous(self, read_gene_profile, isoform_profiles, matched_isoforms):
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

    # resolve when there are 0s  at the ends of read profile
    def match_with_extra_flanking(self, read_id, combined_read_profile):
        read_intron_profile = combined_read_profile.read_intron_profile
        logger.debug("+ + " + str(read_intron_profile.read_profile))
        if read_intron_profile.read_profile[0] == 1 and read_intron_profile.read_profile[-1] == 1:
            logger.warning("+ + Both terminal introns present, odd case")

        assignment = self.match_non_contradictory(read_id, combined_read_profile)

        if not self.params.allow_extra_terminal_introns:
            assignment.set_assignment_type(ReadAssignmentType.contradictory)
        assignment.add_common_event(MatchEventSubtype.extra_intron_out)

        for match in assignment.isoform_matches:
            match.set_match_classification(MatchClassification.nnic)
            match.add_subclassification(MatchEventSubtype.extra_intron_out)

        if not self.params.allow_extra_terminal_introns:
            if assignment.assignment_type == ReadAssignmentType.unique or \
                assignment.assignment_type == ReadAssignmentType.minor:
                assignment.set_assignment_type(ReadAssignmentType.novel)
        else:
            if assignment.assignment_type == ReadAssignmentType.unique:
                assignment.set_assignment_type(ReadAssignmentType.minor)

        return assignment

    # resolve when there are -1s in read profile or when there are no exactly matching isoforms, but no -1s in read profiles
    def match_contradictory(self, read_id, combined_read_profile):
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

        assignment = ReadAssignment(read_id, ReadAssignmentType.contradictory)
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

            matching_events = self.compare_junctions(read_intron_profile.read_features, read_region,
                                                    isoform_introns, isoform_region)
            match_classification = MatchClassification.get_contradiction_classification_from_subtypes(matching_events)
            isoform_match = IsoformMatch(match_classification, isoform_id, self.get_gene_id(isoform_id))
            for m in matching_events:
                isoform_match.add_subclassification(m)
            assignment.add_match(isoform_match)
            logger.debug("+ + Found contradiction for " + isoform_id + ": " + " ".join(matching_events))

        new_assignment_type = None
        # Change assignment from contradictory when contradictions are minor or absent
        if all(m.all_subtypes_are_none() for m in assignment.isoform_matches):
            # No contradiction
            new_assignment_type = ReadAssignmentType.unique if len(best_isoform_ids) == 1 else ReadAssignmentType.ambiguous
        elif self.params.correct_minor_errors and all(m.all_subtypes_are_alignment_artifacts() for m in assignment.isoform_matches):
            # Only alignment artifacts
            new_assignment_type = ReadAssignmentType.minor if len(best_isoform_ids) == 1 else ReadAssignmentType.ambiguous

        if new_assignment_type is not None:
            # Revise all matches as correct
            isoform_matches = self.categorize_multiple_splice_matches(read_intron_profile, best_intron_isoform_ids)
            assignment = ReadAssignment(read_id, new_assignment_type, isoform_matches)
            self.check_for_extra_terminal_seqs(read_split_exon_profile, assignment)
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
        list of detected contradiction events
        """
        if len(read_junctions) == 0:
            if -1 in isoform_junctions:
                if any(contains(read_region, rj) for rj in isoform_junctions):
                    event = MatchEventSubtype.unspliced_intron_retention
                else:
                    event = MatchEventSubtype.unspliced_genic
            else:
                event = MatchEventSubtype.unspliced
            return [event]

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
            return [MatchEventSubtype.extra_intron_out]
        else:
            logger.debug("No contradition detected, odd case")
            return [MatchEventSubtype.none]

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
            contradiction_events.append(copy.deepcopy(event))

        return list(set(contradiction_events))

    def compare_overlapping_contradictional_regions(self, read_junctions, isoform_junctions, read_cregion, isoform_cregion):
        if read_cregion is None:
            return MatchEventSubtype.intron_retention
        elif isoform_cregion is None:
            if self.are_known_introns(read_junctions, read_cregion):
                return MatchEventSubtype.extra_intron_known
            return MatchEventSubtype.extra_intron

        read_intron_total_len = sum(
            [read_junctions[i][1] - read_junctions[i][0] for i in range(read_cregion[0], read_cregion[1] + 1)])
        isoform_intron_total_len = sum(
            [isoform_junctions[i][1] - isoform_junctions[i][0] for i in range(isoform_cregion[0], isoform_cregion[1] + 1)])
        total_intron_len_diff = abs(read_intron_total_len - isoform_intron_total_len)

        read_introns_known = self.are_known_introns(read_junctions, read_cregion)

        if read_cregion[1] == read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            return self.classify_single_intron_alternation(read_junctions, isoform_junctions, read_cregion[0],
                                                           isoform_cregion[0], total_intron_len_diff, read_introns_known)

        elif read_cregion[1] - read_cregion[0] == isoform_cregion[1] - isoform_cregion[0] and \
                total_intron_len_diff < self.params.delta:
            if read_introns_known:
                return MatchEventSubtype.mutually_exclusive_exons_known
            else:
                return MatchEventSubtype.mutually_exclusive_exons_novel

        elif read_cregion[1] == read_cregion[0] and isoform_cregion[1] > isoform_cregion[0]:
            return self.classify_skipped_exons(isoform_junctions, isoform_cregion, total_intron_len_diff, read_introns_known)

        elif read_cregion[1] > read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            if read_introns_known:
                return MatchEventSubtype.exon_gain_known
            else:
                return MatchEventSubtype.exon_gain_novel

        else:
            if read_introns_known:
                return MatchEventSubtype.alternative_structure_known
            else:
                return MatchEventSubtype.alternative_structure_novel

    def classify_skipped_exons(self, isoform_junctions, isoform_cregion,
                               total_intron_len_diff, read_introns_known):
        total_exon_len = sum([isoform_junctions[i + 1][0] - isoform_junctions[i][1] + 1
                              for i in range(isoform_cregion[0], isoform_cregion[1])])

        if total_intron_len_diff < 2 * self.params.delta and total_exon_len <= self.params.max_missed_exon_len:
            return MatchEventSubtype.exon_misallignment
        else:
            if read_introns_known:
                return MatchEventSubtype.exon_skipping_known_intron
            else:
                return MatchEventSubtype.exon_skipping_novel_intron

    def classify_single_intron_alternation(self, read_junctions, isoform_junctions, read_cpos, isoform_cpos,
                                           total_intron_len_diff, read_introns_known):
        if total_intron_len_diff <= 2 * self.params.delta:
            if read_introns_known:
                return MatchEventSubtype.intron_migration
            else:
                if abs(isoform_junctions[isoform_cpos][0] - read_junctions[read_cpos][0]) <= self.params.max_intron_shift:
                    return MatchEventSubtype.intron_shift
                else:
                    return MatchEventSubtype.intron_alternation_novel
        else:
            if abs(isoform_junctions[isoform_cpos][0] - read_junctions[read_cpos][0]) < self.params.delta:
                return MatchEventSubtype.alt_acceptor_site if read_introns_known \
                    else MatchEventSubtype.alt_acceptor_site_novel
            elif abs(isoform_junctions[isoform_cpos][1] - read_junctions[read_cpos][1]) < self.params.delta:
                return MatchEventSubtype.alt_donor_site if read_introns_known \
                    else MatchEventSubtype.alt_donor_site_novel
            else:
                return MatchEventSubtype.intron_alternation_known if read_introns_known \
                    else MatchEventSubtype.intron_alternation_novel

    # === Exon matching ==
    def assign_exons(self, combined_read_profile):
        pass