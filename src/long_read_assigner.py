############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import namedtuple
from enum import Enum, unique
from functools import partial

from .common import (
    contains,
    contains_approx,
    difference_in_present_features,
    equal_profiles_in_range,
    equal_ranges,
    extra_exon_percentage,
    has_overlapping_features,
    jaccard_similarity,
    overlap_intervals,
    overlaps,
    read_coverage_fraction,
    rindex,
)
from .isoform_assignment import (
    IsoformMatch,
    MatchEvent,
    MatchEventSubtype,
    MatchClassification,
    ReadAssignment,
    ReadAssignmentType,
    SupplementaryMatchConstants,
    event_subtype_cost,
    elongation_cost,
)
from .long_read_profiles import OverlappingFeaturesProfileConstructor
from .junction_comparator import JunctionComparator
from .polya_verification import PolyAVerifier


logger = logging.getLogger('IsoQuant')

IsoformDiff = namedtuple("IsoformDiff", ("id", "diff"))


@unique
class AmbiguityResolvingMethod(Enum):
    none = -10
    monoexon_only = 10
    monoexon_and_fsm = 20
    all = 30

    minimal_score = -0.5
    top_scored_factor = 1.5


class LongReadAssigner:
    def __init__(self, gene_info, params, quick_mode=False):
        self.gene_info = gene_info
        self.params = params
        self.quick_mode = quick_mode
        self.intron_comparator = JunctionComparator(params,
                                                    OverlappingFeaturesProfileConstructor
                                                    (self.gene_info.intron_profiles.features,
                                                     (self.gene_info.start, self.gene_info.end),
                                                     comparator=partial(equal_ranges, delta=self.params.delta)))
        self.polya_verifier = PolyAVerifier(self.gene_info, self.params)

    # ======== SUPPORT FUNCTIONS =======
    def get_gene_id(self, transcript_id):
        return self.gene_info.gene_id_map[transcript_id]

    @staticmethod
    def find_distances_to_pos(isoform_ids, transcript_function, pos):
        distances = []
        for isoform_id in isoform_ids:
            distances.append((abs(transcript_function(isoform_id) - pos), isoform_id))
        return distances

    @staticmethod
    def find_overlapping_isoforms(read_exon_split_profile, isoform_profiles, hint=None):
        isoforms = set()
        isoform_set = hint if hint is not None else isoform_profiles.profiles.keys()
        for isoform_id in isoform_set:
            profile_range = overlap_intervals(read_exon_split_profile.gene_profile_range,
                                              isoform_profiles.profile_ranges[isoform_id])
            if has_overlapping_features(isoform_profiles.profiles[isoform_id], read_exon_split_profile.gene_profile,
                                        profile_range=profile_range):
                isoforms.add(isoform_id)
        return isoforms

    def find_containing_isoforms(self, read_exon_split_profile, isoform_profiles, hint=None):
        isoforms = set()
        read_exons = read_exon_split_profile.read_features
        read_region = (read_exons[0][0], read_exons[-1][1])
        isoform_set = hint if hint is not None else isoform_profiles.keys()

        for isoform_id in isoform_set:
            if contains_approx(self.gene_info.transcript_region(isoform_id), read_region,
                               self.params.min_abs_exon_overlap):
                isoforms.add(isoform_id)
        return isoforms

    @staticmethod
    def match_profile(read_gene_profile, isoform_profiles, hint=None, diff_limit=-1, profile_range=None):
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
        isoform_set = hint if hint is not None else isoform_profiles.keys()
        for isoform_id in isoform_set:
            isoform_profile = isoform_profiles[isoform_id]
            diff = difference_in_present_features(isoform_profile, read_gene_profile, diff_limit, profile_range)
            if diff_limit == -1 or diff <= diff_limit:
                isoforms.append(IsoformDiff(isoform_id, diff))
        return sorted(isoforms, key=lambda x: x.diff)

    @staticmethod
    def match_profile_exact_in_range(read_gene_profile, isoform_profiles, profile_range, hint=None):
        isoforms = []
        isoform_set = hint if hint is not None else isoform_profiles.keys()
        for isoform_id in isoform_set:
            isoform_profile = isoform_profiles[isoform_id]
            if equal_profiles_in_range(isoform_profile, read_gene_profile, profile_range):
                isoforms.append(IsoformDiff(isoform_id, 0))
        return isoforms

    def find_matching_isoforms(self, read_profile, isoform_profiles, hint=None):
        """ match read profiles to a known isoform junction profile

        Parameters
        ----------
        read_gene_profile: list of int
        isoform_profiles: dict of str
        hint: set, optional
            potential candidates, all other profiles will be skipped

        Returns
        -------
        result: list(str)
            matching isoforms

        """
        isoforms = self.match_profile_exact_in_range(read_profile.gene_profile, isoform_profiles,
                                                     read_profile.gene_profile_range, hint)
        return [x[0] for x in isoforms]

    # select most similar isoform based on different criteria
    def select_similar_isoforms(self, combined_read_profile):
        read_split_exon_profile = combined_read_profile.read_split_exon_profile
        overlapping_isoforms = self.find_overlapping_isoforms(read_split_exon_profile,
                                                              self.gene_info.split_exon_profiles)
        if not overlapping_isoforms:
            return
        # select isoforms with non-negative nucleotide score
        significantly_overlapping_isoforms = \
            self.resolve_by_nucleotide_score(combined_read_profile, overlapping_isoforms,
                                             similarity_function=self.coverage_based_nucleotide_score,
                                             top_scored_factor=0)
        if not significantly_overlapping_isoforms:
            return

        # find intron matches
        intron_matching_isoforms = self.match_profile(combined_read_profile.read_intron_profile.gene_profile,
                                                      self.gene_info.intron_profiles.profiles,
                                                      profile_range=combined_read_profile.read_intron_profile.gene_profile_range,
                                                      hint=significantly_overlapping_isoforms)

        # add extra terminal bases as potential inconsistency event
        candidates = []
        read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[-1][1])
        for isoform_id, diff_introns in intron_matching_isoforms:
            transcript_start, transcript_end = self.gene_info.transcript_region(isoform_id)
            extra_left = 0.5 if read_region[0] + self.params.delta < transcript_start else 0
            extra_right = 0.5 if read_region[1] - self.params.delta > transcript_end else 0
            candidates.append((isoform_id, diff_introns + extra_right + extra_left))
        # select isoforms that have similar number of potential inconsistencies
        best_diff = min(candidates, key=lambda x: x[1])[1]
        best_candidates = [x[0] for x in filter(lambda x: x[1] <= best_diff + 3, candidates)]
        # logger.debug("+ + Closest matching isoforms " + str(best_candidates))

        return best_candidates

    # detect exon elongation subtyp
    def categorize_exon_elongation_subtype(self, read_split_exon_profile, isoform_id):
        split_exons = self.gene_info.split_exon_profiles.features
        isoform_profile = self.gene_info.split_exon_profiles.profiles[isoform_id]
        profile_range = self.gene_info.split_exon_profiles.profile_ranges[isoform_id]
        read_profile_range = read_split_exon_profile.gene_profile_range

        # find first and last common exons
        common_first_exon = -1
        isoform_first_exon = profile_range[0]
        for i in range(max(isoform_first_exon, read_profile_range[0]), len(split_exons)):
            if isoform_profile[i] == read_split_exon_profile.gene_profile[i] == 1:
                common_first_exon = i
                break

        common_last_exon = -1
        isoform_last_exon = profile_range[1] - 1
        for i in range(min(isoform_last_exon, read_profile_range[1] - 1), -1, -1):
            if isoform_profile[i] == read_split_exon_profile.gene_profile[i] == 1:
                common_last_exon = i
                break
        if common_first_exon == -1 or common_last_exon == -1:
            logger.warning(" + Odd case for exon elongation, no matching exons")

        events = []

        first_read_exon = read_split_exon_profile.read_features[0]
        if overlaps(first_read_exon, split_exons[common_first_exon]):
            extra_left = split_exons[common_first_exon][0] - first_read_exon[0]
            if common_first_exon == isoform_first_exon:
                # logger.debug("+ + Left: %d %d" % (common_first_exon, isoform_first_exon))
                if abs(extra_left) <= self.params.minor_exon_extension:
                    if abs(extra_left) <= self.params.delta:
                        left_match_event = MatchEventSubtype.terminal_site_match_left_precise
                    else:
                        left_match_event = MatchEventSubtype.terminal_site_match_left
                    events.append(MatchEvent(left_match_event, event_info=extra_left))

                if extra_left > self.params.minor_exon_extension:
                    events.append(MatchEvent(MatchEventSubtype.major_exon_elongation_left, event_info=extra_left))
                elif extra_left > self.params.delta:
                    events.append(MatchEvent(MatchEventSubtype.exon_elongation_left, event_info=extra_left))
            elif self.params.minor_exon_extension >= extra_left > self.params.delta:
                # minor exon elongations detected for non-terminal exons too
                events.append(MatchEvent(MatchEventSubtype.exon_elongation_left, event_info=extra_left))

        last_read_exon = read_split_exon_profile.read_features[-1]
        if overlaps(last_read_exon, split_exons[common_last_exon]):
            extra_right = last_read_exon[1] - split_exons[common_last_exon][1]
            if common_last_exon == isoform_last_exon:
                # logger.debug("+ + Right: %d %d" % (common_last_exon, isoform_last_exon))
                if  abs(extra_right) <= self.params.minor_exon_extension:
                    if abs(extra_right) <= self.params.delta:
                        right_match_event = MatchEventSubtype.terminal_site_match_right_precise
                    else:
                        right_match_event = MatchEventSubtype.terminal_site_match_right
                    events.append(MatchEvent(right_match_event, event_info=extra_right))

                if extra_right > self.params.minor_exon_extension:
                    events.append(MatchEvent(MatchEventSubtype.major_exon_elongation_right, event_info=extra_right))
                elif extra_right > self.params.delta:
                    events.append(MatchEvent(MatchEventSubtype.exon_elongation_right, event_info=extra_right))
            elif self.params.minor_exon_extension >= extra_right > self.params.delta:
                events.append(MatchEvent(MatchEventSubtype.exon_elongation_right, event_info=extra_right))

        return events

    # select best assignment based on nucleotide similarity
    # score = similarity_function(read_exons, isoform_exons)
    def resolve_by_nucleotide_score(self, combined_read_profile, matched_isoforms, similarity_function,
                                    min_similarity=AmbiguityResolvingMethod.minimal_score.value,
                                    top_scored_factor=AmbiguityResolvingMethod.top_scored_factor.value):
        if not matched_isoforms:
            return []

        # logger.debug("+ + + Resolving by nucleotide similarity")
        read_exons = combined_read_profile.read_split_exon_profile.read_features

        scores = []
        for isoform_id in matched_isoforms:
            isoform_exons = self.gene_info.all_isoforms_exons[isoform_id]
            scores.append((isoform_id, similarity_function(read_exons, isoform_exons)))

        #logger.debug("Scores: " + str(scores))
        scores = sorted(scores, key=lambda x:x[1], reverse=True)
        best_score = scores[0][1]
        if top_scored_factor == 0:
            top_scored = sorted(filter(lambda x: x[1] >= min_similarity, scores))
        else:
            top_scored = sorted(filter(lambda x: x[1] * top_scored_factor >= best_score and x[1] >= min_similarity,
                                       scores))
        # logger.debug("+ + + Best score = %f, all candidates = %s" % (best_score, str(top_scored)))
        return list(map(lambda x: x[0], top_scored))

    # score = read's covered fraction - flanking fraction (can take values in [-1,1])
    # gives priority to isoforms that contain the read
    def coverage_based_nucleotide_score(self, read_exons, isoform_exons):
        read_coverage = read_coverage_fraction(read_exons, isoform_exons)
        isoform_extended_region = (isoform_exons[0][0] - self.params.minor_exon_extension,
                                   isoform_exons[-1][1] + self.params.minor_exon_extension)
        flanking_percentage = extra_exon_percentage(isoform_extended_region, read_exons)

        return read_coverage - flanking_percentage

    # score = Jaccard similarity - flanking fraction (can take values in [-1,1])
    # gives priority to isoforms that contain the read
    def jaccard_based_nucleotide_score(self, read_exons, isoform_exons):
        js = jaccard_similarity(read_exons, isoform_exons)
        isoform_extended_region = (isoform_exons[0][0] - self.params.minor_exon_extension,
                                   isoform_exons[-1][1] + self.params.minor_exon_extension)
        flanking_percentage = extra_exon_percentage(isoform_extended_region, read_exons)

        return js - flanking_percentage

    # ====== CLASSIFICATION =======

    # check for extra sequences and modify assignment accordingly
    def check_read_ends(self, read_split_exon_profile, assignment):
        for match in assignment.isoform_matches:
            if match.assigned_transcript is None:
                continue

            exon_elongation_types = self.categorize_exon_elongation_subtype(read_split_exon_profile,
                                                                            match.assigned_transcript)
            for e in exon_elongation_types:
                match.add_subclassification(e)

            if any(MatchEventSubtype.is_major_elongation(e.event_type) for e in exon_elongation_types):
                # serious exon elongation
                if not assignment.assignment_type.is_inconsistent():
                    assignment.set_assignment_type(ReadAssignmentType.inconsistent_non_intronic)
            elif any(MatchEventSubtype.is_minor_elongation(e.event_type) for e in exon_elongation_types):
                # minor exon elongation
                if assignment.assignment_type == ReadAssignmentType.unique:
                    assignment.set_assignment_type(ReadAssignmentType.unique_minor_difference)

    # get incompleteness type, assuming introns are matched and read is spliced
    def detect_ism_subtype(self, read_region, isoform_id):
        isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
        isoform_intron_region = (isoform_introns[0][0], isoform_introns[-1][1])
        is_left_truncated = isoform_intron_region[0] < read_region[0]
        is_right_truncated = isoform_intron_region[1] > read_region[1]

        if is_left_truncated and is_right_truncated:
            # logger.debug(" + Internal")
            event_type = MatchEventSubtype.ism_internal
        elif is_left_truncated:
            # logger.debug(" + Truncated on the left")
            event_type = MatchEventSubtype.ism_left
        elif is_right_truncated:
            # logger.debug(" + Truncated on the right")
            event_type = MatchEventSubtype.ism_right
        else:
            # logger.debug(" + No ISM truncation ")
            event_type = MatchEventSubtype.none
        return MatchEvent(event_type)

    # check where it is full splice match, assume all introns are matched and read is spliced
    def is_fsm(self, read_region, isoform_id):
        isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
        # all isoform instrons are within read
        return contains(read_region, (isoform_introns[0][0], isoform_introns[-1][1]))

    # make proper match subtype, assume all introns are matched or read is unspliced
    def categorize_correct_splice_match(self, combined_read_profile, isoform_id):
        read_exons = combined_read_profile.read_split_exon_profile.read_features
        read_region = (read_exons[0][0], read_exons[-1][1])

        if len(combined_read_profile.read_intron_profile.read_profile) == 0 or \
                len(self.gene_info.all_isoforms_introns[isoform_id]) == 0:
            # logger.debug(" + Mono exon")
            match_classification = MatchClassification.mono_exon_match
            match_subclassifications = MatchEvent(MatchEventSubtype.mono_exon_match)
        elif self.is_fsm(read_region, isoform_id):
            # logger.debug("+ + Full splice match " + isoform_id)
            match_classification = MatchClassification.full_splice_match
            match_subclassifications = MatchEvent(MatchEventSubtype.fsm)
        else:
            # logger.debug("+ + Incomplete splice match " + isoform_id)
            match_classification = MatchClassification.incomplete_splice_match
            match_subclassifications = self.detect_ism_subtype(read_region, isoform_id)

        return IsoformMatch(match_classification, self.get_gene_id(isoform_id), isoform_id,
                            match_subclassifications, self.gene_info.isoform_strands[isoform_id])

    # make proper match subtype
    def categorize_multiple_splice_matches(self, combined_read_profile, isoform_ids):
        isoform_matches = []
        for isoform_id in sorted(isoform_ids):
            isoform_matches.append(self.categorize_correct_splice_match(combined_read_profile, isoform_id))
        return isoform_matches

    def categorize_correct_unspliced_match(self, combined_read_profile, isoform_id):
        isoform_exons = self.gene_info.all_isoforms_exons[isoform_id]
        if len(isoform_exons) == 1:
            # both isoform and read are monoexon
            events = [MatchEvent(MatchEventSubtype.mono_exon_match)]
        else:
            # monoexonic read is inside exon
            events = [MatchEvent(MatchEventSubtype.mono_exonic)]

        match_classification = MatchClassification.get_mono_exon_classification(events)
        return IsoformMatch(match_classification, self.get_gene_id(isoform_id), isoform_id, events,
                            self.gene_info.isoform_strands[isoform_id])

    def categorize_multiple_unspliced_matches(self, combined_read_profile, isoform_ids):
        isoform_matches = []
        for isoform_id in sorted(isoform_ids):
            isoform_matches.append(self.categorize_correct_unspliced_match(combined_read_profile, isoform_id))
        return isoform_matches

    # =========== MAIN PART ============

    # === Isoform matching function ===
    def assign_to_isoform(self, read_id, combined_read_profile):
        """ assign read to isoform according to it

        Parameters
        ----------
        read_id: str
        combined_read_profile: CombinedProfile

        Returns
        -------

        """
        read_intron_profile = combined_read_profile.read_intron_profile
        read_split_exon_profile = combined_read_profile.read_split_exon_profile

        # logger.debug("> Read blocks" + str(read_split_exon_profile.read_features))
        # logger.debug("Gene split exons" + str(self.gene_info.split_exon_profiles.features))
        # logger.debug("Read exons profile" + str(read_split_exon_profile.read_profile))
        # logger.debug("Gene exons profile" + str(read_split_exon_profile.gene_profile))

        # logger.debug("< Read introns" + str(read_intron_profile.read_features))
        # logger.debug("Gene introns" + str(self.gene_info.intron_profiles.features))
        # logger.debug("Read intron profile" + str(read_intron_profile.read_profile))
        # logger.debug("Gene intron profile" + str(read_intron_profile.gene_profile))
        if self.gene_info.empty():
            return ReadAssignment(read_id, ReadAssignmentType.intergenic,
                                  IsoformMatch(MatchClassification.intergenic))

        if all(el != 1 for el in read_split_exon_profile.read_profile) \
                or all(el == 0 or el == -2 for el in read_split_exon_profile.gene_profile):
            read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[-1][1])
            gene_region = (self.gene_info.split_exon_profiles.features[0][0],
                              self.gene_info.split_exon_profiles.features[-1][1])
            # none of the blocks matched
            if not overlaps(read_region, gene_region):
                # logger.debug("EMPTY - noninformative")
                assignment = ReadAssignment(read_id, ReadAssignmentType.noninformative, IsoformMatch(MatchClassification.intergenic))
            elif all(el != 1 for el in read_split_exon_profile.gene_profile):
                # logger.debug("EMPTY - intronic")
                assignment = ReadAssignment(read_id, ReadAssignmentType.noninformative, IsoformMatch(MatchClassification.genic_intron))
            else:
                # logger.debug("EMPTY - genic")
                assignment = ReadAssignment(read_id, ReadAssignmentType.noninformative, IsoformMatch(MatchClassification.genic))
            return assignment

        elif any(el == -1 for el in read_intron_profile.read_profile) \
                or any(el == -1 for el in read_split_exon_profile.read_profile):
            # Read has inconsistent exons / introns
            # logger.debug("+ Has inconsistent features (novel intron / exons)")
            assignment = self.match_inconsistent(read_id, combined_read_profile)

        elif any(el == 0 for el in read_intron_profile.read_profile) \
                or any(el == 0 for el in read_split_exon_profile.read_profile):
            # Read has extra flanking exons / introns
            # logger.debug("+ Has extra flanking features")
            assignment = self.match_inconsistent(read_id, combined_read_profile)

        else:
            # logger.debug("+ No contradictory features in read, but inconsistent isoforms may still be detected")
            assignment = self.match_consistent(read_id, combined_read_profile)

            if assignment is None:
                # alternative isoforms made of known introns/exons or intron retention
                # logger.debug("+ + Resolving unmatched ")
                assignment = self.match_inconsistent(read_id, combined_read_profile)

        return assignment

    def match_consistent(self, read_id, combined_read_profile):
        """ match profile when all read features are assigned and read has at least one intron

        Parameters
        ----------
        read_id: str
        combined_read_profile: CombinedReadProfiles

        Returns
        -------

        """
        read_intron_profile = combined_read_profile.read_intron_profile
        read_split_exon_profile = combined_read_profile.read_split_exon_profile
        isoform_split_exon_profiles = self.gene_info.split_exon_profiles.profiles

        containing_isoforms = self.find_containing_isoforms(read_split_exon_profile, isoform_split_exon_profiles)
        if not containing_isoforms:
            return None

        overlapping_isoforms = self.find_overlapping_isoforms(read_split_exon_profile, self.gene_info.split_exon_profiles,
                                                              hint=containing_isoforms)
        if not overlapping_isoforms:
            return None

        consistent_isoforms = self.find_matching_isoforms(read_intron_profile,
                                                          self.gene_info.intron_profiles.profiles,
                                                          hint=overlapping_isoforms)

        if not read_intron_profile.read_profile:
            assignment = self.match_consistent_unspliced(read_id, combined_read_profile, consistent_isoforms)
        else:
            assignment = self.match_consistent_spliced(read_id, combined_read_profile, consistent_isoforms)

        if assignment is not None:
            self.check_read_ends(combined_read_profile.read_split_exon_profile, assignment)
            self.verify_read_ends_for_assignment(combined_read_profile, assignment)

            # if apa detected isoform should be treated as inconsistent and rechecked
            if assignment.assignment_type.is_inconsistent():
                return None

        return assignment

    def match_consistent_spliced(self, read_id, combined_read_profile, consistent_isoforms):
        isoform_split_exon_profiles = self.gene_info.split_exon_profiles.profiles
        matched_isoforms = consistent_isoforms
        if len(consistent_isoforms) > 1:
            exon_matched_isoforms = \
                self.find_matching_isoforms(combined_read_profile.read_split_exon_profile,
                                            isoform_split_exon_profiles, hint=consistent_isoforms)
            if len(exon_matched_isoforms) != 0:
                matched_isoforms = exon_matched_isoforms

        if len(matched_isoforms) > 1:
            read_exons = combined_read_profile.read_split_exon_profile.read_features
            read_region = (read_exons[0][0], read_exons[-1][1])
            if self.params.resolve_ambiguous == AmbiguityResolvingMethod.all or \
               (self.params.resolve_ambiguous == AmbiguityResolvingMethod.monoexon_and_fsm and
                any(self.is_fsm(read_region, isoform_id) for isoform_id in matched_isoforms)):
                # resolve spliced using nucleotide score
                matched_isoforms = \
                    self.resolve_by_nucleotide_score(combined_read_profile, matched_isoforms,
                                                     similarity_function=self.jaccard_based_nucleotide_score)

        read_assignment = None
        if len(matched_isoforms) == 1:
            isoform_id = list(matched_isoforms)[0]
            # logger.debug("+ + UNIQUE intron match found " + isoform_id)
            isoform_match = self.categorize_correct_splice_match(combined_read_profile, isoform_id)
            read_assignment = ReadAssignment(read_id, ReadAssignmentType.unique, isoform_match)

        elif len(matched_isoforms) > 1:
            # logger.debug("+ + Ambiguous read")
            isoform_matches = self.categorize_multiple_splice_matches(combined_read_profile, matched_isoforms)
            return ReadAssignment(read_id, ReadAssignmentType.ambiguous, isoform_matches)

        return read_assignment

    def match_consistent_unspliced(self, read_id, combined_read_profile, consistent_isoforms):
        # logger.debug("+  Resolving monoexonic read")

        matched_isoforms = consistent_isoforms
        if len(consistent_isoforms) > 1  and self.params.resolve_ambiguous != AmbiguityResolvingMethod.none:
            matched_isoforms = self.resolve_by_nucleotide_score(combined_read_profile, sorted(consistent_isoforms),
                                                                similarity_function=self.jaccard_based_nucleotide_score)

        read_assignment = None
        if len(matched_isoforms) == 1:
            isoform_id = matched_isoforms[0]
            # logger.debug("Single monoexonic match: %s" % isoform_id)
            isoform_match = self.categorize_correct_unspliced_match(combined_read_profile, isoform_id)
            read_assignment = ReadAssignment(read_id, ReadAssignmentType.unique, isoform_match)

        elif len(matched_isoforms) > 1:
            # logger.debug("Nucleotide similarity picked multiple isoforms")
            isoform_matches = self.categorize_multiple_unspliced_matches(combined_read_profile, matched_isoforms)
            read_assignment = ReadAssignment(read_id, ReadAssignmentType.ambiguous, isoform_matches)

        return read_assignment

    # resolve when there are -1s in read profile or when there are no exactly matching isoforms, but no -1s in read profiles
    def match_inconsistent(self, read_id, combined_read_profile):
        if self.quick_mode:
            return ReadAssignment(read_id, ReadAssignmentType.inconsistent, IsoformMatch(MatchClassification.genic))

        # select most similar isoforms based on multiple criteria
        best_candidates = self.select_similar_isoforms(combined_read_profile)
        if not best_candidates:
            return ReadAssignment(read_id, ReadAssignmentType.noninformative, IsoformMatch(MatchClassification.genic))

        # logger.debug("* Best candidates for inconsistency detection: " + str(best_candidates))
        # detect inconsistency for each one
        read_matches = self.detect_inconsistensies(read_id, combined_read_profile, best_candidates)
        if not read_matches:
            return ReadAssignment(read_id, ReadAssignmentType.noninformative)
        # logger.debug("* Inconsistencies detected: " + str(read_matches))

        # select ones with the least number of inconsistent events
        best_isoforms, penalty_score = self.select_best_among_inconsistent(combined_read_profile, read_matches)
        if not best_isoforms:
            return ReadAssignment(read_id, ReadAssignmentType.noninformative)
        # logger.debug("* Selected isoforms: " + str(best_isoforms))

        best_isoforms = sorted(best_isoforms)
        assignment_type = self.classify_assignment(best_isoforms, read_matches)
        if not combined_read_profile.read_intron_profile.read_profile:
            isoform_matches = self.create_monoexon_matches(read_matches, best_isoforms)
        elif assignment_type.is_inconsistent():
            isoform_matches = self.create_inconsistent_matches(read_matches, best_isoforms, penalty_score)
        else:
            isoform_matches = self.create_consistent_matches(read_matches, best_isoforms, combined_read_profile)
        return ReadAssignment(read_id, assignment_type, isoform_matches)

    def classify_assignment(self, best_isoforms, read_matches):
        is_abmiguous = len(best_isoforms) > 1
        all_event_types = set()
        for isoform_id in best_isoforms:
            for e in read_matches[isoform_id]:
                all_event_types.add(e.event_type)

        #logger.debug("* All events: " + str(all_event_types))
        if all(MatchEventSubtype.is_consistent(e) for e in all_event_types):
            #logger.debug("* * Assignment seems to be consistent")
            assignment_type = ReadAssignmentType.ambiguous if is_abmiguous else ReadAssignmentType.unique
        elif any(MatchEventSubtype.is_major_inconsistency(e) for e in all_event_types):
            #logger.debug("* * Assignment is inconsistent")
            if any(MatchEventSubtype.is_intronic_inconsistency(e) for e in all_event_types):
                inconsistency_type = ReadAssignmentType.inconsistent
            else:
                inconsistency_type = ReadAssignmentType.inconsistent_non_intronic
            assignment_type = ReadAssignmentType.inconsistent_ambiguous if is_abmiguous else inconsistency_type
        elif any(MatchEventSubtype.is_minor_error(e) for e in all_event_types):
            #logger.debug("* * Assignment has minor errors")
            assignment_type = ReadAssignmentType.ambiguous if is_abmiguous else ReadAssignmentType.unique_minor_difference
        else:
            logger.warning("Unexpected event reported: " + str(all_event_types))
            assignment_type = ReadAssignmentType.noninformative

        return assignment_type

    def create_monoexon_matches(self, read_matches, selected_isoforms):
        matches = []
        for isoform_id in selected_isoforms:
            read_match = read_matches[isoform_id]
            match_classification = MatchClassification.get_mono_exon_classification(read_match)
            isoform_match = IsoformMatch(match_classification, self.get_gene_id(isoform_id), isoform_id,
                                         read_match, self.gene_info.isoform_strands[isoform_id])
            matches.append(isoform_match)
        return matches

    def create_consistent_matches(self, read_matches, selected_isoforms, combined_read_profile):
        matches = []
        for isoform_id in selected_isoforms:
            isoform_match = self.categorize_correct_splice_match(combined_read_profile, isoform_id)
            for e in read_matches[isoform_id]:
                if e.event_type != MatchEventSubtype.none:
                    isoform_match.add_subclassification(e)
            matches.append(isoform_match)
        return matches

    def create_inconsistent_matches(self, read_matches, selected_isoforms, penalty_score=0.0):
        matches = []
        for isoform_id in selected_isoforms:
            read_match = read_matches[isoform_id]
            match_classification = MatchClassification.get_inconsistency_classification(read_match)
            isoform_match = IsoformMatch(match_classification, self.get_gene_id(isoform_id), isoform_id,
                                         read_match, self.gene_info.isoform_strands[isoform_id],
                                         penalty_score=penalty_score)
            matches.append(isoform_match)
        return matches

    # detect inconsistensies for selected isoforms
    def detect_inconsistensies(self, read_id, combined_read_profile, matched_isoforms):
        """ compare read to closest matching isoforms

        Parameters
        ----------
        read_id: str
        matched_isoforms: list of str

        Returns
        -------
        result: list of tuples (isoform_id, list of events)
        """
        # get isoforms that have closes intron and exon profiles
        #logger.debug("Detecting difference for %s, %d matched isoforms" % (read_id, len(matched_isoforms)))
        read_intron_profile = combined_read_profile.read_intron_profile
        read_split_exon_profile = combined_read_profile.read_split_exon_profile

        #logger.debug("+ + Closest matching isoforms " + str(matched_isoforms))
        # isoform_id = best_isoform_ids[0]
        # logger.debug(str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
        # logger.debug(str(self.gene_info.intron_profiles.profiles[isoform_id]))

        read_matches = {}
        for isoform_id in sorted(matched_isoforms):
            # logger.debug("Checking isoform %s" % isoform_id)
            # get intron coordinates
            isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
            # read start-end coordinates
            read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[-1][1])
            # isoform start-end
            isoform_region = self.gene_info.transcript_region(isoform_id)
            # logger.debug("R: " + str(combined_read_profile.read_split_exon_profile.read_features))
            # logger.debug(str(combined_read_profile.read_intron_profile.read_features))
            # logger.debug("I: " + str(self.gene_info.all_isoforms_exons[isoform_id]))

            matching_events = self.intron_comparator.compare_junctions(read_intron_profile.read_features, read_region,
                                                                       isoform_introns, isoform_region)

            if len(matching_events) == 1 and matching_events[0].event_type == MatchEventSubtype.undefined:
                continue

            matching_events += self.categorize_exon_elongation_subtype(read_split_exon_profile, isoform_id)
            matching_events = self.polya_verifier.verify_read_ends(combined_read_profile, isoform_id, matching_events)
            read_matches[isoform_id] = matching_events

            # logger.debug("+ + Found contradiction for " + isoform_id + ": " + " ".join(map(lambda x: x.event_type.name, matching_events)))
        return read_matches

    def select_best_among_inconsistent(self, combined_read_profile, read_matches):
        isoform_scores = []
        for isoform_id, match_events in read_matches.items():
            # logger.debug("* * Scoring inconsistencies for " + isoform_id + ": " + str(match_events))
            penalty_score = 0.0
            for e in match_events:
                event_count = 1
                if e.isoform_region != SupplementaryMatchConstants.undefined_region and \
                        SupplementaryMatchConstants.absent_position not in e.isoform_region and \
                        e.event_type in {MatchEventSubtype.exon_skipping_known, MatchEventSubtype.exon_skipping_novel}:
                    event_count = e.isoform_region[1] - e.isoform_region[0] + 1
                elif e.read_region != SupplementaryMatchConstants.undefined_region and \
                        SupplementaryMatchConstants.absent_position not in e.read_region:
                    event_count = e.read_region[1] - e.read_region[0] + 1
                    if e.event_type in {MatchEventSubtype.exon_gain_novel,
                                        MatchEventSubtype.exon_gain_known,
                                        MatchEventSubtype.mutually_exclusive_exons_novel,
                                        MatchEventSubtype.mutually_exclusive_exons_known,
                                        MatchEventSubtype.exon_detach_known,
                                        MatchEventSubtype.exon_detach_novel}:
                        event_count -= 1
                    elif e.event_type in {MatchEventSubtype.intron_retention,
                                          MatchEventSubtype.unspliced_intron_retention,
                                          MatchEventSubtype.fake_micro_intron_retention,
                                          MatchEventSubtype.incomplete_intron_retention_left,
                                          MatchEventSubtype.incomplete_intron_retention_right}:
                        event_count = 1

                event_cost = event_subtype_cost[e.event_type]
                if e.event_type in {MatchEventSubtype.major_exon_elongation_left,
                                    MatchEventSubtype.major_exon_elongation_right,
                                    MatchEventSubtype.exon_elongation_right,
                                    MatchEventSubtype.exon_elongation_left}:
                    event_cost = elongation_cost(self.params, e.event_info)

                penalty_score +=  event_cost * event_count
                # logger.debug("* * * Event " + str(e.event_type) + ", introns affected " + str(event_count) +
                #             ", event cost " + str(event_cost) +
                #             ". Updated penalty_score: " + str(penalty_score))
            # logger.debug("* * Final penalty_score for isoform " + isoform_id + ": " + str(penalty_score))
            isoform_scores.append((isoform_id, penalty_score))

        min_penalty_score = min(isoform_scores, key=lambda x:x[1])[1]
        # logger.debug("* * Best penalty_score " + str(min_penalty_score))
        best_isoforms = [x[0] for x in filter(lambda x:x[1] == min_penalty_score, isoform_scores)]
        # logger.debug("* * Best isoforms " + str(best_isoforms))

        # if several isoforms are tied select the best according to nucl penalty_score
        if len(best_isoforms) > 1:
            best_isoforms = self.resolve_by_nucleotide_score(combined_read_profile, best_isoforms,
                                                             similarity_function=self.coverage_based_nucleotide_score)
            # logger.debug("* * Resolved by nucl penalty_score " + str(best_isoforms))

        return best_isoforms, min_penalty_score

    # ==== POLYA STUFF ====
    def verify_read_ends_for_assignment(self, combined_read_profile, assignment):
        match_dict = {}
        for match in assignment.isoform_matches:
            if match.assigned_transcript is None:
                continue
            match.match_subclassifications = \
                self.polya_verifier.verify_read_ends(combined_read_profile, match.assigned_transcript,
                                                     match.match_subclassifications)
            match_dict[match.assigned_transcript] = match.match_subclassifications

        assignment.assignment_type = self.classify_assignment(match_dict.keys(), match_dict)



