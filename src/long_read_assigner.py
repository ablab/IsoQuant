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
from src.junction_comparator import *

logger = logging.getLogger('IsoQuant')

IsoformDiff = namedtuple("IsoformDiff", ("id", "diff"))


class ExonAmbiguityResolvingMethod(Enum):
    mono_exonic_only = 1
    full_splice_matches_only = 2
    all = 3

    minimal_jaccard_similarity = 0.3
    top_scored_factor = 1.2


class LongReadAssigner:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params
        self.intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.intron_profiles.features,
                                                  (self.gene_info.start, self.gene_info.end),
                                                  comparator=partial(equal_ranges, delta=self.params.delta))
        self.intron_comparator = JunctionComparator(params, self.intron_profile_constructor)

    # ======== SUPPORT FUNCTIONS =======
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

    def get_gene_id(self, transcript_id):
        return self.gene_info.gene_id_map[transcript_id]

    # check for extra sequences and modify assignment accordingly
    def check_for_extra_terminal_seqs(self, read_split_exon_profile, assignment):
        for match in assignment.isoform_matches:
            exon_elongation_type = self.categorize_exon_elongation_subtype(read_split_exon_profile,
                                                                           match.assigned_transcript)
            if exon_elongation_type == MatchEventSubtype.extra_exon_out:
                # serious exon elongation
                match.add_subclassification(make_event(exon_elongation_type))
                if assignment.assignment_type in (ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference):
                    assignment.set_assignment_type(ReadAssignmentType.contradictory)
            elif exon_elongation_type != MatchEventSubtype.none:
                # minor exon elongation
                match.add_subclassification(make_event(exon_elongation_type))
                if assignment.assignment_type == ReadAssignmentType.unique:
                    assignment.set_assignment_type(ReadAssignmentType.unique_minor_difference)

    # detect exon elongation subtyp
    def categorize_exon_elongation_subtype(self, read_split_exon_profile, isoform_id):
        split_exons = self.gene_info.split_exon_profiles.features
        isoform_profile = self.gene_info.split_exon_profiles.profiles[isoform_id]

        # find first and last common exons
        common_first_exon = -1
        for i in range(len(split_exons)):
            if isoform_profile[i] == 1 and read_split_exon_profile.gene_profile[i] == 1:
                common_first_exon = i
                break
        common_last_exon = -1
        for i in range(len(split_exons)):
            index = len(split_exons) - i - 1
            if isoform_profile[index] == 1 and read_split_exon_profile.gene_profile[index] == 1:
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
        logger.debug("+ + Read: " + str(read_start) + "-" + str(read_end))
        logger.debug("+ + Isoform: " + str(isoform_start) + "-" + str(isoform_end))
        logger.debug("+ + Extra bases: left = %d, right = %d" % (extra_left, extra_right))

        if extra_right <= self.params.delta and extra_left <= self.params.delta:
            logger.debug("+ + + None")
            return MatchEventSubtype.none
        else:
            if extra_right < self.params.max_exon_extension and extra_left < self.params.max_exon_extension:
                logger.debug("+ + + Minor")
                if extra_right > self.params.delta and extra_left > self.params.delta:
                    logger.debug("+ + + Extra sequence on both ends")
                    return MatchEventSubtype.exon_elongation_both
                elif extra_right > self.params.delta:
                    logger.debug("+ + + Extra sequence on right end")
                    if self.gene_info.isoform_strands[isoform_id] == "+":
                        return MatchEventSubtype.exon_elongation3
                    else:
                        return MatchEventSubtype.exon_elongation5
                else:
                    logger.debug("+ + + Extra sequence on left end")
                    if self.gene_info.isoform_strands[isoform_id] == "-":
                        return MatchEventSubtype.exon_elongation3
                    else:
                        return MatchEventSubtype.exon_elongation5
            else:
                logger.debug("+ + + Major")
                return MatchEventSubtype.extra_exon_out

    # get incompleteness type
    def detect_ism_subtype(self, read_intron_profile, isoform_id):
        if len(read_intron_profile.read_profile) == 0:
            logger.debug(" + Mono exon")
            return MatchEventSubtype.mono_exonic

        read_profile = read_intron_profile.gene_profile
        isoform_profile = self.gene_info.intron_profiles.profiles[isoform_id]
        is_left_truncated = left_truncated(read_profile, isoform_profile)
        is_right_truncated = right_truncated(read_profile, isoform_profile)

        if is_left_truncated and is_right_truncated:
            logger.debug(" + Internal")
            event_type = MatchEventSubtype.ism_internal
        elif is_left_truncated:
            logger.debug(" + Truncated on the left")
            if self.gene_info.isoform_strands[isoform_id] == "+":
                event_type = MatchEventSubtype.ism_5
            else:
                event_type = MatchEventSubtype.ism_3
        elif is_right_truncated:
            logger.debug(" + Truncated on the right")
            if self.gene_info.isoform_strands[isoform_id] == "-":
                event_type = MatchEventSubtype.ism_5
            else:
                event_type = MatchEventSubtype.ism_3
        else:
            logger.debug(" + No ISM truncation ")
            event_type = MatchEventSubtype.none
        return make_event(event_type)

    # check where it is full splice match
    def is_fsm(self, read_intron_profile, isoform_id):
        # TODO include check for minor alignment errors
        read_profile = read_intron_profile.gene_profile
        if not read_intron_profile.read_profile or not all(el == 1 for el in read_intron_profile.read_profile):
            return False

        isoform_profile = self.gene_info.intron_profiles.profiles[isoform_id]
        return not left_truncated(read_profile, isoform_profile) \
               and not right_truncated(read_profile, isoform_profile)

    def find_overlapping_isoforms(self, read_exon_split_profile, isoform_profiles):
        isoforms = set()
        for isoform_id, isoform_profile in isoform_profiles.items():
            if has_overlapping_features(isoform_profile, read_exon_split_profile):
                isoforms.add(isoform_id)
        return isoforms

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
        return sorted(isoforms, key=lambda x: x.diff)

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
            isoform_match = IsoformMatch(MatchClassification.full_splice_match, self.get_gene_id(isoform_id),
                                         isoform_id, make_event(MatchEventSubtype.fsm))
        else:
            logger.debug("+ + Incomplete splice match " + isoform_id)
            isoform_match = IsoformMatch(MatchClassification.incomplete_splice_match, self.get_gene_id(isoform_id),
                                         isoform_id, self.detect_ism_subtype(read_intron_profile, isoform_id))
        return isoform_match

    # make proper match subtype
    def categorize_multiple_splice_matches(self, read_intron_profile, isoform_ids):
        isoform_matches = []
        for isoform_id in isoform_ids:
            isoform_matches.append(self.categorize_correct_splice_match(read_intron_profile, isoform_id))
        return isoform_matches

    def categorize_mono_exonic_read_match(self, combined_read_profile, isoform_id):
        read_split_exon_profile = combined_read_profile.read_split_exon_profile
        read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[-1][1])
        isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
        matching_event = self.intron_comparator.get_mono_exon_subtype(read_region, isoform_introns)
        return IsoformMatch(MatchClassification.get_mono_exon_classification_from_subtypes(matching_event),
                            self.get_gene_id(isoform_id), isoform_id, matching_event)

    def categorize_multiple_mono_exon_matches(self, combined_read_profile, isoform_ids):
        isoform_matches = []
        for isoform_id in isoform_ids:
            isoform_matches.append(self.categorize_mono_exonic_read_match(combined_read_profile, isoform_id))
        return isoform_matches
    # =========== END SUPPORT ============

    # === Isoform matching function ===
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

        # logger.debug("> Read blocks" + str(read_split_exon_profile.read_features))
        # logger.debug("Gene split exons" + str(self.gene_info.split_exon_profiles.features))
        # logger.debug("Read exons profile" + str(read_split_exon_profile.read_profile))
        # logger.debug("Gene exons profile" + str(read_split_exon_profile.gene_profile))

        # logger.debug("< Read introns" + str(read_intron_profile.read_features))
        # logger.debug("Gene introns" + str(self.gene_info.intron_profiles.features))
        # logger.debug("Read intron profile" + str(read_intron_profile.read_profile))
        # logger.debug("Gene intron profile" + str(read_intron_profile.gene_profile))

        if all(el == 0 for el in read_split_exon_profile.read_profile) \
                or all(el == 0 for el in read_split_exon_profile.gene_profile):
            # none of the blocks matched
            if all(el == 0 for el in read_intron_profile.gene_profile):
                logger.debug("EMPTY - intergenic")
                assignment = ReadAssignment(read_id, ReadAssignmentType.empty, IsoformMatch(MatchClassification.intergenic))
            else:
                logger.debug("EMPTY - intronic")
                assignment = ReadAssignment(read_id, ReadAssignmentType.empty, IsoformMatch(MatchClassification.genic_intron))
            return assignment

        elif any(el == -1 for el in read_intron_profile.read_profile) \
                or any(el == -1 for el in read_split_exon_profile.read_profile):
            # Read has missing exons / introns
            logger.debug("+ Has contradictory features")
            assignment = self.match_contradictory(read_id, combined_read_profile)

        elif any(el == 0 for el in read_intron_profile.read_profile) \
                or any(el == 0 for el in read_split_exon_profile.read_profile):
            # Read has extra flanking exons / introns
            logger.debug("+ Has extra flanking features")
            assignment = self.match_with_extra_flanking(read_id, combined_read_profile)

        else:
            logger.debug("+ No contradictory features")
            if len(read_intron_profile.read_features) > 0:
                assignment = self.match_non_contradictory_spliced(read_id, combined_read_profile)
            else:
                assignment = self.match_non_contradictory_monoexonic(read_id, combined_read_profile)

            if assignment is None:
                # alternative isoforms made of known introns/exons or intron retention
                logger.debug("+ + Resolving unmatched ")
                assignment = self.match_contradictory(read_id, combined_read_profile)
            else:
                # check for extra flanking sequences
                self.check_for_extra_terminal_seqs(read_split_exon_profile, assignment)

        self.verify_polyA(combined_read_profile, assignment)
        return assignment

    def match_non_contradictory_spliced(self, read_id, combined_read_profile):
        """ match profile when all read features are assigned and read has at least one intron

        Parameters
        ----------
        read_id: str
        combined_read_profile: CombinedReadProfiles

        Returns
        -------

        """
        read_intron_profile = combined_read_profile.read_intron_profile
        assert read_intron_profile
        intron_matched_isoforms = self.find_matching_isoforms(read_intron_profile.gene_profile,
                                                              self.gene_info.intron_profiles.profiles)
        # logger.debug("Intron matched " + str(intron_matched_isoforms))

        read_assignment = None
        if len(intron_matched_isoforms) == 1 and len(read_intron_profile.read_features) > 0:
            isoform_id = list(intron_matched_isoforms)[0]

            logger.debug("+ + UNIQUE intron match found " + isoform_id)
            # logger.debug("Exon profile: " + str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
            # logger.debug("Intron profile: " + str(self.gene_info.intron_profiles.profiles[isoform_id]))

            isoform_match = self.categorize_correct_splice_match(read_intron_profile, isoform_id)
            read_assignment = ReadAssignment(read_id, ReadAssignmentType.unique, isoform_match)

        elif len(intron_matched_isoforms) > 1:
            logger.debug("+ + Multiexonic read, trying exon resolution")
            read_split_exon_profile = combined_read_profile.read_split_exon_profile
            exon_matched_isoforms = self.find_matching_isoforms(read_split_exon_profile.gene_profile,
                                                                self.gene_info.split_exon_profiles.profiles)
            # logger.debug("Exon matched " + str(exon_matched_isoforms))
            matched_isoforms = sorted(intron_matched_isoforms.intersection(exon_matched_isoforms))

            if len(matched_isoforms) == 1:
                isoform_id = list(matched_isoforms)[0]
                logger.debug("+ + Unique exon match found " + isoform_id)
                isoform_match = self.categorize_correct_splice_match(read_intron_profile, isoform_id)
                read_assignment = ReadAssignment(read_id, ReadAssignmentType.unique, isoform_match)
            elif len(matched_isoforms) == 0:
                read_assignment = self.resolve_multiple_assignments_spliced(read_id, combined_read_profile,
                                                                            sorted(intron_matched_isoforms))
            else:
                read_assignment = self.resolve_multiple_assignments_spliced(read_id, combined_read_profile,
                                                                            matched_isoforms)

        return read_assignment

    # resolve case when multiple isoforms are matched
    def resolve_multiple_assignments_spliced(self, read_id, combined_read_profile, matched_isoforms):
        read_intron_profile = combined_read_profile.read_intron_profile
        assert read_intron_profile

        # counter = 0
        # for isoform_id in matched_isoforms:
        #    logger.debug("+ + Detected match " + str(counter) + ": " + isoform_id)
        #    counter += 1
        #    logger.debug("Exon profile: " + str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
        #    logger.debug("Intron profile: " + str(self.gene_info.intron_profiles.profiles[isoform_id]))

        if self.params.resolve_ambiguous in (ExonAmbiguityResolvingMethod.all,
                                             ExonAmbiguityResolvingMethod.full_splice_matches_only):
            jaccard_matched_isoforms = \
                self.resolve_by_nucleotide_jaccard_similarity(combined_read_profile, matched_isoforms)
            jaccard_matched_isoforms = [x[0] for x in jaccard_matched_isoforms]
            if len(jaccard_matched_isoforms) > 1:
                matched_isoforms = jaccard_matched_isoforms
            elif len(jaccard_matched_isoforms) == 1:
                logger.debug("Successfully resolved using Jaccard similarity")
                isoform_id = list(jaccard_matched_isoforms)[0]
                isoform_match = self.categorize_correct_splice_match(read_intron_profile, isoform_id)
                return ReadAssignment(read_id, ReadAssignmentType.unique, isoform_match)

        read_intron_profile = combined_read_profile.read_intron_profile
        isoform_matches = self.categorize_multiple_splice_matches(read_intron_profile, matched_isoforms)
        return ReadAssignment(read_id, ReadAssignmentType.ambiguous, isoform_matches)

    def match_non_contradictory_monoexonic(self, read_id, combined_read_profile):
        logger.debug("+  Resolving monoexonic read")
        read_split_exon_profile = combined_read_profile.read_split_exon_profile
        isoform_split_exon_profiles = self.gene_info.split_exon_profiles.profiles
        overlapping_isoforms = set()
        for isoform_id, isoform_profile in isoform_split_exon_profiles.items():
            if has_overlapping_features(read_split_exon_profile.gene_profile, isoform_profile):
                overlapping_isoforms.add(isoform_id)
        # non-contradictory and non empty profile - some exons must match
        assert overlapping_isoforms

        jaccard_matched_isoforms = self.resolve_by_nucleotide_jaccard_similarity(combined_read_profile,
                                                                                 sorted(overlapping_isoforms),
                                                                                 min_similarity=0.0)
        reliable_isoforms = [
            x[0] for x in filter(lambda x: x[1] >= ExonAmbiguityResolvingMethod.minimal_jaccard_similarity.value,
                                 jaccard_matched_isoforms)
        ]

        if len(reliable_isoforms) == 1:
            isoform_id = list(reliable_isoforms)[0]
            logger.debug("Jaccard similarity picked a single isoform: %s " % isoform_id)
            isoform_match = self.categorize_mono_exonic_read_match(combined_read_profile, isoform_id)
            if isoform_match.all_subtypes_are_none_or_monoexonic():
                assignment_type = ReadAssignmentType.unique
            else:
                assignment_type = ReadAssignmentType.contradictory
            read_assignment = ReadAssignment(read_id, assignment_type, isoform_match)
        elif len(reliable_isoforms) > 0:
            logger.debug("Jaccard similarity picked multiple isoforms")
            isoform_matches = self.categorize_multiple_mono_exon_matches(combined_read_profile, reliable_isoforms)
            if all(im.all_subtypes_are_none_or_monoexonic() for im in isoform_matches):
                assignment_type = ReadAssignmentType.ambiguous
            else:
                assignment_type = ReadAssignmentType.contradictory
            read_assignment = ReadAssignment(read_id, assignment_type, isoform_matches)
        else:
            # unreliable_mono_exon
            logger.debug("Jaccard similarity is not reliable")
            unreliable_isoforms = [x[0] for x in jaccard_matched_isoforms]
            isoform_matches = self.categorize_multiple_mono_exon_matches(combined_read_profile, unreliable_isoforms)
            if all(im.all_subtypes_are_none_or_monoexonic() for im in isoform_matches):
                assignment_type = ReadAssignmentType.unreliable
            else:
                assignment_type = ReadAssignmentType.contradictory
            read_assignment = ReadAssignment(read_id, assignment_type, isoform_matches)

        return read_assignment

    # select best assignment based on nucleotide similarity
    def resolve_by_nucleotide_jaccard_similarity(self, combined_read_profile, matched_isoforms,
                                                 min_similarity=ExonAmbiguityResolvingMethod.minimal_jaccard_similarity.value,
                                                 top_scored_factor=ExonAmbiguityResolvingMethod.top_scored_factor.value):
        assert matched_isoforms
        logger.debug("+ + + Resolving by Jaccard similarity")
        jaccard_similarities = []
        read_exons = combined_read_profile.read_exon_profile.read_features
        for isoform_id in matched_isoforms:
            isoform_exons = self.gene_info.all_isoforms_exons[isoform_id]
            jaccard_similarities.append((isoform_id, jaccard_similarity(read_exons, isoform_exons)))

        # logger.debug(jaccard_similarities)
        best_score = max([x[1] for x in jaccard_similarities])
        top_scored = sorted(filter(lambda x: x[1] * top_scored_factor >= best_score and x[1] >= min_similarity,
                                   jaccard_similarities))
        logger.debug("+ + + Best score = %f, total candidates = %d" % (best_score, len(top_scored)))
        return top_scored

    # resolve when there are 0s  at the ends of read profile
    def match_with_extra_flanking(self, read_id, combined_read_profile):
        read_intron_profile = combined_read_profile.read_intron_profile
        if read_intron_profile.read_profile[0] == 1 and read_intron_profile.read_profile[-1] == 1:
            logger.debug("+ + Both terminal introns present, odd case")

        assignment = self.match_non_contradictory_spliced(read_id, combined_read_profile)
        if assignment is None:
            # alternative isoforms made of known introns/exons or intron retention
            logger.debug("+ + Resolving unmatched with extra flanking ")
            return self.match_contradictory(read_id, combined_read_profile)

        for match in assignment.isoform_matches:
            match.set_classification(MatchClassification.novel_not_in_catalog)
            logger.debug("+ + Found extra terminal introns " + str(read_intron_profile.read_profile))
            isoform_start = self.gene_info.transcript_start(match.assigned_transcript)
            self.intron_comparator.add_extra_out_exon_events(match.match_subclassifications, read_intron_profile.read_profile,
                                                             read_intron_profile.read_features, isoform_start)

        if not self.params.allow_extra_terminal_introns:
            if assignment.assignment_type in (ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference):
                assignment.set_assignment_type(ReadAssignmentType.contradictory)
        else:
            if assignment.assignment_type == ReadAssignmentType.unique:
                assignment.set_assignment_type(ReadAssignmentType.unique_minor_difference)

        return assignment

    # resolve when there are -1s in read profile or when there are no exactly matching isoforms, but no -1s in read profiles
    def match_contradictory(self, read_id, combined_read_profile):
        read_intron_profile = combined_read_profile.read_intron_profile
        read_split_exon_profile = combined_read_profile.read_split_exon_profile

        overlapping_isoforms = self.find_overlapping_isoforms(read_split_exon_profile.gene_profile,
                                                              self.gene_info.split_exon_profiles.profiles)
        intron_matching_isoforms = self.match_profile(read_intron_profile.gene_profile,
                                                      self.gene_info.intron_profiles.profiles,
                                                      hint=overlapping_isoforms)
        exon_matching_isoforms = self.match_profile(read_split_exon_profile.gene_profile,
                                                    self.gene_info.split_exon_profiles.profiles,
                                                    hint=overlapping_isoforms)
        return self.detect_differences(read_id, read_intron_profile, read_split_exon_profile,
                                       intron_matching_isoforms, exon_matching_isoforms)

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
        logger.debug("Detecting difference for %s, %d intron matched and %d exon matched" %
                     (read_id, len(intron_matching_isoforms), len(exon_matching_isoforms)))
        best_intron_isoform_ids = get_first_best_from_sorted(intron_matching_isoforms)
        # TODO add Jaccard similarity instead of simple exon matching?
        best_exon_isoforms = list(filter(lambda x: x[0] in best_intron_isoform_ids, exon_matching_isoforms))
        best_isoform_ids = get_first_best_from_sorted(best_exon_isoforms)
        if not best_isoform_ids:
            return ReadAssignment(read_id, ReadAssignmentType.empty)

        assignment = ReadAssignment(read_id, ReadAssignmentType.contradictory)
        logger.debug("+ + Closest matching isoforms " + str(best_isoform_ids))
        isoform_id = best_isoform_ids[0]
        # logger.debug(str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
        # logger.debug(str(self.gene_info.intron_profiles.profiles[isoform_id]))

        for isoform_id in sorted(best_isoform_ids):
            logger.debug("Checking isoform %s" % isoform_id)
            # get intron coordinates
            isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
            # get_blocks_from_profile(self.gene_info.intron_profiles.features, self.gene_info.intron_profiles.profiles[isoform_id])
            # read start-end coordinates
            read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[-1][1])
            # isoform start-end
            isoform_start = self.gene_info.transcript_start(isoform_id)
            isoform_end = self.gene_info.transcript_end(isoform_id)
            isoform_region = (isoform_start, isoform_end)

            matching_events = self.intron_comparator.compare_junctions(read_intron_profile.read_features, read_region,
                                                                       isoform_introns, isoform_region)
            if len(matching_events) == 1 and matching_events[0].event_type == MatchEventSubtype.undefined:
                continue
            match_classification = MatchClassification.get_contradiction_classification_from_subtypes(matching_events)
            isoform_match = IsoformMatch(match_classification, self.get_gene_id(isoform_id), isoform_id)
            for m in matching_events:
                isoform_match.add_subclassification(m)
            assignment.add_match(isoform_match)
            logger.debug("+ + Found contradiction for " + isoform_id + ": " + " ".join(map(lambda x: x.event_type.name, matching_events)))

        new_assignment_type = None
        if len(assignment.isoform_matches) == 0:
            return ReadAssignment(read_id, ReadAssignmentType.empty)
        # Change assignment from contradictory when contradictions are minor or absent
        if all(m.all_subtypes_are_none_or_monoexonic() for m in assignment.isoform_matches):
            # No contradiction
            new_assignment_type = ReadAssignmentType.unique if len(best_isoform_ids) == 1 else ReadAssignmentType.ambiguous
        elif self.params.correct_minor_errors and all(m.all_subtypes_are_alignment_artifacts() for m in assignment.isoform_matches):
            # Only alignment artifacts
            new_assignment_type = ReadAssignmentType.unique_minor_difference if len(best_isoform_ids) == 1 else ReadAssignmentType.ambiguous

        if new_assignment_type is not None:
            # Revise all matches as correct
            isoform_matches = self.categorize_multiple_splice_matches(read_intron_profile, best_isoform_ids)
            assignment = ReadAssignment(read_id, new_assignment_type, isoform_matches)
            self.check_for_extra_terminal_seqs(read_split_exon_profile, assignment)
        return assignment

    # check consistency with polyA
    def verify_polyA(self, combined_read_profile, read_assignment):
        if not self.params.has_polya:
            return
        logger.debug("+ Validating polyA/T sites")
        if combined_read_profile.polya_pos == -1 and combined_read_profile.polyt_pos == -1:
            logger.debug("+ No sites found, ciao")
            return

        if read_assignment.assignment_type == ReadAssignmentType.ambiguous:
            self.resolve_by_polyA(combined_read_profile, read_assignment)

        apa_found = False
        pos = -1
        for isoform_match in read_assignment.isoform_matches:
            isoform_id = isoform_match.assigned_transcript
            logger.debug("+ Checking isoform %s" % isoform_id)
            if self.gene_info.isoform_strands[isoform_id] == '+' and combined_read_profile.polya_pos != -1:
                isoform_end = self.gene_info.transcript_end(isoform_id)
                pos = combined_read_profile.polya_pos
                dist_to_polya = abs(isoform_end - pos)
                logger.debug("+ Distance to polyA is %d" % dist_to_polya)
                if dist_to_polya > self.params.apa_delta:
                    apa_found = True
                    logger.debug("+ Seems like APA site")
            elif self.gene_info.isoform_strands[isoform_id] == '-' and combined_read_profile.polyt_pos != -1:
                isoform_start = self.gene_info.transcript_start(isoform_id)
                pos = combined_read_profile.polyt_pos
                dist_to_polyt = abs(pos - isoform_start)
                logger.debug("+ Distance to polyT is %d" % dist_to_polyt)
                if dist_to_polyt > self.params.apa_delta:
                    apa_found = True
                    logger.debug("+ Seems like APA site")

            if apa_found:
                isoform_match.add_subclassification(make_event(MatchEventSubtype.alternative_polya_site))
                isoform_match.set_classification(MatchClassification.novel_in_catalog)

        if apa_found and read_assignment.assignment_type == ReadAssignmentType.unique:
            read_assignment.assignment_type = ReadAssignmentType.contradictory

    # try to resolve when polyA position is known
    def resolve_by_polyA(self, combined_read_profile, read_assignment):
        logger.debug("+ + Resolving ambiguous case with polyA: " + str())
        strands = set()
        isoform_ids = []
        for isoform_match in read_assignment.isoform_matches:
            isoform_id = isoform_match.assigned_transcript
            isoform_ids.append(isoform_id)
            strands.add(self.gene_info.isoform_strands[isoform_id])
        logger.debug("+ + Resolving ambiguous case with polyA: " + str(isoform_ids))
        # we need to know strand to check polyA or polyT
        if len(strands) > 1:
            logger.debug("+ + Ambiguous strands")
            return
        strand = list(strands)[0]

        if strand == "+" and combined_read_profile.polya_pos != -1:
            # find distances to polyA site
            distances = self.find_distances_to_pos(isoform_ids, self.gene_info.transcript_end,
                                                   combined_read_profile.polya_pos)
            best_distance = min([x[0] for x in distances])
            # select best isoforms
            closest_isoforms = [x[1] for x in filter(lambda y:y[0] <= 2 * best_distance, distances)]
            logger.debug("+ + PolyA helped: " + str(closest_isoforms))
        else:
            distances = self.find_distances_to_pos(isoform_ids, self.gene_info.transcript_start,
                                                   combined_read_profile.polyt_pos)
            best_distance = min([x[0] for x in distances])
            closest_isoforms = [x[1] for x in filter(lambda y:y[0] <= 2 * best_distance, distances)]
            logger.debug("+ + PolyT helped: " + str(closest_isoforms))

        closest_isoforms = set(closest_isoforms)
        if len(closest_isoforms) < len(isoform_ids):
            # set new matches if polyA did filter something out
            read_assignment.isoform_matches = list(filter(lambda im: im.assigned_transcript in closest_isoforms,
                                                          read_assignment.isoform_matches))
            logger.debug("+ + Assigning new matches")
            if len(closest_isoforms) == 1:
                logger.debug("+ + and it's a unique match")
                read_assignment.set_assignment_type(ReadAssignmentType.unique)

    def find_distances_to_pos(self, isoform_ids, transcript_function, pos):
        distances = []
        for isoform_id in isoform_ids:
            distances.append((abs(transcript_function(isoform_id) - pos), isoform_id))
        return distances
