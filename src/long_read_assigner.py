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

    minimal_score = 0.0
    top_scored_factor = 1.5


class LongReadAssigner:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params
        self.intron_comparator = JunctionComparator(params,
                                                    OverlappingFeaturesProfileConstructor
                                                    (self.gene_info.intron_profiles.features,
                                                     (self.gene_info.start, self.gene_info.end),
                                                     comparator=partial(equal_ranges, delta=self.params.delta)))

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
    def find_overlapping_isoforms(read_exon_split_profile, isoform_profiles):
        isoforms = set()
        for isoform_id, isoform_profile in isoform_profiles.items():
            if has_overlapping_features(isoform_profile, read_exon_split_profile.gene_profile):
                isoforms.add(isoform_id)
        return isoforms

    def find_containing_isoforms(self, read_exon_split_profile, isoform_profiles, hint=None):
        isoforms = set()
        read_exons = read_exon_split_profile.read_features
        for isoform_id, isoform_profile in isoform_profiles.items():
            if hint and isoform_id not in hint:
                continue
            isoform_region = (self.gene_info.transcript_start(isoform_id) - self.params.delta,
                              self.gene_info.transcript_end(isoform_id) + self.params.delta)
            read_region = (read_exons[0][0], read_exons[-1][1])
            if contains(isoform_region, read_region):
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
        return set([x[0] for x in filter(lambda x: x[1] == 0, isoforms)])


    # select most similar isoform based on different criteria
    def select_similar_isoforms(self, read_id, combined_read_profile):
        read_split_exon_profile = combined_read_profile.read_split_exon_profile
        overlapping_isoforms = self.find_overlapping_isoforms(read_split_exon_profile,
                                                              self.gene_info.split_exon_profiles.profiles)
        if not overlapping_isoforms:
            return
        # select isoforms with non-negative nucleotide score
        significantly_overlapping_isoforms = self.resolve_by_nucleotide_score(combined_read_profile,
                                                                               overlapping_isoforms,
                                                                               top_scored_factor=2)
        if not significantly_overlapping_isoforms:
            return

        # find intron matches
        intron_matching_isoforms = self.match_profile(combined_read_profile.read_intron_profile.gene_profile,
                                                      self.gene_info.intron_profiles.profiles,
                                                      hint=significantly_overlapping_isoforms)

        # add extra terminal bases as potential inconsistency event
        candidates = []
        read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[-1][1])
        for isoform_id, diff_introns in intron_matching_isoforms:
            extra_left = 1 if read_region[0] + self.params.delta < self.gene_info.transcript_start(isoform_id) else 0
            extra_right = 1 if read_region[0] - self.params.delta > self.gene_info.transcript_end(isoform_id) else 0
            candidates.append((isoform_id, diff_introns + extra_right + extra_left))
        # select isoforms that have similar number of potential inconsistencies
        best_diff = min(candidates, key=lambda x: x[1])[1]
        best_candidates = [x[0] for x in filter(lambda x: x[1] <= best_diff + 1, candidates)]
        logger.debug("+ + Closest matching isoforms " + str(best_candidates))
        return best_candidates

    # detect exon elongation subtyp
    def categorize_exon_elongation_subtype(self, read_split_exon_profile, isoform_id):
        split_exons = self.gene_info.split_exon_profiles.features
        isoform_profile = self.gene_info.split_exon_profiles.profiles[isoform_id]

        # find first and last common exons
        common_first_exon = -1
        isofrom_first_exon = -1
        for i in range(len(split_exons)):
            if isoform_profile[i] != 1:
                continue
            if isofrom_first_exon == -1:
                isofrom_first_exon = i
            if read_split_exon_profile.gene_profile[i] == 1:
                common_first_exon = i
                break

        common_last_exon = -1
        isofrom_last_exon = -1
        for i in range(len(split_exons)):
            index = len(split_exons) - i - 1
            if isoform_profile[index] != 1:
                continue
            if isofrom_last_exon == -1:
                isofrom_last_exon = index
            if read_split_exon_profile.gene_profile[index] == 1:
                common_last_exon = index
                break

        if common_first_exon == -1 or common_last_exon == -1:
            logger.warning(" + Werid case for exon elongation, no matching exons")

        first_read_exon = read_split_exon_profile.read_features[0]
        if overlaps(first_read_exon, split_exons[common_first_exon]):
            extra_left = split_exons[common_first_exon][0] - first_read_exon[0]
        else:
            # first read exon seems to be to the left of common exon
            extra_left = 0

        last_read_exon = read_split_exon_profile.read_features[-1]
        if overlaps(last_read_exon, split_exons[common_last_exon]):
            extra_right = last_read_exon[1] - split_exons[common_last_exon][1]
        else:
            extra_right = 0

        logger.debug("+ + Checking exon elongation")
        logger.debug("+ + Extra bases: left = %d, right = %d" % (extra_left, extra_right))

        events = []
        # TODO add intron retention position and process it in transcript model constructor
        if extra_left > self.params.max_exon_extension:
            if common_first_exon == isofrom_first_exon:
                events.append(MatchEventSubtype.major_exon_elongation_left)
            else:
                events.append(MatchEventSubtype.incomplete_intron_retention)
        elif extra_left > self.params.delta:
            events.append(MatchEventSubtype.exon_elongation_left)

        if extra_right > self.params.max_exon_extension:
            if common_last_exon == isofrom_last_exon:
                events.append(MatchEventSubtype.major_exon_elongation_right)
            else:
                events.append(MatchEventSubtype.incomplete_intron_retention)
        elif extra_right > self.params.delta:
            events.append(MatchEventSubtype.exon_elongation_right)

        return list(map(make_event, events))

    # select best assignment based on nucleotide similarity
    # score = Jaccard similarity - flanking len / read len ([-1,1])
    # gives priority to isoforms that contain the read
    def resolve_by_nucleotide_score(self, combined_read_profile, matched_isoforms,
                                    min_similarity=ExonAmbiguityResolvingMethod.minimal_score.value,
                                    top_scored_factor=ExonAmbiguityResolvingMethod.top_scored_factor.value):
        if not matched_isoforms:
            return []

        logger.debug("+ + + Resolving by nucleotide similarity")
        read_exons = combined_read_profile.read_exon_profile.read_features
        read_split_exon_profile = combined_read_profile.read_split_exon_profile

        scores = []
        for isoform_id in matched_isoforms:
            isoform_exons = self.gene_info.all_isoforms_exons[isoform_id]
            js = jaccard_similarity(read_exons, isoform_exons)

            isoform_extended_region = (self.gene_info.transcript_start(isoform_id) - self.params.max_exon_extension,
                                       self.gene_info.transcript_end(isoform_id) + self.params.max_exon_extension)
            flanking_percentage = extra_exon_percentage(isoform_extended_region, read_exons)

            scores.append((isoform_id, js - flanking_percentage))

        # logger.debug(jaccard_similarities)
        best_score = max([x[1] for x in scores])
        if top_scored_factor == 0:
            top_scored = sorted(filter(lambda x: x[1] >= min_similarity, scores))
        else:
            top_scored = sorted(filter(lambda x: x[1] * top_scored_factor >= best_score and x[1] >= min_similarity,
                                       scores))
        logger.debug("+ + + Best score = %f, all candidates = %s" % (best_score, str(top_scored)))
        return list(map(lambda x: x[0], top_scored))

    # ====== CLASSIFICATION =======

    # check for extra sequences and modify assignment accordingly
    def check_for_extra_terminal_seqs(self, read_split_exon_profile, assignment):
        for match in assignment.isoform_matches:
            if match.assigned_transcript is None:
                continue

            #if any(s.event_type in {MatchEventSubtype.extra_intron_flanking_left, MatchEventSubtype.extra_intron_flanking_right}
            #       for s in match.match_subclassifications):
            #    continue

            exon_elongation_types = self.categorize_exon_elongation_subtype(read_split_exon_profile,
                                                                            match.assigned_transcript)
            for e in exon_elongation_types:
                match.add_subclassification(make_event(e))
            if  any(MatchEventSubtype.is_major_elongation(e) for e in exon_elongation_types):
                # serious exon elongation
                if assignment.assignment_type in {ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference}:
                    # match.set_classification
                    assignment.set_assignment_type(ReadAssignmentType.inconsistent)
            elif len(exon_elongation_types) > 0:
                # minor exon elongation
                if assignment.assignment_type == ReadAssignmentType.unique:
                    assignment.set_assignment_type(ReadAssignmentType.unique_minor_difference)

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
            event_type = MatchEventSubtype.ism_left
        elif is_right_truncated:
            logger.debug(" + Truncated on the right")
            event_type = MatchEventSubtype.ism_right
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

    # make proper match subtype
    def categorize_correct_splice_match(self, read_intron_profile, isoform_id):
        if self.is_fsm(read_intron_profile, isoform_id):
            logger.debug("+ + Full splice match " + isoform_id)
            isoform_match = IsoformMatch(MatchClassification.full_splice_match,
                                         self.get_gene_id(isoform_id),
                                         isoform_id, make_event(MatchEventSubtype.fsm),
                                         self.gene_info.isoform_strands[isoform_id])
        else:
            logger.debug("+ + Incomplete splice match " + isoform_id)
            isoform_match = IsoformMatch(MatchClassification.incomplete_splice_match,
                                         self.get_gene_id(isoform_id),
                                         isoform_id,
                                         self.detect_ism_subtype(read_intron_profile, isoform_id),
                                         self.gene_info.isoform_strands[isoform_id])
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

        matching_events = self.intron_comparator.get_mono_exon_subtype(read_region, isoform_introns)
        match_classification = MatchClassification.get_mono_exon_classification(matching_events)

        return IsoformMatch(match_classification, self.get_gene_id(isoform_id), isoform_id, matching_events,
                            self.gene_info.isoform_strands[isoform_id])

    def categorize_multiple_mono_exon_matches(self, combined_read_profile, isoform_ids):
        isoform_matches = []
        for isoform_id in isoform_ids:
            isoform_matches.append(self.categorize_mono_exonic_read_match(combined_read_profile, isoform_id))
        return isoform_matches

    # =========== MAIN PART ============

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

            read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[-1][1])
            gene_region = (self.gene_info.split_exon_profiles.features[0][0],
                              self.gene_info.split_exon_profiles.features[-1][1])
            # none of the blocks matched
            if not overlaps(read_region, gene_region):
                assignment = ReadAssignment(read_id, ReadAssignmentType.noninformative, IsoformMatch(MatchClassification.intergenic))
            elif all(el == 0 for el in read_split_exon_profile.gene_profile):
                logger.debug("EMPTY - intronic")
                assignment = ReadAssignment(read_id, ReadAssignmentType.noninformative, IsoformMatch(MatchClassification.genic_intron))
            else:
                logger.debug("EMPTY - genic")
                assignment = ReadAssignment(read_id, ReadAssignmentType.noninformative, IsoformMatch(MatchClassification.genic))
            return assignment

        elif any(el == -1 for el in read_intron_profile.read_profile) \
                or any(el == -1 for el in read_split_exon_profile.read_profile):
            # Read has inconsistent exons / introns
            logger.debug("+ Has inconsistent features (novel intron / exons)")
            assignment = self.match_inconsistent(read_id, combined_read_profile)

        elif any(el == 0 for el in read_intron_profile.read_profile) \
                or any(el == 0 for el in read_split_exon_profile.read_profile):
            # Read has extra flanking exons / introns
            logger.debug("+ Has extra flanking features")
            assignment = self.match_inconsistent(read_id, combined_read_profile)

        else:
            logger.debug("+ No contradictory features in read, but no consistent isoforms still can be found")
            if len(read_intron_profile.read_features) > 0:
                assignment = self.match_consistent_spliced(read_id, combined_read_profile)
            else:
                assignment = self.match_consistent_unspliced(read_id, combined_read_profile)

            if assignment is None:
                # alternative isoforms made of known introns/exons or intron retention
                logger.debug("+ + Resolving unmatched ")
                assignment = self.match_inconsistent(read_id, combined_read_profile)

        # checking polyA
        self.verify_polyA(combined_read_profile, assignment)
        return assignment

    def match_consistent_spliced(self, read_id, combined_read_profile):
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

        read_split_exon_profile = combined_read_profile.read_split_exon_profile
        isoform_split_exon_profiles = self.gene_info.split_exon_profiles.profiles
        overlapping_isoforms = self.find_overlapping_isoforms(read_split_exon_profile, isoform_split_exon_profiles)
        assert overlapping_isoforms
        containing_isoforms = self.find_containing_isoforms(read_split_exon_profile, isoform_split_exon_profiles,
                                                            hint=overlapping_isoforms)
        if not containing_isoforms:
            return None
        intron_matched_isoforms = self.find_matching_isoforms(read_intron_profile.gene_profile,
                                                              self.gene_info.intron_profiles.profiles,
                                                              hint=containing_isoforms)
        # logger.debug("Intron matched " + str(intron_matched_isoforms))

        read_assignment = None
        if len(intron_matched_isoforms) == 1 and len(read_intron_profile.read_features) > 0:
            isoform_id = list(intron_matched_isoforms)[0]

            logger.debug("+ + UNIQUE intron match found " + isoform_id)
            # logger.debug("Exon profile: " + str(isoform_split_exon_profiles[isoform_id]))
            # logger.debug("Intron profile: " + str(self.gene_info.intron_profiles.profiles[isoform_id]))

            isoform_match = self.categorize_correct_splice_match(read_intron_profile, isoform_id)
            read_assignment = ReadAssignment(read_id, ReadAssignmentType.unique, isoform_match)

        elif len(intron_matched_isoforms) > 1:
            logger.debug("+ + Multiexonic read, trying exon resolution")
            exon_matched_isoforms = self.find_matching_isoforms(read_split_exon_profile.gene_profile,
                                                                isoform_split_exon_profiles,
                                                                hint=intron_matched_isoforms)
            # logger.debug("Exon matched " + str(exon_matched_isoforms))
            if len(exon_matched_isoforms) == 1:
                isoform_id = list(exon_matched_isoforms)[0]
                logger.debug("+ + Unique exon match found " + isoform_id)
                isoform_match = self.categorize_correct_splice_match(read_intron_profile, isoform_id)
                read_assignment = ReadAssignment(read_id, ReadAssignmentType.unique, isoform_match)
            elif len(exon_matched_isoforms) == 0:
                read_assignment = self.resolve_multiple_assignments_spliced(read_id, combined_read_profile,
                                                                            sorted(intron_matched_isoforms))
            else:
                read_assignment = self.resolve_multiple_assignments_spliced(read_id, combined_read_profile,
                                                                            exon_matched_isoforms)

        return read_assignment

    # resolve case when multiple isoforms are matched
    def resolve_multiple_assignments_spliced(self, read_id, combined_read_profile, matched_isoforms):
        read_intron_profile = combined_read_profile.read_intron_profile
        assert read_intron_profile

        # FIXME full_splice_matches_only
        if self.params.resolve_ambiguous in (ExonAmbiguityResolvingMethod.all,
                                             ExonAmbiguityResolvingMethod.full_splice_matches_only):
            nucl_matched_isoforms = \
                self.resolve_by_nucleotide_score(combined_read_profile, matched_isoforms)
            nucl_matched_isoforms = [x[0] for x in nucl_matched_isoforms]
            if len(nucl_matched_isoforms) > 1:
                matched_isoforms = nucl_matched_isoforms
            elif len(nucl_matched_isoforms) == 1:
                logger.debug("Successfully resolved using Nucleotide similarity")
                isoform_id = list(nucl_matched_isoforms)[0]
                isoform_match = self.categorize_correct_splice_match(read_intron_profile, isoform_id)
                return ReadAssignment(read_id, ReadAssignmentType.unique, isoform_match)

        read_intron_profile = combined_read_profile.read_intron_profile
        isoform_matches = self.categorize_multiple_splice_matches(read_intron_profile, matched_isoforms)
        return ReadAssignment(read_id, ReadAssignmentType.ambiguous, isoform_matches)

    def match_consistent_unspliced(self, read_id, combined_read_profile):
        logger.debug("+  Resolving monoexonic read")
        read_split_exon_profile = combined_read_profile.read_split_exon_profile
        isoform_split_exon_profiles = self.gene_info.split_exon_profiles.profiles

        overlapping_isoforms = self.find_overlapping_isoforms(read_split_exon_profile, isoform_split_exon_profiles)
        if len(overlapping_isoforms) == 0:
            return ReadAssignment(read_id, ReadAssignmentType.noninformative,
                                  IsoformMatch(MatchClassification.genic_intron))
        containing_isoforms = self.find_containing_isoforms(read_split_exon_profile, isoform_split_exon_profiles,
                                                            hint=overlapping_isoforms)
        if not containing_isoforms:
            return None

        matched_isoforms = self.resolve_by_nucleotide_score(combined_read_profile, sorted(containing_isoforms))
        if len(matched_isoforms) == 1:
            isoform_id = matched_isoforms[0]
            logger.debug("Nucleotide similarity picked a single isoform: %s" % isoform_id)
            isoform_match = self.categorize_mono_exonic_read_match(combined_read_profile, isoform_id)
            if isoform_match.monoexon_is_consistent():
                assignment_type = ReadAssignmentType.unique
            else:
                assignment_type = ReadAssignmentType.inconsistent
            read_assignment = ReadAssignment(read_id, assignment_type, isoform_match)

        elif len(matched_isoforms) > 0:
            logger.debug("Nucleotide similarity picked multiple isoforms")
            isoform_matches = self.categorize_multiple_mono_exon_matches(combined_read_profile, matched_isoforms)
            if all(im.monoexon_is_consistent() for im in isoform_matches):
                assignment_type = ReadAssignmentType.ambiguous
            else:
                assignment_type = ReadAssignmentType.inconsistent
            read_assignment = ReadAssignment(read_id, assignment_type, isoform_matches)

        else:
            logger.debug("Nucleotide similarity is not reliable")
            read_assignment = ReadAssignment(read_id, ReadAssignmentType.noninformative,
                                             IsoformMatch(MatchClassification.genic))

        return read_assignment

    # resolve when there are -1s in read profile or when there are no exactly matching isoforms, but no -1s in read profiles
    def match_inconsistent(self, read_id, combined_read_profile):
        # select most similar isoforms based on multiple criteria
        best_candidates = self.select_similar_isoforms(read_id, combined_read_profile)
        if not best_candidates:
            return ReadAssignment(read_id, ReadAssignmentType.noninformative, IsoformMatch(MatchClassification.genic))

        logger.debug("* Best candidates for inconsistency detection: " + str(best_candidates))
        # detect inconsistency for each one
        read_matches = self.detect_inconsistensies(read_id, combined_read_profile, best_candidates)
        if not read_matches:
            return ReadAssignment(read_id, ReadAssignmentType.noninformative)
        logger.debug("* Inconsistencies detected: " + str(read_matches))

        # select ones with least number of inconsistent events
        best_isoforms = self.select_best_among_inconsistent(read_matches)
        if not best_isoforms:
            return ReadAssignment(read_id, ReadAssignmentType.noninformative)

        is_abmiguous = len(best_isoforms) > 1
        all_event_types = set()
        for isoform_id in best_isoforms:
            for e in read_matches[isoform_id]:
                all_event_types.add(e.event_type)

        logger.debug("* Selected isoforms: " + str(best_isoforms) + ", all events: " + str(all_event_types))
        if all(MatchEventSubtype.is_consistent(e) for e in all_event_types):
            logger.debug("* * Assignment seems to be consistent")
            assignment_type = ReadAssignmentType.ambiguous if is_abmiguous else ReadAssignmentType.unique
        elif any(MatchEventSubtype.is_major_inconsistency(e) for e in all_event_types):
            logger.debug("* * Assignment is INconsistent")
            assignment_type = ReadAssignmentType.inconsistent
        elif any(MatchEventSubtype.is_minor_error(e) for e in all_event_types):
            logger.debug("* * Assignment has minor errors")
            assignment_type = ReadAssignmentType.ambiguous if is_abmiguous else ReadAssignmentType.unique_minor_difference
        else:
            logger.warning("Unexpected event reported: " + str(all_event_types))
            assignment_type = ReadAssignmentType.noninformative

        if not combined_read_profile.read_intron_profile.read_profile:
            isoform_matches = self.create_monoexon_matches(read_matches, best_isoforms)
        elif assignment_type == ReadAssignmentType.inconsistent:
            isoform_matches = self.create_inconsistent_matches(read_matches, best_isoforms)
        else:
            isoform_matches = self.create_consistent_matches(read_matches, best_isoforms, combined_read_profile)
        return ReadAssignment(read_id, assignment_type, isoform_matches)

    def create_monoexon_matches(self, read_matches, selected_isoforms):
        matches = []
        for isoform_id in selected_isoforms:
            match_classification = MatchClassification.get_mono_exon_classification(read_matches[isoform_id])
            isoform_match = IsoformMatch(match_classification, self.get_gene_id(isoform_id), isoform_id,
                                         read_matches[isoform_id], self.gene_info.isoform_strands[isoform_id])
            matches.append(isoform_match)
        return matches

    def create_consistent_matches(self, read_matches, selected_isoforms, combined_read_profile):
        matches = []
        for isoform_id in selected_isoforms:
            isoform_match = self.categorize_correct_splice_match(combined_read_profile.read_intron_profile, isoform_id)
            for e in read_matches[isoform_id]:
                if e.event_type != MatchEventSubtype.none:
                    isoform_match.add_subclassification(e)
            matches.append(isoform_match)
        return matches

    def create_inconsistent_matches(self, read_matches, selected_isoforms):
        matches = []
        for isoform_id in selected_isoforms:
            match_classification = MatchClassification.get_inconsistency_classification(read_matches[isoform_id])
            isoform_match = IsoformMatch(match_classification, self.get_gene_id(isoform_id), isoform_id,
                                         read_matches[isoform_id], self.gene_info.isoform_strands[isoform_id])
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
        logger.debug("Detecting difference for %s, %d matched isoforms" % (read_id, len(matched_isoforms)))
        read_intron_profile = combined_read_profile.read_intron_profile
        read_split_exon_profile = combined_read_profile.read_split_exon_profile

        logger.debug("+ + Closest matching isoforms " + str(matched_isoforms))
        # isoform_id = best_isoform_ids[0]
        # logger.debug(str(self.gene_info.split_exon_profiles.profiles[isoform_id]))
        # logger.debug(str(self.gene_info.intron_profiles.profiles[isoform_id]))

        read_matches = {}
        for isoform_id in sorted(matched_isoforms):
            logger.debug("Checking isoform %s" % isoform_id)
            # get intron coordinates
            isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
            # read start-end coordinates
            read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[-1][1])
            # isoform start-end
            isoform_region = (self.gene_info.transcript_start(isoform_id), self.gene_info.transcript_end(isoform_id))

            matching_events = self.intron_comparator.compare_junctions(read_intron_profile.read_features, read_region,
                                                                       isoform_introns, isoform_region)

            if len(matching_events) == 1 and matching_events[0].event_type == MatchEventSubtype.undefined:
                continue

            matching_events += self.categorize_exon_elongation_subtype(read_split_exon_profile, isoform_id)
            #match_classification = MatchClassification.get_contradiction_classification_from_subtypes(matching_events)
            read_matches[isoform_id] = matching_events

            logger.debug("+ + Found contradiction for " + isoform_id + ": " + " ".join(map(lambda x: x.event_type.name, matching_events)))
        return read_matches

    def select_best_among_inconsistent(self, read_matches):
        isoform_scores = []
        for isoform_id, match_events in read_matches.items():
            logger.debug("* * Scoring inconsistencies for " + isoform_id + ": " + str(match_events))
            score = 0.0
            for e in match_events:
                event_count = 1
                if e.read_region != SupplementaryMatchConstansts.undefined_region:
                    event_count = e.read_region[1] - e.read_region[0] + 1
                    if e.event_type in {MatchEventSubtype.exon_gain_novel,
                                        MatchEventSubtype.exon_gain_known,
                                        MatchEventSubtype.mutually_exclusive_exons_novel,
                                        MatchEventSubtype.mutually_exclusive_exons_known}:
                        event_count -= 1
                    elif e.event_type in {MatchEventSubtype.intron_retention,
                                          MatchEventSubtype.unspliced_intron_retention,
                                          MatchEventSubtype.unspliced_genic,
                                          MatchEventSubtype.incomplete_intron_retention}:
                        event_count = 1
                score += event_subtype_cost[e.event_type] * event_count
                logger.debug("* * * Event " + str(e.event_type) + ", introns affected " + str(event_count) +
                             ", event cost " + str(event_subtype_cost[e.event_type]) +
                             ". Updated score: " + str(score))
            logger.debug("* * Final score for isoform " + isoform_id + ": " + str(score))
            isoform_scores.append((isoform_id, score))

        best_score = min(isoform_scores, key=lambda x:x[1])[1]
        logger.debug("* * Best score " + str(best_score))
        return [x[0] for x in filter(lambda x:x[1] == best_score, isoform_scores)]

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
        for isoform_match in read_assignment.isoform_matches:
            isoform_id = isoform_match.assigned_transcript
            if isoform_id is None:
                continue
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
                isoform_match.set_classification(MatchClassification.novel_not_in_catalog)

        if apa_found and read_assignment.assignment_type in \
                {ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference}:
            read_assignment.assignment_type = ReadAssignmentType.inconsistent

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

