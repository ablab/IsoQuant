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
from src.polya_verification import *
from src.long_read_assigner import IsoformDiff
from src.long_read_assigner import AmbiguityResolvingMethod

logger = logging.getLogger('IsoQuant')


class LongReadSimpleAssigner:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params
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
            isoform_region = (self.gene_info.transcript_start(isoform_id) - self.params.min_abs_exon_overlap,
                              self.gene_info.transcript_end(isoform_id) + self.params.min_abs_exon_overlap)
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
        return [x[0] for x in filter(lambda x: x[1] == 0, isoforms)]


    # select most similar isoform based on different criteria
    def select_similar_isoforms(self, combined_read_profile):
        read_split_exon_profile = combined_read_profile.read_split_exon_profile
        overlapping_isoforms = self.find_overlapping_isoforms(read_split_exon_profile,
                                                              self.gene_info.split_exon_profiles.profiles)
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
        best_candidates = [x[0] for x in filter(lambda x: x[1] <= best_diff + 3, candidates)]
        logger.debug("+ + Closest matching isoforms " + str(best_candidates))

        return best_candidates

    # detect exon elongation subtyp
    def categorize_exon_elongation_subtype(self, read_split_exon_profile, isoform_id):
        split_exons = self.gene_info.split_exon_profiles.features
        isoform_profile = self.gene_info.split_exon_profiles.profiles[isoform_id]

        # find first and last common exons
        common_first_exon = -1
        isofrom_first_exon = isoform_profile.index(1)
        for i in range(len(split_exons)):
            if isoform_profile[i] == read_split_exon_profile.gene_profile[i] == 1:
                common_first_exon = i
                break

        common_last_exon = -1
        isofrom_last_exon = rindex(isoform_profile, 1)
        for i in range(len(split_exons)):
            index = len(split_exons) - i - 1
            if isoform_profile[index] == read_split_exon_profile.gene_profile[index] == 1:
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
        left_event = None
        if extra_left > self.params.minor_exon_extension:
            if common_first_exon == isofrom_first_exon:
                left_event = MatchEventSubtype.major_exon_elongation_left
        elif extra_left > self.params.delta:
            left_event = MatchEventSubtype.exon_elongation_left
        if left_event:
            events.append(MatchEvent(left_event, event_info=extra_left))

        right_event = None
        if extra_right > self.params.minor_exon_extension:
            if common_last_exon == isofrom_last_exon:
                right_event = MatchEventSubtype.major_exon_elongation_right
        elif extra_right > self.params.delta:
                right_event = MatchEventSubtype.exon_elongation_right
        if right_event:
            events.append(MatchEvent(right_event, event_info=extra_right))

        return events

    # select best assignment based on nucleotide similarity
    # score = similarity_function(read_exons, isoform_exons)
    def resolve_by_nucleotide_score(self, combined_read_profile, matched_isoforms, similarity_function,
                                    min_similarity=AmbiguityResolvingMethod.minimal_score.value,
                                    top_scored_factor=AmbiguityResolvingMethod.top_scored_factor.value):
        if not matched_isoforms:
            return []

        logger.debug("+ + + Resolving by nucleotide similarity")
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
        logger.debug("+ + + Best score = %f, all candidates = %s" % (best_score, str(top_scored)))
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
    def check_for_extra_terminal_seqs(self, read_split_exon_profile, assignment):
        for match in assignment.isoform_matches:
            if match.assigned_transcript is None:
                continue

            exon_elongation_types = self.categorize_exon_elongation_subtype(read_split_exon_profile,
                                                                            match.assigned_transcript)
            for e in exon_elongation_types:
                match.add_subclassification(e)

            if  any(MatchEventSubtype.is_major_elongation(e.event_type) for e in exon_elongation_types):
                # serious exon elongation
                assignment.set_assignment_type(ReadAssignmentType.inconsistent)
            elif len(exon_elongation_types) > 0:
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
            logger.debug(" + Mono exon")
            match_classification = MatchClassification.mono_exon_match
            match_subclassifications = MatchEvent(MatchEventSubtype.mono_exon_match)
        elif self.is_fsm(read_region, isoform_id):
            logger.debug("+ + Full splice match " + isoform_id)
            match_classification = MatchClassification.full_splice_match
            match_subclassifications = MatchEvent(MatchEventSubtype.fsm)
        else:
            logger.debug("+ + Incomplete splice match " + isoform_id)
            match_classification = MatchClassification.incomplete_splice_match
            match_subclassifications = self.detect_ism_subtype(read_region, isoform_id)

        return IsoformMatch(match_classification, self.get_gene_id(isoform_id), isoform_id,
                            match_subclassifications, self.gene_info.isoform_strands[isoform_id])

    # make proper match subtype
    def categorize_multiple_splice_matches(self, combined_read_profile, isoform_ids):
        isoform_matches = []
        for isoform_id in isoform_ids:
            isoform_matches.append(self.categorize_correct_splice_match(combined_read_profile, isoform_id))
        return isoform_matches

    def categorize_correct_unspliced_match(self, combined_read_profile, isoform_id):
        isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
        if len(isoform_introns) == 0:
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
        for isoform_id in isoform_ids:
            isoform_matches.append(self.categorize_correct_unspliced_match(combined_read_profile, isoform_id))
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
        if len(read_intron_profile.read_profile) == 0:
                # or any(el == -1 for el in read_split_exon_profile.read_profile):
            # Read has inconsistent exons / introns
            assignment = ReadAssignment(read_id, ReadAssignmentType.mono_exon_skipped)
        elif all(el != 1 for el in read_split_exon_profile.read_profile) \
                or all(el == 0 or el == -2 for el in read_split_exon_profile.gene_profile):
            read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[-1][1])
            gene_region = (self.gene_info.split_exon_profiles.features[0][0],
                              self.gene_info.split_exon_profiles.features[-1][1])
            # none of the blocks matched
            if not overlaps(read_region, gene_region):
                logger.debug("EMPTY - noninformative")
                assignment = ReadAssignment(read_id, ReadAssignmentType.noninformative, IsoformMatch(MatchClassification.intergenic))
            elif all(el != 1 for el in read_split_exon_profile.gene_profile):
                logger.debug("EMPTY - intronic")
                assignment = ReadAssignment(read_id, ReadAssignmentType.noninformative, IsoformMatch(MatchClassification.genic_intron))
            else:
                logger.debug("EMPTY - genic")
                assignment = ReadAssignment(read_id, ReadAssignmentType.noninformative, IsoformMatch(MatchClassification.genic))
            return assignment

        elif any(el == -1 for el in read_intron_profile.read_profile):
                # or any(el == -1 for el in read_split_exon_profile.read_profile):
            # Read has inconsistent exons / introns
            if all(el != 0 for el in read_split_exon_profile.gene_profile):
                logger.debug("EMPTY - intronic")
                assignment = ReadAssignment(read_id, ReadAssignmentType.noninformative,
                                            IsoformMatch(MatchClassification.genic_intron))
            else:
                logger.debug("+ Has inconsistent features (novel intron / exons)")
                assignment = ReadAssignment(read_id, ReadAssignmentType.inconsistent)

        elif any(el == 0 for el in read_intron_profile.read_profile):
                #or any(el == 0 for el in read_split_exon_profile.read_profile):
            # Read has extra flanking exons / introns
            logger.debug("+ Has extra flanking features")
            assignment = ReadAssignment(read_id, ReadAssignmentType.inconsistent)

        else:
            logger.debug("+ No contradictory features in read, but no consistent isoforms still can be found")
            consistent_isoforms = self.find_matching_isoforms(read_intron_profile.gene_profile,
                                                              self.gene_info.intron_profiles.profiles)
            matches = self.categorize_multiple_splice_matches(combined_read_profile, consistent_isoforms)

            if len(consistent_isoforms) == 0:
                assignment = ReadAssignment(read_id, ReadAssignmentType.inconsistent)
            elif len(consistent_isoforms) == 1:
                assignment = ReadAssignment(read_id, ReadAssignmentType.unique, matches)
            else:
                assignment = ReadAssignment(read_id, ReadAssignmentType.ambiguous, matches)


            # assignment = self.match_consistent(read_id, combined_read_profile)

            # if assignment is None:
                # alternative isoforms made of known introns/exons or intron retention
            #    logger.debug("+ + Resolving unmatched ")
            #    assignment = self.match_inconsistent(read_id, combined_read_profile)

        return assignment

    def classify_assignment(self, best_isoforms, read_matches):
        is_abmiguous = len(best_isoforms) > 1
        all_event_types = set()
        for isoform_id in best_isoforms:
            for e in read_matches[isoform_id]:
                all_event_types.add(e.event_type)

        logger.debug("* All events: " + str(all_event_types))
        if all(MatchEventSubtype.is_consistent(e) for e in all_event_types):
            logger.debug("* * Assignment seems to be consistent")
            assignment_type = ReadAssignmentType.ambiguous if is_abmiguous else ReadAssignmentType.unique
        elif any(MatchEventSubtype.is_major_inconsistency(e) for e in all_event_types):
            logger.debug("* * Assignment is inconsistent")
            assignment_type = ReadAssignmentType.inconsistent
        elif any(MatchEventSubtype.is_minor_error(e) for e in all_event_types):
            logger.debug("* * Assignment has minor errors")
            assignment_type = ReadAssignmentType.ambiguous if is_abmiguous else ReadAssignmentType.unique_minor_difference
        else:
            logger.warning("Unexpected event reported: " + str(all_event_types))
            assignment_type = ReadAssignmentType.noninformative

        return assignment_type
