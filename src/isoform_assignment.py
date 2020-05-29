############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import copy

from src.common import *


logger = logging.getLogger('IsoQuant')


class ReadAssignmentType:
    unique = "unique"
    empty = "empty"
    ambiguous = "ambiguous"
    minor = "unique_minor_difference"
    contradictory = "unique_contradictory"
    novel = "novel"


# SQANTI-like
class MatchClassification:
    none = "undefined"
    fsm = "full_splice_match"
    ism = "incomplete_splice_match"
    nic = "novel_in_catalog"
    nnic = "novel_not_in_catalog"
    genic = "genic"
    antisense = "antisense"
    fusion = "fusion"
    intergenic = "intergenic"
    genic_intron = "genic_intron"

    @staticmethod
    def get_contradiction_classification_from_subtypes(match_event_subtypes):
        if any(me in [MatchEventSubtype.alt_donor_site_novel, MatchEventSubtype.alt_acceptor_site_novel,
                      MatchEventSubtype.extra_intron, MatchEventSubtype.extra_intron_out,
                      MatchEventSubtype.mutually_exclusive_exons_novel, MatchEventSubtype.exon_gain_novel,
                      MatchEventSubtype.exon_skipping_novel_intron, MatchEventSubtype.alternative_structure_novel]
               for me in match_event_subtypes):
            return MatchClassification.nnic
        elif any(me in [MatchEventSubtype.unspliced_intron_retention, MatchEventSubtype.intron_retention,
                        MatchEventSubtype.alt_donor_site, MatchEventSubtype.alt_acceptor_site,
                        MatchEventSubtype.extra_intron_known, MatchEventSubtype.intron_migration,
                        MatchEventSubtype.mutually_exclusive_exons_known,
                        MatchEventSubtype.exon_skipping_known_intron, MatchEventSubtype.exon_gain_known,
                        MatchEventSubtype.alternative_structure_known]
                 for me in match_event_subtypes):
            return MatchClassification.nic
        elif any(me in [MatchEventSubtype.unspliced_genic] for me in match_event_subtypes):
            return MatchClassification.genic
        return MatchClassification.none

class MatchEventSubtype:
    none = "none"
    undefined = "undefined"

    # non-contradictory
    unspliced = "mono_exon"
    ism_5 = "ism_5"
    ism_3 = "ism_3"
    ism_internal = "ism_internal"
    # alignment artifacts
    intron_shift = "intron_shift"
    exon_misallignment = "small_exon_misallignment"
    # minor alternations
    exon_elongation5 = "exon_elongation_5prime"
    exon_elongation3 = "exon_elongation_3prime"
    exon_elongation_both = "exon_elongation_both"
    # intron retentions
    intron_retention = "intron_retention"
    unspliced_intron_retention = "mono_exon_intron_retention"
    unspliced_genic = "mono_exon_genic"
    # major alternation
    alt_donor_site = "alt_donor_site"
    alt_acceptor_site = "alt_acceptor_site"
    alt_donor_site_novel = "alt_donor_site_novel"
    alt_acceptor_site_novel = "alt_acceptor_site_novel"
    extra_intron = "additional_novel_intron"
    extra_intron_known = "additional_known_intron"
    extra_intron_out = "additional_terminal_intron"
    extra_exon_out = "additional_terminal_exon"
    intron_migration = "intron_migration"
    intron_alternation_novel = "intron_change_to_novel"
    intron_alternation_known = "intron_change_to_known"
    mutually_exclusive_exons_novel = "mutualy_exclusive_novel_exons"
    mutually_exclusive_exons_known = "mutualy_exclusive_known_exons"
    exon_skipping_known_intron = "exon_skipping_known_intron"
    exon_skipping_novel_intron = "exon_skipping_novel_intron"
    exon_gain_known = "gains_known_exon"
    exon_gain_novel = "gains_novel_exon"
    alternative_structure_novel = "alternative_structure_novel_introns"
    alternative_structure_known = "alternative_structure_known_introns"

    @staticmethod
    def is_alignment_artifact(match_event_subtype):
        return match_event_subtype in [MatchEventSubtype.intron_shift, MatchEventSubtype.exon_misallignment]

    @staticmethod
    def is_minor_error(match_event_subtype):
        return match_event_subtype in [MatchEventSubtype.exon_elongation5, MatchEventSubtype.exon_elongation3,
                                       MatchEventSubtype.exon_elongation_both]


class IsoformMatch:
    def __init__(self, match_classification, assigned_gene = "None", assigned_transcript = "None",
                 match_subclassifications = None):
        self.assigned_gene = assigned_gene
        self.assigned_transcript = assigned_transcript
        self.match_classification = copy.deepcopy(match_classification)
        if match_subclassifications is None:
            self.match_subclassifications = []
        elif isinstance(match_subclassifications, list):
            self.match_subclassifications = match_subclassifications
        else:
            self.match_subclassifications = [copy.deepcopy(match_subclassifications)]
        self.additional_info = {}

    def add_subclassification(self, match_subclassification):
        if len(self.match_subclassifications) == 1 and self.match_subclassifications[0] == MatchClassification.none:
            self.match_subclassifications = [copy.deepcopy(match_subclassification)]
        else:
            self.match_subclassifications.append(copy.deepcopy(match_subclassification))

    def set_classification(self, classification):
        self.match_classification = copy.deepcopy(classification)

    def all_subtypes_are_none(self):
        return all(el == MatchEventSubtype.none for el in self.match_subclassifications)

    def all_subtypes_are_alignment_artifacts(self):
        return all(MatchEventSubtype.is_alignment_artifact(el) for el in self.match_subclassifications)

    def all_subtypes_are_minor_errors(self):
        return all(MatchEventSubtype.is_minor_error(el) for el in self.match_subclassifications)

    def set_additional_info(self, key, value):
        self.additional_info[key] = value


class ReadAssignment:
    def __init__(self, read_id, assignment_type, match = None):
        self.read_id = read_id
        self.assignment_type = copy.deepcopy(assignment_type)
        if match is None:
            self.isoform_matches = []
        elif isinstance(match, list):
            self.isoform_matches = match
        else:
            self.isoform_matches = [match]

    def add_match(self, match):
        self.isoform_matches.append(match)

    def set_assignment_type(self, assignment_type):
        self.assignment_type = copy.deepcopy(assignment_type)

    def merge(self, other):
        pass