############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import copy

from src.common import *
from enum import Enum

logger = logging.getLogger('IsoQuant')


class ReadAssignmentType(Enum):
    unique = 1
    empty = 0
    ambiguous = 10
    unique_minor_difference = 2
    contradictory = 3
    novel = 4


# SQANTI-like
class MatchClassification(Enum):
    undefined = 0
    full_splice_match = 10
    incomplete_splice_match = 11
    novel_in_catalog = 20
    novel_not_in_catalog = 21
    genic = 30
    antisense = 40
    fusion = 50
    intergenic = 32
    genic_intron = 31

    @staticmethod
    def get_contradiction_classification_from_subtypes(match_event_subtypes):
        if any(me in [MatchEventSubtype.alt_donor_site_novel, MatchEventSubtype.alt_acceptor_site_novel,
                      MatchEventSubtype.extra_intron, MatchEventSubtype.extra_intron_out,
                      MatchEventSubtype.mutually_exclusive_exons_novel, MatchEventSubtype.exon_gain_novel,
                      MatchEventSubtype.exon_skipping_novel_intron, MatchEventSubtype.alternative_structure_novel]
               for me in match_event_subtypes):
            return MatchClassification.novel_not_in_catalog
        elif any(me in [MatchEventSubtype.unspliced_intron_retention, MatchEventSubtype.intron_retention,
                        MatchEventSubtype.alt_donor_site, MatchEventSubtype.alt_acceptor_site,
                        MatchEventSubtype.extra_intron_known, MatchEventSubtype.intron_migration,
                        MatchEventSubtype.mutually_exclusive_exons_known,
                        MatchEventSubtype.exon_skipping_known_intron, MatchEventSubtype.exon_gain_known,
                        MatchEventSubtype.alternative_structure_known]
                 for me in match_event_subtypes):
            return MatchClassification.novel_in_catalog
        elif any(me in [MatchEventSubtype.unspliced_genic] for me in match_event_subtypes):
            return MatchClassification.genic
        return MatchClassification.undefined

class MatchEventSubtype(Enum):
    none = 0
    undefined = 1
    # non-contradictory
    mono_exonic = 10
    ism_5 = 15
    ism_3 = 13
    ism_internal = 14
    # alignment artifacts
    intron_shift = 21
    exon_misallignment = 22
    # minor alternations
    exon_elongation5 = 25
    exon_elongation3 = 23
    exon_elongation_both = 24
    # intron retentions
    intron_retention = 31
    unspliced_intron_retention = 32
    unspliced_genic = 33
    # major alternation
    alt_donor_site = 101
    alt_acceptor_site = 102
    alt_donor_site_novel = 103
    alt_acceptor_site_novel = 104
    extra_intron = 112
    extra_intron_known = 111
    extra_intron_out = 113
    extra_exon_out = 120
    intron_migration = 114
    intron_alternation_novel = 115
    intron_alternation_known = 115
    mutually_exclusive_exons_novel = 121
    mutually_exclusive_exons_known = 122
    exon_skipping_known_intron = 123
    exon_skipping_novel_intron = 124
    exon_gain_known = 125
    exon_gain_novel = 126
    alternative_structure_novel = 131
    alternative_structure_known = 132

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
        self.match_classification = match_classification
        if match_subclassifications is None:
            self.match_subclassifications = []
        elif isinstance(match_subclassifications, list):
            self.match_subclassifications = match_subclassifications
        else:
            self.match_subclassifications = [match_subclassifications]
        self.additional_info = {}

    def add_subclassification(self, match_subclassification):
        if len(self.match_subclassifications) == 1 and self.match_subclassifications[0] == MatchClassification.undefined:
            self.match_subclassifications = [match_subclassification]
        else:
            self.match_subclassifications.append(match_subclassification)

    def set_classification(self, classification):
        self.match_classification = classification

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
        self.assignment_type = assignment_type
        if match is None:
            self.isoform_matches = []
        elif isinstance(match, list):
            self.isoform_matches = match
        else:
            self.isoform_matches = [match]

    def add_match(self, match):
        self.isoform_matches.append(match)

    def set_assignment_type(self, assignment_type):
        self.assignment_type = assignment_type

    def merge(self, other):
        pass