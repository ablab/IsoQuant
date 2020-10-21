############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from enum import Enum
from collections import namedtuple

from src.common import *

logger = logging.getLogger('IsoQuant')


class ReadAssignmentType(Enum):
    unique = 1
    empty = 0
    ambiguous = 10
    unique_minor_difference = 2
    contradictory = 3
    novel = 4
    contradictory_monoexon = 11


# SQANTI-like
class MatchClassification(Enum):
    undefined = 0
    full_splice_match = 10
    incomplete_splice_match = 11
    mono_exon_match = 12
    novel_in_catalog = 20
    novel_not_in_catalog = 21
    genic = 30
    antisense = 40
    fusion = 50
    intergenic = 32
    genic_intron = 31

    @staticmethod
    def get_contradiction_classification_from_subtypes(match_event_subtypes):
        if any(me.event_type in nnic_event_types for me in match_event_subtypes):
            return MatchClassification.novel_not_in_catalog
        elif any(me.event_type in nic_event_types for me in match_event_subtypes):
            return MatchClassification.novel_in_catalog
        elif any(me.event_type == MatchEventSubtype.unspliced_genic for me in match_event_subtypes):
            return MatchClassification.genic
        return MatchClassification.undefined

    @staticmethod
    def get_mono_exon_classification_from_subtypes(match_event):
        if match_event.event_type == MatchEventSubtype.unspliced_genic:
            return MatchClassification.genic
        elif match_event.event_type == MatchEventSubtype.unspliced_intron_retention:
            return MatchClassification.novel_in_catalog
        elif match_event.event_type == MatchEventSubtype.mono_exon_match:
            return MatchClassification.mono_exon_match
        elif match_event.event_type == MatchEventSubtype.mono_exonic:
            return MatchClassification.incomplete_splice_match
        else:
            assert False


class MatchEventSubtype(Enum):
    none = 0
    undefined = 1
    # non-contradictory
    mono_exonic = 10
    fsm = 11
    ism_5 = 15
    ism_3 = 13
    ism_internal = 14
    mono_exon_match = 15
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
    # alternative donor/acceptor sites
    alt_donor_site_known = 101
    alt_acceptor_site_known = 102
    alt_donor_site_novel = 103
    alt_acceptor_site_novel = 104
    # additional introns in the middle
    extra_intron = 1012
    extra_intron_known = 1011
    # extra inrons on the sides
    extra_intron_flanking_left = 1013
    extra_intron_flanking_right = 1015
    # significant exon elongation, more than allowed
    major_exon_elongation_both = 1004
    major_exon_elongation_left = 1005
    major_exon_elongation_right = 1003
    # other intron modifications
    intron_migration = 114
    intron_alternation_novel = 115
    intron_alternation_known = 115
    # mutually exclusive
    mutually_exclusive_exons_novel = 121
    mutually_exclusive_exons_known = 122
    # exon skipping
    exon_skipping_known_intron = 123
    exon_skipping_novel_intron = 124
    # exon gain
    exon_gain_known = 125
    exon_gain_novel = 126
    # other
    alternative_structure_novel = 131
    alternative_structure_known = 132
    # TTS and TSS
    alternative_polya_site = 200
    alternative_tss = 201

    def __lt__(self, other):
        return self.value < other.value

    @staticmethod
    def is_alignment_artifact(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.intron_shift, MatchEventSubtype.exon_misallignment}

    @staticmethod
    def is_minor_error(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.exon_elongation5, MatchEventSubtype.exon_elongation3,
                                       MatchEventSubtype.exon_elongation_both}

    @staticmethod
    def is_major_elongation(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.major_exon_elongation_left,
                                       MatchEventSubtype.major_exon_elongation_right,
                                       MatchEventSubtype.major_exon_elongation_both}


nnic_event_types = {
    MatchEventSubtype.alt_donor_site_novel, MatchEventSubtype.alt_acceptor_site_novel,
    MatchEventSubtype.extra_intron, MatchEventSubtype.extra_intron_flanking_left,
    MatchEventSubtype.extra_intron_flanking_right,MatchEventSubtype.mutually_exclusive_exons_novel,
    MatchEventSubtype.exon_gain_novel, MatchEventSubtype.exon_skipping_novel_intron,
    MatchEventSubtype.alternative_structure_novel, MatchEventSubtype.intron_alternation_novel
}

nic_event_types = {
    MatchEventSubtype.unspliced_intron_retention, MatchEventSubtype.intron_retention,
    MatchEventSubtype.alt_donor_site_known, MatchEventSubtype.alt_acceptor_site_known,
    MatchEventSubtype.extra_intron_known, MatchEventSubtype.intron_migration,
    MatchEventSubtype.mutually_exclusive_exons_known, MatchEventSubtype.exon_skipping_known_intron,
    MatchEventSubtype.exon_gain_known, MatchEventSubtype.alternative_structure_known,
    MatchEventSubtype.intron_alternation_known
}


class EventSide(Enum):
    none = 0
    right = 3
    left = 5
    both = 35


elongation_types = {"major": {EventSide.right: MatchEventSubtype.major_exon_elongation_right,
                              EventSide.left: MatchEventSubtype.major_exon_elongation_left,
                              EventSide.both: MatchEventSubtype.major_exon_elongation_both},
                    "minor": {EventSide.right: MatchEventSubtype.exon_elongation3,
                              EventSide.left: MatchEventSubtype.exon_elongation5,
                              EventSide.both: MatchEventSubtype.exon_elongation_both}}

class SupplementaryMatchConstansts:
    extra_left_mod_position = -1000000
    extra_right_mod_position = 1000000
    undefined_position = -2000000


MatchEvent = namedtuple("MatchEvent", ("event_type", "isoform_position", "read_region"))


def make_event(event_type, isoform_position=SupplementaryMatchConstansts.undefined_position,
               read_region=SupplementaryMatchConstansts.undefined_position):
    return MatchEvent(event_type, isoform_position, read_region)


class IsoformMatch:
    def __init__(self, match_classification, assigned_gene=None, assigned_transcript=None,
                 match_subclassification = None):
        self.assigned_gene = assigned_gene
        self.assigned_transcript = assigned_transcript
        self.match_classification = match_classification
        if match_subclassification is None:
            self.match_subclassifications = []
        else:
            self.match_subclassifications = [match_subclassification]

    def add_subclassification(self, match_subclassification):
        if len(self.match_subclassifications) == 1 and \
                self.match_subclassifications[0].event_type == MatchClassification.undefined:
            self.match_subclassifications = [match_subclassification]
        else:
            self.match_subclassifications.append(match_subclassification)

    def set_classification(self, classification):
        self.match_classification = classification

    def all_subtypes_are_none_or_monoexonic(self):
        valid_subtypes = [MatchEventSubtype.none, MatchEventSubtype.mono_exonic, MatchEventSubtype.mono_exon_match]
        return all(el.event_type in valid_subtypes for el in self.match_subclassifications)

    def all_subtypes_are_alignment_artifacts(self):
        return all(MatchEventSubtype.is_alignment_artifact(el.event_type) for el in self.match_subclassifications)

    def all_subtypes_are_minor_errors(self):
        return all(MatchEventSubtype.is_minor_error(el.event_type) for el in self.match_subclassifications)


class ReadAssignment:
    def __init__(self, read_id, assignment_type, match=None):
        self.read_id = read_id
        self.combined_profile = None
        self.gene_info = None
        self.polyA_found = False
        self.cage_found = False
        self.read_group = "."
        self.mapped_strand = "."
        self.assignment_type = assignment_type
        if match is None:
            self.isoform_matches = []
        elif isinstance(match, list):
            self.isoform_matches = match
        else:
            self.isoform_matches = [match]
        self.additional_info = {}

    def add_match(self, match):
        self.isoform_matches.append(match)

    def set_assignment_type(self, assignment_type):
        self.assignment_type = assignment_type

    def chr_id(self):
        return self.gene_info.chr_id

    def start(self):
        return self.combined_profile.read_exon_profile.read_features[0][0]

    def end(self):
        return self.combined_profile.read_exon_profile.read_features[-1][1]

    def length(self):
        return sum([x[1] - x[0] + 1 for x in self.combined_profile.read_exon_profile.read_features])

    def exon_count(self):
        return len(self.combined_profile.read_exon_profile.read_features)

    def set_additional_info(self, key, value):
        self.additional_info[key] = value

    def to_str(self):
        pass

    def from_str(self, string):
        pass


def get_assigned_transcript_id(match):
    return match.assigned_transcript


def get_assigned_gene_id(match):
    return match.assigned_gene
