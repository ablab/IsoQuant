############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from enum import Enum, unique
from collections import namedtuple

from src.common import *

logger = logging.getLogger('IsoQuant')

@unique
class ReadAssignmentType(Enum):
    unique = 1
    noninformative = 0
    ambiguous = 10
    unique_minor_difference = 2
    inconsistent = 3


# SQANTI-like
@unique
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
    def get_inconsistency_classification(match_event_subtypes):
        if any(me.event_type in nnic_event_types for me in match_event_subtypes):
            return MatchClassification.novel_not_in_catalog
        elif any(me.event_type in nic_event_types for me in match_event_subtypes):
            return MatchClassification.novel_in_catalog
        return MatchClassification.undefined

    @staticmethod
    def get_mono_exon_classification(match_event_subtypes):
        # events are not mixed in the list
        if any(e.event_type in {MatchEventSubtype.alternative_polya_site, MatchEventSubtype.fake_polya_site}
               for e in match_event_subtypes):
            return MatchClassification.novel_not_in_catalog
        if any(e.event_type == MatchEventSubtype.unspliced_intron_retention for e in match_event_subtypes):
            return MatchClassification.novel_in_catalog
        elif any(e.event_type in {MatchEventSubtype.incomplete_intron_retention_left,
                                  MatchEventSubtype.incomplete_intron_retention_right} for e in match_event_subtypes):
            return MatchClassification.genic
        elif match_event_subtypes[0].event_type == MatchEventSubtype.mono_exon_match:
            return MatchClassification.mono_exon_match
        elif match_event_subtypes[0].event_type == MatchEventSubtype.mono_exonic:
            return MatchClassification.incomplete_splice_match
        else:
            logger.warning("Unexpected set of monoexonic events: " + str(match_event_subtypes))
            return MatchClassification.undefined


@unique
class MatchEventSubtype(Enum):
    none = 0
    undefined = 1
    # non-contradictory
    mono_exonic = 10
    fsm = 11
    ism_left = 15
    ism_right = 13
    ism_internal = 14
    mono_exon_match = 19
    # alignment artifacts
    intron_shift = 21
    exon_misallignment = 22
    fake_terminal_exon_left = 24
    fake_terminal_exon_right = 26
    # minor alternations
    exon_elongation_left = 25
    exon_elongation_right = 23
    # intron retentions
    intron_retention = 31
    unspliced_intron_retention = 32
    incomplete_intron_retention_left = 38
    incomplete_intron_retention_right = 39
    # major alternation
    # alternative donor/acceptor sites
    alt_left_site_known = 101
    alt_right_site_known = 102
    alt_left_site_novel = 103
    alt_right_site_novel = 104
    # additional introns in the middle
    extra_intron = 1012
    extra_intron_known = 1011
    # extra inrons on the sides
    extra_intron_flanking_left = 1013
    extra_intron_flanking_right = 1015
    # significant exon elongation, more than allowed
    major_exon_elongation_left = 1005
    major_exon_elongation_right = 1003
    # other intron modifications
    intron_migration = 114
    intron_alternation_novel = 115
    intron_alternation_known = 116
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
    fake_polya_site = 201
    alternative_tss = 202

    def __lt__(self, other):
        return self.value < other.value

    @staticmethod
    def is_alignment_artifact(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.intron_shift, MatchEventSubtype.exon_misallignment,
                                       MatchEventSubtype.fake_terminal_exon_left,
                                       MatchEventSubtype.fake_terminal_exon_right}

    @staticmethod
    def is_minor_error(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.exon_elongation_left,
                                       MatchEventSubtype.exon_elongation_right,
                                       MatchEventSubtype.intron_shift,
                                       MatchEventSubtype.exon_misallignment,
                                       MatchEventSubtype.fake_terminal_exon_left,
                                       MatchEventSubtype.fake_terminal_exon_right}

    @staticmethod
    def is_consistent(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.none,
                                       MatchEventSubtype.mono_exonic,
                                       MatchEventSubtype.mono_exon_match}

    @staticmethod
    def is_major_elongation(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.major_exon_elongation_left,
                                       MatchEventSubtype.major_exon_elongation_right,
                                       MatchEventSubtype.incomplete_intron_retention}

    @staticmethod
    def is_major_inconsistency(match_event_subtype):
        return match_event_subtype in nnic_event_types or match_event_subtype in nic_event_types


event_subtype_cost = {
    MatchEventSubtype.none:0,
    MatchEventSubtype.undefined:0,
    MatchEventSubtype.mono_exonic:0,
    MatchEventSubtype.ism_left:0,
    MatchEventSubtype.ism_right:0,
    MatchEventSubtype.ism_internal:0,
    MatchEventSubtype.mono_exon_match:0,
    MatchEventSubtype.intron_shift:0.1,
    MatchEventSubtype.exon_misallignment:0.1,
    MatchEventSubtype.fake_terminal_exon_left:0.2,
    MatchEventSubtype.fake_terminal_exon_right:0.2,
    # minor alternations
    MatchEventSubtype.exon_elongation_left:0.1,
    MatchEventSubtype.exon_elongation_right:0.1,
    # intron retentions
    MatchEventSubtype.intron_retention:0.5,
    MatchEventSubtype.unspliced_intron_retention:0.5,
    MatchEventSubtype.incomplete_intron_retention_left:0.75,
    MatchEventSubtype.incomplete_intron_retention_right:0.75,
    # major alternation
    # alternative donor/acceptor sites
    MatchEventSubtype.alt_left_site_known:1,
    MatchEventSubtype.alt_right_site_known:1,
    MatchEventSubtype.alt_left_site_novel:1,
    MatchEventSubtype.alt_right_site_novel:1,
    # additional introns in the middle
    MatchEventSubtype.extra_intron:1,
    MatchEventSubtype.extra_intron_known:1,
    # extra inrons on the sides
    MatchEventSubtype.extra_intron_flanking_left:1,
    MatchEventSubtype.extra_intron_flanking_right:1,
    # significant exon elongation, more than allowed
    MatchEventSubtype.major_exon_elongation_left:1,
    MatchEventSubtype.major_exon_elongation_right:1,
    # other intron modifications
    MatchEventSubtype.intron_migration:1,
    MatchEventSubtype.intron_alternation_novel:1,
    MatchEventSubtype.intron_alternation_known:1,
    # mutually exclusive
    MatchEventSubtype.mutually_exclusive_exons_novel:1,
    MatchEventSubtype.mutually_exclusive_exons_known:1,
    # exon skipping
    MatchEventSubtype.exon_skipping_known_intron:1,
    MatchEventSubtype.exon_skipping_novel_intron:1,
    # exon gain
    MatchEventSubtype.exon_gain_known:1,
    MatchEventSubtype.exon_gain_novel:1,
    # other
    MatchEventSubtype.alternative_structure_novel:1,
    MatchEventSubtype.alternative_structure_known:1,
    # TTS and TSS
    MatchEventSubtype.alternative_polya_site:0.1,
    MatchEventSubtype.fake_polya_site:0.5,
    MatchEventSubtype.alternative_tss :0.1
}


def elongation_cost(params, elongation_len):
    min_cost = event_subtype_cost[MatchEventSubtype.exon_elongation_left]
    max_cost = event_subtype_cost[MatchEventSubtype.major_exon_elongation_left]
    lower_bound = params.minor_exon_extension
    upper_bound = params.major_exon_extension
    if elongation_len <= lower_bound:
        return min_cost
    elif elongation_len >= upper_bound:
        return max_cost
    else:
        logger.debug(str(elongation_len) + ", " + str(lower_bound) + "," + str(upper_bound) + ", " + str((elongation_len - lower_bound) / (upper_bound - lower_bound)))
        return min_cost + (max_cost - min_cost) * (elongation_len - lower_bound) / (upper_bound - lower_bound)


nnic_event_types = {
    MatchEventSubtype.alt_left_site_novel, MatchEventSubtype.alt_right_site_novel,
    MatchEventSubtype.extra_intron, MatchEventSubtype.extra_intron_flanking_left,
    MatchEventSubtype.extra_intron_flanking_right, MatchEventSubtype.mutually_exclusive_exons_novel,
    MatchEventSubtype.exon_gain_novel, MatchEventSubtype.exon_skipping_novel_intron,
    MatchEventSubtype.alternative_structure_novel, MatchEventSubtype.intron_alternation_novel,
    MatchEventSubtype.alternative_polya_site, MatchEventSubtype.alternative_tss, MatchEventSubtype.fake_polya_site
}

nic_event_types = {
    MatchEventSubtype.unspliced_intron_retention, MatchEventSubtype.intron_retention,
    MatchEventSubtype.alt_left_site_known, MatchEventSubtype.alt_right_site_known,
    MatchEventSubtype.extra_intron_known, MatchEventSubtype.intron_migration,
    MatchEventSubtype.mutually_exclusive_exons_known, MatchEventSubtype.exon_skipping_known_intron,
    MatchEventSubtype.exon_gain_known, MatchEventSubtype.alternative_structure_known,
    MatchEventSubtype.intron_alternation_known, MatchEventSubtype.major_exon_elongation_left,
    MatchEventSubtype.major_exon_elongation_right, MatchEventSubtype.incomplete_intron_retention_left,
    MatchEventSubtype.incomplete_intron_retention_right
}

nonintronic_events = {
    MatchEventSubtype.alternative_polya_site, MatchEventSubtype.alternative_tss, MatchEventSubtype.fake_polya_site,
    MatchEventSubtype.major_exon_elongation_left, MatchEventSubtype.major_exon_elongation_right,
    MatchEventSubtype.exon_elongation_left, MatchEventSubtype.exon_elongation_right,
}


# (side, is_known) -> alternation type
alternative_sites = {("left", True): MatchEventSubtype.alt_left_site_known,
                     ("left", False): MatchEventSubtype.alt_left_site_novel,
                     ("right", True): MatchEventSubtype.alt_right_site_known,
                     ("right", False): MatchEventSubtype.alt_right_site_novel}


class SupplementaryMatchConstansts:
    extra_left_mod_position = -1000000
    extra_right_mod_position = 1000000
    undefined_position = -2000000
    undefined_region = (undefined_position, undefined_position)


MatchEvent = namedtuple("MatchEvent", ("event_type", "isoform_position", "read_region", "event_length"))


def make_event(event_type,
               isoform_position=SupplementaryMatchConstansts.undefined_position,
               read_region=SupplementaryMatchConstansts.undefined_region,
               event_length=0):
    return MatchEvent(event_type, isoform_position, read_region, event_length)


class IsoformMatch:
    def __init__(self, match_classification, assigned_gene=None, assigned_transcript=None,
                 match_subclassification = None, transcript_strand=None):
        self.assigned_gene = assigned_gene
        self.assigned_transcript = assigned_transcript
        self.transcript_strand = transcript_strand
        self.match_classification = match_classification
        if match_subclassification is None:
            self.match_subclassifications = []
        elif isinstance(match_subclassification, list):
            self.match_subclassifications = \
                list(filter(lambda x: x.event_type != MatchEventSubtype.none, match_subclassification))
        else:
            self.match_subclassifications = [match_subclassification]

    def add_subclassification(self, match_subclassification):
        if len(self.match_subclassifications) == 1 and \
                self.match_subclassifications[0].event_type in {MatchEventSubtype.undefined, MatchEventSubtype.none}:
            self.match_subclassifications = [match_subclassification]
        else:
            self.match_subclassifications.append(match_subclassification)

    def set_classification(self, classification):
        self.match_classification = classification

    def monoexon_is_consistent(self):
        valid_subtypes = {MatchEventSubtype.none, MatchEventSubtype.mono_exonic, MatchEventSubtype.mono_exon_match}
        return all(el.event_type in valid_subtypes for el in self.match_subclassifications)

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
