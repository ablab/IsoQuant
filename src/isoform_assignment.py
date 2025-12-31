############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from enum import Enum, unique
from src.common import junctions_from_blocks
from src.id_policy import SimpleIDDistributor
from src.serialization import *
from src.polya_finder import PolyAInfo

logger = logging.getLogger('IsoQuant')


@unique
class ReadAssignmentType(Enum):
    unique = 1
    noninformative = 0
    intergenic = 20
    ambiguous = 10
    unique_minor_difference = 2
    inconsistent = 30
    inconsistent_non_intronic = 31
    inconsistent_ambiguous = 32
    suspended = 255

    def is_inconsistent(self):
        return self in [ReadAssignmentType.inconsistent,
                        ReadAssignmentType.inconsistent_ambiguous,
                        ReadAssignmentType.inconsistent_non_intronic]

    def is_consistent(self):
        return self in [ReadAssignmentType.unique,
                        ReadAssignmentType.unique_minor_difference,
                        ReadAssignmentType.ambiguous]

    def is_unassigned(self):
        return self in [ReadAssignmentType.noninformative,
                        ReadAssignmentType.intergenic]

    def is_unique(self):
        return self in [ReadAssignmentType.unique_minor_difference,
                        ReadAssignmentType.unique]

    def is_ambiguous(self):
        return self in [ReadAssignmentType.ambiguous,
                        ReadAssignmentType.inconsistent_ambiguous]

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
        if any(e.event_type in {MatchEventSubtype.alternative_polya_site_left, MatchEventSubtype.alternative_polya_site_right,
                                MatchEventSubtype.internal_polya_left, MatchEventSubtype.internal_polya_right}
               for e in match_event_subtypes):
            return MatchClassification.novel_not_in_catalog
        if any(e.event_type == MatchEventSubtype.unspliced_intron_retention for e in match_event_subtypes):
            return MatchClassification.novel_in_catalog
        elif any(e.event_type in {MatchEventSubtype.incomplete_intron_retention_left,
                                  MatchEventSubtype.incomplete_intron_retention_right} for e in match_event_subtypes):
            return MatchClassification.genic
        elif match_event_subtypes[0].event_type == MatchEventSubtype.fake_micro_intron_retention:
            return MatchClassification.incomplete_splice_match
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
    antisense = 9
    fsm = 11
    ism_left = 15
    ism_right = 13
    ism_internal = 14
    mono_exon_match = 19
    # alignment artifacts
    intron_shift = 21
    exon_misalignment = 22
    fake_terminal_exon_left = 24
    fake_terminal_exon_right = 26
    terminal_exon_misalignment_right = 28
    terminal_exon_misalignment_left = 29
    # minor alternations
    exon_elongation_left = 25
    exon_elongation_right = 23
    # intron retentions
    intron_retention = 31
    unspliced_intron_retention = 32
    incomplete_intron_retention_left = 38
    incomplete_intron_retention_right = 39
    fake_micro_intron_retention = 34
    # major alternation
    # alternative donor/acceptor sites
    alt_left_site_known = 101
    alt_right_site_known = 102
    alt_left_site_novel = 103
    alt_right_site_novel = 104
    # additional introns in the middle
    extra_intron_novel = 1012
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
    exon_skipping_known = 133
    exon_skipping_novel = 134
    # similar to mutually exclusive exons but when read exon is attached to another exon
    exon_merge_known = 135
    exon_merge_novel = 136
    # exon gain
    exon_gain_known = 145
    exon_gain_novel = 146
    # similar to mutually exclusive exons but when reference exon is attached to another exon
    exon_detach_known = 147
    exon_detach_novel = 148
    # terminal exon shift
    terminal_exon_shift_known = 151
    terminal_exon_shift_novel = 152
    # other
    alternative_structure_novel = 181
    alternative_structure_known = 182
    # TTS and TSS
    alternative_polya_site_right = 200
    alternative_polya_site_left = 201
    internal_polya_left = 202
    internal_polya_right = 203
    alternative_tss_left = 204
    alternative_tss_right = 205
    correct_polya_site_left = 222
    correct_polya_site_right = 223
    aligned_polya_tail = 300

    terminal_site_match_left = 250
    terminal_site_match_right = 251
    terminal_site_match_left_precise = 252
    terminal_site_match_right_precise = 253

    def __lt__(self, other):
        return self.value < other.value

    @staticmethod
    def is_alignment_artifact(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.intron_shift, MatchEventSubtype.exon_misalignment,
                                       MatchEventSubtype.terminal_exon_misalignment_right,
                                       MatchEventSubtype.terminal_exon_misalignment_left,
                                       MatchEventSubtype.fake_terminal_exon_left,
                                       MatchEventSubtype.fake_terminal_exon_right,
                                       MatchEventSubtype.fake_micro_intron_retention}

    @staticmethod
    def is_minor_error(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.exon_elongation_left,
                                       MatchEventSubtype.exon_elongation_right,
                                       MatchEventSubtype.intron_shift,
                                       MatchEventSubtype.exon_misalignment,
                                       MatchEventSubtype.fake_micro_intron_retention,
                                       MatchEventSubtype.fake_terminal_exon_left,
                                       MatchEventSubtype.fake_terminal_exon_right,
                                       MatchEventSubtype.terminal_exon_misalignment_left,
                                       MatchEventSubtype.terminal_exon_misalignment_right}

    @staticmethod
    def is_consistent(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.none,
                                       MatchEventSubtype.mono_exonic,
                                       MatchEventSubtype.mono_exon_match,
                                       MatchEventSubtype.fsm,
                                       MatchEventSubtype.ism_left,
                                       MatchEventSubtype.ism_right,
                                       MatchEventSubtype.ism_internal,
                                       MatchEventSubtype.correct_polya_site_left,
                                       MatchEventSubtype.correct_polya_site_right,
                                       MatchEventSubtype.terminal_site_match_left,
                                       MatchEventSubtype.terminal_site_match_right,
                                       MatchEventSubtype.terminal_site_match_left_precise,
                                       MatchEventSubtype.terminal_site_match_right_precise}

    @staticmethod
    def is_major_elongation(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.major_exon_elongation_left,
                                       MatchEventSubtype.major_exon_elongation_right}

    @staticmethod
    def is_minor_elongation(match_event_subtype):
        return match_event_subtype in {MatchEventSubtype.exon_elongation_left,
                                       MatchEventSubtype.exon_elongation_right}

    @staticmethod
    def is_major_inconsistency(match_event_subtype):
        return match_event_subtype in all_major_events

    @staticmethod
    def is_intronic_inconsistency(match_event_subtype):
        return match_event_subtype in intronic_major_events


event_subtype_cost = {
    MatchEventSubtype.none:0,
    MatchEventSubtype.undefined:0,
    MatchEventSubtype.mono_exonic:0,
    MatchEventSubtype.fsm:0,
    MatchEventSubtype.ism_left:0,
    MatchEventSubtype.ism_right:0,
    MatchEventSubtype.ism_internal:0,
    MatchEventSubtype.mono_exon_match:0,
    MatchEventSubtype.intron_shift:0.1,
    MatchEventSubtype.exon_misalignment:0.1,
    MatchEventSubtype.terminal_exon_misalignment_right:0.1,
    MatchEventSubtype.terminal_exon_misalignment_left:0.1,
    MatchEventSubtype.fake_terminal_exon_left:0.2,
    MatchEventSubtype.fake_terminal_exon_right:0.2,
    # minor alternations
    MatchEventSubtype.exon_elongation_left:0.1,
    MatchEventSubtype.exon_elongation_right:0.1,
    # intron retentions
    MatchEventSubtype.intron_retention:0.6,
    MatchEventSubtype.unspliced_intron_retention:0.6,
    MatchEventSubtype.incomplete_intron_retention_left:0.7,
    MatchEventSubtype.incomplete_intron_retention_right:0.7,
    MatchEventSubtype.fake_micro_intron_retention:0.1,
    # major alternation
    # alternative donor/acceptor sites
    MatchEventSubtype.alt_left_site_known:1,
    MatchEventSubtype.alt_right_site_known:1,
    MatchEventSubtype.alt_left_site_novel:1,
    MatchEventSubtype.alt_right_site_novel:1,
    # additional introns in the middle
    MatchEventSubtype.extra_intron_novel:1,
    MatchEventSubtype.extra_intron_known:1,
    # extra inrons on the sides
    MatchEventSubtype.extra_intron_flanking_left:1,
    MatchEventSubtype.extra_intron_flanking_right:1,
    # significant exon elongation, more than allowed
    MatchEventSubtype.major_exon_elongation_left:0.6,
    MatchEventSubtype.major_exon_elongation_right:0.6,
    # other intron modifications
    MatchEventSubtype.intron_migration:1,
    MatchEventSubtype.intron_alternation_novel:1,
    MatchEventSubtype.intron_alternation_known:1,
    # mutually exclusive
    MatchEventSubtype.mutually_exclusive_exons_novel:0.8,
    MatchEventSubtype.mutually_exclusive_exons_known:0.8,
    # exon skipping
    MatchEventSubtype.exon_skipping_known:1,
    MatchEventSubtype.exon_skipping_novel:1,
    MatchEventSubtype.exon_merge_known:0.75,
    MatchEventSubtype.exon_merge_novel:0.75,
    # exon gain
    MatchEventSubtype.exon_gain_known:1,
    MatchEventSubtype.exon_gain_novel:1,
    MatchEventSubtype.exon_detach_known:0.75,
    MatchEventSubtype.exon_detach_novel:0.75,
    MatchEventSubtype.terminal_exon_shift_known:0.5,
    MatchEventSubtype.terminal_exon_shift_novel:0.5,
    # other
    MatchEventSubtype.alternative_structure_novel:1,
    MatchEventSubtype.alternative_structure_known:1,
    # TTS and TSS
    MatchEventSubtype.alternative_polya_site_left:0.6,
    MatchEventSubtype.alternative_polya_site_right:0.6,
    MatchEventSubtype.internal_polya_left:0.5,
    MatchEventSubtype.internal_polya_right:0.5,
    MatchEventSubtype.alternative_tss_left :0.6,
    MatchEventSubtype.alternative_tss_right :0.6,
    MatchEventSubtype.correct_polya_site_left:0,
    MatchEventSubtype.correct_polya_site_right:0,
    MatchEventSubtype.aligned_polya_tail:0,
    MatchEventSubtype.terminal_site_match_left:0,
    MatchEventSubtype.terminal_site_match_right:0,
    MatchEventSubtype.terminal_site_match_left_precise:0,
    MatchEventSubtype.terminal_site_match_right_precise:0
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
        return min_cost + (max_cost - min_cost) * (elongation_len - lower_bound) / (upper_bound - lower_bound)


nnic_event_types = {
    MatchEventSubtype.alt_left_site_novel, MatchEventSubtype.alt_right_site_novel,
    MatchEventSubtype.extra_intron_novel, MatchEventSubtype.extra_intron_flanking_left,
    MatchEventSubtype.extra_intron_flanking_right, MatchEventSubtype.mutually_exclusive_exons_novel,
    MatchEventSubtype.exon_gain_novel, MatchEventSubtype.exon_skipping_novel,
    MatchEventSubtype.exon_detach_novel, MatchEventSubtype.exon_merge_novel,
    MatchEventSubtype.terminal_exon_shift_novel,
    MatchEventSubtype.alternative_structure_novel, MatchEventSubtype.intron_alternation_novel,
    MatchEventSubtype.alternative_polya_site_left, MatchEventSubtype.alternative_polya_site_right,
    MatchEventSubtype.alternative_tss_right, MatchEventSubtype.alternative_tss_left
}

nic_event_types = {
    MatchEventSubtype.unspliced_intron_retention, MatchEventSubtype.intron_retention,
    MatchEventSubtype.alt_left_site_known, MatchEventSubtype.alt_right_site_known,
    MatchEventSubtype.extra_intron_known, MatchEventSubtype.intron_migration,
    MatchEventSubtype.mutually_exclusive_exons_known, MatchEventSubtype.exon_skipping_known,
    MatchEventSubtype.exon_detach_known, MatchEventSubtype.exon_merge_known,
    MatchEventSubtype.terminal_exon_shift_known,
    MatchEventSubtype.exon_gain_known, MatchEventSubtype.alternative_structure_known,
    MatchEventSubtype.intron_alternation_known, MatchEventSubtype.major_exon_elongation_left,
    MatchEventSubtype.major_exon_elongation_right, MatchEventSubtype.incomplete_intron_retention_left,
    MatchEventSubtype.incomplete_intron_retention_right, MatchEventSubtype.internal_polya_right,
    MatchEventSubtype.internal_polya_left
}

nonintronic_events = {
    MatchEventSubtype.alternative_polya_site_left, MatchEventSubtype.alternative_polya_site_right,
    MatchEventSubtype.alternative_tss_left, MatchEventSubtype.alternative_tss_right,
    MatchEventSubtype.internal_polya_right, MatchEventSubtype.internal_polya_left,
    MatchEventSubtype.major_exon_elongation_left, MatchEventSubtype.major_exon_elongation_right,
    MatchEventSubtype.exon_elongation_left, MatchEventSubtype.exon_elongation_right,
}

all_major_events = nic_event_types.union(nnic_event_types)

intronic_major_events = all_major_events.difference(nonintronic_events)


# (side, is_known) -> alternation type
alternative_sites = {("left", True): MatchEventSubtype.alt_left_site_known,
                     ("left", False): MatchEventSubtype.alt_left_site_novel,
                     ("right", True): MatchEventSubtype.alt_right_site_known,
                     ("right", False): MatchEventSubtype.alt_right_site_novel}


class SupplementaryMatchConstants:
    extra_left_mod_position = (1 << 30) - 1
    extra_right_mod_position = (1 << 30) + 1
    undefined_position = (1 << 31)
    undefined_region = (undefined_position, undefined_position)
    extra_left_region = (extra_left_mod_position, extra_left_mod_position)
    extra_right_region = (extra_right_mod_position, extra_right_mod_position)
    absent_position = (1 << 31) - 1


class MatchEvent:
    def __init__(self, event_type:MatchEventSubtype,
                 isoform_region:tuple=SupplementaryMatchConstants.undefined_region,
                 read_region:tuple=SupplementaryMatchConstants.undefined_region,
                 event_info=0):
        self.event_type = event_type
        self.isoform_region = isoform_region
        self.read_region = read_region
        self.event_info = event_info

    @classmethod
    def deserialize(cls, infile):
        event = cls.__new__(cls)
        event.event_type = MatchEventSubtype(read_int(infile, SHORT_INT_BYTES))
        # TODO: think about 2/3-byte storage
        event.isoform_region = (read_int(infile), read_int(infile))
        event.read_region = (read_int(infile), read_int(infile))
        event.event_info = read_int_neg(infile)
        return event

    def serialize(self, outfile):
        write_int(self.event_type.value, outfile, SHORT_INT_BYTES)
        write_int(self.isoform_region[0], outfile)
        write_int(self.isoform_region[1], outfile)
        write_int(self.read_region[0], outfile)
        write_int(self.read_region[1], outfile)
        write_int_neg(self.event_info, outfile)

    def __repr__(self):
        return "%s:%s,%s,%s" % (self.event_type.name, str(self.isoform_region),
                                str(self.read_region), str(self.event_info))


class IsoformMatch:
    """
    Represents a match between a read and an isoform.

    Memory Optimization:
        Gene and transcript IDs are stored as integers referencing
        shared string pools (string_pools parameter required).
    """
    def __init__(self, match_classification, string_pools, assigned_gene=None, assigned_transcript=None,
                 match_subclassification = None, transcript_strand='.', penalty_score=0):
        # Store string pools reference (required for memory optimization)
        assert string_pools is not None, "string_pools is required"
        self._string_pools = string_pools

        # Store gene/transcript as integers
        self.assigned_gene_id = string_pools.gene_pool.get_int(assigned_gene) if assigned_gene else None
        self.assigned_transcript_id = string_pools.transcript_pool.get_int(assigned_transcript) if assigned_transcript else None

        self.transcript_strand = transcript_strand  # Keep as char (1 byte)
        self.match_classification = match_classification
        self.penalty_score = penalty_score
        if match_subclassification is None:
            self.match_subclassifications = []
        elif isinstance(match_subclassification, list):
            self.match_subclassifications = \
                list(filter(lambda x: x.event_type != MatchEventSubtype.none, match_subclassification))
        else:
            self.match_subclassifications = [match_subclassification]

    @property
    def assigned_gene(self):
        """Return gene string (backward compatibility)"""
        return self._string_pools.gene_pool.get_str(self.assigned_gene_id) if self.assigned_gene_id is not None else None

    @assigned_gene.setter
    def assigned_gene(self, value):
        """Set gene from string"""
        self.assigned_gene_id = self._string_pools.gene_pool.get_int(value) if value is not None else None

    @property
    def assigned_transcript(self):
        """Return transcript string (backward compatibility)"""
        return self._string_pools.transcript_pool.get_str(self.assigned_transcript_id) if self.assigned_transcript_id is not None else None

    @assigned_transcript.setter
    def assigned_transcript(self, value):
        """Set transcript from string"""
        self.assigned_transcript_id = self._string_pools.transcript_pool.get_int(value) if value is not None else None

    @classmethod
    def deserialize(cls, infile, string_pools):
        assert string_pools is not None, "string_pools is required"
        match = cls.__new__(cls)
        match._string_pools = string_pools
        match.assigned_gene_id = read_int_or_none(infile)
        match.assigned_transcript_id = read_int_or_none(infile)
        match.transcript_strand = read_string(infile)
        match.match_classification = MatchClassification(read_short_int(infile))
        match.penalty_score = float(read_int(infile)) / float(SHORT_FLOAT_MULTIPLIER)
        match.match_subclassifications = read_list(infile, MatchEvent.deserialize)
        return match

    def serialize(self, outfile):
        write_int_or_none(self.assigned_gene_id, outfile)
        write_int_or_none(self.assigned_transcript_id, outfile)
        write_string(self.transcript_strand, outfile)
        write_short_int(self.match_classification.value, outfile)
        write_int(int(self.penalty_score * SHORT_FLOAT_MULTIPLIER), outfile)
        write_list(self.match_subclassifications, outfile, MatchEvent.serialize)

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


class BasicReadAssignment:
    """
    Simplified read assignment for serialization and multimap resolution.

    Memory Optimization:
        chr_id, barcode, umi, gene IDs, and isoform IDs are stored as integers
        referencing shared string pools (string_pools parameter required).
    """
    def __init__(self, read_assignment):
        self.assignment_id = read_assignment.assignment_id
        self.read_id = read_assignment.read_id

        # String interning for memory optimization
        assert hasattr(read_assignment, '_string_pools'), "read_assignment must have string_pools"
        self._string_pools = read_assignment._string_pools
        assert self._string_pools is not None, "string_pools is required"

        # Copy interned fields
        self.chr_id_int = read_assignment.chr_id_int
        self.barcode_id = read_assignment.barcode_id
        self.umi_id = read_assignment.umi_id

        self.start = 0
        self.end = 0
        if read_assignment.exons:
            self.start = read_assignment.exons[0][0]
            self.end = read_assignment.exons[-1][1]
        self.genomic_region = read_assignment.genomic_region
        self.multimapper = read_assignment.multimapper
        self.polyA_found = read_assignment.polyA_found
        self.assignment_type = read_assignment.assignment_type
        self.gene_assignment_type = read_assignment.gene_assignment_type
        self.penalty_score = 0.0

        # Store gene/isoform lists as integer IDs
        self.gene_ids = []
        self.isoform_ids = []

        if read_assignment.isoform_matches:
            gene_set = set()
            isoform_set = set()
            for m in read_assignment.isoform_matches:
                self.penalty_score = min(self.penalty_score, read_assignment.isoform_matches[0].penalty_score)
                if m.assigned_gene:
                    gene_set.add(m.assigned_gene_id)
                if m.assigned_transcript:
                    isoform_set.add(m.assigned_transcript_id)

            self.gene_ids = list(gene_set)
            self.isoform_ids = list(isoform_set)

    @property
    def chr_id(self):
        """Return chromosome string"""
        if self.chr_id_int is not None:
            return self._string_pools.chromosome_pool.get_str(self.chr_id_int)
        # Handle unpickled state where chr_id is stored as string
        if hasattr(self, '_chr_id_str'):
            return self._chr_id_str
        return "."

    @chr_id.setter
    def chr_id(self, value):
        """Set chromosome from string"""
        self.chr_id_int = self._string_pools.chromosome_pool.get_int(value) if value else None

    @property
    def barcode(self):
        """Return barcode string"""
        if self.barcode_id is not None:
            return self._string_pools.barcode_pool.get_str(self.barcode_id)
        return None

    @barcode.setter
    def barcode(self, value):
        """Set barcode from string"""
        self.barcode_id = self._string_pools.barcode_pool.get_int(value) if value else None

    @property
    def umi(self):
        """Return UMI string"""
        if self.umi_id is not None:
            return self._string_pools.umi_pool.get_str(self.umi_id)
        return None

    @umi.setter
    def umi(self, value):
        """Set UMI from string"""
        self.umi_id = self._string_pools.umi_pool.get_int(value) if value else None

    @property
    def genes(self):
        """Return list of gene strings"""
        if self.gene_ids:
            return [self._string_pools.gene_pool.get_str(gid) for gid in self.gene_ids]
        # Handle unpickled state where genes are stored as strings
        if hasattr(self, '_genes'):
            return self._genes
        return []

    @genes.setter
    def genes(self, value):
        """Set genes from list of strings"""
        self.gene_ids = [self._string_pools.gene_pool.get_int(g) for g in value]

    @property
    def isoforms(self):
        """Return list of isoform strings"""
        if self.isoform_ids:
            return [self._string_pools.transcript_pool.get_str(tid) for tid in self.isoform_ids]
        # Handle unpickled state where isoforms are stored as strings
        if hasattr(self, '_isoforms'):
            return self._isoforms
        return []

    @isoforms.setter
    def isoforms(self, value):
        """Set isoforms from list of strings"""
        self.isoform_ids = [self._string_pools.transcript_pool.get_int(t) for t in value]

    def __eq__(self, other):
        if isinstance(other, BasicReadAssignment):
            return (self.read_id == other.read_id and
                    self.chr_id == other.chr_id and
                    self.start == other.start and
                    self.end == other.end and
                    self.isoforms == other.isoforms)
        return False

    def __getstate__(self):
        # When pickling, always use string representation for compatibility
        return (self.assignment_id,
                self.read_id,
                self.chr_id,  # Property will convert from int if needed
                self.start,
                self.end,
                self.genomic_region[0],
                self.genomic_region[1],
                self.multimapper,
                self.polyA_found,
                self.assignment_type.value,
                self.gene_assignment_type.value,
                self.penalty_score,
                self.isoforms,  # Property will convert from int list if needed
                self.genes)  # Property will convert from int list if needed

    def __setstate__(self, state):
        # When unpickling, store as strings (no pools available during pickle)
        self._string_pools = None
        self.assignment_id = state[0]
        self.read_id = state[1]
        self._chr_id_str = state[2]
        self.chr_id_int = None
        self.start = state[3]
        self.end = state[4]
        self.genomic_region = (state[5], state[6])
        self.multimapper = state[7]
        self.polyA_found = state[8]
        self.assignment_type = ReadAssignmentType(state[9])
        self.gene_assignment_type = ReadAssignmentType(state[10])
        self.penalty_score = state[11]
        self._isoforms = state[12]
        self._genes = state[13]
        self.barcode_id = None
        self.umi_id = None
        self.gene_ids = []
        self.isoform_ids = []

    @classmethod
    def deserialize(cls, infile, string_pools):
        assert string_pools is not None, "string_pools is required"
        read_assignment = cls.__new__(cls)
        read_assignment._string_pools = string_pools
        read_assignment.assignment_id = read_int(infile)
        read_assignment.read_id = read_string(infile)
        read_assignment.chr_id_int = read_int_or_none(infile)
        read_assignment.start = read_int(infile)
        read_assignment.end = read_int(infile)
        read_assignment.genomic_region = (read_int(infile), read_int(infile))
        bool_arr = read_bool_array(infile, 2)
        read_assignment.multimapper = bool_arr[0]
        read_assignment.polyA_found = bool_arr[1]
        read_assignment.assignment_type = ReadAssignmentType(read_short_int(infile))
        read_assignment.gene_assignment_type = ReadAssignmentType(read_short_int(infile))
        read_assignment.penalty_score = float(read_int(infile)) / float(SHORT_FLOAT_MULTIPLIER)
        read_assignment.gene_ids = read_list(infile, read_int)
        read_assignment.isoform_ids = read_list(infile, read_int)
        return read_assignment

    @classmethod
    def deserialize_from_read_assignment(cls, infile, string_pools):
        assert string_pools is not None, "string_pools is required"
        read_assignment = cls.__new__(cls)
        read_assignment._string_pools = string_pools
        read_assignment.assignment_id = read_int(infile)
        read_assignment.read_id = read_string(infile)
        read_assignment.genomic_region = (read_int(infile), read_int(infile))
        exons = read_list_of_pairs(infile, read_int)
        read_assignment.start = exons[0][0]
        read_assignment.end = exons[-1][1]
        read_list_of_pairs(infile, read_int)
        bool_arr = read_bool_array(infile, 3)
        read_assignment.multimapper = bool_arr[0]
        read_assignment.polyA_found = bool_arr[1]
        read_int_neg(infile)
        read_int_neg(infile)
        read_int_neg(infile)
        read_int_neg(infile)
        # Read group: stored as strings for cross-worker compatibility
        read_group_strings = read_list(infile, read_string)
        read_assignment.read_group_ids = string_pools.read_group_to_ids(read_group_strings) if read_group_strings else []
        # Barcode/UMI pools are built in sorted order, so IDs are deterministic
        read_assignment.barcode_id = read_int_or_none(infile)
        read_assignment.umi_id = read_int_or_none(infile)
        read_string(infile)
        read_string(infile)
        read_assignment.chr_id_int = read_int_or_none(infile)
        read_short_int(infile)
        read_assignment.assignment_type = ReadAssignmentType(read_short_int(infile))
        read_assignment.gene_assignment_type = ReadAssignmentType(read_short_int(infile))
        read_assignment.penalty_score = 0.0
        isoform_matches = read_list(infile, lambda f: IsoformMatch.deserialize(f, string_pools))
        gene_set = set()
        isoform_set = set()
        for m in isoform_matches:
            read_assignment.penalty_score = min(read_assignment.penalty_score, isoform_matches[0].penalty_score)
            if m.assigned_gene:
                gene_set.add(m.assigned_gene_id)
            if m.assigned_transcript:
                isoform_set.add(m.assigned_transcript_id)
        read_assignment.gene_ids = list(gene_set)
        read_assignment.isoform_ids = list(isoform_set)

        read_dict(infile)
        read_dict(infile)
        read_short_int(infile)
        read_list(infile, read_int_neg)
        read_list(infile, read_int_neg)
        return read_assignment

    def serialize(self, outfile):
        write_int(self.assignment_id, outfile)
        write_string(self.read_id, outfile)
        write_int_or_none(self.chr_id_int, outfile)
        write_int(self.start, outfile)
        write_int(self.end, outfile)
        write_int(self.genomic_region[0], outfile)
        write_int(self.genomic_region[1], outfile)
        write_bool_array([self.multimapper, self.polyA_found], outfile)
        write_short_int(self.assignment_type.value, outfile)
        write_short_int(self.gene_assignment_type.value, outfile)
        write_int(int(self.penalty_score * SHORT_FLOAT_MULTIPLIER), outfile)
        write_list(self.gene_ids, outfile, write_int)
        write_list(self.isoform_ids, outfile, write_int)


class ReadAssignment:
    """
    Complete read assignment with isoform matches and additional metadata.

    Memory Optimization:
        chr_id, barcode, and umi are stored as integers referencing
        shared string pools (string_pools parameter required).
    """
    assignment_id_generator = SimpleIDDistributor()

    def __init__(self, read_id, assignment_type, string_pools, match=None):
        assert string_pools is not None, "string_pools is required"
        self.assignment_id = ReadAssignment.assignment_id_generator.increment()
        self.read_id = read_id
        self.genomic_region = (0, 0)
        self.exons = None
        self.corrected_exons = None
        self.corrected_introns = None
        self.gene_info = None
        self.multimapper = False
        self.polyA_found = False
        self.cage_found = False
        self.polya_info = None

        # String interning for memory optimization
        self._string_pools = string_pools
        self.read_group_ids = []  # List of integer IDs for read groups
        self.barcode_id = None  # Integer ID for cell/spatial barcode
        self.umi_id = None  # Integer ID for unique molecular identifier
        self.chr_id_int = None  # Integer ID for chromosome

        self.mapped_strand = "."  # Keep as single char (1 byte)
        self.strand = "."  # Keep as single char (1 byte)
        self.mapping_quality = 0
        self.assignment_type = assignment_type
        if match is None:
            self.isoform_matches = []
        elif isinstance(match, list):
            self.isoform_matches = match
        else:
            self.isoform_matches = [match]
        if self.assignment_type == ReadAssignmentType.ambiguous:
            assigned_genes = set([m.assigned_gene for m in self.isoform_matches])
            self.gene_assignment_type = ReadAssignmentType.ambiguous if len(
                assigned_genes) > 1 else ReadAssignmentType.unique
        elif self.assignment_type == ReadAssignmentType.inconsistent_ambiguous:
            assigned_genes = set([m.assigned_gene for m in self.isoform_matches])
            self.gene_assignment_type = ReadAssignmentType.inconsistent_ambiguous if len(
                assigned_genes) > 1 else ReadAssignmentType.inconsistent
        else:
            self.gene_assignment_type = self.assignment_type

        self.additional_info = {}
        self.additional_attributes = {}
        self.introns_match = False
        self.exon_gene_profile = []
        self.intron_gene_profile = []

    @property
    def chr_id(self):
        """Return chromosome string"""
        if self.chr_id_int is not None:
            return self._string_pools.chromosome_pool.get_str(self.chr_id_int)
        return "."

    @chr_id.setter
    def chr_id(self, value):
        """Set chromosome from string"""
        self.chr_id_int = self._string_pools.chromosome_pool.get_int(value) if value else None

    @property
    def barcode(self):
        """Return barcode string"""
        if self.barcode_id is not None:
            return self._string_pools.barcode_pool.get_str(self.barcode_id)
        return None

    @barcode.setter
    def barcode(self, value):
        """Set barcode from string"""
        self.barcode_id = self._string_pools.barcode_pool.get_int(value) if value else None

    @property
    def umi(self):
        """Return UMI string"""
        if self.umi_id is not None:
            return self._string_pools.umi_pool.get_str(self.umi_id)
        return None

    @umi.setter
    def umi(self, value):
        """Set UMI from string"""
        self.umi_id = self._string_pools.umi_pool.get_int(value) if value else None

    @property
    def read_group(self):
        """Return read group as list of strings"""
        if not self.read_group_ids:
            return []
        return self._string_pools.read_group_from_ids(self.read_group_ids)

    @read_group.setter
    def read_group(self, value):
        """Set read group from list of strings"""
        if not value:
            self.read_group_ids = []
        else:
            self.read_group_ids = self._string_pools.read_group_to_ids(value)

    @classmethod
    def deserialize(cls, infile, gene_info, string_pools):
        assert string_pools is not None, "string_pools is required"
        read_assignment = cls.__new__(cls)
        read_assignment._string_pools = string_pools
        read_assignment.assignment_id = read_int(infile)
        read_assignment.read_id = read_string(infile)
        read_assignment.genomic_region = (read_int(infile), read_int(infile))
        read_assignment.exons = read_list_of_pairs(infile, read_int)
        read_assignment.corrected_exons = read_list_of_pairs(infile, read_int)
        read_assignment.corrected_introns = junctions_from_blocks(read_assignment.corrected_exons)
        read_assignment.gene_info = gene_info
        bool_arr = read_bool_array(infile, 3)
        read_assignment.multimapper = bool_arr[0]
        read_assignment.polyA_found = bool_arr[1]
        read_assignment.cage_found = bool_arr[2]
        read_assignment.polya_info = PolyAInfo(read_int_neg(infile), read_int_neg(infile), read_int_neg(infile), read_int_neg(infile))
        # Read group: stored as strings for cross-worker compatibility
        read_group_strings = read_list(infile, read_string)
        read_assignment.read_group_ids = string_pools.read_group_to_ids(read_group_strings) if read_group_strings else []
        # Barcode/UMI pools are built in sorted order, so IDs are deterministic
        read_assignment.barcode_id = read_int_or_none(infile)
        read_assignment.umi_id = read_int_or_none(infile)
        read_assignment.mapped_strand = read_string(infile)
        read_assignment.strand = read_string(infile)
        read_assignment.chr_id_int = read_int_or_none(infile)
        read_assignment.mapping_quality = read_short_int(infile)
        read_assignment.assignment_type = ReadAssignmentType(read_short_int(infile))
        read_assignment.gene_assignment_type = ReadAssignmentType(read_short_int(infile))
        read_assignment.isoform_matches = read_list(infile, lambda f: IsoformMatch.deserialize(f, string_pools))
        read_assignment.additional_info = read_dict(infile)
        read_assignment.additional_attributes = read_dict(infile)
        read_assignment.introns_match = bool(read_short_int(infile))
        read_assignment.exon_gene_profile = read_list(infile, read_int_neg)
        read_assignment.intron_gene_profile = read_list(infile, read_int_neg)
        return read_assignment

    def serialize(self, outfile):
        write_int(self.assignment_id, outfile)
        write_string(self.read_id, outfile)
        write_int(self.genomic_region[0], outfile)
        write_int(self.genomic_region[1], outfile)
        write_list_of_pairs(self.exons, outfile, write_int)
        write_list_of_pairs(self.corrected_exons, outfile, write_int)
        write_bool_array([self.multimapper, self.polyA_found, self.cage_found], outfile)
        write_int_neg(self.polya_info.external_polya_pos, outfile)
        write_int_neg(self.polya_info.external_polyt_pos, outfile)
        write_int_neg(self.polya_info.internal_polya_pos, outfile)
        write_int_neg(self.polya_info.internal_polyt_pos, outfile)
        # Serialize read_group as strings for cross-worker compatibility
        # (dynamic pools like BAM tags have worker-specific ID mappings)
        write_list(self.read_group, outfile, write_string)
        # Barcode/UMI pools are built in sorted order, so IDs are deterministic
        write_int_or_none(self.barcode_id, outfile)
        write_int_or_none(self.umi_id, outfile)
        write_string(self.mapped_strand, outfile)
        write_string(self.strand, outfile)
        write_int_or_none(self.chr_id_int, outfile)
        write_short_int(self.mapping_quality, outfile)
        write_short_int(self.assignment_type.value, outfile)
        write_short_int(self.gene_assignment_type.value, outfile)
        write_list(self.isoform_matches, outfile, IsoformMatch.serialize)
        write_dict(self.additional_info, outfile)
        write_dict(self.additional_attributes, outfile)
        write_short_int(int(self.introns_match), outfile)
        write_list(self.exon_gene_profile, outfile, write_int_neg)
        write_list(self.intron_gene_profile, outfile, write_int_neg)

    def add_match(self, match):
        self.isoform_matches.append(match)

    def set_assignment_type(self, assignment_type):
        self.assignment_type = assignment_type

    def start(self):
        return self.exons[0][0]

    def end(self):
        return self.exons[-1][1]

    def length(self):
        return sum([x[1] - x[0] + 1 for x in self.exons])

    def exon_count(self):
        return len(self.exons)

    def set_additional_info(self, key, value):
        self.additional_info[key] = value

    def set_additional_attribute(self, key, value):
        self.additional_attributes[key] = value

    def add_match_attribute(self, match_event):
        for m in self.isoform_matches:
            m.add_subclassification(match_event)


match_subtype_printable_names = \
    {MatchEventSubtype.ism_left : ('ism_5', 'ism_3', 'ism'),
     MatchEventSubtype.ism_right : ('ism_3', 'ism_5', 'ism'),
     MatchEventSubtype.exon_elongation_left : ('exon_elongation_5', 'exon_elongation_3', 'exon_elongation'),
     MatchEventSubtype.exon_elongation_right : ('exon_elongation_3', 'exon_elongation_5', 'exon_elongation'),
     MatchEventSubtype.major_exon_elongation_left: ('major_exon_elongation_5', 'major_exon_elongation_3', 'major_exon_elongation'),
     MatchEventSubtype.major_exon_elongation_right: ('major_exon_elongation_3', 'major_exon_elongation_5', 'major_exon_elongation'),
     MatchEventSubtype.fake_terminal_exon_left : ('fake_terminal_exon_5', 'fake_terminal_exon_3', 'fake_terminal_exon'),
     MatchEventSubtype.fake_terminal_exon_right : ('fake_terminal_exon_3', 'fake_terminal_exon_5', 'fake_terminal_exon'),
     MatchEventSubtype.terminal_exon_misalignment_left : ('terminal_exon_misalignment_5', 'terminal_exon_misalignment_3', 'terminal_exon_misalignment'),
     MatchEventSubtype.terminal_exon_misalignment_right : ('terminal_exon_misalignment_3', 'terminal_exon_misalignment_5', 'terminal_exon_misalignment'),
     MatchEventSubtype.incomplete_intron_retention_left: ('incomplete_intron_retention_5', 'incomplete_intron_retention_3', 'incomplete_intron_retention'),
     MatchEventSubtype.incomplete_intron_retention_right: ('incomplete_intron_retention_3', 'incomplete_intron_retention_5', 'incomplete_intron_retention'),
     MatchEventSubtype.extra_intron_flanking_left: ('extra_intron_5', 'extra_intron_3', 'extra_intron'),
     MatchEventSubtype.extra_intron_flanking_right: ('extra_intron_3', 'extra_intron_5', 'extra_intron'),
     MatchEventSubtype.alt_left_site_known: ('alt_donor_site_known', 'alt_acceptor_site_known', 'alternative_splice_site_known'),
     MatchEventSubtype.alt_right_site_known: ('alt_acceptor_site_known', 'alt_donor_site_known', 'alternative_splice_site_known'),
     MatchEventSubtype.alt_left_site_novel: ('alt_donor_site_novel', 'alt_acceptor_site_novel', 'alternative_splice_site_novel'),
     MatchEventSubtype.alt_right_site_novel: ('alt_acceptor_site_novel', 'alt_donor_site_novel', 'alternative_splice_site_novel'),
     MatchEventSubtype.terminal_site_match_left: ('tss_match', 'tes_match', 'terminal_position_match'),
     MatchEventSubtype.terminal_site_match_right: ('tes_match', 'tss_match', 'terminal_position_match'),
     MatchEventSubtype.terminal_site_match_left_precise: ('tss_match_precise', 'tes_match_precise', 'terminal_position_match_precise'),
     MatchEventSubtype.terminal_site_match_right_precise: ('tes_match_precise', 'tss_match_precise', 'terminal_position_match_precise')
     }
     #MatchEventSubtype.alternative_polya_site_left: ('alternative_polya_site_5', 'alternative_polya_site_3'),
     #MatchEventSubtype.alternative_polya_site_right: ('alternative_polya_site_3', 'alternative_polya_site_5'),
     #MatchEventSubtype.internal_polya_left: ('internal_polya_site_5', 'internal_polya_site_3'),
     #MatchEventSubtype.internal_polya_right: ('internal_polya_site_3', 'internal_polya_site_5'),
     #MatchEventSubtype.correct_polya_site_left: ('correct_polya_site_5', 'correct_polya_site_3'),
     #MatchEventSubtype.correct_polya_site_right: ('correct_polya_site_3', 'correct_polya_site_5'),
     #MatchEventSubtype.alternative_tss_left: ('alternative_tss_5', 'alternative_tss_3'),
     #MatchEventSubtype.alternative_tss_right: ('alternative_tss_3', 'alternative_tss_5')}


def match_subtype_to_str(event, strand):
    event_subtype = event.event_type
    if event_subtype in match_subtype_printable_names.keys():
        if strand is None:
            logger.warning("Strand is not set for site-dependent modifications")
        if strand == '-':
            return match_subtype_printable_names[event_subtype][1]
        elif strand == '+':
            return match_subtype_printable_names[event_subtype][0]
        else:
            return match_subtype_printable_names[event_subtype][2]
    return event_subtype.name


def regions_to_str(regions):
    return ",".join([str(x[0]) + "-" + str(x[1]) for x in regions])


def match_subtype_to_str_with_additional_info(event, strand, read_introns, isoform_introns):
    event_subtype = event.event_type
    additional_info = ""
    if event_subtype in {MatchEventSubtype.intron_retention,
                         MatchEventSubtype.unspliced_intron_retention,
                         MatchEventSubtype.incomplete_intron_retention_left,
                         MatchEventSubtype.incomplete_intron_retention_right,
                         MatchEventSubtype.fake_micro_intron_retention}:
        if event.isoform_region != SupplementaryMatchConstants.undefined_region:
            introns = isoform_introns[event.isoform_region[0]:event.isoform_region[1]+1]
            additional_info = ":" + regions_to_str(introns)
    elif event_subtype in {MatchEventSubtype.major_exon_elongation_left, MatchEventSubtype.major_exon_elongation_right,
                           MatchEventSubtype.exon_elongation_left, MatchEventSubtype.exon_elongation_right,
                           MatchEventSubtype.alternative_tss_left, MatchEventSubtype.alternative_tss_right,
                           MatchEventSubtype.alternative_polya_site_left, MatchEventSubtype.alternative_polya_site_right,
                           MatchEventSubtype.correct_polya_site_right, MatchEventSubtype.correct_polya_site_left,
                           MatchEventSubtype.internal_polya_left, MatchEventSubtype.internal_polya_right,
                           MatchEventSubtype.terminal_site_match_left, MatchEventSubtype.terminal_site_match_right,
                           MatchEventSubtype.terminal_site_match_left_precise,
                           MatchEventSubtype.terminal_site_match_right_precise}:
        # elongation events
        additional_info = ":" + str(event.event_info)
    else:
        if event.read_region != SupplementaryMatchConstants.undefined_region and \
                event.read_region[0] >= 0 and event.read_region[1] >= 0:
            introns = read_introns[event.read_region[0]:event.read_region[1]+1]
            additional_info = ":" + regions_to_str(introns)

    return match_subtype_to_str(event, strand) + additional_info


def is_matching_assignment(isoform_assignment):
    if isoform_assignment.assignment_type == ReadAssignmentType.unique:
        return True
    elif isoform_assignment.assignment_type.is_unique():
        allowed_set = {MatchEventSubtype.none,
                       MatchEventSubtype.fsm,
                       MatchEventSubtype.exon_misalignment,
                       MatchEventSubtype.intron_shift,
                       MatchEventSubtype.terminal_site_match_left,
                       MatchEventSubtype.terminal_site_match_right,
                       MatchEventSubtype.terminal_site_match_left_precise,
                       MatchEventSubtype.terminal_site_match_right_precise,
                       MatchEventSubtype.correct_polya_site_right,
                       MatchEventSubtype.correct_polya_site_left,
                       MatchEventSubtype.exon_elongation_left,
                       MatchEventSubtype.exon_elongation_right}

        return all(m.event_type in allowed_set for m in isoform_assignment.isoform_matches[0].match_subclassifications)
    return False
