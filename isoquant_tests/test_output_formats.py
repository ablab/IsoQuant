############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Tests for output format improvements:
- ReadInfoPrinter (read_info.tsv unified format)
- IntronRetentionCounter (event-based intron retention counting)
- convert_read_info (read_info → legacy format conversion)
"""

import os
import tempfile

import pytest

from isoquant_lib.isoform_assignment import (
    ReadAssignment, ReadAssignmentType, IsoformMatch,
    MatchClassification, MatchEvent, MatchEventSubtype,
)
from isoquant_lib.long_read_counter import IntronRetentionCounter, INTRON_RETENTION_EVENT_TYPES
from isoquant_lib.string_pools import StringPoolManager
from isoquant_lib.gene_info import FeatureInfo, FeatureProfiles
from isoquant_lib.convert_read_info import (
    _parse_exons, _exons_to_range_str, _assignment_type_to_read_type,
    convert_to_read_assignments, convert_to_allinfo,
)


# ── helpers ──────────────────────────────────────────────────────────────

def _make_string_pools(genes=("gene1",), transcripts=("tx1",)):
    sp = StringPoolManager()
    for g in genes:
        sp.gene_pool.add(g)
    for t in transcripts:
        sp.transcript_pool.add(t)
    return sp


class _FakeGeneInfo:
    """Minimal GeneInfo stand-in for counter tests."""
    def __init__(self, intron_features, intron_property_map, all_isoforms_introns):
        self.intron_profiles = _FakeProfiles(intron_features)
        self.intron_property_map = intron_property_map
        self.all_isoforms_introns = all_isoforms_introns
        self.reference_region = None


class _FakeProfiles:
    def __init__(self, features):
        self.features = features


class _FakeFeatureProperty:
    """Minimal FeatureInfo stand-in."""
    def __init__(self, fid, name_str):
        self.id = fid
        self._name = name_str
        self.strand = {'+', '-', '.'}

    def to_str(self):
        return self._name


def _make_assignment(string_pools, gene="gene1", transcript="tx1",
                     assignment_type=ReadAssignmentType.inconsistent,
                     events=None):
    """Create a ReadAssignment with isoform matches and events."""
    match = IsoformMatch(
        MatchClassification.full_splice_match,
        string_pools,
        assigned_gene=gene,
        assigned_transcript=transcript,
        match_subclassification=events or [],
    )
    ra = ReadAssignment("read_1", assignment_type, string_pools, match=match)
    ra.exons = [(100, 200), (300, 400), (500, 600)]
    ra.corrected_exons = [(100, 200), (300, 400), (500, 600)]
    ra.exon_gene_profile = [1, 1, 1]
    ra.intron_gene_profile = [1, 1]
    ra.strand = "+"
    ra.chr_id = "chr1"
    return ra


# ── IntronRetentionCounter tests ─────────────────────────────────────────

class TestIntronRetentionEventTypes:
    def test_includes_intron_retention(self):
        assert MatchEventSubtype.intron_retention in INTRON_RETENTION_EVENT_TYPES

    def test_includes_unspliced_intron_retention(self):
        assert MatchEventSubtype.unspliced_intron_retention in INTRON_RETENTION_EVENT_TYPES

    def test_excludes_incomplete(self):
        assert MatchEventSubtype.incomplete_intron_retention_left not in INTRON_RETENTION_EVENT_TYPES
        assert MatchEventSubtype.incomplete_intron_retention_right not in INTRON_RETENTION_EVENT_TYPES

    def test_excludes_fake(self):
        assert MatchEventSubtype.fake_micro_intron_retention not in INTRON_RETENTION_EVENT_TYPES


class TestIntronRetentionCounter:
    def setup_method(self):
        self.tmpdir = tempfile.mkdtemp()
        self.output_prefix = os.path.join(self.tmpdir, "ir_test")
        self.counter = IntronRetentionCounter(self.output_prefix)
        self.string_pools = _make_string_pools()

        # Gene has two introns: (200,300) and (400,500)
        intron_features = [(200, 300), (400, 500)]
        prop0 = _FakeFeatureProperty(0, "chr1\t200\t300\t+\tintron\tgene1")
        prop1 = _FakeFeatureProperty(1, "chr1\t400\t500\t+\tintron\tgene1")
        intron_property_map = [prop0, prop1]
        all_isoforms_introns = {"tx1": [(200, 300), (400, 500)]}

        self.gene_info = _FakeGeneInfo(intron_features, intron_property_map, all_isoforms_introns)

    def _make_ra_with_ir(self, event_type, isoform_region):
        event = MatchEvent(event_type, isoform_region=isoform_region)
        ra = _make_assignment(self.string_pools, events=[event])
        ra.gene_info = self.gene_info
        return ra

    def test_counts_intron_retention(self):
        ra = self._make_ra_with_ir(MatchEventSubtype.intron_retention, (0, 0))
        self.counter.add_read_info(ra)
        assert self.counter.inclusion_feature_counter[0].get(0) == 1

    def test_counts_unspliced_intron_retention(self):
        ra = self._make_ra_with_ir(MatchEventSubtype.unspliced_intron_retention, (1, 1))
        self.counter.add_read_info(ra)
        assert self.counter.inclusion_feature_counter[1].get(0) == 1

    def test_counts_multi_intron_range(self):
        ra = self._make_ra_with_ir(MatchEventSubtype.intron_retention, (0, 1))
        self.counter.add_read_info(ra)
        assert self.counter.inclusion_feature_counter[0].get(0) == 1
        assert self.counter.inclusion_feature_counter[1].get(0) == 1

    def test_ignores_incomplete_retention(self):
        ra = self._make_ra_with_ir(MatchEventSubtype.incomplete_intron_retention_left, (0, 0))
        self.counter.add_read_info(ra)
        assert self.counter.inclusion_feature_counter[0].get(0) == 0

    def test_ignores_fake_micro_retention(self):
        ra = self._make_ra_with_ir(MatchEventSubtype.fake_micro_intron_retention, (0, 0))
        self.counter.add_read_info(ra)
        assert self.counter.inclusion_feature_counter[0].get(0) == 0

    def test_ignores_non_retention_events(self):
        event = MatchEvent(MatchEventSubtype.exon_skipping_novel, isoform_region=(0, 0))
        ra = _make_assignment(self.string_pools, events=[event])
        ra.gene_info = self.gene_info
        self.counter.add_read_info(ra)
        assert self.counter.inclusion_feature_counter[0].get(0) == 0

    def test_skips_none_assignment(self):
        self.counter.add_read_info(None)
        assert len(self.counter.inclusion_feature_counter) == 0

    def test_skips_unassigned(self):
        ra = _make_assignment(self.string_pools, assignment_type=ReadAssignmentType.noninformative)
        ra.gene_info = self.gene_info
        self.counter.add_read_info(ra)
        assert len(self.counter.inclusion_feature_counter) == 0

    def test_skips_ambiguous_gene(self):
        ra = _make_assignment(self.string_pools, assignment_type=ReadAssignmentType.ambiguous)
        ra.gene_info = self.gene_info
        self.counter.add_read_info(ra)
        assert len(self.counter.inclusion_feature_counter) == 0

    def test_skips_out_of_range_index(self):
        """Event with isoform_region beyond intron list should be skipped gracefully."""
        ra = self._make_ra_with_ir(MatchEventSubtype.intron_retention, (5, 5))
        self.counter.add_read_info(ra)
        # Should not crash, no counts added
        assert len(self.counter.inclusion_feature_counter) == 0

    def test_multiple_reads_accumulate(self):
        ra1 = self._make_ra_with_ir(MatchEventSubtype.intron_retention, (0, 0))
        ra2 = self._make_ra_with_ir(MatchEventSubtype.unspliced_intron_retention, (0, 0))
        self.counter.add_read_info(ra1)
        self.counter.add_read_info(ra2)
        assert self.counter.inclusion_feature_counter[0].get(0) == 2

    def test_dump_creates_file(self):
        ra = self._make_ra_with_ir(MatchEventSubtype.intron_retention, (0, 0))
        self.counter.add_read_info(ra)
        self.counter.dump()
        assert os.path.exists(self.output_prefix + "_counts.tsv")


# ── convert_read_info tests ──────────────────────────────────────────────

class TestParseExons:
    def test_parse_normal(self):
        assert _parse_exons("100-200,300-400") == [(100, 200), (300, 400)]

    def test_parse_single(self):
        assert _parse_exons("100-200") == [(100, 200)]

    def test_parse_dot(self):
        assert _parse_exons(".") == []


class TestExonsToRangeStr:
    def test_normal(self):
        assert _exons_to_range_str([(100, 200), (300, 400)]) == "100-200,300-400"

    def test_empty(self):
        assert _exons_to_range_str([]) == "."


class TestAssignmentTypeToReadType:
    def test_unique(self):
        assert _assignment_type_to_read_type("unique") == "known"

    def test_unique_minor(self):
        assert _assignment_type_to_read_type("unique_minor_difference") == "known"

    def test_ambiguous(self):
        assert _assignment_type_to_read_type("ambiguous") == "known_ambiguous"

    def test_inconsistent(self):
        assert _assignment_type_to_read_type("inconsistent") == "novel"

    def test_inconsistent_non_intronic(self):
        assert _assignment_type_to_read_type("inconsistent_non_intronic") == "novel"

    def test_noninformative(self):
        assert _assignment_type_to_read_type("noninformative") == "none"

    def test_intergenic(self):
        assert _assignment_type_to_read_type("intergenic") == "none"


class TestConvertToReadAssignments:
    def test_round_trip(self):
        """Write a minimal read_info file and convert to read_assignments."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("read_id\tchr\tstrand\tgene_id\tgene_assignment_type\tisoform_id\t"
                    "isoform_assignment_type\tassignment_events\tclassification\texons\t"
                    "polyA\tCAGE\tcanonical\tbarcode\tumi\tcell_type\tgroups\tadditional\n")
            f.write("r1\tchr1\t+\tGENE1\tunique\tTX1\tunique\texon_match(0,0)\t"
                    "full_splice_match\t100-200,300-400\tTrue\t.\tTrue\t"
                    "ACGT\tUMI1\t.\t.\t*\n")
            input_path = f.name

        output_path = input_path + ".ra.tsv"
        try:
            convert_to_read_assignments(input_path, output_path)
            with open(output_path) as out:
                lines = out.readlines()

            # Header line
            assert lines[0].startswith("read_id\tchr\tstrand")
            # Find actual data line (skip header-like lines)
            data_lines = [l for l in lines[1:] if not l.startswith("read_id")]
            assert len(data_lines) >= 1
            cols = data_lines[0].strip().split("\t")
            assert cols[0] == "r1"
            assert cols[1] == "chr1"
            assert cols[3] == "TX1"  # isoform_id
            assert cols[4] == "GENE1"  # gene_id
            assert cols[5] == "unique"  # assignment_type
        finally:
            os.unlink(input_path)
            if os.path.exists(output_path):
                os.unlink(output_path)


class TestConvertToAllinfo:
    def test_basic_conversion(self):
        """Write a minimal read_info file and convert to allinfo."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("read_id\tchr\tstrand\tgene_id\tgene_assignment_type\tisoform_id\t"
                    "isoform_assignment_type\tassignment_events\tclassification\texons\t"
                    "polyA\tCAGE\tcanonical\tbarcode\tumi\tcell_type\tgroups\tadditional\n")
            f.write("r1\tchr1\t+\tGENE1\tunique\tTX1\tunique\texon_match(0,0)\t"
                    "full_splice_match\t100-200,300-400\tTrue\t.\tTrue\t"
                    "ACGT\tUMI1\tTypeA\tgrp1\t*\n")
            input_path = f.name

        output_path = input_path + ".allinfo.tsv"
        try:
            convert_to_allinfo(input_path, output_path)
            with open(output_path) as out:
                lines = out.readlines()

            assert len(lines) == 1
            cols = lines[0].strip().split("\t")
            assert cols[0] == "r1"       # read_id
            assert cols[1] == "GENE1"    # gene_id
            assert cols[2] == "TypeA"    # cell_type
            assert cols[3] == "ACGT"     # barcode
            assert cols[4] == "UMI1"     # umi
            assert cols[11] == "TX1"     # isoform_id
        finally:
            os.unlink(input_path)
            if os.path.exists(output_path):
                os.unlink(output_path)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
