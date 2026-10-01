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

from isoquant_lib.assignment.isoform_assignment import (
    ReadAssignment, ReadAssignmentType, IsoformMatch,
    MatchClassification, MatchEvent, MatchEventSubtype,
)
from isoquant_lib.quantification import rna_velocity_counter
from isoquant_lib.quantification.long_read_counter import (
    IntronCounter, IntronRetentionCounter, INTRON_RETENTION_EVENTS,
)
from isoquant_lib.utils.string_pools import StringPoolManager
from isoquant_lib.scripts.convert_read_info import (
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
        assert MatchEventSubtype.intron_retention in INTRON_RETENTION_EVENTS

    def test_includes_unspliced_intron_retention(self):
        assert MatchEventSubtype.unspliced_intron_retention in INTRON_RETENTION_EVENTS

    def test_includes_incomplete(self):
        assert MatchEventSubtype.incomplete_intron_retention_left in INTRON_RETENTION_EVENTS
        assert MatchEventSubtype.incomplete_intron_retention_right in INTRON_RETENTION_EVENTS

    def test_excludes_fake(self):
        assert MatchEventSubtype.fake_micro_intron_retention not in INTRON_RETENTION_EVENTS

    def test_shared_with_velocity(self):
        assert rna_velocity_counter.INTRON_RETENTION_EVENTS is INTRON_RETENTION_EVENTS


class TestIntronRetentionCounter:
    def setup_method(self):
        self.tmpdir = tempfile.mkdtemp()
        self.output_prefix = os.path.join(self.tmpdir, "ir_test")
        self.counter = IntronRetentionCounter(self.output_prefix)
        self.string_pools = _make_string_pools(genes=("gene1", "gene2"), transcripts=("tx1", "tx2"))

        # Gene has two introns: (200,300) and (400,500); tx2 shares the first one
        intron_features = [(200, 300), (400, 500)]
        prop0 = _FakeFeatureProperty(0, "chr1\t200\t300\t+\tintron\tgene1")
        prop1 = _FakeFeatureProperty(1, "chr1\t400\t500\t+\tintron\tgene1")
        intron_property_map = [prop0, prop1]
        all_isoforms_introns = {"tx1": [(200, 300), (400, 500)], "tx2": [(200, 300)]}

        self.gene_info = _FakeGeneInfo(intron_features, intron_property_map, all_isoforms_introns)

    def _make_ra_with_ir(self, event_type, isoform_region, intron_gene_profile=None):
        event = MatchEvent(event_type, isoform_region=isoform_region)
        ra = _make_assignment(self.string_pools, events=[event])
        ra.gene_info = self.gene_info
        # by default the read retains every intron of the event and splices none
        if intron_gene_profile is None:
            intron_gene_profile = [0, 0]
            for idx in range(isoform_region[0], isoform_region[1] + 1):
                if idx < len(intron_gene_profile):
                    intron_gene_profile[idx] = -1
        ra.intron_gene_profile = intron_gene_profile
        return ra

    def _assert_empty(self):
        assert len(self.counter.inclusion_feature_counter) == 0
        assert len(self.counter.exclusion_feature_counter) == 0
        assert len(self.counter.encountered_group_ids) == 0

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

    def test_builds_and_reuses_feature_index_cache(self):
        # coords -> feature index map is built once and reused across reads
        ra = self._make_ra_with_ir(MatchEventSubtype.intron_retention, (0, 0))
        self.counter.add_read_info(ra)
        cached = self.gene_info._intron_feature_index
        assert cached == {(200, 300): 0, (400, 500): 1}
        ra2 = self._make_ra_with_ir(MatchEventSubtype.unspliced_intron_retention, (1, 1))
        self.counter.add_read_info(ra2)
        assert self.gene_info._intron_feature_index is cached  # not rebuilt
        assert self.counter.inclusion_feature_counter[0].get(0) == 1
        assert self.counter.inclusion_feature_counter[1].get(0) == 1

    @pytest.mark.parametrize("event_type", [MatchEventSubtype.incomplete_intron_retention_left,
                                            MatchEventSubtype.incomplete_intron_retention_right])
    def test_counts_incomplete_retention(self, event_type):
        ra = self._make_ra_with_ir(event_type, (0, 0))
        self.counter.add_read_info(ra)
        assert self.counter.inclusion_feature_counter[0].get(0) == 1

    def test_ignores_fake_micro_retention(self):
        ra = self._make_ra_with_ir(MatchEventSubtype.fake_micro_intron_retention, (0, 0))
        self.counter.add_read_info(ra)
        assert self.counter.inclusion_feature_counter[0].get(0) == 0

    def test_ignores_non_retention_events(self):
        event = MatchEvent(MatchEventSubtype.exon_skipping_novel, isoform_region=(0, 0))
        ra = _make_assignment(self.string_pools, events=[event])
        ra.gene_info = self.gene_info
        ra.intron_gene_profile = [-1, 0]
        self.counter.add_read_info(ra)
        assert self.counter.inclusion_feature_counter[0].get(0) == 0

    def test_retention_shared_by_matched_isoforms_counted_once(self):
        # read is inconsistent against tx1 and tx2 of one gene, both sharing the retained intron
        event = MatchEvent(MatchEventSubtype.intron_retention, isoform_region=(0, 0))
        m1 = IsoformMatch(MatchClassification.novel_in_catalog, self.string_pools, assigned_gene="gene1",
                          assigned_transcript="tx1", match_subclassification=[event])
        m2 = IsoformMatch(MatchClassification.novel_in_catalog, self.string_pools, assigned_gene="gene1",
                          assigned_transcript="tx2", match_subclassification=[event])
        ra = ReadAssignment("read_1", ReadAssignmentType.inconsistent_ambiguous, self.string_pools, match=[m1, m2])
        assert ra.gene_assignment_type == ReadAssignmentType.inconsistent
        ra.exon_gene_profile = [1, 1, 1]
        ra.intron_gene_profile = [-1, 0]
        ra.strand = "+"
        ra.gene_info = self.gene_info
        self.counter.add_read_info(ra)
        assert self.counter.inclusion_feature_counter[0].get(0) == 1

    def test_counts_spliced_introns_as_exclusion(self):
        # read retains intron 1 and splices intron 0
        ra = self._make_ra_with_ir(MatchEventSubtype.intron_retention, (1, 1), intron_gene_profile=[1, -1])
        self.counter.add_read_info(ra)
        assert self.counter.exclusion_feature_counter[0].get(0) == 1
        assert self.counter.inclusion_feature_counter[0].get(0) == 0
        assert self.counter.inclusion_feature_counter[1].get(0) == 1
        assert self.counter.exclusion_feature_counter[1].get(0) == 0

    def test_retained_intron_not_counted_as_spliced(self):
        # inconsistent input: profile claims the retained intron is spliced; retention wins
        ra = self._make_ra_with_ir(MatchEventSubtype.intron_retention, (0, 0), intron_gene_profile=[1, 1])
        self.counter.add_read_info(ra)
        assert self.counter.inclusion_feature_counter[0].get(0) == 1
        assert self.counter.exclusion_feature_counter[0].get(0) == 0
        assert self.counter.exclusion_feature_counter[1].get(0) == 1

    def test_spliced_only_read_counted(self):
        ra = _make_assignment(self.string_pools, assignment_type=ReadAssignmentType.unique)
        ra.gene_info = self.gene_info
        self.counter.add_read_info(ra)
        assert self.counter.exclusion_feature_counter[0].get(0) == 1
        assert self.counter.exclusion_feature_counter[1].get(0) == 1
        assert len(self.counter.inclusion_feature_counter) == 0

    def test_exclusion_respects_feature_strand(self):
        self.gene_info.intron_property_map[0].strand = "-"
        ra = _make_assignment(self.string_pools, assignment_type=ReadAssignmentType.unique)
        ra.gene_info = self.gene_info
        self.counter.add_read_info(ra)
        assert self.counter.exclusion_feature_counter[0].get(0) == 0
        assert self.counter.exclusion_feature_counter[1].get(0) == 1

    def test_exclusion_equals_splice_junction_inclusion(self):
        junction_counter = IntronCounter(os.path.join(self.tmpdir, "sj_test"))
        reads = [
            self._make_ra_with_ir(MatchEventSubtype.intron_retention, (1, 1), intron_gene_profile=[1, -1]),
            self._make_ra_with_ir(MatchEventSubtype.unspliced_intron_retention, (0, 0)),
            _make_assignment(self.string_pools, assignment_type=ReadAssignmentType.unique),
        ]
        reads[2].gene_info = self.gene_info
        for ra in reads:
            self.counter.add_read_info(ra)
            junction_counter.add_read_info(ra)
        for feature_id in (0, 1):
            assert self.counter.exclusion_feature_counter[feature_id].get(0) == \
                   junction_counter.inclusion_feature_counter[feature_id].get(0)

    def test_skips_none_assignment(self):
        self.counter.add_read_info(None)
        self._assert_empty()

    def test_skips_unassigned(self):
        # no gene at all: gene assignment stays noninformative
        match = IsoformMatch(MatchClassification.intergenic, self.string_pools)
        ra = ReadAssignment("read_1", ReadAssignmentType.noninformative, self.string_pools, match=match)
        assert ra.gene_assignment_type.is_unassigned()
        ra.exon_gene_profile = [1, 1, 1]
        ra.intron_gene_profile = [1, 1]
        ra.strand = "+"
        ra.gene_info = self.gene_info
        self.counter.add_read_info(ra)
        self._assert_empty()

    def test_skips_ambiguous_gene(self):
        event = MatchEvent(MatchEventSubtype.intron_retention, isoform_region=(0, 0))
        m1 = IsoformMatch(MatchClassification.novel_in_catalog, self.string_pools, assigned_gene="gene1",
                          assigned_transcript="tx1", match_subclassification=[event])
        m2 = IsoformMatch(MatchClassification.novel_in_catalog, self.string_pools, assigned_gene="gene2",
                          assigned_transcript="tx2", match_subclassification=[event])
        ra = ReadAssignment("read_1", ReadAssignmentType.ambiguous, self.string_pools, match=[m1, m2])
        assert ra.gene_assignment_type == ReadAssignmentType.ambiguous
        ra.exon_gene_profile = [1, 1, 1]
        ra.intron_gene_profile = [1, 1]
        ra.strand = "+"
        ra.gene_info = self.gene_info
        self.counter.add_read_info(ra)
        self._assert_empty()

    def test_skips_out_of_range_index(self):
        """Event with isoform_region beyond intron list should be skipped gracefully."""
        ra = self._make_ra_with_ir(MatchEventSubtype.intron_retention, (5, 5))
        self.counter.add_read_info(ra)
        # Should not crash, no counts added
        self._assert_empty()

    def test_multiple_reads_accumulate(self):
        ra1 = self._make_ra_with_ir(MatchEventSubtype.intron_retention, (0, 0))
        ra2 = self._make_ra_with_ir(MatchEventSubtype.unspliced_intron_retention, (0, 0))
        self.counter.add_read_info(ra1)
        self.counter.add_read_info(ra2)
        assert self.counter.inclusion_feature_counter[0].get(0) == 2

    def test_dump_writes_include_and_exclude(self):
        ra = self._make_ra_with_ir(MatchEventSubtype.intron_retention, (1, 1), intron_gene_profile=[1, -1])
        self.counter.add_read_info(ra)
        self.counter.dump()
        with open(self.output_prefix + "_counts.tsv") as f:
            rows = [line.rstrip("\n").split("\t") for line in f]
        assert rows[0][-2:] == ["include_counts", "exclude_counts"]
        counts = {(row[1], row[2]): (int(row[-2]), int(row[-1])) for row in rows[1:]}
        assert counts == {("200", "300"): (0, 1), ("400", "500"): (1, 0)}


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
            # Exactly one data line -- the read_info column header must NOT be
            # emitted as a spurious record (only lines[0] may start with read_id)
            assert all(not line.startswith("read_id") for line in lines[1:])
            data_lines = lines[1:]
            assert len(data_lines) == 1
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

    def test_comment_lines_preserved_and_no_duplicate_header(self):
        """Leading '#' comments pass through; the header is not duplicated as data."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("# Command line: isoquant.py ...\n")
            f.write("# IsoQuant version: 3.x\n")
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

            comment_lines = [line for line in lines if line.startswith("#")]
            assert len(comment_lines) == 2
            # Comments precede the single header line
            header_idx = [i for i, line in enumerate(lines) if line.startswith("read_id")]
            assert header_idx == [2]
            # Exactly one converted data row
            data_lines = [line for line in lines if not line.startswith("#") and not line.startswith("read_id")]
            assert len(data_lines) == 1
            assert data_lines[0].split("\t")[0] == "r1"
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
