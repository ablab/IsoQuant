############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Tests for JointExonCounter: region-based exon quantification."""

import os
import tempfile

import pytest

from isoquant_lib.gene_info import ExonRegion, GeneInfo
from isoquant_lib.long_read_counter import JointExonCounter
from isoquant_lib.isoform_assignment import (
    ReadAssignment, ReadAssignmentType, IsoformMatch, MatchClassification,
)
from isoquant_lib.string_pools import StringPoolManager


# ── helpers ──────────────────────────────────────────────────────────────

def _make_string_pools():
    sp = StringPoolManager()
    sp.gene_pool.add("gene1")
    sp.transcript_pool.add("tx1")
    return sp


class _FakeFeatureInfo:
    def __init__(self, fid, strand, gene_ids=("gene1",)):
        self.id = fid
        # production FeatureInfo.strand is a concatenated string (e.g., "+", "+-")
        self.strand = strand
        self.gene_ids = list(gene_ids)


class _FakeExonProfiles:
    def __init__(self, features):
        self.features = features


class _FakeGeneInfo:
    """Minimal GeneInfo stand-in carrying just what JointExonCounter touches."""
    def __init__(self, exon_features, exon_property_map, exon_overlap_regions):
        self.exon_profiles = _FakeExonProfiles(exon_features)
        self.exon_property_map = exon_property_map
        self.exon_overlap_regions = exon_overlap_regions


def _make_assignment(string_pools, exon_profile, strand="+", gene="gene1"):
    match = IsoformMatch(
        MatchClassification.full_splice_match,
        string_pools,
        assigned_gene=gene,
        assigned_transcript="tx1",
    )
    ra = ReadAssignment("read_1", ReadAssignmentType.unique, string_pools, match=match)
    ra.exon_gene_profile = exon_profile
    ra.intron_gene_profile = [1] * max(0, len(exon_profile) - 1)
    ra.strand = strand
    return ra


def _make_gene_info():
    """Three overlapping exons + one singleton on the same strand.

    overlap region [100, 350]: exons (100,200), (150,250), (200,350)
    singleton region [500, 600]: exon (500,600)
    """
    features = [(100, 200), (150, 250), (200, 350), (500, 600)]
    props = [
        _FakeFeatureInfo(10, '+'),
        _FakeFeatureInfo(11, '+'),
        _FakeFeatureInfo(12, '+'),
        _FakeFeatureInfo(13, '+'),
    ]
    overlap = ExonRegion("chr1", 100, 350, '+')
    overlap.member_exon_indices = [0, 1, 2]
    overlap.gene_ids = {"gene1"}
    single = ExonRegion("chr1", 500, 600, '+')
    single.member_exon_indices = [3]
    single.gene_ids = {"gene1"}
    return _FakeGeneInfo(features, props, [overlap, single])


# ── tests ────────────────────────────────────────────────────────────────

class TestExonOverlapRegions:
    def test_grouping_via_geneinfo_method(self):
        """GeneInfo.build_exon_overlap_regions groups overlapping exons by strand."""
        gi = GeneInfo.__new__(GeneInfo)
        gi.chr_id = "chr1"
        features = [(100, 200), (150, 250), (200, 350), (500, 600)]
        property_map = [_FakeFeatureInfo(0, '+'), _FakeFeatureInfo(1, '+'),
                        _FakeFeatureInfo(2, '+'), _FakeFeatureInfo(3, '+')]
        regions, region_map = gi.build_exon_overlap_regions(features, property_map)
        assert len(regions) == 2
        assert regions[0].start == 100 and regions[0].end == 350
        assert regions[0].member_exon_indices == [0, 1, 2]
        assert regions[1].start == 500 and regions[1].end == 600
        assert regions[1].member_exon_indices == [3]
        assert region_map == [0, 0, 0, 1]

    def test_strand_separation(self):
        """Exons on different strands never share a region even if coords overlap."""
        gi = GeneInfo.__new__(GeneInfo)
        gi.chr_id = "chr1"
        features = [(100, 200), (150, 250)]
        property_map = [_FakeFeatureInfo(0, '+'), _FakeFeatureInfo(1, '-')]
        regions, region_map = gi.build_exon_overlap_regions(features, property_map)
        assert len(regions) == 2
        assert regions[0].strand != regions[1].strand


class TestJointExonCounter:
    def setup_method(self):
        self.tmpdir = tempfile.mkdtemp()
        self.output_prefix = os.path.join(self.tmpdir, "joint_test")
        self.counter = JointExonCounter(self.output_prefix)
        self.string_pools = _make_string_pools()
        self.gene_info = _make_gene_info()

    def _add(self, profile, strand="+"):
        ra = _make_assignment(self.string_pools, profile, strand=strand)
        ra.gene_info = self.gene_info
        self.counter.add_read_info(ra)

    def test_inclusion_in_overlap_region(self):
        # exon 0 included, others excluded; singleton has -1 → region exclusion
        self._add([1, -1, -1, -1])
        assert self.counter.inclusion_counter[(10, "gene1")].get(0) == 1
        assert self.counter.inclusion_counter[(11, "gene1")].get(0) == 0
        assert self.counter.inclusion_counter[(12, "gene1")].get(0) == 0
        # the singleton region's exclusion fired
        singleton_region_id = self.gene_info.exon_overlap_regions[1].id
        assert self.counter.exclusion_counter[(singleton_region_id, "gene1")].get(0) == 1

    def test_all_exclusion_increments_region(self):
        # all three overlap members excluded → +1 region exclusion (not per-exon)
        self._add([-1, -1, -1, 1])
        overlap_region_id = self.gene_info.exon_overlap_regions[0].id
        assert self.counter.exclusion_counter[(overlap_region_id, "gene1")].get(0) == 1
        # singleton inclusion fires
        assert self.counter.inclusion_counter[(13, "gene1")].get(0) == 1

    def test_zero_state_skips(self):
        # any 0 in region's members → skip that region
        self._add([0, -1, -1, 1])
        overlap_region_id = self.gene_info.exon_overlap_regions[0].id
        assert self.counter.exclusion_counter[(overlap_region_id, "gene1")].get(0) == 0
        for fid in (10, 11, 12):
            assert self.counter.inclusion_counter[(fid, "gene1")].get(0) == 0
        # singleton still increments
        assert self.counter.inclusion_counter[(13, "gene1")].get(0) == 1

    def test_strand_filter(self):
        # opposite-strand read → no counts in any region
        self._add([1, -1, -1, 1], strand="-")
        for fid in (10, 11, 12, 13):
            assert self.counter.inclusion_counter[(fid, "gene1")].get(0) == 0
        for region in self.gene_info.exon_overlap_regions:
            assert self.counter.exclusion_counter[(region.id, "gene1")].get(0) == 0

    def test_multiple_inclusions_in_same_region(self):
        # two annotated variants both match (e.g., alternative acceptors share span)
        self._add([1, 1, -1, -1])
        assert self.counter.inclusion_counter[(10, "gene1")].get(0) == 1
        assert self.counter.inclusion_counter[(11, "gene1")].get(0) == 1
        assert self.counter.inclusion_counter[(12, "gene1")].get(0) == 0
        overlap_region_id = self.gene_info.exon_overlap_regions[0].id
        # mixed inclusion + exclusion means region exclusion does NOT fire
        assert self.counter.exclusion_counter[(overlap_region_id, "gene1")].get(0) == 0

    def test_dump_format(self):
        self._add([1, -1, -1, -1])
        self._add([-1, -1, -1, 1])
        self.counter.dump()
        with open(self.counter.output_counts_file_name) as f:
            lines = f.read().strip().split("\n")
        # header + at least 3 data rows: exon 10 inclusion, overlap exclusion, exon 13 inclusion
        assert lines[0] == ("chr\tregion_start\tregion_end\tstrand\texon_start\texon_end"
                            "\tgene_id\tfeature_kind\tgroup_id\tcount")
        data = [l.split("\t") for l in lines[1:]]
        kinds = {(row[0], int(row[1]), int(row[2]), row[7]) for row in data}
        assert ("chr1", 100, 350, "inclusion") in kinds
        assert ("chr1", 100, 350, "exclusion") in kinds
        # singleton inclusion appears with coords 500..600
        assert any(row[0] == "chr1" and row[4] == "500" and row[5] == "600"
                   and row[7] == "inclusion" for row in data)
        # gene_id column populated (non-empty) for every data row
        assert all(row[6] == "gene1" for row in data)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
