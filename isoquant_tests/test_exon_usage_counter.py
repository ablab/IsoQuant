############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Tests for ExonUsageCounter: per-exon full / half / skip / alt quantification."""

import os
import random
import tempfile
from bisect import bisect_left, bisect_right
from functools import partial

import pytest

from isoquant_lib.assignment.long_read_profiles import NonOverlappingFeaturesProfileConstructor
from isoquant_lib.common import overlaps_at_least_when_overlap
from isoquant_lib.gene_info import FeatureInfo, GeneInfo
from isoquant_lib.quantification.convert_grouped_counts import _load_profile_linear
from isoquant_lib.quantification.long_read_counter import (
    ExonUsageCounter, EXON_FULL, EXON_LEFT, EXON_RIGHT, EXON_SKIP, EXON_ALT,
)
from isoquant_lib.assignment.isoform_assignment import (
    ReadAssignment, ReadAssignmentType, IsoformMatch, MatchClassification,
)
from isoquant_lib.utils.string_pools import StringPoolManager


# ── helpers ──────────────────────────────────────────────────────────────

class _Features:
    def __init__(self, features):
        self.features = features


def _make_gene_info(transcripts: dict, strand: str = "+", gene_id: str = "gene1") -> GeneInfo:
    """GeneInfo carrying just what ExonUsageCounter touches, index built by the real GeneInfo code."""
    gene_info = GeneInfo.__new__(GeneInfo)
    gene_info.chr_id = "chr1"
    gene_info.all_isoforms_exons = transcripts
    exons = sorted(set(e for exons in transcripts.values() for e in exons))
    gene_info.exon_profiles = _Features(exons)
    gene_info.split_exon_profiles = _Features(GeneInfo.split_exons(exons))
    gene_info.exon_property_map = [FeatureInfo("chr1", e[0], e[1], strand, "I", [gene_id]) for e in exons]
    return gene_info


# exon2 (300,400) is a cassette exon, (300,450) its alternative 3' variant,
# (300,600) retains the intron between exon2 and exon3, chaining them all into one overlap cluster
MAIN_GENE = {
    "T1": [(100, 200), (300, 400), (500, 600), (700, 800)],
    "T2": [(100, 200), (500, 600)],
    "T3": [(100, 200), (300, 450), (500, 600)],
    "T4": [(100, 200), (300, 600)],
}


def _make_string_pools():
    sp = StringPoolManager()
    sp.gene_pool.add("gene1")
    sp.gene_pool.add("gene2")
    sp.transcript_pool.add("tx1")
    return sp


class TestExonUsageCounter:
    def setup_method(self):
        self.tmpdir = tempfile.mkdtemp()
        self.counter = ExonUsageCounter(os.path.join(self.tmpdir, "exon_usage"), delta=6)
        self.string_pools = _make_string_pools()
        self.gene_info = _make_gene_info(MAIN_GENE)

    def _add(self, blocks, strand="+", gene="gene1", gene_info=None):
        match = IsoformMatch(MatchClassification.full_splice_match, self.string_pools,
                             assigned_gene=gene, assigned_transcript="tx1")
        ra = ReadAssignment("read_1", ReadAssignmentType.unique, self.string_pools, match=match)
        ra.exon_gene_profile = []
        ra.intron_gene_profile = []
        ra.exons = blocks
        ra.corrected_exons = blocks
        ra.strand = strand
        ra.gene_info = gene_info if gene_info is not None else self.gene_info
        self.counter.add_read_info(ra)

    def _state_counts(self, exon, gene="gene1"):
        key = (("chr1", exon[0], exon[1], "+"), gene)
        return self.counter.exon_counts.get(key, {}).get(0, [0] * 5)

    def _states(self):
        result = {}
        for (feature_id, _), group_counts in self.counter.exon_counts.items():
            counts = group_counts[0]
            nonzero = [s for s in range(5) if counts[s]]
            assert len(nonzero) == 1
            result[(feature_id[1], feature_id[2])] = nonzero[0]
        return result

    def test_skip_within_overlap_chain(self):
        # read skips exon2: the cassette exon and its 3' variant are skipped even though
        # the intron-retaining (300,600) is chained with them and overlapped by the read
        self._add([(100, 200), (500, 600)])
        assert self._states() == {
            (100, 200): EXON_FULL,
            (300, 400): EXON_SKIP,
            (300, 450): EXON_SKIP,
            (300, 600): EXON_ALT,
            # last exon of T2 / T3, but internal in T1: the read end does not confirm its right splice site
            (500, 600): EXON_LEFT,
        }

    def test_internal_inclusion_and_alt_variants(self):
        self._add([(100, 200), (300, 400), (500, 600), (700, 800)])
        assert self._states() == {
            (100, 200): EXON_FULL,
            (300, 400): EXON_FULL,
            (300, 450): EXON_ALT,
            (300, 600): EXON_ALT,
            (500, 600): EXON_FULL,
            (700, 800): EXON_FULL,
        }

    def test_alternative_splice_site(self):
        self._add([(100, 200), (302, 448), (500, 600)])
        states = self._states()
        assert states[(300, 450)] == EXON_FULL
        assert states[(300, 400)] == EXON_ALT

    def test_exons_outside_read_are_not_counted(self):
        # read stops (e.g. at its polyA tail) inside exon2: downstream exons are not skipped
        self._add([(100, 200), (300, 380)])
        states = self._states()
        assert (500, 600) not in states
        assert (700, 800) not in states
        assert states[(300, 400)] == EXON_LEFT

    def test_half_inclusion_internal_exon(self):
        # read starts inside the internal cassette exon: only its right splice site is confirmed
        self._add([(350, 400), (500, 600)])
        states = self._states()
        assert states[(300, 400)] == EXON_RIGHT
        assert states[(300, 450)] == EXON_ALT
        assert states[(500, 600)] == EXON_LEFT

    def test_terminal_exons_need_only_internal_splice_site(self):
        self._add([(150, 200), (300, 400), (500, 600), (700, 750)])
        states = self._states()
        assert states[(100, 200)] == EXON_FULL
        assert states[(700, 800)] == EXON_FULL

    def test_terminal_and_internal_exon_is_treated_as_internal(self):
        # (500, 600) is the last exon of T2 / T3 and internal in T1
        self._add([(150, 200), (300, 400), (500, 550)])
        assert self._states()[(500, 600)] == EXON_LEFT

    def test_first_and_last_exon_is_treated_as_internal(self):
        gene_info = _make_gene_info({"T1": [(100, 200), (300, 400)],
                                     "T2": [(10, 50), (100, 200)]})
        self._add([(150, 200), (300, 400)], gene_info=gene_info)
        self._add([(10, 50), (100, 150)], gene_info=gene_info)
        assert self._state_counts((100, 200))[EXON_RIGHT] == 1
        assert self._state_counts((100, 200))[EXON_LEFT] == 1
        assert self._state_counts((100, 200))[EXON_FULL] == 0

    def test_terminal_block_extending_past_exon_is_alt(self):
        self._add([(50, 200), (300, 400)])
        assert self._states()[(100, 200)] == EXON_ALT

    def test_shared_splice_site_counts_for_all_exons(self):
        gene_info = _make_gene_info({"T1": [(100, 200), (300, 400), (500, 600)],
                                     "T2": [(100, 200), (250, 400), (500, 600)]})
        self._add([(350, 400), (500, 600)], gene_info=gene_info)
        assert self._state_counts((300, 400))[EXON_RIGHT] == 1
        assert self._state_counts((250, 400))[EXON_RIGHT] == 1

    def test_counts_accumulate_with_delta(self):
        # the second read ends 2bp past the last exon, within delta
        self._add([(100, 200), (500, 600)])
        self._add([(100, 200), (500, 600), (700, 802)])
        assert self._state_counts((300, 400))[EXON_SKIP] == 2
        assert self._state_counts((700, 800))[EXON_FULL] == 1

    def test_read_with_intron_inside_exon_is_alt(self):
        # both blocks lie within the intron-retaining (300, 600), the read uses shorter exons instead
        self._add([(300, 400), (500, 550)])
        self._add([(350, 400), (500, 600)])
        assert self._state_counts((300, 600))[EXON_ALT] == 2

    def test_uninformative_reads(self):
        self._add([(100, 600)])
        self._add([(100, 200), (500, 600)], strand="-")
        self._add([(100, 200), (500, 600)], gene="gene2")
        assert not self.counter.exon_counts

    def test_dump_format_and_matrix_loading(self):
        self._add([(100, 200), (500, 600)])
        self._add([(100, 200), (300, 400), (500, 600)])
        self._add([(350, 400), (500, 600)])
        self.counter.dump()
        with open(self.counter.output_counts_file_name) as f:
            lines = f.read().strip().split("\n")
        assert lines[0] == ("chr\tstart\tend\tstrand\tflags\tgene_ids\tgroup_id"
                            "\tinclude_counts\texclude_counts\tn_full\tn_left\tn_right\tn_alt")
        rows = {(int(r[1]), int(r[2])): r for r in (line.split("\t") for line in lines[1:])}
        # cassette exon: 1 full + 1 right half included, 1 skipped
        assert rows[(300, 400)][5:] == ["gene1", "NA", "2", "1", "1", "0", "1", "0"]
        assert rows[(500, 600)][7:9] == ["3", "0"]

        df = _load_profile_linear(self.counter.output_counts_file_name)
        assert df is not None
        cassette = df[(df["start"] == 300) & (df["end"] == 400)].iloc[0]
        assert cassette["include_counts"] == 2 and cassette["exclude_counts"] == 1


def _random_blocks(rng: random.Random, max_pos: int) -> list:
    # sorted non-adjacent blocks (at least 1bp introns) within [1, max_pos]
    n = rng.randint(1, 5)
    borders = sorted(rng.sample(range(1, max_pos), 2 * n))
    blocks = [(borders[2 * i], borders[2 * i + 1]) for i in range(n)]
    return [b for i, b in enumerate(blocks) if i == 0 or b[0] > blocks[i - 1][1] + 1]


class TestSlicedSplitProfile:
    """The counter builds the split-exon profile only over segments within the read span;
    it must equal the corresponding part of the profile built over all segments (as the assigner does)."""

    def test_sliced_profile_equals_full_profile(self):
        rng = random.Random(42)
        comparator = partial(overlaps_at_least_when_overlap, delta=5)
        checked = 0
        for _ in range(3000):
            exons = sorted({(s, s + rng.randint(0, 60)) for s in (rng.randint(1, 300) for _ in range(rng.randint(1, 8)))})
            segments = GeneInfo.split_exons(exons)
            blocks = _random_blocks(rng, 400)
            first = bisect_left([seg[1] for seg in segments], blocks[0][0])
            last = bisect_right([seg[0] for seg in segments], blocks[-1][1])
            if first >= last:
                continue
            full = NonOverlappingFeaturesProfileConstructor(
                segments, comparator=comparator).construct_profile(blocks).gene_profile
            sliced = NonOverlappingFeaturesProfileConstructor(
                segments[first:last], comparator=comparator).construct_profile(blocks).gene_profile
            assert sliced == full[first:last], (exons, blocks)
            # segments outside the read span are uninformative
            assert all(v == 0 for v in full[:first] + full[last:]), (exons, blocks)
            checked += 1
        assert checked > 1000


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
