############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Tests for ExonSpliceSiteCounter: region-based exon splice-site counts."""

import importlib.util
import os
import tempfile
from types import SimpleNamespace

import pytest

from isoquant_lib.common import junctions_from_blocks
from isoquant_lib.gene_info import ExonSpliceSiteRegion, GeneInfo
from isoquant_lib.long_read_counter import ExonSpliceSiteCounter
from isoquant_lib.isoform_assignment import MatchEventSubtype
from isoquant_lib.string_pools import StringPoolManager


# ── helpers ──────────────────────────────────────────────────────────────

def _region(start, end, candidates, left_unique, right_unique,
            chr_id="chr1", strand="+", gene="gene1"):
    r = ExonSpliceSiteRegion(chr_id, start, end, strand, gene)
    r.candidates = list(candidates)
    r.left_unique = list(left_unique)
    r.right_unique = list(right_unique)
    return r


def _classify(counter, region, blocks, strand="+", tss=False, polya=False,
              group_id=0, read_id="r"):
    introns = junctions_from_blocks(blocks)
    counter._classify_region(region, blocks, introns, strand, tss, polya, group_id, read_id)


def _rec(counter, region):
    key = (region.chr_id, region.start, region.end, region.strand, region.gene_id)
    return counter.regions.get(key)


def _cnt(bucket):
    """Total count across all groups in a bucket (group_id -> [count, [read_ids]])."""
    return sum(entry[0] for entry in bucket.values())


class _FakeGA:
    """Stand-in for gene_assignment_type: assigned, not ambiguous."""
    def is_unassigned(self):
        return False

    def is_ambiguous(self):
        return False


def _make_ra(gene_info, blocks, strand="+", polyA=False, gene="gene1",
             events=(), read_id="r1", read_group_ids=None):
    """Lightweight ReadAssignment stub carrying only what the counter touches."""
    match = SimpleNamespace(
        assigned_gene=gene,
        match_subclassifications=[SimpleNamespace(event_type=e) for e in events],
    )
    return SimpleNamespace(
        exon_gene_profile=[1],
        intron_gene_profile=[],
        gene_info=gene_info,
        gene_assignment_type=_FakeGA(),
        isoform_matches=[match],
        corrected_exons=blocks,
        exons=blocks,
        strand=strand,
        polyA_found=polyA,
        read_id=read_id,
        read_group_ids=read_group_ids or [],
    )


def _counter(**kwargs):
    tmpdir = tempfile.mkdtemp()
    prefix = os.path.join(tmpdir, "ess_test")
    return ExonSpliceSiteCounter(prefix, **kwargs)


# ── region construction ────────────────────────────────────────────────────

class TestExonSpliceSiteRegions:
    def _gene_info(self):
        gi = GeneInfo.__new__(GeneInfo)
        gi.chr_id = "chr1"
        # tx1/tx2 share donor; alternative acceptors 100 vs 150 into exon ending at 200
        # tx3 is an intron-retention isoform: single exon spanning the (200,300) intron
        gi.all_isoforms_exons = {
            "tx1": [(100, 200), (300, 400)],
            "tx2": [(150, 200), (300, 400)],
            "tx3": [(100, 400)],
        }
        gi.all_isoforms_introns = {
            t: junctions_from_blocks(ex) for t, ex in gi.all_isoforms_exons.items()
        }
        gi.gene_id_map = {"tx1": "gene1", "tx2": "gene1", "tx3": "gene1"}
        gi.gene_strands = {"gene1": "+"}
        return gi

    def test_grouping_and_ir_drop(self):
        regions = self._gene_info().build_exon_splice_site_regions()
        assert len(regions) == 2
        a, b = sorted(regions, key=lambda r: r.start)
        assert (a.start, a.end) == (100, 200)
        assert sorted(a.candidates) == [(100, 200), (150, 200)]
        assert (100, 400) not in a.candidates
        assert (b.start, b.end) == (300, 400)
        assert b.candidates == [(300, 400)]

    def test_edge_uniqueness(self):
        regions = self._gene_info().build_exon_splice_site_regions()
        a = min(regions, key=lambda r: r.start)
        order = {c: i for i, c in enumerate(a.candidates)}
        assert a.left_unique[order[(100, 200)]] is True
        assert a.left_unique[order[(150, 200)]] is True
        assert a.right_unique[order[(100, 200)]] is False
        assert a.right_unique[order[(150, 200)]] is False

    def test_empty_when_no_exons(self):
        gi = GeneInfo.__new__(GeneInfo)
        gi.chr_id = "chr1"
        gi.all_isoforms_exons = {}
        gi.all_isoforms_introns = {}
        gi.gene_id_map = {}
        gi.gene_strands = {}
        assert gi.build_exon_splice_site_regions() == []


# ── per-read classification (via _classify_region) ─────────────────────────

class TestClassification:
    def _shared_right_region(self):
        return _region(100, 200, [(100, 200), (150, 200)], [True, True], [False, False])

    def _shared_left_region(self):
        return _region(100, 250, [(100, 200), (100, 250)], [False, False], [True, True])

    def _single_region(self):
        return _region(100, 200, [(100, 200)], [True], [True])

    def test_full_inclusion_interior_block(self):
        c = _counter()
        region = self._shared_right_region()
        _classify(c, region, [(50, 80), (100, 200), (300, 350)])
        rec = _rec(c, region)
        assert _cnt(rec["cands"][(100, 200)]["full"]) == 1
        assert (150, 200) not in rec["cands"]
        assert _cnt(rec["excl"]) == 0 and _cnt(rec["amb"]) == 0

    def test_left_half_inclusion(self):
        c = _counter()
        region = self._shared_right_region()
        _classify(c, region, [(20, 40), (100, 150), (300, 350)])
        rec = _rec(c, region)
        assert _cnt(rec["cands"][(100, 200)]["left"]) == 1
        assert _cnt(rec["cands"][(100, 200)]["full"]) == 0
        assert _cnt(rec["amb"]) == 0

    def test_right_half_inclusion(self):
        c = _counter()
        region = self._shared_left_region()
        _classify(c, region, [(20, 40), (220, 250), (300, 350)])
        rec = _rec(c, region)
        assert _cnt(rec["cands"][(100, 250)]["right"]) == 1
        assert _cnt(rec["amb"]) == 0

    def test_half_at_shared_edge_is_ambiguous(self):
        c = _counter()
        region = self._shared_right_region()
        _classify(c, region, [(20, 40), (170, 200), (300, 350)])
        rec = _rec(c, region)
        assert _cnt(rec["amb"]) == 1
        assert (100, 200) not in rec["cands"] or _cnt(rec["cands"][(100, 200)]["right"]) == 0

    def test_region_exclusion(self):
        c = _counter(exclusion_margin=50)
        region = self._single_region()
        _classify(c, region, [(20, 40), (360, 400)])
        rec = _rec(c, region)
        assert _cnt(rec["excl"]) == 1
        assert not rec["cands"]

    def test_ambiguous_overlap_no_match(self):
        c = _counter()
        region = self._single_region()
        _classify(c, region, [(120, 180)])
        rec = _rec(c, region)
        assert _cnt(rec["amb"]) == 1

    def test_no_overlap_no_exclusion_skipped(self):
        c = _counter()
        region = self._single_region()
        _classify(c, region, [(500, 600)])
        assert _rec(c, region) is None

    def test_multiple_full_is_ambiguous(self):
        c = _counter()
        region = _region(100, 400, [(100, 200), (300, 400)], [True, True], [True, True])
        _classify(c, region, [(50, 70), (100, 200), (300, 400), (500, 520)])
        rec = _rec(c, region)
        assert _cnt(rec["amb"]) == 1
        assert (100, 200) not in rec["cands"] and (300, 400) not in rec["cands"]

    def test_terminal_full_requires_anchor(self):
        region = self._single_region()
        blocks = [(100, 200), (300, 350)]

        c_no = _counter()
        _classify(c_no, region, blocks, strand="+", tss=False)
        assert _cnt(_rec(c_no, region)["amb"]) == 1

        c_yes = _counter()
        _classify(c_yes, region, blocks, strand="+", tss=True)
        assert _cnt(_rec(c_yes, region)["cands"][(100, 200)]["full"]) == 1


# ── add_read_info: bulk (no barcode) + filters + anchoring ─────────────────

class TestAddReadInfo:
    def _gene_info_single(self):
        gi = SimpleNamespace()
        gi._exon_splice_site_regions = [_region(100, 200, [(100, 200)], [True], [True])]
        return gi

    def test_bulk_read_counted_without_barcode(self):
        c = _counter()
        gi = self._gene_info_single()
        c.add_read_info(_make_ra(gi, [(50, 80), (100, 200), (300, 350)]))
        assert _cnt(_rec(c, gi._exon_splice_site_regions[0])["cands"][(100, 200)]["full"]) == 1

    def test_strand_filter(self):
        c = _counter()
        gi = self._gene_info_single()
        c.add_read_info(_make_ra(gi, [(50, 80), (100, 200), (300, 350)], strand="-"))
        assert not c.regions

    def test_wrong_gene_not_counted(self):
        c = _counter()
        gi = self._gene_info_single()
        c.add_read_info(_make_ra(gi, [(50, 80), (100, 200), (300, 350)], gene="gene2"))
        assert not c.regions

    def test_tss_event_enables_terminal_full(self):
        gi = self._gene_info_single()
        blocks = [(100, 200), (300, 350)]
        c1 = _counter()
        c1.add_read_info(_make_ra(gi, blocks))
        assert _cnt(_rec(c1, gi._exon_splice_site_regions[0])["amb"]) == 1
        c2 = _counter()
        c2.add_read_info(_make_ra(gi, blocks, events=[MatchEventSubtype.terminal_site_match_left]))
        assert _cnt(_rec(c2, gi._exon_splice_site_regions[0])["cands"][(100, 200)]["full"]) == 1


# ── grouped counter: one row per (feature, group) ──────────────────────────

class TestGrouped:
    def test_one_row_per_group(self):
        sp = StringPoolManager()
        pool = sp.get_read_group_pool(0)
        g1, g2 = pool.add("cell1"), pool.add("cell2")
        c = _counter(string_pools=sp, group_index=0)
        region = _region(100, 200, [(100, 200)], [True], [True])
        # two full-inclusion reads from cell1, one from cell2
        _classify(c, region, [(50, 80), (100, 200), (300, 350)], group_id=g1)
        _classify(c, region, [(50, 80), (100, 200), (300, 350)], group_id=g1)
        _classify(c, region, [(50, 80), (100, 200), (300, 350)], group_id=g2)
        c.dump()
        with open(c.output_counts_file_name) as f:
            lines = f.read().strip().split("\n")
        assert lines[0] == "region_gene_candidate\tn_full\tn_left\tn_right\tgroup_id"
        rows = {(r[0], r[4]): r for r in (l.split("\t") for l in lines[1:])}
        feat = "chr1_100_200_+__gene1__chr1_100_200_+"
        assert rows[(feat, "cell1")][1] == "2"
        assert rows[(feat, "cell2")][1] == "1"


# ── output ─────────────────────────────────────────────────────────────────

class TestOutput:
    def test_ungrouped_dump_group_id_is_na(self):
        c = _counter()
        region = _region(100, 200, [(100, 200), (150, 200)], [True, True], [False, False])
        _classify(c, region, [(50, 80), (100, 200), (300, 350)])   # full inclusion
        _classify(c, region, [(20, 40), (360, 400)])               # exclusion
        c.dump()
        with open(c.output_counts_file_name) as f:
            lines = f.read().strip().split("\n")
        assert lines[0] == "region_gene_candidate\tn_full\tn_left\tn_right\tgroup_id"
        rows = {l.split("\t")[0]: l.split("\t") for l in lines[1:]}
        inc = rows["chr1_100_200_+__gene1__chr1_100_200_+"]
        assert inc[1:5] == ["1", "0", "0", "NA"]
        exc = rows["chr1_100_200_+__gene1__exclusion"]
        assert exc[1:5] == ["1", "0", "0", "NA"]

    def test_emit_read_ids(self):
        c = _counter(emit_read_ids=True)
        region = _region(100, 200, [(100, 200)], [True], [True])
        _classify(c, region, [(50, 80), (100, 200), (300, 350)], read_id="read_42")
        c.dump()
        with open(c.output_counts_file_name) as f:
            lines = f.read().strip().split("\n")
        assert lines[0].split("\t") == ["region_gene_candidate", "n_full", "n_left", "n_right",
                                        "group_id", "read_ids_full", "read_ids_left", "read_ids_right"]
        row = lines[1].split("\t")
        assert row[4] == "NA"                       # group_id
        assert row[5:8] == ["read_42", "NA", "NA"]  # read id lists

    def test_finalize_keeps_plain_tsv(self):
        c = _counter()
        region = _region(100, 200, [(100, 200)], [True], [True])
        _classify(c, region, [(50, 80), (100, 200), (300, 350)])
        c.dump()
        c.finalize()
        assert os.path.exists(c.output_counts_file_name)
        assert not os.path.exists(c.output_counts_file_name + ".gz")


# ── conversion script: per-group -> per-molecule group lists ───────────────

def _load_convert():
    path = os.path.join(os.path.dirname(__file__), "..", "misc",
                        "exon_splice_site_to_group_lists.py")
    spec = importlib.util.spec_from_file_location("ess_convert", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestConversionScript:
    def test_per_group_to_group_lists(self):
        conv = _load_convert()
        tmpdir = tempfile.mkdtemp()
        inp = os.path.join(tmpdir, "in.tsv")
        outp = os.path.join(tmpdir, "out.tsv")
        with open(inp, "w") as f:
            f.write("region_gene_candidate\tn_full\tn_left\tn_right\tgroup_id\n")
            f.write("F__inc\t2\t0\t0\tcell1\n")
            f.write("F__inc\t1\t0\t0\tcell2\n")
            f.write("F__exclusion\t3\t0\t0\tcell1\n")
        conv.convert(inp, outp)
        with open(outp) as f:
            lines = f.read().strip().split("\n")
        assert lines[0] == ("region_gene_candidate\tn_full\tn_left\tn_right"
                            "\tgroups_full\tgroups_left\tgroups_right")
        rows = {l.split("\t")[0]: l.split("\t") for l in lines[1:]}
        inc = rows["F__inc"]
        assert inc[1:4] == ["3", "0", "0"]
        assert inc[4] == "cell1;cell1;cell2"
        assert inc[5] == "NA" and inc[6] == "NA"
        exc = rows["F__exclusion"]
        assert exc[1] == "3" and exc[4] == "cell1;cell1;cell1"

    def test_conversion_carries_read_ids(self):
        conv = _load_convert()
        tmpdir = tempfile.mkdtemp()
        inp = os.path.join(tmpdir, "in.tsv")
        outp = os.path.join(tmpdir, "out.tsv")
        with open(inp, "w") as f:
            f.write("region_gene_candidate\tn_full\tn_left\tn_right\tgroup_id"
                    "\tread_ids_full\tread_ids_left\tread_ids_right\n")
            f.write("F__inc\t2\t0\t0\tcell1\tr1;r2\tNA\tNA\n")
        conv.convert(inp, outp)
        with open(outp) as f:
            lines = f.read().strip().split("\n")
        assert lines[0].endswith("read_ids_full\tread_ids_left\tread_ids_right")
        row = lines[1].split("\t")
        assert row[4] == "cell1;cell1"   # groups_full
        assert row[7] == "r1;r2"         # read_ids_full
        assert row[8] == "NA" and row[9] == "NA"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
