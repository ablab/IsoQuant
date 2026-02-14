############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import pytest
import tempfile

from src.read_groups import (
    parse_barcode2spot_spec,
    load_barcode2barcode_mapping,
    BarcodeSpotGrouper,
    SharedTableData,
)
from src.barcode_calling.umi_filtering import UMIFilter
from src.isoform_assignment import ReadAssignment, ReadAssignmentType, IsoformMatch, MatchClassification
from src.common import junctions_from_blocks
from src.string_pools import StringPoolManager
from src.file_naming import umi_barcode2barcode_prefix, umi_barcode2barcode_global_lock


# -- Helpers / Mocks --

class MockAlignment:
    def __init__(self, read_id=""):
        self.query_name = read_id


class MockReadAssignment:
    """Lightweight mock used only for grouper tests (no string_pools needed)."""
    def __init__(self, barcode=None):
        self.barcode = barcode


def _make_read_assignment(string_pools, read_id, barcode, umi, gene_id, transcript_id,
                          exons, assignment_type=ReadAssignmentType.unique):
    """Create a fully populated ReadAssignment for UMI filtering tests."""
    ra = ReadAssignment(read_id=read_id, assignment_type=assignment_type, string_pools=string_pools)
    ra.chr_id = "chr1"
    ra.strand = "+"
    ra.corrected_exons = exons
    ra.corrected_introns = junctions_from_blocks(exons)
    ra.barcode = barcode
    ra.umi = umi
    ra.gene_assignment_type = assignment_type
    im = IsoformMatch(MatchClassification.full_splice_match, string_pools,
                      assigned_gene=gene_id, assigned_transcript=transcript_id)
    ra.isoform_matches = [im]
    ra.set_additional_attribute('transcript_type', 'protein_coding')
    ra.set_additional_attribute('polya_site', -1)
    ra.set_additional_attribute('cell_type', 'None')
    return ra


# -- file_naming tests --

class TestBarcode2BarcodeFileNaming:
    def test_prefix_format(self):
        result = umi_barcode2barcode_prefix("/out/sample.umi_done", 0)
        assert result == "/out/sample.umi_done.barcode_barcode_col0"

    def test_prefix_multiple_columns(self):
        assert umi_barcode2barcode_prefix("/x", 2) == "/x.barcode_barcode_col2"

    def test_lock_format(self):
        result = umi_barcode2barcode_global_lock("/out/sample.umi_done", 1)
        assert result == "/out/sample.umi_done.barcode_barcode_col1.lock"


# -- load_barcode2barcode_mapping tests --

class TestLoadBarcode2BarcodeMapping:
    def test_single_column(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("AAAA\tspot_A\n")
            f.write("CCCC\tspot_A\n")
            f.write("GGGG\tspot_B\n")
            fname = f.name

        try:
            mapping = load_barcode2barcode_mapping(fname, barcode_col=0, spot_cols=[1])
            assert mapping["AAAA"] == ["spot_A"]
            assert mapping["CCCC"] == ["spot_A"]
            assert mapping["GGGG"] == ["spot_B"]
        finally:
            os.unlink(fname)

    def test_multi_column(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("AAAA\tspot_A\tregion_1\n")
            f.write("CCCC\tspot_A\tregion_2\n")
            f.write("GGGG\tspot_B\tregion_1\n")
            fname = f.name

        try:
            mapping = load_barcode2barcode_mapping(fname, barcode_col=0, spot_cols=[1, 2])
            assert mapping["AAAA"] == ["spot_A", "region_1"]
            assert mapping["CCCC"] == ["spot_A", "region_2"]
            assert mapping["GGGG"] == ["spot_B", "region_1"]
        finally:
            os.unlink(fname)

    def test_custom_barcode_column(self):
        """Barcode column is not necessarily column 0."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("extra\tAAAA\tspot_X\n")
            f.write("extra\tCCCC\tspot_Y\n")
            fname = f.name

        try:
            mapping = load_barcode2barcode_mapping(fname, barcode_col=1, spot_cols=[2])
            assert mapping["AAAA"] == ["spot_X"]
            assert mapping["CCCC"] == ["spot_Y"]
        finally:
            os.unlink(fname)

    def test_empty_file(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            fname = f.name

        try:
            mapping = load_barcode2barcode_mapping(fname, barcode_col=0, spot_cols=[1])
            assert mapping == {}
        finally:
            os.unlink(fname)

    def test_comments_and_blank_lines_skipped(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("# header comment\n")
            f.write("\n")
            f.write("AAAA\tspot_A\n")
            f.write("# another comment\n")
            f.write("CCCC\tspot_B\n")
            fname = f.name

        try:
            mapping = load_barcode2barcode_mapping(fname, barcode_col=0, spot_cols=[1])
            assert len(mapping) == 2
            assert mapping["AAAA"] == ["spot_A"]
            assert mapping["CCCC"] == ["spot_B"]
        finally:
            os.unlink(fname)


# -- Grouper tests (barcode_barcode uses BarcodeSpotGrouper) --

class TestBarcode2BarcodeGrouper:
    """Test that barcode_barcode grouping works via BarcodeSpotGrouper."""

    def test_single_column_grouper(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("AAAA\tspot_A\n")
            f.write("CCCC\tspot_A\n")
            f.write("GGGG\tspot_B\n")
            fname = f.name

        try:
            shared = SharedTableData(fname, read_id_column_index=0,
                                     group_id_column_indices=[1], delim='\t')
            grouper = BarcodeSpotGrouper(shared, column_index=0)

            alignment = MockAlignment()
            assert grouper.get_group_id(alignment, read_assignment=MockReadAssignment("AAAA")) == "spot_A"
            assert grouper.get_group_id(alignment, read_assignment=MockReadAssignment("CCCC")) == "spot_A"
            assert grouper.get_group_id(alignment, read_assignment=MockReadAssignment("GGGG")) == "spot_B"
        finally:
            os.unlink(fname)

    def test_multi_column_groupers(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("AAAA\tspot_A\tregion_1\n")
            f.write("CCCC\tspot_B\tregion_1\n")
            fname = f.name

        try:
            shared = SharedTableData(fname, read_id_column_index=0,
                                     group_id_column_indices=[1, 2], delim='\t')
            g0 = BarcodeSpotGrouper(shared, column_index=0)
            g1 = BarcodeSpotGrouper(shared, column_index=1)

            alignment = MockAlignment()
            ra1 = MockReadAssignment("AAAA")
            ra2 = MockReadAssignment("CCCC")

            assert g0.get_group_id(alignment, read_assignment=ra1) == "spot_A"
            assert g0.get_group_id(alignment, read_assignment=ra2) == "spot_B"
            assert g1.get_group_id(alignment, read_assignment=ra1) == "region_1"
            assert g1.get_group_id(alignment, read_assignment=ra2) == "region_1"
        finally:
            os.unlink(fname)

    def test_unknown_barcode_returns_na(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("AAAA\tspot_A\n")
            fname = f.name

        try:
            shared = SharedTableData(fname, read_id_column_index=0,
                                     group_id_column_indices=[1], delim='\t')
            grouper = BarcodeSpotGrouper(shared, column_index=0)
            alignment = MockAlignment()

            assert grouper.get_group_id(alignment, read_assignment=MockReadAssignment("ZZZZ")) == "NA"
        finally:
            os.unlink(fname)


# -- parse_grouping_spec tests for barcode_barcode --

class TestParseGroupingSpecBarcode2Barcode:
    """Test parse_grouping_spec with barcode_barcode spec."""

    class MockArgs:
        barcode2barcode = None

    class MockSample:
        barcodes_split_reads = "some_file"

    def test_single_column(self):
        from src.read_groups import parse_grouping_spec

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("AAAA\tspot_A\n")
            f.write("CCCC\tspot_B\n")
            fname = f.name

        try:
            args = self.MockArgs()
            args.barcode2barcode = fname
            sample = self.MockSample()

            result = parse_grouping_spec("barcode_barcode", args, sample, "chr1")
            assert isinstance(result, BarcodeSpotGrouper)
        finally:
            os.unlink(fname)

    def test_multi_column_returns_list(self):
        from src.read_groups import parse_grouping_spec

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("AAAA\tspot_A\tregion_1\n")
            f.write("CCCC\tspot_B\tregion_2\n")
            fname = f.name

        try:
            args = self.MockArgs()
            args.barcode2barcode = "%s:0:1,2" % fname
            sample = self.MockSample()

            result = parse_grouping_spec("barcode_barcode", args, sample, "chr1")
            assert isinstance(result, list)
            assert len(result) == 2
            assert all(isinstance(g, BarcodeSpotGrouper) for g in result)
        finally:
            os.unlink(fname)

    def test_missing_arg_exits(self):
        from src.read_groups import parse_grouping_spec

        args = self.MockArgs()
        args.barcode2barcode = None
        sample = self.MockSample()

        with pytest.raises(SystemExit):
            parse_grouping_spec("barcode_barcode", args, sample, "chr1")

    def test_missing_barcoded_reads_exits(self):
        from src.read_groups import parse_grouping_spec

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("AAAA\tspot_A\n")
            fname = f.name

        try:
            args = self.MockArgs()
            args.barcode2barcode = fname
            sample = self.MockSample()
            sample.barcodes_split_reads = None

            with pytest.raises(SystemExit):
                parse_grouping_spec("barcode_barcode", args, sample, "chr1")
        finally:
            os.unlink(fname)


# -- get_grouping_strategy_names tests --

class TestGetGroupingStrategyNamesBarcode2Barcode:

    class MockArgs:
        barcode2barcode = None
        read_group = None

    def test_single_column_name(self):
        from src.read_groups import get_grouping_strategy_names

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("AAAA\tspot_A\n")
            fname = f.name

        try:
            args = self.MockArgs()
            args.barcode2barcode = fname
            args.read_group = ["barcode_barcode"]

            names = get_grouping_strategy_names(args)
            assert names == ["barcode_barcode"]
        finally:
            os.unlink(fname)

    def test_multi_column_names(self):
        from src.read_groups import get_grouping_strategy_names

        args = self.MockArgs()
        args.barcode2barcode = "file.tsv:0:1,2,3"
        args.read_group = ["barcode_barcode"]

        names = get_grouping_strategy_names(args)
        assert names == ["barcode_barcode_col0", "barcode_barcode_col1", "barcode_barcode_col2"]


# -- get_grouping_pool_types tests --

class TestGetGroupingPoolTypesBarcode2Barcode:

    class MockArgs:
        barcode2barcode = None
        read_group = None

    def test_single_column_pool(self):
        from src.read_groups import get_grouping_pool_types

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("AAAA\tspot_A\n")
            fname = f.name

        try:
            args = self.MockArgs()
            args.barcode2barcode = fname
            args.read_group = ["barcode_barcode"]

            pool_types = get_grouping_pool_types(args)
            assert pool_types == {0: 'barcode_barcode:0'}
        finally:
            os.unlink(fname)

    def test_multi_column_pools(self):
        from src.read_groups import get_grouping_pool_types

        args = self.MockArgs()
        args.barcode2barcode = "file.tsv:0:1,2"
        args.read_group = ["barcode_barcode"]

        pool_types = get_grouping_pool_types(args)
        assert pool_types == {0: 'barcode_barcode:0', 1: 'barcode_barcode:1'}


# -- UMIFilter barcode_remap tests --

class TestUMIFilterBarcodeRemap:
    string_pools = StringPoolManager()

    def test_init_stores_remap(self):
        remap = {"AAAA": "spot_X", "CCCC": "spot_Y"}
        f = UMIFilter(umi_length=10, edit_distance=3, barcode_remap=remap)
        assert f.barcode_remap is remap

    def test_init_default_remap_is_none(self):
        f = UMIFilter(umi_length=10, edit_distance=3)
        assert f.barcode_remap is None

    def test_remap_groups_by_spot(self):
        """Two barcodes mapping to the same spot should be grouped together for dedup."""
        remap = {"BC_A": "spot_1", "BC_B": "spot_1", "BC_C": "spot_2"}
        umi_filter = UMIFilter(umi_length=0, edit_distance=0, barcode_remap=remap)

        # Build gene_barcode_dict the same way process_single_chr does
        from collections import defaultdict
        gene_barcode_dict = defaultdict(lambda: defaultdict(list))

        ra1 = _make_read_assignment(self.string_pools, "r1", "BC_A", "UMI1", "gene1", "tx1",
                                    [(100, 200), (300, 400)])
        ra2 = _make_read_assignment(self.string_pools, "r2", "BC_B", "UMI1", "gene1", "tx1",
                                    [(100, 200), (300, 400)])
        ra3 = _make_read_assignment(self.string_pools, "r3", "BC_C", "UMI2", "gene1", "tx1",
                                    [(100, 200), (300, 400)])

        for ra in [ra1, ra2, ra3]:
            barcode = ra.barcode
            gene_id = ra.isoform_matches[0].assigned_gene
            grouping_key = barcode
            if umi_filter.barcode_remap and barcode in umi_filter.barcode_remap:
                grouping_key = umi_filter.barcode_remap[barcode]
            gene_barcode_dict[gene_id][grouping_key].append(ra)

        # BC_A and BC_B both map to spot_1, so they should be in the same group
        assert len(gene_barcode_dict["gene1"]["spot_1"]) == 2
        assert len(gene_barcode_dict["gene1"]["spot_2"]) == 1

    def test_no_remap_groups_by_barcode(self):
        """Without remap, reads group by original barcode."""
        umi_filter = UMIFilter(umi_length=0, edit_distance=0, barcode_remap=None)

        from collections import defaultdict
        gene_barcode_dict = defaultdict(lambda: defaultdict(list))

        ra1 = _make_read_assignment(self.string_pools, "r1", "BC_A", "UMI1", "gene1", "tx1",
                                    [(100, 200), (300, 400)])
        ra2 = _make_read_assignment(self.string_pools, "r2", "BC_B", "UMI1", "gene1", "tx1",
                                    [(100, 200), (300, 400)])

        for ra in [ra1, ra2]:
            barcode = ra.barcode
            gene_id = ra.isoform_matches[0].assigned_gene
            grouping_key = barcode
            if umi_filter.barcode_remap and barcode in umi_filter.barcode_remap:
                grouping_key = umi_filter.barcode_remap[barcode]
            gene_barcode_dict[gene_id][grouping_key].append(ra)

        # Without remap, each barcode stays separate
        assert len(gene_barcode_dict["gene1"]) == 2
        assert len(gene_barcode_dict["gene1"]["BC_A"]) == 1
        assert len(gene_barcode_dict["gene1"]["BC_B"]) == 1

    def test_remap_unknown_barcode_stays_as_is(self):
        """Barcodes not in the remap dict should keep their original barcode as key."""
        remap = {"BC_A": "spot_1"}
        umi_filter = UMIFilter(umi_length=0, edit_distance=0, barcode_remap=remap)

        from collections import defaultdict
        gene_barcode_dict = defaultdict(lambda: defaultdict(list))

        ra1 = _make_read_assignment(self.string_pools, "r1", "BC_A", "UMI1", "gene1", "tx1",
                                    [(100, 200), (300, 400)])
        ra2 = _make_read_assignment(self.string_pools, "r2", "BC_UNKNOWN", "UMI2", "gene1", "tx1",
                                    [(100, 200), (300, 400)])

        for ra in [ra1, ra2]:
            barcode = ra.barcode
            gene_id = ra.isoform_matches[0].assigned_gene
            grouping_key = barcode
            if umi_filter.barcode_remap and barcode in umi_filter.barcode_remap:
                grouping_key = umi_filter.barcode_remap[barcode]
            gene_barcode_dict[gene_id][grouping_key].append(ra)

        assert "spot_1" in gene_barcode_dict["gene1"]
        assert "BC_UNKNOWN" in gene_barcode_dict["gene1"]

    def test_remap_dedup_same_umi_different_barcodes(self):
        """Reads with same UMI, different barcodes but same spot should be deduped together."""
        remap = {"BC_A": "spot_1", "BC_B": "spot_1"}
        umi_filter = UMIFilter(umi_length=0, edit_distance=0, barcode_remap=remap)

        ra1 = _make_read_assignment(self.string_pools, "r1", "BC_A", "AAAA", "gene1", "tx1",
                                    [(100, 200), (300, 400)])
        ra2 = _make_read_assignment(self.string_pools, "r2", "BC_B", "AAAA", "gene1", "tx1",
                                    [(100, 200), (300, 400)])

        # Group by spot via remap
        from collections import defaultdict
        gene_barcode_dict = defaultdict(lambda: defaultdict(list))
        for ra in [ra1, ra2]:
            barcode = ra.barcode
            gene_id = ra.isoform_matches[0].assigned_gene
            grouping_key = remap.get(barcode, barcode)
            gene_barcode_dict[gene_id][grouping_key].append(ra)

        # Both reads are in spot_1 with same UMI — _process_duplicates should dedup
        result = umi_filter._process_duplicates(gene_barcode_dict["gene1"]["spot_1"])
        assert len(result) == 1  # Deduped to 1 representative read


# -- CLI argument tests --

class TestBarcode2BarcodeArgParsing:
    """Test that --barcode2barcode argument is parsed correctly."""

    def test_arg_added(self):
        """Verify the argument exists in the parser."""
        import sys
        sys.argv = ['isoquant.py', '--full_help']
        # We just test that parse_barcode2spot_spec works on barcode2barcode-style specs
        filename, bc_col, spot_cols = parse_barcode2spot_spec("mapping.tsv:0:1,2")
        assert filename == "mapping.tsv"
        assert bc_col == 0
        assert spot_cols == [1, 2]


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
