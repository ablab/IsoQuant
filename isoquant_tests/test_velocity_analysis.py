############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""When RNA velocity runs: the velocity analysis, --barcode2cellid and the barcode count gate."""

import argparse
import logging
from types import SimpleNamespace

import pytest

import isoquant
from isoquant_lib.assignment.assignment_aggregator import ReadAssignmentAggregator
from isoquant_lib.assignment.read_groups import (
    BarcodeSpotGrouper,
    get_grouping_pool_types,
    get_grouping_strategy_names,
    parse_grouping_spec,
)
from isoquant_lib.dataset_processor import DatasetProcessor
from isoquant_lib.modes import AnalysisType, IsoQuantMode, MAX_DEFAULT_VELOCITY_CELLS
from isoquant_lib.utils.string_pools import StringPoolManager


def resolve(mode: str, analysis=None, genedb="annotation.db", barcode2cellid=None) -> argparse.Namespace:
    args = argparse.Namespace(mode=IsoQuantMode[mode], analysis=analysis, genedb=genedb,
                              barcode2cellid=barcode2cellid, count_exons=False,
                              count_intron_retentions=False, no_model_construction=False)
    isoquant.resolve_analyses(args)
    return args


class TestResolveVelocityAnalysis:
    @pytest.mark.parametrize("mode", ["tenX_v3", "tenX_v2"])
    def test_on_by_default_for_10x(self, mode):
        args = resolve(mode)
        assert args.rna_velocity
        assert not args.velocity_explicit

    @pytest.mark.parametrize("mode", ["stereoseq", "visium_hd", "curio", "visium_5prime", "custom_sc", "bulk"])
    def test_off_by_default_elsewhere(self, mode):
        assert not resolve(mode).rna_velocity

    @pytest.mark.parametrize("mode", ["stereoseq", "visium_hd", "bulk"])
    def test_barcode2cellid_turns_it_on(self, mode):
        assert resolve(mode, barcode2cellid="cells.tsv").rna_velocity

    def test_explicit_analysis_list_without_velocity_wins(self):
        assert not resolve("tenX_v3", analysis=["quant"]).rna_velocity
        assert not resolve("stereoseq", analysis=["quant"], barcode2cellid="cells.tsv").rna_velocity

    def test_explicit_velocity_is_marked_explicit(self):
        args = resolve("stereoseq", analysis=["quant", "velocity"])
        assert args.rna_velocity
        assert args.velocity_explicit
        assert AnalysisType.velocity in args.analyses

    def test_bulk_without_cells_drops_velocity(self, caplog):
        with caplog.at_level(logging.WARNING, logger="IsoQuant"):
            args = resolve("bulk", analysis=["velocity"])
        assert not args.rna_velocity
        assert "velocity requires" in caplog.text

    def test_requires_genedb(self):
        assert not resolve("tenX_v3", genedb=None).rna_velocity


class TestBarcodeCellGrouping:
    def test_pool_types_and_strategy_names_stay_aligned(self):
        args = SimpleNamespace(read_group=["barcode", "barcode_spot", "barcode_cell"],
                               barcode2spot="spots.tsv:0:1,2", barcode2cellid="cells.tsv",
                               barcode2barcode=None)
        pool_types = get_grouping_pool_types(args)
        names = get_grouping_strategy_names(args)
        assert names == ["barcode", "barcode_spot_col0", "barcode_spot_col1", "barcode_cell"]
        assert pool_types == {0: "barcode", 1: "barcode_spot:0", 2: "barcode_spot:1", 3: "barcode_cell"}

    def test_auto_added_by_barcode2cellid(self):
        args = SimpleNamespace(read_group=None, barcode2cellid="cells.tsv", barcode2spot=None,
                               barcode2barcode=None, mode=IsoQuantMode.stereoseq, large_output=None,
                               input_data=SimpleNamespace(has_replicas=lambda: False),
                               matching_strategy="default", model_construction_strategy="default_ont",
                               splice_correction_strategy="default_ont", fl_data=False)
        isoquant.set_data_dependent_options(args)
        # barcode_cell counts as barcode grouping, so plain barcode grouping is not added on top
        assert args.read_group == ["barcode_cell"]

    def test_grouper_maps_barcode_to_cell(self, tmp_path):
        mapping = tmp_path / "cells.tsv"
        mapping.write_text("AAAA\tcell_1\nCCCC\tcell_2\n")
        args = SimpleNamespace(barcode2cellid=str(mapping), barcoded_bam=True)
        sample = SimpleNamespace(barcodes_split_reads=None)
        grouper = parse_grouping_spec("barcode_cell", args, sample, "chr1")
        assert isinstance(grouper, BarcodeSpotGrouper)
        alignment = SimpleNamespace(query_name="read1")
        assert grouper.get_group_id(alignment, read_assignment=SimpleNamespace(barcode="CCCC")) == "cell_2"

    def test_string_pool_holds_cell_ids(self, tmp_path):
        mapping = tmp_path / "cells.tsv"
        mapping.write_text("AAAA\tcell_1\nCCCC\tcell_2\nGGGG\tcell_1\n")
        pools = StringPoolManager()
        pools.build_barcode_cell_pool(str(mapping))
        pools.set_group_spec_pool_type(0, "barcode_cell")
        pool = pools.get_read_group_pool(0)
        assert pool is pools.barcode_cell_pool
        assert sorted(pool.get_str(i) for i in range(len(pool))) == ["NA", "cell_1", "cell_2"]


def velocity_pool_type(read_group, **flags):
    args = SimpleNamespace(read_group=read_group, barcode2spot=None, barcode2barcode=None,
                           barcode2cellid="cells.tsv", **flags)
    aggregator = ReadAssignmentAggregator.__new__(ReadAssignmentAggregator)
    aggregator.args = args
    return aggregator._velocity_pool_type(get_grouping_pool_types(args))


class TestVelocityGrouping:
    def test_barcode_cell_preferred_over_barcode(self):
        assert velocity_pool_type(["barcode", "barcode_cell"], rna_velocity=True) == "barcode_cell"

    def test_barcode_otherwise(self):
        assert velocity_pool_type(["barcode", "barcode_spot"], rna_velocity=True) == "barcode"

    def test_never_spots(self):
        # 'barcode' pool type never equals barcode_spot / barcode_barcode pool types
        pool_type = velocity_pool_type(["barcode_spot"], rna_velocity=True)
        assert pool_type == "barcode"
        assert pool_type not in get_grouping_pool_types(
            SimpleNamespace(read_group=["barcode_spot"], barcode2spot=None, barcode2barcode=None)).values()

    def test_off_when_analysis_off(self):
        assert velocity_pool_type(["barcode"], rna_velocity=False) is None

    def test_per_sample_decision_overrides_analysis(self):
        assert velocity_pool_type(["barcode"], rna_velocity=True, rna_velocity_active=False) is None


def velocity_is_active(barcode_count, strategies=("barcode",), **flags) -> bool:
    processor = DatasetProcessor.__new__(DatasetProcessor)
    processor.args = SimpleNamespace(**flags)
    processor.grouping_strategy_names = list(strategies)
    return processor.velocity_is_active(barcode_count)


class TestBarcodeCountGate:
    def test_default_under_limit(self):
        assert velocity_is_active(MAX_DEFAULT_VELOCITY_CELLS, rna_velocity=True)

    def test_default_over_limit_is_skipped(self, caplog):
        with caplog.at_level(logging.WARNING, logger="IsoQuant"):
            assert not velocity_is_active(MAX_DEFAULT_VELOCITY_CELLS + 1, rna_velocity=True)
        assert "--barcode2cellid" in caplog.text

    def test_explicit_bypasses_limit(self):
        assert velocity_is_active(MAX_DEFAULT_VELOCITY_CELLS + 1, rna_velocity=True, velocity_explicit=True)

    def test_barcode_cell_bypasses_limit(self):
        assert velocity_is_active(MAX_DEFAULT_VELOCITY_CELLS + 1, strategies=("barcode", "barcode_cell"),
                                  rna_velocity=True)

    def test_analysis_off(self):
        assert not velocity_is_active(10, rna_velocity=False)
