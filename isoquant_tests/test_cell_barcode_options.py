############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Resolution of --n_cells and --barcode_correction."""

import argparse

import pytest

from isoquant_lib.barcode_calling import options
from isoquant_lib.modes import AUTO_BARCODES, BarcodeCorrectionMethod, IsoQuantMode


def make_args(mode="tenX_v3", n_cells=None, barcode_whitelist=None,
              barcoded_reads=None, barcoded_bam=False,
              barcode_correction=BarcodeCorrectionMethod.auto.name):
    return argparse.Namespace(mode=IsoQuantMode[mode], n_cells=n_cells,
                              barcode_whitelist=barcode_whitelist,
                              barcoded_reads=barcoded_reads, barcoded_bam=barcoded_bam,
                              barcode_correction=barcode_correction)


def whitelist_file(tmp_path, count=10):
    path = tmp_path / "whitelist.txt"
    path.write_text("".join("A" * 15 + "%s\n" % base for base in "ACGT" * count))
    return [str(path)]


class TestResolveNCells:
    @pytest.mark.parametrize("value, expected", [(None, None), ("auto", AUTO_BARCODES),
                                                 ("5000", 5000)])
    def test_accepted_values(self, value, expected):
        args = make_args(n_cells=value)
        assert options.resolve_n_cells(args) == expected

    @pytest.mark.parametrize("value", ["0", "-1", "many", "1e3"])
    def test_rejected_values(self, value):
        with pytest.raises(SystemExit) as excinfo:
            options.resolve_n_cells(make_args(n_cells=value))
        assert excinfo.value.code != 0


class TestDetectionIsRequested:
    def test_whitelist_alone_is_the_cell_list(self, tmp_path):
        args = make_args(barcode_whitelist=whitelist_file(tmp_path))
        options.resolve_barcode_correction(args)
        assert args.detect_cell_barcodes is False

    def test_n_cells_turns_the_whitelist_into_a_pool(self, tmp_path):
        args = make_args(n_cells="5000", barcode_whitelist=whitelist_file(tmp_path))
        options.resolve_barcode_correction(args)
        assert args.detect_cell_barcodes is True
        assert args.n_cells == 5000

    def test_auto_whitelist_detects_without_a_pool(self):
        args = make_args(barcode_whitelist=[AUTO_BARCODES])
        options.resolve_barcode_correction(args)
        assert args.detect_cell_barcodes is True
        assert args.n_cells == AUTO_BARCODES

    def test_detect_forces_detection_without_n_cells(self, tmp_path):
        args = make_args(barcode_whitelist=whitelist_file(tmp_path),
                         barcode_correction=BarcodeCorrectionMethod.detect.name)
        options.resolve_barcode_correction(args)
        assert args.detect_cell_barcodes is True
        assert args.n_cells == AUTO_BARCODES

    def test_method_names_describe_what_they_do(self):
        """`graph` was left over from the dropped edit-distance port; it corrects nothing."""
        assert {e.name for e in BarcodeCorrectionMethod} == {"whitelist", "detect", "auto"}

    def test_unsupported_mode_aborts(self, tmp_path):
        args = make_args(mode="curio", n_cells="5000",
                         barcode_whitelist=whitelist_file(tmp_path))
        with pytest.raises(SystemExit) as excinfo:
            options.resolve_barcode_correction(args)
        assert excinfo.value.code != 0


class TestNCellsIgnoredWarnings:
    """Silently dropping --n_cells leaves the user thinking cells were detected."""

    def test_barcoded_bam_warns(self, caplog):
        args = make_args(n_cells="5000", barcoded_bam=True)
        options.resolve_barcode_correction(args)
        assert any("--n_cells is ignored" in r.message for r in caplog.records)
        assert args.detect_cell_barcodes is False

    def test_barcoded_reads_warns(self, caplog):
        args = make_args(n_cells="5000", barcoded_reads=["reads.tsv"])
        options.resolve_barcode_correction(args)
        assert any("--n_cells is ignored" in r.message for r in caplog.records)
        assert args.detect_cell_barcodes is False

    def test_explicit_whitelist_method_warns(self, caplog, tmp_path):
        args = make_args(n_cells="5000", barcode_whitelist=whitelist_file(tmp_path),
                         barcode_correction=BarcodeCorrectionMethod.whitelist.name)
        options.resolve_barcode_correction(args)
        assert any("--n_cells is ignored" in r.message for r in caplog.records)
        assert args.detect_cell_barcodes is False

    def test_no_warning_when_n_cells_is_honoured(self, caplog, tmp_path):
        args = make_args(n_cells="5000", barcode_whitelist=whitelist_file(tmp_path))
        options.resolve_barcode_correction(args)
        assert not any("--n_cells is ignored" in r.message for r in caplog.records)


class TestBothEntryPointsAgree:
    """isoquant.py and isoquant_detect_barcodes.py used to carry separate copies of this."""

    @staticmethod
    def standalone_args(**kwargs):
        """The standalone tool names the whitelist --barcodes and has no --barcoded_*."""
        args = make_args(**kwargs)
        args.barcodes = args.barcode_whitelist
        del args.barcode_whitelist
        del args.barcoded_reads
        del args.barcoded_bam
        return args

    @pytest.mark.parametrize("kwargs", [
        {},
        {"n_cells": "5000"},
        {"n_cells": "auto"},
        {"barcode_whitelist": [AUTO_BARCODES]},
        {"barcode_correction": BarcodeCorrectionMethod.detect.name},
    ])
    def test_same_decision_under_either_option_name(self, kwargs, tmp_path):
        if "barcode_whitelist" not in kwargs:
            kwargs = dict(kwargs, barcode_whitelist=whitelist_file(tmp_path))
        pipeline_args = make_args(**kwargs)
        standalone = self.standalone_args(**kwargs)

        options.resolve_barcode_correction(pipeline_args)
        options.resolve_barcode_correction(standalone, whitelist_option="--barcodes")

        assert pipeline_args.detect_cell_barcodes == standalone.detect_cell_barcodes
        assert pipeline_args.n_cells == standalone.n_cells
