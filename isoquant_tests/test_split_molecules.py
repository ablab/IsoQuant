############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Resolution of --split_molecules and the superseded split mode names."""

import argparse

import pytest

import isoquant
from isoquant_lib.modes import (
    DEPRECATED_MODE_ALIASES,
    IsoQuantMode,
    SPLIT_MOLECULES_AUTO,
    SPLIT_MOLECULES_FALSE,
    SPLIT_MOLECULES_TRUE,
)

SPLITTING_MODES = ["tenX_v3", "tenX_v2", "stereoseq", "visium_5prime"]
NON_SPLITTING_MODES = ["curio", "visium_hd", "custom_sc"]


def resolve(mode, split_molecules=None):
    """Run the same two steps check_input_params does, and return (mode, split flag)."""
    args = argparse.Namespace(mode=mode, split_molecules=split_molecules)
    isoquant._resolve_deprecated_mode(args)
    args.mode = IsoQuantMode[args.mode]
    isoquant._resolve_split_molecules(args)
    return args.mode, args.split_molecules


class TestSplitMoleculesResolution:
    @pytest.mark.parametrize("mode", SPLITTING_MODES)
    def test_auto_splits_where_supported(self, mode):
        assert resolve(mode)[1] is True
        assert resolve(mode, SPLIT_MOLECULES_AUTO)[1] is True

    @pytest.mark.parametrize("mode", NON_SPLITTING_MODES)
    def test_auto_is_silent_where_unsupported(self, mode):
        """auto must not fail for a protocol that cannot split -- it just does not split."""
        assert resolve(mode)[1] is False
        assert resolve(mode, SPLIT_MOLECULES_AUTO)[1] is False

    @pytest.mark.parametrize("mode", SPLITTING_MODES)
    def test_true_splits_where_supported(self, mode):
        assert resolve(mode, SPLIT_MOLECULES_TRUE)[1] is True

    @pytest.mark.parametrize("mode", NON_SPLITTING_MODES + ["bulk"])
    def test_true_aborts_where_unsupported(self, mode):
        """Asking for the impossible must fail rather than quietly produce unsplit output."""
        with pytest.raises(SystemExit) as excinfo:
            resolve(mode, SPLIT_MOLECULES_TRUE)
        assert excinfo.value.code != 0

    @pytest.mark.parametrize("mode", SPLITTING_MODES + NON_SPLITTING_MODES)
    def test_false_never_splits(self, mode):
        assert resolve(mode, SPLIT_MOLECULES_FALSE)[1] is False


class TestDeprecatedModeAliases:
    @pytest.mark.parametrize("alias, expected_mode, expected_split", [
        ("tenX_v3_split", IsoQuantMode.tenX_v3, True),
        ("tenX_v2_split", IsoQuantMode.tenX_v2, True),
        ("stereoseq_nosplit", IsoQuantMode.stereoseq, False),
    ])
    def test_alias_maps_to_mode_and_split(self, alias, expected_mode, expected_split):
        assert resolve(alias) == (expected_mode, expected_split)

    def test_alias_covers_every_removed_name(self):
        assert set(DEPRECATED_MODE_ALIASES) == {"tenX_v3_split", "tenX_v2_split", "stereoseq_nosplit"}

    @pytest.mark.parametrize("alias", list(DEPRECATED_MODE_ALIASES))
    def test_explicit_flag_overrides_the_alias(self, alias):
        """The alias only supplies a default."""
        assert resolve(alias, SPLIT_MOLECULES_FALSE)[1] is False
        assert resolve(alias, SPLIT_MOLECULES_TRUE)[1] is True

    def test_alias_warns(self, caplog):
        resolve("stereoseq_nosplit")
        assert any("deprecated" in r.message for r in caplog.records)

    @pytest.mark.parametrize("alias", list(DEPRECATED_MODE_ALIASES))
    def test_alias_reproduces_the_old_behaviour(self, alias):
        """Each alias must resolve to what the removed mode used to do."""
        old_produced_split_fasta = {"tenX_v3_split": True, "tenX_v2_split": True,
                                    "stereoseq_nosplit": False}
        assert resolve(alias)[1] is old_produced_split_fasta[alias]


class TestDetectorSelection:
    @pytest.mark.parametrize("mode", SPLITTING_MODES)
    def test_registry_returns_a_splitting_detector(self, mode):
        from isoquant_lib.barcode_calling.detect_barcodes import (
            BARCODE_CALLING_MODES, barcode_detector_class)
        enum_mode = IsoQuantMode[mode]
        split_cls = barcode_detector_class(enum_mode, True)
        plain_cls = barcode_detector_class(enum_mode, False)
        assert split_cls is not plain_cls
        assert plain_cls is BARCODE_CALLING_MODES[enum_mode]

    @pytest.mark.parametrize("mode", NON_SPLITTING_MODES)
    def test_non_splitting_modes_keep_their_detector(self, mode):
        from isoquant_lib.barcode_calling.detect_barcodes import (
            BARCODE_CALLING_MODES, barcode_detector_class)
        enum_mode = IsoQuantMode[mode]
        assert barcode_detector_class(enum_mode, False) is BARCODE_CALLING_MODES[enum_mode]

    @pytest.mark.parametrize("mode, split, expected", [
        # the exact detectors the removed mode names used to select
        ("tenX_v3", True, "TenXSplittingBarcodeDetector"),
        ("tenX_v3", False, "TenXBarcodeDetector"),
        ("tenX_v2", True, "TenXv2SplittingBarcodeDetector"),
        ("tenX_v2", False, "TenXv2BarcodeDetector"),
        ("stereoseq", True, "SharedMemoryStereoSplittingBarcodeDetector"),
        ("stereoseq", False, "SharedMemoryStereoBarcodeDetector"),
        # visium_5prime gains splitting by reusing the 10x 3' detector
        ("visium_5prime", True, "TenXSplittingBarcodeDetector"),
        ("visium_5prime", False, "TenXBarcodeDetector"),
    ])
    def test_exact_detector_mapping(self, mode, split, expected):
        from isoquant_lib.barcode_calling.detect_barcodes import barcode_detector_class
        assert barcode_detector_class(IsoQuantMode[mode], split).__name__ == expected
