############################################################################
# Copyright (c) 2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Visium HD barcode calling with --barcode2barcode: barcodes are reported as spot IDs, and a part whose
best candidates tie is still resolved when all tied barcode pairs fall into one spot."""

import argparse
import gzip
import pickle

import pytest

import isoquant
from isoquant_lib.barcode_calling.callers import (
    VisiumHDBarcodeDetector,
    VisiumHDSplittingBarcodeDetector,
)
from isoquant_lib.barcode_calling.common import find_tied_candidates_ssw_var_len, reverese_complement
from isoquant_lib.barcode_calling.options import resolve_spot_barcode_calling, validate_barcode_calling
from isoquant_lib.barcode_calling.spot_mapping import SpotMapping, parse_spot_mapping_spec
from isoquant_lib.modes import IsoQuantMode

# BC1 candidates 0 and 1 differ only in the last base, so a read ending that part with a third base ties
PART1 = ["ACGTACGTACGTACGA", "ACGTACGTACGTACGC", "TTGGCCAATTGGCCAA"]
PART2 = ["CATGCATGCATGCAT", "GGATCCGGATCCGGA"]
AMBIGUOUS_BC1 = "ACGTACGTACGTACGG"

# 2 um barcodes -> coarser spots: PART1[0] and PART1[1] share a spot with PART2[0] but not with PART2[1]
MAPPING_ROWS = [
    (PART1[0] + "|" + PART2[0], "s_002um_00000_00000", "s_008um_00000_00000"),
    (PART1[1] + "|" + PART2[0], "s_002um_00001_00000", "s_008um_00000_00000"),
    (PART1[2] + "|" + PART2[0], "s_002um_00009_00000", "s_008um_00002_00000"),
    (PART1[0] + "|" + PART2[1], "s_002um_00000_00009", "s_008um_00000_00002"),
    (PART1[1] + "|" + PART2[1], "s_002um_00001_00009", "s_008um_00001_00002"),
]

R1 = "CTACACGACGCTCTTCCGATCT"
TSO_TAIL = "ATGTACTCTGCGTTGATACCACTGCTT"
CDNA = "GCTAGCATCGATCGGATCCAGTACGTAGCTAGGCTAACGTAGGCATGCAAGTCGACGTAGCATGCAGTCAGTACGATCGATGCATGCAGT" * 3


def write_mapping(path, rows=MAPPING_ROWS, gzipped=False):
    text = "#barcode\tspot_2um\tspot_8um\n" + "".join("\t".join(row) + "\n" for row in rows)
    if gzipped:
        with gzip.open(str(path), "wt") as f:
            f.write(text)
    else:
        path.write_text(text)
    return str(path)


def molecule(bc1, bc2, umi="ACGTTGCAA"):
    return R1 + umi + bc1 + bc2 + "GC" + "T" * 30 + CDNA + TSO_TAIL


@pytest.fixture
def mapping_file(tmp_path):
    return write_mapping(tmp_path / "barcode2spot.tsv")


@pytest.fixture
def spot_mapping(mapping_file):
    return SpotMapping(PART1, PART2, mapping_file, barcode_col=0, spot_col=2)


class TestSpotMapping:
    def test_spot(self, spot_mapping):
        assert spot_mapping.spot(PART1[0], PART2[0]) == "s_008um_00000_00000"
        assert spot_mapping.spot(PART1[2], PART2[1]) is None  # pair not in the mapping
        assert spot_mapping.spot("AAAA", PART2[0]) is None  # not a whitelist barcode

    def test_unique_spot(self, spot_mapping):
        assert spot_mapping.unique_spot([PART1[0], PART1[1]], [PART2[0]]) == "s_008um_00000_00000"
        assert spot_mapping.unique_spot([PART1[0], PART1[1]], [PART2[1]]) is None  # two spots
        # unmapped pairs do not compete: PART1[2]|PART2[1] is not in the mapping
        assert spot_mapping.unique_spot([PART1[2]], [PART2[0], PART2[1]]) == "s_008um_00002_00000"
        assert spot_mapping.unique_spot([], [PART2[0]]) is None

    def test_column_selection(self, mapping_file):
        fine = SpotMapping(PART1, PART2, mapping_file, barcode_col=0, spot_col=1)
        assert fine.spot(PART1[1], PART2[0]) == "s_002um_00001_00000"
        assert fine.unique_spot([PART1[0], PART1[1]], [PART2[0]]) is None

    def test_gzipped(self, tmp_path):
        mapping = SpotMapping(PART1, PART2, write_mapping(tmp_path / "m.tsv.gz", gzipped=True), 0, 2)
        assert mapping.spot(PART1[0], PART2[0]) == "s_008um_00000_00000"

    def test_skips_foreign_barcodes(self, tmp_path):
        rows = MAPPING_ROWS + [("GGGG|" + PART2[0], "x", "s_008um_09999_09999"), ("NOSEPARATOR", "x", "y")]
        mapping = SpotMapping(PART1, PART2, write_mapping(tmp_path / "m.tsv", rows), 0, 2)
        assert "s_008um_09999_09999" not in mapping.spot_names

    def test_no_matching_barcode_exits(self, tmp_path):
        with pytest.raises(SystemExit):
            SpotMapping(["AAAA"], ["CCCC"], write_mapping(tmp_path / "m.tsv"), 0, 2)

    def test_pickle(self, spot_mapping):
        restored = pickle.loads(pickle.dumps(spot_mapping))
        assert restored.unique_spot([PART1[0], PART1[1]], [PART2[0]]) == "s_008um_00000_00000"


class TestParseSpotMappingSpec:
    def test_default_column(self):
        assert parse_spot_mapping_spec("map.tsv") == ("map.tsv", 0, 1)

    def test_selected_column(self):
        assert parse_spot_mapping_spec("map.tsv:0:2") == ("map.tsv", 0, 2)

    def test_several_columns_exit(self):
        with pytest.raises(SystemExit):
            parse_spot_mapping_spec("map.tsv:0:1,2,3")


class TestTiedCandidates:
    def test_ties(self):
        matches = [(PART1[0], 0, []), (PART1[1], 0, []), (PART1[2], 0, [])]
        score, ties = find_tied_candidates_ssw_var_len(matches, AMBIGUOUS_BC1, min_score=13)
        assert score == 15
        assert sorted(t[0] for t in ties) == sorted(PART1[:2])

    def test_nothing_good_enough(self):
        assert find_tied_candidates_ssw_var_len([(PART1[2], 0, [])], AMBIGUOUS_BC1, min_score=13) == (0, [])


class TestDetectorSpots:
    @pytest.mark.parametrize("detector_class", [VisiumHDBarcodeDetector, VisiumHDSplittingBarcodeDetector])
    def test_tie_within_one_spot_is_assigned(self, detector_class, spot_mapping):
        detector = detector_class([PART1, PART2], spot_mapping=spot_mapping)
        result = detector.find_barcode_umi("read", molecule(AMBIGUOUS_BC1, PART2[0]))
        assert result.get_barcode() == "s_008um_00000_00000"
        assert result.get_umi() == "ACGTTGCAA"

    @pytest.mark.parametrize("detector_class", [VisiumHDBarcodeDetector, VisiumHDSplittingBarcodeDetector])
    def test_tie_across_spots_is_not_assigned(self, detector_class, spot_mapping):
        detector = detector_class([PART1, PART2], spot_mapping=spot_mapping)
        result = detector.find_barcode_umi("read", molecule(AMBIGUOUS_BC1, PART2[1]))
        assert not result.is_valid()

    def test_unique_barcode_reported_as_spot(self, spot_mapping):
        detector = VisiumHDBarcodeDetector([PART1, PART2], spot_mapping=spot_mapping)
        result = detector.find_barcode_umi("read", molecule(PART1[2], PART2[0]))
        assert result.get_barcode() == "s_008um_00002_00000"

    def test_without_mapping_tie_stays_unassigned(self):
        detector = VisiumHDBarcodeDetector([PART1, PART2])
        assert not detector.find_barcode_umi("read", molecule(AMBIGUOUS_BC1, PART2[0])).is_valid()
        result = detector.find_barcode_umi("read", molecule(PART1[2], PART2[0]))
        assert result.get_barcode() == PART1[2] + "|" + PART2[0]

    def test_window_margin_tolerates_shifted_polyt(self, tmp_path):
        """A longer separator before polyT shifts the BC2 window by 3; spot mode searches a wider window."""
        part1, part2 = ["GATCCTAGCTTGACAG"], ["CAGTTCAAGGCTAGC"]
        mapping = SpotMapping(part1, part2, write_mapping(
            tmp_path / "m.tsv", [(part1[0] + "|" + part2[0], "s_002um_00028_00028", "s_008um_00007_00007")]), 0, 2)
        read = R1 + "ACGTTGCAA" + part1[0] + part2[0] + "GCAGC" + "T" * 30 + CDNA + TSO_TAIL
        assert not VisiumHDBarcodeDetector([part1, part2]).find_barcode_umi("read", read).is_valid()
        result = VisiumHDBarcodeDetector([part1, part2], spot_mapping=mapping).find_barcode_umi("read", read)
        assert result.get_barcode() == "s_008um_00007_00007"
        assert result.get_umi() == "ACGTTGCAA"

    def test_no_polyt_path(self, spot_mapping):
        """R1-anchored detection without polyT resolves ties the same way."""
        detector = VisiumHDBarcodeDetector([PART1, PART2], spot_mapping=spot_mapping)
        read = R1 + "ACGTTGCAA" + AMBIGUOUS_BC1 + PART2[0] + "GCAGCAGCAGCAGCAG"
        assert detector.find_barcode_umi("read", read).get_barcode() == "s_008um_00000_00000"

    def test_concatenated_read(self, spot_mapping):
        detector = VisiumHDSplittingBarcodeDetector([PART1, PART2], spot_mapping=spot_mapping)
        read = molecule(AMBIGUOUS_BC1, PART2[0]) + reverese_complement(molecule(PART1[2], PART2[0]))
        result = detector.find_barcode_umi("read", read)
        assert sorted((p.get_barcode(), p.strand) for p in result.detected_patterns) == \
            [("s_008um_00000_00000", "+"), ("s_008um_00002_00000", "-")]

    def test_pickle(self, spot_mapping):
        detector = pickle.loads(pickle.dumps(VisiumHDSplittingBarcodeDetector([PART1, PART2], spot_mapping)))
        result = detector.find_barcode_umi("read", molecule(AMBIGUOUS_BC1, PART2[0]))
        assert result.get_barcode() == "s_008um_00000_00000"


def calling_args(**kwargs):
    defaults = dict(mode=IsoQuantMode.visium_hd, barcode_whitelist=["p1.tsv", "p2.tsv"],
                    barcode2barcode="map.tsv:0:2", read_group=None)
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


class TestResolveSpotBarcodeCalling:
    def test_visium_hd_with_whitelist(self):
        args = calling_args()
        resolve_spot_barcode_calling(args)
        assert args.spot_barcode_calling

    @pytest.mark.parametrize("kwargs", [
        dict(barcode2barcode=None),
        dict(barcode_whitelist=None),  # barcodes from --barcoded_reads / --barcoded_bam
        dict(mode=IsoQuantMode.tenX_v3),
    ])
    def test_not_applied(self, kwargs):
        args = calling_args(**kwargs)
        resolve_spot_barcode_calling(args)
        assert not args.spot_barcode_calling

    def test_several_spot_columns_exit(self):
        with pytest.raises(SystemExit):
            resolve_spot_barcode_calling(calling_args(barcode2barcode="map.tsv:0:1,2,3"))

    def test_barcode_barcode_group_exits(self):
        with pytest.raises(SystemExit):
            resolve_spot_barcode_calling(calling_args(read_group=["barcode_barcode"]))

    def test_bulk_resets_flag(self):
        args = argparse.Namespace(mode=IsoQuantMode.bulk)
        validate_barcode_calling(args)
        assert args.spot_barcode_calling is False


class _InputData:
    @staticmethod
    def has_replicas():
        return False


def data_dependent_args(spot_barcode_calling):
    return argparse.Namespace(
        matching_strategy="default", data_type=isoquant.NANOPORE_DATA, model_construction_strategy="default_ont",
        fl_data=False, splice_correction_strategy="default_ont", read_group=None, input_data=_InputData(),
        barcode2spot=None, barcode2barcode="map.tsv:0:2", large_output=None, mode=IsoQuantMode.visium_hd,
        spot_barcode_calling=spot_barcode_calling)


class TestReadGroupsWithSpotCalling:
    def test_spot_calling_groups_by_barcode(self):
        args = data_dependent_args(True)
        isoquant.set_data_dependent_options(args)
        assert args.read_group == ["barcode"]
        assert "allinfo" in args.large_output

    def test_barcode2barcode_at_umi_stage_keeps_its_group(self):
        args = data_dependent_args(False)
        isoquant.set_data_dependent_options(args)
        assert args.read_group == ["barcode_barcode"]
