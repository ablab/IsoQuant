############################################################################
# Copyright (c) 2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Tests for ground truth parsing and split-mode scoring in misc/assess_barcode_quality.py."""

import importlib.util
import os

import pytest


def _load_qa():
    path = os.path.join(os.path.dirname(__file__), "..", "misc",
                        "assess_barcode_quality.py")
    spec = importlib.util.spec_from_file_location("assess_barcode_qa", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


QA = _load_qa()

# simulated Visium HD read names: per molecule BC1|BC2 (variable part lengths), UMI, strand; then NanoSim fields
VISIUM_BARCODES = ["GTCACAATATGCATG|TCACAATATGCATG", "GTCAACAATATGCATG|TCAACAATATGCATG",
                   "GTCACTATACTGCAAG|AAGATCATGTGCAAG", "GTTGAGCATCCCCAG|ATATGTCATGCACTG"]
NANOSIM_SUFFIX = "_0_aligned_17_R_0_3412_0"


def visium_read_name(n_molecules):
    name = "READ_5_ENSMUST00000207341.3"
    for i in range(n_molecules):
        name += "_%s_%s_%s" % (VISIUM_BARCODES[i], "AAATACCAC", "+-"[i % 2])
    return name + NANOSIM_SUFFIX


class TestGroundTruthSplit:
    @pytest.mark.parametrize("n_molecules", [1, 2, 3, 4])
    def test_visium_joint_barcodes(self, n_molecules):
        truth = QA.extract_ground_truth_split(visium_read_name(n_molecules), 31, joint_barcodes=True)
        assert truth == set(VISIUM_BARCODES[:n_molecules])

    def test_fixed_length_ignores_visium_barcodes(self):
        # without joint_barcodes only fields of exactly the barcode length count
        assert QA.extract_ground_truth_split(visium_read_name(2), 31) == set()

    def test_fixed_length_modes_unchanged(self):
        read_id = "READ_8_ENST00000367.4_CAACAGTAGGTTCCAT_GGACTCCGTCTT_-_TTTGGTTCAGGTCAGA_ACGTACGTACGT_+_0_aligned_1_F_0_900_0"
        assert QA.extract_ground_truth_split(read_id, 16) == {"CAACAGTAGGTTCCAT", "TTTGGTTCAGGTCAGA"}

    def test_mode_tables(self):
        assert "visium_hd_split" in QA.SPLIT_MODES
        assert "visium_hd_split" in QA.MODE_BARCODE_LENGTHS
        assert "visium_hd_split" in QA.JOINT_BARCODE_MODES


class TestAssessVisiumSplit:
    def test_split_suffixes_grouped_and_scored(self, tmp_path):
        read = visium_read_name(3)
        rows = [
            # two molecules found, named with the split segment suffix
            "%s_10_900_+\t%s\tAAATACCAC\t30" % (read, VISIUM_BARCODES[0]),
            "%s_15_870_-\t%s\tAAATACCAC\t29" % (read, VISIUM_BARCODES[1]),
            # a wrong barcode, and an undetected read
            "%s_20_500_+\t%s\tAAATACCAC\t28" % (read, VISIUM_BARCODES[3]),
            "%s\t*\t*\t-1" % visium_read_name(1),
        ]
        tsv = tmp_path / "barcodes.tsv"
        tsv.write_text("#read_id\tbarcode\tUMI\tBC_score\n" + "\n".join(rows) + "\n")
        metrics = QA.assess_split_barcode_mode(str(tsv), 31, joint_barcodes=True)
        assert metrics["total_reads"] == 2
        assert metrics["total_true_barcodes"] == 4
        assert metrics["correct_assignments"] == 2
        assert metrics["incorrect_assignments"] == 1
        assert metrics["no_barcodes_assigned"] == 1


# barcodes 0 and 1 share a spot, 2 and 3 have their own
SPOTS = {VISIUM_BARCODES[0]: "s_008um_00000_00000", VISIUM_BARCODES[1]: "s_008um_00000_00000",
         VISIUM_BARCODES[2]: "s_008um_00001_00000", VISIUM_BARCODES[3]: "s_008um_00002_00000"}


class TestSpotLevelScoring:
    def test_load_spot_map_column(self, tmp_path):
        mapping = tmp_path / "map.tsv"
        mapping.write_text("#barcode\t2um\t8um\n" + "".join("%s\tfine\t%s\n" % kv for kv in SPOTS.items()))
        assert QA.load_spot_map(str(mapping) + ":0:2") == SPOTS

    def test_to_spot(self):
        assert QA.to_spot(VISIUM_BARCODES[1], SPOTS) == "s_008um_00000_00000"
        assert QA.to_spot("s_008um_00000_00000", SPOTS) == "s_008um_00000_00000"  # caller reported spots
        assert QA.to_spot(VISIUM_BARCODES[1], None) == VISIUM_BARCODES[1]

    def test_single_mode_counts_neighbouring_barcode_as_correct(self, tmp_path):
        read = visium_read_name(1)  # true barcode VISIUM_BARCODES[0]
        tsv = tmp_path / "barcodes.tsv"
        tsv.write_text("#read_id\tbarcode\tUMI\tBC_score\n%s\t%s\tAAATACCAC\t30\n" % (read, VISIUM_BARCODES[1]))
        assert QA.assess_single_barcode_mode(str(tsv), 31)["correct_barcodes"] == 0
        assert QA.assess_single_barcode_mode(str(tsv), 31, spot_map=SPOTS)["correct_barcodes"] == 1

    def test_split_mode_with_spot_ids(self, tmp_path):
        read = visium_read_name(3)
        rows = ["%s_10_900_+\ts_008um_00000_00000\tAAATACCAC\t30" % read,
                "%s_15_870_-\t%s\tAAATACCAC\t29" % (read, VISIUM_BARCODES[2])]
        tsv = tmp_path / "barcodes.tsv"
        tsv.write_text("#read_id\tbarcode\tUMI\tBC_score\n" + "\n".join(rows) + "\n")
        metrics = QA.assess_split_barcode_mode(str(tsv), 31, joint_barcodes=True, spot_map=SPOTS)
        # truth: barcodes 0 and 1 collapse into one spot, plus barcode 2
        assert metrics["total_true_barcodes"] == 2
        assert metrics["correct_assignments"] == 2
