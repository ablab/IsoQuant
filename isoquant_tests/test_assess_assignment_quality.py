############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Tests for AssignmentData header-based column detection in
misc/assess_assignment_quality.py: the default --tool isoquant must parse both
the new read_info.tsv and the legacy read_assignments.tsv correctly."""

import importlib.util
import os
import tempfile

import pytest


def _load_qa():
    path = os.path.join(os.path.dirname(__file__), "..", "misc",
                        "assess_assignment_quality.py")
    spec = importlib.util.spec_from_file_location("assess_qa", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# read_info.tsv layout: isoform_id at col 5, isoform_assignment_type at col 6
READ_INFO_HEADER = ("read_id\tchr\tstrand\tgene_id\tgene_assignment_type\tisoform_id\t"
                    "isoform_assignment_type\tassignment_events\tclassification\texons\t"
                    "polyA\tCAGE\tcanonical\tbarcode\tumi\tcell_type\tgroups\tadditional\n")
READ_INFO_ROW = ("r1\tchr1\t+\tGENE1\tunique\tENST1\tunique\texon_match\t"
                 "full_splice_match\t100-200\tTrue\t.\tTrue\tACGT\tUMI1\t.\t.\t*\n")

# legacy read_assignments.tsv layout: isoform_id at col 3, assignment_type at col 5
READ_ASSIGN_HEADER = ("read_id\tchr\tstrand\tisoform_id\tgene_id\tassignment_type"
                      "\tassignment_events\texons\tadditional_info\tgroups\n")
READ_ASSIGN_ROW = "r1\tchr1\t+\tENST2\tGENE1\tunique\texon_match\t100-200\t*\tgrp\n"


def _write(tmp, content):
    path = os.path.join(tmp, "in.tsv")
    with open(path, "w") as f:
        f.write(content)
    return path


class TestHeaderDetection:
    def test_isoquant_preset_parses_read_info(self):
        qa = _load_qa()
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, READ_INFO_HEADER + READ_INFO_ROW)
            data = qa.AssignmentData(path, preset="isoquant")
        assert data.assigned_isoforms["r1"] == "ENST1"

    def test_isoquant_preset_parses_legacy_read_assignments(self):
        # The whole point of the fix: default --tool isoquant must still parse a
        # legacy read_assignments.tsv (detection overrides the read_info preset).
        qa = _load_qa()
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, READ_ASSIGN_HEADER + READ_ASSIGN_ROW)
            data = qa.AssignmentData(path, preset="isoquant")
        assert data.assigned_isoforms["r1"] == "ENST2"

    def test_legacy_preset_parses_read_info(self):
        qa = _load_qa()
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, READ_INFO_HEADER + READ_INFO_ROW)
            data = qa.AssignmentData(path, preset="isoquant_legacy")
        assert data.assigned_isoforms["r1"] == "ENST1"

    def test_headerless_falls_back_to_preset(self):
        # No header row -> keep the given preset (read_info column layout).
        qa = _load_qa()
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, READ_INFO_ROW)
            data = qa.AssignmentData(path, preset="isoquant")
        assert data.assigned_isoforms["r1"] == "ENST1"

    def test_detection_scoped_to_isoquant_presets(self):
        # talon/sqanti presets must never be overridden by a read_id header.
        qa = _load_qa()
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, READ_INFO_HEADER + READ_INFO_ROW)
            for preset in ("talon", "sqanti"):
                data = qa.AssignmentData.__new__(qa.AssignmentData)
                data.preset = preset
                tokens = READ_INFO_HEADER.strip().split("\t")
                assert data._detect_params_from_header(tokens) is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
