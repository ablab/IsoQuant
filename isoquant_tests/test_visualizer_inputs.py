############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Tests for locating and parsing IsoQuant output in the visualizer.

Covers the 4.0.0 output layout: the per-read file is SAMPLE.read_info.tsv[.gz]
(the old SAMPLE.read_assignments.tsv is still accepted), and grouped counts carry
the grouping strategy in their name.
"""

import gzip
import os
import pickle
from argparse import Namespace

import pytest

from isoquant_lib.scripts.convert_read_info import convert_to_read_assignments
from isoquant_lib.visualizer.post_process import DictionaryBuilder, OutputConfig

PREFIX = "S"

READ_INFO_HEADER = ("read_id\tchr\tstrand\tgene_id\tgene_assignment_type\tisoform_id\t"
                    "isoform_assignment_type\tassignment_events\tclassification\texons\t"
                    "polyA\tCAGE\tcanonical\tbarcode\tumi\tcell_type\tgroups\tadditional\n")

READ_INFO_ROWS = [
    ("r1\tchr1\t+\tGENE1\tunique\tTX1\tunique\texon_match\tfull_splice_match\t100-200\t"
     "True\t.\tTrue\tACGT\tUMI1\t.\tbc1\t*\n"),
    ("r2\tchr1\t+\tGENE1\tunique\tTX1\tunique_minor_difference\texon_elongation\t"
     "incomplete_splice_match\t100-190\tTrue\t.\tTrue\tACGT\tUMI2\t.\tbc1\t*\n"),
    ("r3\tchr1\t-\t.\tnoninformative\t.\tnoninformative\t.\t.\t500-600\tFalse\t.\t.\t"
     ".\t.\t.\tbc2\t*\n"),
]


def _write_read_info(path, gzipped=False):
    content = "# Command line: isoquant.py\n# IsoQuant version: 4.0.0\n" + READ_INFO_HEADER
    content += "".join(READ_INFO_ROWS)
    if gzipped:
        with gzip.open(path, "wt") as f:
            f.write(content)
    else:
        with open(path, "w") as f:
            f.write(content)
    return path


def _make_output(tmp_path, file_names=(), prefix=PREFIX):
    """Build a minimal IsoQuant output directory: <out>/.params + <out>/<prefix>/."""
    out_dir = tmp_path / "out"
    sample_dir = out_dir / prefix
    sample_dir.mkdir(parents=True)
    gtf = tmp_path / "ref.gtf"
    gtf.write_text("")
    params = Namespace(genedb=str(gtf), genedb_filename=None, fastq=None,
                       prefix=prefix, yaml=None)
    with open(str(out_dir / ".params"), "wb") as f:
        pickle.dump(params, f)
    for name in file_names:
        (sample_dir / name).write_text("")
    return out_dir, sample_dir


class TestReadFileDiscovery:
    def test_finds_read_info(self, tmp_path):
        out_dir, sample_dir = _make_output(tmp_path, [PREFIX + ".read_info.tsv"])
        config = OutputConfig(str(out_dir))
        assert config.read_assignments == str(sample_dir / (PREFIX + ".read_info.tsv"))

    def test_prefers_read_info_over_legacy(self, tmp_path):
        out_dir, sample_dir = _make_output(
            tmp_path, [PREFIX + ".read_info.tsv.gz", PREFIX + ".read_assignments.tsv"])
        config = OutputConfig(str(out_dir))
        assert config.read_assignments == str(sample_dir / (PREFIX + ".read_info.tsv.gz"))

    def test_falls_back_to_legacy(self, tmp_path):
        out_dir, sample_dir = _make_output(tmp_path, [PREFIX + ".read_assignments.tsv.gz"])
        config = OutputConfig(str(out_dir))
        assert config.read_assignments == str(sample_dir / (PREFIX + ".read_assignments.tsv.gz"))

    def test_gzipped_file_is_not_decompressed(self, tmp_path):
        out_dir, sample_dir = _make_output(tmp_path, [PREFIX + ".read_info.tsv.gz"])
        config = OutputConfig(str(out_dir))
        assert config.read_assignments.endswith(".gz")
        # A single-cell read_info is huge; unpacking it next to the original would
        # fill the user's output directory.
        assert not os.path.exists(str(sample_dir / (PREFIX + ".read_info.tsv")))

    def test_ignores_other_per_read_files(self, tmp_path):
        # Same format, different file - and the SQANTI-like output shares the
        # read_assignments stem.
        out_dir, _ = _make_output(tmp_path, [PREFIX + ".transcript_model_reads.tsv.gz",
                                             PREFIX + ".read_assignments.SQANTI-like.tsv"])
        config = OutputConfig(str(out_dir))
        assert config.read_assignments is None

    def test_missing_file_message_names_both_formats(self, tmp_path):
        out_dir, _ = _make_output(tmp_path)
        config = OutputConfig(str(out_dir))
        with pytest.raises(FileNotFoundError) as error:
            DictionaryBuilder(config).build_read_assignment_and_classification_dictionaries()
        assert "read_info" in str(error.value) and "read_assignments" in str(error.value)


class TestGroupedCountsDiscovery:
    def test_finds_strategy_named_files(self, tmp_path):
        out_dir, sample_dir = _make_output(tmp_path, [
            PREFIX + ".gene_grouped_barcode_counts.tsv",
            PREFIX + ".gene_grouped_barcode_tpm.tsv",
            PREFIX + ".transcript_grouped_barcode_counts.tsv",
            PREFIX + ".discovered_transcript_grouped_barcode_counts.tsv",
        ])
        config = OutputConfig(str(out_dir))
        assert config.group_strategies == ["barcode"]
        assert config.conditions
        assert config.gene_grouped_counts == str(sample_dir / (PREFIX + ".gene_grouped_barcode_counts.tsv"))
        assert config.transcript_grouped_counts is not None
        assert config.transcript_model_grouped_counts is not None

    def test_strategy_with_underscores(self, tmp_path):
        out_dir, _ = _make_output(tmp_path, [PREFIX + ".gene_grouped_file0_col1_counts.tsv"])
        config = OutputConfig(str(out_dir))
        assert config.group_strategies == ["file0_col1"]

    def test_discovered_files_are_not_read_as_plain_ones(self, tmp_path):
        out_dir, _ = _make_output(tmp_path, [
            PREFIX + ".discovered_transcript_grouped_barcode_counts.tsv"])
        config = OutputConfig(str(out_dir))
        assert config.transcript_grouped_counts is None
        assert config.transcript_model_grouped_counts is not None

    def test_selects_requested_strategy(self, tmp_path):
        out_dir, sample_dir = _make_output(tmp_path, [
            PREFIX + ".gene_grouped_barcode_counts.tsv",
            PREFIX + ".gene_grouped_file_name_counts.tsv",
        ])
        config = OutputConfig(str(out_dir), read_group_strategy="file_name")
        assert config.group_strategies == ["barcode", "file_name"]
        assert config.gene_grouped_counts == str(sample_dir / (PREFIX + ".gene_grouped_file_name_counts.tsv"))

    def test_unknown_strategy_is_rejected(self, tmp_path):
        out_dir, _ = _make_output(tmp_path, [PREFIX + ".gene_grouped_barcode_counts.tsv"])
        with pytest.raises(ValueError) as error:
            OutputConfig(str(out_dir), read_group_strategy="cell_type")
        assert "barcode" in str(error.value)

    def test_mtx_only_falls_back_to_ungrouped(self, tmp_path):
        out_dir, _ = _make_output(tmp_path, [
            PREFIX + ".gene_grouped_barcode_counts.matrix.mtx",
            PREFIX + ".gene_grouped_barcode_counts.barcodes.tsv",
            PREFIX + ".gene_grouped_barcode_counts.features.tsv",
            PREFIX + ".gene_tpm.tsv",
        ])
        config = OutputConfig(str(out_dir))
        assert config.group_strategies == ["barcode"]
        assert not config.conditions  # no matrix to plot, ungrouped files are used
        assert config.gene_grouped_counts is None

    def test_too_many_groups_are_not_pivoted(self, tmp_path):
        out_dir, sample_dir = _make_output(tmp_path)
        linear = sample_dir / (PREFIX + ".gene_grouped_barcode_counts.linear.tsv")
        rows = ["feature_id\tgroup_id\tcount\n"]
        rows += ["GENE1\tbc%d\t1.00\n" % i for i in range(150)]
        linear.write_text("".join(rows))
        config = OutputConfig(str(out_dir))
        assert config.group_strategies == ["barcode"]
        assert config.gene_grouped_counts is None
        assert not config.conditions

    def test_small_linear_file_is_converted(self, tmp_path):
        out_dir, sample_dir = _make_output(tmp_path)
        linear = sample_dir / (PREFIX + ".gene_grouped_barcode_counts.linear.tsv")
        linear.write_text("feature_id\tgroup_id\tcount\nGENE1\tbc1\t2.00\nGENE1\tbc2\t3.00\n")
        config = OutputConfig(str(out_dir))
        assert config.conditions
        assert config.gene_grouped_counts is not None
        assert os.path.exists(config.gene_grouped_counts)


class TestReadFileParsing:
    def _counts(self, config, path):
        config.read_assignments = path
        return DictionaryBuilder(config).build_read_assignment_and_classification_dictionaries()

    def test_read_info_counts(self, tmp_path):
        out_dir, sample_dir = _make_output(tmp_path)
        path = _write_read_info(str(sample_dir / (PREFIX + ".read_info.tsv")))
        config = OutputConfig(str(out_dir))
        classifications, assignment_types = self._counts(config, path)
        assert classifications == {"full_splice_match": 1, "incomplete_splice_match": 1, ".": 1}
        assert assignment_types == {"unique": 1, "unique_minor_difference": 1, "noninformative": 1}

    def test_gzipped_read_info_counts(self, tmp_path):
        out_dir, sample_dir = _make_output(tmp_path)
        plain = _write_read_info(str(sample_dir / (PREFIX + ".read_info.tsv")))
        gzipped = _write_read_info(str(sample_dir / (PREFIX + ".read_info.tsv.gz")), gzipped=True)
        config = OutputConfig(str(out_dir))
        assert self._counts(config, plain) == self._counts(config, gzipped)

    def test_legacy_format_gives_the_same_counts(self, tmp_path):
        # The shipped converter is the oracle: the same reads in the deprecated
        # layout must produce the same two dictionaries.
        out_dir, sample_dir = _make_output(tmp_path)
        read_info = _write_read_info(str(sample_dir / (PREFIX + ".read_info.tsv")))
        legacy = str(sample_dir / (PREFIX + ".read_assignments.tsv"))
        convert_to_read_assignments(read_info, legacy)
        config = OutputConfig(str(out_dir))
        assert self._counts(config, read_info) == self._counts(config, legacy)

    def test_legacy_classification_is_not_taken_from_groups(self, tmp_path):
        # additional_info is the 9th column, the last one is groups; reading the last
        # field returned the group id (or "NA") instead of the classification.
        out_dir, sample_dir = _make_output(tmp_path)
        legacy = sample_dir / (PREFIX + ".read_assignments.tsv")
        legacy.write_text(
            "# Command line: isoquant.py\n# IsoQuant version: 4.0.0\n"
            "read_id\tchr\tstrand\tisoform_id\tgene_id\tassignment_type\tassignment_events\t"
            "exons\tadditional_info\tgroups\n"
            "r1\tchr1\t+\tTX1\tGENE1\tunique\t.\t100-200\t"
            "PolyA=True; Classification=full_splice_match;\tNA\n")
        config = OutputConfig(str(out_dir))
        classifications, _ = self._counts(config, str(legacy))
        assert classifications == {"full_splice_match": 1}
