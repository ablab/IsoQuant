############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Tests for the per-experiment run summary (SAMPLE.summary_stats.json / .summary.html)."""

import json
import os

from isoquant_lib.report.html_report import render_html
from isoquant_lib.report.run_summary import RunSummary, read_stats_tsv


class FakeEnum:
    """Stand-in for the enum members EnumStats counts by."""

    def __init__(self, name):
        self.name = name

    def __lt__(self, other):
        return self.name < other.name


def _enum_stats(**counts):
    return {FakeEnum(name): count for name, count in counts.items()}


class FakeSample:
    def __init__(self, directory, prefix="S"):
        self.prefix = prefix
        base = os.path.join(directory, prefix)
        self.out_gene_counts_tsv = base + ".gene"
        self.out_transcript_counts_tsv = base + ".transcript"
        self.out_gene_grouped_counts_tsv = base + ".gene_grouped"
        self.out_transcript_grouped_counts_tsv = base + ".transcript_grouped"
        self.barcodes_tsv = base + ".barcoded_reads"
        self.out_cell_barcodes_stats = base + ".cell_barcodes.stats"
        self.out_umi_filtered = base + ".UMI_filtered"


def _write(path, content):
    with open(path, "w") as f:
        f.write(content)


def _populated_sample(tmp_path):
    """Write the stat files a single-cell run leaves behind."""
    directory = str(tmp_path)
    sample = FakeSample(directory)
    _write(sample.out_gene_counts_tsv + "_counts.tsv",
           "feature_id\tcount\nGENE1\t60.00\nGENE2\t40.00\n"
           "__ambiguous\t5\n__no_feature\t3\n__not_aligned\t10\n")
    _write(sample.out_transcript_counts_tsv + "_counts.tsv",
           "feature_id\tcount\nTX1\t30.00\nTX2\t20.00\n"
           "__ambiguous\t9\n__no_feature\t4\n__not_aligned\t10\n")
    _write(sample.barcodes_tsv + "_0.tsv.stats",
           "Total reads\t100\nBarcode detected\t80\nPolyT detected\t95\n")
    _write(sample.barcodes_tsv + "_1.tsv.stats",
           "Total reads\t100\nBarcode detected\t90\nPolyT detected\t99\n")
    _write(sample.out_cell_barcodes_stats, "Cell barcodes detected\t3\n")
    _write(sample.out_umi_filtered + ".ED4.stats.tsv",
           "Total reads saved\t60\nTotal assignments processed\t150\n"
           "Assigned to any gene and barcoded\t120\n")
    _write(sample.out_gene_grouped_counts_tsv + "_barcode_counts.linear.tsv",
           "feature_id\tgroup_id\tcount\n"
           "GENE1\tbc1\t10.00\nGENE2\tbc1\t5.00\nGENE1\tbc2\t4.00\nGENE1\tbc3\t1.00\n")
    _write(sample.out_transcript_grouped_counts_tsv + "_barcode_counts.linear.tsv",
           "feature_id\tgroup_id\tcount\nTX1\tbc1\t8.00\nTX1\tbc2\t2.00\n")
    return sample


def _summary(tmp_path):
    summary = RunSummary("S", isoquant_version="4.0.0", command_line="isoquant.py -o out",
                         mode="tenX_v3")
    summary.set_alignment_stats(_enum_stats(primary=190, secondary=50, supplementary=20,
                                            unaligned=10))
    summary.set_assignment_stats(_enum_stats(unique=100, unique_minor_difference=20,
                                             ambiguous=30, inconsistent=25,
                                             inconsistent_non_intronic=10, noninformative=5,
                                             intergenic=5))
    summary.set_transcript_model_stats(_enum_stats(known=7, novel_in_catalog=2))
    summary.set_polya_stats(195, 156)
    summary.collect_output_files(_populated_sample(tmp_path), ["barcode"])
    return summary


class TestReadStatsTsv:
    def test_missing_file_is_empty(self):
        assert read_stats_tsv("/no/such/file.stats") == {}
        assert read_stats_tsv(None) == {}

    def test_skips_malformed_lines(self, tmp_path):
        path = str(tmp_path / "x.stats")
        _write(path, "Total reads\t10\nnot a stat line\nBad count\tNaN\nBarcode detected\t7\n")
        assert read_stats_tsv(path) == {"Total reads": 10, "Barcode detected": 7}


class TestRates:
    def test_input_reads_exclude_extra_alignment_records(self, tmp_path):
        summary = _summary(tmp_path)
        # secondary and supplementary records belong to reads already counted as primary
        assert summary.total_reads == 200
        assert summary.mapping_rate == 0.95

    def test_barcode_rate_sums_all_input_files(self, tmp_path):
        summary = _summary(tmp_path)
        assert summary.barcodes["Total reads"] == 200
        assert summary.barcodes["Barcode detected"] == 170
        assert summary.barcode_rate == 0.85

    def test_assignment_rollup(self, tmp_path):
        rollup = _summary(tmp_path).assignment_rollup()
        assert rollup["total"] == 195
        assert rollup["unique"] == 120
        assert rollup["ambiguous"] == 30
        assert rollup["inconsistent"] == 35  # inconsistent + inconsistent_non_intronic
        assert rollup["unassigned"] == 10  # noninformative + intergenic

    def test_counted_reads_are_not_turned_into_a_rate_of_their_own(self, tmp_path):
        # Reads dropped by the counting strategy are in neither the features nor the
        # __ambiguous / __no_feature tail, so the share is taken over input reads.
        summary = _summary(tmp_path)
        gene = summary.to_dict()["quantification"]["gene"]
        assert gene["counted"] == 100.0
        assert gene["share_of_reads"] == 0.5
        assert "assignment_rate" not in gene

    def test_group_rollup(self, tmp_path):
        groups = _summary(tmp_path).group_rollup("barcode")
        assert groups["gene"]["groups"] == 3
        assert groups["gene"]["reads"] == 20.0
        assert groups["gene"]["share_of_reads"] == 0.1
        assert groups["gene"]["median_reads_per_group"] == 4.0
        assert groups["gene"]["median_features_per_group"] == 1.0
        assert groups["transcript"]["groups"] == 2

    def test_umi_duplication_rate(self, tmp_path):
        umi = _summary(tmp_path).to_dict()["umi_filtering"]
        assert umi["edit_distance"] == 4
        assert umi["molecules"] == 60
        assert umi["duplication_rate"] == 0.6


class TestJson:
    def test_written_json_round_trips(self, tmp_path):
        summary = _summary(tmp_path)
        path = str(tmp_path / "S.summary_stats.json")
        summary.write_json(path)
        with open(path) as f:
            data = json.load(f)
        assert data["sample"] == "S"
        assert data["mode"] == "tenX_v3"
        assert set(data) >= {"alignment", "assignment", "quantification", "transcript_models",
                             "barcodes", "umi_filtering", "groups"}
        assert data["assignment"]["polya_rate"] == 0.8
        assert data["transcript_models"] == {"known": 7, "novel_in_catalog": 2}

    def test_per_group_depths_stay_out_of_the_json(self, tmp_path):
        # The rank curve behind the figure would be one entry per barcode.
        summary = _summary(tmp_path)
        path = str(tmp_path / "S.summary_stats.json")
        summary.write_json(path)
        with open(path) as f:
            assert "ranked_reads" not in f.read()

    def test_empty_summary_is_still_valid(self, tmp_path):
        summary = RunSummary("S")
        path = str(tmp_path / "empty.json")
        summary.write_json(path)
        with open(path) as f:
            assert json.load(f)["sample"] == "S"


class TestHtml:
    def _render(self, tmp_path, summary=None):
        path = str(tmp_path / "S.summary.html")
        render_html(summary or _summary(tmp_path), path)
        with open(path) as f:
            return f.read()

    def test_page_is_self_contained(self, tmp_path):
        page = self._render(tmp_path)
        assert page.startswith("<!DOCTYPE html>")
        # no external stylesheets, scripts or images: figures are inlined as SVG
        assert "src=" not in page
        assert "<link" not in page and "<script" not in page

    def test_headline_metrics_are_shown(self, tmp_path):
        page = self._render(tmp_path)
        assert "95.0%" in page  # mapping rate
        assert "85.0%" in page  # valid barcode rate
        assert "Valid barcodes" in page
        assert "Barcode calling" in page
        assert "UMI deduplication (edit distance 4)" in page

    def test_sections_are_skipped_when_empty(self, tmp_path):
        summary = RunSummary("S", isoquant_version="4.0.0")
        summary.set_alignment_stats(_enum_stats(primary=10, unaligned=0))
        page = self._render(tmp_path, summary)
        assert "Alignment" in page
        assert "Barcode calling" not in page
        assert "Grouped counts" not in page
