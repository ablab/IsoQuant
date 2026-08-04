############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import random

import pytest

from isoquant_lib.barcode_calling.correct_barcodes import (
    NOSEQ,
    correct_barcodes,
    count_barcodes,
    write_corrected_table,
)
from isoquant_lib.barcode_calling.barcode_graph import BarcodeGraph

HEADER = "read_id\tbarcode\tUMI\tBC_score\tvalid_UMI\tstrand\tpolyT_start\tR1_end"


def random_barcode(rnd, length=16):
    return "".join(rnd.choice("ACGT") for _ in range(length))


def write_raw_table(path, rows):
    with open(path, "w") as handle:
        handle.write(HEADER + "\n")
        for read_id, barcode in rows:
            handle.write("%s\t%s\tACGTACGTACGT\t0\tTrue\t+\t50\t30\n" % (read_id, barcode))


def read_table(path):
    with open(path) as handle:
        lines = [line.rstrip("\n").split("\t") for line in handle]
    return lines[0], lines[1:]


class TestCounting:
    def test_header_is_not_counted_as_a_read(self, tmp_path):
        """The raw table header is not comment-prefixed, so it must be detected by content."""
        raw = tmp_path / "raw.tsv"
        write_raw_table(raw, [("r1", "ACGTACGTACGTACGT"), ("r2", "ACGTACGTACGTACGT")])
        graph = BarcodeGraph(barcode_length=16)
        rows = count_barcodes([str(raw)], graph)
        assert rows == 2
        assert graph.malformed_count == 0
        assert dict(graph.counts) == {"ACGTACGTACGTACGT": 2}

    def test_missing_barcodes_are_skipped(self, tmp_path):
        raw = tmp_path / "raw.tsv"
        write_raw_table(raw, [("r1", "ACGTACGTACGTACGT"), ("r2", NOSEQ)])
        graph = BarcodeGraph(barcode_length=16)
        assert count_barcodes([str(raw)], graph) == 2
        assert dict(graph.counts) == {"ACGTACGTACGTACGT": 1}
        assert graph.malformed_count == 0

    def test_counts_accumulate_across_files(self, tmp_path):
        first, second = tmp_path / "a.tsv", tmp_path / "b.tsv"
        write_raw_table(first, [("r1", "ACGTACGTACGTACGT")])
        write_raw_table(second, [("r2", "ACGTACGTACGTACGT"), ("r3", "TTTTGGGGCCCCAAAA")])
        graph = BarcodeGraph(barcode_length=16)
        assert count_barcodes([str(first), str(second)], graph) == 3
        assert graph.counts["ACGTACGTACGTACGT"] == 2


class TestRewrite:
    def test_header_gains_the_raw_column(self, tmp_path):
        raw, out = tmp_path / "raw.tsv", tmp_path / "out.tsv"
        write_raw_table(raw, [("r1", "ACGTACGTACGTACGT")])
        write_corrected_table(str(raw), str(out), {"ACGTACGTACGTACGT": "ACGTACGTACGTACGT"}, 16)
        header, rows = read_table(out)
        assert header == HEADER.split("\t") + ["raw_barcode"]
        assert len(rows) == 1

    def test_corrected_barcode_replaces_the_raw_one(self, tmp_path):
        raw, out = tmp_path / "raw.tsv", tmp_path / "out.tsv"
        write_raw_table(raw, [("r1", "ACGTACGTACGTACGA")])
        stats = write_corrected_table(str(raw), str(out),
                                      {"ACGTACGTACGTACGA": "ACGTACGTACGTACGT"}, 16)
        _, rows = read_table(out)
        assert rows[0][1] == "ACGTACGTACGTACGT"
        assert rows[0][-1] == "ACGTACGTACGTACGA"
        assert stats["corrected"] == 1

    def test_unassigned_barcode_becomes_noseq(self, tmp_path):
        raw, out = tmp_path / "raw.tsv", tmp_path / "out.tsv"
        write_raw_table(raw, [("r1", "TTTTGGGGCCCCAAAA")])
        stats = write_corrected_table(str(raw), str(out), {}, 16)
        _, rows = read_table(out)
        assert rows[0][1] == NOSEQ
        assert rows[0][-1] == "TTTTGGGGCCCCAAAA"  # the raw window is preserved
        assert stats["dropped"] == 1

    def test_reads_without_a_barcode_pass_through(self, tmp_path):
        raw, out = tmp_path / "raw.tsv", tmp_path / "out.tsv"
        write_raw_table(raw, [("r1", NOSEQ)])
        stats = write_corrected_table(str(raw), str(out), {}, 16)
        _, rows = read_table(out)
        assert rows[0][1] == NOSEQ
        assert stats["no_barcode"] == 1

    def test_overlong_window_is_looked_up_trimmed(self, tmp_path):
        raw, out = tmp_path / "raw.tsv", tmp_path / "out.tsv"
        write_raw_table(raw, [("r1", "ACGTACGTACGTACGTA")])
        write_corrected_table(str(raw), str(out), {"ACGTACGTACGTACGT": "ACGTACGTACGTACGT"}, 16)
        _, rows = read_table(out)
        assert rows[0][1] == "ACGTACGTACGTACGT"

    def test_column_order_is_preserved(self, tmp_path):
        """Everything downstream reads read_id, barcode, UMI from columns 0, 1, 2."""
        raw, out = tmp_path / "raw.tsv", tmp_path / "out.tsv"
        write_raw_table(raw, [("r1", "ACGTACGTACGTACGT")])
        write_corrected_table(str(raw), str(out), {"ACGTACGTACGTACGT": "ACGTACGTACGTACGT"}, 16)
        _, rows = read_table(out)
        assert rows[0][0] == "r1"
        assert rows[0][2] == "ACGTACGTACGT"


class TestEndToEnd:
    def _build_sample(self, tmp_path, rnd, cell_count=30, reads_per_cell=40):
        cells = [random_barcode(rnd) for _ in range(cell_count)]
        whitelist = tmp_path / "whitelist.txt"
        with open(whitelist, "w") as handle:
            for barcode in sorted(set(cells) | {random_barcode(rnd) for _ in range(2000)}):
                handle.write(barcode + "\n")

        rows, truth = [], {}
        for index, barcode in enumerate(cells):
            for i in range(reads_per_cell):
                read_id = "r%d_%d" % (index, i)
                observed = barcode
                if i % 4 == 3:  # a quarter of the reads carry one substitution
                    position = rnd.randrange(16)
                    observed = (barcode[:position] +
                                rnd.choice([c for c in "ACGT" if c != barcode[position]]) +
                                barcode[position + 1:])
                rows.append((read_id, observed))
                truth[read_id] = barcode
        rnd.shuffle(rows)

        raw = tmp_path / "raw_0.tsv"
        write_raw_table(raw, rows)
        return str(raw), str(whitelist), truth

    @pytest.mark.parametrize("implementation", ["centers", "full"])
    def test_recovers_error_bearing_barcodes(self, tmp_path, implementation):
        rnd = random.Random(41)
        raw, whitelist, truth = self._build_sample(tmp_path, rnd)
        out = str(tmp_path / ("corrected_%s.tsv" % implementation))

        correct_barcodes([raw], [out], [whitelist], n_cells=30, implementation=implementation)

        _, rows = read_table(out)
        assert len(rows) == len(truth)
        assert all(row[1] == truth[row[0]] for row in rows)

    def test_both_implementations_produce_the_same_file(self, tmp_path):
        rnd = random.Random(42)
        raw, whitelist, _ = self._build_sample(tmp_path, rnd)
        outputs = []
        for implementation in ("centers", "full"):
            out = str(tmp_path / ("out_%s.tsv" % implementation))
            correct_barcodes([raw], [out], [whitelist], n_cells=30, implementation=implementation)
            outputs.append(open(out).read())
        assert outputs[0] == outputs[1]

    def test_stats_file_is_written(self, tmp_path):
        rnd = random.Random(43)
        raw, whitelist, _ = self._build_sample(tmp_path, rnd, cell_count=5, reads_per_cell=20)
        stats_file = tmp_path / "stats.tsv"
        stats = correct_barcodes([raw], [str(tmp_path / "out.tsv")], [whitelist],
                                 n_cells=5, stats_file=str(stats_file))
        written = dict(line.split("\t") for line in stats_file.read_text().splitlines())
        assert int(written["Cluster centers"]) == stats["Cluster centers"] == 5

    def test_cell_number_is_estimated_when_not_given(self, tmp_path):
        rnd = random.Random(44)
        raw, whitelist, truth = self._build_sample(tmp_path, rnd, cell_count=50)
        out = str(tmp_path / "out.tsv")
        stats = correct_barcodes([raw], [out], [whitelist])
        assert stats["Cluster centers"] == 50

    def test_no_whitelist_still_works(self, tmp_path):
        rnd = random.Random(45)
        raw, _, truth = self._build_sample(tmp_path, rnd, cell_count=10)
        out = str(tmp_path / "out.tsv")
        correct_barcodes([raw], [out], None, n_cells=10)
        _, rows = read_table(out)
        assert all(row[1] == truth[row[0]] for row in rows)

    def test_empty_input(self, tmp_path):
        raw, out = tmp_path / "raw.tsv", tmp_path / "out.tsv"
        write_raw_table(raw, [])
        stats = correct_barcodes([str(raw)], [str(out)], None, n_cells=10)
        assert stats["Distinct barcodes observed"] == 0
        header, rows = read_table(out)
        assert header[-1] == "raw_barcode"
        assert rows == []
