############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import random

import pytest

from isoquant_lib.barcode_calling.barcode_graph import BarcodeGraph
from isoquant_lib.barcode_calling.correct_barcodes import (
    NOSEQ,
    count_barcodes,
    load_whitelist,
    select_cell_barcodes,
)

HEADER = "read_id\tbarcode\tUMI\tBC_score\tvalid_UMI\tstrand\tpolyT_start\tR1_end"


def random_barcode(rnd, length=16):
    return "".join(rnd.choice("ACGT") for _ in range(length))


def write_raw_table(path, rows):
    with open(path, "w") as handle:
        handle.write(HEADER + "\n")
        for read_id, barcode in rows:
            handle.write("%s\t%s\tACGTACGTACGT\t0\tTrue\t+\t50\t30\n" % (read_id, barcode))


def read_barcodes(path):
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


class TestCounting:
    def test_header_is_not_counted_as_a_read(self, tmp_path):
        """Raw tables carry a plain, non-commented header, so it must be spotted by content."""
        raw = tmp_path / "raw.tsv"
        write_raw_table(raw, [("r1", "ACGTACGTACGTACGT"), ("r2", "ACGTACGTACGTACGT")])
        graph = BarcodeGraph(barcode_length=16)
        assert count_barcodes([str(raw)], graph) == 2
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


class TestLoadWhitelist:
    def test_auto_means_no_pool(self):
        assert load_whitelist(["auto"]) is None

    def test_missing_means_no_pool(self):
        assert load_whitelist(None) is None
        assert load_whitelist([]) is None

    def test_files_are_unioned(self, tmp_path):
        first, second = tmp_path / "a.txt", tmp_path / "b.txt"
        first.write_text("AAAAAAAAAAAAAAAA\nCCCCCCCCCCCCCCCC\n")
        second.write_text("CCCCCCCCCCCCCCCC\nGGGGGGGGGGGGGGGG\n")
        assert load_whitelist([str(first), str(second)]) == {
            "AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC", "GGGGGGGGGGGGGGGG"}


class TestSelectCellBarcodes:
    def _sample(self, tmp_path, rnd, cell_count=30, noise=400):
        """Cells with a spread of read counts, plus a tail of one-off background barcodes."""
        cells = [random_barcode(rnd) for _ in range(cell_count)]
        rows = []
        for index, barcode in enumerate(cells):
            for i in range(rnd.randrange(40, 200)):
                rows.append(("r%d_%d" % (index, i), barcode))
        for i in range(noise):
            rows.append(("n%d" % i, random_barcode(rnd)))
        rnd.shuffle(rows)
        raw = tmp_path / "raw_0.tsv"
        write_raw_table(raw, rows)
        return str(raw), cells

    def test_detects_the_cell_barcodes(self, tmp_path):
        rnd = random.Random(1)
        raw, cells = self._sample(tmp_path, rnd)
        out = str(tmp_path / "cells.tsv")
        select_cell_barcodes([raw], out, None, n_cells=30)
        assert sorted(read_barcodes(out)) == sorted(cells)

    def test_cell_number_is_estimated(self, tmp_path):
        """--n_cells auto: the knee of the count distribution stands in for a given count."""
        rnd = random.Random(2)
        raw, cells = self._sample(tmp_path, rnd)
        out = str(tmp_path / "cells.tsv")
        select_cell_barcodes([raw], out, None, n_cells="auto")
        assert sorted(read_barcodes(out)) == sorted(cells)

    def test_whitelist_filters_candidates(self, tmp_path):
        """A frequent barcode outside the pool must not be taken for a cell."""
        rnd = random.Random(3)
        raw, cells = self._sample(tmp_path, rnd)
        impostor = "TTTTTTTTTTTTTTTT"
        with open(raw, "a") as handle:
            for i in range(5000):  # more abundant than any real cell
                handle.write("imp%d\t%s\tACGTACGTACGT\t0\tTrue\t+\t50\t30\n" % (i, impostor))
        pool = tmp_path / "pool.txt"
        pool.write_text("\n".join(cells) + "\n")
        out = str(tmp_path / "cells.tsv")
        select_cell_barcodes([raw], out, [str(pool)], n_cells=30)
        detected = read_barcodes(out)
        assert impostor not in detected
        assert sorted(detected) == sorted(cells)

    def test_returns_the_output_path(self, tmp_path):
        rnd = random.Random(4)
        raw, _ = self._sample(tmp_path, rnd)
        out = str(tmp_path / "cells.tsv")
        assert select_cell_barcodes([raw], out, None, n_cells=30) == out

    def test_stats_file_is_written(self, tmp_path):
        rnd = random.Random(5)
        raw, cells = self._sample(tmp_path, rnd)
        stats_file = tmp_path / "cells.stats"
        select_cell_barcodes([raw], str(tmp_path / "cells.tsv"), None,
                             n_cells=30, stats_file=str(stats_file))
        stats = dict(line.split("\t") for line in stats_file.read_text().splitlines())
        assert int(stats["Cell barcodes detected"]) == len(cells)
        assert int(stats["Distinct barcodes extracted"]) > len(cells)

    def test_empty_input(self, tmp_path):
        raw, out = tmp_path / "raw.tsv", tmp_path / "cells.tsv"
        write_raw_table(raw, [])
        select_cell_barcodes([str(raw)], str(out), None, n_cells="auto")
        assert read_barcodes(out) == []
