############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import random

import pytest

from isoquant_lib.barcode_calling.cell_selection import (
    NOSEQ,
    CellBarcodeSelector,
    count_barcodes,
    estimate_cell_number,
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


class TestEstimateCellNumber:
    def test_too_few_barcodes(self):
        assert estimate_cell_number([9, 8, 7]) == 3

    def test_finds_knee(self):
        # 500 real cells at high counts, then a long tail of noise
        counts = [1000] * 500 + [3] * 20000
        assert 400 <= estimate_cell_number(counts) <= 600

    def test_ignores_singletons(self):
        assert estimate_cell_number([1] * 1000) == 0

    def test_flat_distribution_has_no_knee(self):
        """Equally abundant barcodes are all equally cell-like, so take all of them.

        The chord-distance argmax is 0 everywhere here; reading that as a knee would
        silently report a single cell.
        """
        assert estimate_cell_number([60] * 30 + [1] * 400) == 30

    def test_knee_survives_a_noise_tail(self):
        counts = [1000 - i for i in range(200)] + [2] * 5000
        assert 150 <= estimate_cell_number(counts) <= 260


class TestBarcodeTally:
    def test_rejects_malformed(self):
        selector = CellBarcodeSelector(barcode_length=16)
        assert selector.add_barcode("ACGTACGTACGTACGT")
        assert not selector.add_barcode("ACGTNCGTACGTACGT")      # N is not a nucleotide
        assert not selector.add_barcode("ACGT")                  # too short
        assert not selector.add_barcode("ACGTACGTACGTACGTACGT")  # too long
        assert selector.malformed_count == 3
        assert dict(selector.counts) == {"ACGTACGTACGTACGT": 1}

    def test_trailing_base_is_trimmed(self):
        """Windows one base longer than the barcode are truncated, as Badger did."""
        selector = CellBarcodeSelector(barcode_length=16)
        assert selector.add_barcode("ACGTACGTACGTACGTA")
        assert dict(selector.counts) == {"ACGTACGTACGTACGT": 1}

    def test_merge_counts(self):
        selector = CellBarcodeSelector(barcode_length=16)
        selector.add_barcode("ACGTACGTACGTACGT")
        selector.merge_counts({"ACGTACGTACGTACGT": 4, "TTTTGGGGCCCCAAAA": 2}, malformed=7)
        assert selector.counts["ACGTACGTACGTACGT"] == 5
        assert selector.counts["TTTTGGGGCCCCAAAA"] == 2
        assert selector.malformed_count == 7


class TestSelect:
    def _selector(self, counts):
        selector = CellBarcodeSelector(barcode_length=16)
        selector.counts = dict(counts)
        return selector

    def test_whitelist_filters_candidates(self):
        rnd = random.Random(1)
        real = [random_barcode(rnd) for _ in range(5)]
        impostor = random_barcode(rnd)
        counts = {bc: 100 for bc in real}
        counts[impostor] = 1000  # most abundant, but not whitelisted
        centers = self._selector(counts).select(n_cells=5, whitelist=set(real))
        assert impostor not in centers
        assert sorted(centers) == sorted(real)

    def test_abundant_artifact_does_not_raise_the_cutoff(self):
        """The cutoff comes from the candidates, not from everything that was observed.

        A frequent non-whitelisted artifact used to inflate mean(top n_cells) and push
        genuine cells below the cutoff.
        """
        rnd = random.Random(11)
        real = [random_barcode(rnd) for _ in range(20)]
        counts = {bc: 60 for bc in real}
        counts[random_barcode(rnd)] = 100000  # artifact, orders of magnitude more abundant
        centers = self._selector(counts).select(n_cells=20, whitelist=set(real))
        assert sorted(centers) == sorted(real)

    def test_fewer_barcodes_than_requested_cells(self):
        """Badger raised IndexError here; asking for more cells than exist must be safe."""
        rnd = random.Random(2)
        counts = {random_barcode(rnd): 50 for _ in range(3)}
        assert len(self._selector(counts).select(n_cells=1000)) == 3

    def test_overshooting_n_cells_does_not_scrape_the_noise_floor(self):
        """Padding towards n_cells must stop above the singletons."""
        rnd = random.Random(12)
        real = [random_barcode(rnd) for _ in range(50)]
        counts = {bc: 100 for bc in real}
        counts.update({random_barcode(rnd): 1 for _ in range(5000)})
        centers = self._selector(counts).select(n_cells=2000)
        assert sorted(centers) == sorted(real)

    def test_no_barcodes(self):
        assert self._selector({}).select(n_cells=10) == []

    def test_cutoff_uses_the_most_abundant_barcodes(self):
        """The cutoff is derived from the top n_cells, not an arbitrary dict slice."""
        rnd = random.Random(3)
        real = [random_barcode(rnd) for _ in range(10)]
        noise = [random_barcode(rnd) for _ in range(5000)]
        counts = {bc: 10000 for bc in real}
        counts.update({bc: 1 for bc in noise})
        centers = self._selector(counts).select(n_cells=10, whitelist=set(real) | set(noise))
        # cutoff is mean(top 10)/5 = 2000, so none of the singletons may be picked up
        assert sorted(centers) == sorted(real)

    def test_stats(self):
        rnd = random.Random(4)
        real = [random_barcode(rnd) for _ in range(5)]
        selector = self._selector({bc: 100 for bc in real})
        selector.select(n_cells=5)
        stats = selector.stats()
        assert stats["Cell barcodes detected"] == 5
        assert stats["Reads exactly matching a cell barcode"] == 500


class TestCounting:
    def test_header_is_not_counted_as_a_read(self, tmp_path):
        """Raw tables carry a plain, non-commented header, so it must be spotted by content."""
        raw = tmp_path / "raw.tsv"
        write_raw_table(raw, [("r1", "ACGTACGTACGTACGT"), ("r2", "ACGTACGTACGTACGT")])
        graph = CellBarcodeSelector(barcode_length=16)
        assert count_barcodes([str(raw)], graph) == 2
        assert graph.malformed_count == 0
        assert dict(graph.counts) == {"ACGTACGTACGTACGT": 2}

    def test_missing_barcodes_are_skipped(self, tmp_path):
        raw = tmp_path / "raw.tsv"
        write_raw_table(raw, [("r1", "ACGTACGTACGTACGT"), ("r2", NOSEQ)])
        graph = CellBarcodeSelector(barcode_length=16)
        assert count_barcodes([str(raw)], graph) == 2
        assert dict(graph.counts) == {"ACGTACGTACGTACGT": 1}
        assert graph.malformed_count == 0

    def test_counts_accumulate_across_files(self, tmp_path):
        first, second = tmp_path / "a.tsv", tmp_path / "b.tsv"
        write_raw_table(first, [("r1", "ACGTACGTACGTACGT")])
        write_raw_table(second, [("r2", "ACGTACGTACGTACGT"), ("r3", "TTTTGGGGCCCCAAAA")])
        graph = CellBarcodeSelector(barcode_length=16)
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
