############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import random

import pytest

from isoquant_lib.barcode_calling.barcode_graph import (
    BarcodeGraph,
    bounded_distance,
    estimate_cell_number,
    optimal_kmer_size,
)


def random_barcode(rnd, length=16):
    return "".join(rnd.choice("ACGT") for _ in range(length))


def substitute(rnd, barcode, position=None):
    position = rnd.randrange(len(barcode)) if position is None else position
    replacement = rnd.choice([c for c in "ACGT" if c != barcode[position]])
    return barcode[:position] + replacement + barcode[position + 1:]


class TestBoundedDistance:
    def test_identical(self):
        assert bounded_distance("ACGTACGTACGTACGT", "ACGTACGTACGTACGT", 1) == 0

    def test_single_substitution(self):
        assert bounded_distance("ACGTACGTACGTACGT", "ACGTACGTACGTACGA", 1) == 1

    def test_two_substitutions_exceed_threshold(self):
        assert bounded_distance("ACGTACGTACGTACGT", "ACGTACGTACGTACAA", 1) > 1

    def test_shifted_window_counts_as_one_edit(self):
        """A deletion inside the barcode shifts the window and appends a spurious base."""
        assert bounded_distance("CGTACGTACGTACGTA", "ACGTACGTACGTACGT", 1) == 1

    def test_symmetric(self):
        rnd = random.Random(5)
        for _ in range(50):
            seq1, seq2 = random_barcode(rnd), random_barcode(rnd)
            assert bounded_distance(seq1, seq2, 2) == bounded_distance(seq2, seq1, 2)


class TestEstimateCellNumber:
    def test_too_few_barcodes(self):
        assert estimate_cell_number([9, 8, 7]) == 3

    def test_finds_knee(self):
        # 500 real cells at high counts, then a long tail of noise
        counts = [1000] * 500 + [3] * 20000
        estimated = estimate_cell_number(counts)
        assert 400 <= estimated <= 600

    def test_ignores_singletons(self):
        assert estimate_cell_number([1] * 1000) == 0


class TestKmerSizeSelection:
    """k is picked as large as the q-gram bound allows, which is what makes lookups cheap."""

    @pytest.mark.parametrize("length, threshold, expected", [
        (16, 1, 8), (16, 2, 5), (16, 3, 4), (25, 1, 12), (14, 1, 7),
    ])
    def test_expected_sizes(self, length, threshold, expected):
        assert optimal_kmer_size(length, threshold) == expected

    def test_never_below_the_floor(self):
        assert optimal_kmer_size(8, 5) == 4

    @pytest.mark.parametrize("threshold", [1, 2, 3])
    def test_bound_stays_usable(self, threshold):
        assert BarcodeGraph(threshold=threshold, barcode_length=16)._qgram_bound() >= 1

    @pytest.mark.parametrize("threshold", [1, 2, 3])
    def test_bound_never_rejects_a_close_pair(self, threshold):
        """Every pair within `threshold` edits must share at least `bound` k-mers.

        If this ever fails the k-mer filter silently drops real neighbours, which is the
        exact failure mode that makes whitelist matching miss error-bearing barcodes.
        """
        rnd = random.Random(77)
        graph = BarcodeGraph(threshold=threshold, barcode_length=16)
        k, bound = graph.kmer_size, graph._qgram_bound()
        for _ in range(300):
            seq1 = random_barcode(rnd)
            seq2 = seq1
            for _ in range(threshold):
                seq2 = substitute(rnd, seq2)
            if bounded_distance(seq1, seq2, threshold) > threshold:
                continue
            kmers1 = [seq1[i:i + k] for i in range(len(seq1) - k + 1)]
            kmers2 = [seq2[i:i + k] for i in range(len(seq2) - k + 1)]
            assert sum(kmers2.count(x) for x in kmers1) >= bound

    def test_explicit_size_is_respected(self):
        assert BarcodeGraph(threshold=1, barcode_length=16, kmer_size=6).kmer_size == 6

    def test_clustering_is_independent_of_kmer_size(self):
        rnd = random.Random(78)
        cells = [random_barcode(rnd) for _ in range(30)]
        counts = {}
        for barcode in cells:
            counts[barcode] = 60
            for _ in range(5):
                counts[substitute(rnd, barcode)] = 2
        for _ in range(300):
            counts[random_barcode(rnd)] = 1

        results = []
        for kmer_size in (6, 8):
            graph = BarcodeGraph(threshold=1, kmer_size=kmer_size)
            graph.counts = dict(counts)
            graph.select_cluster_centers(n_cells=30, whitelist=set(cells))
            graph.cluster_from_centers(rounds=2)
            results.append(graph.clustering)
        assert results[0] == results[1]


class TestCounting:
    def test_rejects_malformed(self):
        graph = BarcodeGraph(barcode_length=16)
        assert graph.add_barcode("ACGTACGTACGTACGT")
        assert not graph.add_barcode("ACGTNCGTACGTACGT")  # N is not a nucleotide
        assert not graph.add_barcode("ACGT")              # too short
        assert not graph.add_barcode("ACGTACGTACGTACGTACGT")  # too long
        assert graph.malformed_count == 3
        assert dict(graph.counts) == {"ACGTACGTACGTACGT": 1}

    def test_trailing_base_is_trimmed(self):
        """Windows one base longer than the barcode are truncated, as Badger does."""
        graph = BarcodeGraph(barcode_length=16)
        assert graph.add_barcode("ACGTACGTACGTACGTA")
        assert dict(graph.counts) == {"ACGTACGTACGTACGT": 1}

    def test_merge_counts(self):
        graph = BarcodeGraph(barcode_length=16)
        graph.add_barcode("ACGTACGTACGTACGT")
        graph.merge_counts({"ACGTACGTACGTACGT": 4, "TTTTGGGGCCCCAAAA": 2}, malformed=7)
        assert graph.counts["ACGTACGTACGTACGT"] == 5
        assert graph.counts["TTTTGGGGCCCCAAAA"] == 2
        assert graph.malformed_count == 7


class TestCenterSelection:
    def _graph(self, counts):
        graph = BarcodeGraph(barcode_length=16)
        graph.counts = dict(counts)
        return graph

    def test_whitelist_filters_candidates(self):
        rnd = random.Random(1)
        real = [random_barcode(rnd) for _ in range(5)]
        impostor = random_barcode(rnd)
        counts = {bc: 100 for bc in real}
        counts[impostor] = 1000  # most abundant, but not whitelisted
        graph = self._graph(counts)
        centers = graph.select_cluster_centers(n_cells=5, whitelist=set(real))
        assert impostor not in centers
        assert sorted(centers) == sorted(real)

    def test_fewer_barcodes_than_requested_cells(self):
        """Badger raised IndexError here; asking for more cells than exist must be safe."""
        rnd = random.Random(2)
        counts = {random_barcode(rnd): 50 for _ in range(3)}
        graph = self._graph(counts)
        centers = graph.select_cluster_centers(n_cells=1000)
        assert len(centers) == 3

    def test_no_barcodes(self):
        assert self._graph({}).select_cluster_centers(n_cells=10) == []

    def test_cutoff_uses_the_most_abundant_barcodes(self):
        """The count cutoff is derived from the top n_cells, not an arbitrary dict slice."""
        rnd = random.Random(3)
        real = [random_barcode(rnd) for _ in range(10)]
        noise = [random_barcode(rnd) for _ in range(5000)]
        counts = {bc: 10000 for bc in real}
        counts.update({bc: 1 for bc in noise})
        graph = self._graph(counts)
        centers = graph.select_cluster_centers(n_cells=10, whitelist=set(real) | set(noise))
        # cutoff is mean(top 10)/5 = 2000, so none of the singletons may be picked up
        assert sorted(centers) == sorted(real)

    def test_supplied_cell_barcodes_take_precedence(self):
        rnd = random.Random(4)
        real = [random_barcode(rnd) for _ in range(3)]
        graph = self._graph({bc: 1 for bc in real})
        unobserved = random_barcode(rnd)
        centers = graph.select_cluster_centers(true_barcodes=set(real) | {unobserved})
        assert sorted(centers) == sorted(real)


class TestClustering:
    def test_error_variants_are_corrected(self):
        rnd = random.Random(11)
        cells = [random_barcode(rnd) for _ in range(20)]
        graph = BarcodeGraph(threshold=1)
        variants = {}
        for barcode in cells:
            graph.counts[barcode] = 100
            for position in range(0, 16, 4):
                variant = substitute(rnd, barcode, position)
                graph.counts[variant] = 2
                variants[variant] = barcode

        graph.select_cluster_centers(n_cells=20, whitelist=set(cells))
        graph.cluster_from_centers(rounds=2)
        assignments = graph.get_assignments()
        for variant, origin in variants.items():
            assert assignments[variant] == origin

    def test_equidistant_barcode_stays_unassigned(self):
        center1 = "ACGTACGTACGTACGT"
        center2 = "ACGTACGTACGTACAA"
        middle = "ACGTACGTACGTACGA"  # one edit from either center
        graph = BarcodeGraph(threshold=1)
        graph.counts = {center1: 100, center2: 100, middle: 3}

        graph.select_cluster_centers(n_cells=2, whitelist={center1, center2})
        graph.cluster_from_centers(rounds=2)
        assert middle not in graph.get_assignments()
        assert graph.clustering[middle][0] is None

    def test_unreachable_barcode_is_not_assigned(self):
        rnd = random.Random(12)
        center = random_barcode(rnd)
        unrelated = random_barcode(rnd)
        graph = BarcodeGraph(threshold=1)
        graph.counts = {center: 100, unrelated: 1}

        graph.select_cluster_centers(n_cells=1, whitelist={center})
        graph.cluster_from_centers(rounds=2)
        assert unrelated not in graph.get_assignments()

    def test_second_round_reaches_via_an_observed_intermediate(self):
        center = "ACGTACGTACGTACGT"
        hop1 = "ACGTACGTACGTACGA"
        hop2 = "ACGTACGTACGTAAGA"
        graph = BarcodeGraph(threshold=1)
        graph.counts = {center: 100, hop1: 5, hop2: 2}
        graph.select_cluster_centers(n_cells=1, whitelist={center})

        graph.cluster_from_centers(rounds=1)
        assert hop2 not in graph.get_assignments()

        graph.cluster_from_centers(rounds=2)
        assert graph.get_assignments()[hop2] == center


class TestImplementationsAgree:
    """cluster_from_centers must reproduce construct_graph + cluster exactly."""

    @pytest.mark.parametrize("seed", [21, 22, 23])
    def test_same_clustering(self, seed):
        rnd = random.Random(seed)
        cells = [random_barcode(rnd) for _ in range(40)]
        whitelist = set(cells) | {random_barcode(rnd) for _ in range(500)}
        counts = {}
        for barcode in cells:
            counts[barcode] = rnd.randrange(50, 200)
            for _ in range(4):
                variant = substitute(rnd, barcode)
                counts[variant] = counts.get(variant, 0) + 2
                counts[substitute(rnd, variant)] = 1
        for _ in range(200):
            counts[random_barcode(rnd)] = 1

        results = []
        for use_full_graph in (True, False):
            graph = BarcodeGraph(threshold=1)
            graph.counts = dict(counts)
            if use_full_graph:
                graph.construct_graph(threads=1)
                graph.select_cluster_centers(n_cells=40, whitelist=whitelist)
                graph.cluster(rounds=2)
            else:
                graph.select_cluster_centers(n_cells=40, whitelist=whitelist)
                graph.cluster_from_centers(rounds=2)
            results.append(graph.clustering)

        assert results[0] == results[1]

    def test_low_complexity_barcodes_are_not_pruned(self):
        """Repetitive barcodes inflate shared k-mer counts; hit pruning must stay disabled."""
        center = "CAACGAACCCCGCCCC"
        variant = "CAACGAACCCCCCCCC"
        counts = {center: 100, variant: 2}
        # plenty of poly-C decoys, which is what used to push the true neighbour out
        for i in range(200):
            counts["C" * 12 + "ACGT"[i % 4] + "CCC"] = 1
            counts["CCCC" + "ACGT"[i % 4] + "C" * 11] = 1

        graph = BarcodeGraph(threshold=1)
        graph.counts = dict(counts)
        graph.select_cluster_centers(n_cells=1, whitelist={center})
        graph.cluster_from_centers(rounds=2)
        assert graph.get_assignments()[variant] == center


class TestParallelGraphConstruction:
    def test_matches_single_thread(self):
        rnd = random.Random(31)
        cells = [random_barcode(rnd) for _ in range(30)]
        counts = {}
        for barcode in cells:
            counts[barcode] = 40
            for position in range(0, 16, 3):
                counts[substitute(rnd, barcode, position)] = 2

        graphs = []
        for threads in (1, 4):
            graph = BarcodeGraph(threshold=1)
            graph.counts = dict(counts)
            graph.construct_graph(threads=threads)
            graphs.append(graph)

        assert graphs[0].dists == graphs[1].dists
        assert ({k: sorted(v) for k, v in graphs[0].edges.items()} ==
                {k: sorted(v) for k, v in graphs[1].edges.items()})
