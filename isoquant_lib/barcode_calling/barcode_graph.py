############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Edit-distance graph over observed barcodes.

Ported from Badger (https://github.com/algbio/Badger). The idea: with a stock 10x
whitelist, per-read whitelist matching degenerates into exact matching, so instead of
searching millions of candidate barcodes per read we

  1. extract barcode windows verbatim (see TenXBarcodeDetector, whitelist_matching=False),
  2. count how often each distinct window is observed over the whole run,
  3. pick the cluster centers -- the barcodes that were actually sequenced -- from the
     count distribution, keeping only those present in the whitelist,
  4. correct every other observed barcode onto a center by walking an edit-distance graph.

The whitelist is only a membership filter on candidate centers, never a per-read search
space. See .claude/BARCODE_GRAPH_CORRECTION.md.
"""

import logging
import math
from collections import defaultdict
from statistics import mean
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import editdistance
import numpy

from .indexers import ArrayKmerIndexer

logger = logging.getLogger('IsoQuant')

# Cluster centers must be observed more often than mean(top n_cells counts) / this
CENTER_COUNT_FRACTION = 5.0
# ... and never less often than this, however sparse the data
MIN_CENTER_COUNT = 5
# Barcodes rarer than this are not considered when estimating the cell number
MIN_COUNT_FOR_ESTIMATION = 2
# Number of distinct barcodes handed to one graph construction worker
BC_CHUNK_SIZE = 10000

NUCLEOTIDES = frozenset("ACGT")


def bounded_distance(seq1: str, seq2: str, threshold: int) -> int:
    """Edit distance between two barcodes, tolerating a shifted last base.

    Barcode windows are cut at a fixed offset from the R1 primer, so an indel inside the
    barcode shifts the window and leaves a spurious trailing base. Comparing the truncated
    variants (as Badger does) absorbs that. Returns threshold + 1 when the true distance
    exceeds threshold -- callers only ever test against the threshold.
    """
    # Fast path: equal-length strings differing in at most `threshold` positions are within
    # `threshold` edits, no alignment needed. This covers plain substitutions, the common case.
    mismatches = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 != c2:
            mismatches += 1
            if mismatches > threshold:
                break
    if mismatches <= threshold:
        return mismatches

    return min(editdistance.eval(seq1, seq2),
               editdistance.eval(seq1[:-1], seq2),
               editdistance.eval(seq1, seq2[:-1]))


def estimate_cell_number(sorted_counts: Sequence[int]) -> int:
    """Estimate the number of cell-associated barcodes from the count distribution.

    Knee of the log-log count-versus-rank curve, located as the point furthest below the
    chord joining its ends. Used when --n_cells is not given.
    """
    counts = [c for c in sorted_counts if c >= MIN_COUNT_FOR_ESTIMATION]
    if len(counts) < 10:
        return len(counts)

    x = numpy.log10(numpy.arange(1, len(counts) + 1, dtype=numpy.float64))
    y = numpy.log10(numpy.asarray(counts, dtype=numpy.float64))
    # distance from each point to the chord, signed so that points below it are positive
    dx, dy = x[-1] - x[0], y[-1] - y[0]
    norm = math.hypot(dx, dy)
    if norm == 0.0:
        return len(counts)
    distance = (dx * (y[0] - y) - (x[0] - x) * dy) / norm
    return int(numpy.argmax(distance)) + 1


# Populated once per worker process by _init_compare_worker; keeping the index here rather
# than passing it to every task means it is pickled once per worker, not once per chunk.
_WORKER_STATE: Dict[str, object] = {}


def _compare_chunk(barcode_chunk: List[str]) -> List[Tuple[str, str, int]]:
    """Find graph edges for one chunk of barcodes. Runs in a worker process."""
    index = _WORKER_STATE["index"]
    threshold = _WORKER_STATE["threshold"]
    min_kmers = _WORKER_STATE["min_kmers"]
    hits_delta = _WORKER_STATE["hits_delta"]

    edges = []
    for barcode in barcode_chunk:
        for candidate, _, _ in index.get_occurrences(barcode, max_hits=0, min_kmers=min_kmers,
                                                     hits_delta=hits_delta):
            # every pair is reported from both sides, keep it once
            if candidate <= barcode:
                continue
            dist = bounded_distance(barcode, candidate, threshold)
            if dist <= threshold:
                edges.append((barcode, candidate, dist))
    return edges


def _set_worker_state(index, threshold: int, min_kmers: int, hits_delta: int) -> None:
    _WORKER_STATE["index"] = index
    _WORKER_STATE["threshold"] = threshold
    _WORKER_STATE["min_kmers"] = min_kmers
    _WORKER_STATE["hits_delta"] = hits_delta


def _init_compare_worker(log_file, log_level, index, threshold: int,
                         min_kmers: int, hits_delta: int) -> None:
    from ..common import setup_worker_logging
    setup_worker_logging(log_file, log_level)
    _set_worker_state(index, threshold, min_kmers, hits_delta)


class BarcodeGraph:
    """Undirected graph connecting observed barcodes within `threshold` edits."""

    def __init__(self, threshold: int = 1, barcode_length: int = 16, kmer_size: int = 6):
        self.threshold: int = threshold
        self.barcode_length: int = barcode_length
        self.kmer_size: int = kmer_size
        self.counts: Dict[str, int] = defaultdict(int)
        self.edges: Dict[str, List[str]] = defaultdict(list)
        self.dists: Dict[Tuple[str, str], int] = {}
        # observed barcode -> (center, hop); center is None for ambiguous barcodes
        self.clustering: Dict[str, Tuple[Optional[str], int]] = {}
        self.centers: List[str] = []
        self.malformed_count: int = 0

    # ------------------------------------------------------------------ counting

    def add_barcode(self, barcode: str) -> bool:
        """Tally one observed barcode. Returns False if it was unusable."""
        if len(barcode) == self.barcode_length + 1:
            barcode = barcode[:-1]
        if len(barcode) != self.barcode_length or not set(barcode) <= NUCLEOTIDES:
            self.malformed_count += 1
            return False
        self.counts[barcode] += 1
        return True

    def add_barcodes(self, barcodes: Iterable[str]) -> None:
        for barcode in barcodes:
            self.add_barcode(barcode)

    def merge_counts(self, counts: Dict[str, int], malformed: int = 0) -> None:
        """Merge a partial count table, e.g. produced by another process."""
        for barcode, count in counts.items():
            self.counts[barcode] += count
        self.malformed_count += malformed

    # ------------------------------------------------- graph construction

    def _kmers_per_barcode(self) -> int:
        return self.barcode_length - self.kmer_size + 1

    def _qgram_bound(self) -> int:
        """Minimum shared k-mer occurrences two barcodes within `threshold` edits must have.

        Standard q-gram lemma: each edit destroys at most `kmer_size` of the
        `barcode_length - kmer_size + 1` k-mers.
        """
        return max(self._kmers_per_barcode() - self.kmer_size * self.threshold, 1)

    def _no_hit_pruning(self) -> int:
        """hits_delta that provably keeps every candidate the q-gram bound admits.

        The indexers drop candidates more than hits_delta shared k-mers behind the best one,
        which is exactly the pruning that makes whitelist matching miss error-bearing
        barcodes. The graph needs every candidate, so the delta has to exceed the largest
        possible count: each of the query's k-mers can match at most that many occurrences
        in one indexed barcode, and low-complexity barcodes really do reach it.
        """
        return self._kmers_per_barcode() ** 2

    def build_index(self, barcodes: Optional[Sequence[str]] = None) -> ArrayKmerIndexer:
        if barcodes is None:
            barcodes = list(self.counts.keys())
        return ArrayKmerIndexer(barcodes, kmer_size=self.kmer_size)

    def add_edge(self, seq1: str, seq2: str, dist: int) -> None:
        self.edges[seq1].append(seq2)
        self.edges[seq2].append(seq1)
        self.dists[(seq1, seq2)] = dist
        self.dists[(seq2, seq1)] = dist

    def construct_graph(self, threads: int = 1) -> None:
        """Connect every pair of observed barcodes within `threshold` edits."""
        barcodes = list(self.counts.keys())
        logger.info("Constructing barcode graph over %d distinct barcodes" % len(barcodes))
        index = self.build_index(barcodes)
        min_kmers = self._qgram_bound()
        hits_delta = self._no_hit_pruning()

        if threads <= 1:
            _set_worker_state(index, self.threshold, min_kmers, hits_delta)
            for seq1, seq2, dist in _compare_chunk(barcodes):
                self.add_edge(seq1, seq2, dist)
        else:
            for seq1, seq2, dist in self._construct_graph_in_parallel(barcodes, index,
                                                                     min_kmers, hits_delta, threads):
                self.add_edge(seq1, seq2, dist)

        logger.info("Barcode graph has %d nodes with at least one edge, %d edges" %
                    (len(self.edges), len(self.dists) // 2))

    def _construct_graph_in_parallel(self, barcodes: List[str], index: ArrayKmerIndexer,
                                     min_kmers: int, hits_delta: int,
                                     threads: int) -> Iterable[Tuple[str, str, int]]:
        import concurrent.futures
        import multiprocessing
        from ..common import _get_log_params

        log_file, log_level = _get_log_params()
        chunks = [barcodes[i:i + BC_CHUNK_SIZE] for i in range(0, len(barcodes), BC_CHUNK_SIZE)]
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=threads,
                mp_context=multiprocessing.get_context('spawn'),
                initializer=_init_compare_worker,
                initargs=(log_file, log_level, index, self.threshold, min_kmers, hits_delta)) as proc:
            for chunk_edges in proc.map(_compare_chunk, chunks, chunksize=1):
                for edge in chunk_edges:
                    yield edge

    # ------------------------------------------------- center selection

    def sorted_barcodes(self) -> List[str]:
        """Distinct observed barcodes, most frequent first."""
        return sorted(self.counts, key=lambda bc: (self.counts[bc], bc), reverse=True)

    def select_cluster_centers(self, n_cells: Optional[int] = None,
                               whitelist: Optional[Set[str]] = None,
                               true_barcodes: Optional[Set[str]] = None,
                               interval: int = 25) -> List[str]:
        """Pick the barcodes that represent actually sequenced cells."""
        if true_barcodes:
            centers = [bc for bc in true_barcodes if bc in self.counts]
            logger.info("Using %d of %d supplied cell barcodes as cluster centers "
                        "(%d were never observed)" %
                        (len(centers), len(true_barcodes), len(true_barcodes) - len(centers)))
            self.centers = centers
            return centers

        by_count = self.sorted_barcodes()
        if not by_count:
            self.centers = []
            return []

        if n_cells is None:
            n_cells = estimate_cell_number([self.counts[bc] for bc in by_count])
            logger.info("Estimated number of cell-associated barcodes: %d" % n_cells)
        n_cells = max(1, min(n_cells, len(by_count)))

        cutoff = max(mean(self.counts[bc] for bc in by_count[:n_cells]) / CENTER_COUNT_FRACTION,
                     MIN_CENTER_COUNT)
        upper = n_cells + n_cells * interval / 100.0
        lower = n_cells - n_cells * interval / 100.0
        logger.info("Selecting up to %d cluster centers, minimal barcode count %.1f" % (upper, cutoff))

        centers: List[str] = []
        i = 0
        while i < len(by_count) and self.counts[by_count[i]] > cutoff and len(centers) <= upper:
            if whitelist is None or by_count[i] in whitelist:
                centers.append(by_count[i])
            i += 1
        # if the count cutoff was too aggressive, keep going until the lower bound is reached
        while i < len(by_count) and len(centers) < lower:
            if whitelist is None or by_count[i] in whitelist:
                centers.append(by_count[i])
            i += 1

        logger.info("Selected %d cluster centers, rarest one observed %d times" %
                    (len(centers), self.counts[centers[-1]] if centers else 0))
        self.centers = centers
        return centers

    # ------------------------------------------------- clustering

    def cluster(self, centers: Optional[Sequence[str]] = None, rounds: int = 2) -> None:
        """Grow clusters outwards from the centers, `rounds` graph hops at a time.

        A barcode reached from two different centers within the same hop is genuinely
        ambiguous and stays unassigned; a barcode reached at an earlier hop keeps its
        first (closer) center.
        """
        if centers is None:
            centers = self.centers
        clustering: Dict[str, Tuple[Optional[str], int]] = {bc: (bc, 0) for bc in centers}
        frontier = list(centers)

        for hop in range(1, rounds + 1):
            claims: Dict[str, Set[str]] = defaultdict(set)
            for node in frontier:
                center = clustering[node][0]
                for neighbor in self.edges.get(node, ()):
                    if neighbor not in clustering:
                        claims[neighbor].add(center)

            next_frontier = []
            ambiguous = 0
            for neighbor, claiming_centers in claims.items():
                if len(claiming_centers) == 1:
                    clustering[neighbor] = (next(iter(claiming_centers)), hop)
                    next_frontier.append(neighbor)
                else:
                    # claimed by several centers at the same distance, cannot be resolved
                    clustering[neighbor] = (None, hop)
                    ambiguous += 1
            logger.info("Clustering round %d: %d barcodes assigned, %d ambiguous" %
                        (hop, len(next_frontier), ambiguous))
            frontier = next_frontier
            if not frontier:
                break

        self.clustering = clustering

    def cluster_from_centers(self, centers: Optional[Sequence[str]] = None, rounds: int = 2) -> None:
        """Cluster without ever materialising the graph.

        The graph built by construct_graph() is only ever consumed by a walk outwards from
        the centers, so the edges nobody reaches need not be found. Each round indexes just
        the previous round's frontier -- the centers first (a few thousand sequences), then
        whatever they claimed -- and queries the still-unassigned barcodes against it.
        Produces exactly the same clustering as construct_graph() + cluster(), but the
        index never has to hold all observed barcodes at once.
        """
        if centers is None:
            centers = self.centers
        clustering: Dict[str, Tuple[Optional[str], int]] = {bc: (bc, 0) for bc in centers}
        min_kmers = self._qgram_bound()
        hits_delta = self._no_hit_pruning()
        frontier = list(centers)

        for hop in range(1, rounds + 1):
            if not frontier:
                break
            index = self.build_index(frontier)
            next_frontier = []
            ambiguous = 0
            for barcode in self.counts:
                if barcode in clustering:
                    continue
                claiming_centers: Set[str] = set()
                for candidate, _, _ in index.get_occurrences(barcode, max_hits=0,
                                                             min_kmers=min_kmers,
                                                             hits_delta=hits_delta):
                    if bounded_distance(barcode, candidate, self.threshold) <= self.threshold:
                        claiming_centers.add(clustering[candidate][0])
                        if len(claiming_centers) > 1:
                            break
                if not claiming_centers:
                    continue
                if len(claiming_centers) == 1:
                    clustering[barcode] = (next(iter(claiming_centers)), hop)
                    next_frontier.append(barcode)
                else:
                    # reached from several centers at the same distance, cannot be resolved
                    clustering[barcode] = (None, hop)
                    ambiguous += 1
            logger.info("Clustering round %d: %d barcodes assigned, %d ambiguous" %
                        (hop, len(next_frontier), ambiguous))
            frontier = next_frontier

        self.clustering = clustering

    def get_assignments(self) -> Dict[str, str]:
        """Observed barcode -> corrected barcode, for every barcode that could be resolved."""
        return {bc: center for bc, (center, _) in self.clustering.items() if center is not None}

    def assignment_stats(self) -> Dict[str, int]:
        assigned_reads = 0
        ambiguous_reads = 0
        for bc, count in self.counts.items():
            center = self.clustering.get(bc, (None, -1))[0]
            if center is not None:
                assigned_reads += count
            elif bc in self.clustering:
                ambiguous_reads += count
        return {
            "Distinct barcodes observed": len(self.counts),
            "Malformed barcodes skipped": self.malformed_count,
            "Cluster centers": len(self.centers),
            "Distinct barcodes corrected": len(self.get_assignments()),
            "Distinct barcodes ambiguous": sum(1 for c, _ in self.clustering.values() if c is None),
            "Reads with corrected barcode": assigned_reads,
            "Reads with ambiguous barcode": ambiguous_reads,
        }
