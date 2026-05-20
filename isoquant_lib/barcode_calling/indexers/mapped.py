############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Disk-backed (mmap) sparse 3-anchor k-mer indexer.

A drop-in successor to :class:`SharedMemorySparseAnchorIndexer` whose large
buffers live in a directory of memory-mapped files instead of POSIX
shared-memory segments. Resident-set size becomes whatever the page cache
chooses to retain, capping RAM usage for 10B-50B-scale whitelists at the
cost of disk-backed random reads at query time.

File layout (versioned via ``meta.json``):

    {cache_dir}/
      meta.json           -- versioned metadata (LAYOUT_VERSION 1)
      known_bin_seq.bin   -- N * 8 bytes (uint64 2-bit-packed barcodes)
      index_ranges.bin    -- (4^k + 1) * 8 bytes (uint64 cumulative offsets)
      inverted_list.bin   -- index_size * 5 bytes (uint40 packed barcode IDs)

The build path mirrors :class:`SharedMemorySparseAnchorIndexer` (low-mem
two-pass anchor extraction) but writes directly into ``numpy.memmap``-backed
arrays so the inverted list never has to be resident as a single
contiguous allocation.
"""

import json
import logging
import math
import mmap
import os
import shutil
import tempfile
from typing import List, Optional, Tuple

import numpy

from ..common import bit_to_str, str_to_2bit
from .shared_memory import (
    INDEX_BYTES_PER_ENTRY,
    MAX_INDEXED_BARCODES,
    MAX_KMER_SIZE_FOR_UINT32_ANCHORS,
    NUM_ANCHORS,
    _compute_and_count_anchors,
    _compute_and_populate_anchors,
    _compute_anchors,
    _count_anchors,
    _hamming_2bit_python,
    _populate_anchors,
    _splitmix64_python,
)

logger = logging.getLogger('IsoQuant')


def _madvise_random(memmap_array: numpy.memmap) -> None:
    """Apply ``MADV_RANDOM`` to a numpy.memmap-backed array.

    Best-effort: silently no-ops on platforms or interpreters where
    ``mmap.mmap.madvise`` or ``MADV_RANDOM`` are unavailable. Only the
    underlying ``mmap.mmap`` is reachable via the private ``_mmap`` attribute
    of ``numpy.memmap``; that attribute has been stable since NumPy 1.x.
    """
    underlying = getattr(memmap_array, "_mmap", None)
    if underlying is None:
        return
    madvise = getattr(underlying, "madvise", None)
    advice = getattr(mmap, "MADV_RANDOM", None)
    if madvise is None or advice is None:
        return
    try:
        madvise(advice)
    except (OSError, ValueError):
        pass


def _madvise_sequential(memmap_array: numpy.memmap) -> None:
    """Apply ``MADV_SEQUENTIAL`` for the linear-scan build phase."""
    underlying = getattr(memmap_array, "_mmap", None)
    if underlying is None:
        return
    madvise = getattr(underlying, "madvise", None)
    advice = getattr(mmap, "MADV_SEQUENTIAL", None)
    if madvise is None or advice is None:
        return
    try:
        madvise(advice)
    except (OSError, ValueError):
        pass


class MappedIndexInfo:
    """
    Metadata container for handing a mmap-based index to worker processes.

    Workers reconstruct the indexer by re-opening the files in the named
    cache directory. The directory is owned by the main process; workers
    open files read-only and share the kernel's page cache.
    """

    def __init__(
        self,
        cache_dir: str,
        barcode_count: int,
        kmer_size: int,
        seq_len: int,
        index_size: int,
        layout_version: int,
    ):
        self.cache_dir = cache_dir
        self.barcode_count = barcode_count
        self.kmer_size = kmer_size
        self.seq_len = seq_len
        self.index_size = index_size
        self.layout_version = layout_version

    def __getstate__(self):
        return (
            self.cache_dir,
            self.barcode_count,
            self.kmer_size,
            self.seq_len,
            self.index_size,
            self.layout_version,
        )

    def __setstate__(self, state):
        (
            self.cache_dir,
            self.barcode_count,
            self.kmer_size,
            self.seq_len,
            self.index_size,
            self.layout_version,
        ) = state


class MappedSparseAnchorIndexer:
    """
    Disk-backed sparse 3-anchor indexer.

    Construction allocates three files under ``cache_dir`` and writes the
    anchor postings into them via ``numpy.memmap``. Worker processes
    reconstruct the indexer with :meth:`from_sharable_info` and open the
    same files read-only. Query semantics match
    :class:`SharedMemorySparseAnchorIndexer` exactly.

    The main-process instance owns the cache directory and removes it on
    deletion unless ``persistent=True``. Workers never delete files.
    """

    SEQ_DTYPE = numpy.uint64
    INDEX_BYTE_DTYPE = numpy.uint8
    RANGES_DTYPE = numpy.uint64

    LAYOUT_VERSION: int = 1
    META_FILENAME: str = "meta.json"
    BARCODES_FILENAME: str = "known_bin_seq.bin"
    RANGES_FILENAME: str = "index_ranges.bin"
    INDEX_FILENAME: str = "inverted_list.bin"

    def __init__(
        self,
        known_bin_seq,
        kmer_size: int = 14,
        seq_len: int = 25,
        *,
        cache_dir: Optional[str] = None,
        persistent: bool = False,
        low_mem: bool = True,
    ):
        """
        Build a mmap-backed sparse 3-anchor index.

        Args:
            known_bin_seq: Pre-encoded sequences as integers (use
                :func:`str_to_2bit`). Accepts a ``numpy.uint64`` array
                (preferred — used without copy) or any sequence convertible
                by ``numpy.asarray``.
            kmer_size: Length of k-mers (anchor length). Default 14.
            seq_len: Length of sequences. Default 25.
            cache_dir: Directory in which to allocate the index files. If
                ``None``, an empty directory is created via
                :func:`tempfile.mkdtemp` honoring ``$TMPDIR``. The caller
                is encouraged to pass an explicit fast-local-disk path
                (e.g., NVMe scratch) for large whitelists.
            persistent: If ``False`` (default), the cache directory is
                removed when the main-process instance is destroyed. If
                ``True``, the files survive — useful for build-once,
                query-many workflows.
            low_mem: If ``True`` (default), build using the fused two-pass
                anchor-extraction path (no transient ``(N, 3)`` uint32
                buffer). Set ``False`` only if profiling shows the extra
                CPU pass dominates and RAM headroom is available.
        """
        if seq_len < kmer_size + 2:
            raise ValueError(
                f"MappedSparseAnchorIndexer requires seq_len >= "
                f"kmer_size + 2 (got seq_len={seq_len}, kmer_size={kmer_size}); "
                f"need at least one internal window for the minimizer."
            )
        if kmer_size > MAX_KMER_SIZE_FOR_UINT32_ANCHORS:
            raise ValueError(
                f"MappedSparseAnchorIndexer stores anchor k-mer values "
                f"as uint32, so kmer_size must be <= "
                f"{MAX_KMER_SIZE_FOR_UINT32_ANCHORS} (got kmer_size={kmer_size})."
            )

        self.main_instance: bool = True
        self.persistent: bool = persistent
        self.k: int = kmer_size
        self.seq_len: int = seq_len
        self.mask: int = (1 << (2 * self.k)) - 1
        self.seq_mask: int = (1 << (2 * self.seq_len)) - 1
        self.total_sequences: int = len(known_bin_seq)
        total_kmers = int(math.pow(4, self.k))

        self._first_shift: int = (seq_len - kmer_size) * 2
        self._last_shift: int = 0
        self._mid_shifts_py: Tuple[int, ...] = tuple(
            (seq_len - kmer_size - p) * 2
            for p in range(1, seq_len - kmer_size)
        )

        self.cache_dir: str = (
            cache_dir
            if cache_dir is not None
            else tempfile.mkdtemp(prefix="barcode_index_")
        )
        os.makedirs(self.cache_dir, exist_ok=True)

        self._barcodes_path: str = os.path.join(
            self.cache_dir, self.BARCODES_FILENAME)
        self._ranges_path: str = os.path.join(
            self.cache_dir, self.RANGES_FILENAME)
        self._index_path: str = os.path.join(
            self.cache_dir, self.INDEX_FILENAME)
        self._meta_path: str = os.path.join(
            self.cache_dir, self.META_FILENAME)

        if self.total_sequences == 0:
            self.index_size: int = 0
            self.known_bin_seq = numpy.array([], dtype=self.SEQ_DTYPE)
            self.index = numpy.array([], dtype=self.INDEX_BYTE_DTYPE)
            self.index_ranges = numpy.zeros(
                total_kmers + 1, dtype=self.RANGES_DTYPE)
            self._persist_meta()
            return

        if self.total_sequences > MAX_INDEXED_BARCODES:
            raise ValueError(
                f"MappedSparseAnchorIndexer supports at most "
                f"{MAX_INDEXED_BARCODES} sequences (uint40 packed entries); "
                f"got {self.total_sequences}"
            )

        seqs = numpy.asarray(known_bin_seq, dtype=self.SEQ_DTYPE)

        first_shift = numpy.uint64(self._first_shift)
        last_shift = numpy.uint64(self._last_shift)
        k_mask = numpy.uint64(self.mask)
        mid_shifts = numpy.array(self._mid_shifts_py, dtype=numpy.uint64)

        logger.info(
            "Building mmap-backed sparse 3-anchor index for %d sequences "
            "(low_mem: %s, cache_dir: %s)",
            self.total_sequences, low_mem, self.cache_dir,
        )

        # --- Pass 1: count anchors per bucket (no large allocation) -------
        kmer_counts = numpy.zeros(total_kmers, dtype=numpy.uint64)
        if low_mem:
            _compute_and_count_anchors(
                seqs, first_shift, last_shift, k_mask, mid_shifts, kmer_counts,
            )
            anchors = None
            logger.info("Anchor counting complete (low-mem)")
        else:
            anchors = numpy.empty(
                (self.total_sequences, NUM_ANCHORS), dtype=numpy.uint32,
            )
            _compute_anchors(
                seqs, first_shift, last_shift, k_mask, mid_shifts, anchors,
            )
            logger.info("Anchor extraction complete")
            _count_anchors(anchors, kmer_counts)
            logger.info("Anchor counting complete")

        self.index_size = int(NUM_ANCHORS * self.total_sequences)

        # --- Allocate barcodes memmap and copy in ------------------------
        self.known_bin_seq = numpy.memmap(
            self._barcodes_path,
            dtype=self.SEQ_DTYPE,
            mode="w+",
            shape=(self.total_sequences,),
        )
        _madvise_sequential(self.known_bin_seq)
        self.known_bin_seq[:] = seqs[:]
        del seqs

        # --- Allocate index_ranges memmap and write cumsum ---------------
        self.index_ranges = numpy.memmap(
            self._ranges_path,
            dtype=self.RANGES_DTYPE,
            mode="w+",
            shape=(total_kmers + 1,),
        )
        _madvise_sequential(self.index_ranges)
        self.index_ranges[0] = 0
        self.index_ranges[1:] = numpy.cumsum(kmer_counts)
        del kmer_counts

        # --- Allocate inverted-list memmap as the populate target --------
        self.index = numpy.memmap(
            self._index_path,
            dtype=self.INDEX_BYTE_DTYPE,
            mode="w+",
            shape=(self.index_size * INDEX_BYTES_PER_ENTRY,),
        )
        # Random scatter writes follow.
        _madvise_random(self.index)

        if low_mem:
            _compute_and_populate_anchors(
                self.known_bin_seq, first_shift, last_shift, k_mask,
                mid_shifts, self.index_ranges, self.index,
            )
            logger.info(
                "Sparse anchor index population complete (low-mem, mmap)")
        else:
            _populate_anchors(anchors, self.index_ranges, self.index)
            del anchors
            logger.info("Sparse anchor index population complete (mmap)")

        # Flush dirty pages so workers opening the file read consistent data.
        self.known_bin_seq.flush()
        self.index_ranges.flush()
        self.index.flush()

        # Mark for query phase: read pattern is random.
        _madvise_random(self.known_bin_seq)
        _madvise_random(self.index_ranges)

        del mid_shifts

        self._persist_meta()

    # --- Lifecycle ------------------------------------------------------

    def __del__(self):
        # Release file handles by dropping memmap refs explicitly so the
        # directory removal below has no busy maps to fight (matters on
        # filesystems where unlink-while-mapped is awkward to debug).
        for attr in ("known_bin_seq", "index_ranges", "index"):
            buf = getattr(self, attr, None)
            if isinstance(buf, numpy.memmap):
                underlying = getattr(buf, "_mmap", None)
                if underlying is not None:
                    try:
                        underlying.close()
                    except (OSError, ValueError):
                        pass
            setattr(self, attr, None)

        if not getattr(self, "main_instance", False):
            return
        if getattr(self, "persistent", False):
            return
        cache_dir = getattr(self, "cache_dir", None)
        if cache_dir and os.path.isdir(cache_dir):
            try:
                shutil.rmtree(cache_dir)
            except OSError:
                pass

    # --- Metadata persistence -------------------------------------------

    def _persist_meta(self) -> None:
        meta = {
            "layout_version": self.LAYOUT_VERSION,
            "barcode_count": int(self.total_sequences),
            "kmer_size": int(self.k),
            "seq_len": int(self.seq_len),
            "index_size": int(self.index_size),
        }
        tmp_path = self._meta_path + ".tmp"
        with open(tmp_path, "w") as fh:
            json.dump(meta, fh)
        os.replace(tmp_path, self._meta_path)

    @classmethod
    def _load_meta(cls, cache_dir: str) -> dict:
        meta_path = os.path.join(cache_dir, cls.META_FILENAME)
        with open(meta_path, "r") as fh:
            meta = json.load(fh)
        version = int(meta.get("layout_version", -1))
        if version != cls.LAYOUT_VERSION:
            raise ValueError(
                f"MappedSparseAnchorIndexer cache at {cache_dir} has layout "
                f"version {version}, expected {cls.LAYOUT_VERSION}. Rebuild."
            )
        return meta

    # --- Cross-process handoff ------------------------------------------

    def get_sharable_info(self) -> MappedIndexInfo:
        return MappedIndexInfo(
            cache_dir=self.cache_dir,
            barcode_count=self.total_sequences,
            kmer_size=self.k,
            seq_len=self.seq_len,
            index_size=self.index_size,
            layout_version=self.LAYOUT_VERSION,
        )

    @classmethod
    def from_sharable_info(
        cls, shared_mem_index_info: MappedIndexInfo,
    ) -> "MappedSparseAnchorIndexer":
        if shared_mem_index_info.layout_version != cls.LAYOUT_VERSION:
            raise ValueError(
                f"MappedIndexInfo layout version "
                f"{shared_mem_index_info.layout_version} != "
                f"{cls.LAYOUT_VERSION}; rebuild the index."
            )

        kmer_index = cls.__new__(cls)
        kmer_index.main_instance = False
        kmer_index.persistent = True  # worker must not unlink files
        kmer_index.k = shared_mem_index_info.kmer_size
        kmer_index.mask = (1 << (2 * kmer_index.k)) - 1
        kmer_index.seq_len = shared_mem_index_info.seq_len
        kmer_index.seq_mask = (1 << (2 * kmer_index.seq_len)) - 1
        kmer_index.index_size = shared_mem_index_info.index_size
        kmer_index.total_sequences = shared_mem_index_info.barcode_count
        total_kmers = int(math.pow(4, kmer_index.k))

        kmer_index._first_shift = (kmer_index.seq_len - kmer_index.k) * 2
        kmer_index._last_shift = 0
        kmer_index._mid_shifts_py = tuple(
            (kmer_index.seq_len - kmer_index.k - p) * 2
            for p in range(1, kmer_index.seq_len - kmer_index.k)
        )

        kmer_index.cache_dir = shared_mem_index_info.cache_dir
        kmer_index._barcodes_path = os.path.join(
            kmer_index.cache_dir, cls.BARCODES_FILENAME)
        kmer_index._ranges_path = os.path.join(
            kmer_index.cache_dir, cls.RANGES_FILENAME)
        kmer_index._index_path = os.path.join(
            kmer_index.cache_dir, cls.INDEX_FILENAME)
        kmer_index._meta_path = os.path.join(
            kmer_index.cache_dir, cls.META_FILENAME)

        if kmer_index.total_sequences == 0:
            kmer_index.known_bin_seq = numpy.array([], dtype=cls.SEQ_DTYPE)
            kmer_index.index = numpy.array([], dtype=cls.INDEX_BYTE_DTYPE)
            kmer_index.index_ranges = numpy.zeros(
                total_kmers + 1, dtype=cls.RANGES_DTYPE)
            return kmer_index

        kmer_index.known_bin_seq = numpy.memmap(
            kmer_index._barcodes_path,
            dtype=cls.SEQ_DTYPE,
            mode="r",
            shape=(kmer_index.total_sequences,),
        )
        kmer_index.index_ranges = numpy.memmap(
            kmer_index._ranges_path,
            dtype=cls.RANGES_DTYPE,
            mode="r",
            shape=(total_kmers + 1,),
        )
        kmer_index.index = numpy.memmap(
            kmer_index._index_path,
            dtype=cls.INDEX_BYTE_DTYPE,
            mode="r",
            shape=(kmer_index.index_size * INDEX_BYTES_PER_ENTRY,),
        )
        _madvise_random(kmer_index.known_bin_seq)
        _madvise_random(kmer_index.index_ranges)
        _madvise_random(kmer_index.index)
        return kmer_index

    # --- Public API -----------------------------------------------------

    def empty(self) -> bool:
        return len(self.known_bin_seq) == 0

    def _query_anchors(self, seq: int) -> Tuple[int, int, int]:
        first = (seq >> self._first_shift) & self.mask
        last = (seq >> self._last_shift) & self.mask
        best_hash = 0xFFFFFFFFFFFFFFFF
        best_kmer = 0
        for shift in self._mid_shifts_py:
            kmer = (seq >> shift) & self.mask
            h = _splitmix64_python(kmer)
            if h < best_hash:
                best_hash = h
                best_kmer = kmer
        return first, last, best_kmer

    def get_occurrences(
        self,
        sequence: str,
        max_hits: int = 0,
        min_kmers: int = 1,
        hits_delta: int = 1,
        ignore_equal: bool = False,
    ) -> List[Tuple[str, int, List[int]]]:
        """Sparse-anchor query — semantics identical to the SHM variant."""
        if len(self.known_bin_seq) == 0:
            return []

        seq = str_to_2bit(sequence)
        if seq > self.seq_mask:
            return []

        q_first, q_last, q_min = self._query_anchors(seq)

        slot_masks: dict[int, int] = {}
        index_buf = self.index
        for slot_bit, hash_val in (
                (0b001, q_first), (0b010, q_last), (0b100, q_min)):
            start_index = int(self.index_ranges[hash_val])
            end_index = int(self.index_ranges[hash_val + 1])
            for entry_idx in range(start_index, end_index):
                base = entry_idx * INDEX_BYTES_PER_ENTRY
                seq_idx = (int(index_buf[base])
                           | (int(index_buf[base + 1]) << 8)
                           | (int(index_buf[base + 2]) << 16)
                           | (int(index_buf[base + 3]) << 24)
                           | (int(index_buf[base + 4]) << 32))
                barcode = int(self.known_bin_seq[seq_idx])
                slot_masks[barcode] = slot_masks.get(barcode, 0) | slot_bit

        if not slot_masks:
            return []

        tier_a: List[Tuple[int, int, int]] = []
        tier_b: List[Tuple[int, int, int, int]] = []
        tier_c: List[Tuple[int, int, int, int]] = []
        first_pos = 0
        last_pos = self.seq_len - self.k
        for barcode, mask in slot_masks.items():
            if ignore_equal and barcode == seq:
                continue
            anchor_count = (mask & 1) + ((mask >> 1) & 1) + ((mask >> 2) & 1)
            if anchor_count < min_kmers:
                continue
            if (mask & 0b011) == 0b011:
                tier_a.append((barcode, anchor_count, mask))
            elif anchor_count == 2:
                ham = _hamming_2bit_python(seq, barcode)
                tier_b.append((barcode, anchor_count, mask, ham))
            else:
                ham = _hamming_2bit_python(seq, barcode)
                tier_c.append((barcode, anchor_count, mask, ham))

        tier_a.sort(key=lambda x: -x[1])
        tier_b.sort(key=lambda x: x[3])
        tier_c.sort(key=lambda x: x[3])

        result: List[Tuple[str, int, List[int]]] = []

        def _emit(barcode: int, anchor_count: int, mask: int, ham: int) -> None:
            score = 1000 * anchor_count - ham
            positions: List[int] = []
            if mask & 0b001:
                positions.append(first_pos)
            if mask & 0b010:
                positions.append(last_pos)
            if mask & 0b100:
                positions.append(-1)
            result.append((bit_to_str(barcode, self.seq_len), score, positions))

        for barcode, anchor_count, mask in tier_a:
            _emit(barcode, anchor_count, mask, 0)
            if max_hits and len(result) >= max_hits:
                return result
        for barcode, anchor_count, mask, ham in tier_b:
            _emit(barcode, anchor_count, mask, ham)
            if max_hits and len(result) >= max_hits:
                return result
        for barcode, anchor_count, mask, ham in tier_c:
            _emit(barcode, anchor_count, mask, ham)
            if max_hits and len(result) >= max_hits:
                return result

        return result
