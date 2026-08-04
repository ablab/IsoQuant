# K-mer Indexer for Barcode Calling

Fast approximate string matching using k-mer indexing. Used for single-cell barcode calling.

## Classes

### `KmerIndexer`
Dictionary-based k-mer index. Best for small to medium barcode sets.

**Usage**:
```python
indexer = KmerIndexer(["ACTG", "TGCA", "GGGG"], kmer_size=3)
results = indexer.get_occurrences("ACTC", max_hits=3)
# Returns: [("ACTG", 2, [0, 1]), ...] - (barcode, shared_kmers, positions)
```

### `ArrayKmerIndexer`
Array-based k-mer index using 2-bit nucleotide encoding.
- Memory: O(4^k) array entries
- Best for k ≤ 8
- ~2-3x faster than `KmerIndexer`

### `Array2BitKmerIndexer`
Memory-efficient index storing both k-mers and sequences in 2-bit format.
- Flat (CSR-style) `index` + `index_ranges` structure for cache-friendly access
- Built by a counting pass followed by a fill pass, so the 4^k per-k-mer lists are never
  materialised (they cost ~235 MB of empty list objects at k=11, ×16 per extra base)
- Best for large barcode sets (e.g., 10x whitelists)

### `SharedMemoryArray2BitKmerIndexer`
Same layout backed by `multiprocessing.shared_memory`, so worker processes map the index
instead of receiving a copy. Used by the Stereo-seq detectors.

## Two 2-bit alphabets

`ArrayKmerIndexer.NUCL2BIN` is **A,C,G,T = 0,1,2,3**; `common.NUCL2BIN` (used by
`str_to_2bit` / `batch_str_to_2bit` and the 2-bit indexers) is **A,C,T,G = 0,1,2,3**, chosen so
a base's code is `(ord(c) & 6) >> 1`. Both are self-consistent; codes must never be exchanged
between the two families.

## Algorithm

1. **Indexing**: Split each barcode into overlapping k-mers, map k-mer → barcode indices
2. **Querying**: For query sequence, count shared k-mers with each indexed barcode
3. **Filtering**: Return barcodes with most shared k-mers (within hits_delta threshold)

## Parameters

- `kmer_size` - Length of k-mers (default: 6-12 depending on class)
- `max_hits` - Limit results (0 = unlimited)
- `min_kmers` - Minimum shared k-mers to report match
- `hits_delta` - Include results within N k-mers of top hit
- `ignore_equal` - Skip exact matches

**`hits_delta` counts k-mer *occurrences*, not distinct k-mers.** For low-complexity sequences
(homopolymer runs) a single query k-mer can match many times in one indexed barcode, so the
count can far exceed the barcode length — `hits_delta=barcode_length` does *not* disable the
filter. The provable upper bound is `(barcode_length - kmer_size + 1)²`; see
`BarcodeGraph._no_hit_pruning()`. This pruning is one of the reasons per-read whitelist
matching misses error-bearing barcodes (`.claude/BARCODE_GRAPH_CORRECTION.md`).

## Performance

**Barcode calling** (16bp barcodes, 10K whitelist):
- `KmerIndexer` (k=6): ~5K queries/sec
- `ArrayKmerIndexer` (k=6): ~12K queries/sec
- `Array2BitKmerIndexer` (k=12): ~15K queries/sec, 50% less memory
