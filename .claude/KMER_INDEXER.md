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
- Flat array structure for cache-friendly access
- Best for large barcode sets (e.g., 10x whitelists)
- Lowest memory usage

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

## Performance

**Barcode calling** (16bp barcodes, 10K whitelist):
- `KmerIndexer` (k=6): ~5K queries/sec
- `ArrayKmerIndexer` (k=6): ~12K queries/sec
- `Array2BitKmerIndexer` (k=12): ~15K queries/sec, 50% less memory
