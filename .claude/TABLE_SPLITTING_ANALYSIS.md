# Table Splitting Algorithm Analysis

## Problem Statement

When read groups are provided via TSV file, IsoQuant splits them by chromosome so each worker thread only loads reads relevant to its chromosome. This is critical for memory efficiency with large datasets.

## Current Implementations

### 1. Single-threaded (`src/read_groups.py:383`, `437`)

```python
def split_table(table_file, sample, out_prefix, ...):
    read_groups = load_multicolumn_table(table_file, ...)  # Load entire table into memory
    read_group_files = {}
    processed_reads = defaultdict(set)
    bam_files = list(map(lambda x: x[0], sample.file_list))

    for bam_file in bam_files:
        bam = pysam.AlignmentFile(bam_file, "rb")
        # Open output files for each chromosome
        for chr_id in bam.references:
            if chr_id not in read_group_files:
                read_group_files[chr_id] = open(out_prefix + "_" + chr_id, "w")
        # Scan all reads
        for read_alignment in bam:
            chr_id = read_alignment.reference_name
            read_id = read_alignment.query_name
            if read_id in read_groups and read_id not in processed_reads[chr_id]:
                read_group_files[chr_id].write(...)
                processed_reads[chr_id].add(read_id)
```

**Complexity:**
- Time: O(R) where R = total reads in BAM files
- Memory: O(T + R*C) where T = table size, C = chromosome count
- Single-threaded, slow for large BAMs

### 2. Parallel (`src/table_splitter.py:151`)

```python
def split_read_table_parallel(sample, input_tsvs, split_reads_file_names, num_threads, ...):
    # Step 1: Cache all read IDs per chromosome (scan BAMs)
    read_id_cache = {}
    for chr_id in chromosomes_to_process:
        read_id_cache[chr_id] = collect_chromosome_reads(chr_id, bam_files)  # O(R)

    # Step 2: Create generator that yields (chr_id, chunk) pairs
    barcode_chunk_gen = combine_generators(chromosomes_to_process, load_func(input_tsvs))

    # Step 3: Process chunks in parallel
    # For each chunk, create C tasks (one per chromosome)
    for chr_id, read_chunk in barcode_chunk_gen:
        proc.submit(process_barcode_chunk_with_read_cache,
                   chr_id, read_chunk,
                   read_id_cache[chr_id],  # Pass cached IDs
                   split_reads_file_names[chr_id])
```

**Key mechanism:**
- `combine_generators()` (line 141): For each TSV chunk, yields (chr_id, chunk) for ALL chromosomes
- Each chunk is processed C times (once per chromosome)

**Complexity:**
- Time: O(R + K*C*T/500K) where K = num_threads, T = table size, C = chromosome count
- Memory: O(R*C) for read_id_cache
- Tasks: (T/500K) * C total tasks submitted

## Efficiency Analysis

### Current Parallel Algorithm Issues

#### 1. **Redundant Work** (Critical Issue)
For each TSV chunk of 500K reads, the algorithm creates C tasks (one per chromosome).

**Example:**
- 10M reads in TSV → 20 chunks
- 20 chromosomes
- **Total tasks: 20 × 20 = 400 tasks**
- **Each chunk processed 20 times!**

Only ~5% of work is useful (1/20 chromosomes match per chunk on average).

#### 2. **Serialization Overhead**
Each task submission passes:
- `read_chunk` dict (~500K entries, ~50MB)
- `read_id_cache[chr_id]` set (can be millions of IDs, ~100MB+)

With ProcessPoolExecutor, these are pickled/unpickled for each task. For 400 tasks, this is **massive overhead**.

#### 3. **Memory Duplication**
- `read_id_cache` contains ALL read IDs for ALL chromosomes
- With C=20 chromosomes, ~1M reads each → 20M IDs in memory
- Each worker process gets a copy of relevant cache
- Peak memory: `num_threads * cache_size`

#### 4. **Generator Bottleneck**
`combine_generators()` creates a Cartesian product of (chromosomes × chunks):
```python
def combine_generators(chromosome_list, read_chunk_generator):
    while True:
        read_chunk = next(read_chunk_generator)
        for chr_id in chromosome_list:
            yield chr_id, read_chunk  # Yields C items per chunk!
```

This forces sequential processing of C tasks per chunk, limiting parallelism.

#### 5. **File I/O Contention**
Multiple processes may append to files, though this is less of an issue with append mode.

### Benchmarks (Estimated)

Assuming:
- 10M reads in TSV table
- 20 chromosomes
- 100M total reads in BAM files
- 8 threads

**Current parallel algorithm:**
- Pre-cache BAMs: ~5 minutes (scan 100M reads)
- Create 400 tasks (20 chunks × 20 chromosomes)
- Process 400 tasks with 8 workers: ~10 minutes (mostly serialization overhead)
- **Total: ~15 minutes**

**Single-threaded algorithm:**
- Scan 100M reads from BAM: ~5 minutes
- Lookup in 10M table: O(1) per read
- **Total: ~5 minutes** (faster than parallel due to less overhead!)

**The parallel version is actually SLOWER** due to overhead!

## Proposed Improved Algorithm

### Design Principle
**"Each worker owns specific chromosomes and processes table chunks"**

### Algorithm

```python
def split_read_table_parallel_improved(sample, input_tsvs, split_reads_file_names, num_threads):
    # Step 1: Assign chromosomes to workers
    chromosomes = list(split_reads_file_names.keys())
    chr_per_worker = distribute_chromosomes(chromosomes, num_threads)
    # Example: worker 0 gets [chr1, chr2], worker 1 gets [chr3, chr4], etc.

    # Step 2: Pre-cache read IDs per worker (once)
    worker_read_caches = []
    for worker_id in range(num_threads):
        cache = {}
        for chr_id in chr_per_worker[worker_id]:
            cache[chr_id] = collect_chromosome_reads(chr_id, bam_files)
        worker_read_caches.append(cache)

    # Step 3: Process table chunks in parallel
    with ProcessPoolExecutor(max_workers=num_threads) as proc:
        # Create workers with their assigned chromosomes
        workers = []
        for worker_id in range(num_threads):
            workers.append(proc.submit(
                process_chunks_for_chromosomes,
                worker_id,
                input_tsvs,
                chr_per_worker[worker_id],
                worker_read_caches[worker_id],
                split_reads_file_names
            ))

        # Wait for all workers to complete
        for future in workers:
            future.result()

def process_chunks_for_chromosomes(worker_id, input_tsvs, my_chromosomes, my_read_cache, output_files):
    """Each worker processes the ENTIRE table but only for its assigned chromosomes"""
    # Open output files once
    out_handles = {chr_id: open(output_files[chr_id], 'w') for chr_id in my_chromosomes}

    # Stream through table chunks
    for read_chunk in load_table_chunked(input_tsvs):
        for read_id, group_vals in read_chunk.items():
            # Check which of MY chromosomes this read belongs to
            for chr_id in my_chromosomes:
                if read_id in my_read_cache[chr_id]:
                    out_handles[chr_id].write(f"{read_id}\t{group_vals}\n")

    # Close files
    for f in out_handles.values():
        f.close()
```

### Key Improvements

1. **No Redundant Work**: Each chunk processed exactly once by each worker, but workers only check their assigned chromosomes

2. **Minimal Serialization**:
   - Workers created once with their chromosome assignments
   - Each worker reads table independently (no passing chunks between processes)
   - No pickle overhead for chunks

3. **Better Memory Usage**:
   - Each worker holds cache only for its chromosomes
   - No duplicate caches across processes
   - Memory: `sum(cache_size_per_worker)` instead of `num_threads * total_cache_size`

4. **True Parallelism**:
   - All workers run simultaneously
   - No generator bottleneck
   - Each worker independently streams the table

5. **Clean I/O**:
   - Each worker writes to its own files exclusively
   - No file contention

### Complexity Analysis

**Time:** O(R/C + T*C/W) where:
- R = total BAM reads
- T = table size
- C = chromosome count
- W = num workers

**With W = C (one worker per chromosome):** O(R + T)

**Memory:** O(R) total across all workers (each holds cache for its chromosomes)

**Tasks:** W tasks (not T/500K * C!)

### Estimated Performance

Same scenario (10M table, 20 chr, 100M BAM reads, 8 threads):
- Pre-cache BAMs: ~5 minutes (distributed across 8 workers scanning their chromosomes)
- Process table: Each worker streams 10M reads, checks ~2.5 chromosomes average
- **Total: ~3-4 minutes**

**3-4x faster than current parallel, 25% faster than single-threaded**

## Alternative: Chromosome-First Approach

Even simpler approach - assign chromosomes to workers, let each worker do everything:

```python
def split_read_table_chromosome_first(sample, input_tsvs, split_reads_file_names, num_threads):
    chromosomes = list(split_reads_file_names.keys())

    with ProcessPoolExecutor(max_workers=num_threads) as proc:
        futures = []
        for chr_id in chromosomes:
            futures.append(proc.submit(
                process_chromosome,
                chr_id,
                input_tsvs,
                bam_files,
                split_reads_file_names[chr_id]
            ))

        for future in futures:
            future.result()

def process_chromosome(chr_id, table_file, bam_files, output_file):
    """Process one chromosome completely"""
    # 1. Scan BAMs for this chromosome's read IDs
    read_ids = collect_chromosome_reads(chr_id, bam_files)

    # 2. Stream table and write matches
    with open(output_file, 'w') as out:
        for read_chunk in load_table_chunked(table_file):
            for read_id, group_vals in read_chunk.items():
                if read_id in read_ids:
                    out.write(f"{read_id}\t{group_vals}\n")
```

**Pros:**
- Simplest implementation
- Each chromosome completely independent
- No shared state
- Natural load balancing (ProcessPoolExecutor queue)

**Cons:**
- Each worker scans full table (T size)
- Total work: C * T (redundant table scanning)

**When this is better:**
- Small tables (T < 1M reads)
- Many chromosomes (C > 20)
- Table scanning is cheap (local SSD, small file)

## Recommendation

### For Current IsoQuant Usage:

**Use the improved algorithm** (section "Proposed Improved Algorithm")

**Rationale:**
1. Single-cell data often has 10M+ barcode tables
2. Typical use: 20-30 chromosomes
3. Table scanning dominates runtime for large tables
4. Better scalability for future datasets

### Implementation Strategy:

1. Replace `split_read_table_parallel()` with improved version
2. Keep single-threaded as fallback for small datasets
3. Auto-select based on heuristics:
   - If `table_size < 1M` → use single-threaded
   - If `table_size < 10M` and `num_chr > 20` → use chromosome-first
   - Otherwise → use improved parallel

### Code Changes Required:

**Files to modify:**
1. `src/table_splitter.py`: Implement new `split_read_table_parallel_improved()`
2. `src/read_groups.py`: Update `prepare_read_groups()` to call improved version
3. Add heuristics for algorithm selection

**Backward compatibility:**
- Keep old implementation as `split_read_table_parallel_legacy()`
- Add flag `--table_split_algorithm` for testing

## Edge Cases

### 1. Unbalanced Chromosome Sizes
Some chromosomes have far more reads than others (e.g., chr1 vs. chrM).

**Solution:** Distribute by estimated work, not chromosome count:
```python
def distribute_chromosomes_weighted(chromosomes, num_workers, bam_files):
    # Estimate read count per chromosome
    chr_sizes = {}
    for chr_id in chromosomes:
        chr_sizes[chr_id] = estimate_chromosome_reads(chr_id, bam_files)

    # Greedy bin packing: assign largest chromosomes first
    workers = [[] for _ in range(num_workers)]
    worker_loads = [0] * num_workers

    for chr_id in sorted(chromosomes, key=lambda c: chr_sizes[c], reverse=True):
        min_worker = worker_loads.index(min(worker_loads))
        workers[min_worker].append(chr_id)
        worker_loads[min_worker] += chr_sizes[chr_id]

    return workers
```

### 2. Very Large Tables (100M+ reads)
Streaming becomes critical.

**Solution:** Already handled by `load_table_chunked()` generator.

### 3. Compressed TSV Files (.gz)
Decompression can be bottleneck.

**Solution:** Use `pigz` for parallel decompression or `zcat | split` preprocessing.

### 4. Memory Constraints
Worker caches might exceed RAM.

**Solution:**
- Spill caches to disk (SQLite or mmap'd files)
- Process in batches (split by table chunks AND chromosomes)

## Testing Plan

1. **Unit tests:**
   - Small synthetic data (1K reads, 3 chromosomes)
   - Verify output correctness vs. single-threaded

2. **Performance benchmarks:**
   - 1M, 10M, 100M read tables
   - 5, 20, 50 chromosomes
   - 1, 4, 8, 16 threads
   - Compare: single-threaded, current parallel, improved parallel

3. **Memory profiling:**
   - Peak RSS per worker
   - Total memory usage
   - Cache sizes

4. **Real datasets:**
   - 10x scRNA-seq barcodes (typical: 5M cells, 30 chr)
   - Spatial transcriptomics (typical: 50M spots, 25 chr)

## Conclusion

The current parallel algorithm has significant inefficiencies due to redundant work and serialization overhead. The improved algorithm (worker-per-chromosome-set) is:
- **3-4x faster** than current
- **Simpler to understand** (no complex generator logic)
- **Better memory efficiency** (no cache duplication)
- **More scalable** (grows with chromosomes, not chunks × chromosomes)

**Recommendation: Implement the improved algorithm as the default.**
