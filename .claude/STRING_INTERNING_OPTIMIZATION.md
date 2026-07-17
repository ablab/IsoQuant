# String Interning Optimization for Memory Efficiency

## Problem Analysis

### Current Memory Usage

With **1 billion reads**, the current string storage in ReadAssignment and IsoformMatch causes massive memory duplication:

**Per ReadAssignment object** (~400-800 bytes of string data):
- `chr_id`: ~10 bytes (e.g., "chr1", "chrX")
- `read_group`: list of strings, e.g., ["CB_ATCG", "file_sample1"] = ~30-50 bytes
- `barcode`: ~16-20 bytes (cell barcode)
- `umi`: ~10-12 bytes (UMI sequence)
- `strand`, `mapped_strand`: 1-2 bytes each

**Per IsoformMatch object** (1-10 per read, ~50-100 bytes of string data each):
- `assigned_gene`: ~15-30 bytes (e.g., "ENSG00000123456")
- `assigned_transcript`: ~15-30 bytes (e.g., "ENST00000654321")
- `transcript_strand`: 1 byte

**Memory Impact Example**:
- 1 billion reads
- Average 3 IsoformMatches per read
- ~600 bytes of duplicated string data per read
- **Total: ~600 GB just for strings** (most of which are duplicates!)

### String Cardinality Analysis

Strings fall into different categories based on cardinality:

1. **Tiny finite set** (< 10 values):
   - `strand`, `mapped_strand`: {'+', '-', '.'}
   - `transcript_strand`: {'+', '-', '.'}

2. **Small finite set** (< 100 values):
   - `chr_id`: ~25-100 chromosomes per organism

3. **Medium finite set** (known in advance):
   - `assigned_gene`: ~20,000-60,000 genes (from GTF)
   - `assigned_transcript`: ~100,000-300,000 transcripts (from GTF)

4. **Large finite set** (known in advance):
   - `barcode`: ~10,000-100,000 barcodes (from barcode calling, split per chromosome)
   - `read_group` from TSV files: ~1,000s-100,000s values (split per chromosome)
   - `read_group` from file_name: ~10-100 filenames (known globally)
   - `read_group` from barcode_spot: ~10-1,000 spots/cell types (known globally)

5. **Large finite set** (discovered during processing):
   - `read_group` values from BAM tags (tag:CB, tag:UB, etc.)
   - `read_group` values from read ID suffixes (read_id:DELIM)
   - Could be 1,000s to 100,000s of unique values

6. **Very large or infinite set**:
   - `umi`: Could be millions of unique UMIs
   - `read_id`: Unique per read (cannot optimize)

### Serialization Impact

Current serialization uses:
- `write_string()` / `read_string()`: Stores full string with length prefix
- Each duplicate string is written completely to disk
- Disk usage: Same ~600 GB problem for temporary files

## Optimization Strategy

### High-Level Approach: String Interning with Integer References

Replace strings with integer indices that reference shared string pools (dictionaries).

**Core concept**: Store each unique string once, reference it many times by index.

### Architecture Design

#### 1. String Pool Manager

Create a centralized `StringPoolManager` that maintains multiple string pools:

```python
class StringPoolManager:
    def __init__(self):
        # Global known-in-advance pools (same across all chromosomes)
        self.gene_pool = StringPool()          # gene_id -> int (from GTF)
        self.transcript_pool = StringPool()    # transcript_id -> int (from GTF)
        self.chromosome_pool = StringPool()    # chr_id -> int (from GTF/reference)
        self.strand_pool = StringPool()        # strand values -> int (tiny: +/-/.)

        # Per-chromosome known-in-advance pools (loaded from split files)
        self.barcode_pool = StringPool()       # barcode -> int (from split barcode file)
        self.read_group_tsv_pool = StringPool()  # read_group from TSV -> int (from split read_group file)

        # Global known-in-advance read_group pools (for non-TSV grouping types)
        self.file_name_pool = StringPool()     # file_name -> int (from sample file list)
        self.barcode_spot_pool = StringPool()  # spot/cell_type -> int (from barcode2spot file)

        # Discovered-during-processing pools (only for tag/read_id grouping)
        self.read_group_dynamic_pool = StringPool()  # tag/read_id values -> int

        # UMI pool (optional, might be too large)
        self.umi_pool = StringPool()           # umi -> int (optional)
```

#### 2. StringPool Implementation

```python
class StringPool:
    """Bidirectional mapping between strings and integers"""
    def __init__(self, known_values=None):
        self.str_to_int = {}   # string -> int
        self.int_to_str = []   # list[string], index = int
        self.next_id = 0

        if known_values:
            # Pre-populate with known values
            for value in known_values:
                self.add(value)

    def add(self, string):
        """Add string to pool, return its integer ID"""
        if string in self.str_to_int:
            return self.str_to_int[string]

        int_id = self.next_id
        self.str_to_int[string] = int_id
        self.int_to_str.append(string)
        self.next_id += 1
        return int_id

    def get_int(self, string):
        """Get integer ID for string (adds if not present)"""
        return self.add(string)

    def get_str(self, int_id):
        """Get string for integer ID"""
        return self.int_to_str[int_id]
```

#### 3. Modified Data Classes

**IsoformMatch**:
```python
class IsoformMatch:
    def __init__(self, match_classification, assigned_gene=None,
                 assigned_transcript=None, transcript_strand='.',
                 penalty_score=0, string_pools=None):
        # Store integers instead of strings
        self.assigned_gene_id = None      # int
        self.assigned_transcript_id = None  # int
        self.transcript_strand_id = None  # int

        if string_pools:
            if assigned_gene:
                self.assigned_gene_id = string_pools.gene_pool.get_int(assigned_gene)
            if assigned_transcript:
                self.assigned_transcript_id = string_pools.transcript_pool.get_int(assigned_transcript)
            self.transcript_strand_id = string_pools.strand_pool.get_int(transcript_strand)

        # ... rest unchanged

    # Properties for backward compatibility
    @property
    def assigned_gene(self):
        if self.assigned_gene_id is None:
            return None
        return self._string_pools.gene_pool.get_str(self.assigned_gene_id)

    @property
    def assigned_transcript(self):
        if self.assigned_transcript_id is None:
            return None
        return self._string_pools.transcript_pool.get_str(self.assigned_transcript_id)
```

**ReadAssignment**:
```python
class ReadAssignment:
    def __init__(self, read_id, assignment_type, match=None, string_pools=None):
        self.read_id = read_id  # Keep as string (unique)

        # Replace strings with integers
        self.chr_id_int = None           # int
        self.strand_id = None            # int
        self.mapped_strand_id = None     # int
        self.barcode_id = None           # int
        self.umi_id = None               # int (optional)
        self.read_group_ids = []         # list[int]

        # ... rest unchanged
```

### Parallel Processing Strategy

**Challenge**: Each worker process in ProcessPoolExecutor needs access to string pools.

**Solution Options**:

#### Option A: Per-Worker Pools with Pre-Loading (RECOMMENDED)

1. **Before parallel processing** (main process):
   - Load global pools from external sources:
     - Gene/transcript pools from GTF
     - Chromosome pool from GTF/reference
     - File_name pool from sample file list
     - Barcode_spot pool from barcode2spot files
     - Strand pool with {'+', '-', '.'}
   - These pools are read-only and identical across all workers

2. **During parallel read collection** (per chromosome worker):
   - Each worker receives a copy of global pools (via pickling)
   - Load chromosome-specific pools from split files:
     - Barcode pool from `{sample.barcodes_split_reads}_{chr_id}`
     - Read_group_tsv pool from `{sample.read_group_file}_{chr_id}` (if TSV grouping used)
   - Build dynamic pool only if needed:
     - Read_group_dynamic pool for tag/read_id grouping (discovered during processing)
   - All pools except read_group_dynamic are **read-only** during processing

3. **Serialization** (per chromosome):
   - Serialize only the pools needed for this chromosome:
     - Global pools (gene, transcript, chromosome, strand, file_name, barcode_spot)
     - Chromosome-specific pools (barcode, read_group_tsv)
     - Dynamic pool (read_group_dynamic, if applicable)
   - Format: `[string_pool_header][pool_data][read_assignments]`
   - Each chromosome file is self-contained

4. **Deserialization** (when loading):
   - Load string pools first
   - Use them to reconstruct ReadAssignment objects

5. **Pool merging** (rarely needed):
   - Only read_group_dynamic pools differ across chromosomes
   - If cross-chromosome operations needed, merge dynamic pools and remap IDs
   - Most operations are per-chromosome, so merging rarely needed

#### Option B: Shared Pools via Manager (More Complex)

Use multiprocessing.Manager to share pools across processes, but:
- ❌ Slower due to IPC overhead
- ❌ More complex synchronization
- ❌ Not compatible with ProcessPoolExecutor's pickling

**Recommendation**: Use Option A.

### Read Group Pool Strategy by Type

Different `--read_group` types have different pool strategies:

| Group Type | Example | Pool Type | When Loaded | Scope |
|------------|---------|-----------|-------------|-------|
| `file:TSV` | `file:groups.tsv:0:1` | `read_group_tsv_pool` | Per-chromosome (from split TSV) | Read-only |
| `file_name` | `file_name` | `file_name_pool` | Global (from file list) | Read-only |
| `barcode_spot` | `barcode_spot` | `barcode_spot_pool` | Global (from barcode2spot) | Read-only |
| `tag:TAG` | `tag:CB` | `read_group_dynamic_pool` | During processing | Dynamic |
| `read_id:DELIM` | `read_id:_` | `read_group_dynamic_pool` | During processing | Dynamic |

**Benefits of this approach**:
1. **Maximize pre-loading**: 3 out of 5 group types use read-only pools
2. **Consistent with barcode handling**: Split files → pre-loaded pools
3. **Minimal dynamic building**: Only tag/read_id grouping needs dynamic pools
4. **Simple implementation**: ReadTableGrouper already loads the TSV, just extract unique values
5. **Memory efficient**: Pre-loaded pools shared across workers (via pickle), no per-worker overhead

**Implementation detail**:
When loading split read_group TSV files in `ReadTableGrouper.__init__()`:
```python
# Current code loads read_id -> group_value mapping
self.read_map = load_table(table_tsv_file, read_id_column_index,
                           group_id_column_index, delim)

# Enhancement: also extract unique group values for pool
unique_groups = set(self.read_map.values())
for group in unique_groups:
    string_pools.read_group_tsv_pool.add(group)
```

### Implementation Phases

#### Phase 1: Infrastructure (No Behavior Change)

1. Create `StringPool` and `StringPoolManager` classes
2. Add `string_pools` parameter to IsoformMatch and ReadAssignment constructors
3. Add properties for backward compatibility (`.assigned_gene` returns string from pool)
4. Update serialization to handle both old (string) and new (int + pool) formats
5. Add pool serialization/deserialization

**Goal**: Code works with both string and interned modes, controlled by flag.

#### Phase 2: Pre-Known Pools

**Global pools** (loaded once, shared across all workers):
1. Load gene/transcript pools from GTF during initialization
2. Load chromosome pool from reference/GTF
3. Pre-populate strand pool with {'+', '-', '.'}
4. Load file_name pool from sample file lists
5. Load barcode_spot pool from barcode2spot files (if provided)

**Per-chromosome pools** (loaded by each worker for its chromosome):
6. Load barcode pool from split barcode file: `{sample.barcodes_split_reads}_{chr_id}`
7. Load read_group_tsv pool from split read_group file: `{sample.read_group_file}_{chr_id}` (if TSV grouping)

**Goal**: All annotation-based and file-based strings are interned with read-only pools.

#### Phase 3: Dynamic Pools (Only for tag/read_id grouping)

1. Update AlignmentTagReadGrouper to use read_group_dynamic_pool
2. Update ReadIdSplitReadGrouper to use read_group_dynamic_pool
3. Ensure proper serialization of dynamic pools
4. Decide on UMI pool strategy (intern vs keep as string, measure cardinality first)

**Goal**: Tag/read_id based grouping strings are interned via dynamic pools.

**Note**: With Phase 2 complete, most common use cases (TSV grouping, file_name grouping, barcode_spot grouping) are already fully optimized with read-only pools. Phase 3 only adds value for less common tag/read_id grouping scenarios.

#### Phase 4: Serialization Optimization

1. Modify serialization format:
   ```
   [MAGIC_HEADER]
   [VERSION]
   [STRING_POOLS_HEADER]
     [gene_pool: int_count, list[string]]
     [transcript_pool: int_count, list[string]]
     [chr_pool: int_count, list[string]]
     [read_group_pool: int_count, list[string]]
     ...
   [READ_ASSIGNMENTS]
     [assignment_id, read_id_str, chr_id_int, gene_ids_int[], ...]
   ```

2. Update BasicReadAssignment to use integer lists
3. Optimize pool storage (compress, deduplicate)

**Goal**: Disk usage reduced by ~80-90%.

### Memory Savings Estimation

**Before optimization** (1B reads, 3 IsoformMatches avg):
- Strings: ~600 GB
- Other data: ~100 GB
- **Total: ~700 GB**

**After optimization**:

String pools:
- 50K genes × 25 bytes = 1.25 MB
- 250K transcripts × 25 bytes = 6.25 MB
- 50 chromosomes × 10 bytes = 500 bytes
- 50K barcodes × 18 bytes = 900 KB
- 10K read_groups × 30 bytes = 300 KB
- **Pool total: ~10 MB** (negligible!)

Integer references:
- 1B reads × 4 bytes (chr_id) = 4 GB
- 1B reads × 4 bytes (barcode_id) = 4 GB
- 1B reads × 2 ints × 4 bytes (strand IDs) = 8 GB
- 3B IsoformMatches × 2 ints × 4 bytes (gene/transcript) = 24 GB
- Read groups: ~1B × 2 ints × 4 bytes = 8 GB
- **References total: ~48 GB**

Other data: ~100 GB

**Total after: ~150 GB**
**Savings: ~550 GB (78% reduction)**

### Special Considerations

#### UMI Handling

UMIs might have millions of unique values. Options:

1. **Don't intern UMIs**: Keep as strings if cardinality too high
2. **Conditional interning**: Intern only if unique count < threshold (e.g., 1M)
3. **Compress instead**: Use shorter representation (e.g., base64 encoding)

**Recommendation**: Start without UMI interning, measure cardinality, decide later.

#### Backward Compatibility

- Keep property-based access to maintain API compatibility
- Support loading old serialization format (pure strings)
- Add version header to serialized files
- Gradual migration: workers can output interned format, old code can still read

#### Resume Functionality

- String pools must be saved/loaded during resume
- Pool IDs must be consistent across resume sessions
- Solution: Serialize pools at the start, load them back on resume

#### Read Group Discovery

Some read_group values are discovered during processing:
- BAM tag values
- Read ID prefixes

Strategy:
- Each worker maintains its own read_group pool
- Serialize pool per chromosome
- When combining results, pools remain separate (one per chromosome)
- Only merge if cross-chromosome operations needed (rare)

### Testing Strategy

1. **Unit tests**: StringPool add/get operations
2. **Serialization tests**: Round-trip with pools
3. **Memory benchmarks**: Measure before/after on test dataset
4. **Compatibility tests**: Old format → new format conversion
5. **Parallel processing tests**: Verify pool independence across workers
6. **Large-scale test**: Run on chromosome with 100M+ reads

### Rollout Plan

1. Implement Phase 1 (infrastructure) as optional feature (`--intern_strings`)
2. Test on small datasets
3. Implement Phase 2 (pre-known pools)
4. Benchmark memory usage on large dataset
5. Implement Phase 3 (dynamic pools) if Phase 2 successful
6. Implement Phase 4 (serialization) after runtime optimization proven
7. Make interning default after thorough testing
8. Eventually remove old string-based code

## Alternative Approaches Considered

### 1. Python's `sys.intern()`

**Pros**: Built-in, automatic deduplication in memory
**Cons**:
- No control over pool size
- No serialization benefit
- Doesn't work across processes
- Less memory savings (still stores full strings)

### 2. Flyweight Pattern with Weak References

**Pros**: Automatic garbage collection of unused strings
**Cons**:
- Complex to implement
- Overhead of weak reference management
- Doesn't help with serialization

### 3. String Compression (e.g., zlib)

**Pros**: Reduces disk usage
**Cons**:
- CPU overhead for compression/decompression
- No memory savings (must decompress to use)
- Doesn't address duplication problem

### 4. Custom String Encoding

For specific fields (e.g., barcodes with fixed alphabet):
- Encode ATCG as 2 bits per base
- 16-base barcode: 4 bytes instead of 16

**Pros**: Even more memory savings for specific fields
**Cons**: Complex, type-specific, harder to maintain

**Recommendation**: Consider as Phase 5 optimization for specific high-impact fields.

## Risks and Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| Integer overflow (> 2^31 values) | Critical | Use 64-bit ints for pools that might exceed 2B values |
| Pool memory overhead | Medium | Monitor pool sizes, add size limits, fallback to strings |
| Serialization compatibility | High | Version headers, support old format reading indefinitely |
| Performance regression | Medium | Benchmark thoroughly, optimize hot paths, use caching |
| Code complexity | Medium | Clear documentation, gradual rollout, property-based API |

## Success Metrics

- Memory usage reduction: Target 70-80%
- Serialization size reduction: Target 80-90%
- Performance impact: < 5% slowdown acceptable
- Code coverage: > 95% for new string pool code
- Production ready: No regressions in existing tests

## Conclusion

String interning via integer reference pools offers:
- **Massive memory savings** (550 GB → negligible for 1B reads)
- **Disk space savings** (proportional to memory savings)
- **Minimal code changes** (property-based backward compatibility)
- **Incremental rollout** (optional flag → default → required)

The parallel processing architecture fits well with per-worker pools, and the serialization format can efficiently store pools once per chromosome.

**Recommendation**: Proceed with implementation, starting with Phase 1 infrastructure and Phase 2 pre-known pools.
