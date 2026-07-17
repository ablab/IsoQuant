# Barcode Calling Architecture

This document describes the barcode calling subsystem for single-cell and spatial transcriptomics.

## Directory Structure

```
src/barcode_calling/
├── __init__.py                    # Re-exports all public classes
├── common.py                      # Utilities (str_to_2bit, bit_to_str, find_polyt_start, etc.)
│
├── indexers/                      # K-mer indexers
│   ├── __init__.py
│   ├── base.py                    # KmerIndexer, ArrayKmerIndexer
│   ├── two_bit.py                 # Dict2BitKmerIndexer, Array2BitKmerIndexer
│   └── shared_memory.py           # SharedMemoryArray2BitKmerIndexer, Numba JIT functions
│
├── callers/                       # Barcode detectors
│   ├── __init__.py
│   ├── base.py                    # Result classes (BarcodeDetectionResult, etc.)
│   ├── curio.py                   # CurioBarcodeDetector, CurioBruteForceDetector, CurioIlluminaDetector
│   ├── stereo.py                  # StereoBarcodeDetector, StereoSplittingBarcodeDetector, SharedMemory variants
│   └── tenx.py                    # TenXBarcodeDetector, VisiumHDBarcodeDetector
│
└── umi_filtering.py               # UMI-based deduplication
```

## K-mer Indexers

### Selection Guide

| Indexer | Use Case | Memory | Performance |
|---------|----------|--------|-------------|
| KmerIndexer | Small sets (<10K), variable length | Dict overhead | Good |
| ArrayKmerIndexer | Small sets, k≤8 | O(4^k) | Fast lookups |
| Dict2BitKmerIndexer | Medium-large sets (10K-1M) | ~40% less than KmerIndexer | Good |
| Array2BitKmerIndexer | Large fixed-length sets | Flat array | Fast |
| SharedMemoryArray2BitKmerIndexer | Very large sets (>1M) with multiprocessing | Shared memory | Numba JIT optimized |

### 2-Bit DNA Encoding

All 2-bit indexers use this encoding:
- A = 00
- C = 01
- G = 10
- T = 11

25bp sequence → 50 bits → stored in uint64

### SharedMemoryArray2BitKmerIndexer

Optimized for large barcode sets (>1M) with parallel processing:

**Features:**
- Two-pass construction: count k-mers first, then populate index
- Numba JIT compilation for 15-94x speedup (with fallback to pure Python)
- Shared memory allows worker processes to share index without copying
- Automatic cleanup via `__del__` method

**Key Methods:**
```python
# Create in main process
indexer = SharedMemoryArray2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=25)

# Get info for passing to workers
info = indexer.get_sharable_info()

# Reconstruct in worker process
worker_indexer = SharedMemoryArray2BitKmerIndexer.from_sharable_info(info)

# Query
matches = indexer.get_occurrences(sequence, max_hits=10, min_kmers=2)
```

## Barcode Detectors

### Platform-Specific Classes

**Curio Platform (formerly DoubleBarcodeDetector):**
- `CurioBarcodeDetector` - Main detector
- `CurioBruteForceDetector` - Exact linker matching (faster, less tolerant)
- `CurioIlluminaDetector` - Optimized for Illumina short reads

**Stereo-seq Platform:**
- `StereoBarcodeDetector` - Standard single-barcode mode
- `SharedMemoryStereoBarcodeDetector` - With shared memory for large whitelists
- `StereoSplittingBarcodeDetector` - Read splitting mode (concatenated reads)
- `SharedMemoryStereoSplittingBarcodeDetector` - Splitting with shared memory

**10x Genomics:**
- `TenXBarcodeDetector` - 10x v3 single-cell
- `TenXv2BarcodeDetector` - 10x v2 single-cell (UMI_LEN=10)
- `TenXSplittingBarcodeDetector` - 10x v3 split mode (concatenated reads → multiple FASTA); see algorithm docs below
- `TenXv2SplittingBarcodeDetector` - 10x v2 split mode (inherits from TenXSplittingBarcodeDetector, different TSO sequence)
- `VisiumHDBarcodeDetector` - Visium HD spatial (two-part barcode)

### Detection Flow (Standard Single-Barcode Modes)

1. Find polyT tail (indicates 3' end)
2. Find linker/primer sequences
3. Extract barcode region
4. Match against whitelist using k-mer index + Smith-Waterman alignment
5. Extract and validate UMI

### 10x Split Mode Algorithm (`TenXSplittingBarcodeDetector`)

Detects multiple barcodes in concatenated ONT reads where several cDNA molecules are ligated end-to-end during library preparation.

**Molecule template structure** (10x v3, one molecule on + strand):
```
R1(22bp) → Barcode(16bp) → UMI(12bp) → polyT → cDNA → TSO(32bp)
```
In concatenated reads, multiple molecules are joined via TSO boundaries. Molecules can alternate orientation (+ and - strand), so the read must be scanned in both directions.

**Entry point**: `find_barcode_umi(read_id, sequence) → SplittingBarcodeDetectionResult`

**Top-level flow**:
1. Scan forward strand for all molecule patterns via `_scan_strand(read_id, sequence, "+", result)`
2. Scan reverse complement for all molecule patterns via `_scan_strand(read_id, rev_seq, "-", result)`
3. If no patterns found on either strand, attempt single-molecule detection as fallback
4. Filter results: keep only valid barcoded detections (or best non-barcoded if none valid)

**`_scan_strand` loop** — iterates over one strand finding successive molecules:
```
current_start = 0
while remaining sequence >= MIN_REMAINING_SEQ (50bp):
    1. Find one barcode pattern in sequence[current_start:]
    2. If no polyT found → break (no more molecules on this strand)
    3. Check consistency (R1 and polyT within MAX_R1_POLYT_DISTANCE = 50bp)
       - Consistent: append result, advance to TSO end (or polyT + 100 if no TSO)
       - Inconsistent: skip to R1 position (polyT was spurious, R1 belongs to next molecule)
    4. Enforce MIN_SPLIT_STEP (150bp) to prevent overlapping detections
```

**Per-molecule detection** (`_find_barcode_umi_split_fwd`):

Step 1 — `_find_barcode_umi_fwd_local` (local-only barcode search):
1. Find polyT via sliding window (`find_polyt_start`)
2. Search for R1 adapter in `sequence[0:polyT+1]` only (NOT the full read)
   - Uses relaxed k-mer matching (min_score=11, R1 k-mer indexer)
   - No fallback to full-read search — eliminates O(n) search on long reads
3. Extract 16bp barcode region immediately after R1
4. Match barcode against whitelist (k-mer index + Smith-Waterman, min_score=14 for 9K whitelist)
5. Extract UMI between barcode end and polyT start

Step 2 — TSO search (only if valid barcode found):
1. Search from polyT to end of read (no distance cap — cDNA can be 5-6kb+)
2. Uses `detect_first_exact_positions` — takes the first acceptable TSO match (nearest to polyT), not best-scoring
3. TSO k-mer indexer: k=9, min_score=18, start_delta=3, end_delta=3 (relaxed vs standard detection)

**Consistency check** (`_is_consistent_detection`):
- Validates that R1 end is before polyT and within 50bp (expected ~28bp: BC=16 + UMI=12)
- Rejects spurious polyT detections in cDNA that would match R1 from a different molecule
- When rejected, the scan skips to `R1 - len(R1) - 10` so the next iteration starts before the real R1

**Output pipeline** (`_process_read_split` in `detect_barcodes.py`):
- Each detected molecule becomes a separate FASTA record
- Segment coordinates: `R1 - 10` to `TSO + 35` (or read end if no TSO)
- Read ID suffix: `_{start}_{end}_{strand}` for provenance tracking
- If multiple patterns detected, only molecules with TSO are emitted (`require_tso` flag)

**Performance optimizations**:
- `_find_barcode_umi_fwd_local` eliminates the expensive full-read R1 fallback search from the base `_find_barcode_umi_fwd` (7.4x speedup on 1K reads)
- TSO search gated on `is_valid()` — only searches when a barcode was actually found, avoiding unnecessary k-mer extraction on spurious polyT hits
- `detect_first_exact_positions` for TSO — early termination on first acceptable match instead of scanning all k-mer hits for best score

**Key parameters** (`TenXSplittingBarcodeDetector`):
| Parameter | Value | Purpose |
|-----------|-------|---------|
| `MAX_R1_POLYT_DISTANCE` | 50 | Max bp between R1 end and polyT for consistency |
| `MIN_REMAINING_SEQ` | 50 | Stop scanning when less than this remains |
| `DEFAULT_POLYT_STEP` | 100 | Step past polyT when no TSO found |
| `MIN_SPLIT_STEP` | 150 | Minimum advance between scan iterations |
| TSO k-mer size | 9 | K-mer size for TSO indexer |
| TSO min_score | 18 | Minimum SSW score for TSO detection |
| TSO start/end_delta | 3 | Allowed terminal mismatch in TSO alignment |

**CI tests**: `Barcode.Mouse.10x_concat.yml` — two datasets:
- `Mouse.10x_concat.3M.800K` — 800K simulated reads, 3M full whitelist (stress test)
- `Mouse.10x_concat.realistic.800K` — 800K simulated reads, 9K realistic whitelist

### Result Classes

```
BarcodeDetectionResult (base)
├── LinkerBarcodeDetectionResult (adds polyT, primer, linker positions)
│   └── TSOBarcodeDetectionResult (adds TSO5 position — Stereo-seq)
└── TenXBarcodeDetectionResult (adds R1, polyT positions)
    └── TenXSplitBarcodeDetectionResult (adds tso position, get_fasta_segment_start/end)

SplittingBarcodeDetectionResult (list of TSOBarcodeDetectionResult or TenXSplitBarcodeDetectionResult)
```

**Generic segment interface** (for `_process_read_split` in `detect_barcodes.py`):
- `get_fasta_segment_start() -> int` — start of FASTA segment to extract
- `get_fasta_segment_end(seq_len: int) -> int` — end of FASTA segment
- `get_tso_position() -> int` — TSO position (-1 if not found); used for `require_tso` check

Both `TSOBarcodeDetectionResult` and `TenXSplitBarcodeDetectionResult` implement this interface, making `_process_read_split` platform-agnostic.

## detect_barcodes.py Integration

The `BARCODE_CALLING_MODES` dictionary maps IsoQuantMode to detector class:

```python
BARCODE_CALLING_MODES = {
    IsoQuantMode.tenX_v3: TenXBarcodeDetector,
    IsoQuantMode.tenX_v2: TenXv2BarcodeDetector,
    IsoQuantMode.tenX_v3_split: TenXSplittingBarcodeDetector,
    IsoQuantMode.tenX_v2_split: TenXv2SplittingBarcodeDetector,
    IsoQuantMode.curio: CurioBarcodeDetector,
    IsoQuantMode.stereoseq_nosplit: SharedMemoryStereoBarcodeDetector,
    IsoQuantMode.stereoseq: SharedMemoryStereoSplittingBarcodeDetector,
    IsoQuantMode.visium_5prime: TenXBarcodeDetector,
    IsoQuantMode.visium_hd: VisiumHDBarcodeDetector,
}
```

## Common Utilities (common.py)

| Function | Purpose |
|----------|---------|
| `str_to_2bit(seq)` | Convert DNA string to 2-bit integer |
| `bit_to_str(val, length)` | Convert 2-bit integer back to string |
| `reverese_complement(seq)` | Reverse complement (note: preserved typo) |
| `find_polyt_start(seq)` | Detect polyT tail using sliding window |
| `align_pattern_ssw(pattern, text, ...)` | Smith-Waterman alignment |
| `find_candidate_with_max_score_ssw(...)` | Select best barcode from candidates |
| `detect_exact_positions(...)` | Locate pattern using k-mer index hits |

## Multiprocessing Notes

- Uses spawn context (not fork) for compatibility
- SharedMemory indexers serialize/deserialize via `__getstate__`/`__setstate__`
- Worker processes attach to existing shared memory blocks
- Main process is responsible for cleanup (`unlink()`)

## Adding New Platforms

1. Create result class in `callers/base.py` if needed
2. Create detector class in appropriate file (or new file)
3. Add to `callers/__init__.py` exports
4. Add to `src/barcode_calling/__init__.py` exports
5. Add to `detect_barcodes.py` BARCODE_CALLING_MODES
6. Update BARCODE_FILES_REQUIRED if whitelist format differs

## Barcode/UMI Pipeline Integration

Barcode/UMI data is integrated into `ReadAssignment` objects (loaded once during read collection, persisted through pipeline):
- `ReadAssignment.barcode` and `.umi` fields in `src/isoform_assignment.py`
- Loaded from split barcode files during `collect_reads_in_parallel()` in `src/dataset_processor.py`
- `AlignmentCollector` (`src/alignment_processor.py`) receives `barcode_dict`, populates fields
- `UMIFilter` (`src/barcode_calling/umi_filtering.py`) reads barcode/umi directly from ReadAssignment
- Barcode table splitting at pipeline start (`split_read_barcode_table()`), each file loaded once

## Barcode-Spot Grouping

Groups read counts by cell type/spatial region instead of individual barcodes.
- `BarcodeSpotGrouper` in `src/read_groups.py` — maps read_id -> barcode -> spot/cell type
- CLI: `--read_group barcode_spot` or `--read_group barcode_spot:FILE`
- Works with multi-grouping: `--read_group tag:CB file_name barcode_spot`

## Barcode-to-Barcode Mapping (--barcode2barcode)

Enables spot-level UMI deduplication: reads with different barcodes but the same spatial spot ID are grouped together for dedup.

**CLI**: `--barcode2barcode file.tsv` or `--barcode2barcode file.tsv:barcode_col:spot_cols`

**File format**: Same TSV spec as `--barcode2spot`:
```
ATGCATGC    spot_A    region_1
GCTAGCTA    spot_A    region_2
TTTTAAAA    spot_B    region_1
```

**UMI dedup behavior** (in `DatasetProcessor.filter_umis()`):
- Round 0 (standard): Original barcodes, produces `.allinfo` AND `.filtered_*` read IDs
- Rounds 1-N (per spot column): Grouped by spot_id via `barcode_remap`, produces `.allinfo` only (no `.filtered_*`)
- Output naming: `{sample}.UMI_filtered.barcode_barcode_col{K}.ED{N}.allinfo`

**Implementation across files**:
- `isoquant.py`: CLI argument + file existence check
- `src/read_groups.py`: `barcode_barcode` handling in `parse_grouping_spec()`, `get_grouping_strategy_names()`, `get_grouping_pool_types()`. Reuses `BarcodeSpotGrouper` and `SharedTableData`. Also `load_barcode2barcode_mapping()` helper
- `src/barcode_calling/umi_filtering.py`: `UMIFilter` accepts optional `barcode_remap: Dict[str, str]`. When set, `process_single_chr()` uses remapped key instead of barcode for grouping in `gene_barcode_dict`
- `src/parallel_workers.py`: `filter_umis_in_parallel()` accepts `barcode_remap` and `output_prefix_override` params. Forces `output_filtered_reads=False` when remap is active
- `src/dataset_processor.py`: After standard dedup loop, iterates spot columns calling `filter_umis_in_parallel()` with column-specific `barcode_remap` and prefix
- `src/file_naming.py`: `umi_barcode2barcode_prefix()`, `umi_barcode2barcode_global_lock()` helpers

**Read grouping**: `--read_group barcode_barcode` works identically to `barcode_spot` but reads from `--barcode2barcode`. Auto-added when `--barcode2barcode` is set (same pattern as `--barcode2spot` auto-adds `barcode_spot`).

**Utility script**: `misc/prepare_visium_spot_ids.py` generates the barcode-to-spot TSV from two per-part coordinate files (part1→Y, part2→X) via cartesian product. Produces composite barcodes (`BC1|BC2`) with spot IDs formatted as `prefix + Y_coord + delimiter + X_coord + suffix` for each resolution column. Default prefixes match Visium HD format (`s_002um_`, `s_008um_`, `s_016um_`). Supports custom prefix/suffix/delimiter. Output is directly usable as `--barcode2barcode` input.

**Tests**: `tests/test_barcode2barcode.py` — 26 tests covering file naming, mapping loading, grouper behavior, parse_grouping_spec, strategy names, pool types, UMI remap grouping, and dedup logic.

## Barcoded BAM (--barcoded_bam)

Reads barcodes and UMIs directly from BAM tags (e.g. cellranger output) instead of running barcode calling.

**CLI**: `--barcoded_bam` (flag), optionally with `--barcode_tag CB --umi_tag UB` (defaults) and `--strip_barcode_suffix` (removes trailing `-1` etc.)

**Behavior**:
- Bypasses barcode calling entirely — no FASTQ, no whitelist needed
- `AlignmentCollector` (`src/alignment_processor.py`) reads `CB` and `UB` tags from each BAM record
- If `--strip_barcode_suffix` is set, removes the portion after the last `-` in the barcode (e.g. `ATCG-1` → `ATCG`)
- Barcode/UMI populated into `ReadAssignment` objects, same as normal barcode calling path
- Compatible with all downstream features: UMI dedup, read grouping, barcode2barcode
- Skips `detect_barcodes` and `split_read_barcode_table` steps (`DatasetProcessor` checks `barcoded_bam` flag)

**Implementation across files**:
- `isoquant.py`: CLI arguments (`--barcoded_bam`, `--barcode_tag`, `--umi_tag`, `--strip_barcode_suffix`), validation
- `src/alignment_processor.py`: `AlignmentCollector.__init__` stores tag names and suffix stripping flag; `_collect_alignments_with_barcodes()` and `_process_alignment_group()` read tags from `pysam.AlignedSegment`
- `src/dataset_processor.py`: Skips barcode table loading when `barcoded_bam` is set
- `src/read_groups.py`: `barcode_spot` and `barcode_barcode` grouping accept `barcoded_bam` as valid barcode source

**Data prep utility**: `misc/inject_bam_tags.py` — injects CB/UB tags from a `barcoded_reads.tsv` file into an existing BAM. Used to create test data from IsoQuant's own barcode calling output. Supports `--suffix` to add cellranger-style `-1` suffixes.

## Universal Barcode Calling

Custom molecule structures via Molecule Definition Files (MDF), any sequencing platform without code changes.

**CLI**: `isoquant.py --molecule molecule_definition.mdf ...`

**Location**: `src/barcode_calling/callers/`

**Key Files**:
- `molecule_structure.py` - MDF parsing, ElementType enum, MoleculeElement/MoleculeStructure
- `universal_extraction.py` - UniversalSingleMoleculeExtractor for pattern detection
- `extraction_result.py` - ExtractionResult dict-based result, DetectedElement container
- `protocol.py` - BarcodeResult protocol interface

### MDF Format

```
# Line 1: Element order (colon-separated)
R1:Barcode:UMI:PolyT:cDNA:TSO

# Subsequent lines: element_name  type  value [optional_length]
R1       CONST      CTACACGACGCTCTTCCGATCT
Barcode  VAR_FILE   barcodes.tsv
UMI      VAR_ANY    12
TSO      CONST      CCCATGTACTCTGCGTTGATACCACTGCTT
```

### ElementType Enum

| Type | Description | Value field |
|------|-------------|-------------|
| `PolyT` | PolyT/PolyA tract | (none) |
| `cDNA` | cDNA region | (none) |
| `CONST` | Constant sequence | Sequence string |
| `VAR_ANY` | Variable (extract as-is) | Length (int) |
| `VAR_LIST` | Variable from inline list | Comma-separated |
| `VAR_FILE` | Variable from file | TSV file path |
| `VAR_ANY_SEPARATOR` | Separator | Length (int) |
| `VAR_ANY_NON_T_SEPARATOR` | Non-T separator | Length (int) |

**Key methods**: `needs_correction()` (VAR_LIST/FILE), `needs_sequence_extraction()` (VAR_ANY), `is_constant()` (CONST), `is_variable()` (VAR_ANY/LIST/FILE), `is_base_separator()`

### Barcode/UMI Identification

By naming convention (case-insensitive prefix):
- "barcode" prefix -> barcode elements (e.g., `Barcode`, `barcode1`)
- "umi" prefix -> UMI elements (e.g., `UMI`, `umi1`)

### Linked Elements

Two types of linked elements for multi-part barcodes/UMIs. Both require:
- All parts variable type (VAR_ANY, VAR_LIST, VAR_FILE)
- Same type and values across parts (same whitelist)
- 1-indexed, consecutive (1, 2, 3, ...)

**Concatenated** (`|`): Parts of a single barcode separated on the molecule.
- MDF: `Barcode|1:Barcode|2:PolyT:cDNA` with `Barcode|1 VAR_LIST whitelist 8`
- `elements_to_concatenate`: `"Barcode|1" -> ("Barcode", 1)`
- `concatenated_elements_counts`: `"Barcode" -> 2`
- Parts extracted independently, concatenated, corrected as one sequence
- `__init__`: ONE index per base name using total_length for k-mer size
- `process_concatenated_elements(storage, results, sequence)`
- Without correction (`correct_sequences=False`): raw concatenation, no index

**Duplicated** (`/`): Multiple copies of the same barcode (redundancy).
- MDF: `Barcode/1:Barcode/2:PolyT:cDNA`
- `duplicated_elements`: `"Barcode/1" -> ("Barcode", 1)`
- `duplicated_elements_counts`: `"Barcode" -> 2`
- Each copy corrected independently, then majority vote
- `__init__`: ONE shared index per base name using individual part length
- `process_duplicated_elements(storage, results, sequence)`:
  - With correction: both agree -> report; one succeeds -> report; disagree -> NOSEQ; both fail -> NOSEQ
  - Without correction: all match -> report; differ -> comma-separated; one detected -> report

**Shared infrastructure**:
- `_identify_barcode_umi_elements()`: uses base name for both types (deduplicates)
- `_backup_concatenated_element()` / `_backup_duplicated_element()`: save coordinates during extraction
- `_detect_variable_element()`: early return for linked elements -> backup only
- `ExtractionResult.__str__()` / `.header()`: skip individual parts (data under base name)

### Extraction Flow

`UniversalSingleMoleculeExtractor.find_barcode_umi`:
1. Try forward and reverse complement
2. `_find_patterns_fwd()`:
   a. Detect polyT (`find_polyt`)
   b. Detect constant elements (`detect_const_elements`) — k-mer index + SSW alignment
   c. Extract variable elements from 5' side (`extract_variable_elements5`)
   d. Extract variable elements from 3' side (`extract_variable_elements3`)
   e. `process_concatenated_elements` — combine parts, correct against full whitelist
   f. `process_duplicated_elements` — majority vote across copies
3. Choose strand by polyT presence, or by more detected elements

### Test Data

**Unit tests**: `tests/universal_data/`: `10x.mdf`, `simple.mdf`, `barcodes_small.tsv`, `test_reads.fq`

Read name convention: `@READ_1_ATCCTTAGTGTTTGTC_AACCTTAGACCG_forward`

**Unit tests**: `tests/test_molecule_structure.py`, `tests/test_universal_extraction.py`

### CI Integration Tests (Universal Barcode Calling)

CI tests validate that universal barcode calling (`custom_sc` mode with MDF files) produces equivalent results to protocol-specific callers. Tests are in `/abga/work/andreyp/ci_isoquant/data/barcodes/` and run via `tests/github/run_barcode_test.py` (which supports a `molecule` config key for MDF files).

**Test runner**: `tests/github/run_barcode_test.py` — `molecule` config key adds `--molecule <mdf_file>` to detect_barcodes command.

**Protocols tested** (8 tests total, 2 per protocol):

| Protocol | Small variant | Large variant | Small indexer | Large indexer |
|----------|--------------|---------------|---------------|---------------|
| 10x | `Mouse.10x.custom_sc.small` | `Mouse.10x.custom_sc.large` | ArrayKmerIndexer | Array2BitKmerIndexer |
| Curio | `Human.Curio.custom_sc.small` | `Human.Curio.custom_sc.large` | ArrayKmerIndexer | Array2BitKmerIndexer |
| Visium HD | `Mouse.VisiumHD.custom_sc.small` | `Mouse.VisiumHD.custom_sc.large` | ArrayKmerIndexer | Array2BitKmerIndexer |
| Stereo-seq | `Mouse.StereoSeq.custom_sc.small` | `Mouse.StereoSeq.custom_sc.large` | KmerIndexer | Dict2BitKmerIndexer |

**Indexer path coverage**: The `prepare_barcode_index()` function in `universal_extraction.py` selects indexer type based on whitelist size (threshold: 100K) and kmer_size (threshold: 8). Each protocol has a "small" variant (<100K barcodes) and a "large" variant (>=100K barcodes) to exercise both code paths.

**Large barcode sets**: Created by inflating original whitelists with random dummy DNA sequences (script: `misc/create_test_barcodes.py`). Dummy barcodes don't match any reads, so detection results stay identical — only the indexer type changes.

**Stereo test data**: Barcode subsets (80K and 1M) sampled from the full 10M whitelist, with reads filtered to matching barcodes (script: `misc/create_stereo_test_subsets.py`).

**MDF files**: 8 files (`10x.small.mdf`, `10x.large.mdf`, `curio.small.mdf`, `curio.large.mdf`, `visium.small.mdf`, `visium.large.mdf`, `stereo.small.mdf`, `stereo.large.mdf`) — all use absolute paths for VAR_FILE references since MoleculeStructure resolves paths relative to CWD.

**GitHub Actions workflows**:
- `Barcode.Custom.10x.yml` — Monday 4:00 AM
- `Barcode.Custom.Curio.yml` — Wednesday 3:00 AM
- `Barcode.Custom.VisiumHD.yml` — Thursday 2:00 AM
- `Barcode.Custom.StereoSeq.yml` — Friday 1:00 AM

**Note**: Visium HD results are poor (low recall for large variant) because universal barcode calling is not designed for variable-length similar barcodes. Tests are kept for regression detection.

## BarcodeResult Protocol and Class Hierarchy

**Protocol** (`protocol.py`): `BarcodeResult` — PEP 544 structural typing, `runtime_checkable`

```
BarcodeResult (Protocol)
+-- BarcodeDetectionResult (base class, attribute-based)
|   +-- LinkerBarcodeDetectionResult (double barcodes: Curio, Stereo-seq)
|   |   +-- TSOBarcodeDetectionResult (Stereo-seq with TSO)
|   +-- TenXBarcodeDetectionResult (10x Genomics)
+-- SplittingBarcodeDetectionResult (container for multiple patterns)
+-- ExtractionResult (dict-based, universal extraction)
```

**Common interface**: `get_barcode()`, `get_umi()`, `is_valid()`, `has_barcode()`, `has_umi()`, `set_strand()`, `update_coordinates()`, `more_informative_than()`, `get_additional_attributes()`, `header()`, `__str__()`

Prefer `get_barcode()`/`get_umi()` over `.barcode`/`.UMI` attributes in new code.

## PolyT-Anchored Barcode Fallback

### Problem
When R1/linker is partially present or error-corrupted, the k-mer seeding step produces zero seeds, so SSW alignment is never attempted and the barcode is not detected. This affects ~450K reads in 10x datasets where scisorseqR (exact matching) finds barcodes but IsoQuant does not.

### Root Cause
R1 k-mer indexer uses k=7 for seeding. Partial R1 (<12bp) or error-corrupted R1 can have zero matching 7-mers, so the alignment is never attempted. Without R1 detection, no barcode search is triggered.

### Solution: PolyT-Anchored Search
When R1/linker is not found but polyT is detected, search for the barcode directly in the region before polyT using the barcode whitelist. The polyT anchor constrains the search region since the molecule structure is known: `...BC → UMI → polyT → cDNA...`

### Implementation (10x — committed)

**`TenXBarcodeDetector._find_barcode_near_polyt()`** (`callers/tenx.py`):
- Searches region `[polyT - BC_LEN - UMI_LEN - MARGIN, polyT)` for barcode
- Uses existing `barcode_indexer.get_occurrences()` + `find_candidate_with_max_score_ssw()`
- Extracts UMI from between barcode end and polyT
- Returns `TenXBarcodeDetectionResult` with barcode, UMI, score

**Non-split mode** (`TenXBarcodeDetector._find_barcode_umi_fwd`):
- Called inline when R1 not found but polyT exists — safe because called once per strand

**Split mode** (`TenXSplittingBarcodeDetector.find_barcode_umi`):
- Called as **per-read fallback** after scanning both strands, only when no valid barcodes found
- MUST NOT be called inline in `_find_barcode_umi_fwd_local` — doing so causes massive false positives in concatenated reads (spurious polyT in cDNA triggers random barcode matches at every scan position)

### Test Results (10x concat, realistic 9K whitelist)
- **Precision**: 99.80% (baseline 99.89%)
- **Recall**: 86.80% (baseline 85.81%, +1%)
- **excessive_assignments**: 190 (baseline 188, no regression)
- **incorrect_assignments**: 2,286 (baseline 1,197, +1,089 tradeoff for +13,417 correct)
- Etalon needs update to reflect improved recall

### Stereo-seq: Not Implemented (deferred)

PolyT-anchored fallback was prototyped for Stereo-seq but **reverted** due to high false positive rate:
- With 452M barcode whitelist (25bp barcodes), random matches dominate the search region
- Stereo concat test: 4,178 incorrect (baseline 651) with ~0 new correct assignments
- Stereo non-split test: marginal improvement (+1,146 reads, precision 99.81%)

**Potential approaches for future work**:
1. Require higher `min_score` (e.g., 24 out of 25) for the fallback with large whitelists
2. Require both barcode AND valid UMI length for fallback results
3. Use linker fragment detection as additional confirmation
4. Only enable fallback when whitelist is small (<1M barcodes)
