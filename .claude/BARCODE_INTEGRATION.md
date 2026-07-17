# Barcode/UMI Integration Optimization

## Overview

This document describes the optimization implemented in December 2024 to integrate barcode and UMI data into ReadAssignment objects, eliminating redundant I/O operations during UMI filtering.

## Problem Statement

Prior to this change, barcode tables were being loaded twice in the pipeline:
1. **During read collection**: Loaded to potentially use during alignment processing
2. **During UMI filtering**: Loaded again to associate barcodes/UMIs with reads for PCR deduplication

This caused:
- Redundant I/O operations (reading the same files twice)
- Increased memory usage (barcode dicts held in memory during both stages)
- Slower pipeline execution

## Solution Design

### Core Concept

Store barcode and UMI information as properties of ReadAssignment objects so they:
1. Are loaded once during read collection
2. Persist through serialization/deserialization
3. Are available during UMI filtering without reloading files

### Implementation Strategy

1. Add `barcode` and `umi` fields to ReadAssignment
2. Load barcodes during read collection and populate these fields
3. Serialize/deserialize these fields with ReadAssignment
4. Modify UMI filtering to read from ReadAssignment instead of loading files

## Pipeline Architecture

### Barcode Table Splitting (Early in Pipeline)

Barcode tables are split by chromosome **before** read collection begins, in `process_assigned_reads()`:

**File**: `src/dataset_processor.py` (lines 484-498)

```python
if self.args.mode.needs_pcr_deduplication():
    split_barcodes_dict = {}
    for chr_id in self.get_chr_list():
        split_barcodes_dict[chr_id] = sample.barcodes_split_reads + "_" + chr_id

    barcode_split_done = split_barcodes_lock_filename(sample)
    if self.args.resume and os.path.exists(barcode_split_done):
        logger.info("Barcode table was split during the previous run")
    else:
        self.split_read_barcode_table(sample, split_barcodes_dict)
        open(barcode_split_done, "w").close()
```

This ensures split barcode files are available when `collect_reads_in_parallel()` runs.

**Cleanup** happens after `load_read_info()` (lines 516-520) to ensure files are available throughout the pipeline.

## Implementation Details

### Phase 1: ReadAssignment Enhancement

**File**: `src/isoform_assignment.py`

Added fields to `ReadAssignment` class:
```python
# Line 638-639
self.barcode = None  # Cell/spatial barcode
self.umi = None      # Unique molecular identifier
```

Added fields to `BasicReadAssignment` for serialization:
```python
# Line 489-490
self.barcode = read_assignment.barcode if hasattr(read_assignment, 'barcode') else None
self.umi = read_assignment.umi if hasattr(read_assignment, 'umi') else None
```

**Serialization/Deserialization**:

Added to `ReadAssignment.serialize()` (lines 715-716):
```python
write_string_or_none(self.barcode, outfile)
write_string_or_none(self.umi, outfile)
```

Added to `ReadAssignment.deserialize()` (lines 686-687):
```python
read_assignment.barcode = read_string_or_none(infile)
read_assignment.umi = read_string_or_none(infile)
```

**Backward Compatibility**:
- Old serialization format (without barcode/umi) can still be read
- `BasicReadAssignment` old-format deserialization includes dummy reads for barcode/umi fields (lines 580-581)
- This ensures existing saved data can be loaded without errors

### Phase 2: Barcode Loading During Read Collection

**File**: `src/dataset_processor.py`

Modified `collect_reads_in_parallel()` function (lines 121-137):

```python
# Load barcode dict for this chromosome if available
barcode_dict = {}
if sample.barcodes_split_reads:
    barcode_file = sample.barcodes_split_reads + "_" + chr_id
    if os.path.exists(barcode_file):
        logger.debug(f"Loading barcodes from {barcode_file}")
        for line in open(barcode_file):
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 3:
                barcode_dict[parts[0]] = (parts[1], parts[2])

# Pass to AlignmentCollector
alignment_collector = AlignmentCollector(
    chr_id, bam_file_pairs, args, illumina_bam, gffutils_db,
    current_chr_record, read_grouper, barcode_dict,
    args.max_coverage_small_chr, args.max_coverage_normal_chr)
```

**Barcode File Format**:
- TSV with 3 columns: `read_id`, `barcode`, `umi`
- Created by table_splitter.py during barcode table splitting
- One file per chromosome: `{sample.barcodes_split_reads}_{chr_id}`

### Phase 3: AlignmentCollector Update

**File**: `src/alignment_processor.py`

Added barcode_dict parameter and population logic:

```python
# Line 239 - __init__ parameter
def __init__(self, chr_id, bam_pairs, params, illumina_bam,
             genedb=None, chr_record=None, read_groupper=DefaultReadGrouper(),
             barcode_dict=None,
             small_chr_max_coverage=1000000,
             usual_gene_max_coverage=-1):
    ...
    self.barcode_dict = barcode_dict if barcode_dict else {}

# After line 390 and 453 - populate when creating ReadAssignment
if read_id in self.barcode_dict:
    read_assignment.barcode, read_assignment.umi = self.barcode_dict[read_id]
```

### Phase 4: UMIFilter Simplification

**File**: `src/barcode_calling/umi_filtering.py`

Removed split_barcodes_dict dependency:

```python
# Line 159 - Removed split_barcodes_dict parameter
def __init__(self, umi_length=0, edit_distance=3, disregard_length_diff=True,
             only_unique_assignments=False, only_spliced_reads=False):
    # No longer stores self.split_barcodes_dict

# Lines 521-522 - Read from ReadAssignment in process_single_chr()
# No hasattr checks needed - barcode/umi always exist after serialization
barcode = read_assignment.barcode
umi = read_assignment.umi

# Lines 589-590 - Read from ReadAssignment in process_from_raw_assignments()
# No hasattr checks needed - barcode/umi always exist after serialization
barcode = read_assignment.barcode
umi = read_assignment.umi
```

Removed barcode file loading:
- Removed `barcode_dict = self.load_barcodes_simple(self.split_barcodes_dict[chr_id])` from both processing methods

### Phase 5: Function Signature Updates

**File**: `src/dataset_processor.py`

Updated `filter_umis_in_parallel()`:

```python
# Line 257 - Removed split_barcodes_dict parameter
def filter_umis_in_parallel(sample, chr_id, args, edit_distance, output_filtered_reads=False):
    ...
    # Line 274 - Updated UMIFilter instantiation
    umi_filter = UMIFilter(args.umi_length, edit_distance)
    ...

# Lines 798-805 - Updated parallel call generator
umi_gen = (
    filter_umis_in_parallel,
    itertools.repeat(sample),
    self.get_chr_list(),
    itertools.repeat(self.args),           # No split_barcodes_dict
    itertools.repeat(edit_distance),
    itertools.repeat(output_filtered_reads),
)
```

## Benefits

1. **Reduced I/O**: Each barcode file is read once per chromosome instead of twice
2. **Lower Memory Usage**: Barcode dict not held in memory during UMI filtering stage
3. **Better Resume Functionality**: Barcode/UMI data persists through serialization
4. **Cleaner Architecture**: Barcode data flows naturally through the pipeline
5. **Faster Execution**: Eliminated redundant file parsing

## Backward Compatibility

The implementation uses `hasattr()` checks when reading barcode/umi from ReadAssignment objects, ensuring compatibility with:
- Old serialized data that doesn't have these fields
- Non-barcoded runs where these fields are None

## Testing Considerations

When testing this implementation:
1. Verify barcode files are loaded during read collection
2. Verify barcode/umi populate ReadAssignment objects correctly
3. Verify UMI filtering works without loading barcode files
4. Test serialization/deserialization of barcode/umi fields
5. Test resume functionality with barcoded data
6. Test non-barcoded runs still work correctly

## Related Files

- `src/isoform_assignment.py` - ReadAssignment data structures
- `src/alignment_processor.py` - Alignment collection and ReadAssignment creation
- `src/dataset_processor.py` - Pipeline orchestration
- `src/barcode_calling/umi_filtering.py` - UMI deduplication
- `src/table_splitter.py` - Barcode table splitting (unchanged)

## Future Enhancements

Potential improvements:
1. Add barcode validation during read collection
2. Add metrics for barcode coverage per chromosome
3. Consider streaming barcode loading for very large files
4. Optimize barcode dict memory usage with more compact representations
