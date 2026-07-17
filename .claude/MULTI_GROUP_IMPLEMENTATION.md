# Multi-Group Read Grouping Implementation

This document describes the multi-group read grouping functionality implemented in IsoQuant.

## Overview

IsoQuant now supports **multiple simultaneous grouping strategies** for reads. Users can group reads by multiple criteria at once (e.g., by cell barcode AND by UMI AND by file name).

## Implementation Timeline

Based on git history (commits from Dec 17-18, 2025):

1. **4d82cad** - "Support for multiple read groups" - Initial implementation
2. **194c4ce** - Fix deserialization for list read_group
3. **df2f38b** - Temporary fix using read_group[0] in counters
4. **0c1b37d** - Join read_group list in TSV output
5. **ba7e0b3** - Use first group in model construction temporarily
6. **385e174** - **Implement proper multi-group counters with strategy-based naming**
7. **ebbab3e** - Fix references to grouped counters (now lists)
8. **de3ed9d** - Update test to expect new grouped file naming format
9. **df3a542** - Fix multi-group support in GraphBasedModelConstructor
10. **d51f06b** - **Fix technical replicas check to only use file_name group**

## Key Concepts

### 1. Read Group Specification

Users specify grouping strategies via `--read_group` (nargs='+'):

```bash
# Single grouping strategy
--read_group tag:CB

# Multiple grouping strategies
--read_group tag:CB tag:UB file_name

# File-based grouping with multiple columns
--read_group file:barcodes.tsv:0:1,2,3
```

### 2. Data Structure Changes

**Before (single group):**
- `read_assignment.read_group` was a single string
- Counters were simple dictionaries

**After (multi-group):**
- `read_assignment.read_group` is a **list** of group IDs (one per grouping strategy)
- Counters are **lists of dictionaries** (one dict per grouping strategy)

### 3. Key Components

#### A. Read Groupers (`src/read_groups.py`)

**Single Groupers:**
- `DefaultReadGrouper` - No grouping (returns "NA")
- `AlignmentTagReadGrouper` - Groups by BAM tag (e.g., CB for cell barcode)
- `ReadIdSplitReadGrouper` - Groups by read ID suffix
- `FileNameGrouper` - Groups by input file name
- `ReadTableGrouper` - Groups by TSV file lookup
- `MultiColumnReadTableGrouper` - Multiple group columns from one TSV file

**Multi-Grouper:**
- `MultiReadGrouper` - Manages multiple groupers, returns list of group IDs

**Key Functions:**
- `create_read_grouper(args, sample, chr_id)` - Creates grouper(s) from command-line args
- `get_grouping_strategy_names(args)` - Returns list of strategy names (e.g., ["tag_CB", "tag_UB", "file_name"])
- `prepare_read_groups(args, sample)` - Prepares file-based grouping tables

#### B. Counter Architecture (`src/long_read_counter.py`)

**AbstractCounter:**
- `group_index` - Which group in the read_group list to use
- `read_groups` - Set of group IDs for this specific grouping strategy

**GroupedFeatureCounter:**
- When `ignore_read_groups=False` and multiple groups exist:
  - Creates one counter per grouping strategy
  - Each counter tracks its own `group_index`
  - Counters are stored as lists when multiple groups are active

**Key Changes:**
- Line 262: `group_id = read_assignment.read_group[self.group_index]` - Extract specific group from list
- Counters now iterate over `group_index` positions

#### C. File Naming (`src/file_naming.py`)

**Strategy-Based Naming:**
- When multiple grouping strategies exist, files are named using strategy names
- Example: `gene_grouped.tag_CB_counts.tsv`, `gene_grouped.file_name_counts.tsv`
- Uses `get_grouping_strategy_names()` to generate consistent names

#### D. Dataset Processor (`src/dataset_processor.py`)

**Multi-Group Support:**
- Line 454: `self.args.use_technical_replicas = self.args.read_group == "file_name" and len(sample.file_list) > 1`
- Creates multiple grouped counters (one per strategy)
- Passes `grouping_strategy_names` to components

#### E. Graph-Based Model Construction (`src/graph_based_model_construction.py`)

**Technical Replicas Check:**
- Line 63: `self.file_name_group_idx = self.grouping_strategy_names.index("file_name") if "file_name" in self.grouping_strategy_names else -1`
- Line 510-519: Uses `file_name_group_idx` to check if reads come from same file
- Line 514: `file_groups = set([a.read_group[self.file_name_group_idx] for a in read_assignments])`
- **Critical Fix (d51f06b)**: Only uses file_name group index, not just read_group[0]

## Usage Examples

### Example 1: Single-cell data with cell barcodes and file names

```bash
isoquant.py \
  --reference genome.fa \
  --genedb annotation.gtf \
  --bam sample1.bam sample2.bam \
  --data_type nanopore \
  --read_group tag:CB file_name \
  -o output/
```

**Result:**
- Reads grouped by both cell barcode (CB tag) and file name
- Output files:
  - `gene_grouped.tag_CB_counts.tsv` (cells × genes)
  - `gene_grouped.file_name_counts.tsv` (files × genes)

### Example 2: Multiple tags

```bash
isoquant.py \
  --reference genome.fa \
  --genedb annotation.gtf \
  --bam sample.bam \
  --data_type nanopore \
  --read_group tag:CB tag:UB \
  -o output/
```

**Result:**
- Reads grouped by cell barcode (CB) and UMI (UB)
- Output files:
  - `gene_grouped.tag_CB_counts.tsv`
  - `gene_grouped.tag_UB_counts.tsv`

### Example 3: Multi-column TSV file

Given `barcodes.tsv`:
```
read1   cell1   umi1    batch1
read2   cell2   umi2    batch1
```

```bash
isoquant.py \
  --reference genome.fa \
  --genedb annotation.gtf \
  --fastq sample.fq \
  --data_type nanopore \
  --read_group file:barcodes.tsv:0:1,2,3 \
  -o output/
```

**Result:**
- Reads grouped by 3 columns (cell, UMI, batch)
- Output files:
  - `gene_grouped.file_col1_counts.tsv` (cells)
  - `gene_grouped.file_col2_counts.tsv` (UMIs)
  - `gene_grouped.file_col3_counts.tsv` (batches)

## Automatic file_name Addition

**New Feature (current commit):**

When multiple input files are provided and `--read_group` is specified, `file_name` is **automatically added** if not already present:

```bash
# User runs:
isoquant.py --bam file1.bam file2.bam --read_group tag:CB

# Internally becomes:
args.read_group = ["tag:CB", "file_name"]
```

This ensures technical replicas check works correctly even when users specify custom grouping.

## Technical Replicas Check

**Purpose:** Avoid calling novel transcripts that appear in only one technical replicate (file).

**Implementation:**
1. `args.use_technical_replicas = args.read_group is not None and "file_name" in args.read_group`
2. In GraphBasedModelConstructor:
   - Find index of "file_name" in grouping strategies
   - For each novel transcript candidate, check if supporting reads come from multiple files
   - Skip transcript if all reads come from single file

## Testing

**Key Tests:**
- `tests/console_test.py::test_clean_start` - Tests with `--read_group file:...:0:1`
- `tests/test_file_naming.py` - Tests grouped file naming
- Various GitHub Actions workflows test multi-group functionality

## Important Notes

1. **List Structure:** `args.read_group` is ALWAYS a list (nargs='+') or None, never a string
2. **Backward Compatibility:** Code in `read_groups.py` still handles semicolon-separated strings for internal compatibility
3. **Counter Lists:** When multiple groups exist, many counters become lists instead of single objects
4. **Strategy Names:** File naming depends on `get_grouping_strategy_names()` returning consistent names

## Common Pitfalls

1. **Don't assume read_group is a string** - It's always a list or None
2. **Don't access read_group[0] blindly** - Use the appropriate group_index for your use case
3. **file_name group index** - For technical replicas, must specifically use file_name group, not arbitrary group[0]
4. **Counter structure changes** - Check if counters are lists when multiple groups exist

## Files Modified

Key files involved in multi-group implementation:
- `isoquant.py` - Argument handling, auto-addition of file_name
- `src/read_groups.py` - Grouper classes, multi-grouper, strategy names
- `src/long_read_counter.py` - Multi-group counter support
- `src/dataset_processor.py` - Counter creation for multiple groups
- `src/graph_based_model_construction.py` - Technical replicas check with file_name_group_idx
- `src/file_naming.py` - Strategy-based file naming
- `src/isoform_assignment.py` - read_group as list in ReadAssignment
