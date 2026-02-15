# Single-cell and spatial transcriptomics

**NB! This feature is experimental and is not part of the official IsoQuant release.**

IsoQuant supports single-cell and spatial transcriptomics data from multiple platforms.
When a single-cell or spatial mode is selected, IsoQuant automatically performs
barcode calling and UMI-based PCR deduplication as part of the pipeline.

## Overview

The single-cell/spatial pipeline extends the standard bulk pipeline with these additional steps:

1. **Barcode calling** -- extract cell/spot barcodes and UMIs from raw reads
2. **Standard IsoQuant processing** -- alignment, read-to-isoform assignment
3. **UMI deduplication** -- remove PCR/RT duplicates within each barcode group
4. **Grouped quantification** -- produce per-cell/per-spot count matrices

Note: UMI deduplication relies on read-to-gene assignment. Reads that are not assigned to any gene are discarded.
Hence, novel gene discovery will not be performed in single-cell/spatial mode. 
We recommend using `bulk` mode for novel gene and transcript discovery.

Barcode calling is handled by the built-in [barcode calling module](barcode_calling.md),
which can also be used as a standalone tool.

## Supported platforms

| Mode | Platform | Barcode whitelist files | UMI length | Notes                                                                                |
|------|----------|--------------------|------------|--------------------------------------------------------------------------------------|
| `tenX_v3` | 10x Genomics 3' v3 | 1 | 12 | Single 16bp barcode                                                                  |
| `visium_5prime` | 10x Genomics Visium 5' | 1 | 12 | Same detector as tenX_v3                                                             |
| `visium_hd` | 10x Genomics Visium HD | 2 | 9 | Two barcodes (15/16bp + 14/15bp)                                                     |
| `curio` | Curio Bioscience | 1  | 9 | Double barcode (8bp + 6bp) with linker                                               |
| `stereoseq` | Stereo-seq | 1 | 10 | 25bp barcode, read splitting mode                                                    |
| `stereoseq_nosplit` | Stereo-seq | 1 | 10 | 25bp barcode, no read splitting                                                      |
| `custom_sc` | Any platform | 0 (uses MDF) | varies | User-defined molecule structure via [MDF file](#molecule-definition-file-mdf-format) |

## Quick start examples

10x Genomics single-cell:
```bash
isoquant.py --reference genome.fa --genedb genes.gtf --complete_genedb \
  --fastq reads.fastq.gz --data_type nanopore \
  --mode tenX_v3 --barcode_whitelist 3M-february-2018.txt.gz \
  -o sc_output
```

Stereo-seq spatial:
```bash
isoquant.py --reference genome.fa --genedb genes.gtf --complete_genedb \
  --fastq reads.fastq.gz --data_type nanopore \
  --mode stereoseq --barcode_whitelist barcodes.txt \
  -o stereo_output
```

Custom platform with molecule definition file:
```bash
isoquant.py --reference genome.fa --genedb genes.gtf --complete_genedb \
  --fastq reads.fastq.gz --data_type nanopore \
  --mode custom_sc --molecule molecule_definition.mdf \
  -o custom_output
```

Pre-called barcodes (skip barcode calling):
```bash
isoquant.py --reference genome.fa --genedb genes.gtf --complete_genedb \
  --bam aligned.bam --data_type nanopore \
  --mode tenX_v3 --barcoded_reads barcodes.tsv \
  -o sc_output
```

## Command line options

`--mode` or `-m`

IsoQuant processing mode. Available modes:

* `bulk` -- standard bulk RNA-seq mode (default)
* `tenX_v3` -- 10x Genomics single-cell 3' gene expression
* `curio` -- Curio Bioscience single-cell
* `visium_hd` -- 10x Genomics Visium HD spatial transcriptomics
* `visium_5prime` -- 10x Genomics Visium 5' spatial transcriptomics
* `stereoseq` -- Stereo-seq spatial transcriptomics (BGI), with read splitting
* `stereoseq_nosplit` -- Stereo-seq without read splitting
* `custom_sc` -- custom single-cell/spatial mode using a molecule definition file (MDF)

All modes except `bulk` enable automatic barcode calling and UMI-based deduplication.


`--barcode_whitelist`
Path to file(s) with barcode whitelist(s) for barcode calling.
Required for single-cell/spatial modes unless `--barcoded_reads` is provided.

File should contain one barcode sequence per line. 
More than 1 tab-separated column is allowed, but only the first will be used.
Supports plain text and gzipped files.

Note: barcode calling is performed much better if the whitelist contains a small number of barcodes. 
If you have a subset of barcodes from short-read data, provide them instead of the full whitelist.

The number of whitelist files depends on the mode:

* 1 file: `tenX_v3`, `visium_5prime`, `stereoseq`, `stereoseq_nosplit`, `curio` (combined 14bp barcodes)
* 2 files: `visium_hd` (part 1 and part 2 barcode lists)
* Not needed: `custom_sc` (barcode lists are specified inside the MDF file)

`--barcoded_reads`
Path to TSV file(s) with pre-called barcoded reads.
Format: `read_id<TAB>barcode<TAB>umi` (one read per line).
If provided, IsoQuant skips barcode calling and uses these assignments directly.
More than 3 columns are allowed, but only the first 3 will be used.

Note! IsoQuant does not read barcodes or UMIs from BAM file tags.

`--barcode2spot`
Path to TSV file mapping barcodes to cell types, spatial spots, or other barcode properties.
By default, barcode is in the first column, cell type in the second.
However, you can specify one or more columns via colon symbol (similar to `--read_group`): 
`file.tsv:barcode_column:spot_column(s)` (e.g., `cell_types.tsv:0:1,2,3` for multiple barcode properties).

When `--barcode2spot` is set, `--read_group barcode_spot` will be set automatically
to group counts by cell type, spatial regions, or other provided properties.

`--molecule`
Path to a molecule description format (MDF) file for `custom_sc` mode.
This file defines the structure of the sequencing molecule (barcodes, UMIs, linkers, polyT, cDNA)
and allows IsoQuant to process reads from any single-cell or spatial platform.
See the [MDF format](#molecule-definition-file-mdf-format) section below for details.

## Molecule description format (MDF)

The MDF format allows users to describe the structure of their sequencing molecule
so that IsoQuant can extract barcodes and UMIs from any platform.
The molecule is described in the 3' to 5' direction (primer end first, cDNA last).

An MDF file has two parts:

1. **Header line**: colon-separated list of element names defining the order of elements on the molecule (3' to 5').
2. **Element definitions**: one line per element, tab-separated: `name  type  value  [length]`.

### Element types

| Type | Description                                                                   | Value field              |
|------|-------------------------------------------------------------------------------|--------------------------|
| `CONST` | Constant/known sequence (primer, linker, TSO)                                 | Sequence                 |
| `VAR_FILE` | Variable sequence matched against a whitelist file                            | Path to TSV file         |
| `VAR_LIST` | Variable sequence matched against an inline list                              | Comma-separated sequences |
| `VAR_ANY` | Variable fixed-length sequence extracted as-is                                | Length (integer)         |
| `VAR_ANY_SEPARATOR` | Fixed-length variable separator sequence (not extracted)                      | Length (integer)        |
| `VAR_ANY_NON_T_SEPARATOR` | Fixed-length variable separator sequence without T nucleoties (not extracted) | Length (integer)         |
| `PolyT` | PolyT tail                                                                    | (none)                   |
| `cDNA` | cDNA region                                                                   | (none)                   |

At the moment, only a single `cDNA` and a single `PolyT` are supported. 

Variable elements are expected to have a fixed length (`VAR_FILE` and `VAR_LIST`).
Using variable-length barcodes may result in suboptimal performance.

### Barcode and UMI identification

Elements are identified as barcodes or UMIs by their name prefix (case-insensitive):
- Names starting with `barcode` are treated as barcode elements (e.g., `Barcode`, `barcode1`)
- Names starting with `umi` are treated as UMI elements (e.g., `UMI`, `umi1`)

When multiple barcode (or UMI) elements are present, their sequences are concatenated
in the order they appear.

### Examples

**10x 3' single-cell v3**:
```
R1:Barcode:UMI:PolyT:cDNA:TSO
R1        CONST      CTACACGACGCTCTTCCGATCT
Barcode   VAR_FILE   barcodes.tsv
UMI       VAR_ANY    12
TSO       CONST      CCCATGTACTCTGCGTTGATACCACTGCTT
```

**10x 3' single-cell v3** (inline barcodes):
```
R1:Barcode:UMI:PolyT:cDNA:TSO
R1        CONST      CTACACGACGCTCTTCCGATCT
Barcode   VAR_FILE   AAACCCGGGTTTAAAC,TTTGGGCCCAAATTTG,GGGGAAAACCCCTTTT
UMI       VAR_ANY    12
TSO       CONST      CCCATGTACTCTGCGTTGATACCACTGCTT
```

### Linked elements for multi-part barcodes

When a barcode is split into multiple parts in the molecule, use linked element notation.
All linked parts must be of the same variable type and reference the same whitelist.
Parts are numbered consecutively starting from 1.

**Concatenated** (`|`): parts of the barcode separated by other elements.
Parts are extracted independently, then concatenated and corrected as one sequence against the full whitelist.
Such structure appears in, for example, Curio Bioscience single-cell protocol:

```
PCR_PRIMER:Barcode|1:Linker:Barcode|2:UMI:PolyT:cDNA
PCR_PRIMER  TACACGACGCTCTTCCGATCT
Barcode|1   VAR_FILE   barcodes.tsv   8
Linker      CONST      TCTTCAGCGTTCCCGAGA
Barcode|2   VAR_FILE   barcodes.tsv   6
UMI         VAR_ANY    9
```

**Duplicated** (`/`): multiple redundant copies of the same barcode in the molecule.
Each copy is corrected independently against the whitelist;
a majority vote determines the final barcode.
For example:
```
PCR_PRIMER:Barcode/1:UMI:PolyT:cDNA:Barcode/2:TSO
PCR_PRIMER  TACACGACGCTCTTCCGATCT
Barcode/1   VAR_FILE   barcodes.tsv   16
Barcode/2   VAR_FILE   barcodes.tsv   16
UMI         VAR_ANY    12
TSO       CONST      CCCATGTACTCTGCGTTGATACCACTGCTT
```

## UMI deduplication

All single-cell and spatial modes perform UMI-based PCR deduplication after isoform assignment.
Within each cell barcode and gene, reads with similar UMIs (similarity criteria depends on the UMI length)
are collapsed into a single representative read.

The representative read is selected based on:
1. Unique isoform assignment over ambiguous
2. More exons
3. Longer transcript alignment

### Spot-level UMI deduplication

For spatial transcriptomics at higher resolution than the capture spots
(e.g. Visium HD 2um barcodes mapping to 8um or 16um spots), use
`--barcode2barcode` to deduplicate UMIs at the spot level:

```bash
isoquant.py --reference genome.fa --genedb genes.gtf --complete_genedb \
  --fastq reads.fastq.gz --data_type nanopore \
  --mode visium_hd --barcode_whitelist part1.txt part2.txt \
  --barcode2barcode barcode2spot.tsv:0:1,2 \
  -o visium_output
```

For Visium HD, generate the mapping file from per-part coordinate files:

```bash
python misc/prepare_visium_spot_ids.py part1_to_y.tsv part2_to_x.tsv -o barcode2spot.tsv
```

See `python misc/prepare_visium_spot_ids.py --help` for custom prefix/delimiter options.

## Output

### Count matrices

Single-cell/spatial modes produce grouped count matrices in addition to the standard IsoQuant output.
Use `--counts_format` to control the output format:

* `default` -- automatic selection: matrix format for small numbers of groups (<=100), MTX for larger datasets
* `matrix` -- standard matrix format with genes/transcripts as rows and barcodes as columns
* `mtx` -- Matrix Market (MTX) format compatible with Seurat and Scanpy
* `none` -- no conversion (only internal linear format is produced)

Grouped counts can also be converted after the run using `src/convert_grouped_counts.py`.

### Grouping reads

Use `--read_group` to control how reads are grouped for quantification.
Multiple grouping strategies can be combined (space-separated), producing separate count tables for each.

The most common use-cases for single-cell/spatial data are grouping by barcode property (cell type, spot):

`--read_group barcode_spot` (requires `--barcode2spot`)

or grouping by individual barcode (not recommended for datasets with many barcodes):

`--read_group barcode file_name`


See the [read grouping options](cmd.md#other-input-options) for the full list of grouping strategies.
