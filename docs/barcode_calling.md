# Barcode calling

IsoQuant includes a built-in barcode calling module (`detect_barcodes.py`) that extracts
cell barcodes and UMIs from raw long reads. This module supports multiple single-cell and
spatial transcriptomics platforms and can also be used as a standalone tool.

See [single-cell and spatial transcriptomics](single_cell.md) for supported platforms,
pipeline integration, MDF format, and UMI deduplication.

## How it works

The barcode calling module processes each read by:

1. Searching for known constant sequences (linkers, primers, TSO) on the read
2. Locating the polyT tail to determine read orientation
3. Extracting the barcode region based on its expected position relative to the anchoring constant sequences
4. Matching the extracted barcode against a whitelist
5. Extracting the UMI sequence from its expected position

Each read is assigned:
- A **barcode** (cell identity) corrected against the whitelist, or `*` if no match is found
- A **UMI** (unique molecular identifier), or `*` if not detected
- A **strand** orientation (`+`, `-`, or `.` if unknown)
- Platform-specific features (polyT position, linker positions, etc.)

## Integration with IsoQuant

When running IsoQuant with a single-cell or spatial mode (`--mode tenX_v3`, etc.),
barcode calling is performed automatically as the first pipeline step.
The barcode detection output is a TSV file with one line per read, which is then used
throughout the rest of the pipeline for grouping and UMI deduplication.

If barcodes have already been called externally (e.g., by Cell Ranger),
use `--barcoded_reads` to skip the barcode calling step.

## Standalone usage

`detect_barcodes.py` can be run independently from the IsoQuant pipeline:

```bash
python detect_barcodes.py \
  --input reads.fastq.gz \
  --barcodes barcode_whitelist.txt \
  --mode tenX_v3 \
  --output output_prefix \
  --threads 16
```

### Command line options

`--input` or `-i` (required)

One or more input read files in FASTQ, FASTA, BAM, or SAM format. Gzipped FASTQ/FASTA files are supported.

`--output` or `-o` (required)

Output prefix. The barcode calling results will be written to `<prefix>.barcoded_reads.tsv`.
When multiple input files are provided, outputs are numbered: `<prefix>_0.barcoded_reads.tsv`, etc.

`--barcodes` or `-b`

One or more barcode whitelist files. The number of files depends on the mode:

* 1 file: `tenX_v3`, `visium_5prime`, `stereoseq`, `stereoseq_nosplit`
* 1 or 2 files: `curio`
* 2 files: `visium_hd`
* Not needed for `custom_sc` (barcodes defined in MDF file)

`--mode`

Barcode calling mode. Available modes: `tenX_v3`, `curio`, `stereoseq`, `stereoseq_nosplit`,
`visium_5prime`, `visium_hd`, `custom_sc`. Default: `stereoseq`.

`--molecule`

Path to a molecule definition file (MDF) for `custom_sc` mode.
See [MDF format](single_cell.md#molecule-definition-file-mdf-format) for the format specification.

`--threads` or `-t`

Number of threads for parallel processing (default: 16).

`--min_score`

Override the minimum alignment score for barcode matching.
The scoring system uses +1 for match, -1 for mismatch, -1 for gap open, -1 for gap extension.
By default, the minimum score is set automatically based on barcode length.

`--tmp_dir`

Folder for temporary files during parallel processing.

### Example commands

10x Genomics:
```bash
python detect_barcodes.py \
  -i reads.fastq.gz \
  -b 3M-february-2018.txt.gz \
  --mode tenX_v3 \
  -o results/my_sample \
  -t 16
```

Stereo-seq:
```bash
python detect_barcodes.py \
  -i reads.fastq.gz \
  -b stereo_barcodes.txt \
  --mode stereoseq \
  -o results/stereo_sample \
  -t 16
```

Curio with two barcode files:
```bash
python detect_barcodes.py \
  -i reads.fastq.gz \
  -b left_barcodes.txt right_barcodes.txt \
  --mode curio \
  -o results/curio_sample
```

Custom molecule:
```bash
python detect_barcodes.py \
  -i reads.fastq.gz \
  --mode custom_sc \
  --molecule my_platform.mdf \
  -o results/custom_sample
```

## Output format

The main output is a TSV file (`*.barcoded_reads.tsv`) with a header line followed by one line per read.
The columns depend on the platform mode.

### Common columns (all modes)

| Column | Description |
|--------|-------------|
| read_id | Read identifier |
| barcode | Corrected barcode sequence, or `*` if not found |
| UMI | UMI sequence, or `*` if not detected |
| BC_score | Barcode alignment score (-1 if not matched) |
| strand | Read orientation: `+`, `-`, or `.` |

### Additional columns by mode

**10x Genomics** (`tenX_v3`, `visium_5prime`, `visium_hd`):

| Column | Description |
|--------|-------------|
| polyT | Position of polyT tail start (-1 if not found) |
| R1 | Position of R1 primer end (-1 if not found) |

**Curio** (`curio`):

| Column | Description |
|--------|-------------|
| polyT | PolyT tail position |
| primer | PCR primer position |
| linker_start | Linker start position |
| linker_end | Linker end position |

**Stereo-seq** (`stereoseq_nosplit`, `stereoseq`):

| Column | Description |
|--------|-------------|
| polyT | PolyT tail position |
| primer | Primer position |
| linker_start | Linker start position |
| linker_end | Linker end position |
| TSO | TSO position (-1 if not found) |

**Custom** (`custom_sc`): columns depend on the elements defined in the MDF file.
Each variable element produces a column with its detected sequence and score.

### Statistics file

A statistics file (`*.barcoded_reads.tsv.stats`) is generated alongside the TSV output,
reporting the total number of reads processed, barcodes detected, and UMIs found.

## Supported platforms

### 10x Genomics (`tenX_v3`, `visium_5prime`)

Molecule structure (3' to 5'):
```
[R1 primer]---[Barcode (16bp)]---[UMI (12bp)]---[PolyT]---[cDNA]---[TSO]
```

Requires 1 barcode whitelist file (e.g., the 10x `3M-february-2018.txt.gz`).

### Visium HD (`visium_hd`)

Molecule structure (3' to 5'):
```
[R1 primer]---[BC_part1 (16bp)]---[separator]---[BC_part2 (15bp)]---[UMI (9bp)]---[PolyT]---[cDNA]---[TSO]
```

Requires 2 barcode whitelist files (part 1 and part 2).
The final barcode is the concatenation of both parts.

### Curio Bioscience (`curio`)

Molecule structure (3' to 5'):
```
[PCR primer]---[Left BC (8bp)]---[Linker]---[Right BC (6bp)]---[UMI (9bp)]---[PolyT]---[cDNA]
```

The linker sequence anchors the barcode detection.
Accepts either 1 combined whitelist file (14bp barcodes) or 2 files (8bp left + 6bp right).

### Stereo-seq (`stereoseq`, `stereoseq_nosplit`)

Molecule structure (3' to 5'):
```
[Primer]---[Barcode (25bp)]---[Linker]---[UMI (10bp)]---[PolyT]---[cDNA]---[TSO]
```

- `stereoseq` mode: splits concatenated reads at TSO boundaries, producing a new FASTA with individual subreads
- `stereoseq_nosplit` mode: processes reads without splitting

### Custom platform (`custom_sc`)

Uses a molecule definition file (MDF) to describe arbitrary molecule structures.
See [MDF format](single_cell.md#molecule-definition-file-mdf-format).

The universal extractor:
1. Tries both forward and reverse complement orientations
2. Detects polyT tail to determine strand
3. Finds constant elements (primers, linkers) on the read
4. Extracts variable elements (barcodes, UMIs) from their expected positions
5. Corrects barcodes against whitelists when provided (`VAR_FILE` or `VAR_LIST`)
6. Handles linked elements: concatenated (`|`) parts are joined then corrected; duplicated (`/`) parts use majority vote
