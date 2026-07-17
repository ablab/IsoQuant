# IsoQuant Testing System

This document describes the comprehensive testing system for IsoQuant, including configuration, execution, evaluation, and CI/CD integration.

## Overview

The testing system validates IsoQuant's functionality across multiple dimensions:
- **Read assignment quality** - Accuracy of assigning reads to isoforms
- **Transcript construction** - Quality of novel transcript discovery
- **Quantification accuracy** - TPM/count accuracy for known/novel transcripts and genes
- **Computational performance** - Memory usage, CPU time, and wall-clock time
- **Resume functionality** - Checkpoint/restart capabilities
- **Read grouping** - Multi-group and single-cell/spatial modes

## System Architecture

### Components

1. **Main Test Runner**: `tests/github/run_pipeline.py`
2. **Barcode Test Runner**: `tests/github/run_barcode_test.py`
3. **Baseline Updater**: `tests/github/update_defaults.py`
4. **Config Converter**: `tests/github/cfg2yaml.py` (legacy `.cfg` → `.yaml`)
5. **Config Files**: `/abga/work/andreyp/ci_isoquant/data/*.yaml` (YAML format with embedded baselines)
6. **Evaluation Scripts**: `misc/assess_*.py`, `misc/*_stats.py`, `misc/*_gffcompare.py`
7. **GitHub Actions**: `.github/workflows/*.yml`
8. **Test Data**: `/abga/work/andreyp/ci_isoquant/data/` and `/abga/work/andreyp/data/`
9. **Unit Tests**: `tests/test_ci_config.py` (38 tests for config loading, validation, conversion)

## Test Runner: `tests/github/run_pipeline.py`

### Config Loading

The runner uses extension-based dispatch to load configs:

```python
def load_config(config_file: str) -> tuple[dict, dict]:
    ext = os.path.splitext(config_file)[1].lower()
    if ext in (".yaml", ".yml"):
        return load_yaml_config(config_file)  # returns (config_dict, baselines)
    else:
        return load_tsv_config(config_file), {}  # legacy, no embedded baselines
```

For YAML configs, `load_yaml_config()` extracts the `baselines:` section and returns config values as strings (matching TSV behavior). Assessment functions accept an optional `baselines` parameter and fall back to external files when not provided.

### Main Workflow

The test runner orchestrates the entire testing pipeline:

```python
1. Parse config file (.yaml or legacy .cfg) — auto-detected by extension
2. Run IsoQuant with specified parameters
3. Run quality assessment scripts based on run_type
4. Compare metrics against baselines (embedded in YAML or from external files)
5. Check for expected output files
6. Return exit code (0 = success, negative = failure)
```

### Config File Format

Test configs use YAML format (`.yaml`) with embedded baselines. Legacy TSV format (`.cfg`) is still supported but deprecated. The runner auto-detects format by file extension.

#### YAML Format (current)

```yaml
name: Mouse.ONT_simulated.R10.reduced_db
output: /abga/work/andreyp/ci_isoquant/output
run_type: transcripts,quantification_known,quantification_novel,quantification_genes,performance
genome: /path/to/genome.fa
genedb: /path/to/annotation.gtf
bam: /path/to/aligned.bam
datatype: ont
isoquant_options: '"-t 12 --complete_genedb --force"'
reduced_db: /path/to/reduced_annotation_prefix
reference_tpm: /path/to/ground_truth_tpm.tsv

baselines:
  transcripts:
    full_recall: 84.70
    full_precision: 96.60
    known_recall: 88.60
    known_precision: 98.00
  quantification_ref:
    reported_overlap: 29084
    corr_overlap: 0.8878
  performance:
    max_rss: 8615989248.00
    cpu_time: 20479.00
    clock_time: 2622.00
```

Baselines are embedded under a `baselines:` section with sub-keys:
- `transcripts` — transcript quality (was `.etl` via `etalon` key)
- `assignment` — read assignment quality (was `.etl` via `etalon_assignment` key)
- `quantification_ref`, `quantification_novel`, `quantification_gene` — quantification accuracy (was `.qnt` files)
- `performance` — computational performance (was `.prf` file)
- `allinfo` — UMI filtering quality (was `.etl` via `etalon_allinfo` key)
- `barcode` — barcode calling accuracy (barcode tests only)
- `polya_prediction` — polyA-site prediction accuracy vs a synthetic ground-truth GTF (`misc/assess_polya_prediction.py`)

#### Legacy TSV Format (.cfg) — deprecated

```tsv
name	<test_name>
output	<output_directory>
run_type	<comma_separated_run_types>
genome	<path_to_genome_fasta>
genedb	<path_to_annotation_gtf>
bam	<path_to_aligned_bam>
datatype	<ont|pacbio>
isoquant_options	"<additional_flags>"
```

Old `.cfg` files and their associated `.etl`/`.qnt`/`.prf` files are preserved in `_depr/` subdirectories.

#### Run Types

- `void` - Only run IsoQuant, no evaluation
- `assignment` - Evaluate read assignment quality
- `transcripts` - Evaluate transcript construction (GTF output)
- `quantification_known` - Evaluate known transcript quantification
- `quantification_novel` - Evaluate novel transcript quantification
- `quantification_genes` - Evaluate gene-level quantification
- `performance` - Track computational performance (memory, CPU, time)
- `allinfo` - Evaluate UMI filtering quality (single-cell/spatial modes)
- `polya_prediction` - Evaluate polyA-site prediction accuracy against a ground-truth GTF (requires `reference_polya_gtf` config key)

Multiple run types can be combined: `run_type	transcripts,quantification_known,performance`

#### Special Config Keys

- `resume` - Path to existing IsoQuant output for resume testing
- `label` - Override sample label (required for resume tests)
- `reduced_db` - Prefix to reduced gene database for transcript evaluation
- `reference_tpm` - Ground truth TPM values for transcripts
- `reference_gene_tpm` - Ground truth TPM values for genes
- `check_files` - Space-separated list of expected output file suffixes
- `check_input_files` - Space-separated list of expected output file suffixes (checked as `{label}.{suffix}`)
- `qa_options` - Additional options for quality assessment scripts
- `edit_distance` - Edit distance for UMI filtering (default 3; used by `allinfo` run type)
- `tolerance` - Tolerance fraction for metric comparison (default 0.01 = 1%)

**Legacy-only keys** (not needed in YAML, baselines are inline):
- `etalon` - Path to reference values for transcript quality (.etl file)
- `etalon_assignment` - Path to reference values for assignment quality
- `etalon_quantification_ref`, `etalon_quantification_novel`, `etalon_quantification_gene` - Quantification baselines (.qnt)
- `performance_baseline` - Performance baselines (.prf)
- `etalon_allinfo` - Allinfo quality baselines (.etl file)

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| -2 | Missing required config key |
| -3 | Config file not found |
| -4 | Missing required parameter |
| -5 | Output files missing |
| -6 | One or more quality checks failed |
| -10 | BAM file not found |
| -11 | IsoQuant or QA script failed |
| -12 | Metric check failed |
| -20 | Output GTF not found |
| -21 | Transcript evaluation failed |
| -30 | Output TPM file not found |
| -31 | Reference TPM file not found |
| -33 | Tracking file not found |
| -34 | Quantification evaluation failed |
| -35 | Etalon file not found |
| -36 | Metric not found in output |

## Evaluation Scripts

### 1. Read Assignment Quality: `misc/assess_assignment_quality.py`

**Purpose**: Evaluates how accurately IsoQuant assigns reads to isoforms.

**Inputs**:
- `--tsv` - IsoQuant read assignments (`*.read_assignments.tsv`)
- `--gene_db` - Gene annotation (GTF)
- `--mapping` - Aligned BAM file
- `--fasta` - Original reads (FASTQ/FASTA)

**Outputs**:
- `report.tsv` - TSV with assignment quality metrics

**How it works**:
1. Parses simulated read names to extract ground truth isoform IDs
2. Loads IsoQuant assignments from TSV
3. Compares assigned isoforms vs ground truth
4. Calculates metrics (see below)

**Key Metrics**:
- `correct_assignments` - Reads assigned to correct isoform
- `incorrect_same_gene` - Wrong isoform, correct gene
- `incorrect_other_gene` - Wrong gene
- `not_assigned` - Unassigned reads
- `precision`, `recall`, `f1_score` - Standard classification metrics

### 2. Transcript Quality: `misc/reduced_db_gffcompare.py`

**Purpose**: Evaluates quality of constructed transcript models using gffcompare.

**Inputs**:
- `--genedb` - Prefix to reduced gene database (contains `.expressed.gtf`, `.expressed_kept.gtf`, `.excluded.gtf`)
- `--gtf` - IsoQuant output GTF (`*.transcript_models.gtf`)
- `--tool` - Tool name (isoquant, talon, etc.)

**Outputs**:
- `<tool>.full.stats` - gffcompare stats for all transcripts
- `<tool>.known.stats` - gffcompare stats for known transcripts
- `<tool>.novel.stats` - gffcompare stats for novel transcripts

**How it works**:
1. Separates GTF into full/known/novel categories
2. Runs gffcompare against appropriate reference for each category
3. Parses gffcompare output for recall/precision

**Key Metrics** (from `baselines.transcripts` in YAML or legacy `.etl` files):
- `full_recall` - Recall for all transcripts
- `full_precision` - Precision for all transcripts
- `known_recall` - Recall for known transcripts
- `known_precision` - Precision for known transcripts
- `novel_recall` - Recall for novel transcripts
- `novel_precision` - Precision for novel transcripts

### 3. Quantification Accuracy: `misc/quantification_stats.py`

**Purpose**: Evaluates TPM/count accuracy compared to ground truth.

**Inputs**:
- `--ref_expr` - Ground truth TPM values
- `--tpm` - IsoQuant output TPM (`*.transcript_tpm.tsv`, `*.gene_tpm.tsv`, or `*.discovered_transcript_tpm.tsv`)
- `--tracking` - Optional gffcompare tracking file (for novel transcripts)

**Outputs**:
- `<mode>.quantification.tsv` - TSV with quantification metrics

**How it works**:
1. Loads reference and output TPM values
2. If tracking file provided, maps IsoQuant transcript IDs to reference IDs
3. Compares TPM values in two modes:
   - **overlap** - Only transcripts present in both files
   - **full** - All transcripts (missing ones = 0 TPM)
4. Calculates correlation and match statistics

**Key Metrics** (from `baselines.quantification_*` in YAML or legacy `.qnt` files):
- `reported_overlap` - Number of transcripts reported (overlap mode)
- `corr_overlap` - Pearson correlation (overlap mode)
- `good_matches_overlap` - TPM within ±10% (overlap mode)
- `fair_matches_overlap` - TPM within ±20% (overlap mode)
- `false_reports_overlap` - False positives (overlap mode)
- `missed_overlap` - False negatives (overlap mode)
- Same metrics with `_full` suffix for full mode

### 4. Performance Tracking: `tests/github/performance_counter.py`

**Purpose**: Tracks computational resource usage during IsoQuant execution.

**Inputs**:
- `--cmd_file` - Shell script with IsoQuant command

**Outputs**:
- `stats.tsv` - Performance metrics
- `performance.tsv` - Time-series data (RSS, VMS, CPU% over time)

**How it works**:
1. Launches IsoQuant command as subprocess
2. Monitors process and all children every interval (default 1s)
3. Tracks memory (RSS, VMS) and CPU usage
4. Records max RSS, total CPU time, and wall-clock time

**Key Metrics** (from `baselines.performance` in YAML or legacy `.prf` files):
- `max_rss` - Maximum resident set size (bytes)
- `cpu_time` - Total CPU seconds across all processes
- `clock_time` - Wall-clock time (seconds)

## Baseline Values

### YAML Embedded Baselines (current)

Baselines are embedded directly in the YAML config under the `baselines:` section. Each baseline type is a sub-dictionary of metric name → numeric value:

```yaml
baselines:
  transcripts:
    full_recall: 84.70
    full_precision: 96.60
  quantification_ref:
    reported_overlap: 29084
    corr_overlap: 0.8878
  performance:
    max_rss: 8615989248.00
    cpu_time: 20479.00
    clock_time: 2622.00
  barcode:  # barcode tests only
    precision: 99.89
    recall: 86.82
```

When loaded, all values are converted to strings (matching the TSV behavior) for consistent comparison.

### Legacy Separate Files (deprecated, in `_depr/` folders)

Previously baselines were stored in separate TSV files (`.etl`, `.qnt`, `.prf`) with metric-value pairs, referenced by config keys like `etalon`, `etalon_quantification_ref`, `performance_baseline`, etc.

## Metric Validation

### Tolerance

The `check_value()` function in `run_pipeline.py` validates metrics with configurable tolerance:

```python
def check_value(etalon_value, output_value, name, percent=0.01):
    lower_bound = etalon_value * (1-percent)
    upper_bound = etalon_value * (1+percent)
    # Check if output_value is within bounds
```

**Default tolerance**: ±1% for most metrics
**Performance metrics**: ±20% tolerance (more variability expected)

**Special handling for negative values**: Requires exact match (used for counts that should be exactly 0 or -1)

## GitHub Actions Integration

### Workflow Structure

Each `.github/workflows/*.yml` file:
1. Defines when to run (schedule cron, workflow_dispatch for manual)
2. Sets environment variables:
   - `RUN_NAME` - Test name
   - `LAUNCHER` - Path to `run_pipeline.py`
   - `CFG_DIR` - Path to config directory
   - `BIN_PATH` - Path to external tools (minimap2, samtools, gffcompare)
   - `OUTPUT_BASE` - Base output directory
3. Runs on self-hosted runner with label `isoquant`
4. Executes one or more test steps

### Typical Workflow Step

```yaml
- name: 'Mouse ONT R10 simulated'
  if: always()  # Run even if previous steps failed
  shell: bash
  env:
    STEP_NAME: Mouse.ONT_simulated.R10.reduced_db
  run: |
    export PATH=$PATH:${{env.BIN_PATH}}
    python3 ${{env.LAUNCHER}} ${{env.CFG_DIR}}/${{env.STEP_NAME}}.yaml -o ${{env.OUTPUT_BASE}}
```

### Scheduling

Tests run on different schedules via cron:
- **Weekly**: Sunday 2:00 AM - Read grouping tests
- **Bi-weekly**: Wednesday 3:00 AM - Various simulated data tests
- **Twice weekly**: Tuesday/Friday 6:00 AM - Resume tests
- **Weekly barcode tests**: Protocol-specific (various days) and universal calling custom_sc tests (Mon/Wed/Thu/Fri)
- **Weekly SC pipeline tests**: Saturday — SC.Mouse.10x (5:00), SC.Mouse.VisiumHD.barcode2barcode (5:30), SC.Mouse.10x.barcoded_bam (6:00); Sunday — SC.Human.Curio.custom_sc (6:00)

All barcode tests check for commits in last 7 days before running (skip if no changes).

### Self-Hosted Runners

Tests run on servers with label `isoquant`:
- `dx3` - dx3-523-28534.ad.helsinki.fi
- `dx4` - Available for SSH access

Data and configs are shared via `/abga/work/andreyp/ci_isoquant/`

## Test Categories

### 1. Simulated Data Tests

**Purpose**: Test core functionality with known ground truth

**Examples**:
- `Mouse.ONT_simulated.R10.reduced_db.yaml`
- `Human.PB_simulated.uniform_cov.yaml`
- `SIRVs.simulated_R10.yaml`

**Data characteristics**:
- Simulated reads with isoform IDs embedded in read names
- Known transcript structures
- Known expression levels
- Various sequencing technologies (ONT, PacBio)
- Various annotation completeness (reduced, full, no GTF)

### 2. Real Data Tests

**Purpose**: Test on actual experimental data (SIRV spike-ins)

**Examples**:
- `SIRVs.Set4.R10.opts1.yaml`
- `SIRVs.Set4.R10.ubam.yaml`
- `SIRVs.Set4.real.yaml`

**Data**: Lexogen SIRV Set 4 sequenced with ONT

### 3. Read Grouping Tests

**Purpose**: Test multi-group functionality

**Examples**:
- `GROUP1.TAG.SIRVs.R10.yaml` - Group by BAM tag (RG)
- `GROUP2.TSV.SIRVs.R10.yaml` - Group by TSV file
- `GROUP3.READID.SIRVs.R10.yaml` - Group by read ID pattern
- `GROUP4.FILES.SIRVs.R10.yaml` - Group by file name
- `GROUP5.CSV.SIRVs.R10.yaml` - Group by CSV file
- `GROUP6-9.MTX*.yaml` - Matrix Market format output

**Run type**: Usually `transcripts` with `check_files` to verify grouped outputs

### 4. Resume Tests

**Purpose**: Test checkpoint/resume functionality

**Examples**:
- `RESUME1.SIRVs.R10.yaml` - Resume from assignment step
- `RESUME2.SIRVs.R10.yaml` - Resume from construction step
- `RESUME3.SIRVs.R10.yaml` - Resume from multimapper step
- `RESUME_GR1.SIRVs.R10.yaml` - Resume with read groups

**Special config key**: `resume` points to directory with previous run

### 5. Performance Tests

**Purpose**: Track computational efficiency

**Run type**: `performance`

**Baseline**: `baselines.performance` section in YAML (or legacy `.prf` files)

### 6. Barcode Calling Tests

**Purpose**: Validate barcode detection accuracy for single-cell and spatial modes

**Test runner**: `tests/github/run_barcode_test.py` (separate from main `run_pipeline.py`)

**Config and data location**: `/abga/work/andreyp/ci_isoquant/data/barcodes/`

**Two categories**:

#### Protocol-specific barcode calling tests
Test built-in barcode callers (TenXBarcodeDetector, CurioBarcodeDetector, etc.):
- `Mouse.10x.5k.no_truncation.ONT_Stereo.yaml` — 10x v3 barcode detection
- `Human.Curio.70K.yaml` — Curio barcode detection
- `Mouse.StereoSeq.D04620C2.10M.800K.yaml` — Stereo-seq barcode detection
- `Mouse.VisiumHD.v1.s1.sep.10M.correct_struct.800K.yaml` — Visium HD barcode detection

#### Universal barcode calling tests (custom_sc mode)
Test `UniversalSingleMoleculeExtractor` with MDF files, 2 variants per protocol (small/large whitelist to exercise different k-mer indexers):
- `Mouse.10x.custom_sc.{small,large}.yaml` — 10x via universal calling
- `Human.Curio.custom_sc.{small,large}.yaml` — Curio via universal calling
- `Mouse.VisiumHD.custom_sc.{small,large}.yaml` — Visium HD via universal calling
- `Mouse.StereoSeq.custom_sc.{small,large}.yaml` — Stereo-seq via universal calling

**Config format**: YAML with barcode-specific keys (`mode`, `molecule`, `qa_mode`, `run_type`) and `baselines.barcode` section containing precision/recall/UMI metrics

**QA scripts**: `misc/assess_barcode_quality.py` — computes precision, recall, UMI accuracy against ground truth barcodes embedded in read IDs

See `.claude/BARCODE_CALLING.md` for detailed documentation of universal calling tests and indexer coverage.

### 7. Single-Cell/Spatial Pipeline Tests (allinfo)

**Purpose**: End-to-end pipeline tests for single-cell/spatial modes with UMI filtering quality assessment

**Test runner**: `tests/github/run_pipeline.py` with `run_type allinfo`

**QA script**: `misc/assess_allinfo_quality.py` — validates UMI filtering output (unique gene-barcode pairs, read counts, spliced reads, etc.)

**Tests**:

| Config | Workflow | Mode | Features tested |
|--------|----------|------|-----------------|
| `SC.Mouse.10x.allinfo.yaml` | `SC.Mouse.10x.yml` | tenX_v3 | Standard 10x pipeline with whitelist barcode calling |
| `SC.Mouse.10x.barcoded_bam.allinfo.yaml` | `SC.Mouse.10x.barcoded_bam.yml` | tenX_v3 | `--barcoded_bam --strip_barcode_suffix` (barcodes from BAM CB/UB tags) |
| `SC.Mouse.VisiumHD.barcode2barcode.allinfo.yaml` | `SC.Mouse.VisiumHD.barcode2barcode.yml` | visium_hd | `--barcode2barcode` spot-level UMI dedup + `check_input_files` for barcode_barcode output |
| `SC.Human.Curio.custom_sc.allinfo.yaml` | `SC.Human.Curio.custom_sc.yml` | custom_sc | Universal barcode calling via MDF |

**Data locations** (`/abga/work/andreyp/ci_isoquant/data/barcodes/`):
- `Mouse.10x.5k.ONT_cDNA.R10.4.no_trunc.bam` — Original 10x BAM
- `Mouse.10x.5k.ONT_cDNA.R10.4.no_trunc.tagged.bam` — Same BAM with injected CB/UB tags (created by `misc/inject_bam_tags.py`)
- `Mouse.VisiumHD.800K.aligned.bam` — Pre-aligned Visium HD BAM (800K reads)
- `Mouse.VisiumHD.800K.barcoded_reads.tsv` — Pre-computed barcodes for Visium HD

**Baseline format**: `baselines.allinfo` section in YAML with metrics from `misc/assess_allinfo_quality.py`. Key metrics: `Unique gene-barcodes pairs`, `Total reads saved`, `Spliced reads saved`, `duplicate_triplets`, `allinfo_stats_match`.

**Initial setup for new allinfo test**:
1. Create `.yaml` without a `baselines.allinfo` section
2. Run pipeline manually to verify success
3. Run `update_defaults.py` to populate `baselines.allinfo` from output
4. Verify the baselines look correct

## Test Data Organization

### Server Locations

**Config files** (with embedded baselines):
```
/abga/work/andreyp/ci_isoquant/data/
  *.yaml                   # Test configurations with embedded baselines
  _depr/                   # Deprecated .cfg/.etl/.qnt/.prf files
  ref/                     # Reference genomes and annotations
  mouse/                   # Mouse-specific data
  groups/                  # Grouped read test data
  resume/                  # Resume test checkpoints
  barcodes/
    *.yaml                 # Barcode test configurations
    _depr/                 # Deprecated barcode .cfg/.etl files
```

**Input data**:
```
/abga/work/andreyp/data/
  reference/               # Genome assemblies
  sirvs/                   # SIRV datasets
```

**Output**:
```
/abga/work/andreyp/ci_isoquant/output/<branch_name>/
  <test_name>/             # Per-test output
    <label>/               # IsoQuant output
    report.tsv             # Assignment quality report
    gffcompare/            # Transcript quality results
    *.quantification.tsv   # Quantification results
    stats.tsv              # Performance results
```

## Running Tests Locally

### Single Test

```bash
# On dx3 or dx4 server
cd /path/to/IsoQuant2
export PATH=$PATH:/abga/work/andreyp/ci_isoquant/bin/

# Run a specific pipeline test
python3 tests/github/run_pipeline.py \
  /abga/work/andreyp/ci_isoquant/data/Mouse.ONT_simulated.R10.reduced_db.yaml \
  -o /abga/work/andreyp/ci_isoquant/output/test_run/

# Run a specific barcode test
python3 tests/github/run_barcode_test.py \
  /abga/work/andreyp/ci_isoquant/data/barcodes/Mouse.10x.custom_sc.small.yaml \
  -o /abga/work/andreyp/ci_isoquant/output/test_run/
```

### With Additional Options

```bash
python3 tests/github/run_pipeline.py \
  /abga/work/andreyp/ci_isoquant/data/SIRVs.Set4.R10.opts1.yaml \
  -o /tmp/test_output \
  -a "--threads 4 --debug"
```

### SSH Access

```bash
# From another machine
ssh dx3  # or ssh dx4
cd /home/andreyp/IsoQuant2
# Run tests as above
```

## Updating Baselines

### Using update_defaults.py

IsoQuant includes a utility script to automate baseline updates: `tests/github/update_defaults.py`

**Purpose**: Updates baselines from test run output. For YAML configs, baselines are updated **in-place** in the YAML file. For legacy `.cfg` files, baseline files are copied.

**Usage**:

```bash
# Update baselines from master branch output
python3 tests/github/update_defaults.py \
  /abga/work/andreyp/ci_isoquant/data/Mouse.ONT_simulated.R10.reduced_db.yaml

# Update baselines from specific branch
python3 tests/github/update_defaults.py \
  --branch sc_3.9 \
  /abga/work/andreyp/ci_isoquant/data/Mouse.ONT_simulated.R10.reduced_db.yaml

# Update multiple tests at once
python3 tests/github/update_defaults.py \
  --branch master \
  /abga/work/andreyp/ci_isoquant/data/Mouse.*.yaml

# Update barcode test baselines
python3 tests/github/update_defaults.py \
  --branch master \
  /abga/work/andreyp/ci_isoquant/data/barcodes/Mouse.10x.custom_sc.small.yaml

# Use custom output folder
python3 tests/github/update_defaults.py \
  --output /custom/output/path \
  /abga/work/andreyp/ci_isoquant/data/SIRVs.Set4.R10.yaml
```

**How it works (YAML mode)**:

1. Reads YAML config to find test name and existing baselines
2. Locates output directory: `{output}/{branch}/{test_name}/`
3. For each baseline section present in the YAML, loads the corresponding output file:
   - `baselines.transcripts` ← `gffcompare/new_gtf_etalon.tsv`
   - `baselines.quantification_ref` ← `ref.quantification.tsv`
   - `baselines.quantification_novel` ← `novel.quantification.tsv`
   - `baselines.quantification_gene` ← `gene.quantification.tsv`
   - `baselines.assignment` ← `new_assignment_etalon.tsv`
   - `baselines.performance` ← `new_performance_baseline.tsv`
   - `baselines.allinfo` ← `new_allinfo_etalon.tsv`
   - `baselines.barcode` ← `{run_name}.new_etalon.tsv` (at `{output}/{branch}/` level)
4. Updates the `baselines:` section in the YAML file in-place (preserving all other config keys)

**Workflow for updating baselines**:

```bash
# 1. Run test on current branch
python3 tests/github/run_pipeline.py \
  /abga/work/andreyp/ci_isoquant/data/MyTest.yaml

# 2. If changes are acceptable, update baselines in-place
python3 tests/github/update_defaults.py \
  --branch sc_3.9 \
  /abga/work/andreyp/ci_isoquant/data/MyTest.yaml

# 3. Review the updated YAML
cat /abga/work/andreyp/ci_isoquant/data/MyTest.yaml

# 4. Re-run the test to verify it passes
python3 tests/github/run_pipeline.py \
  /abga/work/andreyp/ci_isoquant/data/MyTest.yaml
```

**Best practices**:

1. **Always review before committing**: Check updated YAML baselines to ensure changes are expected
2. **Document changes**: Commit message should explain why baselines changed
3. **Test after update**: Re-run the test to verify it passes with new baselines
4. **Batch updates carefully**: When updating multiple tests, review each individually
5. **Branch awareness**: Make sure you're using the correct branch's output

## Creating New Tests

### Step 1: Prepare Test Data

1. Create or identify input data (reads, genome, annotation)
2. If using simulated data, ensure read IDs contain ground truth isoform IDs
3. For quantification tests, prepare ground truth TPM files

### Step 2: Create Config File

```bash
# Create .yaml file in /abga/work/andreyp/ci_isoquant/data/
nano /abga/work/andreyp/ci_isoquant/data/MyTest.yaml
```

Example minimal config:
```yaml
name: MyTest
output: /abga/work/andreyp/ci_isoquant/output
run_type: void
genome: /path/to/genome.fa
genedb: /path/to/annotation.gtf
reads: /path/to/reads.fastq
datatype: ont
isoquant_options: '"-t 8 --complete_genedb"'
```

### Step 3: Run Initial Test

```bash
python3 tests/github/run_pipeline.py \
  /abga/work/andreyp/ci_isoquant/data/MyTest.yaml
```

### Step 4: Add Baselines

After successful run, the test generates `new_*_etalon.tsv` files in the output directory. Add empty baseline sections to your YAML, then use `update_defaults.py` to populate them:

```bash
# First, add empty baseline sections to YAML for the run types you want to track.
# For example, add:
#   baselines:
#     transcripts: {}
#     quantification_ref: {}
#     performance: {}

# Then update baselines from output
python3 tests/github/update_defaults.py \
  --branch <branch> \
  /abga/work/andreyp/ci_isoquant/data/MyTest.yaml
```

This updates the `baselines:` section in-place in the YAML file.

### Step 5: Verify Updated Config

The YAML file should now contain embedded baselines:

```yaml
name: MyTest
output: /abga/work/andreyp/ci_isoquant/output
run_type: transcripts,quantification_known,quantification_novel,performance
genome: /path/to/genome.fa
genedb: /path/to/annotation.gtf
reads: /path/to/reads.fastq
datatype: ont
isoquant_options: '"-t 8 --complete_genedb"'
reduced_db: /path/to/reduced_annotation_prefix
reference_tpm: /path/to/ground_truth_tpm.tsv

baselines:
  transcripts:
    full_recall: 84.70
    full_precision: 96.60
    known_recall: 88.60
    known_precision: 98.00
  quantification_ref:
    reported_overlap: 29084
    corr_overlap: 0.8878
  performance:
    max_rss: 8615989248.00
    cpu_time: 20479.00
    clock_time: 2622.00
```

### Step 6: Create GitHub Action Workflow

```bash
nano .github/workflows/MyTest.yml
```

Example:
```yaml
name: My Custom Test

on:
  workflow_dispatch:
  schedule:
  - cron: '0 4 * * 1'  # Monday 4:00 AM

env:
  RUN_NAME: MyTest
  LAUNCHER: ${{github.workspace}}/tests/github/run_pipeline.py
  CFG_DIR: /abga/work/andreyp/ci_isoquant/data
  BIN_PATH: /abga/work/andreyp/ci_isoquant/bin/
  OUTPUT_BASE: /abga/work/andreyp/ci_isoquant/output/${{github.ref_name}}/

concurrency:
  group: ${{github.workflow}}
  cancel-in-progress: false

jobs:
  launch-runner:
    runs-on:
      labels: [isoquant]
    name: 'Running IsoQuant and QC'

    steps:
      - name: 'Cleanup'
        run: >
          set -e &&
          shopt -s dotglob &&
          rm -rf *

      - name: 'Checkout'
        uses: actions/checkout@v3
        with:
          fetch-depth: 1

      - name: 'MyTest'
        if: always()
        shell: bash
        env:
          STEP_NAME: MyTest
        run: |
          export PATH=$PATH:${{env.BIN_PATH}}
          python3 ${{env.LAUNCHER}} ${{env.CFG_DIR}}/${{env.STEP_NAME}}.yaml -o ${{env.OUTPUT_BASE}}
```

### Step 7: Test and Commit

```bash
# Test workflow manually via GitHub Actions UI
# or trigger via workflow_dispatch

# Once verified, commit
git add .github/workflows/MyTest.yml
git commit -m "Add MyTest to CI pipeline"
git push
```

## Troubleshooting

### Test Failures

1. **Check logs**: Look at output in GitHub Actions or local terminal
2. **Review metrics**: Compare output metrics vs baselines in output directory
3. **Inspect outputs**: Check if expected files were generated
4. **Verify data**: Ensure input data is accessible at specified paths
5. **Check tools**: Verify minimap2, samtools, gffcompare are in PATH

### Common Issues

**Metric slightly out of tolerance**:
- Review if legitimate change or regression
- If intentional improvement, update baseline files
- If regression, investigate and fix

**Missing output files**:
- Check IsoQuant logs for errors
- Verify isoquant_options are correct
- Check disk space on server

**Performance degradation**:
- Review code changes for algorithmic changes
- Check if parallelization settings changed
- Consider if acceptable trade-off for accuracy

**Path issues**:
- Ensure BIN_PATH includes all required tools
- Verify data paths are absolute, not relative
- Check file permissions on server

## Best Practices

1. **Always test locally before pushing**: Run tests with new configs locally first
2. **Use descriptive names**: Test names should indicate what's being tested
3. **Document baselines**: When updating baselines, note why in commit message
4. **Conservative tolerances**: Keep tolerances tight to catch regressions early
5. **Comprehensive coverage**: Test different modes, data types, and options
6. **Regular updates**: Update baselines when intentional improvements are made
7. **Monitor trends**: Watch for gradual performance degradation over time

## Areas for Improvement

### Current Limitations and Pain Points

#### 1. Configuration Management (partially addressed)

**Addressed by YAML migration:**
- Replaced ad-hoc TSV `.cfg` format with standard YAML — better readability, tooling support
- Eliminated separate baseline files (`.etl`, `.qnt`, `.prf`) — baselines embedded in YAML
- Single file per test instead of multiple — easier to track and manage
- Conversion tool `tests/github/cfg2yaml.py` available for any remaining migrations
- Unit tests in `tests/test_ci_config.py` validate config loading and conversion

**Remaining issues:**
- Config files still stored outside repository in `/abga/work/andreyp/ci_isoquant/data/`
- No schema validation for YAML configs
- Could benefit from moving configs into repo for version control

#### 2. Baseline Management (partially addressed)

**Addressed by YAML migration:**
- Baselines now embedded in YAML — single file per test, no sync issues between config and baselines
- `update_defaults.py` now updates baselines in-place in YAML files (not file copies)
- Supports both pipeline and barcode baseline updates
- Simpler workflow: run test → run update_defaults.py → review YAML diff

**Remaining issues:**
- Still requires manual review to decide when to update baselines
- No diff preview before updating
- No validation that metrics actually improved
- Baselines still outside repo (no git history)

**Potential improvements:**
- Add interactive mode or diff preview to `update_defaults.py`
- Store YAML configs in repo for version-controlled baseline changes
- Add validation: warn if precision/recall decreased while updating

#### 3. Test Execution and Parallelization

**Issues:**
- GitHub Actions workflows run tests sequentially, even when independent
- Each workflow defined separately with lots of duplication
- No easy way to run subset of tests locally
- Tests take long time to complete (some hours)
- Resource contention when multiple workflows run

**Potential Improvements:**
- Implement parallel test execution in workflows (GitHub Actions matrix strategy)
- Create test groups by dependency and runtime
- Develop local test runner that can run selected tests
- Add "smoke test" suite with small/fast tests for quick validation
- Implement test result caching to skip unchanged tests

#### 4. Error Reporting and Debugging

**Issues:**
- Exit codes are generic (-2, -11, etc.) and don't provide context
- Have to dig through logs to find actual failure reason
- No summary of which metrics failed and by how much
- Difficult to compare failed run against baseline visually
- No automatic notification system for failures

**Potential Improvements:**
- Add detailed error messages with context (which metric, expected vs actual, file paths)
- Generate HTML report showing side-by-side comparison of metrics
- Implement Slack/email notifications for test failures with summary
- Add "test summary" page aggregating all test results
- Create failure categorization (data issue, code regression, baseline drift, etc.)

#### 5. Metrics and Monitoring

**Issues:**
- No historical tracking of metrics over time
- Can't see if performance is gradually degrading
- No visualization of metric trends
- Tolerance is hardcoded (1% for most, 20% for performance)
- No statistical analysis of metric variability

**Potential Improvements:**
- Create metrics database to track all test runs over time
- Build dashboard showing metric trends (line charts, regression detection)
- Implement adaptive tolerance based on historical variability
- Add statistical tests (e.g., detect significant changes)
- Track "flaky" tests that occasionally fail due to variability
- Generate performance regression reports

#### 6. Test Coverage

**Issues:**
- No clear mapping of which features are tested by which tests
- Hard to know if new feature needs new test
- No coverage metrics for different IsoQuant modes/options
- Some edge cases may not be covered

**Potential Improvements:**
- Create test coverage matrix (features × tests)
- Add code coverage tracking to CI (pytest-cov integration)
- Document which options are tested and which aren't
- Implement automatic detection of untested code paths
- Add property-based testing for complex algorithms

#### 7. Data Management

**Issues:**
- Test data scattered across server in different locations
- Large data files not version controlled
- Hard to reproduce tests without server access
- Data duplication across tests
- No checksums to verify data integrity

**Potential Improvements:**
- Centralize test data with clear organization
- Implement data versioning (e.g., DVC, git-lfs)
- Create "test data catalog" documenting all datasets
- Add data validation checks (checksums, format validation)
- Develop minimal test datasets that can run locally without server
- Share common reference files across tests

#### 8. Documentation and Discoverability

**Issues:**
- Hard to find which test covers specific scenario
- No index of all available tests with descriptions
- Config format not well documented in code
- Learning curve for adding new tests is steep

**Potential Improvements:**
- Auto-generate test catalog from configs with descriptions
- Add docstrings to evaluation scripts
- Create interactive guide for adding new tests
- Implement `--list-tests` option to discover available tests
- Add examples for common test patterns

#### 9. Integration with Development Workflow

**Issues:**
- Tests run on schedule, not on PRs
- No fast feedback loop for developers
- Can't easily test changes before merge
- Unit tests (pytest) separate from integration tests (CI)

**Potential Improvements:**
- Add PR-triggered tests (subset of fast/critical tests)
- Create pre-commit hooks to run local validation
- Implement test selection based on changed files
- Link unit test failures to relevant integration tests
- Add "test plan" generation based on code changes

#### 10. Tool Dependencies

**Issues:**
- External tools (minimap2, samtools, gffcompare) managed manually
- Version inconsistencies across runners
- No automatic installation or verification
- PATH manipulation required

**Potential Improvements:**
- Use containerization (Docker) for consistent environments
- Implement dependency management (conda environment)
- Add version checking and automatic installation
- Create setup script for test environment
- Pin tool versions in requirements

### Priority Recommendations

**High Priority (Quick Wins):**
1. Move config files into repository for version control
2. Add detailed error messages with metric comparisons
3. Implement parallel test execution in GitHub Actions
4. ~~Standardize config format~~ — **Done**: migrated to YAML with embedded baselines

**Medium Priority (Important but Complex):**
1. Build metrics tracking database and visualization
2. Implement baseline approval workflow
3. Add PR-triggered smoke tests
4. Containerize test environment

**Low Priority (Nice to Have):**
1. Property-based testing for algorithms
2. Auto-generated test documentation
3. Advanced statistical analysis of metrics
4. Local minimal test datasets

### Technical Debt

1. **Hardcoded paths**: Many absolute paths in configs should use environment variables
2. **Code duplication**: Workflow YAML files have lots of repeated boilerplate
3. **Error handling**: Many silent failures that should raise warnings
4. ~~**Legacy formats**: Some configs use old format (should migrate to consistent schema)~~ — **Addressed**: All configs migrated to YAML; legacy `.cfg` files in `_depr/` folders
5. **Performance counter**: Could use more robust process tracking
6. **Magic numbers**: Tolerance percentages (0.01, 0.2) should be configurable per metric
7. **Code duplication across runners**: `load_tsv_config()`, `fix_path()`, `load_yaml_config()`, `load_config()` are duplicated in `run_pipeline.py`, `run_barcode_test.py`, and `update_defaults.py` — could be extracted into a shared module

## Config Format Migration Tools

### cfg2yaml.py Converter

`tests/github/cfg2yaml.py` converts legacy `.cfg` files (with associated `.etl`/`.qnt`/`.prf` baseline files) into single YAML files.

**Usage:**
```bash
# Convert pipeline configs
python3 tests/github/cfg2yaml.py /abga/work/andreyp/ci_isoquant/data/MyTest.cfg

# Convert barcode configs (auto-detected or forced)
python3 tests/github/cfg2yaml.py /abga/work/andreyp/ci_isoquant/data/barcodes/Mouse.10x.cfg
python3 tests/github/cfg2yaml.py --barcode /path/to/config.cfg

# Dry run (preview without writing)
python3 tests/github/cfg2yaml.py --dry-run /path/to/config.cfg

# Convert multiple configs at once
python3 tests/github/cfg2yaml.py /abga/work/andreyp/ci_isoquant/data/*.cfg
```

**Auto-detection**: Barcode configs are auto-detected by presence of `mode`/`molecule` keys without `genome`/`genedb` keys. Use `--barcode` to force barcode mode.

**Baseline key mapping:**
- `etalon` → `baselines.transcripts` (pipeline) or `baselines.barcode` (barcode)
- `etalon_assignment` → `baselines.assignment`
- `etalon_quantification_ref` → `baselines.quantification_ref`
- `etalon_quantification_novel` → `baselines.quantification_novel`
- `etalon_quantification_gene` → `baselines.quantification_gene`
- `performance_baseline` → `baselines.performance`
- `etalon_allinfo` → `baselines.allinfo`

## Unit Tests for CI Config System

`tests/test_ci_config.py` contains 38 unit tests across 9 test classes:

| Class | Tests | What it covers |
|-------|-------|----------------|
| `TestLoadTsvConfig` | 4 | TSV config loading (basic, comments, blank lines, missing file) |
| `TestLoadYamlConfig` | 6 | YAML config loading (baselines extraction, string conversion, edge cases) |
| `TestLoadConfigDispatch` | 4 | Extension-based dispatch (.yaml/.yml/.cfg) |
| `TestCheckValue` | 7 | Tolerance checking in metric validation |
| `TestDetectBarcodeMode` | 4 | Pipeline vs barcode config auto-detection |
| `TestCfg2Yaml` | 4 | Config conversion (pipeline, barcode, force mode, missing etalon) |
| `TestLoadBaselineFile` | 3 | Baseline TSV file loading |
| `TestUpdateYamlBaselines` | 4 | In-place YAML baseline updating (pipeline, barcode, missing output) |
| `TestYamlTsvRoundTrip` | 2 | Equivalence between converted YAML and original TSV configs |

Run with:
```bash
pytest tests/test_ci_config.py -v
```
