# Simulating Concatenated 10x Reads for tenX_v3_split / tenX_v2_split Testing

## Overview

The simulation pipeline has two stages:

1. **Template creation** — `simulate_barcoded.py` generates ideal (error-free) read templates in FASTA format, one template per simulated read, with barcodes/UMIs embedded in the correct positions
2. **Error simulation** — NanoSim reads the templates and applies ONT-specific error profiles to produce realistic FASTQ reads

This pipeline is used to validate that `TenXSplittingBarcodeDetector` (and `TenXv2SplittingBarcodeDetector`) correctly recovers barcodes from concatenated molecules.

---

## Stage 1: Template Creation

**Script**: `~/ngs-analysis/barcode_detection/simulate_barcoded.py`

**Template structure** (10x v3 single molecule, forward strand):
```
[R1 (22bp)] [BC (16bp)] [UMI (12bp)] [polyT (30bp)] [RC(cDNA)] [TSO (30bp)]
```

For v2: UMI is 10bp, TSO RC is 25bp.

The TSO sequence embedded at the end of each molecule segment acts as the molecule boundary marker that `TenXSplittingBarcodeDetector` uses to locate the next molecule's start.

**Concatenation**: `--concatenate_templates` chains multiple molecules end-to-end in a single FASTA entry:
```
[mol1: R1+BC+UMI+polyT+RC(cDNA)+TSO] [mol2: R1+BC+UMI+polyT+RC(cDNA)+TSO] ...
```

Concatenation probabilities (hardcoded in script):
- 1 molecule: 50%
- 2 molecules: 35%
- 3 molecules: 10%
- 4 molecules: 5%

Each sub-molecule after the first has ~40% probability of being on the reverse strand (RC of the cDNA segment only — the R1/BC/UMI/polyT/TSO flanks remain in forward orientation relative to the read).

**Command**:
```bash
python ~/ngs-analysis/barcode_detection/simulate_barcoded.py \
  --transcriptome /abga/work/andreyp/simulation/ref/<species>.fa \
  --counts /abga/work/andreyp/simulation/ref/<species>_counts.tsv \
  --template_count 50000 \
  --mode 10x \
  --barcodes /path/to/3M-february-2018.txt \
  -o /abga/work/andreyp/simulation/<run_name>/ref/sim_10x_concat \
  --concatenate_templates
```

For v2 (10bp UMI, 25bp TSO RC):
```bash
  --mode 10x_v2
```

**Outputs**:
- `<prefix>.templates.fasta` — template sequences (one per simulated read)
- `<prefix>.template_counts.tsv` — expression table (template_id → count) required by NanoSim

---

## Stage 2: ONT Error Simulation (NanoSim)

**Environment**: `nanosim` conda env

**Command**:
```bash
conda run -n nanosim python simulator.py transcriptome \
  --ref_t /abga/work/andreyp/simulation/<run_name>/ref/sim_10x_concat.templates.fasta \
  --exp  /abga/work/andreyp/simulation/<run_name>/ref/sim_10x_concat.template_counts.tsv \
  --model_prefix /abga/work/andreyp/simulation/models/<model_name>/training \
  -n 50000 \
  -b guppy \
  -r cDNA_1D \
  --no_model_ir \
  --aligned_only \
  --truncation_mode none \
  -o /abga/work/andreyp/simulation/<run_name>/data/sim_10x_concat
```

Key flags:
- `--aligned_only` — only emit reads that align to the template (no random noise reads)
- `--truncation_mode none` — preserve full-length templates (truncation would cut molecules in the middle and confuse barcode detection)
- `--no_model_ir` — disable intron retention model (not relevant for cDNA)

Pre-trained models at: `/abga/work/andreyp/simulation/models/`

**Outputs**:
- `<prefix>_aligned_reads.fasta` — simulated reads with ONT errors
- `<prefix>_aligned_error_profile` — per-read error statistics

---

## Validation

After simulation, run `detect_barcodes.py` directly on the FASTA to check split detection:

```bash
python isoquant_lib/barcode_calling/detect_barcodes.py \
  --fastq /abga/work/andreyp/simulation/<run_name>/data/sim_10x_concat_aligned_reads.fasta \
  --mode tenX_v3_split \
  --barcode_whitelist /path/to/3M-february-2018.txt \
  --output_sequences \
  -o /tmp/split_test
```

Expected TSV output: multi-line entries (one line per detected molecule per read) for concatenated reads, with correct barcodes matching those embedded in stage 1.

**Key things to check**:
- Fraction of reads yielding ≥2 detected patterns (should match ~50% concatenation rate)
- Barcode recovery rate vs. ground truth barcodes embedded in template names
- TSO detection rate (tso_start column != -1 for interior molecules)

---

## Existing Data

Pre-simulated data for reference/comparison:
- `/abga/work/andreyp/simulation/` — root directory
- `/abga/work/andreyp/simulation/ref/` — reference transcriptomes and barcode whitelists
- `/abga/work/andreyp/simulation/models/` — pre-trained NanoSim models
- `/abga/work/andreyp/simulation/data/` (or per-experiment subdirs) — simulated reads
- Per-experiment `info/info.txt` — command lines used for that run

---

## Notes on Strand Handling

The detector handles strand-switching transparently:
- `TenXSplittingBarcodeDetector.find_barcode_umi()` tries both forward and reverse complement of the full read
- For concatenated reads where sub-molecules switch strand mid-read, the current implementation picks the globally better strand rather than per-molecule strand detection
- This means strand-switched sub-molecules may be missed — a known limitation for the first implementation
- Template names from `simulate_barcoded.py` encode the strand of each sub-molecule (check `info.txt` or script source for exact format)
