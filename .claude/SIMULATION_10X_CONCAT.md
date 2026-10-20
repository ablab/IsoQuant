# Simulating Concatenated 10x Reads for read splitting (--split_molecules) testing

## Overview

The simulation pipeline has two stages:

1. **Template creation** — `simulate_barcoded.py` generates ideal (error-free) read templates in FASTA format, one template per simulated read, with barcodes/UMIs embedded in the correct positions
2. **Error simulation** — NanoSim reads the templates and applies ONT-specific error profiles to produce realistic FASTQ reads

This pipeline is used to validate that `TenXSplittingBarcodeDetector` (and `TenXv2SplittingBarcodeDetector`) correctly recovers barcodes from concatenated molecules.

---

## Stage 1: Template Creation

**Script**: `~/ngs-analysis/barcode_detection/simulate_barcoded.py`

**Template structure** (10x v3 single molecule, in barcode-first "detector" orientation):
```
[R1 (22bp)] [BC (16bp)] [UMI (12bp)] [polyT (30bp)] [RC(cDNA)] [TSO (30bp)]
```
The script builds this and then reverse-complements it, so a template labelled `+` is in mRNA
sense (`RC(TSO) cDNA polyA RC(UMI) RC(BC) RC(R1)`).

For v2: UMI is 10bp; the TSO is the same 30bp sequence as v3.

Visium HD (`--mode visium_hd`): `[R1 (22bp)] [UMI (9bp)] [BC1 BC2 (29-31bp)] [VV] [polyT (30bp)] [RC(cDNA)] [TSO tail (27bp, 10x TSO without CCC)]`;
barcodes are taken from a `BC1|BC2` list and keep the `|` in read names.

The TSO sequence embedded at the end of each molecule segment acts as the molecule boundary marker that `TenXSplittingBarcodeDetector` uses to locate the next molecule's start.

**Concatenation**: `--concatenate_templates` chains multiple molecules end-to-end in a single FASTA entry, with no linker:
```
[mol1] [mol2] ...
```
- Every sub-molecule is a copy of the **same transcript**; only barcode and UMI differ.
- `--concat_probabilities` (default `0.5,0.35,0.1,0.05`): probabilities of 1, 2, 3, ... molecules per template.
- `--strand_switch_prob` (default 0.4): each sub-molecule, **including the first**, is independently
  reverse-complemented **as a whole** (adapters included).
- `--seed`: makes template generation reproducible.

**Ground truth in names**: `READ_<i>_<transcript>_<BC>_<UMI>_<strand>[_<BC>_<UMI>_<strand>...]`; `i` restarts per transcript.

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

For v2 (10bp UMI):
```bash
  --mode 10x_v2
```

**Outputs**:
- `<prefix>.templates.fasta` — template sequences (one per simulated read)
- `<prefix>.template_counts.tsv` — expression table (template_id → count) required by NanoSim. It starts with a
  `#` header line because NanoSim skips the first line of the expression file (files made before this was added
  silently lose their first template).

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
- `--seed` does not make runs reproducible with `-t > 1`: worker processes reseed Python's `random` after fork.
- NanoSim additionally reverse-complements whole reads (the `_R` tag in read names), so both strands occur
  regardless of `--strand_switch_prob`.

Pre-trained models at: `/abga/work/andreyp/simulation/models/`. Low-error Stereo-seq ONT model (3.8% error):
`/abga/work/andreyp/simulation/models/Stereo/StereoSeq.NanoSim.Model/StereoSeq_S2_PC1-CID_L_Exome_LW`
(the older prefix `models/Stereo/StereoSeq_S2_PC1-CID_L_Exome_LW` no longer exists).

**Outputs**:
- `<prefix>_aligned_reads.fasta` — simulated reads with ONT errors
- `<prefix>_aligned_error_profile` — per-read error statistics

---

## Validation

After simulation, run `detect_barcodes.py` directly on the FASTA to check split detection:

```bash
python isoquant_lib/barcode_calling/detect_barcodes.py \
  --fastq /abga/work/andreyp/simulation/<run_name>/data/sim_10x_concat_aligned_reads.fasta \
  --mode tenX_v3 --split_molecules true \
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

The splitting detectors scan both the read and its reverse complement and collect molecules from both,
so concatenated reads whose sub-molecules have mixed orientations are handled per molecule.
Template names from `simulate_barcoded.py` encode the strand of each sub-molecule (see "Ground truth in names").
