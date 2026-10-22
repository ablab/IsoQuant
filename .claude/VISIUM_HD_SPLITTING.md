# Visium HD read splitting (de-concatenation)

Branch `visium_hd_split`. Concatenated ONT Visium HD reads (several cDNA molecules ligated end to end)
are split into one FASTA record per molecule, as already done for 10x v3/v2/visium_5prime and
Stereo-seq. Detector architecture is in `.claude/BARCODE_CALLING.md` § "Visium HD Split Mode"; this
file records the decisions, the simulated data and the measurements behind them.

## Decisions

- `visium_hd` is in `IsoQuantMode.supports_molecule_splitting()`, so the default
  `--split_molecules auto` splits it, like 10x.
- Aligned input + splitting still aborts (`isoquant.py` `_reject_splitting_aligned_input`), also
  with `--barcoded_bam` / `--barcoded_reads` and under `auto`. Deliberate: a BAM means "use my
  alignments", the protocol implies "split", and splitting needs remapping, so neither request is
  dropped silently. BAM-input CI configs pin `--split_molecules false`.
- No `CHECKPOINT_FORMAT` bump: stage names are unchanged, and `--resume` keeps the resolved bool
  stored in `.params`, so an old Visium HD run resumes unsplit.

## Molecule

Barcode-first ("detector") orientation:
`R1(22) UMI(9) BC1(15-16) BC2(14-15) VV polyT(30) RC(cDNA) TSO`; polyT starts 41–43 bp after
the R1 end. The read carries the 10x v3 TSO; the simulator writes its last 27 bp
(`ATGTACTCTGCGTTGATACCACTGCTT`), which is also what the detector searches, so reads with and
without the leading `CCC` match.

## Simulated data

Simulator: `~/ngs-analysis/barcode_detection/simulate_barcoded.py` (ngs-analysis `77bf07a`):
Visium HD UMI 9 bp (was 10), `--concat_probabilities`, `--strand_switch_prob`, `--seed`, a `#`
header in `template_counts.tsv` (NanoSim skips the first line, so the first template used to be
lost), `--umis` fix. All sub-molecules of a template copy one transcript; each, including the first,
is reverse-complemented as a whole with `--strand_switch_prob`.

```
python3 ~/ngs-analysis/barcode_detection/simulate_barcoded.py \
  --transcriptome /abga/work/andreyp/data/reference/mouse/gencode.vM26.basic.transcripts.fasta \
  --counts /abga/work/andreyp/simulation/ref/Mouse.ideal_counts.tsv --template_count 2000000 \
  --mode visium_hd --barcodes /abga/work/andreyp/data/visium/spaceranger/visium_v1_hd.slide1.joint_barcodes.sep.tsv \
  --concatenate_templates --concat_probabilities 0.3,0.35,0.2,0.15 --seed 42 \
  -o /abga/work/andreyp/simulation/ref/visium/Mouse.VisiumHD_concat.v1.s1.2M
~/miniconda3/envs/nanosim/bin/python /abga/work/andreyp/simulation/lrgasp-simulation/src/NanoSim/src/simulator.py transcriptome \
  --ref_t .../Mouse.VisiumHD_concat.v1.s1.2M.templates.fasta --exp .../Mouse.VisiumHD_concat.v1.s1.2M.template_counts.tsv \
  --model_prefix /abga/work/andreyp/simulation/models/Stereo/StereoSeq.NanoSim.Model/StereoSeq_S2_PC1-CID_L_Exome_LW \
  -n 2000000 -b guppy -r cDNA_1D --no_model_ir --aligned_only --truncation_mode none --seed 99 -t 16 --fastq \
  -o /abga/work/andreyp/simulation/data/visium/Mouse.VisiumHD_concat.v1.s1.2M
```
- 2,078,310 templates, mean 2.2 molecules (0.30/0.35/0.20/0.15 for 1–4), 40% reverse-strand
  molecules; Stereo-seq ONT model (3.8% error). NanoSim `--seed` is not reproducible with `-t > 1`.
- CI subset: first 800K reads → `ci_isoquant/data/barcodes/Mouse.VisiumHD_concat.v1.s1.800K.fastq.gz`
  (1,759,901 molecules). Dev subset: 100K later reads sharing no template with the CI subset
  (`/abga/work/andreyp/visium_hd_split/dev.100K.fastq.gz`); smoke: 50K dev reads with short names.
- Read names reach 239 characters (4 molecules), and the split suffix can push them past SAM's 254:
  use this dataset for barcode calling, rename reads before running the full pipeline.
- Ground truth: barcode fields at name positions 3, 6, 9, 12 contain `|`; QA mode `visium_hd_split`
  in `misc/assess_barcode_quality.py` recognises them by the separator.

## Detector design and variant measurements

Three choices were measured on the dev set (100K concatenated reads) and on 100K non-concatenated
reads (`correct_struct` CI data), 16 threads:

| Variant | concat precision | concat recall | non-concat precision | non-concat recall |
|---|---|---|---|---|
| no splitting | 97.74 | 35.24 | 97.59 | 74.81 |
| (a) R1 required (10x-style) | 97.73 | 72.99 | 97.59 | 74.83 |
| (a)+(c) polyT rescue near R1 | 97.72 | 73.03 | 97.59 | 74.83 |
| (b) R1 optional | 97.68 | 74.11 | 97.58 | 74.83 |
| **(b)+(c), chosen** | **97.68** | **74.15** | **97.58** | **74.83** |

- (b) reads both barcode parts backwards from polyT even without R1 — the two-part whitelist match
  keeps random hits rare (added detections 94.8% correct, 3 excessive per 100K). 10x cannot
  do this: a single 16 bp barcode next to a spurious polyT matches too often.
- (c) relaxed polyT search right after the barcodes when R1 is found but polyT is > 50 bp away;
  small but free. Without it `_scan_strand` jumps past such a molecule.
- Per-read fallback to the non-split detector keeps single-molecule reads identical to non-split.

## Final numbers (CI datasets, 800K reads)

| Data | Mode | Precision | Recall |
|---|---|---|---|
| concatenated (2.2 molecules/read) | `--split_molecules false` | 97.71 | 35.09 |
| concatenated | split | 97.69 | 73.97 |
| non-concatenated `correct_struct` | `--split_molecules false` (visium_hd QA) | 97.58 | 74.64 |
| non-concatenated `correct_struct` | split (default `auto`) | 97.58 | 74.66 |

Concatenated recall with splitting is 0.99x the single-molecule recall. TSO found for 95.5% of
detected molecules. Splitting costs 2–3x barcode-calling time (dev, 16 threads: 73 → 160 s
non-concatenated, 125 → 407 s concatenated). Full pipeline smoke test (50K concatenated reads,
GRCm39 + vM36): split FASTA realigned, 85,910 molecules counted, report says "molecules".

Refactor safety: extracting `SplittingScanMixin` and the Visium polyT helper left 10x split
(realistic and 3M whitelists) and Visium non-split output byte-identical to master; Visium
non-split on the full 800K CI reads is identical to master too (its slight drift from the old stored
baselines, e.g. 611,917 vs 611,984 barcoded reads, predates this branch; baselines are now refreshed
together with the UMI fix below).

## UMI without R1 (fixed)

Without R1 the UMI used to be read as the 10 bases before BC1 (`_extract_barcodes_near_polyt`:
`max(0, end - UMI_LEN_VIS)` .. `end`, inclusive), one more than the 9 bp UMI; now `end - UMI_LEN_VIS + 1`.
This also changes non-split output. Measured on the concatenated CI reads (9 bp UMIs), correctly
barcoded molecules without R1 (25,710 of 1,301,799): exact UMI 0.09% → 50.31%, edit distance ≤ 1
55.35% → 63.60%; molecules with R1 unchanged. Barcode metrics are byte-identical. On the older
`correct_struct` simulation, whose UMIs are 10 bp, `umi_within_1` drops slightly (572,352 → 572,013),
because the extra base used to match the extra simulated base; baselines refreshed.

## Known limitations

- No real concatenated Visium HD ONT data was available; false-positive behaviour of the R1-optional
  path on real reads is unmeasured.
- Recall is bounded by barcode matching (~75% on this simulation even without concatenation), not by
  splitting.
