# Spot-level Visium HD barcode calling (--barcode2barcode during calling)

Branch `visium_hd_spot_calling` (on top of `visium_hd_split`). A Visium HD barcode `PART1|PART2` identifies a
2 µm spot; part barcodes of neighbouring rows/columns differ by 1–2 edits, so a read with an error in its
barcode can match two neighbours equally well and stay unassigned. Analyses usually use 8 or 16 µm spots,
where such ties often fall into one spot. With `--barcode2barcode` given while IsoQuant calls Visium HD
barcodes, the detector reports spot IDs and assigns a tie when all its candidates share a spot.

## Behaviour and decisions

- Applies when mode is `visium_hd`, `--barcode_whitelist` is given and `--barcode2barcode` is set
  (`resolve_spot_barcode_calling`, `barcode_calling/options.py`). The result is `args.spot_barcode_calling`,
  set before `.params` is saved and read with `getattr(args, "spot_barcode_calling", False)`, so runs started
  before this change resume with the old behaviour; no `CHECKPOINT_FORMAT` bump (stage names unchanged).
- Exactly one spot column (`FILE:0:2`); several columns abort (`parse_spot_mapping_spec`). An explicit
  `--read_group barcode_barcode` aborts: barcodes are spot IDs already, every read would become `NA`.
- Downstream: `barcode_barcode` group not auto-added (the default `barcode` group counts per spot),
  `allinfo` still auto-added, `filter_umis` skips the barcode2barcode rounds (round 0 is spot-level).
- `--barcoded_reads` / `--barcoded_bam` (e.g. CI `SC.Mouse.VisiumHD.barcode2barcode.allinfo`) and other modes:
  the mapping is used at UMI deduplication as before.
- Standalone `isoquant_detect_barcodes.py --barcode2barcode FILE:0:COL` (visium_hd only).
- `--barcode2spot` files must be keyed by spot ID in this mode.

## Implementation

- `SpotMapping` (`barcode_calling/spot_mapping.py`): part → index dicts plus an `int32` matrix of spot indices
  (n1·n2, −1 unmapped) and the spot names. The 11.2M-row slide mapping (1 GB TSV) loads in ~12 s and pickles
  into each spawn worker cheaply; a composite-keyed dict would be ~2.5 GB per worker. Rows whose parts are not
  in the whitelists are skipped with a warning.
- `find_tied_candidates_ssw_var_len` (`common.py`): all candidates with the best score ≥ min_score. The
  existing `find_candidate_with_max_score_ssw_var_len` is untouched (its tie handling is order-dependent,
  ties broken by edit distance), so uniquely called parts are byte-identical; ties are computed only when it
  returns None and a mapping is set.
- `VisiumHDBarcodeDetector._match_part` / `_resolve_barcode` (`callers/tenx.py`), used by both the
  polyT-anchored and the R1-anchored no-polyT paths, hence by the splitting detector and its fallback.
  Ties' positions are used only when they agree, otherwise the ambiguous-part default (expected position).
  Without a mapping `_resolve_barcode` returns `PART1|PART2` exactly as before.
- QA: `misc/assess_barcode_quality.py --barcode2spot FILE[:bc:col]` maps true and detected `PART1|PART2`
  barcodes to spots (detected spot IDs pass through), single and split modes; `run_barcode_test.py` config keys
  `spot_mapping` + `spot_column` feed both the detect command and the QA; `spot_mapping` is a manifest path key.
- CI: barcode level `Mouse.VisiumHD.v1.s1.sep.10M.correct_struct.{8um,16um}` and
  `Mouse.VisiumHD_concat.v1.s1.800K.8um`; pipeline level `SC.Mouse.VisiumHD.spot_calling.allinfo` (simulated
  aligned BAM with `--barcode_whitelist ... --barcode2barcode ...:0:2 --split_molecules false`, i.e. the whole
  route: flag resolution, spot IDs in allinfo, no extra UMI rounds), a step of `SC.Mouse.VisiumHD.barcode2barcode.yml`.

## Measurements (simulated ONT, Stereo-seq NanoSim model, 800K-read CI sets, default `auto` splitting)

Baseline = today's 2 µm calls mapped to spots afterwards (what `--barcode2barcode` at UMI dedup sees);
scored with the same QA. Scoring the baseline at 2 µm via the spot column reproduces the plain barcode QA
exactly (97.58/74.66 non-concatenated, 97.69/73.97 concatenated).

| data | resolution | baseline precision / recall | ties only (margin 0) | default (margin 3, CI baselines) |
|---|---|---|---|---|
| non-concatenated (`correct_struct`) | 8 µm | 99.29 / 75.98 | 99.28 / 76.42 | 99.16 / 78.82 |
| non-concatenated | 16 µm | 99.60 / 76.21 | 99.59 / 76.88 | 99.49 / 79.30 |
| concatenated (2.2 molecules/read) | 8 µm | 99.32 / 75.21 | 99.32 / 75.63 | 99.21 / 77.95 |
| concatenated | 16 µm | 99.64 / 75.45 | 99.63 / 76.07 | — |

Without a mapping, Visium HD output (split and non-split, 100K reads) is byte-identical to `visium_hd_split`.
Cost (800K non-concatenated, 16 threads): 3:59.7 → 4:05 wall clock (mapping load 12 s), peak RSS
2.31 → 2.57 GB. Full-pipeline smoke test (50K concatenated reads, 8 µm): 83,677 molecules with `s_008um_*`
IDs, count matrices keyed by spot, no `barcode_barcode_col*` outputs; `:0:1,2,3` and explicit
`--read_group barcode_barcode` abort.

## Why the gain is small: failure breakdown

`spot/diagnose.py` (in `/abga/work/andreyp/visium_hd_split/`) re-runs the detector steps with the truth on the
9,332 unassigned reads among the first 40K non-concatenated reads:

| category | BC1 | BC2 |
|---|---|---|
| correct | 3,082 | 2,339 |
| true barcode in candidates but scores < min_score (13) | 3,089 | 3,814 |
| true barcode not among the 20 k-mer candidates (mostly score < 13) | 2,620 | 2,885 |
| tie including the truth | 367 | 197 |
| wrong unique barcode | 170 | 94 |

- Ties are ~6% of the failures; where they occur the tied barcodes are 1–2 coordinates apart (BC1 spans:
  1 → 217, 2 → 136, ≥3 → 15), i.e. genuine neighbours, which 8/16 µm spots absorb as designed.
- The dominant loss is the true barcode scoring below min_score in its window. The windows are exactly the
  barcode length; the true BC2 starts 1–3 bases off the window in ~40% of failed reads (polyT start estimate),
  and two sequencing errors in a 15 bp part already drop the score below 13.
- Spot-level calling makes looser matching safer (extra ties between neighbours resolve to one spot), so
  wider windows and/or a lower min_score in spot mode are the natural next step — see below.

## Looser matching in spot mode (window margin, min score)

Spot mode only (`VisiumHDBarcodeDetector.SPOT_WINDOW_MARGIN`, `SPOT_MIN_SCORE`; without a mapping nothing
changes): the BC2 window is widened by the margin on both sides (the bases after BC2 are T-rich and cannot
match a barcode), the BC1 window to the left only (BC2 follows, and part 2 barcodes are part 1 barcodes minus
their first base). When tied candidates disagree on their position the part is assumed at the expected offset.

Sweep on template-disjoint dev sets (100K non-concatenated reads 800,001–900,000 of the 10M `correct_struct`
simulation; 100K concatenated `dev.100K`), scored at 8 µm; baseline = 2 µm calls mapped afterwards:

| margin, min score | non-concat precision / recall | concat precision / recall |
|---|---|---|
| baseline (2 µm) | 99.29 / 75.87 | 99.32 / 75.40 |
| 0, 13 | 99.28 / 76.28 | 99.31 / 75.82 |
| 2, 13 | 99.20 / 78.41 | 99.24 / 77.94 |
| **3, 13 (default)** | **99.18 / 78.69** | **99.21 / 78.20** |
| 5, 13 | 99.15 / 78.86 | 99.16 / 78.38 |
| 0, 12 | 98.70 / 81.06 | 98.72 / 80.32 |
| 3, 12 | 98.30 / 82.89 | 98.34 / 82.11 |
| 3, 11 | 95.76 / 84.60 | 95.88 / 83.64 |

At 16 µm: baseline 99.61 / 76.11, (0,13) 99.60 / 76.73, **(3,13) 99.50 / 79.14**, (5,13) 99.48 / 79.32,
(0,12) 99.11 / 81.81, (3,12) 98.71 / 83.62 (non-concatenated; concatenated alike).

Marginal precision of what each step adds (non-concatenated, 8 µm): margin 3 vs 0: +2,409 correct /
+99 wrong (96%); min score 12 on top: +4,201 / +776 (84%); 11 on top of 12: +1,712 / +2,313, i.e. more wrong
than right, and excessive molecules jump from 4 to 445 — rejected. Margin 3 is the knee (5 adds 0.2% recall for
more errors) and stays within the ≤0.2% precision cost agreed for this feature, so it is the default;
min score 12 (+4.2% recall for −0.9% precision) is left as a decision for the user. Runtime is unchanged.
