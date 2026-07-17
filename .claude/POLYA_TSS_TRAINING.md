# Training the polyA / TSS peak filter

Procedure for regenerating `isoquant_lib/data/model_polya.json` or
`model_tss.json`. The runtime side is documented in
`.claude/POLYA_TSS_DETECTION.md`.

## When you need to retrain

- You changed any of `FEATURE_COLUMNS`, `PEAK_DISTANCE`,
  `PEAK_REL_HEIGHT`, `HISTOGRAM_PAD`, or `ANNOTATION_TOLERANCE` in
  `isoquant_lib/terminal_counter.py`. The classifier's input vector
  has to match what the runtime produces, and there is no schema
  check at inference time.
- You want a per-organism / per-protocol model (e.g. ONT direct-RNA
  vs PacBio IsoSeq cleavage geometry differs enough to benefit from
  a dedicated classifier).
- You're chasing accuracy regressions on a new dataset.

## Two-command workflow

```bash
# 1. Collect labelled per-peak features from a real IsoQuant run.
python isoquant.py \
    --reference REFERENCE.fa.gz \
    --genedb ANNOTATION.gtf.gz --complete_genedb \
    --bam ALIGNED_READS.bam \
    --data_type {nanopore|pacbio_ccs} \
    --fl_data \                              # only when collecting TSS
    -o train_run -t N --prefix run \
    --collect_polya_training polya_peaks.csv \
    --collect_tss_training tss_peaks.csv     # only when collecting TSS

# 2. Train a fresh model.
python misc/train_polya_tss_model.py \
    --features polya_peaks.csv \
    --output isoquant_lib/data/model_polya.json
```

That's it. No source edits. Repeat the second command with
`tss_peaks.csv` → `model_tss.json` for the TSS classifier.

## The hidden flags

`--collect_polya_training` and `--collect_tss_training` are
intentionally hidden from `--help` and `--full_help` (`add_hidden_option`
in `isoquant.py`). They're developer-only and not part of the user
contract. Brief usage comment in `isoquant.py` next to the flag
definitions for future maintainers.

When either flag is given, IsoQuant logs a multi-line WARNING block
at startup ("DEVELOPER MODE: polyA/TSS training-data collection
enabled.") so it's obvious from the log that the run is not a normal
inference run.

## What changes inside IsoQuant in collection mode

- `PolyACounter.TRAINING_ARG = "collect_polya_training"` /
  `TSSCounter.TRAINING_ARG = "collect_tss_training"`. The counter
  reads its flag in `__init__` and sets `self._collecting_training`.
- In collection mode, `dump()` runs the same peak-detection +
  feature-extraction pipeline as inference but skips the XGBoost
  `predict()` call. Every detected peak (and every zero-peak mode
  fallback) becomes a row in the per-chr training CSV next to the
  per-chr prediction TSV, with columns
  `FEATURE_COLUMNS + ['chromosome', 'true_peak']`.
- `true_peak` is computed inline:
  `1 if |peak_location + start − annotated| ≤ ANNOTATION_TOLERANCE`,
  else `0`. Same threshold the runtime uses for Known/Novel
  classification.
- An empty prediction TSV is still emitted (header only) so the rest
  of the pipeline sees the file it expects.
- Grouped counters are skipped entirely in collection mode (only the
  ungrouped counter contributes; grouped output is meaningless when
  the goal is feature collection).
- `dataset_processor.merge_assignments` concatenates per-chr
  fragments into the user-supplied CSV path and removes the per-chr
  files (`merge_files` with `header_lines=1`).

## Training CSV schema

One row per detected peak:

```
var, skew, peak_count, peak_width, entropy, mean_height,
peak_heights, relative_height, chromosome, true_peak
```

- The first eight columns are exactly `FEATURE_COLUMNS` in their
  declared order — the trainer reads them by name, so reordering is
  safe but renaming is not.
- `chromosome` drives the train/test split inside the trainer.
- `true_peak ∈ {0, 1}` is the label.

`skew` can be missing for tiny clusters (the counter sets it to
`None` when fewer than 4 reads or near-zero variance). XGBoost
tolerates NaN; you do not need to filter these rows.

## Sanity-check the dataset

- **Class balance.** The shipped model was trained with
  `scale_pos_weight ≈ 0.297` in
  `misc/train_polya_tss_model.py:DEFAULT_PARAMS`, implying roughly
  3-4 negatives per positive in the original training set. Print
  `df.true_peak.value_counts()` after collection; if your balance is
  very different, tune `scale_pos_weight` to roughly `n_neg /
  n_pos`.
- **Chromosome coverage.** The trainer splits by chromosome and
  refuses to run with fewer than two distinct chromosomes. Aim for
  ≥4 — the more, the more stable the test F1.
- **Dataset size.** A few thousand labelled peaks is usually enough
  to start; tens of thousands produces a noticeably better model.
  The original shipped models came from a ~50 K-row CSV.

## Running the trainer

```
python misc/train_polya_tss_model.py \
    --features polya_peaks.csv \
    --output isoquant_lib/data/model_polya.json
```

Optional flags:
- `--seed INT` (default 42): seed for the chromosome-level split.
- `--test-fraction FLOAT` (default 0.5): fraction of chromosomes
  withheld for testing.

What it does:
- Loads the CSV; validates that `true_peak` and every
  `FEATURE_COLUMNS` column are present.
- Splits chromosomes deterministically into train/test buckets.
- Fits an `XGBClassifier` with `DEFAULT_PARAMS` (same hyperparameters
  used for the shipped models; kept in-file so re-runs are
  reproducible).
- Prints test-set F1 and ROC-AUC.
- Writes the JSON model.

Healthy ballpark on a quality dataset (similar to what produced the
shipped model): F1 ≈ 0.75-0.85, ROC-AUC ≈ 0.90-0.95. Much lower
than that usually means feature columns drifted from what the
runtime now produces, or the labelling threshold is misaligned.

## Hyperparameters

`DEFAULT_PARAMS` in `misc/train_polya_tss_model.py`:

```python
{
    'booster': 'dart',
    'learning_rate': 0.25989,
    'gamma': 4.0875,
    'max_depth': 6,
    'min_child_weight': 5.7230,
    'reg_lambda': 4.8476,
    'reg_alpha': 0.0380,
    'subsample': 0.8810,
    'colsample_bytree': 0.7483,
    'n_estimators': 111,
    'scale_pos_weight': 0.2971,
}
```

These came out of an Optuna sweep on the original training data and
were carried over verbatim. To re-tune, edit `DEFAULT_PARAMS` or
plug in a separate Optuna script.

## Replacing the shipped model

1. Overwrite `isoquant_lib/data/model_polya.json` (or `model_tss.json`)
   with the trainer's output.
2. Re-run the unit suite — `python -m pytest
   isoquant_tests/test_polya_prediction.py -v`. The mocked
   `_StubModel` keeps these tests independent of the file content,
   so they should pass without modification (sanity check the
   counter wiring, not the model accuracy).
3. Run the integration tests — `python -m pytest
   isoquant_tests/console_test.py::test_with_bam_and_polya
   isoquant_tests/console_test.py::test_with_bam_polya_and_fl_data`.
   These only assert non-emptiness, but they at least confirm the
   new JSON loads without an XGBoost version mismatch.
4. Smoke against `isoquant_tests/simple_data/` from a non-repo cwd
   to catch any accidental cwd-dependent regressions:
   ```bash
   mkdir -p /tmp/polya_smoke && cd /tmp/polya_smoke
   python $REPO/isoquant.py \
     --reference $REPO/isoquant_tests/simple_data/chr9.4M.fa.gz \
     --genedb $REPO/isoquant_tests/simple_data/chr9.4M.gtf.gz --complete_genedb \
     --bam $REPO/isoquant_tests/simple_data/chr9.4M.ont.sim.polya.bam \
     --data_type nanopore --fl_data -t 2 -o out --prefix smoke
   wc -l out/smoke/smoke.polyA_prediction.tsv out/smoke/smoke.TSS_prediction.tsv
   ```
   Reference numbers on the currently shipped models: 35 polyA + 11
   TSS predictions. A retrained model will land near those counts
   but not exactly.
5. Commit the new JSON alongside whatever code change motivated the
   retraining.

## Versioning

There is no in-file model version. If you ship a retrained model:
- Bump the `.json` filename if the change is breaking (e.g.
  `model_polya_v2.json`) and update `POLYA_MODEL_PATH` in
  `isoquant_lib/terminal_counter.py`, **OR**
- Leave the filename and note the retraining in `changelog.md` /
  release notes. The previous file is recoverable from git history.

XGBoost JSON is forward-compatible across minor versions (1.x → 2.x
works); cross-major-version migration may need a re-save through
the matching XGBoost release.

## Quick reference: where the moving parts live

| Concern                         | File                                                  |
| ------------------------------- | ----------------------------------------------------- |
| Hidden CLI flags                | `isoquant.py` (search `--collect_polya_training`)     |
| Warning block                   | `isoquant.py:main`                                    |
| Training-dump branch in counter | `isoquant_lib/terminal_counter.py:_dump_training_features` |
| Skipping grouped counters       | `isoquant_lib/assignment_aggregator.py`               |
| Per-chr CSV merge               | `isoquant_lib/dataset_processor.py:merge_assignments` |
| Trainer                         | `misc/train_polya_tss_model.py`                       |
| Per-chr fragment suffix         | `TRAINING_SUFFIX = ".training.csv"` in `terminal_counter.py` |
| Training-CSV column list        | `TRAINING_COLUMNS` in `terminal_counter.py`           |

## Future work

- Move `FEATURE_COLUMNS`-based feature engineering into a shared
  helper so the counter and the trainer cannot drift.
- Optuna sweep helper as a sibling script under `misc/`.
- Optional `--n_jobs 1` on the XGBoost model so the lazy-load
  fork-deadlock workaround can be removed.
- Hash-check the shipped JSON on first load and warn (not fail)
  when it does not match a known good hash — guards against
  accidental retraining without a docs/changelog update.
