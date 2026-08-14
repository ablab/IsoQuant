# CI test configs (snapshot)

Version-controlled copy of the IsoQuant CI test configurations that normally live on the
self-hosted runners at `/abga/work/andreyp/ci_isoquant/data/`.

These files are **not** GitHub Actions workflows. GitHub Actions only reads workflow files
placed directly in `.github/workflows/`; it does not descend into subdirectories, so nothing
here is parsed, validated, or executed as a workflow.

## Why they are here

The configs carry the embedded `baselines:` sections that every CI test is checked against.
They used to exist only on the runner filesystem, with no history — a baseline update left no
record of what changed or why. This directory gives them git history.

## What is here

| Path | Count | Contents |
|------|-------|----------|
| `*.yaml` | 62 | Pipeline test configs, consumed by `isoquant_tests/github/run_pipeline.py` |
| `barcodes/*.yaml` | 24 | Barcode calling test configs, consumed by `isoquant_tests/github/run_barcode_test.py` |
| `fusion/*.yaml` | 2 | Fusion detection test configs |
| `groups/*.yaml` | 2 | IsoQuant `--yaml` dataset descriptors used as input by the grouped tests |

See `.claude/TESTING_SYSTEM.md` for the config format, run types, baseline sections and the
`update_defaults.py` workflow.

## Important: this is a snapshot, not the live copy

CI reads the configs from `CFG_DIR: /abga/work/andreyp/ci_isoquant/data` (set in each
workflow), **not** from this directory. Editing a file here changes nothing about what CI
runs.

That means the two copies drift apart whenever baselines are updated on the runner. After
running `update_defaults.py`, re-sync with:

```bash
SRC=/abga/work/andreyp/ci_isoquant/data
DST=.github/workflows/.internal
cp -p "$SRC"/*.yaml "$DST"/
cp -p "$SRC"/barcodes/*.yaml "$DST"/barcodes/
cp -p "$SRC"/fusion/*.yaml "$DST"/fusion/
cp -p "$SRC"/groups/*.yaml "$DST"/groups/
```

and commit the diff, so the reason for each baseline change is recorded.

Making CI read the configs from the repository directly would remove the drift entirely;
`.claude/TESTING_SYSTEM.md` lists that as a high-priority improvement.

Snapshot taken 2026-08-14 from `/abga/work/andreyp/ci_isoquant/data/`.
