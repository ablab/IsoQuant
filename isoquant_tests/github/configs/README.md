# CI test configs

The configurations every IsoQuant CI test runs from. **CI reads these files**, not a
copy on the runners: each workflow sets

```yaml
CFG_DIR: ${{github.workspace}}/isoquant_tests/github/configs
```

so the checked-out repository is the single source of truth. Editing a config here
changes what CI does on the next run.

## What is here

| Path | Count | Consumed by |
|------|-------|-------------|
| `*.yaml` | 62 | `run_pipeline.py` |
| `barcodes/*.yaml` | 23 | `run_barcode_test.py` |
| `fusion/*.yaml` | 2 | `run_pipeline.py` (`run_type: fusion`) |
| `groups/*.yaml` | 2 | IsoQuant `--yaml` dataset descriptors, referenced by the grouped tests |
| `data_manifest.tsv` | 1 | `data_manifest.py` - size + SHA-256 of every input file |

See `.claude/TESTING_SYSTEM.md` for the config format, run types and baseline sections.

## Updating baselines

Baselines live in the `baselines:` section of each config. After a run whose results
you accept:

```bash
python3 isoquant_tests/github/update_defaults.py --branch master \
    isoquant_tests/github/configs/MyTest.yaml
git add isoquant_tests/github/configs/MyTest.yaml
git commit -m "update MyTest baselines: <why they moved>"
```

Because the configs are version-controlled and CI reads them directly, every baseline
change is a reviewable diff with a reason attached. This replaced an earlier
arrangement where the configs lived only on the runners and drifted from the
repository copy, which happened three times in the first four days.

## Paths must stay absolute

All path values are absolute. The runners resolve relative paths against the config
file's own directory (`fix_path()`), so a relative path would break the moment the
config directory moves. Keep new entries absolute.

## Input data

The large inputs (genomes, BAMs, annotations) stay on shared storage and are *not* in
git. `data_manifest.tsv` records each one's size and SHA-256 so a baseline can be tied
to specific input bytes:

```bash
# after adding or replacing an input file
python3 isoquant_tests/github/data_manifest.py generate

# check inputs are unchanged (size only; --deep re-hashes)
python3 isoquant_tests/github/data_manifest.py verify --deep
```

`run_pipeline.py` runs the fast check automatically before each test; disable with
`--verify_inputs off`.
