# The `--analysis` option

Single primary CLI option that selects which pipeline stages IsoQuant runs,
replacing four scattered per-stage flags. Added on branch `model-read-info-output`
(commit "Add --analysis option and deprecate per-stage flags").

## Values (with short aliases)

| value | alias | drives internal flag(s) | notes |
| --- | --- | --- | --- |
| `quantification` | `quant` | `run_quantification` | gene/transcript counts; also enables polyA/TSS prediction. Needs `--genedb`. |
| `transcript_discovery` | `td` | `no_model_construction=False` | novel transcript models; also enables polyA/TSS prediction. |
| `exon_quantification` | `ex_quant` | `count_exons=True`, `count_intron_retentions=True` | exon + splice-junction + intron-retention counts. Needs `--genedb`. |
| `fusion` | — | `fusion=True` | fusion detection. Needs `--genedb` + a BAM. |
| `velocity` | — | `rna_velocity=True`, `velocity_explicit` | per-cell RNA velocity counts + loom. Needs `--genedb` and cells (non-bulk mode or `--barcode2cellid`). See `.claude/RNA_VELOCITY.md`. |

`--analysis` takes a space-separated list (`nargs='+'`, argparse `choices=ANALYSIS_CHOICES`).
The enum + alias table live in `isoquant_lib/modes.py` (`AnalysisType`, `ANALYSIS_ALIASES`,
`ANALYSIS_CHOICES`).

## Defaults (when `--analysis` is omitted)

- `--genedb` + bulk mode → `{quantification, transcript_discovery}` (unchanged from before).
- no `--genedb` → `{transcript_discovery}` only (de-novo; unchanged).
- single-cell / spatial modes (`mode.needs_barcode_calling()`) → `{quantification}` only.
  **Behavior change**: single-cell no longer builds models by default.
- `velocity` is added on top for `tenX_v3` / `tenX_v2` (`mode.velocity_by_default()`) and
  whenever `--barcode2cellid` is set. In the default case it is skipped per sample when the
  `barcode` grouping exceeds `MAX_DEFAULT_VELOCITY_CELLS` (100k) barcodes.

## Feasibility rules (warn, never abort)

- No `--genedb`: `quantification` / `exon_quantification` / `fusion` / `velocity` are dropped with a
  warning; if nothing feasible remains, the run falls back to `transcript_discovery`.
- `velocity` in `bulk` mode without `--barcode2cellid`: dropped with a warning.
- Single-cell + `transcript_discovery` requested: warns that model construction runs after
  UMI dedup (no novel genes, possibly incomplete — prefer a pseudo-bulk run), then proceeds.

## Resolution

`resolve_analyses(args)` in `isoquant.py` (called from `check_input_params`, right after
`args.mode` is resolved and before `save_params`) translates `--analysis` + the deprecated
flags into canonical booleans on `args`:

```
run_quantification      = quantification in analyses
count_exons             = exon_quantification in analyses  or legacy --count_exons
count_intron_retentions = exon_quantification in analyses  or legacy --count_intron_retentions
fusion                  = fusion in analyses
no_model_construction   = transcript_discovery not in analyses
rna_velocity            = velocity in analyses
velocity_explicit       = rna_velocity and --analysis was given
predict_terminal_sites  = bool(genedb) and (run_quantification or not no_model_construction)
analyses                = the resolved set of AnalysisType
```

Downstream code was **not** renamed — it keeps reading `count_exons` /
`count_intron_retentions` / `fusion` / `no_model_construction`; only the new
`run_quantification` and `predict_terminal_sites` gates were added.

## Deprecated flags (still functional, hidden from `--help`, shown under `--full_help`)

`--count_exons`, `--count_intron_retentions`, `--fusion`, `--no_model_construction`.
Each is captured *before* being overwritten in `resolve_analyses`, applied additively
onto the chosen analyses, and logs a one-line deprecation warning. Legacy
`--count_intron_retentions` keeps its granular effect (IR only), whereas the new
`exon_quantification` value turns on both exon and IR counting.

## Where the gates are consumed

- `isoquant_lib/assignment/assignment_aggregator.py`: ungrouped + grouped gene/transcript
  counters gated on `run_quantification`; polyA/TSS counters gated on `predict_terminal_sites`
  (TSS also `--fl_data`); exon/IR counters on `count_exons` / `count_intron_retentions`.
- `isoquant_lib/dataset_processor.py` `finalize()`: gene/transcript "counts stored in …" log
  lines guarded by `run_quantification`.
- `isoquant.py` `run_pipeline`: `combine_counts` (multi-sample) gated on `run_quantification`;
  fusion trigger block on `args.fusion`.

## Verification (smoke)

```bash
# genedb defaults: counts + models + polyA
python3 isoquant.py -r REF -g ANN.gtf --complete_genedb --bam B.bam -d nanopore -o out
# discovery only: models + polyA, no gene/transcript counts
python3 isoquant.py ... --analysis td
# quant only: counts + polyA, no models
python3 isoquant.py ... --analysis quant
# exon only: exon/SJ/IR counts, nothing else
python3 isoquant.py ... --analysis ex_quant
# no genedb + infeasible request: warns, falls back to transcript_discovery
python3 isoquant.py -r REF --bam B.bam -d nanopore --analysis quant fusion -o out
```

`console_test.py` still exercises the deprecated flags (`--count_exons`,
`--count_intron_retentions`, `--sqanti_output`) and remains byte-compatible.
