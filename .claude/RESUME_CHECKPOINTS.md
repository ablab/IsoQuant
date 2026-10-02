# Resume checkpoints (`--resume`)

Branch `resume_checkpoints`. Replaces ~15 ad-hoc lock/done files (`_lock`, `_collected`, `_processed`,
`.UMI_filtered.done*`, `.barcodes_done_<i>.tsv`, `.raw_barcodes_done`, `_tagged_bam_done`, ...) with one small
framework, `isoquant_lib/utils/checkpoints.py`. The full design history (three independent reviews, decisions
D1-D9) is in `.claude/RESUME_CHECKPOINTS_PLAN.md`.

## Model

A run is a **linear** sequence of stages; each finished stage leaves a marker file, and `--resume` skips every
stage whose marker exists. Plain markers: no fingerprints (`--resume` forbids changing options anyway).

- `CheckpointStore(root, resume)`: a directory of markers. Holds only a path + flag, so it is picklable and
  passed to `ProcessPoolExecutor` workers. Each marker is its own JSON file
  (`{"isoquant_version", "finished_at", "payload"}`), written to `<name>.done.tmp.<pid>` and `os.replace`d, so
  per-chromosome workers mark concurrently and a marker is never half-written.
- `run_stage(store, name, run, restore=None, cleanup=None, enabled=True, keep_tmp=False)`:
  disabled -> nothing; done -> log `"<name>: done in the previous run, skipping"`, `restore(payload)`;
  otherwise `payload = run()`, then the marker. `cleanup` runs after the marker **and again on skip**
  (removes fragments left by a crash between marker and deletion), never with `--keep_tmp`.
- `mark_stage_done(store, name, payload)`: the marker plus the debug hook; per-chr workers use it directly.
- `marker_name(*parts)`: `/`-joined, every part sanitised like chromosome file names.
- No `invalidate`: a fresh run wipes the stores; on resume a missing marker already means "not done", and
  per-chr child markers (`collect/<chr>`) must survive when the parent stage reruns.

Stores:
- run store `<output>/checkpoints/` (`run_checkpoint_dir`): `cell_barcodes/<prefix>`, `barcodes/<prefix>`,
  `sample/<prefix>` (experiment complete), `combine_counts`, `fusion/<prefix>`;
- sample store `<sample>/aux/checkpoints/` (`SampleData.checkpoint_dir`): every stage below.

Fresh vs resume (`isoquant.py`):
- fresh run: `reset_checkpoints` in `check_and_load_args`, **before `save_params`** -- removes every sample store
  and resets the run store, so a `.params` written by this version always comes with an (empty) run store and no
  old sample marker survives;
- `--resume` without `<output>/checkpoints/` -> exit `RESUME_INCOMPATIBLE` (26): the run predates checkpoints;
- a sample whose directory would be `<output>/checkpoints` is rejected (`_check_sample_dirs_avoid_checkpoints`).

## Stages of one experiment (`DatasetProcessor.process_sample`)

- **Prologue (always runs)**: new `RunSummary`, reset `alignment_stat_counter`, `chr_ids`,
  `use_technical_replicas`, `sample.barcoded_reads = args.barcoded_reads`, tagged-BAM reference list.
1. `read_groups_split`
2. `polya_split` (`--polya_trimmed list:/flnc:`, not `--read_assignments`)
3. `barcode_split` (PCR dedup, not `--barcoded_bam`)
4. `tagged_bam`
5. `collect` + per-chr `collect/<chr>`. A done chromosome is **reloaded**, not skipped (multimapper resolution
   needs all). Payload: alignment stats + `covers_all_reads`, restored by `restore_collection_stats`. `_info` stays
   the source of totals / `all_read_groups`. "To keep these intermediate files for" is logged after the marker
   (CI resume3 kills on it).
6. `umi/ED<d>` (+ per-chr), `umi_bc2bc/col<c>/ED<d>` (+ per-chr). Per-chr allinfo/stats under
   `out_umi_filtered_tmp` (aux), removed by the stage cleanup.
7. `dedup_bam`
- **Derive (always runs)**: `load_read_info`, barcode warning, polyA stats, polyA flags. D1: the flags are derived
  from `original_require_mono{intronic,exonic}_polya` saved in `__init__`, not from the previous experiment's.
8. `construct` + per-chr `construct/<chr>`. Workers close every output before the marker.
- **Post-construct (always runs)**: `report_construction_stats` sums per-chr `_read_stat` / `_transcript_stat`
  (under the saves prefix), prints "Read assignment statistics" / "Transcript model statistics" (TestMode greps
  them), fills the run summary.
9. `merge` (parent marker) with units `merge/<unit>`:
   `read_info`, `assignments`, `corrected_bed`, `sqanti_t2t`, `polya_training`, `tss_training`, `gtf`,
   `extended_gtf`, `model_reads`, and **two units per counter**: `merge/counter/<basename>` (truncate + merge counts
   and stats rows; cleanup keeps `.usable`) and `merge/counter_finalize/<basename>` (`load_usable_fragments` +
   `finalize`: TPM/matrix/MTX/loom; cleanup removes `.usable`). Merging is mechanical; conversion depends on
   options and group count.
   - The merge aggregator is built with `truncate_outputs=False`: no printers, no counter truncates
     (`truncate_output` threaded through every counter class incl. `TerminalCounter`, `RNAVelocityCounter`,
     grouped model counters). Each unit opens/truncates only its own output, so a rerun cannot empty a finished one.
   - Units close outputs before their marker (`close()` on printers, `with` in `merge_counts`).
   - `merge_files` / `merge_counts` take `remove_inputs` (units pass False) and `fragment_dir`.
   - Unit names come from counter file names, which depend on the string pools, hence the parent marker.
10. **Finish**: run summary, compress barcode tables; then `process_all_samples` marks `sample/<prefix>` and
    **after the marker** deletes shared intermediates by glob (`remove_sample_intermediates`, also for skipped
    samples): split barcode / polyA tables, normalized polyA table, `aux/per_chr/` (unless `--keep_tmp`), and
    `out_raw_file_*`, `read_group_file*` (unless `--keep_tmp` or `--read_assignments`).

Per-chr merge fragments live in `<sample>/aux/per_chr/` (`SampleData.chr_fragment_dir`), not the output dir.
`--keep_tmp` keeps every temporary file **except** per-chr tagged/dedup BAM fragments (`merge_bam_files`
deletes them as before, D9).

## Adding a stage

- One function, wrapped in `run_stage` under a unique name (per sample: the sample store).
- In-memory state later code needs -> small JSON payload returned by `run`, re-applied in `restore`. No paths
  (CI copies checkpoint dirs elsewhere), nothing large.
- Fragments only this stage consumes -> idempotent `cleanup`.
- Intermediates shared with later stages -> `remove_sample_intermediates` glob list.
- State needed whether or not stages are skipped -> prologue / derive / post-construct.
- Close every output before returning: the marker is written right after.

## Debugging and testing resume

`ISOQUANT_DEBUG_FAIL=<stage>[:before|:after]` (default `after`) raises `DebugFailure` at that marker, in
`run_stage` and in per-chr workers (`collect/<chr>`, `construct/<chr>`, `umi/ED3/<chr>`). Parsed with the last
`:` only when the suffix is `before`/`after` (chromosome names may contain `:`). An env var, not a CLI flag, so it
is never pickled into `.params`. Typical check: run with the variable, `--resume` without it, diff against a clean
run ignoring `# Command line:` headers.

Unit tests: `isoquant_tests/test_checkpoints.py`, `isoquant_tests/test_resume_merge.py`.

CI resume data (`isoquant_tests/github/generate_resume_test.sh`) must be regenerated on the CI host: old
checkpoints exit with 26. Its stop strings still occur: "Finished processing chromosome SIRV1/SIRV5" (collect and
construct workers), "Gene and transcript records have unequal strands: ..." (construct), "To keep these
intermediate files for" (after collect).

## Output changes in a clean run

- D1: multi-sample runs no longer leak `require_mono*_polya` from one experiment to the next.
- Developer training mode no longer leaves an empty `<p>.polyA_prediction.tsv.training.csv` in the output dir.
- Per-chr fragments never appear in the output dir; chromosome ids containing `/` now merge correctly.
- Log order changes (summary stats are printed after construction, merge logs per unit).
