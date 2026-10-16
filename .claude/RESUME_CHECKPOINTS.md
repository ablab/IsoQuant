# Resume checkpoints (`--resume`)

Branch `resume_checkpoints`. Replaces ~15 ad-hoc lock/done files (`_lock`, `_collected`, `_processed`,
`.UMI_filtered.done*`, `.barcodes_done_<i>.tsv`, `.raw_barcodes_done`, `_tagged_bam_done`, ...) with one small
framework, `isoquant_lib/utils/checkpoints.py`. The full design history (three independent reviews, decisions
D1-D9) is in `.claude/RESUME_CHECKPOINTS_PLAN.md`.

## Model

A run is a **linear** sequence of stages; each finished stage leaves a checkpoint file, and `--resume` skips every
stage whose checkpoint exists. Plain checkpoints: no fingerprints (`--resume` forbids changing options anyway).

- `CheckpointStore(root, resume)`: a directory of checkpoints. Holds only a path + flag, so it is picklable and
  passed to `ProcessPoolExecutor` workers. Each checkpoint is its own JSON file (`{"finished_at", "state"}`,
  `checkpoint_path(name)`), written to `<name>.done.tmp.<pid>` and `os.replace`d, so per-chromosome workers mark
  concurrently and a killed process never leaves a half-written checkpoint. A missing or unreadable checkpoint
  (e.g. damaged by a power loss -- nothing is fsynced) means "not done": the stage just runs again.
- `run_stage(store, name, run, restore=None, cleanup=None, enabled=True, keep_tmp=False)`:
  disabled -> nothing; done -> log `"<name>: done in the previous run, skipping"`, `restore(state)` with the
  stage's saved `StageState` (`load_state`); otherwise `state = run()`, then the checkpoint. `cleanup` runs after
  the checkpoint **and again on skip** (removes fragments left by a crash between checkpoint and deletion), never
  with `--keep_tmp`.
- `mark_done_with_debug_hook(checkpoints, name, state)`: `CheckpointStore.mark_done` plus the debug hook;
  per-chr workers use it directly.
- `checkpoint_name(*parts)`: `/`-joined, every part sanitised like chromosome file names.
- Naming: a `CheckpointStore` argument is called `checkpoints` everywhere (`run_checkpoints` for the run store,
  `DatasetProcessor.sample_checkpoints(sample)` for an experiment's).
- No `invalidate`: a fresh run wipes the stores; on resume a missing checkpoint already means "not done", and
  per-chr child checkpoints (`collect/<chr>`) must survive when the parent stage reruns.

Stores:
- run store `<output>/checkpoints/` (`run_checkpoint_dir`): `cell_barcodes/<prefix>`, `barcodes/<prefix>`,
  `sample/<prefix>` (experiment complete), `combine_counts`, `fusion/<prefix>`;
- sample store `<sample>/aux/checkpoints/` (`SampleData.checkpoint_dir`): every stage below.

Fresh vs resume (`isoquant.py`):
- fresh run: `reset_checkpoints` in `check_and_load_args`, **before `save_params`** -- removes every sample store,
  resets the run store and writes the run record `<output>/checkpoints/run.json`
  (`{"checkpoint_format", "run_id", "isoquant_version", "started_at"}`); the run id is also saved in `.params`
  (`args.checkpoint_run_id`). So no old sample checkpoint survives a fresh run, and checkpoints can be tied to
  `.params`;
- `--resume` checks the record (`get_resume_incompatibility_reason`) and exits `RESUME_INCOMPATIBLE` (26) when it is missing or
  unreadable (a run of a version without checkpoints), has another `CHECKPOINT_FORMAT`, or its run id differs
  from `.params` (e.g. an older IsoQuant version reused the folder after a newer one: its `.params` has no run id,
  while the newer run's checkpoints are still there). Another IsoQuant version with the same format only warns;
- **bump `CHECKPOINT_FORMAT`** whenever stage names, saved stage states or the files a stage leaves for later stages
  change, and regenerate the CI resume data;
- a sample whose directory would be `<output>/checkpoints` is rejected (`_reject_reserved_sample_names`).

## Stages of one experiment (`DatasetProcessor.process_sample`)

- **Prologue (always runs)**: new `RunSummary`, reset `alignment_stat_counter`, `chr_ids`,
  `use_technical_replicas`, `sample.barcoded_reads = args.barcoded_reads`, tagged-BAM reference list.
1. `read_group_split`
2. `polya_split` (`--polya_trimmed list:/flnc:`, not `--read_assignments`)
3. `barcode_split` (PCR dedup, not `--barcoded_bam`)
4. `tagged_bam`
5. `collect` + per-chr `collect/<chr>`. A done chromosome is **reloaded**, not skipped (multimapper resolution
   needs all), and the reload calls `processed_reads_manager.finalize(chr_id)`: `ProcessedReadsManagerNoSecondary`
   (the default) computes the per-chromosome totals there, so without it a resumed run reported 0 assignments for
   every reloaded chromosome (and could flip `requires_polya_for_construction`). State: alignment stats +
   `covers_all_reads`, restored by `restore_collection_stats`. `_info` stays the source of totals /
   `all_read_groups`. "To keep these intermediate files for" is logged after the checkpoint (CI resume3 kills on it).
6. `umi/ED<d>` (+ per-chr), `umi_barcode2barcode/col<c>/ED<d>` (+ per-chr). Per-chr allinfo/stats under
   `out_umi_filtered_tmp` (aux), removed by the stage cleanup.
7. `dedup_bam`
- **Derive (always runs)**: `load_read_info`, barcode warning, polyA stats, polyA flags. D1: the flags are derived
  from `user_defined_require_mono{intronic,exonic}_polya` saved in `__init__`, not from the previous experiment's.
8. `process_assigned_reads` + per-chr `process_assigned_reads/<chr>` (`DatasetProcessor.process_assigned_reads`):
   assignment processing and counting in every mode, plus model construction when it is on. Workers close every
   output before the checkpoint.
- **After assignment processing (always runs)**: `report_assignment_stats` sums per-chr `_read_stat` / `_transcript_stat`
  (under the saves prefix), prints "Read assignment statistics" / "Transcript model statistics" (TestMode greps
  them), fills the run summary.
9. `merge` (parent checkpoint) with one stage per output, `merge/<output>` (built by `merge_printer_output`,
   `merge_gtf`, `merge_counter`, `merge_training_csv` in `DatasetProcessor.merge_outputs`):
   `read_info`, `assignments`, `corrected_bed`, `sqanti_t2t`, `polya_training`, `tss_training`, `gtf`,
   `extended_gtf`, `model_reads`, and **two stages per counter**: `merge/counter/<basename>` (truncate + merge counts
   and stats rows; cleanup keeps `.usable`) and `merge/counter_finalize/<basename>` (`load_usable_fragments` +
   `finalize`: TPM/matrix/MTX/loom; cleanup removes `.usable`). Merging is mechanical; conversion depends on
   options and group count.
   - **Counter constructors never write anything.** Each counter lists the files it appends to in
     `files_to_empty_before_counting()` (default `[output_file]`; `TerminalCounter` adds its training CSV), and the counting
     (per-chromosome) aggregator empties them all once at the end of its `__init__` (`empty_counter_outputs`).
     The merge aggregator is built with `for_merge=True`: no printers, nothing emptied. Each merge stage
     empties only its own final output, so a rerun cannot empty a finished one -- and a new counter is safe by
     default (forgetting an extra file in `files_to_empty_before_counting()` can only leave stale rows in a chromosome piece rerun
     after a crash, never empty a finished output).
   - Merge stages close outputs before their checkpoint (`close()` on printers, `with` in `merge_counts`).
   - `merge_files` / `merge_counts` take `remove_inputs` (merge stages pass False) and `fragment_dir`.
   - Merge stage names come from counter file names, which depend on the string pools, hence the parent checkpoint.
10. **Finish**: run summary, compress barcode tables; then `process_all_samples` marks `sample/<prefix>` and
    **after the checkpoint** deletes shared intermediates by glob (`remove_sample_intermediates`, also for skipped
    samples): split barcode / polyA tables, normalized polyA table, `aux/chr_fragments/` (unless `--keep_tmp`), and
    `out_raw_file_*`, `read_group_file*` (unless `--keep_tmp` or `--read_assignments`).

Per-chr merge fragments (one output's piece for one chromosome) live in `<sample>/aux/chr_fragments/`
(`SampleData.chr_fragment_dir`), not the output dir.
`--keep_tmp` keeps every temporary file **except** per-chr tagged/dedup BAM fragments (`merge_bam_files`
deletes them as before, D9).

## Adding a stage

- One function, wrapped in `run_stage` under a unique name (per sample: the sample store).
- In-memory state later code needs -> small JSON-able dict (`StageState`) returned by `run`, re-applied in `restore`. No paths
  (CI copies checkpoint dirs elsewhere), nothing large.
- Fragments only this stage consumes -> idempotent `cleanup`.
- Intermediates shared with later stages -> `remove_sample_intermediates` glob list.
- State needed whether or not stages are skipped -> prologue / derive / after-assignment-processing steps.
- Close every output before returning: the checkpoint is written right after.
- New counter: write nothing in `__init__`; if it appends to more files than `output_file`, list them in
  `files_to_empty_before_counting()`. New printer in `ReadAssignmentAggregator._init_printers`: None/Void when
  `for_merge=True` (printers open their file when built), plus a merge stage that opens, merges and
  closes only its own final file.

Per-chromosome completion lines ("Finished processing chromosome", "PCR duplicates filtered for chromosome") are
logged **after** the chromosome's checkpoint, so a run killed on such a line has that chromosome marked done; skipped
chromosomes are logged too ("... detected", "... during the previous run").

## Debugging and testing resume

`ISOQUANT_DEBUG_FAIL=<stage>[:before|:after]` (default `after`) raises `DebugFailure` at that checkpoint, in
`run_stage` and in per-chr workers (`collect/<chr>`, `process_assigned_reads/<chr>`, `umi/ED3/<chr>`). Parsed with the last
`:` only when the suffix is `before`/`after` (chromosome names may contain `:`). An env var, not a CLI flag, so it
is never pickled into `.params`. Typical check: run with the variable, `--resume` without it, diff against a clean
run ignoring `# Command line:` headers.

Unit tests: `isoquant_tests/test_checkpoints.py`, `isoquant_tests/test_resume_merge.py`.

CI resume data (`isoquant_tests/github/generate_resume_test.sh`) must be regenerated on the CI host: old
checkpoints exit with 26. Its stop strings still occur: "Finished processing chromosome" (collect and
process_assigned_reads workers), "Gene and transcript records have unequal strands: ..." (model construction),
"To keep these intermediate files for" (after collect); it checks the expected checkpoints after each partial run.

## Output changes in a clean run

- D1: multi-sample runs no longer leak `require_mono*_polya` from one experiment to the next.
- Developer training mode no longer leaves an empty `<p>.polyA_prediction.tsv.training.csv` in the output dir.
- Per-chr fragments never appear in the output dir; chromosome ids containing `/` are now made file-name safe
  everywhere (merges, UMI filtered reads / allinfo fragments, read-group splits).
- Log order changes (summary stats are printed after construction, merge logs per stage).
