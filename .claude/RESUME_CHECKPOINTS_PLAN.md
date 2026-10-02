# Stage/checkpoint framework for `--resume`

## Context

Resume checkpoints use about 15 different ad-hoc marker kinds (`_lock`, `_collected`, `_processed`, `.lock`, `.done*`,
`_tagged_bam_done`, `.barcodes_done_<i>.tsv`, `.raw_barcodes_done`). Each is created, checked and removed by hand in
`dataset_processor.py`, `parallel_workers.py`, `barcode_calling/pipeline.py` and `file_naming.py`, and temp-file
deletion is mixed in with the stages. Two independent reviews checked everything below against the code (file:line).

Defects (step 1 reproduces 1-3 before any change):
1. The split polyA/barcode tables and their locks are deleted right after `process_assigned_reads`
   (`dataset_processor.py:299-307`, `913-919`). Any crash after that point makes `--resume` split them again.
2. Resume destroys a completed sample. Sample 1's per-chr construct locks survive until the end-of-run `clean_up`
   (`:151-156`), so after a crash in sample 2 the merge runs again. Building the merge aggregator truncates every
   final output (`long_read_counter.py:195`, `terminal_counter.py:92,113`, `rna_velocity_counter.py:119`,
   `assignment_io.py:75`, `transcript_printer.py:38`). `merge_files` skips the per-chr files that are already deleted
   (`file_utils.py:153,166-168`). `merge_counts` then opens a deleted per-chr `.usable` file
   (`file_utils.py:177-178` → `long_read_counter.py:459`) and raises FileNotFoundError. Counters without a usable file
   (exon/intron) are left empty without any error. A crash in the middle of a merge has the same effect.
3. A skipped `collect_reads` (`:401-404`) leaves `alignment_stat_counter` empty. The `__not_aligned` row then reports 0
   instead of the unmapped count (`file_utils.py:190`; the per-chr stats always carry `__not_aligned 0`,
   `long_read_counter.py:434`), and the run summary has no alignment stats. The barcode skip restores sample state by
   hand (`pipeline.py:129-134`).
4. `filter_umis` `return`s when it finds the per-ED marker (`:657-660`). That skips the barcode2barcode rounds and the
   global lock (`:806`), and the per-ED and barcode2barcode locks are never removed. `clean_up` then runs `os.remove`
   on the missing global lock (`:158`) and raises FileNotFoundError.
5. The deduplicated BAM, the merges, `combine_counts` and fusion have no checkpoint.
6. An existing cross-sample leak: `process_sample` ORs `args.require_mono{intronic,exonic}_polya` with their previous
   values (`:288-295`). Under the auto strategy, a flag set by sample N-1 therefore stays set for sample N.

Agreed decisions:
- **Scope:** all lock-file stages, from barcode calling through UMI filtering, plus the stages that have no
  checkpoint today.
- **Model:** linear; a stage is skipped if its marker exists.
- **Shared intermediates:** deleted after the sample marker; completed samples are skipped entirely on resume.
- **Stage-private fragments:** deleted only after that stage's own marker.
- **Markers:** plain markers with an optional small JSON payload, no fingerprints.
- **Implementation:** in-house, no new dependencies.
- **Location:** `<sample>/aux/checkpoints/` and `<output>/checkpoints/`.
- **Compatibility:** none with old checkpoints; the CI resume data is regenerated.

Out of scope:
- The mtime caches and the gunzipped reference. Its cleanup in `clean_up` is dead code, because
  `args.gunzipped_reference` is always None (`dataset_processor.py:117`).
- `--read_assignments` appears to be broken already: `get_chromosome_ids` opens the `.save` path as a BAM (`:346-347`),
  and every sample uses `read_assignments[0]` (`:249`). Its code paths must keep compiling, but it isn't tested here (D4).
- `--clean_start` re-maps every sample on resume. That is the mtime-cache domain, so it is unchanged.

## Design

### `isoquant_lib/utils/checkpoints.py`

```python
from __future__ import annotations
from typing import Any, Callable, Dict, Optional

Payload = Optional[Dict[str, Any]]      # runtime alias: no builtin generics or `|` at module level (3.8)

class CheckpointStore:
    """A marker directory. Holds only (root, resume); picklable for ProcessPoolExecutor."""
    def __init__(self, root: str, resume: bool) -> None
    def reset(self) -> None                     # rmtree + makedirs; driver only, fresh runs only
    def ensure_dir(self, name: str) -> None     # optimisation only; mark_done creates parents itself
    def is_done(self, name: str) -> bool        # resume and the marker file exists
    def load(self, name: str) -> Payload
    def mark_done(self, name: str, payload: Payload = None) -> None
        # makedirs(parent, exist_ok=True); JSON {"isoquant_version","finished_at","payload"}
        # -> "<name>.done.tmp.<pid>" -> os.replace
def marker_name(*parts: str) -> str             # "/"-join; EVERY part sanitised by convert_chr_id_to_file_name_str

def run_stage(store: CheckpointStore, name: str, run: Callable[[], Payload],
              restore: Callable[[Payload], None] | None = None,      # annotations are not evaluated: fine on 3.8
              cleanup: Callable[[], None] | None = None,
              enabled: bool = True, keep_tmp: bool = False) -> None
    # not enabled -> return
    # is_done -> log "<name>: done in the previous run, skipping"; restore(payload); cleanup() unless keep_tmp
    # else payload = run(); [debug hook before]; mark_done(payload); [debug hook after]; cleanup() unless keep_tmp
```
- **No `invalidate`.** A fresh run calls `reset()`. On resume a missing marker already means "not done", and per-chr
  child markers (`collect/<chr>`) must survive when the parent stage runs again.
- **`cleanup` is idempotent.** It also runs on skip, which removes fragments left by a crash between the marker and
  the deletion.
- **Payloads are small JSON only:** counts and enum-name→int dicts.
  - Never paths: CI copies checkpoint dirs to other locations.
  - Never `all_read_groups`.
  - Helpers reuse `report.run_summary.enum_stats_to_dict`, and the reverse maps names back via `EnumClass[name]`.
- **Debug hook:** the env var `ISOQUANT_DEBUG_FAIL=<stage>[:before|:after]` (default `after`).
  - It raises in `run_stage` and in the per-chr workers.
  - It is parsed with `rsplit(":", 1)`, and the suffix counts only if it is `before`/`after`, because chromosome names
    can contain `:`.
  - It is not a CLI flag, because a flag would be pickled into `.params`.
- **Stores are built from paths** (`sample.checkpoint_dir`, `args.output`) and never stored on `SampleData`, which
  `map_reads` rebuilds.

### Fresh vs resume
- **Fresh run:** in `check_and_load_args`, right before `save_params`:
  - first `shutil.rmtree(sample.checkpoint_dir, ignore_errors=True)` for every sample (no makedirs; the sample dirs may
    not exist yet), then `reset()` the run store;
  - a `.params` file written by this version then implies that `<output>/checkpoints/` exists, and that no sample marker
    from an earlier run in the same directory survives (no window where new `.params` sit next to old sample markers);
  - `create_output_dirs` creates `aux/checkpoints/` and `aux/per_chr/` with `exist_ok=True`.
- **Resume:** if `<output>/checkpoints/` is missing, the run predates checkpoints. Exit with
  `IsoQuantExitCode.RESUME_INCOMPATIBLE = 26` and ask to restart without `--resume`. No marker globs.
- `check_input_params` rejects any sample with
  `os.path.normpath(sample.out_dir) == os.path.normpath(os.path.join(args.output, "checkpoints"))` (INVALID_PARAMETER).
  That sample directory would be the run store, and `reset()` would delete it. Comparing paths also catches YAML names
  like `checkpoints/`; CLI `--prefix`, YAML `name` and the `<prefix><i>` defaults all end up in `sample.out_dir`.

### Stage list

Run store `<output>/checkpoints/`:
- `cell_barcodes/<prefix>` (pass 1) and `barcodes/<prefix>`.
  - `restore` recomputes `barcoded_reads` (and `file_list` with `--split_molecules`) from `SampleData` paths; no
    payload.
  - The cleanup may also remove `aux/barcode_calling_*` dirs left by a crash (`detect_barcodes.py:578-632`).
- (unchanged: gunzip, GTF→DB, mapping)
- `sample/<prefix>`: the sample is complete, so `process_sample` is skipped. Only the glob-based deletion of step 10
  runs again; it needs no prologue state.
- `combine_counts`; `fusion/<prefix>`, marked only when that sample succeeds.

Sample store `<sample>/aux/checkpoints/`. `process_sample` runs these in order:
- **Prologue (always runs):**
  - new `RunSummary`; reset `alignment_stat_counter`;
  - `chr_ids`, `use_technical_replicas`;
  - `sample.barcoded_reads = args.barcoded_reads` if set;
  - the tagged-BAM reference list.
1. `read_groups_split`
2. `polya_split`; enabled for `polya_trimmed.uses_read_table()` and not `read_assignments`
3. `barcode_split`; enabled for PCR dedup and not `barcoded_bam`
4. `tagged_bam`; BAM fragments handled as today
5. `collect` + per-chr `collect/<chr>`; enabled unless `read_assignments`
   - The per-chr skip keeps the existing reload path and its fallback, because multimapper resolution needs every
     chromosome.
   - Payload: alignment stats + `covers_all_reads`. `restore` refills `alignment_stat_counter` and the run summary.
   - `_info` is still written and stays the only source of the totals and `all_read_groups`.
   - "To keep these intermediate files for" is logged **after** the `collect` marker (CI resume3 kills on it).
6. `umi/ED<d>` (+ per-chr), then `umi_bc2bc/col<c>/ED<d>` (+ per-chr); enabled for PCR dedup unless `read_assignments`
   - The per-chr files use the prefix `out_umi_filtered_tmp` (it exists and is unused today).
   - `out_umi_filtered_done` is removed; `umi_barcode2barcode_prefix` is kept and applied to the tmp prefix.
7. `dedup_bam`; enabled for the same conditions as `umi`, plus `--large_output deduplicated_bam`
- **Derive (always runs):**
  - `load_read_info`;
  - the barcode warning;
  - the polyA log + `set_polya_stats` ('polyA tail detected in N');
  - the polyA flags (D1).
8. `construct` + per-chr `construct/<chr>`
   - The worker **closes every output it opened before `mark_done`**: the aggregator printers (new
     `ReadAssignmentCompositePrinter.close()`), `sqanti_t2t_printer`, `model_reads_printer` and both GFF printers.
     Today the per-chr read_info / read_assignments / corrected_bed / t2t fragments are flushed only by `__del__`,
     after the lock (`parallel_workers.py:369`), so a kill in between leaves a truncated fragment marked done.
- **Post-construct (always runs):**
  - sum the per-chr `_read_stat` files, and with model construction the `_transcript_stat` files too, under the
    `saves_file` prefix (`sample.out_raw_file`; `read_assignments[0]` under D4);
  - print "Read assignment statistics" (genedb only) and "Transcript model statistics", and call
    `set_assignment_stats` / `set_transcript_model_stats`;
  - TestMode needs 'unique: 1' and 'known: 2' from these prints;
  - the stat files live under `out_raw_file` until the sample finishes, so no payload is needed.
9. `merge/<unit>` (D2)
   - **Lazy setup via a parent marker `merge`:**
     - `merge` is written after the last unit and the catch-all cleanup. If it exists, nothing in step 9 runs.
     - Otherwise the merge aggregator and its string pools are built **once**, before the units, with
       `truncate_outputs=False`.
     - They can't be built per unit, because counter unit names come from `output_counts_file_name`. That depends on
       the string pools (`ignore_read_groups = string_pools is None`, `long_read_counter.py:311,581,857` →
       `.linear.tsv` vs `.tsv`, `:192`).
     - Otherwise they are built exactly as today: `chr_id=None`, `load_barcode_pool=False`, `load_tsv_pools=False`.
   - **Text units:** `read_info`, `assignments`, `corrected_bed`, `sqanti_t2t`, `gtf`, `extended_gtf`, `model_reads`,
     `polya_training`, `tss_training`. Units for outputs that aren't enabled are disabled.
   - **Two units per counter**, for every counter in `global_counter`, `transcript_model_global_counter` and
     `gene_model_global_counter` (D7). Merging is a mechanical concatenation of fragments. Converting is a separate
     concern that depends on `--counts_format`, `--normalization_method` and the number of groups. As a side effect,
     peak disk stays as today:
     - `counter/<counts basename>`: truncate `counter.output_file`, merge the count fragments and the stats rows, and
       `close()`. Its cleanup deletes the count and `.stats` fragments and keeps `.usable`. `merge_counts` no longer
       loads `.usable`: that would load it twice in a clean run (`load_usable` adds with `+=`,
       `long_read_counter.py:461,469`), and so change TPM under `--normalization_method usable_reads`.
     - `counter_finalize/<counts basename>`: the new `load_usable_fragments(counter, label, chr_ids, fragment_dir)`,
       then `finalize`, which writes the TPM/matrix/MTX/loom files. Its cleanup deletes `.usable`. This works in a
       fresh process:
       - the merge-time pools are built only from static inputs, and `load_usable` never adds to a pool;
       - `TerminalCounter` and `ExonSpliceSiteCounter` finalize are no-ops;
       - `RNAVelocityCounter.finalize` reads only the merged TSV.
   - **Unaligned count:** the restored `AlignmentType.unaligned` for the global and transcript-model counters, and 0
     for the gene-model counters (as today).
   - **Finalize gating:** the global-counter finalize only with genedb; the model finalizes only with model
     construction.
   - **No unit truncates another unit's output:**
     - printers are built inside the unit's run;
     - `truncate_output=False` is threaded through every counter constructor (see File changes), including the three
       classes that bypass `AbstractCounter.__init__`: `TerminalCounter` (both the TSV and the training CSV),
       `PolyACounter`/`TSSCounter`, and `RNAVelocityCounter`;
     - a counter unit truncates its own file itself.
   - **Units close their outputs before returning.** There is a new `close()` on `TextFileAssignmentPrinter` (and so
     SqantiTSVPrinter) and on `GFFPrinter`, and `merge_counts` uses `with`. A marker therefore never precedes an
     unflushed gzip stream.
   - **Units can be rerun:** fragments are deleted only after the unit's marker. `merge_files`/`merge_counts` get
     `remove_inputs: bool = True`, and the units pass False.
   - **Every finalize is idempotent:** the converters open with "w", and the velocity finalize removes a stale loom
     (`rna_velocity_counter.py:301`).
   - After all units: `_cleanup_per_chr_temp_files` as a catch-all, now globbing the fragment dir.
   - **Fragment location (D8):** the per-chr merge fragments move from `sample.out_dir` to a new
     `sample.chr_fragment_dir = <sample>/aux/per_chr/`.
     - All `SampleData.get_*_file(chr_id)` getters built from `get_chr_prefix` (`input_data_storage.py:183-253`)
       use this dir.
     - So does the per-chr `GFFPrinter` in `construct_models_in_parallel` (`parallel_workers.py:279-300`, today
       `sample.out_dir`).
     - `merge_file_list` / `merge_files` / `merge_counts` get `fragment_dir: Optional[str] = None`. When it is set, the
       per-chr names are `os.path.join(fragment_dir, f"{label}_{sanitised chr}{rest}")` and not the final file's dir.
       - Every merge caller passes `fragment_dir=sample.chr_fragment_dir`: read_info, assignments, corrected_bed,
         sqanti_t2t, gtf, extended_gtf, model_reads, the three counter loops, and the two training-CSV merges
         (`dataset_processor.py:967,973`). The training-CSV *writer* follows automatically
         (`output_prefix + TRAINING_SUFFIX`), but its merge doesn't.
       - Chromosome ids are sanitised with `convert_chr_id_to_file_name_str`, like `get_chr_prefix`. This fixes an
         existing mismatch for chr ids containing `/`, since `merge_file_list` used the raw id (`file_utils.py:143`).
       - With `fragment_dir` set and `label` not a basename prefix, it raises `ValueError` instead of using the
         `rreplace` fallback.
     - Everything else that writes per-chr files derives its path from these getters or the two GFF printers. The
       readers of `out_dir` (run summary, visualizer, `combine_counts`) only ever see merged files.
     - With `--keep_tmp` (D3) the fragments therefore stay under `aux/`, not next to the outputs. The visualizer's
       per-chr filter (`post_process.py:315-320`) becomes redundant but stays.
10. **Finish (not a stage):**
    - `compress_barcode_tables`, `write_run_summary`, mark `sample/<prefix>`;
    - then delete by glob:
      - `barcodes_split_reads + "_*"`, `polya_split_reads + "_*"`, `polya_reads_normalized`: unless `keep_tmp`;
      - `out_raw_file + "_*"`, `read_group_file + "*"`: unless `keep_tmp` or `read_assignments`.

`DatasetProcessor.clean_up` and `clean_locks` are removed. The BAM stages need no per-fragment handling, since a
rerun rewrites all fragments.

### File changes
- **New:** `isoquant_lib/utils/checkpoints.py`, `isoquant_tests/test_checkpoints.py`.
- **`dataset_processor.py`:**
  - restructure `process_sample` into prologue / stages / derive / construct / post-construct / merge / finish;
  - split `process_assigned_reads` into `construct_models` + `merge_outputs`;
  - make `filter_umis` one `run_stage` per ED/column;
  - remove `clean_up` and `clean_up_external_polya`;
  - D1: `__init__` saves the original `require_monointronic_polya` / `require_monoexonic_polya`.
- **`parallel_workers.py`:**
  - the three per-chr workers take a `CheckpointStore` + stage name instead of lock paths;
  - each closes every output it opened before `mark_done`. Collect and UMI already close theirs (`:208`,
    `umi_filtering.py`); construct needs the fix described in stage 8;
  - the per-chr GFF printers write to `sample.chr_fragment_dir`.
- **`barcode_calling/pipeline.py`:** `call_barcodes` and `detect_cell_barcodes` take the run store.
- **`isoquant.py`:**
  - create the stores; reset the run store before `save_params`;
  - exit 26 on resume without `checkpoints/`;
  - reject a sample named `checkpoints`;
  - add the `combine_counts`/fusion markers.
- **`utils/input_data_storage.py`:**
  - in `_init_paths` (so the `SampleData` that `read_mapper.map_reads` rebuilds has them too), add
    `checkpoint_dir = <aux>/checkpoints` and `chr_fragment_dir = <aux>/per_chr`;
  - all `get_chr_prefix`-based getters (`:181-253`) use `chr_fragment_dir`;
  - remove `barcodes_done`, `raw_barcodes_done`, `out_umi_filtered_done` and the `get_*_lock_file` methods.
- **`utils/file_naming.py`:** remove the lock/done helpers and `clean_locks`; keep `umi_barcode2barcode_prefix`.
- **`utils/file_utils.py`:**
  - `merge_file_list(fname, label, chr_ids, fragment_dir=None)`: sanitised chr ids, fragment dir, and `ValueError`
    instead of the fallback when `fragment_dir` is set;
  - `merge_files` gains `remove_inputs: bool = True` and `fragment_dir`;
  - `merge_counts(counter, label, chr_ids, unaligned_reads=0, remove_inputs=True, fragment_dir=None)` merges counts +
    stats rows only, through `with`;
  - its `.usable` loop moves to the new `load_usable_fragments(counter, label, chr_ids, fragment_dir)`.
- **`truncate_output: bool = True`** goes on the following, so the merge driver can pass False:
  - **`quantification/long_read_counter.py`:**
    - `AbstractCounter`, `AssignedFeatureCounter`;
    - `ProfileFeatureCounter` + ExonCounter/ExonUsageCounter/IntronCounter/IntronRetentionCounter;
    - `ExonSpliceSiteCounter`;
    - `create_gene_counter`/`create_transcript_counter`.
  - **`terminal_prediction/terminal_counter.py`:** `TerminalCounter` (gates `:92` and `:113`), `PolyACounter`,
    `TSSCounter`.
  - **`quantification/rna_velocity_counter.py`:** gates `:119`.
- **`assignment/assignment_aggregator.py`:** `__init__(…, truncate_outputs: bool = True)`.
  - It stores `self.truncate_outputs` and passes it to every counter it builds: all `_add_*` methods **and
    `_init_grouped_model_counters`** (`:318-341`), which builds the grouped discovered_transcript/discovered_gene
    counters directly.
  - With False, the four printers are None (`t2t_sqanti_printer` is `VoidTranscriptPrinter()`), `global_printer` is
    an empty composite, and construction has no filesystem side effects.
- **`assignment/assignment_io.py`, `model_construction/transcript_printer.py`:** add `close()` on
  `TextFileAssignmentPrinter`, `ReadAssignmentCompositePrinter`, `GFFPrinter` and the Void printers.
- **`utils/bam_utils.py`:** unchanged. The per-chr tagged/dedup BAM fragments are deleted even with `--keep_tmp`, as
  today (D9).
- **`utils/error_codes.py` + `isoquant_tests/github/error_codes.py`:** add `RESUME_INCOMPATIBLE = 26`.
- **Tests:**
  - drop the `umi_barcode2barcode_global_lock` test in `test_barcode2barcode.py` (lines 21, 69; keep the prefix
    usage);
  - add `remove_inputs=False` and `fragment_dir` cases to `test_file_utils.py`, including a chr id containing `/`;
  - an aggregator test: building with `truncate_outputs=False` over pre-filled outputs leaves them byte-identical
    (includes the grouped model counters and the terminal/velocity counters);
  - the checkpoint tests also cover the `Payload` import on 3.8, debug-var parsing with `:` in names, and the
    `truncate_output=False` constructors.
- **Docs:**
  - new `.claude/RESUME_CHECKPOINTS.md`: stage list, prologue/derive/post-construct rules, how to add a stage,
    cleanup rules, the debug env var; this plan is also saved there;
  - update `CLAUDE.md`, `docs/cmd.md` (`--resume`, exit 26), `SC_IO_OUTPUTS.md:30-45,129-132` + §Resume,
    `CELL_BARCODE_SELECTION.md`, `BARCODE_INTEGRATION.md`, `BARCODE_CALLING.md:308`, `BARCODE2BARCODE_READ_INFO.md:77`,
    `INTRON_GRAPH_FLOW_DUMP.md:346`, `TESTING_SYSTEM.md:159-161`, `POLYA_TSS_TRAINING.md:75,218` and
    `POLYA_TSS_DETECTION.md:277` (both cite `merge_assignments`);
  - `docs/output.md:243`: describe `aux/per_chr/`, `aux/checkpoints/` and `<output>/checkpoints/`.
- **Clean-run output changes besides D1, for the PR:**
  - developer training mode no longer leaves an empty `<p>.polyA_prediction.tsv.training.csv` in the output dir;
  - per-chr fragments are no longer in the output dir;
  - chr ids containing `/` now merge correctly;
  - the log order changes.
- **Nits:** restore the alignment stats through `EnumStats.add`, and use `glob.escape(prefix)` in the new globs.
- **Conventions:** type hints throughout; `requires-python` stays `>=3.8`.
  - New modules start with `from __future__ import annotations` and may write `X | None` / `list[str]` **in
    annotations only**. Module-level aliases and runtime code use `typing` (`Optional`, `Dict`).
  - Edits to existing modules follow each file's current style; `str.removeprefix` is not used.

**Rule for adding a stage:**
- Write one function and wrap it in `run_stage` under a unique name.
- Return in-memory state that later code needs as a small JSON payload, and re-apply it in `restore`. No paths.
- Remove fragments that only this stage consumes in an idempotent `cleanup`.
- Put intermediates shared with later stages on the sample finish glob list.
- State needed whether or not stages are skipped goes into the prologue / derive / post-construct steps.
- Close every output before returning (the marker follows immediately).

## Implementation order
1. Reproduce defects 1-3 on master with `run_until.py`. Use a 2-sample YAML written under `/abga/work/andreyp`, not the
   repo or `/tmp` (`chr9.4M.yaml` has one sample) and kill it on "Processing experiment <sample2>". Record the failures.
2. `checkpoints.py` + unit tests.
3. The per-chr workers + the `collect`/`construct`/`umi` stages.
4. The `truncate_output` threading and printer `close()`; moving fragments to `aux/per_chr/` (`fragment_dir`); the
   per-output `merge/<unit>` stages; the
   prologue/derive/post-construct split; the D1 fix.
5. The split/tagged/dedup stages, the barcode stages, the finish step, the completed-sample skip, the run-level
   markers, exit 26 and the `checkpoints` name check.
6. Docs. The CI resume data is then regenerated on the CI host with `generate_resume_test.sh`. Check that its 4 stop
   strings (across 6 runs) still occur:
   - 'Finished processing chromosome SIRV1/SIRV5' (from both the collect and construct workers);
   - 'Gene and transcript records have unequal strands: SIRV5: -, SIRV504: +';
   - 'To keep these intermediate files for'.

## Verification (only when asked, per CLAUDE.md)
- `pytest isoquant_tests/test_checkpoints.py`:
  - atomic write; payload round-trip; reset;
  - skip + restore; cleanup after the marker and again on skip;
  - disabled stages; picklable store; sanitised names; debug-hook parsing.
- For each stage S, on these runs:
  - `simple_data` bulk BAM;
  - a 2-sample YAML;
  - a 10x run with UMI + `--large_output tagged_bam deduplicated_bam` + `--polya_trimmed list:`. The repo has no 10x
    UMI toy data, so the source needs to be named; the CI SC data is the candidate;
  - a barcode2barcode run.

  Each check:
  - run with `ISOQUANT_DEBUG_FAIL=S`, then `--resume`; include `construct/<chr>:after` with `-t 1`;
  - also `S:before` for `merge/<unit>` and `collect/<chr>`, and a crash after the first counter unit (finished units
    must not be truncated);
  - the outputs must equal a clean run, ignoring the `# Command line:` header lines and the summary `command_line`,
    which resume rewrites;
  - the log must show no rerun of anything already marked done.
- Specific checks:
  - crash in sample 2 → sample 1's files are byte-identical;
  - crash after construct → no "Splitting read barcode table" in the log;
  - the resumed `__not_aligned` value and the summary alignment stats equal the clean run;
  - resuming a completed run does nothing (except re-mapping under `--clean_start`, out of scope);
  - multi-sample resume gives the same polyA flags as a clean run;
  - after a clean run, no `<prefix>_<chr>.*` files are left in the output dir; with `--keep_tmp` they are all under
    `aux/per_chr/`;
  - `--test` log contains 'total assignments 4', 'polyA tail detected in 2', 'unique: 1', 'known: 2',
    'Processed 1 experiment'.
- `./isoquant.py --test`, full `pytest`, `console_test.py`, and the CI `Resume_tests` with the regenerated data.

## Decisions (answered)
- **D1:** fix the polyA-flag leak. Save the original `require_monointronic_polya` / `require_monoexonic_polya`
  (set by `set_model_construction_options`, `isoquant.py:1191-1192`) in `DatasetProcessor.__init__` and derive each
  sample's flags from them. `requires_polya_for_construction` needs no saved original, because it is recomputed per
  sample. This can change clean multi-sample results for sample 2 onwards; it will be noted in the changelog and PR.
- **D2:** per-output merge markers, so peak disk stays as today.
- **D3:** `--keep_tmp` keeps every temporary file, including split tables and merge fragments (except the BAM
  fragments, D9).
- **D4:** `--read_assignments` is fixed in a separate task; here its code paths only have to keep compiling.
- **D5:** the CI resume data is regenerated on the CI host after merging.
- **D6:** `--resume` with no `<output>/checkpoints/` exits with 26, including for runs that completed under an old
  version. No marker globs.
- **D7:** two merge units per counter (`counter/<b>`, `counter_finalize/<b>`). This keeps merging separate from the
  parameter- and group-count-dependent conversion; that peak disk stays as today is secondary.
- **D8:** the per-chr merge fragments live in `<sample>/aux/per_chr/` instead of the output dir.
- **D9:** the per-chr tagged/dedup BAM fragments are an explicit exception to D3. `merge_bam_files` deletes them even
  with `--keep_tmp`, as today: they are as large as the final BAM and rarely useful for debugging.
