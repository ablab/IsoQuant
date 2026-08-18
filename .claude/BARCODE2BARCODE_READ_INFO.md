# Per-column read_info for `--barcode2barcode` UMI rounds (deferred plan)

> **Status: DEFERRED — post-release enhancement.** Captured for a later
> iteration; too large for a pre-release change. Not yet implemented.

## Context / problem

In single-cell/spatial modes the default per-read output is the unified
`read_info.tsv.gz`; the legacy `allinfo` format is opt-in (`--large_output allinfo`).
When `--barcode2barcode` (spot-level UMI dedup) is set, IsoQuant runs an **extra
dedup round per spot column** (`dataset_processor.py:619-707`). Each extra round
does the full dedup compute but, by default, emits **only** a `.stats.tsv` — its
per-read result is exposed only as `allinfo`, which is off by default. So without
`--large_output allinfo` these rounds "run dry": compute is spent, no consumable
per-read output is produced.

Goal: make each extra round emit a per-column **read_info** file by default —
`SAMPLE_ID.UMI_filtered.barcode_barcode_col{C}.ED{N}.read_info.tsv[.gz]` — so the
spot-level deduplicated read set is available in the modern format.

## Why not the obvious shortcuts

- **Filter the main `read_info` by the column's survivor IDs** — rejected: dedup
  keeps a chromosome-global `selected_reads` set with an order-dependent
  "skip-if-already-selected" rule (`umi_filtering.py:367-382`), so spot-level
  survivors are **not** guaranteed ⊆ main survivors; filtering could silently drop
  reads.
- **Emit read_info inline during the dedup pass** — rejected as the primary path:
  the dedup loader `MergingSimpleReadAssignmentLoader` (`assignment_loader.py:93-124`)
  attaches **no `gene_info`**, but `ReadInfoPrinter.add_read_info` needs
  `gene_info.all_isoforms_introns` (`assignment_io.py:439`) and resolved
  `read_group`. Making it work means swapping in the genedb-backed loader and
  re-verifying group/cell_type resolution — more surface, weaker format guarantee.

## Recommended approach (reuse the exact main-read_info path, filtered per column)

Run the **same** loader+printer the main read_info uses, but with the column's own
survivor set as the filter. Output is then byte-for-byte the main read_info format,
restricted to that column's survivors (loaded from the raw saves, so correctness
does not depend on any subset property).

### Changes

1. **Parametrize the read filter path** — `assignment_loader.py`
   - `prepare_read_filter(chr_id, saves_prefix, use_filtered_reads, filtered_reads_path=None)`
     and `create_assignment_loader(..., filtered_reads_path=None)`: default to the
     current `filtered_reads_file_name(saves_prefix, chr_id)`; when given, use the
     supplied path. No behavior change for existing callers.

2. **Write a column-specific survivor-ID file during each extra round** —
   `dataset_processor.py:653-662` + `parallel_workers.py:344-392` +
   `umi_filtering.py:process_single_chr`
   - The extra-round `filter_umis_in_parallel` call currently passes
     `output_filtered_reads=False`. Give it a **column-specific** per-chr filtered-reads
     path derived from `bc2bc_prefix` (new helper in `file_naming.py`, e.g.
     `umi_filtered_reads_file_name(bc2bc_prefix, chr_id, edit_distance) + ".ids"`).
     `_process_chunk` already writes survivor read-IDs when `read_ids_outf` is set
     (`umi_filtering.py:379-380`); just route it to the column path.

3. **New read_info-only parallel worker** — `parallel_workers.py`
   - `write_read_info_in_parallel(sample, chr_id, chr_ids, args, saves_prefix,
     filtered_reads_path, output_read_info_path)`: mirror the read_info-producing
     subset of `construct_models_in_parallel` (`parallel_workers.py:204-223, 278-290`)
     — build string pools, `create_assignment_loader(..., use_filtered_reads=True,
     filtered_reads_path=<column file>)`, construct one `ReadInfoPrinter`
     (`additional_header=common_header`, `gzipped=args.gzipped`), iterate the loader
     and `add_read_info`. **No** model construction, counting, or GFF. Returns the
     per-chr file. Reuse the aggregator's `common_header` builder
     (`assignment_aggregator.py:60-79`) — factor it into a small helper if needed.

4. **Orchestrate + merge per column** — `dataset_processor.py:619-707`
   - After each column's dedup + stats, gate on `large_output_enabled(args, "read_info")`
     (default on), run `write_read_info_in_parallel` across `chr_ids` (same
     `ProcessPoolExecutor` pattern as `filter_umis`), and `merge_files(...)` the per-chr
     files into `output_prefix + ".read_info.tsv"` (+`.gz`), matching the main merge
     (`dataset_processor.py:738-739`, `copy_header=False, header_lines=3`). Add a lock
     file for resume and remove per-chr read_info + column filtered-reads temps after merge.

5. **Docs** — `docs/formats.md` (UMI section): document
   `SAMPLE_ID.UMI_filtered.barcode_barcode_col{C}.ED{N}.read_info.tsv[.gz]` as the
   default per-column output; keep allinfo as the opt-in legacy. Note it in
   `output.md` SC section and `.claude/TESTING_SYSTEM.md`/`BARCODE_CALLING.md` as
   appropriate.

### Minimal fallback (if the above is deemed too large)
Ungate the already-computed allinfo for the extra rounds (always merge it), giving a
per-column `.allinfo` by default. Smaller change, but keeps output in the legacy
format the user is moving away from — not recommended.

## Verification

- Unit: extend `isoquant_tests/test_umi_filtering.py` — assert the extra round writes
  a column survivor-ID file and that `write_read_info_in_parallel` produces a
  header + one line per survivor with the 18-column `read_info` layout.
- End-to-end (SC data on a self-hosted runner): run a `visium_hd` sample with
  `--barcode2barcode file.tsv:0:1,2` **without** `--large_output allinfo`; confirm
  `…UMI_filtered.barcode_barcode_col0.ED4.read_info.tsv.gz` (and col1) exist,
  are non-empty, share the main read_info header, and that line counts match the
  round's `.stats.tsv` survivor totals. Re-run with `--large_output allinfo` and
  confirm both read_info and allinfo appear.
- CI: add the new `read_info.tsv[.gz]` names to `check_input_files` in
  `SC.Mouse.VisiumHD.barcode2barcode.allinfo.yaml`.
- Regression: a run **without** `--barcode2barcode` is unchanged (new code path is
  gated on `args.barcode2barcode`); `create_assignment_loader` default path
  unchanged for all existing callers.
