# RNA Velocity Counts

Per-`(cell barcode, gene)` **spliced / unspliced** read tallies plus a
velocyto-style `.loom` export, for downstream RNA velocity analysis (scVelo).

Origin: branch `rna_velocity` (`isoquant_lib/rna_velocity_counter.py`,
commits `5278c85c` … `c58da75f`, by @arina-c and @andrewprzh), rebased onto
master as `isoquant_lib/quantification/rna_velocity_counter.py`. The branch
also carried the polyA/TSS prediction work, which reached master separately —
see `.claude/POLYA_TSS_DETECTION.md`. Only the velocity part was replayed.

## What it produces

Two files per cell-barcode grouping strategy:

| File | Contents |
|---|---|
| `SAMPLE.RNA_velocity_grouped_<strategy>` | TSV: `#cell_id`, `gene_id`, `spliced`, `unspliced` |
| `SAMPLE.RNA_velocity_grouped_<strategy>.loom` | Same counts; genes × cells, layers `''` (= spliced), `spliced`, `unspliced` |

The unnamed main layer repeats the spliced matrix, matching velocyto's
convention, so scVelo finds the layers it expects by name.

## When it runs

No CLI flag. `_add_grouped_velocity_counter`
(`isoquant_lib/assignment/assignment_aggregator.py`) adds the counter when
**all** of the following hold:

1. `args.read_group` and `args.genedb` are set (the enclosing
   `_init_grouped_counters` gate), **and**
2. `args.mode != IsoQuantMode.bulk`, **and**
3. the grouping strategy's pool type starts with `barcode` — i.e. `barcode`,
   `barcode_spot`, or `barcode_barcode`.

Condition 3 is what keeps velocity off `file_name` / `tag:` / `file:` groupings,
where a "cell" would be a whole sample. It reads `get_grouping_pool_types(args)`
(`isoquant_lib/assignment/read_groups.py`), computed **once** in
`_init_grouped_counters` rather than per strategy.

> **Index alignment.** `get_grouping_pool_types` is keyed by *grouper* index,
> not spec index (its docstring says "spec_index" and is stale — it increments a
> separate `grouper_index` that expands for multi-column file specs and
> multi-column `barcode2spot`). `get_grouping_strategy_names` expands in exactly
> the same order, so `grouping_pool_types[group_idx]` lines up with
> `grouping_strategy_names[group_idx]`. Any change to the expansion rules must
> keep the two functions in lockstep.

Each strategy writes to its own path. Sharing one path across strategies makes
the merge stage delete the same per-chr fragment twice (`FileNotFoundError`).

## Architecture

`RNAVelocityCounter` is an ordinary `AbstractCounter` living in
`ReadAssignmentAggregator.global_counter`, so it rides the existing counter
machinery with no special-casing:

```
worker (per chromosome)                   main process (per sample)
──────────────────────                    ─────────────────────────
aggregator.global_counter
  .add_read_info(ra)   → tally in dicts
  .dump()              → per-chr TSV fragment
                                          merge_counts(counter, ...)
                                            → concatenates fragments,
                                              header_lines=1
                                          global_counter.finalize(args)
                                            → create_loom()
```

It bypasses `AbstractCounter.__init__` (same trick as `TerminalCounter`): the
output path is already a full name, so the `counts_file_name` suffix machinery
would only get in the way. `output_counts_file_name` is set to the same path
because `merge_counts` reads that attribute; `output_stats_file_name` and
`usable_file_name` are `None` so the stats/usable merge branches are skipped.

`add_confirmed_features` / `add_unassigned` / `add_unaligned` are explicit
no-ops. `CompositeCounter` fans those out to every counter, and the base class
raises `NotImplementedError`; nothing dispatches them on `global_counter` today,
but inheriting the raise is a trap for whoever adds the first call.

### Header handling

Every per-chr fragment carries a one-line `#cell_id ...` header, because
`merge_counts` passes `header_lines=1` — it keeps the first fragment's header
and strips exactly one line from each subsequent fragment. A headerless
fragment would lose its first *data* row in the merge.

### Interned keys

Tallies are keyed on `(int group id, int gene id)`, resolved to strings once in
`dump()`. A chromosome's worth of `(cell, gene)` pairs runs to millions of keys
in a real single-cell run; two ints per key instead of two strings is the
difference between a workable and an unworkable worker footprint, and it keeps
`pool.get_str()` off the per-read path. This matches how the feature counters in
`long_read_counter.py` store their group ids.

`dump()` clears both dicts. It is called once per chromosome today, but a second
call must not re-emit the same counts.

## Classification — and its limits

```python
SPLICED   = {unique, unique_minor_difference, ambiguous}
UNSPLICED = {inconsistent, inconsistent_ambiguous}
```

Everything else is dropped, as are reads with no barcode and reads with no gene.

**This uses assignment consistency as a proxy for splicing status, which is not
what RNA velocity actually needs.** Three consequences, in rough order of impact:

1. **Intronic reads — the strongest unspliced signal — are dropped entirely.**
   A read lying inside an intron goes through `assign_to_overlapping_genes`,
   which leaves `assignment_type == noninformative` and records the gene only in
   `gene_assignment_type` (`inconsistent_genic` / `inconsistent_multigenic`).
   Neither type is in `ACCEPTED_ASSIGNMENT_TYPES`, and `noninformative` is not
   either, so pre-mRNA reads never reach the counter. See
   `.claude/READ_ASSIGNMENT_LIFECYCLE.md` § "Genic reads without a transcript".
2. **`inconsistent` over-collects.** A read with a novel exon, a shifted splice
   site, or an alternative TSS is fully spliced but lands in `inconsistent` and
   is counted as unspliced.
3. **`inconsistent_non_intronic` is excluded**, which is right — it explicitly
   marks non-intronic disagreement — but it is excluded by omission rather than
   by intent, so it is easy to "fix" wrongly.

The principled signal is the intron-retention match events already computed by
the assigner (`MatchEventSubtype.intron_retention`,
`unspliced_intron_retention`, `incomplete_intron_retention_{left,right}`; note
`fake_micro_intron_retention` is an artifact class and must be excluded), plus
the `genic_intron` reads from tier 1. Moving to those would mean reading
`gene_assignment_type` alongside `assignment_type` and inspecting
`isoform_matches[*].match_subclassifications`. Not done — it changes the numbers
and wants validation against velocyto on a real dataset.

## Gene attribution

`_get_gene_id` returns the first `isoform_matches[*].assigned_gene_id` that is
not `None`, mirroring the feature counters. Taking `isoform_matches[0]`
unconditionally is wrong: `match_inconsistent`'s `quick_mode` exit
(`long_read_assigner.py:618`) builds an `IsoformMatch(MatchClassification.genic)`
with **no gene**, which produced rows with an empty `gene_id` column.

For ambiguous reads (`ambiguous`, `inconsistent_ambiguous`) the first named gene
wins outright — no 1/N split across candidate genes, unlike the gene counter.

## Barcode attribution

`_get_group_id` drops a read when:

- `ignore_read_groups` (no string pools), or
- `read_group_ids` is empty — no barcode was detected. The feature counters fold
  these into group `0`; velocity output is per cell, so inventing a cell is worse
  than dropping the read. Indexing unconditionally here raised `IndexError`
  inside the worker.
- `read_group_ids[group_index] < 0`. `read_group_to_ids`
  (`isoquant_lib/utils/string_pools.py`) stores `-1` when a strategy produced no
  value for the read. `StringPool.get_str` indexes a plain list, so `-1` resolves
  to *the last barcode in the pool* — a silent mis-attribution to a real cell.
  **This sentinel is unhandled elsewhere in the codebase too** (e.g.
  `long_read_counter.py:393`, `:595`, `:761`, `:986` all call `get_str(group_id)`
  on a possibly-`-1` id); that is pre-existing and out of scope here.

## Known limitations

- **No CLI gate.** Any non-bulk run with a barcode grouping and a genedb
  produces these files. There is no way to opt out, and `loompy` is a hard
  dependency in `requirements.txt` / `pyproject.toml` purely for this export.
  The natural fit is a new `--analysis rna_velocity` value (see
  `.claude/ANALYSIS_OPTION.md`) resolving to an internal `args.rna_velocity`
  flag. Not implemented — it changes default behaviour for every single-cell
  user.
- **`finalize()` loads the whole merged TSV into a pandas DataFrame**, then
  builds two matrices. Loom is a *dense* HDF5 format: 100k cells × 30k genes is
  ~12 GB per layer before compression. Fine at current scales, a wall at
  billion-read ones. The TSV is always written first, so the counts survive even
  if the loom step is skipped.
- **`loompy` import failure is caught** and degraded to a warning — it happens at
  the very end of a long run and the TSV already holds the same numbers.
- **No end-to-end CI test.** Unit coverage is
  `isoquant_tests/test_rna_velocity_counter.py` (21 tests: acceptance rules per
  assignment type, fragment layout, dump-clears-state, the three drop cases, loom
  layers, duplicate-row summing). A CI workflow on a 10x sample would need to
  check that the loom opens and the layer sums match the TSV.
- **No UMI deduplication.** Counts are per read, not per molecule, even though
  single-cell modes run UMI filtering. Velocyto counts molecules.

## Files

- `isoquant_lib/quantification/rna_velocity_counter.py` — the counter
- `isoquant_lib/assignment/assignment_aggregator.py` — `_add_grouped_velocity_counter`
- `isoquant_lib/utils/input_data_storage.py` — `SampleData.out_rna_velocity_grouped`
- `isoquant_tests/test_rna_velocity_counter.py` — unit tests
- `docs/output.md` — user-facing output description
- `requirements.txt`, `requirements_tests.txt`, `pyproject.toml` — `loompy>=3.0.7`
