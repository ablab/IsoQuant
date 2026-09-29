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
| `SAMPLE.RNA_velocity_grouped_<strategy>` | TSV: `#cell_id`, `gene_id`, `spliced`, `unspliced`, `ambiguous` |
| `SAMPLE.RNA_velocity_grouped_<strategy>.loom` | Same counts; genes × cells, layers `''` (= spliced), `spliced`, `unspliced`, `ambiguous` |

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

### Why no CLI flag

Deliberate: the output is cheap enough that a flag would cost more than it
saves. Conditions 2 and 3 already restrict it to runs where the output is
meaningful, and guarding a small extra output behind an option makes the
`--analysis` surface larger for no benefit.

The per-read cost is the part that has to be negligible, and it is:
`add_read_info` walks `isoform_matches` and their `match_subclassifications`
once. `IntronRetentionCounter`, already in `global_counter` under
`--analysis exon_quantification`, walks the same matches *and* an inner loop over
isoform intron indices, so velocity is strictly less per-read work than a counter
the pipeline already runs.

The finalize-time loom export is the part that needed fixing rather than gating —
see `create_loom` below. After the windowed write, peak RSS is flat in the cell
count (0.79 GB at 5k genes x 20k cells, 0.65 GB at 80k cells, where a dense
int32 grid would be 0.4 GB and 1.6 GB). What remains is *time*, linear in
genes x cells because loom is dense: measured ~65 s and ~250 MB of file per 1e9
grid entries (4.5 s / 12.8 MB at 5k x 10k, scaling cleanly to 25.9 s / 98.6 MB at
5k x 80k). Extrapolated, a 30k-gene 10x run with 100k cells spends ~3 min and
~750 MB on the export; a 500k-spot spatial run spends ~15 min and ~4 GB. That is
tolerable at the end of a multi-hour run but is the number to revisit if spatial
cell counts grow — the fix then is a sparsity threshold on the loom step, not a
CLI flag, since the TSV is written either way.

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

Tallies are keyed on `(int group id, int gene id)` in three dicts (one per
`SplicingStatus`, dispatched through `self._tallies`), resolved to strings once in
`dump()`. A chromosome's worth of `(cell, gene)` pairs runs to millions of keys
in a real single-cell run; two ints per key instead of two strings is the
difference between a workable and an unworkable worker footprint, and it keeps
`pool.get_str()` off the per-read path. This matches how the feature counters in
`long_read_counter.py` store their group ids.

`dump()` clears all three dicts. It is called once per chromosome today, but a second
call must not re-emit the same counts.

### The loom write (`create_loom`)

Layers are filled **one column block at a time** via `loompy.new()` +
`add_columns()`, not by assigning a sparse matrix to `ds.layers[name]`.

The reason is a bug in loompy itself. Loom is a dense HDF5 format, and
`LayerManager.__setitem__` does intend to fill a sparse layer in windows — but it
computes the window as

```python
window = max(1, 1024**3 // 8 * m.shape[0])   # loompy 3.0.8, layer_manager.py:150
```

which **multiplies where it means to divide** (`1024**3 // (8 * m.shape[0])`).
The result is ~1.3e8 x n_genes, so the subsequent `min(window, remaining)` always
clamps to every remaining column and the whole matrix is densified in one
`toarray()`. Measured at 5k genes x 20k cells: +1.4 GB peak RSS for a single
int64 layer, and the old code wrote four layers. At 30k genes x 500k spots it
would ask for ~60 GB.

`add_columns()` is loompy's streaming API and does not go through that path.
`LOOM_WINDOW_BYTES` (64 MB) caps the dense block per layer, so the write costs
~256 MB across the four layers regardless of cell count. Counts are cast to
`int32` — per-(cell, gene) read counts cannot approach 2^31, and it halves both
the block and the file.

Do not "simplify" this back to `loompy.create(...)` + `ds.layers[x] = matrix`
unless loompy has fixed that expression; the unit test
`test_loom_is_written_in_several_column_blocks` shrinks `LOOM_WINDOW_BYTES` so
the multi-block path is actually exercised, but it will not catch the memory
regression.

## Classification

Splicing status comes from **intronic evidence**, not from whether the read
agrees with the annotation.

### Per-match verdict (`_match_status`)

| Condition | Verdict |
|---|---|
| `match_classification == genic_intron` | unspliced — the read lies inside an intron |
| any event in `INTRON_RETENTION_EVENTS` | unspliced — the read covers intronic sequence |
| `match_classification in {genic, intergenic}` | `None` — resembles no isoform, undecidable |
| otherwise | spliced |

```python
INTRON_RETENTION_EVENTS = {
    intron_retention, unspliced_intron_retention,
    incomplete_intron_retention_left, incomplete_intron_retention_right,
}
```

`fake_micro_intron_retention` is excluded: `is_alignment_artifact()` classes it
as an artifact, and it is the only IR-named event treated as a *minor* error
rather than a major inconsistency. `MatchClassification.undefined` is **not**
undecidable — it means the SQANTI-style classifier had nothing to say, not that
the read resembles no isoform, and such a read still carries events the verdict
can be read off. Note `long_read_counter.INTRON_RETENTION_EVENT_TYPES` holds only
the two complete variants, so velocity defines its own set rather than importing.

**The incomplete variants already carry a 50 bp floor**, so the counter adds no
length check of its own. Both detection sites in `junction_comparator.py` (line
198 for spliced reads, 415 for mono-exonic ones) gate them on
`overlaps_at_least(read_region, intron, params.minor_exon_extension)`, and
`args.minor_exon_extension = 50` (`isoquant.py:1032` — hardcoded, not data-type
dependent, not a CLI option). Verified numerically: the predicate flips at exactly
50 bases of intronic overlap, symmetrically on both sides. The counter could not
apply its own threshold anyway — `MatchEvent` stores `isoform_region` /
`read_region` as **index** pairs plus an `event_info` int, with no coordinates and
no overlap length; recomputing retained bases would mean intersecting
`gene_info.all_isoforms_introns` with the read's exons per read, on the hot path.
Changing the cutoff therefore means changing `minor_exon_extension`, which is
shared with exon-extension and APA logic across the assigner.

(`overlaps_at_least` short-circuits its length test when one interval contains
the other. Read-contains-intron is caught earlier as *full* IR, and
read-inside-intron is excluded at the mono-exon site by
`not contains(intron, read_region)`; a spliced read whose whole span sits inside
an isoform intron can still reach the incomplete-IR branch at the first site, but
such a read is pre-mRNA anyway, which is the verdict it gets.)

### Per-read aggregation (`_read_status`)

Verdicts are collected across matches that name a gene; matches returning `None`
are skipped, so an undecidable match cannot dilute a decided read.

- all unspliced → **unspliced**
- all spliced → **spliced**
- mixed → **ambiguous** — the read is retained against one isoform and mature
  against another, i.e. genuinely compatible with either model. This is what
  velocyto's third category is for, and it keeps such reads out of unspliced.
- none → read dropped

### Accepted assignment types

`MATCHED_ASSIGNMENT_TYPES` (`unique`, `unique_minor_difference`, `ambiguous`,
`inconsistent`, `inconsistent_ambiguous`, `inconsistent_non_intronic`) go through
the verdict above. `noninformative` is accepted **only** when some match is
`genic_intron` — that is the pre-mRNA path, and it is the reason intronic reads
now reach the counter at all (`assign_to_overlapping_genes` sets `assigned_gene`
on every match it builds, so attribution works). `intergenic`, `discarded` and
`suspended` are dropped, as are reads with no barcode and reads with no gene.

Reads with real intron retention are always `inconsistent`/`inconsistent_ambiguous`
in practice (the IR events are in `intronic_major_events`, and
`is_major_inconsistency` is tested before `is_minor_error`), but the verdict is
computed from events for every accepted type rather than inferred from the type.

### History: why not assignment consistency

The original rule was `{unique, unique_minor_difference, ambiguous}` → spliced,
`{inconsistent, inconsistent_ambiguous}` → unspliced, everything else dropped.
It failed in both directions: it dropped intronic reads entirely (they are
`noninformative`), and it counted every annotation disagreement — novel exon,
shifted splice site, alternative TSS — as unspliced.

**Measured on the same SIRV run** (`GROUP12.SC.SIRVs.R10` inputs, two BAMs,
`--mode tenX_v3 --read_group barcode --barcode2spot ...`, reduced annotation,
simulated from mature transcripts so the true unspliced fraction is ~0):

| rule | spliced | unspliced | ambiguous | unspliced % |
|---|---|---|---|---|
| assignment consistency (old) | 1402 | 3013 | — | **68.2** |
| intron-retention events (now) | 3830 | 585 | 0 | **13.3** |

The old accounting was exact — the run reports `unique: 1359`,
`unique_minor_difference: 43`, `inconsistent: 3013`, so every `inconsistent` read
became an unspliced count, driven by the reduced annotation rather than by
biology. The new 585 map exactly onto the 585 `intron_retention` events in that
run; the residue is genuine annotation ambiguity (a read matching a transcript
whose intron is absent from the reduced GTF looks retained against the isoform
that keeps it), not a rule artifact.

### Real-data check (Mouse 10x, chr19)

`Mouse.10x.5k.ONT_cDNA.R10.4.no_trunc.bam` subset to chr19 (115,841 alignments →
49,756 barcoded → 25,056 assigned), gencode vM36 basic, `-m tenX_v3`:

```
spliced=24992  unspliced=64  ambiguous=0   unspliced%=0.3   cells=4948  genes=616
IR events present: unspliced_intron_retention 55, intron_retention 5,
                   incomplete_intron_retention_3 2, incomplete_intron_retention_5 1
classifications:   full_splice_match 19636, mono_exon_match 5062,
                   novel_not_in_catalog 307, novel_in_catalog 59,
                   incomplete_splice_match 47, genic_intron 3
```

What this does and does not establish:

- **Every path fires on real data.** All four IR event types appear, and the
  `genic_intron`/`noninformative` branch fires (3 reads). Loom layers equal the
  TSV column sums.
- **The accounting matches the design exactly.** All 25,056 assigned reads are
  counted, against 24,997 under the old rule — a difference of 59 = the 56
  `inconsistent_non_intronic` reads (now spliced, previously dropped) plus the 3
  `genic_intron` reads (now unspliced, previously dropped).
- **It is not a biological validation.** 98.6% of the assigned reads are
  `full_splice_match` or `mono_exon_match`: this is a curated full-length CI
  benchmark BAM, not raw 10x, so 0.3% unspliced says more about the input than
  about the rule. A real check needs unfiltered data and a velocyto run to
  compare against — **still outstanding**.

## Gene attribution

`_get_gene_id` returns the first `isoform_matches[*].assigned_gene_id` that is
not `None`, mirroring the feature counters. Taking `isoform_matches[0]`
unconditionally is wrong: `match_inconsistent`'s `quick_mode` exit
(`long_read_assigner.py:618`) builds an `IsoformMatch(MatchClassification.genic)`
with **no gene**, which produced rows with an empty `gene_id` column.

For ambiguous reads (`ambiguous`, `inconsistent_ambiguous`) the first named gene
wins outright — no 1/N split across candidate genes, unlike the gene counter.
The same applies to a multigenic intronic read: it is counted once, against the
first gene that names it, even though it sits in an intron of several.

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

- **`finalize()` loads the whole merged TSV into a pandas DataFrame** before
  factorizing it. That is proportional to the number of populated (cell, gene)
  pairs, not to the grid, and it is the one part of the export that is not
  bounded. The TSV is always written first, so the counts survive even if the
  loom step fails.
- **`loompy` import failure is caught** and degraded to a warning — it happens at
  the very end of a long run and the TSV already holds the same numbers.
- **No end-to-end CI test.** Verified manually once (see the SIRV run above:
  both `barcode` and `barcode_spot` looms open, layer sums equal the TSV sums,
  `''` layer equals `spliced`, and no velocity file is produced for the
  `file_name` strategy — plus the Mouse chr19 run above). Unit coverage is
  `isoquant_tests/test_rna_velocity_counter.py` (32 tests: the per-match verdict
  for every IR event type and classification, mixed-verdict ambiguity,
  `genic_intron` acceptance, fragment layout, dump-clears-state, the drop cases,
  loom layers, duplicate-row summing). A CI workflow on a 10x sample would need
  to check that the loom opens and the layer sums match the TSV.
- **No UMI deduplication.** Counts are per read, not per molecule, even though
  single-cell modes run UMI filtering. Velocyto counts molecules.
- **Never compared against velocyto.** Neither validation dataset can settle
  whether the unspliced fraction is biologically right: SIRV is simulated from
  mature transcripts, and the Mouse benchmark BAM is 98.6% full-splice-match. The
  open question is whether including `incomplete_intron_retention_*` (a read end
  running >=50 bp into an intron, which could equally be an unannotated 3'UTR
  extension) helps or hurts on unfiltered data. Dropping those two event types
  from `INTRON_RETENTION_EVENTS` is a one-line change if it turns out to hurt.

## Files

- `isoquant_lib/quantification/rna_velocity_counter.py` — the counter
- `isoquant_lib/assignment/assignment_aggregator.py` — `_add_grouped_velocity_counter`
- `isoquant_lib/utils/input_data_storage.py` — `SampleData.out_rna_velocity_grouped`
- `isoquant_tests/test_rna_velocity_counter.py` — unit tests
- `docs/output.md` — user-facing output description
- `requirements.txt`, `requirements_tests.txt`, `pyproject.toml` — `loompy>=3.0.7`
