# Joint Exon Counts

Region-based exon quantification that groups overlapping annotated exons into
**exon-overlap regions** and emits **N+1 features per region** (N inclusion
variants + 1 region-level exclusion) instead of the per-exon
inclusion/exclusion pair produced by the classic `ExonCounter`.

This counter runs **alongside** the classic `ExonCounter` under the same
switch — both outputs are written. Exon quantification is enabled by
`--analysis exon_quantification` (short alias `ex_quant`); the old `--count_exons`
flag still works but is deprecated. Either path sets the internal `args.count_exons`
flag (see `resolve_analyses()` in `isoquant.py` and `.claude/ANALYSIS_OPTION.md`).
`--analysis exon_quantification` also turns on intron-retention counting
(`args.count_intron_retentions`), which the deprecated flags toggle separately.

## Motivation

The classic `ExonCounter` (`isoquant_lib/long_read_counter.py:644`) treats
every annotated exon independently:

- `profile[i] == +1` → +1 inclusion for exon `i`
- `profile[i] == -1` → +1 exclusion for exon `i`
- `profile[i] == 0`  → ignore exon `i`

For genes with alternative 3'/5' acceptors, cassette variants, or other
overlapping isoform-specific exons, this over-counts exclusions: a read that
clearly picks variant A of an overlap cluster {A, B, C} contributes +1
inclusion to A **and** +1 exclusion to B and C — even though B and C are not
really "excluded" by the read; they're just unselected alternatives.

Joint counts solve this by treating the cluster as a single quantification
unit. The read either selects a variant (+1 to that variant), skips the whole
cluster (+1 to a region-level "exclusion" feature), or is ignored.

## Algorithm

### Region construction (at `GeneInfo` build time)

`GeneInfo.build_exon_overlap_regions(exon_features, exon_property_map)` runs a
strand-aware sweep over the already-sorted `exon_profiles.features`:

```
for each strand bucket (keyed by exon_property_map[i].strand):
    current = None
    for i in indices (sorted):
        exon = features[i]
        if current is None or exon.start > current.end:
            current = new ExonRegion(...)
            regions.append(current)
        else:
            current.end = max(current.end, exon.end)
        current.member_exon_indices.append(i)
        region_map[i] = index_of(current)
```

Singletons (exons with no overlap) form one-member regions.

Result stored on `GeneInfo`:
- `exon_overlap_regions: list[ExonRegion]`
- `exon_overlap_region_map: list[int]` — `region_map[i]` is the index into
  `exon_overlap_regions` for the exon at `exon_profiles.features[i]`

Wired into every `GeneInfo` construction path that builds
`exon_property_map`: `__init__` (line 248) and `deserialize` (line 456).
The alt constructors (`from_models`, `from_model`, `from_region`) set
`exon_overlap_regions = []` / `exon_overlap_region_map = None`.

### Per-read classification (in `JointExonCounter.add_read_info`)

First, pull `read_gene` from the read's first non-null
`isoform_matches[*].assigned_gene`. Reads without an assigned gene are
dropped. `read_gene` becomes part of the counter key so the same physical exon
shared by two overlapping genes gets attributed only to the gene the read was
actually assigned to.

For each `ExonRegion`:

1. **Strand filter** — `if read.strand not in property_map[members[0]].strand: skip`.
2. **Read region states** — `states = [profile[i] for i in members]`.
3. **Any `0` in states** → skip region (insufficient information).
4. **Any `+1` in states** → for every `i` with `states[i] == +1`, skip if
   `read_gene not in property_map[i].gene_ids`; otherwise
   `inclusion_counter[(exon_FeatureInfo.id, read_gene)].inc(group_id)`.
5. **All states `-1`** → skip if `read_gene not in region.gene_ids`;
   otherwise `exclusion_counter[(region.id, read_gene)].inc(group_id)`.

Note: mixed `+1`/`-1` (without any `0`) is treated as case 4 — region
exclusion does **not** fire. That's by design: the read picked at least one
variant, so it isn't "skipping the region".

`ExonRegion.gene_ids` is the union of its member exons' `FeatureInfo.gene_ids`,
pre-computed in `build_exon_overlap_regions` so the exclusion gate is O(1).

## Output

**File names:**
- Ungrouped: `<prefix>.joint_exon.counts.tsv`
- Per group strategy: `<prefix>.joint_exon_grouped_<strategy>.counts.tsv`
- Per-chromosome shards: `get_joint_exon_counts_file(chr_id)` /
  `get_grouped_counts_file(chr_id, "joint_exon", strategy_name)`

**Schema (10 columns):**

```
chr  region_start  region_end  strand  exon_start  exon_end  gene_id  feature_kind  group_id  count
```

- **Inclusion row:** `feature_kind = "inclusion"`, `exon_start`/`exon_end` =
  specific annotated exon coordinates, `gene_id` = the read's assigned gene.
- **Exclusion row:** `feature_kind = "exclusion"`, `exon_start = exon_end = "."`,
  `gene_id` = the read's assigned gene.

The same physical exon or region may emit multiple rows — one per distinct
`(feature_id, read_gene)` pair seen — when overlapping genes share the same
annotated exons. Counter keys are `(feature_id, gene_id)` tuples; rows are
emitted in first-seen order.

Zero-count rows are suppressed (same convention as `ExonCounter`).

**No matrix conversion** — `finalize()` is a no-op. The schema doesn't match
`convert_profile_to_matrix`'s include/exclude pair layout. Linear TSV is the
only output for now; if a matrix view is needed later, write a dedicated
converter rather than coercing this schema into the existing one.

## Files

| File | Role |
|------|------|
| `isoquant_lib/gene_info.py` | `ExonRegion` class; `build_exon_overlap_regions` method; wiring in `__init__` / `deserialize`; safe defaults in `from_models` / `from_model` / `from_region`. |
| `isoquant_lib/long_read_counter.py` | `JointExonCounter(AbstractCounter)` — classifier, dump, no-op finalize. |
| `isoquant_lib/assignment_aggregator.py` | Instantiates `JointExonCounter` ungrouped (next to `ExonCounter`) and grouped per strategy when `args.count_exons` is set (via `--analysis exon_quantification` or deprecated `--count_exons`). |
| `isoquant_lib/input_data_storage.py` | `out_joint_exon_counts_tsv`, `out_joint_exon_grouped_counts_tsv`, `get_joint_exon_counts_file`. |
| `isoquant_tests/test_joint_exon_counter.py` | Unit tests: region grouping, strand separation, inclusion / all-exclusion / zero-state / strand-filter / multi-inclusion / dump format. |

## CLI

No dedicated flag. Exon quantification (`--analysis exon_quantification` / `ex_quant`,
or the deprecated `--count_exons`) sets `args.count_exons`, which toggles both
`ExonCounter` and `JointExonCounter` together. The joint counter is cheap to
run alongside the classic one (single extra pass over
`exon_overlap_regions` per read, no extra I/O).

If the dual output proves noisy in practice, add `--count_joint_exons` and
gate the wiring on it. Deferred until requested.

## Verification

```bash
pytest isoquant_tests/test_joint_exon_counter.py -v
```

End-to-end spot check:

```bash
python3 isoquant.py \
    --reference tests/simple_data/chr9.4M.fa.gz \
    --genedb tests/simple_data/chr9.4M.gtf.gz --complete_genedb \
    --fastq tests/simple_data/chr9.4M.ont.sim.fq.gz \
    --data_type nanopore --analysis quant td exon_quantification \
    -o test_joint_exon -t 2
```

- `*.exon.counts.tsv` must remain byte-identical to a baseline run on the same
  data (regression guard — classic `ExonCounter` is untouched).
- `*.joint_exon.counts.tsv` must exist, be non-empty, and have the 10-column
  header. For genes with no overlapping exons, joint inclusion counts should
  match per-exon inclusion counts in `*.exon.counts.tsv` (summed across
  `gene_id` rows if an exon belongs to multiple genes).

## Edge cases

- **Empty `exon_overlap_regions`** (e.g., `GeneInfo` built via
  `from_region`/`from_model`/`from_models` with `prepare_profiles=False`):
  `add_read_info` returns early.
- **`is_valid` / `is_assigned_to_gene` checks** are reused from
  `ProfileFeatureCounter` so unassigned, ambiguous, or profile-less reads are
  filtered the same way as in `ExonCounter`.
- **Strand stored as concatenated string** on `FeatureInfo` (e.g., `"+-"` for
  exons that appear on both strands across isoforms). Region construction
  uses the strand string as a dict key so `"+"` and `"+-"` form separate
  buckets — intentional, mirrors `set_feature_properties`' encoding.

## Out of scope

- Joint counters for introns / splice junctions / intron retentions.
- Matrix / MTX conversion for the joint output.
- Reporting which annotated variant a read matched in the per-read TSV.
- Region merging across opposite-strand exons on the same coordinates.
