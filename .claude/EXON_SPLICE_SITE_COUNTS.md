# Exon Splice-Site Counts

Region-based exon **splice-site** quantification: for each read, records whether it
demonstrates a candidate exon's splice sites **fully**, **only on the left**, or
**only on the right**, plus per-region **exclusion** and **ambiguous** outcomes.
Runs under `--count_exons` alongside `JointExonCounter` (default `exon_counts`) and
`IntronCounter` (`splice_junction_counts`). Originally specced (single-cell only) as
`AlternativeExonCounter` / `alternative_exon`; generalized and renamed.

## Output

Plain TSV (no gzip — `finalize()` is a no-op; per-group rows scale like the other
grouped counters). Ungrouped (bulk included) + one file per grouping strategy.
- Ungrouped: `SAMPLE.exon_splice_site_counts.tsv`
- Per grouping strategy: `SAMPLE.exon_splice_site_grouped_<strategy>_counts.linear.tsv`
- Per-chr shards: `get_exon_splice_site_counts_file(chr_id)` /
  `get_grouped_counts_file(chr_id, "exon_splice_site", strategy)`

**Schema (one row per feature and group):**
`region_gene_candidate  n_full  n_left  n_right  group_id
[ read_ids_full read_ids_left read_ids_right ]`
- `region_gene_candidate = {chr}_{rstart}_{rend}_{strand}__{gene}__{candidate}`;
  candidate = `{chr}_{estart}_{eend}_{strand}` | `exclusion` | `ambiguous`.
- exclusion/ambiguous rows: count in `n_full`, `n_left=n_right=0`.
- `group_id`: read group name; `NA` when ungrouped/bulk. With the `barcode`
  strategy this is the cell barcode.
- `read_ids_*`: optional per (feature, group, bucket) read-id lists, only with
  `--emit_read_ids`.

Bucket = dict `group_id -> [count, [read_ids]]`. `misc/exon_splice_site_to_group_lists.py`
converts this per-group format to the per-molecule group-list format (`groups_*`
columns, one entry per read) for downstream cell/cell-type aggregation.

## Region construction (`GeneInfo.build_exon_splice_site_regions`)

`ExonSpliceSiteRegion` per gene, built **lazily** from the counter and cached on
`gene_info._exon_splice_site_regions` (zero cost when `--count_exons` off, nothing
serialized). Per gene: union of `all_isoforms_exons`; **drop intron-retention
exons** (exon spanning an annotated gene intron end-to-end); overlap-group on the
same strand → connected components; precompute per-candidate left/right
edge-uniqueness (needed for half-inclusion).

## Classification (`ExonSpliceSiteCounter._classify_region`)

Per read (blocks = `corrected_exons`, introns = `junctions_from_blocks`), strand
must match region strand. Anchoring: 5'/TSS from `terminal_site_match_left/right`
events (strand-aware), 3'/polyA from `read_assignment.polyA_found`. Precedence:
1. **full** — a single block matches both edges of exactly one candidate within
   `delta`, with ≥1 splice edge and both edges demonstrated (interior splice site,
   or anchored terminal); >1 candidates → ambiguous.
2. else **exclusion** — a read intron interior `(start+margin, end-margin)`
   contains the whole region envelope (`DEFAULT_EXCLUSION_MARGIN = 50`, code-only).
3. else (read overlaps region) **unique half-inclusion** — exactly one candidate
   with a demonstrated, unique left (or right) edge and the other edge not reached.
4. else **ambiguous**.

No re-dedup: counting is post-UMI-dedup, so each read is one molecule. Group id is
`read_group_ids[group_index]` (0/`NA` when ungrouped).

## Files / CLI

| File | Role |
|------|------|
| `isoquant_lib/gene_info.py` | `ExonSpliceSiteRegion`, `build_exon_splice_site_regions` |
| `isoquant_lib/long_read_counter.py` | `ExonSpliceSiteCounter(AbstractCounter)`, `DEFAULT_EXCLUSION_MARGIN` |
| `isoquant_lib/assignment_aggregator.py` | ungrouped + per-strategy grouped instances under `--count_exons` |
| `isoquant_lib/input_data_storage.py` | `out_exon_splice_site_counts_tsv`, grouped path, `get_exon_splice_site_counts_file` |
| `isoquant.py` | `--emit_read_ids` (read-id columns); no `--exclusion_margin` (code-only) |
| `isoquant_tests/test_exon_splice_site_counter.py` | unit tests |

## Verification

```bash
pytest isoquant_tests/test_exon_splice_site_counter.py -v
```
