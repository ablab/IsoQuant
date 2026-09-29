# Exon Counts (per-exon usage)

`exon_counts.tsv` is produced by `ExonUsageCounter`
(`isoquant_lib/quantification/long_read_counter.py`) when exon quantification is on
(`--analysis exon_quantification` / `ex_quant`, deprecated `--count_exons`). It replaced
the region-based `JointExonCounter` shipped in 4.0.0, which was removed after
discussion #426 showed it was unusable for exon skipping:
- people care about individual exons, so PSI values per region were meaningless;
- a read counted towards a region only if it was informative (non-0) for **every**
  member exon. In long overlap chains (e.g. linked by intron-retention exons) most reads
  were dropped and exclusion was badly underestimated;
- "exclusion" also covered shifted splice sites and intron retention, and exons past a
  polyA tail (profile −2) were counted as excluded.

Runs alongside `ExonSpliceSiteCounter` (`exon_splice_site_counts.tsv`, ScisorseqR format,
see `EXON_SPLICE_SITE_COUNTS.md`) and the legacy `ExonCounter` (`--old_exon_count_format`).

## Per-(read, exon) states

Every annotated exon (`gene_info.exon_profiles.features`) is classified independently
per read. States are mutually exclusive for one (read, exon) pair, but one read contributes
to many exons.

| state | rule |
|---|---|
| full | an internal read block matches both exon borders within `delta`; **or** the first/last read block matches the exon's internal splice site, its free end lies inside the exon (`delta` outward tolerance) and the exon is the leftmost/rightmost exon of some transcript (no TSS/polyA anchoring required) |
| left / right | same as the terminal case, but the exon is not transcript-terminal on that side: only its left/right splice site is confirmed |
| skip | all split-exon segments of the exon are −1 (exon lies within a read intron) |
| alt | a block overlaps the exon (some segment +1) but none of the above holds: other splice site, block running past the exon, retained intron |
| — | otherwise (exon not reached by the read) nothing is counted |

Filters: `ProfileFeatureCounter.is_valid` / `is_assigned_to_gene` (unique gene only),
read's `assigned_gene` must be in the exon's `gene_ids`, read strand in exon strand,
reads with a single block are ignored (no splice sites). Blocks are `corrected_exons`
(fallback `exons`). A splice site shared by several exons gives each of them a half-inclusion.

## Engine: split-exon segments

`GeneInfo.get_exon_usage_index()` lazily builds and caches an `ExonUsageIndex`
(`isoquant_lib/gene_info.py`) from `split_exon_profiles.features` (non-overlapping atomic
segments, monotone starts and ends): exon → segment range, segment → exons, and
per-exon left/right transcript-terminal flags (from `all_isoforms_exons`).

Per read: bisect the segments to the read span, recompute the split profile on that
slice with `NonOverlappingFeaturesProfileConstructor` (comparator
`overlaps_at_least_when_overlap`, `minimal_exon_overlap`, **no polyA positions → no −2**),
then visit only exons touching informative segments. Skip requires the exon to be fully
inside the slice with all segments −1; any +1 triggers the coordinate check
(`_inclusion_state`). The split profile is built by the assigner too, but it is not stored
on `ReadAssignment` and counting runs on deserialized assignments, so it is recomputed
(cheap: O(log S + blocks + features under the read)); nothing is serialized.

Note: a segment only slightly overlapped by a block (< `minimal_exon_overlap`) stays 0 in
the split profile, not −1, so such an exon is neither skipped nor included.

## Output

```
chr start end strand flags gene_ids group_id include_counts exclude_counts n_full n_left n_right n_alt
```
- first 9 columns = the profile-counter layout (`FeatureInfo` + group + include/exclude),
  `include = full + left + right`, `exclude = skip`; `gene_ids` = the read's assigned gene
  (one row per (exon, gene); key uses coordinate-based `FeatureInfo.id`, so rows merge
  across GeneInfo builds);
- all-zero rows suppressed;
- `finalize()` is inherited from `ProfileFeatureCounter`: grouped matrix / MTX conversion
  works. `_load_profile_linear` (`convert_grouped_counts.py`) uses the first 9 columns and
  ignores the extra ones.

PSI = include / (include + exclude); strict: n_full / (n_full + exclude); usage fraction:
add n_alt to the denominator.

## Tests

`isoquant_tests/test_exon_usage_counter.py` — skip inside an overlap chain, alt 3′ variant,
intron-retaining exon as alt, half inclusion, terminal exons, shared splice site,
exons beyond the read end not counted, uninformative reads, dump + matrix loading.
