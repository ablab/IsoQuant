# Fusion Detection

Optional fusion-gene discovery stage that runs **after** the normal isoform
pipeline when the user requests it via `--analysis fusion` (deprecated alias:
`--fusion`). Lives on the `fusion` branch (PR ablab/IsoQuant#392) and is gated
behind a single analysis value — no impact on the default pipeline. Both the
new value and the deprecated flag set the internal `args.fusion` flag in
`resolve_analyses()` (see `.claude/ANALYSIS_OPTION.md`), which the trigger block
reads; nothing downstream changed.

## Quick run

```bash
python isoquant.py \
  --reference ref.fa --genedb ann.gtf --complete_genedb \
  --bam sample.bam --data_type nanopore \
  -o out --analysis fusion
```

Output: one TSV per input BAM named `out/fusion_<bam_basename>.tsv`
(see [Output format](#output-format) below).

Requires `mappy>=2.24` and `intervaltree>=3.0` (already in `requirements.txt`).

## Files

| File | What it holds |
| --- | --- |
| `isoquant_lib/fusion_detector.py` | `FusionDetector` — main class. BAM scan, breakpoint estimation, soft-clip realignment, clustering, validation orchestration, TSV reporter. |
| `isoquant_lib/fusion_validator.py` | `FusionValidator` — biotype/multicopy/frequency filters, confidence scoring, duplicate merging. |
| `isoquant_lib/fusion_metadata.py` | `FusionMetadata` — gating, gene assignment from canonical key, per-side biotype computation, key remapping. |
| `isoquant_lib/genomic_interval_index.py` | `GenomicIntervalIndex` — IntervalTree wrapper around gffutils for O(log n) gene/exon coordinate lookups. |
| `isoquant_tests/test_fusion_*.py` + `test_genomic_interval_index.py` | 161 unit tests (1 legitimately skipped). |
| `isoquant.py` | `--analysis fusion` / deprecated `--fusion` flag (`parse_args`) resolved to `args.fusion` in `resolve_analyses`, `get_bam_files_from_samples`, `run_fusion_detection_on_bams`, the post-isoform trigger block in `run_pipeline`. |

## Pipeline integration

`run_pipeline(args)` in `isoquant.py`:

1. Normal stages: barcode calling → reference prep → GTF→DB → mapping → isoform
   assignment → count aggregation.
2. **If `args.fusion` is set:**
   - Collect BAMs across all samples via `get_bam_files_from_samples`.
   - Create **one** `FusionDetector(bam_files[0], args.genedb,
     reference_fasta=args.reference)` — the interval index and module-level
     caches are reused across BAMs.
   - For each BAM: `fd.bam_path = bam`, `fd.clear_state()`, `fd.detect_fusions()`,
     `fd.report(output_path=...)`. Errors per-BAM are logged but don't abort the run.

No per-BAM parallelism — `--threads` is ignored by fusion detection itself.
(Multi-sample throughput is sample-by-sample.)

## Algorithm

### 1. Breakpoint discovery — `detect_fusions()`

Iterate primary alignments in the BAM (skips secondary, supplementary,
unmapped). Per read:

- **Quality filters** (`_passes_read_filters`): aligned length
  ≥ `min_al_len_primary=50`, MAPQ ≥ `min_sa_mapq=10`. Also detects a single
  terminal soft-clip ≥ 50 bp (`detect_softclip`).
- **SA-tag path** (`_process_sa_entries`): parse the `SA` tag
  (`parse_sa_entries`), filter SAs by MAPQ and aligned length, keep top 3 by
  (mapq, aligned_len). For each SA mate:
  - Estimate left/right breakpoints (`estimate_breakpoint`) using the
    primary's soft-clip side and the SA's CIGAR to anchor the right side
    at start vs end.
  - Deduplicate (chr1, pos1//jitter, chr2, pos2//jitter) at
    `jitter_window=50`.
  - Assign genes via `assign_fusion_gene_cached` for both sides.
  - If genes differ, `record_fusion(...)` adds the read to the candidate set
    keyed by alphabetically-sorted `"GENE_A--GENE_B"`.
- **Soft-clip realignment path** (`realign_softclip`): if the primary has a
  strong terminal clip and `mappy.Aligner` is available, realign the clipped
  subsequence against the reference, accept the first hit with mapq ≥ 20 and
  aligned length ≥ 50, then `record_fusion` on the implied breakpoint pair.

### 2. State after BAM scan

Per `FusionDetector`:

```
fusion_candidates       Dict[fusion_key, Set[read_name]]
fusion_breakpoints      Dict[fusion_key, Dict[(c1,p1,c2,p2), int]]   # observation counts
fusion_assigned_pairs   Dict[fusion_key, Dict[read_name, (left_gene, right_gene)]]
fusion_read_scores      Dict[fusion_key, Dict[read_name, (left_score, right_score)]]
fusion_metadata         Dict[fusion_key, dict]                        # filled by build_metadata
```

`fusion_key` is the alphabetically-sorted `"GENE_A--GENE_B"`. Sorting avoids
duplicate entries for left/right orientation; the original (left, right)
pairs are preserved per-read in `fusion_assigned_pairs` for later salvage.

### 3. Metadata build — `FusionMetadata.process_all`

For each fusion key (`_process_fusion_candidate`):

- **`_passes_all_gates`**: require `support ≥ min_support` (default 1 in
  build phase, raised to 2 later); breakpoints exist; clustering returns a
  consensus.
- **`cluster_breakpoints(bp_counts, window=2000)`** in `FusionDetector`:
  brute-O(n²) per chromosome pair — for each candidate center `i`, sum
  weights of items within `±window` on both `p1` and `p2`; pick the densest
  square; consensus is the weighted centroid. *n* per fusion is typically tens.
  (Previous two-pointer implementation was broken — see [history](#history--recent-fixes).)
- **`_assign_genes`**: take left/right gene names from the alphabetically-sorted
  key directly (avoiding majority-voting corruption from inconsistent per-read
  assignments). Per-side scores = mean of per-read scores recorded earlier.
- **`_update_fusion_key_mappings`**: if normalization changes a partner name,
  re-key metadata + candidates + breakpoints + assigned pairs.
- **`_compute_biotypes`**: per side, look up biotype via
  `get_gene_biotype(symbol, chrom, pos)`.

### 4. Validation pipeline — `validate_candidates`

Wrapper inside `FusionDetector.validate_candidates(min_support=2, ...)`:

1. `build_metadata(min_support=1)` → fills `fusion_metadata` (see §3).
2. `_clear_read_level_data()` → drop `fusion_assigned_pairs` and
   `fusion_read_scores` (no longer needed; saves RAM on big datasets).
3. `FusionValidator.filter_based_on_biotype()` — drop fusions where either
   partner's biotype is not in `{protein_coding, IG_{C,V,D,J}_gene,
   TR_{V,D,J,C}_gene}`.
4. `FusionValidator.filter_multicopy_artifact_pairs()` — mark invalid when
   **both** partners are multicopy artifacts (RPL/RPS/MRPL/MRPS, HIST/H1-/H2A/H2B/H3-/H4-,
   HLA-, MICA/B, KRT, OR, ZNF, DEFA/B, IGH/IGK/IGL, TRAV/TRBV/TRGV/TRDV,
   MUC, AMY, CYP).
5. `FusionValidator.apply_frequency_filters()` — mark invalid when a partner
   gene shows up across many fusions in this sample: threshold 2 for
   ribosomal/histone, 4 otherwise.
6. Per-surviving fusion gates 1–5:
   - **Gate 1**: support ≥ `min_support=2`.
   - **Gate 2**: consensus breakpoint exists.
   - **Classify** (`_apply_classification_and_filters`): `intragenic`
     (antisense/divergent on same locus → invalid), `cis-SAGe`
     (same chromosome, distance ≤ `max_intra_chr_distance`),
     `intergenic` (either side unassigned), else `canonical`.
   - **Gate 3**: mitochondrial filter (`_is_mitochondrial_candidate` — chrM/MT
     contigs, gene names starting with `MT-`). If passing, run
     `_attempt_reconstruction_and_realignment`: fetch exons within ±100 bp of
     each breakpoint, concatenate ref sequence, realign with mappy. Failure
     to reconstruct or realign cleanly invalidates the fusion.
   - **Gate 4**: cis-SAGe policy enforcement (`allow_cis_sage=True` by default).
   - **Gate 5**: confidence threshold `0.30` with `support < 2` cutoff.

### 5. Confidence score — `FusionValidator.confidence`

```
support_norm = min(support, 10) / 10
recon        = 0.2 if reconstruction_ok else 0
realign      = min(hits, 3) * 0.1 + min(best_hit_mapq, 30) / 300
priors       = -0.30 if either side multicopy artifact, +0.20 if either side driver gene
confidence   = clamp(support_norm + recon + realign + priors, 0, 1)
```

Driver whitelist (`FusionValidator.is_driver_gene`): BCR, ABL1, RUNX1,
RUNX1T1, PML, RARA, ETV6, NUP98, NUP214, KMT2A, MLLT3, EWSR1, ERG, FGFR1,
RET, ROS1, ALK, TMPRSS2, FLI1, MYH11, CBFB.

### 6. Gene assignment — `assign_fusion_gene(chrom, pos, window=1000)`

Fast path (any gene whose exon contains `pos`) → `_score_exonic_genes` picks
the best by tuple `(is_protein_coding, score, -boundary_dist, -exon_dist,
-body_dist, -min(glen, 500k), name)`. Slow path (no exonic hit) →
`_score_intronic_genes` with the same tuple over gene-body proximity.
`_compute_gene_score` is binary: 1.0 inside gene bounds, 0.0 otherwise.

### 7. Gene name normalization

- `resolve_gene_name`: ENS* or other IDs → gene symbol via cached
  `db[name_or_id]` lookups, with `_fallback_gene_lookup` if direct lookup fails.
- `canonical_locus_name` + `ANTISENSE_SUFFIX_RE` (`-(AS\d+|DT|DIVERGENT|NAT)$`):
  collapse antisense/divergent transcript names to the canonical locus.
- `normalize_gene_label`: resolves → collapses → if still unresolved
  (`None`, starts with `ENSG`/`RP11`) returns `"intergenic"`.

## Caches

Module-level dicts (shared across `FusionDetector` instances within a Python
process) with FIFO eviction at fixed caps:

| Cache | Purpose | Cap | Evicts |
| --- | --- | ---: | ---: |
| `_CONTEXT_CACHE` | `(chrom, pos) → (gene_name, region_type)` for `_context_query` | 200 000 | 20 000 |
| `_GENE_ASSIGNMENT_CACHE` | `(chrom, pos) → (gene_name, score)` for `assign_fusion_gene_cached` | 500 000 | 50 000 |
| `_ALIGNER_MAP_CACHE` | `seq → tuple(hits)` for mappy realignment | 10 000 | 1 000 |
| `_CIGAR_CACHE` | `cigar_string → aligned_length` | unbounded | — |

Per-instance caches:

| Cache | Purpose | Population |
| --- | --- | --- |
| `self.exon_cache` | `gene_id → [(start, end), ...]` for `_get_cached_exons` | **lazy** in `_get_cached_exons`. `_build_exon_cache` exists as an eager pre-populate (call is commented in `__init__` — uncomment to switch back to eager). |
| `self._resolved_name_cache` | per-detector `name_or_id → symbol` | lazy in `resolve_gene_name` |
| `self._symbol_biotype_cache` | per-detector `symbol → biotype` | lazy via `_build_symbol_biotype_index` (one-shot) |
| `self.interval_index` | `GenomicIntervalIndex` — interval trees of genes + exons by chromosome | built once in `__init__` |

## Output format

`out/fusion_<bam_basename>.tsv`, tab-separated, columns:

```
LeftGene  LeftBiotype  LeftScore  LeftChromosome  LeftBreakpoint
RightGene RightBiotype RightScore RightChromosome RightBreakpoint
SupportingReads  FusionName  Class  Valid  Confidence  Reasons
```

Default report filters (in `report()`): `confidence ≥ 0.3`,
`class in ("canonical", "cis-SAGe")`, `only_valid=False` (invalid entries are
included for debugging, with the failure reason in the last column).

`_merge_fully_identical` runs first to collapse exact duplicate fusions
(same genes + same consensus breakpoint).

## Tests

`pytest isoquant_tests/test_fusion_*.py isoquant_tests/test_genomic_interval_index.py`
→ 161 passed, 1 skipped.

Breakdown:

- `test_fusion_detector.py` — 40 tests: init, CIGAR parsing, softclip
  detection, gene/biotype caches, breakpoint estimation, SA parsing,
  `record_fusion`, mitochondrial filter, mocked `gffutils.FeatureDB` +
  `GenomicIntervalIndex`.
- `test_fusion_metadata.py` — 35 tests: gating, gene assignment from sorted
  key, biotype computation, key remapping, supporting-reads merging.
- `test_fusion_validator.py` — 56 tests: biotype whitelist, multicopy
  patterns, frequency filters, driver genes, confidence formula,
  `_merge_fully_identical`, `_remove_discarded_fusions_internal`.
- `test_genomic_interval_index.py` — 30 tests + 1 skipped (the fixture has no
  position with multiple overlapping genes).

End-to-end: `python3 isoquant.py --test` exercises the full pipeline on the
chr9 toy data (no fusion, since `--test` doesn't request `--analysis fusion`).

## History / recent fixes

Recent local commits on the `fusion` branch (unpushed at time of writing —
file:line references in this doc may shift after rebase):

```
39f3be5 cosmetic cleanup: confidence float format, dangling comments, trailing blanks
eb660e1 fix cluster_breakpoints: brute O(n^2) sweep instead of broken two-pointer
c70b89e remove unused fusion graph hooks from intron_graph
8ba96e0 remove dead IntervalTree None-guards in genomic_interval_index
6d86c45 remove dead mappy guard with stale requirements_fusion.txt reference
63be9e9 make exon cache lazy by default, keep eager build code commented
13d65b4 fix _build_exon_cache to actually populate self.exon_cache
5745e3c gate fusion detection behind --fusion flag
0b4ddaf fix README test command to use repo-relative path
```

Notable behavior changes from those commits:

- **Fusion must be requested**. Previously fusion ran on every pipeline
  invocation (auto-trigger after isoform detection). Now: not requested → no
  fusion. The fusion-only early-return path was deleted; requesting fusion runs
  isoform *and* fusion in order. Originally gated behind `--fusion`; that flag is
  now the deprecated alias for `--analysis fusion` (both set `args.fusion`).
- **Exon cache is lazy by default**. Eager `_build_exon_cache` walks every
  gene in the DB; for GENCODE that's tens of thousands of unused queries when
  zero fusions surface. Now `_get_cached_exons` populates on demand. The
  eager function is still present and correct — uncomment
  `self._build_exon_cache()` in `__init__` to re-enable.
- **`cluster_breakpoints` is brute O(n²)**. The previous two-pointer
  implementation had two bugs: (a) the `k` pointer only excluded items with
  *smaller* p2, never larger p2 outliers; (b) `k` was monotonic, so removals
  were permanent even when a later `i` should have brought items back.
  Symptoms: under-counted cluster support → false-negative fusions; weighted
  centroid smeared by outliers → wrong gene assignment downstream. n per
  fusion is small (tens), so the O(n²) cost is negligible.

## Known issues / TODO

Open from the PR #392 review (not yet fixed in any commit):

- **Deferred dead code in `fusion_detector.py`**: ~5 `self.interval_index is
  not None` branches with a `gffutils` fallback `else:` arm that's unreachable
  (intervaltree is hard-imported at the module top). Cleanup would delete
  ~30 lines.
- **No multithreading**. Per-BAM serial. `--threads` only flows into mapping,
  not fusion. For multi-sample runs this is wasted parallelism.
- **`_compute_gene_score` is binary {0, 1}**. Could distinguish "inside exon"
  vs "inside intron" vs "near boundary" vs "outside" with continuous scoring.
- **`realign_softclip` takes the first good hit**. No ranking, no
  multi-mapping handling.
- **Per-side gene assignment uses majority-voted reads only at metadata
  build**. `_assign_genes` derives gene names from the sorted fusion_key —
  if upstream `record_fusion` got the wrong gene, that error persists.
- **No chromosome blacklist**. Only chrM/MT is hardcoded; users can't pass
  `--discard_chr` to fusion detection.
- **cis-SAGe distance threshold defaults to `None`** (no filter).
- **Per-BAM `FusionDetector` reuse skips the interval-index rebuild**, but
  the module-level caches keep growing across BAMs (until FIFO eviction
  kicks in). Hot keys from BAM 1 may evict valid entries before BAM 2 needs
  them. Worth measuring before optimizing.

Roadmap from the original author (per PR discussion): the eventual plan is to
expose a larger `--fusion_algorithm` flag selecting different detection
strategies; today's `--analysis fusion` is the first implementation.
