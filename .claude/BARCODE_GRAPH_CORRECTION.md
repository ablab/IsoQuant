# Graph-based barcode correction

Branch `badger`. Ported from [Badger](https://github.com/algbio/Badger) (student project,
University of Helsinki). Local reference checkout: `/home/andreyp/Badger`.

## Why

IsoQuant's original barcode calling seeds k-mers against the whitelist and picks the best SSW
alignment. That works for a short whitelist derived from short-read data, but with the stock
10x whitelists (737k for v2, ~3M for v3) it **collapses into exact matching**:

| # | Where | Mechanism |
|---|---|---|
| 1 | `callers/tenx.py` (`bc_count > 100000` branch) | `min_score = 16 == BARCODE_LEN_10X` with SSW `match_score=1, mismatch_penalty=1` — only a perfect 16 bp match can reach that score. Plus `min_matching_kmers=2`, `max_barcodes_hits=10`, `score_diff=1`. |
| 2 | `common.py find_optimal_kmer_size` | k grows with whitelist size; 16 bp × 3M barcodes lands at k≈11, leaving only 6 k-mers per barcode. A substitution at positions 5–10 destroys all six, so the true barcode is never even a candidate. |
| 3 | all indexers, `hits_delta=1` | candidates more than one shared k-mer behind the best are dropped before alignment. |

## Benchmarks

Measured on the real CI datasets. "stock" whitelists are the 7.4M `3M-3pgex-may-2023` list
with the dataset's true cell barcodes prepended (`10x.stocklike_7.4M.txt`,
`wl_stocklike_hipp.txt`) — i.e. what a user gets by passing the stock 10x whitelist instead
of a short-read-derived one.

| dataset | cells | whitelist | whitelist matching | graph |
|---|---|---|---|---|
| Mouse 10x ONT StereoQ, 600k reads | 5000 | matched (5K) | 87.35% / 99.73% | 86.55% / **99.84%** |
| | | stock (7.4M) | 74.36% / 99.60% | **86.55%** / **99.84%** |
| Mouse 10x ONT cDNA R10.4, 304k reads | 5000 | matched (5K) | 89.63% / 99.80% | 88.64% / **99.90%** |
| | | stock (7.4M) | 76.27% / 99.58% | **88.64%** / **99.90%** |
| Mouse 10x concat, 800k reads (`tenX_v3_split`) | 8975 | matched (8975) | 86.80% / 99.80% | **87.24%** / 99.80% |
| | | stock (7.4M) | 74.31% / 99.63% | **87.24%** / **99.80%** |

(recall / precision)

Three things this establishes:

1. **Graph correction is exactly insensitive to whitelist size.** The matched and stock rows
   are identical down to the raw counts, because the whitelist only filters candidate
   centers. Whitelist matching loses 12–13 points across the same span.
2. **With a stock whitelist graph wins by ~12 points of recall**, and gains precision too.
3. **With a correctly sized whitelist graph costs ~1 point of recall** on the single-molecule
   datasets (and *gains* 0.44 in split mode). Hence `auto`: whitelist matching below
   `LARGE_WHITELIST_SIZE`, graph above.

Split mode was the open risk — `SplittingBarcodeDetectionResult.filter()` keeps every valid
detection, and in raw mode every window is valid, so wrong-strand duplicates could have been
emitted. They are not: 1,343,140 molecules emitted against 1,359,897 true ones, versus
1,231,545 for whitelist matching, i.e. closer to truth rather than inflated. The structural
`_is_consistent_detection` check is what holds the line.

### Datasets this does *not* apply to

`Mouse.10x.v3.3M.full` (3,982,878 distinct true barcodes over 10M reads) and
`Mouse.10x.v2.737K.full` (552,857) simulate reads from the *entire* whitelist at ~2.5 reads
per barcode. A flat count distribution has no knee, and `MIN_CENTER_COUNT` excludes nearly
everything — an accidental run against a mismatched whitelist picked 1038 centers out of 5000
requested. These are whitelist-matching stress tests, not single-cell experiments; their ~58%
recall baselines are not a target this feature can move.

Note also that the 5000 cells in `10xMultiome_5K.tsv` come from the 10x ARC/Multiome
whitelist, which is not among the stock lists under `simulation/ref/10x/` (only 458/5000 are
in `3M-3pgex-may-2023.txt`) — hence the synthesised stock-like whitelist. If the real
`737K-arc-v1` list is ever added, point the CI configs at it instead.

`custom_sc` is excluded by `supports_graph_correction()`, so `Mouse.10x.custom_sc.{large,small}`
(the same 5k reads through an MDF) cannot use it yet. Given the numbers above, extending it
there looks worthwhile.

## The idea

Stop searching the whitelist per read. Instead:

1. **Extract** the barcode window verbatim (no matching).
2. **Count** each distinct window over the whole run.
3. **Select centers** — the barcodes that were actually sequenced — from the count
   distribution, keeping only those present in the whitelist.
4. **Correct** every other observed barcode onto a center by edit distance.

The whitelist becomes a membership filter over a few thousand candidate centers, never a
per-read search space. Errors are resolved by evidence from the data rather than by
alignment score against millions of decoys.

## Pipeline flow

```
reads ──> TenXBarcodeDetector(whitelist_matching=False)
              │  raw 16 bp window, no whitelist index built at all
              ▼
    aux/<prefix>.raw_barcodes_<i>.tsv        (read_id, barcode, UMI, ...)
              │
              ▼  correct_barcodes.py
    count ──> select centers ──> cluster ──> rewrite
              ▼
    <prefix>.barcoded_reads_<i>.tsv          (+ trailing raw_barcode column)
    <prefix>.barcode_correction.stats
```

`barcoded_reads_<i>.tsv` keeps its meaning (corrected barcode in column 1), so
`dataset_processor.split_read_barcode_table` (`read_column=0, group_columns=(1,2)`) and
everything downstream is untouched. The extra `raw_barcode` column is ignored downstream.

## Code map

| File | Contents |
|---|---|
| `isoquant_lib/barcode_calling/barcode_graph.py` | `BarcodeGraph`, `bounded_distance`, `estimate_cell_number` |
| `isoquant_lib/barcode_calling/correct_barcodes.py` | the pipeline stage: streaming count pass, rewrite pass, stats |
| `callers/tenx.py` | `whitelist_matching` flag, `_extract_raw_barcode`, `_make_umi_result` |
| `isoquant.py` | `--barcode_correction` resolution (`_resolve_barcode_correction`), `correct_sample_barcodes` |
| `modes.py` | `BarcodeCorrectionMethod`, `LARGE_WHITELIST_SIZE`, `IsoQuantMode.supports_graph_correction()` |

## Options

| Option | Default | Meaning |
|---|---|---|
| `--barcode_correction {auto,whitelist,graph}` | `auto` | `auto` → `graph` when the whitelist exceeds `LARGE_WHITELIST_SIZE` (100k), `whitelist` otherwise |
| `--n_cells` | estimated | expected number of cell-associated barcodes |
| `--n_cells_interval` | 25 | % by which the selected count may differ from `--n_cells` |
| `--barcode_graph_threshold` | 1 | maximal edit distance for a graph edge |
| `--barcode_graph_rounds` | 2 | graph hops barcodes are corrected over |
| `--barcode_graph_impl {centers,full}` | `centers` | see below; results are identical |

Supported modes: `tenX_v2`, `tenX_v3`, `tenX_v2_split`, `tenX_v3_split`, `visium_5prime`
(`IsoQuantMode.supports_graph_correction()`). Stereo-seq is deliberately excluded — every
spot barcode there is genuine, so count-based cell selection does not apply.

## Algorithm details

**Raw extraction.** `whitelist_matching=False` skips index construction entirely and returns
`sequence[r1_end+1 : r1_end+1+16]`. Two things differ from the whitelist path and matter:

- The window is **exactly** `BARCODE_LEN_10X`; the whitelist path deliberately takes 17 bases
  for alignment slack. The graph needs fixed-length input.
- `find_barcode_umi` must **not** short-circuit on the forward strand. In raw mode any
  R1 hit yields a "valid" result, so the early return would never let the reverse strand be
  considered.
- **Strand selection is structural, not score-based** — see below. This is the single most
  important detail in the raw path.
- R1 anchoring is relaxed (`RAW_R1_MIN_SCORE=9`, `RAW_TERMINAL_MATCH_DELTA=4` vs 11/2), as in
  Badger: more reads yield a window, errors are fixed downstream rather than rejected.
  Measured: tightening these back to 11/2 *lowers* recall (64.31% vs 66.49% at the time of
  the measurement), so the relaxation is doing real work.
- `_find_barcode_near_polyt` is unavailable — without a whitelist there is nothing to anchor
  against, so "no R1" means "no barcode".

### Gotcha: strand selection dominates raw-mode accuracy

Whitelist mode picks the orientation on which the barcode actually matched — strong evidence.
Raw mode has no such signal: every window is accepted and `BC_score` is always 0, which leaves
`more_informative_than` comparing polyT positions, i.e. deciding close to at random. A wrong
strand yields a reverse-complemented window that no amount of graph correction can rescue.

This was measured, not guessed. Before the fix, of the 102,696 reads that whitelist matching
got right and graph correction did not, **102,590 were wrong-strand picks** — and their raw
windows sat 7–11 edits from the true barcode, i.e. unrelated sequence rather than a
correctable error.

`_raw_strand_score` / `_pick_raw_strand` rank an orientation by molecule structure instead:
R1 and polyT both present, `0 < polyT - r1 <= MAX_R1_POLYT_DISTANCE`, and the span as close as
possible to `BARCODE_LEN_10X + UMI_LEN`. Recall went 66.49% → **86.55%** on the real dataset.
Regression tests: `TestTenXRawExtraction::test_strand_is_chosen_by_molecule_structure`,
`::test_structural_score_prefers_a_consistent_layout`.

**Center selection** (`select_cluster_centers`). Sort observed barcodes by count; take
`cutoff = max(mean(top n_cells counts) / 5, 5)`; walk down keeping whitelisted barcodes until
the count drops below the cutoff or `n_cells * (1 + interval/100)` is reached; then pad to
`n_cells * (1 - interval/100)`. Without `--n_cells`, `estimate_cell_number` finds the knee of
the log-log count-versus-rank curve (the point furthest below the chord joining its ends).

**Clustering.** Breadth-first from the centers. A barcode reached from exactly one center
keeps it; a barcode reached from two different centers **within the same hop** is genuinely
ambiguous and stays unassigned; an earlier hop always wins over a later one.

**Distance.** `bounded_distance` first counts mismatches (equal-length strings differing in
≤ threshold positions are within threshold edits — the common substitution case, no alignment
needed), then falls back to `min(ed(a,b), ed(a[:-1],b), ed(a,b[:-1]))`. The truncated variants
absorb the spurious trailing base left when an indel inside the barcode shifts the window.

**Two implementations, identical results.**

- `full` — materialise the whole graph (`construct_graph` + `cluster`), Badger's structure.
  Kept for cross-checking.
- `centers` (default) — the graph is only ever consumed by a walk outwards from the centers,
  so the edges nobody reaches need not be found. Each round indexes just the previous round's
  frontier (the centers first, ~n_cells sequences; then whatever they claimed) and queries the
  still-unassigned barcodes against it. Round 2 therefore reaches distance 2 *only through an
  actually-observed intermediate*, which is exactly what the full graph does.

`isoquant_tests/test_barcode_graph.py::TestImplementationsAgree` asserts they match. Measured
on 90k distinct observed barcodes: **4.9× faster**, and the index never has to hold every
observed barcode at once.

### k-mer size drives clustering cost

`optimal_kmer_size(barcode_length, threshold)` picks the **largest** k for which the q-gram
bound stays usable: two sequences within `threshold` edits share at least
`L - k + 1 - k*threshold` k-mer occurrences, so the filter is lossless while that is >= 1,
i.e. while `k <= L / (threshold + 1)`. For 16 bp at threshold 1 that is **k=8**, not the 6
Badger used.

This is not a micro-optimisation. k-mers land in `4^k` buckets, so each extra base quarters
the candidates a query walks. Measured on 140,187 distinct observed barcodes: clustering went
**284 s → 31 s (9.1x)** with byte-identical assignments, and the whole correction stage
431 s → 176 s, putting it level with whitelist matching. Regression test:
`TestKmerSizeSelection::test_bound_never_rejects_a_close_pair`.

### Round 2 is the remaining bottleneck

Round 1 indexes the centers (a few thousand sequences); round 2 indexes everything round 1
claimed, which is where cost scales with the observed barcode set rather than the center set.
On the concat dataset round 2 took 35 of 63 minutes. It buys **+1.30 points of recall**
(85.25% → 86.55%) for roughly 2/3 of the runtime, and costs 0.12 of precision.

`--barcode_graph_rounds 1` is therefore a reasonable trade on very large inputs. A count gate
on the round-2 frontier (most of it is singleton error variants) would be the principled fix
and is not implemented.

### Gotcha: hit pruning must stay disabled

`get_occurrences` counts shared k-mer **occurrences**, not distinct k-mers, so for
low-complexity barcodes (long homopolymer runs) the count can far exceed
`barcode_length`. Passing `hits_delta=barcode_length` therefore does **not** disable the
top-hit filter, and the true neighbour of a poly-C barcode gets pruned — the same failure
mode as mechanism 3 above. `BarcodeGraph._no_hit_pruning()` returns
`(barcode_length - kmer_size + 1)²`, the provable maximum. Regression test:
`TestImplementationsAgree::test_low_complexity_barcodes_are_not_pruned`.

## Badger defects deliberately not ported

Kept here so nobody reintroduces them from the upstream source.

- `get_cluster_centers:314` — `cutoff = mean(list(self.counts.values())[:n_cells])` reads the
  **unsorted** dict rather than `sorted_counts`, so the cutoff reflects arbitrary barcodes
  (mostly count 1) and `max(x/5, 5)` collapses to 5 almost always. Silent quality loss.
- `get_cluster_centers:329-337` — `bc_by_counts[i]` unguarded → `IndexError` when the data has
  fewer distinct barcodes than `n_cells`; the pad loop also ignores the whitelist.
- `postprocessing:440-455` — O(|observed| × |centers|) full edit distance (2·10¹¹ calls on a
  real run). Also accepts distance < 3 while `threshold` is 1, and its `defaultdict(str)`
  lookups insert an empty string per unseen barcode.
- `compare_chunk` / `index_chunk` submitted to `ProcessPoolExecutor` as **bound methods** —
  pickles the entire `BarcodeGraph`, index included, once per 10k-barcode chunk. Our workers
  take explicit state via an `initializer` instead.
- `QGramIndex` — 4096 buckets of per-barcode dicts, ~220M boxed ints (20+ GB) for 20M observed
  barcodes. Replaced by the existing `ArrayKmerIndexer`.
- `barcodes.py:95-105` — `pd.read_csv` of the entire read table plus a parallel list of tuples:
  tens of GB on a real run. Our passes stream.
- `barcodes.py:65-71` — `--data_type` choices are `["10x", "Visium"]` but the code tests for
  `"Double"`, so `Visium` always hits the error branch.
- `get_assignments:366` / `isoquant_output:394` reference `self.barcodes`, commented out in
  `__init__` → `AttributeError`; the latter also calls `assign_by_cluster()` without `bc_len`
  and writes hardcoded `Human_R9_V5.3_*.tsv` filenames.
- `index.py:89` — `get_closest` returns mid-function, rest unreachable.
  `index.py:152 KMerIndex.add_to_index(num, barcode)` is missing `self`.
- `barcode_extraction/barcode_callers.py:164` — builds an `Array2BitKmerIndexer` over the whole
  whitelist and never queries it.
- `rank()` raises `KeyError` on `N`; extracted windows can contain `N`. We drop malformed
  windows at count time and report them in the stats.
- `Chunk`/`Edgedist` use class-level mutables, constructing a `QGramIndex(1,16,6)` at import.
- Unused heavy imports (`networkx`, `igraph`, `matplotlib`, `edlib`, `Levenshtein`, `pandas`)
  and stray `print()` calls.

## Testing

```bash
pytest isoquant_tests/test_barcode_graph.py isoquant_tests/test_correct_barcodes.py -v
pytest isoquant_tests/test_barcode_detectors.py -k RawExtraction -v
```

End-to-end:

```bash
python isoquant.py --reference ref.fa --genedb ann.gtf --complete_genedb \
  --fastq reads.fastq --data_type nanopore \
  --mode tenX_v3 --barcode_whitelist 3M-february-2018.txt \
  --barcode_correction graph --n_cells 5000 -o out -t 16
```

Accuracy against simulated data (ground truth in the read IDs):
`misc/assess_barcode_quality.py`.

CI: `.github/workflows/Barcode.Mouse.10x.GraphCorrection.yml` runs the two configs
`Mouse.10x.5k.stock_whitelist{,.graph}.yaml` under `/abga/work/andreyp/ci_isoquant/data/barcodes`
— the same reads and whitelist through both correction methods, so a regression in either
shows up as a diverging pair. Config template: `barcode_test_templates/tenX_v3_graph.cfg`.
