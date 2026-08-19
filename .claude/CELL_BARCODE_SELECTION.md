# Cell barcode detection

Branch `badger`. The counting and selection here are what remain of a port of
[Badger](https://github.com/algbio/Badger) (student project, University of Helsinki; local
reference checkout `/home/andreyp/Badger`). Its edit-distance graph correction was
implemented, benchmarked against the existing SSW matcher, and then dropped — see
"What was dropped and why".

## Why

Matching every read against a whitelist works when the whitelist is the list of cells. With
the stock 10x whitelists (737k for v2, ~3M for v3) it **collapses into exact matching**:

| # | Where | Mechanism |
|---|---|---|
| 1 | `callers/tenx.py` (`bc_count > 100000` branch) | `min_score = 16 == BARCODE_LEN_10X` with SSW `match_score=1, mismatch_penalty=1` — only a perfect 16 bp match reaches that score |
| 2 | `common.py find_optimal_kmer_size` | k grows with whitelist size; 16 bp × 3M barcodes gives k≈12, leaving 5 k-mers per barcode, and `min_matching_kmers=2`. An error in the middle destroys all of them, so the true barcode is never even a candidate |
| 3 | all indexers, `hits_delta=1` | candidates more than one shared k-mer behind the best are dropped before alignment |

A run has a few thousand cells, not millions. Counting the extracted barcodes first turns the
whitelist into a **filter over a few thousand candidates** instead of a per-read search space,
and then ordinary SSW matching against that short list does the correction.

## Design

`--barcode_whitelist` is mandatory for single-cell modes. `--n_cells` decides what it *means*:
unset, it is the set of cell barcodes; set, it is a pool to select cell barcodes from.

| `--barcode_whitelist` | `--n_cells` | behaviour | read passes |
|---|---|---|---|
| file | unset | the file **is** the cell set → match against it | 1 |
| file | `N` or `auto` | file is a candidate pool → detect cells → match against those | 2 |
| `auto` | unset → implies `auto` | detect cells from counts alone → match against those | 2 |
| `auto` | `N` or `auto` | detect N cells from counts → match against those | 2 |

`--barcode_whitelist auto` cannot work without a cell count, so it implies `--n_cells auto`.
Correction is **always** the existing SSW matcher; only the list it matches against varies.
`--barcode_correction` is retained as a hidden override (visible under `--full_help`).

Supplying a large whitelist with `--n_cells` unset warns, because that is the degenerate case
above.

## Flow

```
                     ┌─ --n_cells unset ─────────────────────────┐
reads ───────────────┤                                            ├─> barcoded_reads_<i>.tsv
                     └─ --n_cells set ──> pass 1: extract and count
                                          (in memory, nothing written)
                                                  │
                                          select cell barcodes
                                          <prefix>.cell_barcodes.tsv   (+ .stats)
                                                  │
                                          pass 2: ordinary barcode calling
                                                  against the detected cells
```

Pass 2 is the unmodified detector, so `barcoded_reads_<i>.tsv` keeps its exact format and
everything downstream (`dataset_processor.split_read_barcode_table`) is untouched. In
splitting modes only pass 2 writes the FASTA — pass 1 is purely for counting.

Pass 1 writes nothing. Workers tally barcodes into a `CellBarcodeSelector` and return count
dicts, which the parent merges with `merge_counts`; only the counts are needed, and a barcode
per read would cost ~126 bytes/read of intermediate — 12 GB on a 100M-read run — written and
re-read for nothing. Counting and selection run as one child process
(`detect_cell_barcode_list`) so the count table never reaches the main process.

Resume markers: `aux/<prefix>.raw_barcodes_done` for pass 1, the existing
`aux/<prefix>.barcodes_done_<i>.tsv` for pass 2, so a re-run reuses a detected cell list
without repeating the extraction.

## Code map

| File | Contents |
|---|---|
| `isoquant_lib/barcode_calling/cell_selection.py` | `CellBarcodeSelector`, `estimate_cell_number`, `select_cell_barcodes` |
| `detect_barcodes.py` | `count_barcodes_in_reads`, `count_chunk`, `detect_cell_barcode_list`, `run_chunks_in_parallel` |
| `callers/tenx.py` | `whitelist_matching` flag, `_extract_raw_barcode`, `_raw_strand_score` |
| `isoquant.py` | `_resolve_barcode_correction`, `detect_cell_barcodes`, `_run_barcode_calling` |
| `modes.py` | `AUTO_BARCODES`, `BarcodeCorrectionMethod`, `IsoQuantMode.supports_cell_barcode_detection()` |

Supported modes: `tenX_v2`, `tenX_v3`, `visium_5prime` (with or without `--split_molecules`).
Stereo-seq is excluded — every spot barcode there is genuine, so count-based selection does
not apply. `custom_sc` is excluded because it is MDF-driven; extending it looks worthwhile.

## Benchmarks

recall / precision, measured with `misc/assess_barcode_quality.py`. "stock" is the 7.4M
`3M-3pgex-may-2023` list with the dataset's true cell barcodes prepended.

| dataset | cells | supplied cell list | stock list only | stock + `--n_cells N` | stock + `--n_cells auto` | `-b auto` |
|---|---|---|---|---|---|---|
| Mouse 10x ONT StereoQ, 600k | 5000 | 87.33 / 99.75 | 74.36 / 99.60 | **87.33 / 99.75** | **87.33 / 99.75** | 87.33 / 99.75 |
| Mouse 10x ONT cDNA R10.4, 304k | 5000 | 89.61 / 99.83 | 76.27 / 99.58 | **89.61 / 99.83** | **89.61 / 99.83** | 85.61 / 94.82 |
| Mouse 10x concat, 800k (split) | 8975 | 86.78 / 99.81 | 74.31 / 99.63 | **86.78 / 99.81** | **86.78 / 99.81** | 86.72 / 99.34 |

Two results matter:

1. **With a whitelist as candidate pool the outcome is identical to supplying the true cell
   list** — on all three datasets, to two decimal places. That is the whole point: the user no
   longer has to know which barcodes are the cells. Against taking the stock list at face
   value this is +13 points of recall and +0.2 precision.
2. **Cell number accuracy barely matters.** `--n_cells auto` equals an exact `--n_cells`, and
   selection tolerates being told 4000 or 6000 for a true 5000 (the count cutoff does the
   work; `n_cells` only bounds the search). Estimates: 5315 for a true 5000, 10389 for 8975.

Cost: two passes, ~2.3x the single-pass time (100s → 240s on StereoQ, 667s → 2101s on concat).
Pass 1 writes nothing, so the two passes cost read time but no extra disk.

### `-b auto` is measurably worse — prefer a pool

With no whitelist at all, R10.4 drops to **85.61 / 94.82** (−4.0 recall, −5.0 precision). It
detected 5008 barcodes: all 5000 true cells, plus 8 spurious ones. Those 8 are **not** error
variants — they sit 4–6 edits from any real cell barcode. They are recurring extraction
artifacts (a mis-anchored barcode window repeats across reads and looks exactly like an
abundant cell), and counts alone cannot distinguish them from cells. A whitelist rejects them
because an artifact is not a valid protocol barcode. The concat dataset shows the same effect
more mildly (9031 detected, precision 99.81 → 99.34).

`select_cell_barcodes` warns when it runs without a pool.

## What was dropped and why

Badger corrects barcodes by building an edit-distance graph over all observed barcodes and
walking outwards from the selected cells. That was implemented, and then measured against the
existing SSW matcher **given the same 5000 barcodes**:

| | recall | precision |
|---|---|---|
| SSW matcher | **87.35%** | 99.73% |
| graph correction | 86.55% | 99.84% |

They disagree on 3.6% of reads. Of the 13,105 reads only SSW gets right, 8,373 (64%) have an
identical extracted window — a local alignment slides to find a barcode that shifted inside
its window, a fixed-window edit distance cannot. The rest come from whitelist-driven strand
selection (3,555) and stricter R1 anchoring (1,177), neither of which is available without a
whitelist during extraction.

Graph multi-hop (reaching a barcode two edits away through an *observed* intermediate) does
recover 8,322 reads SSW misses, but inspection showed these are mostly frame-shifted
windows — exactly what SSW handles natively. The union ceiling is 88.7%; keeping both
mechanisms was judged not worth a second correction path and the round-2 bottleneck (2/3 of
runtime on the concat dataset).

Defects found in the Badger source and deliberately not carried over, kept here so nobody
reintroduces them from upstream:

- `get_cluster_centers` computed its count cutoff from the **unsorted** count dict, so the
  cutoff reflected arbitrary barcodes and `max(x/5, 5)` collapsed to 5 almost always.
- `bc_by_counts[i]` unguarded → `IndexError` when the data has fewer distinct barcodes than
  `n_cells`; the pad loop also ignored the whitelist.
- `postprocessing` was O(|observed| × |centers|) full edit distance — 2·10¹¹ calls on a real
  run.
- `compare_chunk` / `index_chunk` were submitted to `ProcessPoolExecutor` as bound methods,
  pickling the entire graph once per 10k-barcode chunk.
- `pd.read_csv` of the whole read table — tens of GB on a real run. Our passes stream.
- `--data_type` choices were `["10x", "Visium"]` while the code tested for `"Double"`, so
  `Visium` always hit the error branch.
- `rank()` raised `KeyError` on `N`. Malformed windows are dropped at count time and reported.

## Gotchas

### Strand selection dominates raw-mode accuracy

Whitelist mode picks the orientation on which the barcode matched — strong evidence. Raw
extraction has no such signal: every window is accepted and `BC_score` is always 0, leaving
`more_informative_than` comparing polyT positions, i.e. deciding close to at random. A wrong
strand yields a reverse-complemented window nothing can rescue.

Measured, not guessed: of 102,696 reads that whitelist matching got right and the raw path did
not, **102,590 were wrong-strand picks**, their windows sitting 7–11 edits from truth.
`_raw_strand_score` / `_pick_raw_strand` rank an orientation by molecule structure instead —
R1 and polyT both present, `0 < polyT - r1 <= MAX_R1_POLYT_DISTANCE`, span closest to
`BARCODE_LEN_10X + UMI_LEN`. Recall went 66.49% → 86.55%. Note the splitting detectors do
**not** use this: they scan both strands and keep molecules from each, with
`_is_consistent_detection` doing the structural filtering.

### `score_diff` was broken

`find_candidate_with_max_score_ssw` is supposed to reject a read when the winning barcode has
no margin over the runner-up. It never did:

- `second_best_score` was only assigned in the branches that *replaced* the best, so a genuine
  lower-scoring runner-up never updated it; it stayed `0` and `best - 0 >= score_diff` always
  passed.
- The tie branch was guarded by `alignment.reference_start < best_match[1]`, so two candidates
  tying at the same offset — the normal case — never registered.
- It was `0` on the ≤100k branch, which had no ambiguity check at all.

Now both are tracked and `score_diff=1` on both 10x branches. Curio and Stereo call this
function with the default `score_diff=0` and are unaffected by design.

**The effect scales with how many ties the whitelist creates**, so it is negligible for 10x
and large for the universal extractor (`custom_sc`), which sets `SCORE_DIFF = 1` and had the
broken check. Measured across the full barcode CI suite:

| test | whitelist | recall | precision | errors removed |
|---|---|---|---|---|
| tenX_v3 StereoQ | 5k cells | 87.35 → 87.33 | 99.73 → 99.75 | — |
| tenX_v3 R10.4 | 5k cells | 89.63 → 89.61 | 99.80 → 99.83 | — |
| 10x concat (split) | 8975 | 86.80 → 86.78 | 99.80 → 99.81 | 2286 → 2127 |
| custom_sc 10x | 5k | 85.48 → 84.85 | 99.73 → 99.76 | 1407 → 1231 |
| custom_sc VisiumHD large | 644k | 3.64 → 3.49 | 96.79 → **98.95** | 966 → 297 (−69%) |
| custom_sc VisiumHD small | 644k | 60.02 → **54.93** | 90.88 → **98.95** | 48199 → 4679 (−90%) |

The VisiumHD `custom_sc` cases trade real recall for large precision gains. That was accepted
deliberately: a misassigned barcode puts a read in the wrong cell and silently corrupts that
cell's counts, whereas a dropped read only costs depth. Baselines for the five affected
configs were refreshed from measurement. Note the native `visium_hd` mode is unaffected — it
uses `find_candidate_with_max_score_ssw_var_len`, a separate function.

### Cell number estimation on a flat distribution

`estimate_cell_number` finds the knee as the point furthest below the chord of the log-log
count-rank curve. When the curve never bends, every chord distance is 0 and `argmax` returns
index 0 — silently reporting a single cell. It now falls back to "everything above the noise
floor". Regression test: `test_flat_distribution_has_no_knee`.

## Testing

```bash
pytest isoquant_tests/test_cell_selection.py isoquant_tests/test_barcode_common.py -v
pytest isoquant_tests/test_barcode_detectors.py -k RawExtraction -v
```

End-to-end:

```bash
python isoquant.py --reference ref.fa --genedb ann.gtf --complete_genedb \
  --fastq reads.fastq --data_type nanopore \
  --mode tenX_v3 --barcode_whitelist 3M-february-2018.txt --n_cells auto \
  -o out -t 16
```

CI: `.github/workflows/Barcode.Mouse.10x.CellSelection.yml` runs the same reads and whitelist
through both the supplied-list and detected-cells paths, so a regression in either shows up as
a diverging pair. Config template: `barcode_test_templates/tenX_v3_cell_selection.cfg`.
