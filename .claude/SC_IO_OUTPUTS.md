# Single-cell I/O outputs

Branch `sc_outputs` (based on `badger`, PR #420). Four independent changes to what the
single-cell pipeline writes to disk: two reclaim space on existing outputs, two add BAM
outputs carrying information the pipeline already computes but previously only emitted as text.

## 1. The split-reads FASTA is gzipped

`<prefix>.split_reads_<i>.fa.gz` unless `--no_gzip` (dest `args.gzipped`, default `True`).
The name is built in `call_barcodes` (`isoquant.py`), which picks the suffix from `args.gzipped`.

Compression happens **inside the chunk workers**, not serially in the parent:
`detect_barcodes.py` writes each per-chunk temp through `open_text_write`, and the merge that
follows is a plain byte copy. That works because concatenated gzip members form a valid gzip
stream — the merge does not decompress and recompress. `numbered_chunk_name()` keeps `.gz`
last so the temps are recognised as compressed.

Nothing downstream slows down: the only consumer is minimap2, which gets a path on the command
line and decompresses natively. No Python code opens this file.

Two things the extension change touched, both fixed:

- `string_pools.py` and `read_groups.py` derive file-name read groups with a single
  `os.path.splitext`, so `.fa.gz` would have left a trailing `.fa` in the group name. Both now
  call `strip_compression_suffix()` first. `test_file_compression.py` pins this.
- `read_mapper.py` derives the BAM name and already handled `.fa.gz` correctly (verified, not
  changed).

## 2. The barcoded read tables are gzipped at the end of the run

`<prefix>.barcoded_reads_<i>.tsv` (`sample.barcodes_tsv`) stays **plain for the whole run** —
`split_read_barcode_table` reads it, and so does the tagged-BAM build. `compress_barcode_tables`
runs in `process_all_samples` after every sample is done and before `clean_up()`.

The per-worker split tables under `aux/` (`sample.barcodes_split_reads + "_<chr>"`) are
deliberately **left uncompressed**: short-lived temporaries read back by Python, so compressing
them would cost CPU in every chromosome worker for a transient saving.

Resume hazard, handled: `call_barcodes` short-circuits on the `barcodes_done` markers and
re-populates `sample.barcoded_reads` with the *uncompressed* names, which no longer exist after
compression. Readers go through `resolve_optionally_gzipped()`, which returns whichever of
`<path>` / `<path>.gz` exists.

Helpers live in `isoquant_lib/utils/file_utils.py`: `open_text_write`, `open_text_read`,
`resolve_optionally_gzipped`, `gzip_file_in_place`, `strip_compression_suffix`.

Both writers pass `compresslevel=GZIP_LEVEL` (6). Python's `gzip` defaults to 9. Measured on a
real 25 MB barcode table from 411K ONT reads (UUID read ids, so realistic entropy — an earlier
synthetic table was near-incompressible and badly overstated the gap):

| level | barcode table | vs L6 size | split FASTA | vs L6 size |
|-------|---------------|-----------|-------------|-----------|
| 1 | 94.6 MB/s | +8.5% | 92.5 MB/s | +13.8% |
| 4 | 55.8 MB/s | +2.0% | 65.3 MB/s | +5.7% |
| 5 | 44.2 MB/s | +0.4% | 32.5 MB/s | +3.4% |
| 6 | 39.2 MB/s | — | 12.5 MB/s | — |
| 9 | 16.0 MB/s | −3.0% | 2.7 MB/s | −3.9% |

So the two kinds of output take different levels, picked by `gzip_level_for(name)` from the
file's own extension:

- **tables** (TSV, BED, MTX, allinfo) — `GZIP_LEVEL = 6`. Going below is not worth it:
  `compress_barcode_tables` is the one serial compressor, but at 39 MB/s a ~100 GB table
  (roughly a billion reads) costs ~45 min against ~105 min at level 9, and 6→1 would save
  another ~25 min for 8.5% more disk while 4-5 save almost nothing.
- **sequences** (FASTA/FASTQ) — `GZIP_LEVEL_SEQUENCES = 4`. Nucleotide data sits near gzip's
  entropy floor so the high levels grind for nothing: 5.2x the throughput for 5.7% more output.

Every gzip *writer* in the project now goes through `open_text_write` / `gzip_file_in_place`,
so both levels apply everywhere: the allinfo writers in `dataset_processor`, the four in
`convert_grouped_counts`, `TextFileAssignmentPrinter` (read_info / read_assignments /
corrected_bed / read2transcripts), and the two `scripts/` converters. None of them call
`gzip.open` any more.

### Barcode calling output order

`run_chunks_in_parallel` waits on `FIRST_COMPLETED` and handed results to `handle_result` in
**completion order**, and `_process_single_file_in_parallel` appended the chunk temp names to a
list in that order — so which chunk finished first decided the row order of
`barcoded_reads.tsv` and the record order of the split FASTA. Two runs of identical code on
identical input differed (measured: at line 200003, a chunk boundary). Row sets were always the
same and every consumer keys by read id, so nothing was wrong; the outputs simply were not
reproducible.

`handle_result` now takes `(result, chunk_index)`. The merging caller keys a dict on the index
and sorts once at merge time; the counting caller ignores it, being order-insensitive already
(`CellBarcodeSelector.sorted_barcodes` sorts by `(count, barcode)`, a total order). Scheduling
is untouched — waiting for chunks *in order* instead would idle the pool behind one slow chunk,
which is the property `run_chunks_in_parallel` exists to provide.

Verified: two runs now produce a byte-identical barcode table and a byte-identical decompressed
FASTA, with the same row and read-id sets as before the change. The compressed FASTA still
differs in exactly 10 bytes, all of them gzip header `mtime` fields — the deflate streams are
identical. Passing `mtime=0` would close that too, if byte-comparable `.gz` outputs are ever
wanted.

Only barcode calling had this pattern; every other parallel stage uses `proc.map`, which
preserves input order.

**Gotcha, caught only by measuring the output size.** The split FASTA is compressed in the
chunk workers, and those temps were called `subreads.gz` — no `.fa`, so the inference gave them
the *table* level while the single-threaded path, which passes the real output name, gave them
the sequence level. The level differed by thread count. `numbered_chunk_name` now keeps the
whole extension chain last (`subreads.fa.gz` → `subreads_3.fa.gz`) and the temp is named for
what it holds. Verified end to end: the FASTA went 5.28 MB → 5.58 MB, exactly the level-4
number, with identical read-id sets.

The table has to stay plain during the run for a second reason worth recording:
`split_read_table_parallel` streams it line by line in *every* worker, so a gzipped table would
be decompressed once per thread rather than read once.

## 3. `--large_output tagged_bam`

`<prefix>.tagged.bam` — a copy of the input alignments with `--barcode_tag` / `--umi_tag`
(CB/UB) added. Off by default.

Named `tagged_bam` rather than `barcoded_bam` because `--barcoded_bam` is already an *input*
flag meaning the opposite direction.

**Purely a side output.** The barcode-table split happens regardless, and nothing downstream
reads the result — building this BAM needs those split tables in the first place, so there was
never a split to save by reusing its tags.

Built in `DatasetProcessor.write_tagged_bam`, right after the split-table block in
`process_sample` while the tables still exist. One fragment per chromosome via
`map_over_chromosomes(write_tagged_bam_in_parallel, ...)`, then merged and indexed.

Guarded by its own resume marker, `tagged_bam_lock_filename(sample)` (in `aux/`, so it outlives
`clean_up`) plus an existence check on the BAM itself. This output is a full copy of the input,
the most expensive thing on the branch, and it is written *before* read collection — without
the marker every `--resume` after a crash in the long stages re-copied the whole BAM. The
marker follows the `barcodes_done` precedent and is deliberately not deleted at the end of
`process_sample`.

The reference list is computed **once** in `process_sample` and passed into `write_tagged_bam`.
It has to be the same list that drove the barcode-table split, or a fragment finds no table and
its reads come out untagged; computing it twice left that invariant implicit.

Every alignment is kept — primary, secondary, supplementary — plus a separate
`write_unmapped_bam` pass, because `fetch(chr)` never returns unmapped reads.

The chromosome list comes from `references_with_alignments(bam_files)`, **not** from
`get_chr_list()`. This was a bug found only on the full CI dataset: IsoQuant analyses the
22 assembled mouse chromosomes, so iterating the analysed list silently dropped 1177 records
sitting on unplaced scaffolds (`GL456382.1`, `JH584299.1`, …) — in a file documented as a copy
of the input. The chr19 subset used during development contained no scaffolds and could not
surface it. Empty references are filtered out via the BAM index so a fragmented assembly does
not spawn a task per contig.

Restoring the records is only half of it: they came back *untagged*, because the
per-chromosome barcode split is driven by `split_barcodes_dict`, which was keyed on the
analysed chromosomes too. Barcode calling runs over the whole input before any chromosome
filtering, so the tags exist in `<prefix>.barcoded_reads_<i>.tsv` — they simply never reached a
split file the tagging worker could read. `process_sample` now widens that dict to
`references_with_alignments(...)` when `tagged_bam` is enabled (and only then, since it costs
an extra pass over the barcode table).

### Unmapped reads were never tagged

A separate bug, on the same output, found the same way. `write_unmapped_bam` copied unplaced
records verbatim with **no tags at all**, so on the CI dataset 3193 of the 3577 unmapped reads
lost barcodes they genuinely had. A barcode is called from the read sequence and does not
depend on the read aligning anywhere, so those tags are meaningful — arguably more so, since
recovering unmapped reads per cell is a reason to want this file.

They belong to no chromosome and so appear in no split table. `write_tagged_bam` therefore
collects the unmapped read ids first (`collect_unmapped_read_ids`, a few thousand) and scans
the whole-sample barcode table for just those (`load_barcode_umi_tags(..., read_ids=...)`),
which keeps memory at the size of the unmapped set rather than the 1.8M-row table.

Worth recording how this was nearly missed: after the scaffold fix the tag-mismatch count
stayed at exactly 3193, and the first reading was that the scaffolds explained it — the
arithmetic `1794571 - 1791378 = 3193` matched. It was a coincidence of two unrelated numbers.
Checking which read ids were actually untagged showed all 3193 were unmapped, and that the
scaffold records had contributed no new barcoded read ids at all.

With `--barcoded_bam` as input there are no split tables and the input already carries the
tags, so IsoQuant warns and skips; likewise in modes with no barcodes at all.

## 4. `--large_output deduplicated_bam`

`<prefix>.deduplicated.bam` — **primary alignments only**, restricted to the reads that survived
UMI filtering, tagged with barcode, UMI, `GX` (gene) and `TX` (transcript). Off by default.

All four values are in hand exactly where survivors are chosen, so no extra pass over the
assignments is needed: `UMIFilter._process_chunk` already holds the `ReadAssignment`.
`UMIFilter._survivor_record` writes them as extra tab-separated columns on the existing
per-chromosome survivors file (`<aux>/<SAMPLE>.save_filtered_<chr>`) when `output_read_tags`
is set — which `parallel_workers.filter_umis_in_parallel` derives from
`large_output_enabled(args, "deduplicated_bam")`.

Its one other consumer, `prepare_read_filter` (`assignment_loader.py`), now takes
`line.split("\t")[0]`, which makes it tolerant of both the bare-id and the tagged format.

Two ordering constraints, both respected in `process_sample`:

- built **after** `filter_umis` and **before** `clean_up`, which deletes `out_raw_file + "_*"`
  including the survivors files;
- only the first edit distance writes those files, and the `barcode2barcode` rounds never do,
  so the subset is defined by the primary dedup round.

### Deliberately *not* wired into fusion detection

An earlier version of this branch auto-enabled `deduplicated_bam` for fusion runs and fed it to
`FusionDetector` in place of the original BAM, reasoning that PCR duplicates would otherwise
inflate breakpoint support. **That was wrong and has been reverted**: fusion evidence lives in precisely the reads UMI filtering removes. The filter keeps
one read per assigned (gene, barcode, UMI) molecule and requires a gene assignment, so chimeric
reads spanning two genes — the inconsistent reads fusion calling is built on — are collapsed or
dropped outright.

(The secondary/supplementary dimension was in fact fine: `fusion_detector.py:470` skips both and
takes breakpoints from the SA tag on the primary record. That is not what makes the subset
unsuitable — the read-level filtering is.)

Deduplicating for fusion, if wanted at all, belongs **inside** the fusion algorithm where it can
see the chimeric reads before they are filtered. Fusion detection therefore reads the original
BAMs, exactly as it did before this branch — `get_bam_files_from_samples` is untouched.

## Shared code

`isoquant_lib/utils/bam_utils.py`:

- `index_bam` — BAI with a CSI fallback for references too long for BAI
- `merge_bam_files` — merge per-chromosome fragments and index; single fragment is moved, not
  merged; `None` and missing fragments are skipped (chromosomes with nothing to write).
  Merging goes through `_merge_in_rounds` in batches of `BAM_MERGE_BATCH` (500): samtools opens
  every input at once and gives up just above a thousand handles regardless of `ulimit -n`
  (measured: 1200 fragments fail at fragment 1019 with `RLIMIT_NOFILE` at 1048576). One
  fragment per non-empty reference means GRCh38's full analysis set, at 3366 contigs, would
  have crashed the merge. A lone leftover batch is carried into the next round rather than
  copied; the intermediates are deleted once the final merge succeeds.
- `unplaced_reads` — the unmapped records with no coordinates. `fetch("*")` seeks straight to
  them on an indexed BAM, with a `fetch(until_eof=True)` filter as the fallback. The two
  callers (`collect_unmapped_read_ids` and `write_unmapped_bam`) previously each scanned the
  entire BAM linearly just to reach the tail, both serially in the parent
- `write_tagged_chromosome_bam` — the one copy loop, parameterised by `primary_only` and
  `keep_untagged`; feature 3 uses `(False, True)`, feature 4 uses `(True, False)`
- `references_with_alignments` — every reference carrying reads, empty ones dropped via the
  index; what the tagged BAM iterates instead of the analysed chromosome list
- `write_unmapped_bam`
- `load_survivor_tags` / `load_barcode_umi_tags` — the two tag sources
- `PLACEHOLDERS = {"*", ".", "None"}` — a tag is **omitted** rather than carrying a
  placeholder. `"None"` is in there because IsoQuant stores the literal string `"None"` as
  `assigned_transcript` for novel and ambiguous reads (visible in `allinfo` too), and
  `TX:Z:None` would be a trap for anything reading the tag.

`DatasetProcessor.map_over_chromosomes(worker, sample, *extra)` factors out the
`ProcessPoolExecutor` boilerplate the per-chromosome stages all repeat.

## Regression fixed along the way

`--split_molecules auto` (default since the `badger` branch) made every splitting-capable mode
split, including when input was an **aligned BAM**. Splitting rewrites the reads, so the pieces
need re-aligning, but `--bam` input skips the mapping stage — `call_barcodes` replaced
`sample.file_list` with the FASTA and `get_chromosome_ids` then tried to open it as a BAM
(`ValueError: file has no sequences defined`). This broke `-m tenX_v3 --bam ...`, which the
`SC.Mouse.10x.allinfo` CI config uses.

Fixed by `_reject_splitting_aligned_input` in `isoquant.py`, run right after
`resolve_split_molecules`: splitting with aligned input aborts, under `auto` as well as `true`.
The check lives in the pipeline rather than in `options.resolve_split_molecules` because that
module is barcode-calling helpers, shared with the standalone `isoquant_detect_barcodes.py`,
which has no `input_data` at all.

Aborting on `auto` rather than quietly not splitting is deliberate: passing a BAM says "do not
map", asking for splitting says "rewrite the reads", and either way of guessing hands the user
something they did not ask for. Four CI configs (`SC.Mouse.10x.allinfo`,
`SC.Mouse.10x.barcoded_bam.allinfo`, `GROUP12/13.SC.SIRVs.R10`) pair `--bam` with a splitting
mode and now carry an explicit `--split_molecules false`; all four were already broken by the
crash before this.

## Verification performed

Unit: `isoquant_tests/test_bam_utils.py` (35 tests), `isoquant_tests/test_file_compression.py`
(18 tests). The merge test monkeypatches `BAM_MERGE_BATCH` to 3 and runs 11/12/13 fragments,
straddling the batch boundary where a lone leftover has to be carried forward; a separate test
pins that `unplaced_reads` returns the same records with and without an index.

End-to-end on chr19 of `Mouse.10x.5k.ONT_cDNA.R10.4.no_trunc.bam` (115841 records, 55392
primary), `-m tenX_v3` with the 5K whitelist — and then on the **full** CI dataset via
`SC.Mouse.10x.allinfo`, which is what caught the scaffold bug above:

- `deduplicated_bam`: 25056 records, 0 secondary/supplementary; read-id set is **exactly** the
  25056 rows of `allinfo`; CB/UB/GX all match `allinfo` on every row; records byte-identical to
  the input apart from tags; SA tags preserved.
- `tagged_bam`: 115841 records — identical to the input count, secondary and supplementary
  included; CB/UB match the barcode table exactly on all 49756 barcoded reads.
- Both are **pure side outputs**: every other file in the run is byte-identical to a run
  without the flag, apart from two `.gz` files whose embedded mtime differs (decompressed
  content identical).
- Fusion: verified untouched — `--analysis fusion` neither enables nor reads either BAM, and
  the fusion code path is byte-identical to the branch base.

On the full CI dataset (2M reads, 3660591 records) via `SC.Mouse.10x.allinfo`:

- `deduplicated_bam`: 903570 records, all primary, exactly matching that run's own
  `Total reads saved`.
- `tagged_bam`, after all three fixes: 3660591 records (exact), 1655717 secondary and 4874
  supplementary preserved, 3577 unmapped preserved of which 3193 barcoded, and CB/UB matching
  the barcode table on **all 1794571** barcoded reads with 0 differences.
- Tag *values* were never wrong at any stage — 0 conflicting and 0 extra tags throughout; every
  defect was a read or a tag going missing, never a wrong one.
- The `allinfo` baselines pass within tolerance on every run.

### CI coverage

`SC.Mouse.10x.allinfo` requests `tagged_bam deduplicated_bam`; `SC.Mouse.10x.barcoded_bam.allinfo`
requests `deduplicated_bam` only — `tagged_bam` is deliberately omitted there because
`--barcoded_bam` skips the split-table block it reads from, so it would warn and skip. That
second config is worth having because it exercises the dedup BAM on the tags-read-from-BAM
path, where the tags never pass through a barcode table at all.

Note that CI *produces* these BAMs but does not assert anything about them — the baselines only
cover `allinfo`. A crash or a knock-on regression would be caught; a wrong tag value would not.

### File-name read groups

`strip_compression_suffix` was added for the `.fa.gz` split reads, but it changes
`--read_group file_name` for **any** gzipped input: `reads.fq.gz` used to group as `reads.fq`
(one `splitext` off a two-part extension) and now groups as `reads`. That is the intended
name, and `FileNameGrouper` and `StringPoolManager` were changed together so they still agree,
but it is a visible change to group labels for existing bulk runs. No CI baseline moves:
`STEREO.TOY` is the only config pairing `file_name` with a gzipped input and it is
`run_type: void` (checks that the grouped count files exist, not their contents).

### Resume

`--resume` restores `large_output` from the pickled `.params` (the resume parser rejects it on
the command line), so the survivors-file format is always consistent within a run and
`load_survivor_tags` can never meet a format the run did not write.

Resuming an **interrupted** run works, and the tagged BAM is now skipped rather than rebuilt
(verified: the skip is logged and the file's mtime does not move).

Resuming a **completed** run does not, and did not
before this branch either: `clean_up()` deletes `out_raw_file + "_*"`, which includes the
survivors files, and `prepare_read_filter` opens them without an existence guard. Verified on
the branch base (b59bfbcc), where the same scenario fails even earlier. `clean_up` and that
guard are untouched here.

The mechanism, for whoever fixes it: `clean_up` removes the *global* UMI lock but not the
per-edit-distance one (`umi_filtered_lock_file_name`, in `aux/`), so a resumed `filter_umis`
returns early without rewriting the survivors files it just deleted. `write_deduplicated_bam`
then finds nothing and warns "No reads survived UMI filtering", which is a misdiagnosis — but
the run dies seconds later in `prepare_read_filter` for the same underlying reason, so nothing
was added here to paper over it.
