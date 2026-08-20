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

Every alignment is kept — primary, secondary, supplementary — plus a separate
`write_unmapped_bam` pass, because `fetch(chr)` never returns unmapped reads. Chromosomes with
no split table (not processed by this run) are still copied, just untagged, so the result stays
a complete copy of what the run covered.

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
  merged; `None` and missing fragments are skipped (chromosomes with nothing to write)
- `write_tagged_chromosome_bam` — the one copy loop, parameterised by `primary_only` and
  `keep_untagged`; feature 3 uses `(False, True)`, feature 4 uses `(True, False)`
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

`_resolve_split_molecules` now also checks `args.input_data.input_type.needs_mapping()`:
`auto` silently skips splitting for aligned input, `true` is a hard error.

## Verification performed

Unit: `isoquant_tests/test_bam_utils.py` (22 tests), `isoquant_tests/test_file_compression.py`
(18 tests).

End-to-end on chr19 of `Mouse.10x.5k.ONT_cDNA.R10.4.no_trunc.bam` (115841 records, 55392
primary), `-m tenX_v3` with the 5K whitelist:

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

### Resume

`--resume` restores `large_output` from the pickled `.params` (the resume parser rejects it on
the command line), so the survivors-file format is always consistent within a run and
`load_survivor_tags` can never meet a format the run did not write.

Resuming an **interrupted** run works. Resuming a **completed** run does not, and did not
before this branch either: `clean_up()` deletes `out_raw_file + "_*"`, which includes the
survivors files, and `prepare_read_filter` opens them without an existence guard. Verified on
the branch base (b59bfbcc), where the same scenario fails even earlier. `clean_up` and that
guard are untouched here.
