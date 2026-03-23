# Output file formats

Although most output files include headers that describe the data, a brief explanation of the output files is provided below.

## Read info (default output)

The default per-read output format. Tab-separated values, the columns are:

* `read_id` - read id;
* `chr` - chromosome id;
* `strand` - strand of the assigned isoform (not to be confused with read mapping strand);
* `gene_id` - gene id to which the read was assigned;
* `gene_assignment_type` - gene-level assignment type;
* `isoform_id` - isoform id to which the read was assigned;
* `isoform_assignment_type` - transcript-level assignment type (same values as `assignment_type` below);
* `assignment_events` - list of detected events (same format as below);
* `classification` - SQANTI-like classification: `full_splice_match`, `incomplete_splice_match`, `novel_in_catalog`, `novel_not_in_catalog`, `genic`, `antisense`, `fusion`, `intergenic`, `genic_intron`;
* `exons` - corrected exon coordinates (1-based, Illumina-corrected when available);
* `polyA` - `True` if poly-A tail is detected, `False` otherwise;
* `CAGE` - `True`/`False` if CAGE data is provided, `.` otherwise;
* `canonical` - `True` if all introns are canonical, `Unspliced` for mono-exon reads, `.` if `--check_canonical` is not set;
* `barcode` - cell barcode (`.` in bulk mode);
* `umi` - UMI sequence (`.` in bulk mode);
* `cell_type` - cell type or spot ID (`.` in bulk mode);
* `groups` - comma-separated read group IDs;
* `additional` - remaining key=value pairs (BAM tags, transcript_type, etc.).

A single read may occur more than once if assigned ambiguously.

To convert `read_info.tsv` to legacy formats, use the conversion script:
```
python -m isoquant_lib.convert_read_info --read_info SAMPLE.read_info.tsv.gz --format read_assignments --output SAMPLE.read_assignments.tsv
python -m isoquant_lib.convert_read_info --read_info SAMPLE.read_info.tsv.gz --format allinfo --output SAMPLE.allinfo.tsv
```

## Read to isoform assignment (legacy)

Tab-separated values, the columns are:

* `read_id` - read id;
* `chr` - chromosome id;
* `strand` - strand of the assigned isoform (not to be confused with read mapping strand);
* `isoform_id` - isoform id to which the read was assigned;
* `gene_id` - gene id to which the read was assigned;
* `assignment_type` - assignment type, can be:
    - `unique` - reads was unambiguously assigned to a single known isoform;
    - `unique_minor_difference` - read was assigned uniquely but has alignment artifacts;
    - `inconsistent` - read was matched with inconsistencies, closest match(es) are reported;
    - `inconsistent_non_intronic` - read was matched with inconsistencies, which do not affect intron chain (e.g. olly TSS/TES);
    - `inconsistent_ambiguous`  - read was matched with inconsistencies equally well to two or more isoforms;
    - `ambiguous` - read was assigned to multiple isoforms equally well;
    - `noninfomative` - reads is intronic or has an insignificant overlap with a known gene;
    - `intergenic` - read is intergenic.
* `assignment_events` - list of detected inconsistencies; for each assigned isoform a list of detected inconsistencies relative to the respective isoform is stored; values in each list are separated by `+` symbol, lists are separated by comma, the number of lists equals to the number of assigned isoforms; possible events are (see graphical representation below):
    - consistent events:
        - `none` / `.` / `undefined` - no special event detected;
        - `mono_exon_match` mono-exonic read matched to mono-exonic transcript;
        - `fsm` - full splice match;
        - `ism_5/3` - incomplete splice match, truncated on 5'/3' side;
        - `ism_internal` - incomplete splice match, truncated on both sides;
        - `mono_exonic` - mono-exonic read matching spliced isoform;
        - `tss_match` / `tss_match_precise` - 5' read is located less than 50 / `delta` bases from the TSS of the assigned isoform
        - `tes_match` / `tes_match_precise` - 3' read is located less than 50 / `delta` bases from the TES of the assigned isoform (can be reported without detecting polyA sites)
    - alignment artifacts:
        - `intron_shift` - intron that seems to be shifted due to misalignment (typical for Nanopores);
        - `exon_misalignment` - short exon that seems to be missed due to misalignment  (typical for Nanopores);
        - `fake_terminal_exon_5/3` - short terminal exon at 5'/3' end that looks like an alignment artifact (typical for Nanopores);  
        - `terminal_exon_misalignment_5/3` - missed reference short terminal exon;
        - `exon_elongation_5/3` - minor exon extension at 5'/3' end (not exceeding 30bp);
        - `fake_micro_intron_retention` - short annotated introns are often missed by the aligners and thus are not considered as intron retention;
    - intron retentions:
        - `intron_retention` - intron retention;
        - `unspliced_intron_retention`  - intron retention by mono-exonic read;
        - `incomplete_intron_retention_5/3` - terminal exon at 5'/3' end partially covers adjacent intron;
    - significant inconsistencies (each type end with `_known` if _all_ resulting read introns are annotated and `_novel` otherwise):
        - `major_exon_elongation_5/3` - significant exon extension at 5'/3' end (exceeding 30bp);
        - `extra_intron_5/3` - additional intron on the 5'/3' end of the isoform;
        - `extra_intron` - read contains additional intron in the middle of exon;
        - `alt_donor_site` - read contains alternative donor site;
        - `alt_acceptor_site` - read contains alternative annotated acceptor site;
        - `intron_migration` - read contains alternative annotated intron of approximately the same length as in the isoform;
        - `intron_alternation` - read contains alternative intron, which doesn't fall intro any of the categories above;
        - `mutually_exclusive_exons` - read contains different exon(s) of the same total length comparing to the isoform;
        - `exon_skipping` - read skips exon(s) comparing to the isoform;
        - `exon_merge` - read skips exon(s) comparing to the isoform, but a sequence of a similar length is attached to a neighboring exon;
        - `exon_gain` - read contains additional exon(s) comparing to the isoform;
        - `exon_detach` - read contains additional exon(s) comparing to the isoform, but a neighboring exon looses a sequnce of a similar length;
        - `terminal_exon_shift` - read has alternative terminal exon;   
        - `alternative_structure` - reads has different intron chain that does not fall into any of categories above;
    - alternative transcription start / end (reported when poly-A tails are present):
        - `alternative_polya_site` - read has alternative polyadenylation site;
        - `internal_polya_site` - poly-A tail detected but seems to be originated from A-rich intronic region;
        - `correct_polya_site` - poly-A site matches reference transcript end;
        - `aligned_polya_tail` - poly-A tail aligns to the reference;  
        - `alternative_tss` - alternative transcription start site.
* `exons` - list of coordinates for normalized read exons (1-based, indels and polyA exons are excluded);
* `additional` - field for supplementary information, which may include:
    - `gene_assignment` - Gene assignment classification; possible values are the same as for transcript classification.
    - `PolyA` - True if poly-A tail is detected;
    - `Canonical` - True if all read introns are canonical, Unspliced is used for mono-exon reads; (use `--check_canonical`);
    - `Classification` - SQANTI-like assignment classification.

Note, that a single read may occur more than once if assigned ambiguously.

## Expression table format

Non-grouped counts/TPM values are stored in a simple TSV file, the columns are:

* `feature_id` - genomic feature ID;
* `TPM` or `count` - expression value (float).

Grouped counts are stored in linear format by default - a TSV file with 3 columns:

* `feature_id` - genomic feature ID;
* `group_id` - group name;
* `count` - read count of the feature in this group. 

By default, IsoQuant converts grouped counts with small number of groups/samples (<=100) to standard matrix format; 
larger matrices (e.g. for single-cell experiments) will be saved to MTX format, which is compatible with the Seurat package.
In standard matrix rows represent features, columns represent groups. While being more human-readable, 
this file make take substantial disk space when the number of groups is large.
In MTX format, `.matrix.mtx` represents counts, `.features.tsv` contains the feature list and 
`.barcodes.tsv` contains group list (typically, barcodes are used as groups in single-cell and spatial experiments).
See [options](cmd.md#specific-output-options) to tune your output.


## UMI filtering allinfo format

**Important: this is an internal format subject to change without notice.**

Produced in single-cell/spatial modes when `--large_output allinfo` is enabled (enabled by default).
The file `SAMPLE_ID.UMI_filtered.ED{N}.allinfo[.gz]` contains one line per read that survived UMI deduplication.
Tab-separated values, the columns are:

* `read_id` - original read identifier;
* `gene_id` - assigned gene ID;
* `cell_type` - cell type or spot ID from barcode-to-spot mapping (`None` if not available);
* `barcode` - cell barcode sequence;
* `umi` - UMI sequence;
* `introns` - semicolon-separated list of intron coordinates in `chr_start_end_strand` format;
* `TSS` - transcription start site coordinate in `chr_pos_pos_strand` format (`NoTSS` if not detected);
* `polyA` - poly-A site coordinate in `chr_pos_pos_strand` format (`NoPolyA` if not detected);
* `exons` - semicolon-separated list of exon coordinates in `chr_start_end_strand` format;
* `read_type` - assignment category: `known` (unique), `known_ambiguous` (ambiguous), `novel` (inconsistent), or `none`;
* `intron_count` - number of introns in the read;
* `transcript_id` - assigned reference transcript ID;
* `transcript_type` - biotype of the assigned transcript (e.g. `protein_coding`, `lncRNA`).

Coordinate lists use `;%;` as the element separator.

An accompanying `SAMPLE_ID.UMI_filtered.ED{N}.stats.tsv` file contains summary statistics of the UMI filtering process.

## Exon and intron count format

Tab-separated values, the columns are:

* `chr` - chromosome ID;
* `start` - feature leftmost 1-based positions;
* `end` - feature rightmost 1-based positions;
* `strand` - feature strand;
* `flags` - symbolic feature flags, can contain the following characters:
    - `X` - terminal feature;
    - `I` - internal feature;
    - `T` - feature appears as both terminal and internal in different isoforms;
    - `S` - feature has similar positions to some other feature;
    - `C` - feature is contained in another feature;
    - `U` - unique feature, appears only in a single known isoform;
    - `M` - feature appears in multiple different genes.
* `gene_ids` - list if gene ids feature belong to;
* `group_id` - read group if provided (NA by default);
* `include_counts` - number of reads that include this feature;
* `exclude_counts` - number of reads that span, but do not include this feature;

## Transcript models format

Constructed transcript models are stored in usual [GTF format](https://www.ensembl.org/info/website/upload/gff.html).
Contains `exon`, `transcript` and `gene` features.

Known genes and transcripts are reposted with their reference IDs. 
Novel genes IDs have format `novel_gene_XXX_###` and novel transcript IDs are formatted as `transcript###.XXX.TYPE`,
where `###` is the unique number (not necessarily consecutive), `XXX` is the chromosome name and TYPE can be one of the following:

* nic - novel in catalog, new transcript that contains only annotated introns;
* nnic - novel not in catalog, new transcript that contains unannotated introns.

Each exon also has a unique ID stored in `exon_id` attribute.

In addition, each transcript contains `canonical` property if `--check_canonical` is set.

If `--sqanti_output` option is set, each novel transcript also has a `similar_reference_id` field containing ID of
a most similar reference isoform and `alternatives` attribute, which indicates the exact differences between
this novel transcript and the similar reference transcript.