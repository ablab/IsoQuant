# IsoQuant changelog

## IsoQuant 3.4.1, 9 May 2024

- Fixes `IndexError: list index out of range` when `--sqanti_output` is set ([#186](https://github.com/ablab/IsoQuant/issues/186)).

- Fixes `IndexError: list index out of range` in printing grouped transcript models TPMs ([#187](https://github.com/ablab/IsoQuant/issues/187)).

- Reduced running time when `--sqanti_output` is set.

## IsoQuant 3.4.0, 9 May 2024

Major novelties and improvements:

- Significant speed-up on datasets containing regions with extremely high coverage,
    often encountered on mitochondrial chromosomes ([#97](https://github.com/ablab/IsoQuant/issues/97)).

- Added support for Illumina reads for spliced alignment correction (thanks to [@rkpfeil](https://github.com/rkpfeil)).

- Added support YAML files (thanks to [@rkpfeil](https://github.com/rkpfeil)). Old options `--bam_list` and `--fastq_list` are still availble, but deprecated since this version.

Transcript discovery and GTF processing:

- Fixed missing genes in extended GTF ([#140](https://github.com/ablab/IsoQuant/issues/140),
    [#147](https://github.com/ablab/IsoQuant/issues/147), [#151](https://github.com/ablab/IsoQuant/issues/151),
    [#175](https://github.com/ablab/IsoQuant/issues/175)).

- Fixed strand detection and output of transcripts with `.` strand ([#107](https://github.com/ablab/IsoQuant/issues/107)).

- Added `--report_canonical` and `--polya_requirement` options that allows to control level of
    filtering of output transcripts based on canonical splice sites and the presence of poly-A tails. ([#128](https://github.com/ablab/IsoQuant/issues/128))

- Added check for input GTFs ([#155](https://github.com/ablab/IsoQuant/issues/155)).

- Extract CDS, other features and attributes from reference GTF to the output GTFs ([#176](https://github.com/ablab/IsoQuant/issues/176)).

- Reworked novel gene merging procedure ([#164](https://github.com/ablab/IsoQuant/issues/164)).

- Revamped algorithm for assigning reads to novel transcripts and their quantification ([#127](https://github.com/ablab/IsoQuant/issues/127)).

Read assignment and quantification:

- Optimized read-to-isoform assignment algorithm.

- Added `gene_assignment_type` attribute to read assignments.

- Fixed duplicated records in `read_assignments.tsv` ([#168](https://github.com/ablab/IsoQuant/issues/168)).

- Improved gene and transcript quantification. Only unique assignments are now used for transcript quantification.
Added more options for quantification strategies (`--gene_quantification` and `--transcript_quantification`). 

- Improved consistency between `trascript_counts.tsv` and `transcript_model_counts.tsv` ([#137](https://github.com/ablab/IsoQuant/issues/137)).

- Introduced mapping quality filtering: `--min_mapq`, `--inconsistent_mapq_cutoff` and
    `--simple_alignments_mapq_cutoff` ([#110](https://github.com/ablab/IsoQuant/issues/110)).

Minor fixes and improvements:

- Added `--bam_tags` option to import additional information from BAM files to read assignments output.

- Large output files are now gzipped by default, `--no_gzip` can be used to keep uncompressed output
    ([#154](https://github.com/ablab/IsoQuant/issues/154)).

- BAM stats are now printed to the log ([#139](https://github.com/ablab/IsoQuant/issues/139)).

- Various minor fixes and requests ([#106](https://github.com/ablab/IsoQuant/issues/106),
    [#141](https://github.com/ablab/IsoQuant/issues/141), [#143](https://github.com/ablab/IsoQuant/issues/143),
    [#146](https://github.com/ablab/IsoQuant/issues/146), [#179](https://github.com/ablab/IsoQuant/issues/179)).


Special acknowledgement to [@almiheenko](https://github.com/almiheenko) for testing and reviewing PRs, 
and to [@alexandrutomescu](https://github.com/alexandrutomescu) for supporting the project. 


## IsoQuant 3.3.1, 26 July 2023

- Fixed `UnboundLocalError: local variable 'match' referenced before assignment` error in SQANTI-like output.

## IsoQuant 3.3, 13 June 2023

- Fixed read to novel models assignment.

- Improved command line options for providing multiple files, added `--prefix` option.

- Additional checks for various unusual cases in input GTFs.

- Do not output empty files when no GTF is provided.

## IsoQuant 3.2, 27 March 2023

- Unspliced novel transcripts are not reported by the default for ONT data, use `--report_novel_unspliced` to generate them.

- When multiple BAM/FASTQ files are provided via `--bam` / `--fastq`,
    they are treated as different replicas/samples of the same experiment;
    a single GTF and per-sample counts are generated automatically.

- 10-15 times lower RAM consumption with the same running time.

- ~5 times lower disk consumption for temporary files.

- `--low_memory` option has no effect (used by default); `--high_memory` mimics old behavior by storing alignments in RAM.

- Read assignment reports transcript start and end (TSS/TES) matches.

- `--sqanti_output` generates SQANTI-like output for novel vs reference transcripts.

- Resulting annotation contains exon ids.

- Supplementary gene attributes are copied from the reference annotation to the output annotations.

- Improved `--resume` and `--force` behaviour.

- `--model_construction_strategy sensitive_pacbio` is now more sensitive.

## IsoQuant 3.1.2, 7 February 2023

- Fixed strand detection that caused lower precision for novel transcripts.

- Fixed known transcript filtering that caused lower recall.

- Fixed duplicate transcript entries in the output annotation.

- Fixed duplicate canonical attribute in extended annotation.

- Fix `--resume` option when relative paths were provided.

## IsoQuant 3.1.1, 16 January 2023

- Fixed error caused by introns of length 0 (strange corner case, but it does happen).

- Fixed error when using a read grouping file.

## IsoQuant 3.1.0, 3 January 2023

- Implement `--resume` option for resuming failed runs.

- Fix SQANTI-like output for raw reads.

- Fix read strand detection, improves transcript discovery as well.

## IsoQuant 3.0.3, 17 October 2022

- Simplify transcript naming, IDs of known transcripts are preserved in the output.

- More information about novel transcripts in GTF

## IsoQuant 3.0.2, 10 October 2022

- Fix GTF attributes, thanks to @rsalz.

## IsoQuant 3.0.1, 31 August 2022

- Fix `--check_canonical` option.

## IsoQuant 3.0.0, 31 August 2022

- Annotation-free mode for de novo transcript discovery.

- Significant speed-up.

- Extended annotation (all reference + novel transcripts) is now part of the output.

- Intermediate BAM files have nicer names.

- Proper single-thread mode without thread pool usage.

## IsoQuant 2.3.0, 27 May 2022

- New options for controlling quantification strategies. Default behaviour is changed as well. 

- New option `--genedb_output` for providing a separate folder for gene database in case the
    output directory is located on a shared disk.

- Possibility to provide read group tables in gzipped format.

## IsoQuant 2.2.2, 18 May 2022

- Fixed `--check_canonical` option.

- Improved running time for the read assignment step (noticeable only for genes with > 100 exons).

## IsoQuant 2.2.1, 28 Apr 2022

- Minor fixes and improvement in output files. Note, that GTFs and some other files have now multiline headers.

## IsoQuant 2.2.0, 5 Mar 2022

- Parallel processing of transcript model construction phase.

- Minor improvements in quantification of reference transcripts.

## IsoQuant 2.1.1, 7 Feb 2022

- Fixed counts/TPM for novel transcript models.

- Fixed processing of BAM records without sequence data (e.g. secondary alignment).

- Fixed `list index out of range` bug in long read counter.

## IsoQuant 2.1.0, 19 Jan 2022

- Improved recall by introducing relative coverage cutoffs.

- More careful handling of transcript terminal positions.

- Fixed GTF to BED conversion.

## IsoQuant 2.0.0, 7 Oct 2021

- Completely new transcript discovery algorithm with significantly higher recall.

- Algorithm for read alignment correction.

- Support for technical replicas within a single sample.

## IsoQuant 1.3.0, 19 Jun 2021

- Significantly improved running time and RAM consumption;

- Annotation is now fed into minimap2;

- Extended output format.

## IsoQuant 1.2.2, 21 May 2021

- Support for GFF3 mRNA features. 

## IsoQuant 1.2.1, 21 Apr 2021

- Support for BAM files with =/X in CIGAR strings; 

- Fixed canonical splice site detection.

## IsoQuant 1.2.0, 28 Mar 2021

- Multi-threading;

- Intermediate results are saved to disc to enable quick restart via --read_assignments option;

- Significantly improved precision for novel transcript detection;

- Secondary alignments are now used by default;

- Fixed several bugs in inconsistency detection algorithm;

- Reworked polyA detection and reporting once again;

- Slightly modified read assignment output format;

- More informative GTF output;

- Removed --has_polya option, --polya_trimmed is now used as the opposite;

- Added --check_canonical option.

## IsoQuant 1.1.0, 11 Dec 2020

- Significantly reworked polyA detection and reporting;

- Improved detection of inconsistencies, added several new event types;

- Better recall and precision for read assignment algorithm;

- Fixed several bug and flaws;

- Added script for counting simple stats for GTF files (srt/gtf_stats.py).

## IsoQuant 1.0.0, 12 Jul 2020

- Initial release.

- [IsoQuant GitHub page](https://github.com/ablab/IsoQuant)

