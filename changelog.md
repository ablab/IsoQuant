# IsoQuant changelog

## IsoQuant 3.10.0, 21 October 2025

- New option `--polya_trimmed` to indicate polyA-trimmed reads (thanks @hmutpw for the suggestion [#342](https://github.com/ablab/IsoQuant/issues/342)).
- 
- New option `--process_only_chr` to process a specific list of chromosomes.


## IsoQuant 3.9.0, 1 October 2025

- Secondary alignments are not used by default from now on. It significantly improves running time and RAM consumption, but barely affects the results' quality.
Use `--use_secondary` to process secondary alignments.

- New options that force IsoQuant to use only a faction of reads in high-coverage loci.
Significantly improves running time and RAM consumption, but affects gene/isoform counts. 
New default behaviour only affects small chromosomes and scaffolds (<500kbp).

    In some cases, high-coverage regions take too much time to process due to extreme number of mapped reads,
especially `chrM` (up to 10x longer compared to normal chromosomes). However, using only a fraction of these
reads is enough to obtain reliable results.

    These options allow to process only up to given number of reads mapping to a hich-coverage loci on short and normal chromosomes:
  - `--max_coverage_small_chr` (default value is 1 million);
  - `--max_coverage_normal_chr` (default value is infinity, so usual chromosomes are not affected by default even if some genes have extreme coverage).

- New option `--discard_chr` to discard a list chromosomes from the analysis.

## IsoQuant 3.8.0, 8 September 2025

- Fixed `--report_canonical` preset ([#332](https://github.com/ablab/IsoQuant/issues/332), thanks to @wwliao).

- Fixed counts for novel genes in `discovered_gene_counts.tsv` and `discovered_gene_tpm.tsv`
([#337](https://github.com/ablab/IsoQuant/issues/337), thanks to @yjliuhub).

- Fixed `--genedb_output` option ([#335](https://github.com/ablab/IsoQuant/issues/335), thanks to @YalanBi).

## IsoQuant 3.7.1, 5 August 2025

- Support for indexing BAMs with large chromosomes, fixes [#327](https://github.com/ablab/IsoQuant/issues/327). Thanks to @maol-corteva.

- `CDS` features are now used when `exon` features are absent, fixes [#309](https://github.com/ablab/IsoQuant/issues/309).
 
- Chromosome names are now checked for consistency between reference genome, annotation and BAM file provided. 
Only overlapping chromosome names are used if inconsistent.

- Fixed SQANTI-like output headers, fixes [#318](https://github.com/ablab/IsoQuant/issues/318).

- Some minor cosmetics.

## IsoQuant 3.7.0, 16 May 2025

- Optimized grouped counts output. By default, all counts are stored in linear format, which saves time and disk space.
Matrices with small number of columns are automatically converted to usual matrix in TSV format, 
larger matrices typical for single-cell and spatial data are converted to MTX format.
See `--counts_format` paramter for options. It is also possible to convert counts after IsoQuant is finished using `src/convert_grouped_counts.py`.
Fixes issues mentioned in [#248](https://github.com/ablab/IsoQuant/issues/248).

- Renamed counts related to discovered transcripts and genes to avoid confusion.

- New options `--indexing_options` and `--mapping_options` that allow to pass options to the indexing and mapping commands.
Fixes [#284](https://github.com/ablab/IsoQuant/issues/284) and [#259](https://github.com/ablab/IsoQuant/issues/259).

- STARlong is now an alternative options for aligning, can be set via `--aligner starlong` (not recommended for ONT reads).
Fixes [#284](https://github.com/ablab/IsoQuant/issues/284)

- Exon/splice junction counts now only come from reads assigned to the same strand, fixes [#253](https://github.com/ablab/IsoQuant/issues/253).

- Use only gene-assigned reads for exon counting, fixes [#283](https://github.com/ablab/IsoQuant/issues/283).

- Fixed rare serialization bug [#304](https://github.com/ablab/IsoQuant/issues/304).

## IsoQuant 3.6.3, 13 January 2025

- Fix penalty score for terminal exon elongation when selecting similar isoforms for inconsistent reads [#270](https://github.com/ablab/IsoQuant/issues/270).

- Fix `transcript_model_grouped_counts` output format [#275](https://github.com/ablab/IsoQuant/issues/275).

## IsoQuant 3.6.2, 15 November 2024

- Fixes linear grouped counts output [#258](https://github.com/ablab/IsoQuant/issues/258).

## IsoQuant 3.6.1, 25 September 2024

- Import exon attributes from the reference annotation [#175](https://github.com/ablab/IsoQuant/issues/175).

- Fixed annotation checks for GFF3 [#240](https://github.com/ablab/IsoQuant/issues/240).


## IsoQuant 3.6.0, 13 September 2024

- Fixed duplicated `noninformative` and `intergenic` reads assignments.
As a results, fixed duplicated novel transcripts [#236](https://github.com/ablab/IsoQuant/issues/236).

## IsoQuant 3.5.2, 3 September 2024

- Fixes exon counting algorithm [#229](https://github.com/ablab/IsoQuant/issues/229).

## IsoQuant 3.5.1, 26 August 2024

- Fixed YAML support in visualization [#222](https://github.com/ablab/IsoQuant/issues/222).

- Fixed transcript naming when IsoQuant-generated GTF is provided as input [#219](https://github.com/ablab/IsoQuant/issues/219).

- Fixed `exons` attribute duplication [#219](https://github.com/ablab/IsoQuant/issues/219).

- Exon ids are now consistent between input and output annotations if present.

- New `--count_format` option for setting desired grouped counts format (matrix/linear/both), fixes [#223](https://github.com/ablab/IsoQuant/issues/223).


## IsoQuant 3.5.0, 2 August 2024

- New visualization software developed by [@jackfreeman88](https://github.com/jackfreeman88). See more [here](https://ablab.github.io/IsoQuant/visualization.html).

- Dramatically reduced RAM consumption for grouped counts, about 10-20x decrease on datasets with large number of groups.
  Important fix for single-cell data processing. Should fix [#189](https://github.com/ablab/IsoQuant/issues/189).

- Fixed [#195](https://github.com/ablab/IsoQuant/issues/195): output GTF contained very similar isoforms and estimated their expression as 0.

- New documentation is now available at [ablab.github.io/IsoQuant](https://ablab.github.io/IsoQuant/).

## IsoQuant 3.4.2, 13 July 2024

- Dramatically reduced RAM consumption.
Should fix [#209](https://github.com/ablab/IsoQuant/issues/209).

IsoQuant 3.4.2 was tested on a simulated ONT dataset with 30M reads using 12 threads.
In the default mode RAM consumption decreased from 280GB to 12GB when using
the reference annotation and from 230GB down to 6GB in the reference-free mode.
Running time in the default mode increased by approximately 20-25%.
When using `--high_memory` option, running time remains the same as in 3.4.1,
RAM consumption in the reference-based mode is 46GB, and 36GB in the reference-free mode.
Note, that in general RAM consumption depends on the particular data being used and the number of threads.

In brief, in 3.4.0 and 3.4.1 inadequate RAM consumption was caused by
[this commit](https://github.com/ablab/IsoQuant/commit/557e5834d0503587b918a0eedf3ff5cee3253141).
Apparently, adding a couple of `int` fields to the `BasicReadAssignment` class made the default pickle serialization
not to clean used memory (possibly, a leak). Since some large lists of `BasicReadAssignment` were sent between
processes, this caused the main process to consume unnecessary RAM. When later new processes were created
for GTF construction, total RAM consumption exploded thanks to the way Python multiprocessing works.
This release implements two ways fixing the issue: sending objects via disk (default) and
using custom pickle serialization (when `--high_memory` is used).

- Transcript and exon ids are now identical between runs, including ones with different number of threads.

## IsoQuant 3.4.1, 9 May 2024

- Fixed `IndexError: list index out of range` when `--sqanti_output` is set ([#186](https://github.com/ablab/IsoQuant/issues/186)).

- Fixed `IndexError: list index out of range` in printing grouped transcript models TPMs ([#187](https://github.com/ablab/IsoQuant/issues/187)).

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

- Fixed `--resume` option when relative paths were provided.

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

