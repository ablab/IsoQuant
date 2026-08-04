# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## About IsoQuant

IsoQuant is a bioinformatics tool for genome-based analysis of long RNA reads (PacBio, Oxford Nanopore). It reconstructs and quantifies transcript models, assigns reads to annotated isoforms based on intron/exon structure, and performs gene/isoform/exon/intron quantification. Supports bulk, single-cell (10x), and spatial transcriptomics modes.

## Working Guidelines

**Commit Message Style**:
- Use single-line commit messages only (no multi-line body)
- Keep messages short and concise
- Use lowercase for fixes: "fix serialization bug"
- Capitalize for features: "Add barcode_spot grouping"
- No periods at the end, use present tense
- Examples: "fix serialization and barcode split", "add unit tests for string_pools"

**Code Style**:
- Add Python type hints to all new and refactored code
- Use proper type annotations for function signatures and class attributes

**Testing**:
- Do NOT run tests automatically unless explicitly requested
- Tests should only be executed when the developer asks for them

## Development Commands

### Testing

Run all tests with coverage:
```bash
tox
# or directly:
pytest --cov --cov-branch
```

Run a specific test file:
```bash
pytest tests/test_long_read_assigner.py -v
```

Run a specific test function:
```bash
pytest tests/console_test.py::test_clean_start -v
```

Run integration tests (end-to-end pipeline):
```bash
pytest tests/console_test.py -v
```

### Running IsoQuant

Basic test run with toy data:
```bash
./isoquant.py --test
```

Typical development run with simple test data:
```bash
python isoquant.py \
  --reference tests/simple_data/chr9.4M.fa.gz \
  --genedb tests/simple_data/chr9.4M.gtf.gz --complete_genedb \
  --fastq tests/simple_data/chr9.4M.ont.sim.fq.gz \
  --data_type nanopore \
  -o test_output \
  -t 2 \
  --prefix test_sample
```

Run with pre-aligned BAM:
```bash
python isoquant.py \
  --reference tests/simple_data/chr9.4M.fa.gz \
  --genedb tests/simple_data/chr9.4M.gtf.gz \
  --bam tests/simple_data/chr9.4M.ont.sim.polya.bam \
  --data_type nanopore \
  -o test_output
```

### Dependencies

Install development dependencies:
```bash
pip install -r requirements.txt
pip install -r requirements_tests.txt
```

External tools required (must be in PATH):
- minimap2 (for read alignment)
- samtools (for BAM manipulation)

## Code Architecture

### Pipeline Flow

Entry point: `isoquant.py` → `main()` → `run_pipeline()`

Pipeline stages execute sequentially:
1. **Barcode Calling** (single-cell/spatial modes only) - `isoquant_lib/barcode_calling/detect_barcodes.py`
2. **Reference Preparation** - Decompress genome if needed
3. **GTF/GFF Conversion** - Convert annotation to gffutils database format
4. **Read Mapping** - Map reads using minimap2 or STAR (via `isoquant_lib/utils/read_mapper.py`)
5. **Isoform Assignment** - Core analysis (via `isoquant_lib/dataset_processor.py`)
6. **Count Aggregation** - Combine counts across samples

### Core Processing Pipeline

The main processing happens in `DatasetProcessor` (`isoquant_lib/dataset_processor.py`), which orchestrates:

1. **Per-chromosome parallel processing**:
   - `AlignmentProcessor` (`isoquant_lib/alignment/alignment_processor.py`) - Collects alignments per gene, handles primary/secondary/supplementary alignments, merges BAM files
   - `LongReadAssigner` (`isoquant_lib/assignment/long_read_assigner.py`) - Assigns reads to isoforms using junction comparison and profile matching
   - `GraphBasedModelConstructor` (`isoquant_lib/model_construction/model_construction.py`) - Discovers novel transcript models from intron graphs
   - `LongReadCounter` (`isoquant_lib/quantification/long_read_counter.py`) - Counts reads/transcripts/exons/introns with multiple strategies

2. **Data structures**:
   - `GeneInfo` (`isoquant_lib/gene_info.py`) - Reference gene and transcript models
   - `ReadAssignment` (`isoquant_lib/assignment/isoform_assignment.py`) - Links reads to isoforms with match quality
   - `TranscriptModel` - Novel transcript representations

3. **File organization**:
   - `InputDataStorage` (`isoquant_lib/utils/input_data_storage.py`) - Manages input files and metadata
   - `file_naming.py` - Systematic naming for per-chromosome intermediates, lock files for crash recovery

### Key Enums and Modes

**IsoQuantMode** (`isoquant_lib/modes.py`):
- `bulk` - Standard bulk RNA-seq
- `tenX_v3`, `tenX_v2`, `curio` - Single-cell modes
- `tenX_v3_split`, `tenX_v2_split` - 10x split modes for concatenated ONT reads where multiple cDNA molecules are ligated end-to-end. Detects multiple barcode/UMI/TSO patterns per read, splits into individual molecule FASTA records. Uses `TenXSplittingBarcodeDetector` (see `.claude/BARCODE_CALLING.md` for algorithm details)
- `stereoseq`, `stereoseq_nosplit`, `visium_hd`, `visium_5prime` - Spatial transcriptomics

Different modes trigger different processing pipelines (barcode calling, UMI deduplication, etc.)

**ReadAssignmentType** (`isoquant_lib/assignment/isoform_assignment.py`):
- `unique` - Assigned to single isoform
- `ambiguous` - Multiple equally good matches
- `inconsistent` - Mismatches with reference
- `intergenic`, `noninformative` - Unassigned reads
- `inconsistent_genic` / `inconsistent_multigenic` - transcript-unassigned reads that
  overlap gene body/bodies but resemble no isoform. Only ever set as
  `gene_assignment_type` (transcript `assignment_type` stays `noninformative`); the
  `_genic` variant mirrors `inconsistent` (one gene), `_multigenic` mirrors
  `inconsistent_ambiguous` (>1 gene, so it's in both `is_inconsistent()` and
  `is_ambiguous()`). See `.claude/READ_ASSIGNMENT_LIFECYCLE.md`.

**MatchClassification** (SQANTI-like categories):
- `full_splice_match` (FSM), `incomplete_splice_match` (ISM)
- `novel_in_catalog` (NIC), `novel_not_in_catalog` (NNIC)

**CountingStrategy** (`isoquant_lib/quantification/long_read_counter.py`):
- `unique_only` - Only uniquely assigned reads
- `with_ambiguous` - Include ambiguous assignments
- `fsm_only` - Only full splice matches (subset of unique_only): reads with an `fsm` or `mono_exon_match` event
- `all` - Include inconsistent reads

**AmbiguityResolvingMethod** (`isoquant_lib/assignment/long_read_assigner.py`):
- `none`, `monoexon_only`, `monoexon_and_fsm`, `all`

### Parallelization

Chromosome-level parallelization using `ProcessPoolExecutor`:
- `collect_reads_in_parallel()`
- `construct_models_in_parallel()`
- `filter_umis_in_parallel()`

Lock file mechanism prevents race conditions and enables crash recovery/resume.

### Output Formats

**Assignment Printers** (`isoquant_lib/assignment/assignment_io.py`):
- `BEDPrinter` - BED format for read assignments
- `TSVPrinter` - Tab-separated read assignments
- `SQANTIPrinter` - SQANTI-compatible classification

**Grouped Counts** (for single-cell/spatial):
- Default TSV format
- Linear format (optimized)
- Matrix Market (MTX) format for sparse matrices

### Serialization and Resume

`isoquant_lib/utils/serialization.py` handles binary serialization of:
- Gene info databases
- Read assignments
- Intermediate processing results

Enables checkpoint/resume functionality for long-running analyses.

## Important Implementation Details

### Multi-Group Support

**Architecture Overview**:

IsoQuant supports multiple simultaneous grouping strategies (e.g., `--read_group tag:CB file_name barcode_spot`). Each strategy creates separate counters for independent quantification.

**Key Implementation Details**:

1. **Data Structures**:
   - `all_read_groups`: `list[set]` - One set of group IDs per grouping strategy
   - `ReadAssignment.read_group`: `list[str]` - One group ID per strategy for multi-group, single `str` for single group
   - Each counter knows its `group_index` to extract the appropriate group ID from the list

2. **CompositeCounter Architecture** (`isoquant_lib/quantification/long_read_counter.py`):
   - All counters (ungrouped and grouped) are contained within global `CompositeCounter` objects
   - Three global composite counters in `DatasetProcessor`:
     - `global_counter` - Contains gene, transcript, exon, intron counters
     - `transcript_model_global_counter` - Contains transcript model counters
     - `gene_model_global_counter` - Contains gene model counters
   - When methods like `add_read_info_raw(read_id, feature_ids, group_ids)` are called on a `CompositeCounter`, it automatically forwards to all internal counters
   - Each internal `AssignedFeatureCounter` extracts its own group ID using `group_ids[self.group_index]`

3. **Counter Initialization** (`isoquant_lib/dataset_processor.py`, lines 337-407):
   - Ungrouped counters created first and added to composite counter
   - For each grouping strategy, grouped counters created with `group_index` parameter
   - Each grouped counter added to the appropriate global composite counter
   - Example:
     ```python
     # Ungrouped counters
     self.global_counter = CompositeCounter()
     self.global_counter.add_counters([gene_counter, transcript_counter])

     # Grouped counters (one per strategy)
     for group_idx, strategy_name in enumerate(grouping_strategy_names):
         gene_counter = create_gene_counter(..., group_index=group_idx)
         transcript_counter = create_transcript_counter(..., group_index=group_idx)
         self.global_counter.add_counters([gene_counter, transcript_counter])
     ```

4. **GraphBasedModelConstructor Integration** (`isoquant_lib/model_construction/model_construction.py`):
   - Receives only the global `CompositeCounter` objects (NOT separate lists of grouped counters)
   - Calls `add_read_info_raw(read_id, feature_ids, read_assignment.read_group)` with the full group ID list
   - `CompositeCounter` automatically distributes to all internal counters
   - Each counter extracts its appropriate group using `self.group_index`
   - This eliminates manual distribution logic - the composite counter handles everything

5. **Other Implementation Notes**:
   - File naming uses strategy-based approach (see `isoquant_lib/utils/file_naming.py`)
   - Technical replicas check uses `file_name` group only (identified by `grouping_strategy_names.index("file_name")`)
   - Serialization saves/loads `all_read_groups` as list of sets
   - Per-chromosome group files save number of strategies + semicolon-separated groups per strategy

**Benefits of CompositeCounter Architecture**:
- Cleaner separation of concerns - callers don't need to know about multiple counters
- Automatic distribution to all counters (ungrouped and grouped)
- Easy to add new grouping strategies without changing caller code
- Each counter is independent and self-contained

### Illumina Integration

IsoQuant can use short Illumina reads (`--illumina_bam`) to correct long-read exon boundaries via `isoquant_lib/assignment/illumina_exon_corrector.py`. Short reads are NOT used for transcript discovery or abundance estimation.

### PolyA Handling

`PolyAUsageStrategies` enum controls polyA tail usage:
- Tool expects reads to contain polyA tails
- Don't trim polyA tails for better transcript model construction
- PolyA verification integrated into assignment process

## Testing Strategy

**Unit Tests** (`tests/test_*.py`):
- Cover individual modules: alignment info, junction comparison, intron graphs, serialization, etc.
- Run with: `pytest tests/test_<module>.py`

**Integration Tests** (`tests/console_test.py`):
- End-to-end pipeline validation
- `test_clean_start()` - Fresh run with full options
- `test_usual_start()` - Standard run
- `test_with_bam_and_polya()` - BAM input with polyA
- `test_with_illumina()` - Illumina correction
- `test_with_yaml()` - YAML config input

**Test Data**:
- `tests/simple_data/` - Small chr9 region for quick tests
- `tests/toy_data/` - MAPT.Mouse reference for examples

## Git Workflow

Main branch: `master`
Current development: `isoquant_4.0.0` branch

CI/CD workflows in `.github/workflows/`:
- `Unit_tests.yml` - Primary test suite (Python 3.8, Ubuntu)
- `Group_tests.yml` - Multi-group functionality
- `YAML_tests.yml` - Configuration parsing
- Various performance and resume tests
- `Barcode.*.yml` - Protocol-specific barcode calling tests (10x, Curio, Stereo-seq, Visium HD)
- `Barcode.Custom.*.yml` - Universal barcode calling tests via custom_sc mode with MDF files
- `SC.Mouse.10x.yml` - 10x single-cell pipeline with allinfo QA
- `SC.Mouse.10x.barcoded_bam.yml` - 10x pipeline using `--barcoded_bam` with CB/UB tags from BAM
- `SC.Mouse.VisiumHD.barcode2barcode.yml` - Visium HD pipeline with `--barcode2barcode` spot-level UMI dedup
- `SC.Human.Curio.custom_sc.yml` - Curio custom_sc pipeline with allinfo QA

## Code Organization by Function

**Input Processing**:
- `isoquant_lib/utils/input_data_storage.py` - Input file organization
- `isoquant_lib/utils/read_mapper.py` - Mapping orchestration
- `isoquant_lib/assignment/read_groups.py` - Grouping/barcode handling

**Alignment Analysis**:
- `isoquant_lib/alignment/alignment_processor.py` - Alignment collection
- `isoquant_lib/alignment/alignment_info.py` - Alignment data structures
- `isoquant_lib/assignment/multimap_resolver.py` - Multi-mapping read handling

**Isoform Assignment**:
- `isoquant_lib/assignment/long_read_assigner.py` - Main assignment logic
- `isoquant_lib/assignment/long_read_profiles.py` - Read profile structures
- `isoquant_lib/assignment/junction_comparator.py` - Splice junction comparison
- `isoquant_lib/assignment/isoform_assignment.py` - Assignment classifications

**Model Construction**:
- `isoquant_lib/model_construction/model_construction.py` - Novel transcript discovery
- `isoquant_lib/model_construction/intron_graph.py` - Intron graph structures
- `isoquant_lib/visualizer/gene_model.py` - Gene model representations

**Quantification**:
- `isoquant_lib/quantification/long_read_counter.py` - Read counting with multiple strategies
- `isoquant_lib/utils/stats.py` - Count aggregation and statistics
- `isoquant_lib/quantification/convert_grouped_counts.py` - Format conversion for grouped data

**Specialized Features**:
- `isoquant_lib/terminal_prediction/polya_finder.py`, `isoquant_lib/terminal_prediction/polya_verification.py` - PolyA detection
- `isoquant_lib/terminal_prediction/cage_finder.py` - CAGE peak integration
- `isoquant_lib/assignment/illumina_exon_corrector.py` - Short-read correction
- `isoquant_lib/barcode_calling/umi_filtering.py` - UMI-based deduplication

**Output**:
- `isoquant_lib/assignment/assignment_io.py` - Assignment output writers
- `isoquant_lib/model_construction/transcript_printer.py` - GTF/GFF output
- `isoquant_lib/visualizer/plot_output.py` - Visualization generation

**Utilities**:
- `isoquant_lib/common.py` - Shared utility functions
- `isoquant_lib/utils/file_utils.py` - File I/O utilities
- `isoquant_lib/utils/serialization.py` - Binary serialization

## Recent Implementation Changes

See `.claude/BARCODE_CALLING.md` for barcode calling architecture (pipeline integration, barcode-spot grouping, barcode2barcode spot-level UMI dedup, barcoded_bam on-the-fly tag reading, universal barcode calling, MDF format, linked elements, result class hierarchy).

See `.claude/BARCODE_GRAPH_CORRECTION.md` for graph-based barcode correction (branch `badger`, ported from https://github.com/algbio/Badger) — the alternative to per-read whitelist matching for stock 10x whitelists, where `min_score == BARCODE_LEN_10X` makes matching degenerate into exact matching (measured on real 600k-read Mouse 10x data against a 7.4M whitelist: recall 74.4% → 86.6%, precision 99.6% → 99.8%, versus 87.35% for the matched 5000-barcode whitelist). Extract barcode windows verbatim (`TenXBarcodeDetector(whitelist_matching=False)`), count them over the run, select cluster centers from the count distribution ∩ whitelist, and correct the rest by edit distance. Selected via `--barcode_correction {auto,whitelist,graph}` (auto → graph above 100k whitelist barcodes) plus `--n_cells`. Architecture: `BarcodeGraph` (`isoquant_lib/barcode_calling/barcode_graph.py`) + the `correct_barcodes` stage (`correct_barcodes.py`), writing `aux/<prefix>.raw_barcodes_<i>.tsv` → `<prefix>.barcoded_reads_<i>.tsv`. Two provably equivalent implementations, `--barcode_graph_impl {centers,full}`; `centers` is the default and 4.9x faster. Also available from the standalone `isoquant_detect_barcodes.py`. Two gotchas documented there, both found by measurement: strand selection must be structural (R1-to-polyT span), because raw mode has no match score to pick an orientation with — this was worth 20 points of recall; and `hits_delta` counts k-mer *occurrences*, so it must be set to `(L-k+1)²` to actually disable hit pruning.

See `.claude/ANALYSIS_OPTION.md` for the `--analysis` interface — the single option (values `quantification`/`quant`, `transcript_discovery`/`td`, `exon_quantification`/`ex_quant`, `fusion`) that selects pipeline stages, with context-aware defaults. It supersedes the now-hidden/deprecated `--count_exons`, `--count_intron_retentions`, `--fusion`, `--no_model_construction`, which still work. Resolution lives in `resolve_analyses()` (`isoquant.py`), producing internal flags `run_quantification` / `count_exons` / `count_intron_retentions` / `fusion` / `no_model_construction` / `predict_terminal_sites`.

See `.claude/JOINT_EXON_COUNTS.md` for joint exon counts — region-based exon quantification that groups overlapping annotated exons into regions and emits N+1 features per region (N inclusion variants + 1 region-level exclusion). Runs alongside the classic `ExonCounter` when exon quantification is enabled (`--analysis exon_quantification`, aka the deprecated `--count_exons`).

See `.claude/FUSION_DETECTION.md` for fusion gene discovery (PR #392, `fusion` branch) — enabled via `--analysis fusion` (deprecated alias: `--fusion`), triggers post-isoform fusion calling via SA-tag breakpoints + soft-clip realignment, with biotype/multicopy/frequency filters, mappy-based reconstruction realignment, and a per-sample TSV report (`<experiment>/<prefix>.fusions.tsv`). Architecture: `FusionDetector` + `FusionValidator` + `FusionMetadata` + `GenomicIntervalIndex`.

See `.claude/POLYA_TSS_DETECTION.md` for polyA / TSS site prediction — per-transcript read-end histograms + scipy peak detection + XGBoost peak filter (shipped models in `isoquant_lib/data/`), gated by `args.predict_terminal_sites` (set in `resolve_analyses`: requires `--genedb` and either `quantification` or `transcript_discovery` active; TSS additionally needs `--fl_data`). Prediction is part of quantification and is also consumed by model construction. Outputs `SAMPLE.polyA_prediction.tsv` / `SAMPLE.TSS_prediction.tsv`. Architecture: `TerminalCounter` + `PolyACounter` / `TSSCounter` in `isoquant_lib/terminal_counter.py`, integrated as ordinary counters in `ReadAssignmentAggregator.global_counter`. The shared peak detector lives in `isoquant_lib/terminal_peaks.py` (`detect_peaks`), also reused by transcript discovery; a per-gene `flush()` exposes `last_gene_predictions` for intron-graph terminal-vertex refinement, while the ungrouped output is predicted whole-chromosome (chromosome-wide `_all_transcripts` buffer, so a transcript whose reads span gene blocks is not emitted as duplicate rows). Companion training doc: `.claude/POLYA_TSS_TRAINING.md` — two-command workflow via the hidden `--collect_polya_training` / `--collect_tss_training` dev-only CLI flags + `misc/train_polya_tss_model.py`.

See `.claude/GFFCOMPARE.md` for how transcript-discovery accuracy is scored with gffcompare (`misc/reduced_db_gffcompare.py` full/known/novel splits vs `expressed`/`expressed_kept`/`excluded` references) and why its Transcript/Intron-chain levels are **end-agnostic** for multi-exon transcripts (so polyA/TSS end refinements are invisible there). Documents the local gffcompare fork at `/home/andreyp/gffcompare_term_delta/` adding `--terminal-delta <d>` (end-sensitive transcript matching, independent of `-e`), and the end-sensitivity benchmark showing the end refinement helps novel transcripts at tight `d`.

See `.claude/REFINE_TRANSCRIPT_ENDS.md` for alternative-polyA/TSS **NIC** discovery (branch `transcript_model_ends`) — **DEFAULT ON whenever `--genedb` is given**; the `--refine_transcript_ends` flag was dropped (commit 518ae62e). Gates (`__init__`): `create_nics=bool(genedb)` (alt-end NIC creation), `use_tss_model=bool(fl_data)` (5' side), `novel_apa=bool(flag)` (hidden); end polishing (`correct_novel_transcript_ends`) is always on — the `polish_ends` flag + the old `_correct_novel_transcript_ends_simple` were removed. Side-appropriate vertex snapping (`intron_graph._refine_positions`) refines 3' terminal vertices toward polyA predictions, 5' toward TSS. Two mechanisms: (A) graph-level — `construct_fl_isoforms` emits a NIC when a goal-1-refined FL-path terminal vertex disagrees with the matched reference end, keeping the known (union; 3' polyA `VERTEX_polya`/`polyt` always, 5' TSS only with `--fl_data`); (B) post-construction per-transcript pass — known-model alt-end NICs from per-model polyA-confirmed-read peak calling, with a relative-support gate + intron-chain/ends dedup. Big win on TSS (Mouse novel 0→~59%, Human →~53%, full_recall@td10 up, known preserved); Human polyA 0→37–41% @ ~94%; Mouse polyA / SIRVs neutral. De-novo (no genedb): only end *polishing* runs (no NIC creation) — proven byte-identical to master on Mouse no_gtf; `--novel_apa` extends to novel chains (inert on these sims). New de-novo CI test `Human.ONT_simulated.polyA_2.no_gtf` (full 47.9/95.3). Part-3 diagnostic (`tme_no_postpass`): reduced_db precision drop is mostly intrinsic to the graph-level NIC (the win-bearing step), so the current design is kept. c5950439 is local-only; parked early post-pass on branch `goal2_end_refinement_wip`.

## Planned Optimizations

### String Interning for Memory Efficiency (Planned)

**Motivation**: With billion-read datasets, ReadAssignment and IsoformMatch objects consume ~600GB of memory, mostly from duplicated strings (gene IDs, transcript IDs, barcodes, read groups).

**Strategy**: Replace strings with integer indices referencing shared string pools.

See detailed analysis and implementation plan in: `.claude/STRING_INTERNING_OPTIMIZATION.md`

**Expected Impact**:
- Memory reduction: 78% (700 GB → 150 GB for 1B reads)
- Disk usage reduction: 80-90% for serialized files
- Minimal performance impact: < 5% slowdown

**Key Design Points**:
- Per-worker string pools during parallel processing
- Pre-populate known pools (genes, transcripts, chromosomes) from GTF
- Discover dynamic pools (read groups) during processing
- Property-based API for backward compatibility
- Incremental rollout with optional flag
