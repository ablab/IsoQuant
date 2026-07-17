# Transcript Model Construction

Genome-based discovery of expressed transcript models from long-read assignments.
Entry point: `GraphBasedModelConstructor.process(read_assignment_storage)` in
`isoquant_lib/graph_based_model_construction.py` — a thin orchestrator over the
`isoquant_lib/model_construction/` package. A new constructor (and `ModelStore`) is
built **per gene block**; `process()` runs once per instance.

## Package layout (`isoquant_lib/model_construction/`)

| Module | Class | Responsibility |
|---|---|---|
| `context.py` | `ModelConstructionContext` | Read-only inputs (gene_info, args, string_pools, id_distributor, polyA/TSS predictions) + the intron graph (via `build_graph`) + shared helpers (`get_transcript_id`, `select_reference_gene`, `get_known_spliced_isoforms`). Owns `.store`. |
| `context.py` | `ModelStore` | Mutable model list + per-read construction bookkeeping + read-binding primitives. Single source of truth. |
| `intron_graph.py` | `IntronGraph`, `TerminalVertex`, `IntronCollector` | Intron graph over corrected read introns, with terminal polyA/TSS/read-end vertices. |
| `intron_path.py` | `IntronPathStorage`, `IntronPathProcessor` | Thread a read's introns into a graph path + bucket reads by path (`fl_paths`, `paths_to_reads`). |
| `fl_graph_constructor.py` | `FLGraphConstructor` | Stage 1 — full-length isoforms from graph FL paths. |
| `assignment_based_constructor.py` | `AssignmentBasedConstructor` | Stage 2 — known spliced non-FL, known/novel mono-exon models. |
| `end_processor.py` | `TranscriptEndProcessor` | Stage 3 — polyA/TSS end refinement + alternative-end NIC discovery. |
| `model_filter.py` | `ModelFilter` | Stage 4a — structural filtering + splice-site correction (delegates end refinement to stage 3). |
| `read_assigner.py` | `ConstructionReadAssigner` | Stage 4b — construction-phase read (re)assignment + ownership resolution. |
| `model_read_counter.py` | `ModelReadCounter` | Stage 4c — final honest per-read assignment → counters + model read_info. |
| `transcript_printer.py` | `GFFPrinter`, … | GTF/GFF output for discovered models. |
| `transcript_to_gene_joiner.py` | `TranscriptToGeneJoiner` | Group discovered transcripts into (novel) genes. |
| `transcript_splice_site_corrector.py` | — | Shift novel splice sites onto canonical motifs from read-deletion clusters. |

`graph_based_model_construction.py` keeps `GraphBasedModelConstructor` (the facade,
same public interface used by `parallel_workers.py`) and re-exports
`StrandnessReportingLevel` for `isoquant.py`.

## Workflow (`process()` — exact, load-bearing order)

```
GraphBasedModelConstructor.process(read_assignment_storage)
│
├─ ctx.build_graph(reads)                    ModelConstructionContext
│     IntronGraph → IntronPathProcessor → IntronPathStorage.fill()
│     known_isoforms_in_graph / known_introns
│
│  ── CONSTRUCTION (each appends models to store + binds reads) ──
├─ 1. construct_fl_isoforms()                FLGraphConstructor
│        graph FL paths → known / alt-end NIC / novel (NIC/NNIC)
│        (_reference_isoform_for_path gates known-vs-alt-end reclassification)
├─ 2. construct_assignment_based_isoforms()  AssignmentBasedConstructor
│        partitions unique reads, then:
│          • construct_monoexon_isoforms()     known mono-exon
│          • construct_nonfl_isoforms()        known spliced non-FL (fills nonfl_known_reads)
│          • construct_monoexon_novel()        novel mono-exon from graph terminal exons
│
│  ── POST-PROCESSING ──
├─ 3. pre_filter_transcripts()               ModelFilter (4a)   coverage/mapq prefilter
├─ 4. assign_reads_to_models()  (pass 1)     ConstructionReadAssigner (4b)
│                                              bind leftover reads (unique-known gate + reprofile)
├─ 5. filter_transcripts()                   ModelFilter (4a) ⇄ TranscriptEndProcessor (3)
│        detect_similar_isoforms()
│        per model: correct_novel_transcript_ends() + _mark_terminal_confirmation()
│        detect_similar_isoforms(remove_unconfirmed=True)   drop ISM truncation artifacts
│        _add_known_alternative_end_models() + _drop_duplicate_alt_end_models()
├─ 6. correct_transcript_splice_sites()      ModelFilter (4a)
├─ 7. assign_reads_to_models()  (pass 2)     ConstructionReadAssigner (4b)
│                                              re-bind reads freed by filtering
├─ 8. resolve_read_ownership_conflicts()     ConstructionReadAssigner (4b)
│                                              hand shared reads back to FL models; drop empties
├─ 9. TranscriptToGeneJoiner.join_transcripts()
├─ 10. build_model_read_assignments()        ModelReadCounter (4c)
│         FINAL honest per-read assignment → discovered transcript/gene counters
│         + model_read_assignments (read_info output; see model-read-info-output branch)
└─ 11. compare_models_with_known()           facade (sqanti only)
```

## Read-binding model (`ModelStore`)

```
save_assigned_read(read, tid):
    transcript_read_ids[tid].append(read)      # reads backing model tid
    internal_counter[tid]        += 1          # coverage used by filters
    read_assignment_counts[read] += 1          # how many models claim this read
    construction_assignment[read] = tid        # the model this read built (last-wins)
```

- Construction (steps 1–2) and both `assign_reads_to_models` passes mutate these dicts.
- `move_read_to_model` / `delete_from_storage` keep them consistent when models are
  reshaped (alt-end NICs) or dropped (filtering).
- **Final counts are re-derived in step 10** (per read, against the surviving model set),
  so the construction bookkeeping only decides which models *survive*, not the reported
  numbers. `ReadAssignmentType.discarded` marks a read whose known isoform was filtered
  but whose gene survives (counts toward the gene only).
- Cross-gene-block dedup of already-emitted known isoforms lives in
  `ModelStore.detected_known_isoforms` — a worker-shared class var, never reset.

## Refactor note

This package was extracted from a former ~1500-line `GraphBasedModelConstructor` god
class in 11 behavior-preserving commits (byte-identical outputs verified per step vs.
`--genedb` / no-genedb / `-t 1` runs on `tests/simple_data`). Extraction pattern: each
stage takes the shared `ctx` (and `ctx.store`), mirrors the frequently-used read-only
inputs onto `self` so method bodies stay `self.args` / `self.gene_info`, and reads
graph-derived state live via `self.ctx`. See `.claude/REFINE_TRANSCRIPT_ENDS.md` and
`.claude/POLYA_TSS_DETECTION.md` for the stage-3 terminal logic.
