# Read assignment & classification lifecycle

How a read goes from alignment to (a) a reference assignment, (b) membership in a
constructed transcript model (which drives the model counts), and (c) a row in the
`transcript_model_reads.tsv` read_info output. Also: when transcript models are
created and filtered.

All code is in `isoquant_lib/graph_based_model_construction.py` unless noted.

## The three "assignments" (do not conflate)

| | Where | Assigner mode | Purpose |
|---|---|---|---|
| **A** reference assignment | collect phase (`alignment_processor.py`, `LongReadAssigner`) | full (`quick_mode=False`) | read vs the **annotation** (`--genedb`) → the "source" `ReadAssignment` |
| **B** model membership | `assign_reads_to_models` + construction | fast (`quick_mode=True`) | which model(s) a read belongs to → **counts** (`transcript_read_ids`) |
| **C** read_info classification | `populate_model_read_assignments` | full (`quick_mode=False`) or reuse of A | how each bound read matches its model → `transcript_model_reads.tsv` |

`transcript_read_ids` (membership, set in phase 2) is the single source of truth for
counts. Phase 3 never recomputes membership — it only *describes* each bound read.

## Phase 1 — Reference assignment (before any model exists)

`LongReadAssigner.assign_to_isoform` (built at `alignment_processor.py`, `quick_mode=False`)
compares each read's profile to the reference `gene_info`:
- type: `unique | unique_minor_difference | ambiguous | inconsistent | noninformative`
- class: `FSM | ISM | NIC | NNIC | genic | intergenic` (+ match events)

The resulting `ReadAssignment` ("source") is serialized (carries `corrected_exons`,
`polya_info`, barcode/UMI, read groups, …) and is the input to construction.

### Genic reads without a transcript — the two `gene_assignment_type` tiers

A read can overlap a gene body yet resemble no isoform: the `genic` / `genic_intron`
branches of `assign_to_isoform`, plus the three empty exits of `match_inconsistent`
(`select_similar_isoforms` / `detect_inconsistensies` / `select_best_among_inconsistent`
all empty). These go through `assign_to_overlapping_genes`, which keeps the transcript
**unassigned** (`assignment_type = noninformative`, no `assigned_transcript`) but attaches
one `IsoformMatch` per overlapping gene (via `gene_info.get_gene_regions()`), classed
`genic` / `genic_intron`.

`ReadAssignment.__init__` then derives `gene_assignment_type` from the distinct gene count:
- 1 gene  → `inconsistent_genic`     (mirrors `inconsistent`)
- >1 gene → `inconsistent_multigenic` (mirrors `inconsistent_ambiguous`; in `is_ambiguous()`)
- 0 genes (true intergenic) → stays `noninformative` / `intergenic`

Both new types are in `is_inconsistent()`, so the **gene** counter counts them under the
`use_inconsistent` strategies (`all`, `unique_inconsistent`): 1.0 for one gene, 1/N split
across N genes. This works because `AssignedFeatureCounter.add_read_info` short-circuits on
the **extractor's** assignment type (`GeneAssignmentExtractor.get_assignment_type` →
`gene_assignment_type`), not the raw transcript `assignment_type`. Transcript counting and
model construction key on `assignment_type` (still `noninformative`), so they are unaffected;
default `--gene_quantification` (`unique_splicing_consistent`) leaves per-gene counts
unchanged (these reads weigh 0 there — only the `__no_feature` summary row shrinks, since
they now carry a gene). History: gene attribution added in b2895a7d; the two dedicated
`inconsistent_*genic` types + gene counting in 5ea2b554.

## Phase 2 — Model construction & filtering (`GraphBasedModelConstructor.process`, per gene block)

```
 1. IntronGraph(reads)                    collapse splice sites within tolerance into vertices;
                                          refine terminal (polyA/TSS) vertices
 2. IntronPathStorage.fill(reads)         each read threaded onto EXACTLY ONE intron path;
                                          full-length (FL) paths flagged
 3. get_known_spliced_isoforms            which annotated isoforms appear in the graph

 4. construct_fl_isoforms()               CREATE models from FL paths (ordered by coverage):
      • path == a known isoform  -> KNOWN model (transcript_from_reference), if cov >= min_known_count
      • else -> NOVEL model: NIC if all introns known else NNIC,
                if cov >= min_novel_count and strand/polyA gates pass
      • save_assigned_read binds the path's reads   [B, construction: each read in ONE model]
 5. construct_assignment_based_isoforms() CREATE mono-exon / by-reference-isoform models from
                                          reads not yet consumed; binds them (skips read if
                                          read_assignment_counts>0)

 6. pre_filter_transcripts()              DROP short (<=2 exon) NOVEL models below coverage/mapq
                                          (delete_from_storage unbinds their reads)

 7. assign_reads_to_models()  1st pass    [B, re-profile] for reads NOT bound in construction,
                                          FULL model set, quick_mode=True:
        reference-unique to a KNOWN iso?
            iso survived -> keep read on it                 (rule 1)
            iso dropped  -> unassign (read_assignment_counts=0)  (rule 2)
        else re-profile: is_consistent() -> bind to matched model(s) (may be several)
                         else            -> unassign             (rule 4-membership)

 8. filter_transcripts()                  DROP models:
        • detect_similar_isoforms -> substitute/drop redundant novel models (two passes)
        • novel cov < cutoff(component coverage) or mapq < cutoff -> drop
        • correct_novel_transcript_ends (polyA/TSS end refinement; always on)
        • (+ _add_known_alternative_end_models if create_nics/novel_apa)
                                          (every drop unbinds that model's reads)

 9. assign_reads_to_models()  2nd pass    re-bind the reads freed by step 8 (same logic as 7)

10. join_transcripts()                    group surviving models into genes
11. forward_counts()                      membership (transcript_read_ids) -> feed the
                                          transcript_model & gene_model COUNTERS
```

Create vs filter vs bind timeline:
```
CREATE:  construct_fl_isoforms -> construct_assignment_based_isoforms -> (+alt-end NIC models)
FILTER:  pre_filter_transcripts -> filter_transcripts (dedup x2 + coverage/mapq)
BIND:    construction -> assign_reads_to_models(1) -> [filter frees reads] -> assign_reads_to_models(2)
```

## Phase 3 — read_info output (`populate_model_read_assignments`)

Only with `--large_output read2transcripts`. Membership is read from `transcript_read_ids`
(NOT recomputed), so this output is byte-consistent with the counts and never changes them.
Written per chromosome by a `ReadInfoPrinter` (`assignment_io.py`) and merged into
`SAMPLE.transcript_model_reads.tsv[.gz]` — same 18-column format as `SAMPLE.read_info.tsv`,
so `misc`/`convert_read_info.py` conversions work on it. `isoform_id` holds the model.

Per bound read on model M, the read_info classification/events are obtained as:

1. **Reuse** — M is a **known** model and the read's source (phase-1) assignment was
   `unique`/`unique_minor_difference` on that same transcript → reuse the source
   `IsoformMatch` verbatim (real class + events; no re-profiling can spuriously drop it).
2. **Reprofile** — otherwise → classify the read against a **single-model** `GeneInfo` with
   `LongReadAssigner(..., quick_mode=False)` (same setting as phase 1), giving real
   FSM/ISM/NIC/NNIC + events, always assigned to M.

Reads bound to several models → `ambiguous`, one row per model. Reads on no model
(dropped or never bound) → `noninformative` (`isoform_id="."`).

Wiring: `parallel_workers.construct_models_in_parallel` creates the per-chr
`ReadInfoPrinter` and calls `populate_model_read_assignments` (gated on the flag);
`dataset_processor` creates the final printer and merges per-chr files (header_lines=3).
The old 2-column `read_id\ttranscript_id` format and `GFFPrinter.dump_read_assignments`
were removed.

### The four rules (what phases 2–3 guarantee)

1. Read that contributes to a **kept known** model → keep information (reuse the source match).
2. Read reference-uniquely assigned to a **dropped known** model → dropped / noninformative.
3. Read that contributes to a **novel** model → proper reassignment via the single-model assigner.
4. All other reads → re-assigned "honestly" (`quick_mode=False`), never blanked.

## `quick_mode` — the subtlety

`match_inconsistent` (`long_read_assigner.py`) short-circuits to `inconsistent`/`genic`
with **no assigned transcript** when `quick_mode=True`; `quick_mode=False` runs the full
inconsistency detection (NIC/NNIC + events, with the assigned transcript).

- Phase 1 (A) and phase 3 (C) use `quick_mode=False` → real classifications.
- The membership re-profile in `assign_reads_to_models` (B) uses `quick_mode=True`. This is
  the only place quick mode still drives a *binding* decision: a read whose sole relation to
  the models is a minor difference that trips `match_inconsistent` (beyond `delta`) is
  **dropped** here, whereas `quick_mode=False` could classify it `unique_minor_difference`
  and bind it. Switching B to `quick_mode=False` would bind more reads → **changes counts**;
  left as-is by default.
- The genic/blank-`isoform_id` bug in an earlier read_info draft was exactly this: phase 3
  had inherited `quick_mode=True` from B. Fixed by using `quick_mode=False` in phase 3.

## Key functions

- `process` — orchestrates phase 2 (order above).
- `construct_fl_isoforms`, `construct_assignment_based_isoforms` — create models + bind (`save_assigned_read`).
- `pre_filter_transcripts`, `filter_transcripts`, `detect_similar_isoforms`, `delete_from_storage` — filtering (unbind on drop).
- `assign_reads_to_models` — membership (B); `_reference_unique_known_isoform` is the rule-1/2 gate.
- `forward_counts` — membership → counters (counts).
- `populate_model_read_assignments`, `_copy_read_aux` — read_info output (C).
