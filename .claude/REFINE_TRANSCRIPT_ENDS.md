# Alternative polyA/TSS NIC discovery (DEFAULT ON with --genedb)

Branch `transcript_model_ends` (off `polya_tss_on_master_jun`). Commits:
- **66f0b3d** graph-level FL-path NIC;
- **1eda1aa** post-construction per-transcript pass;
- **6f6e4d2f** relative-support gate + polyA-confirmed histograms;
- **518ae62e** made default — dropped the `--refine_transcript_ends` flag (origin HEAD);
- **c5950439** de-novo polishing + de-novo novel generation; keeps the genedb path
  byte-identical (gate decoupling only). **Local-only — not pushed** (origin is at 518ae62e).

There is **no `--refine_transcript_ends` flag** anymore. The behavior is controlled
entirely by whether `--genedb` / `--fl_data` / `--novel_apa` are present.

## Gate flags (`__init__`, graph_based_model_construction.py)
```python
self.create_nics   = bool(getattr(args, 'genedb', None))   # alt-end NIC *creation* (graph + post-pass)
self.use_tss_model = bool(getattr(args, 'fl_data', False)) # use the TSS model (polish + NIC, 5' side)
self.novel_apa     = bool(getattr(args, 'novel_apa', False))  # hidden opt-in; novel-chain siblings
```
- Constructed novel-model ends are **always** refined (no flag): `correct_novel_transcript_ends`
  (the renamed `_refine_novel_transcript_ends`) snaps each unsupported terminal to the dominant
  trained peak. The old `polish_ends` flag + dispatcher and `_correct_novel_transcript_ends_simple`
  were removed; the positional fallback survives only as `_closest_inward`. The polyA model is
  always available in `_terminal_model`.
- `create_nics` → graph-level NIC (`construct_fl_isoforms` / `_reference_isoform_for_path`)
  AND the post-pass known→NIC.
- `use_tss_model` → 5' (TSS) side everywhere (polishing and NIC); needs `--fl_data` because
  read starts are unreliable otherwise.
- Post-pass call site in `filter_transcripts`: runs when `create_nics or novel_apa`.

### Behavior matrix
| Input | Graph NIC | Polishing | Post-pass | Novel siblings |
|---|---|---|---|---|
| **--genedb** (any --fl_data) | yes (3' polyA always; 5' TSS only w/ --fl_data) | yes | yes (Part 2; Part 3 if --novel_apa) | only w/ --novel_apa |
| no genedb, no --novel_apa | — | yes (Part 1 only) | — | — |
| no genedb, **--novel_apa** | — | yes | yes (Part 3 only) | yes |

**genedb path is provably unchanged from 518ae62e:** with `--genedb`, `create_nics=True`,
so `create_nics or novel_apa`, `polish_ends or create_nics`, etc. all evaluate exactly as
the old single `refine_ends` did; `use_tss_model`==old `refine_tss`; `novel_apa`==old value.
c5950439 only changed which gate each call site reads. Verified chr9 genedb output identical.

**De-novo (no genedb):** the goal-1 intron-graph vertex *snapping* needs the per-gene
polyA/TSS predictions, which require `--genedb`, so it is inactive de-novo. De-novo
"polishing" is purely the per-model detector refinement (Part 1), which falls back to the
simple positional correction when there is no confident peak. Proven **byte-identical to
master** at all 3 deltas on Mouse no_gtf → safe, adds nothing measurable. Human polyA_2
de-novo discovers full **47.9/95.3** (def), 45.2/90.0 (td50), 37.5/74.6 (td10). See the
new CI test below.

## Post-construction pass (per-transcript, after assign_reads_to_models)
In `filter_transcripts`, gated by `create_nics or novel_apa`:
- **Part 1 — refine novel ends** (`correct_novel_transcript_ends`, the renamed
  `_refine_novel_transcript_ends`): snap each unsupported novel-model terminal
  to the dominant trained peak; both ends in one pass (concordant). Verdict:
  **neutral** (goal-1 already refined novel ends; de-novo it falls back to the
  `_closest_inward` positional correction).
- **Part 2 — known alt-end NICs** (`_add_known_alternative_end_models` →
  `derive_alternative_end_models`): per known model, peak-call its own reads, emit a
  `.nic` per confident alternative end (multi-site, union). Verdict: **the win** —
  TSS novel_recall +3–5% (Mouse 53.9→58.8, Human 50.1→53.2), known preserved;
  polyA neutral. Catches alt-TSS on knowns the FL-path missed.
- **Part 3 — `--novel_apa`** (hidden, default off): also spin off siblings for
  novel chains (works with or without --genedb). Verdict: **inert** on these sims
  (alt-ends live on known chains); kept as opt-in.
- **Dedup** (`_drop_duplicate_alt_end_models`): drop any non-known model whose
  intron chain + both ends (within `apa_delta`) match a reference isoform or an
  already-kept model — catches graph-level + post-pass NICs; 0 duplicate structures.
- **Precision gating (6f6e4d2f), critical for SIRVs:**
  (1) **relative-support gate** — an alt-end peak must clear `min_novel_count` AND
  `terminal_position_rel × dominant_terminal_peak` (reuse the intron-graph terminal
  threshold; self-normalizing, not tuned);
  (2) **polyA-confirmed histograms** (`_terminal_histograms` / `_confirmed_polya_pos`)
  — the 3' polyA histogram is built from polyA-CONFIRMED reads at the detected
  cleavage site (`external_polya_pos`/`_polyt_pos`), matching `PolyACounter`'s
  trained domain; 5' TSS uses all stranded reads (matches `TSSCounter`). Root cause
  of SIRVs over-calling was feeding the model out-of-domain all-read histograms.
  Result: SIRVs novel_precision 75.0→**85.7** (back to neutral), TSS gains held,
  known preserved, polyA slightly improved.

## What the graph-level mechanism does
In `construct_fl_isoforms`, when an FL path's intron chain matches a reference,
`_reference_isoform_for_path` keeps the annotated **known** UNLESS a *refined terminal
vertex* disagrees with the annotated end by > `apa_delta` → then it returns None and the
path falls through to the novel branch, emitting a **NIC** from `novel_exons` (= the
path's refined ends). Union emerges at the population level: base-end reads and alt-end
reads thread to *different* terminal vertices → separate FL paths → base stays known, alt
becomes NIC. So the known is never lost (fixes the prior reclassify-in-place regression).

### The gate (critical — ungated regresses)
`_reference_isoform_for_path` (gated by `create_nics`):
- **3' polyA side**: NIC only when a *detected polyA vertex* differs — `VERTEX_polyt`
  on genomic-left (3' of `-`), `VERTEX_polya` on genomic-right (3' of `+`). Bare
  `read_end`/`read_start` never trigger here (degraded-ONT-end safety). Needs only `--genedb`.
- **5' TSS side**: NIC only when `use_tss_model` (= `--fl_data`), strand-aware:
  5' = genomic-left `read_start` for `+`, genomic-right `read_end` for `-`.
- Threading (`thread_ends`/`thread_starts`) is UNCHANGED (user kept "closest vertex").
- Side-appropriate snapping (`intron_graph._refine_positions` / `_attach_side`): the
  genomic 3' side (`read_end=True`) snaps terminal vertices toward the **polyA**
  predictions, the 5' side toward the **TSS** predictions. (Earlier the 3' read-end
  vertices were snapped toward TSS predictions — a bug, now fixed.)

Why the gate: ungated (compare both ends, any vertex) lost known_recall −6.9 on normal
reduced_db and −4.9 on Mouse polyA (degraded read-ends → spurious NICs). The
trusted-vertex gate restored known fully. Reuses goal-1's already-refined FL-path terminal
vertices (snapped to the trained polyA/TSS predictions) directly, instead of re-deriving
ends from noisy per-model histograms — this is why graph-level beats the parked post-pass.

## Results (default-on, vs master/no-flag baseline; 3-delta gffcompare)
- **Mouse TSS (--fl_data) — big win**: novel 0→~**59%** (graph 53.9 + post-pass to 58.8) @
  ~94% prec; full_recall@td10 +~22; known preserved 96.5, known_prec +4.1.
- **Human TSS — win**: novel 0→~**53%** @ ~92–95% prec; known preserved.
- **Human polyA_1 / polyA_2 — graph-level win**: novel 0→**37–41%** @ ~94% prec (baseline
  novel ~0 → full headroom); known preserved, known_prec +2.7..+4.5, full_recall@td10 +12..+19.
- **Mouse polyA — ~neutral** (novel 18.9→19.0): the assigner is already end-sensitive for
  polyA so the baseline already recovers alt-polyA; little headroom. (Data only 49.6%
  polyA → `requires_polya_for_construction=False`.)
- **SIRVs — neutral** (novel_precision 85.7, restored by the precision gates).
- **Normal reduced_db (reduced-CHAIN truth, no alt-ends)** — small recall + (graph-level
  +1.3–1.5), modest precision cost (see Part-3 diagnostic). This is the price of the
  discovery win, not a bug.
- Only systematic cost: full_precision at the DEFAULT (end-agnostic) delta dips — an
  artifact (alt-end NICs share the chain → gffcompare reads them as duplicates); it
  improves/flips at td10.

## Part-3 diagnostic — is the reduced_db precision drop avoidable?
Throwaway branch **`tme_no_postpass`** (commit 1652e98f, local-only): graph-level NIC +
gate + polishing kept, only the post-pass Part-2 known→NIC disabled. Non-polyA reduced_db
tests, novel recall/precision, master → graph-only → graph+Part2:

| config | master | graph-only | graph+Part2 |
|---|---|---|---|
| Mouse R10.reduced_db | 61.1 / 85.9 | 62.6 / 83.4 | 62.6 / 81.1 |
| Mouse PB.reduced_db | 76.2 / 94.3 | 77.5 / 91.5 | 77.5 / 90.9 |
| Mouse R9.trunc.reduced_db | 62.0 / 88.5 | 63.3 / 83.4 | 63.3 / 83.3 |

**Conclusion:** the precision drop is MOSTLY INTRINSIC to the graph-level NIC (P −2.5 /
−2.8 / −5.1) — the same mechanism that yields the polyA/TSS wins and the +1.3–1.5 recall.
The post-pass Part 2 adds **0 recall** and only a small/variable extra precision cost
(R10 −2.3, PB −0.6, R9.trunc −0.1), in exchange for the +5% TSS gain on fl_data configs.
**Decision: keep the current design** (graph + Part 2, default-on with --genedb). The
precision sacrifice is unavoidable + partly a reduced-CHAIN truth artifact (alt-end NICs
are penalized as FP against a held-out-chain reference); the wins justify it.

## Benchmarks / how to run
6 alt-end sim configs in `/abga/work/andreyp/ci_isoquant/data/`: Mouse.ONT_simulated.polyA,
.TSS, SIRVs.simulated_perfect_polya, Human.ONT_simulated.{polyA_1,polyA_2,TSS}. Each:
`run_type: transcripts,<polya|tss>_prediction`, `reduced_db` = alt-end split (built by
`misc/prepare_simulated_reduced_db.py`; excluded = `<base>_<coord>` alt-end variants = the
"novel" bucket). The refinement is now DEFAULT, so the configs carry NO extra flag.

**New de-novo CI test (added 2026-06-25): `Human.ONT_simulated.polyA_2.no_gtf`**
- Config `/abga/work/andreyp/ci_isoquant/data/Human.ONT_simulated.polyA_2.no_gtf.yaml`:
  no `genedb` (de-novo discovery), `isoquant_options: "-t 16 --force "`, `run_type:
  transcripts`, `reduced_db` kept only for scoring.
- Workflow `.github/workflows/Human.ONT_simulated.polyA_2.no_gtf.yml` (mirrors the Mouse
  no_gtf workflows; PREPENDS the CI bin so the gffcompare fork is used for td50/td10).
- Baselines: only `full_*` is meaningful de-novo (vs the complete expressed reference) —
  full 47.9/95.3 (def), 45.2/90.0 (td50), 37.5/74.6 (td10); `known_*: -1.0` (empty query →
  parse returns -1, sentinel matches); `novel_*` omitted (for a no-genedb run the "novel"
  split scores all discovered models against only the held-out subset → precision is a
  scoring artifact, so it is intentionally not gated).

Run + score (NO refinement flag — it is default):
```
PATH=/abga/work/andreyp/ci_isoquant/bin:$PATH   # gffcompare fork FIRST
python3 isoquant_tests/github/run_pipeline.py -o <outdir> <config>.yaml
```
3-delta scoring (default/td50/td10) via `misc/reduced_db_gffcompare.py` + the gffcompare
fork (`~/bin/gffcompare` IS the fork now; stock at `~/bin/gffcompare.stock.bak`). Stats
land in `<outdir>/<name>/gffcompare/isoquant.{full,known,novel}{,.td50,.td10}.stats`; parse
"Transcript level" Sn/Sp. A metric is checked only if its key is present in the config's
`baselines` block (±1% band; `-1.0` requires exact match). See `.claude/GFFCOMPARE.md`.

## `requires_polya_for_construction` interaction
auto = (polyA fraction ≥ 0.7). Gates FL-path formation on polyA presence
(`IntronPathStorage.fill`). When True, FL paths need a polyA vertex → alt-polyA NICs are
inherently polyA-confirmed. Force with `--polya_requirement required`. Open lever for the
weak Mouse polyA result; the real polyA knob is the assigner's FSM end-tolerance
(`long_read_assigner.py`), not the graph.

## Branch / parking state
- `transcript_model_ends` — the real branch. HEAD **c5950439** is local-only; origin is at
  **518ae62e** (decide whether to push c5950439).
- `tme_no_postpass` (1652e98f) — throwaway Part-3 diagnostic branch; keep local, do not merge.
- `goal2_end_refinement_wip` (4794f76) — parked early/standalone post-pass; superseded.
