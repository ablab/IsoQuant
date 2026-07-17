# gffcompare internals + IsoQuant transcript-discovery testing

How IsoQuant's transcript-model accuracy is scored with `gffcompare`, what
that score does and does **not** measure, and a local fork that adds an
independent terminal-end tolerance so end accuracy becomes visible.

All claims below were verified against the gffcompare v0.12.6 source + gclib.
The `--terminal-delta` fork now lives in its own repo
(**github.com/andrewprzh/gffcompare**, `master`) at **v0.12.10** — same gating
logic, just at shifted line numbers (`findRefMatch` refactored, `cds-match`
added upstream); clone it to `~/gffcompare` and `make release`.

## TL;DR

- gffcompare's **Transcript level** and **Intron-chain level** Sn/Sp are
  **end-agnostic for multi-exon transcripts**: a query matches (`=`) iff its
  *intron chain* is identical to a reference. The 5′ start (TSS) and 3′ end
  (polyA) — the coordinates polyA/TSS prediction refines — are ignored.
- So the standard reduced-db benchmark can confirm an end-refinement does
  **no harm**, but structurally **cannot show it helping**.
- The fork adds `--terminal-delta <d>`: enforce a terminal tolerance `d` at
  the **transcript level only** (independent of `-e`/exon-level), making
  Transcript-level Sn/Sp end-sensitive while Intron-chain stays end-agnostic.
- With it, the polyA/TSS end refinement shows a real benefit on novel
  transcripts at tight `d` (see results below).

## gffcompare matching internals (verified)

### Class codes that count as a match
- `=` : full match. For **multi-exon**: identical intron chain *and* both
  terminal ends within `terminalMatchRange`. For **single-exon**: boundaries
  within range.
- `~` : intron-chain match for multi-exon **with terminal ends differing**
  beyond `terminalMatchRange` (or ≥80% overlap for single-exon).

### `transcriptMatch` (gclib `gff.cpp`, the live matcher)
`tMatch` in `gtf_tracking.cpp` is dead code (inside `/* */`); the real
function is `transcriptMatch(a, b, ovlen, trange, cdsMatch)`. For multi-exon:
```cpp
for (int i=1;i<=imax;i++)                       // INTRON CHAIN only
  if (a.exons[i-1]->end != b.exons[i-1]->end || // donor
      a.exons[i]->start != b.exons[i]->start)   // acceptor
      return 0;                                  // intron mismatch
//--- full intron chain match
if (abs(a.exons[0]->start   - b.exons[0]->start)   <= trange &&  // 5' TSS
    abs(a.exons.Last()->end - b.exons.Last()->end) <= trange)    // 3' polyA
     result = '=';
else result = '~';
```
The outer terminal coordinates (`exons[0]->start`, `exons.Last()->end`) are
**never** part of the chain comparison — they only choose `=` vs `~`.
Single-exon goes through `singleExonTMatch` (overlap-based, ≥80%), so for
mono-exon transcripts the ends *do* affect the match.

### The `~`→`=` promotion (why default is end-agnostic)
`terminalMatchRange` defaults to **0** (`gtf_tracking.cpp`), and is raised to
`exonEndRange` (= the `-e` value, default 100) **only** with `--strict-match`
(`gffcompare.cpp` ~line 345). But in non-strict mode `~` is reported as `=`:
- `gffcompare.cpp` ~836 (stats tally): `if (!stricterMatching) tmatch='=';`
- `gffcompare.cpp` ~1795 (class-code output): `if (eqcode=='~' && !stricterMatching) { eqcode='='; }`

So in **default mode** any intron-chain match → `=` regardless of ends.

### How the levels are tallied (`gffcompare.cpp` ~836–856)
```cpp
char tmatch = transcriptMatch(q, r, ovlen, terminalMatchRange);
if (tmatch) {
  if (!stricterMatching) tmatch='=';        // promotion
  if (q->exons.Count()>1) super->ichainTP++; // Intron-chain TP: ANY match
  if (tmatch=='=')        super->mrnaTP++;    // Transcript TP: only '='
}
```
- **Intron-chain level** Sn/Sp ← `ichainTP` (counts `=` *and* `~`) → always
  end-agnostic.
- **Transcript level** Sn/Sp ← `mrnaTP` (only `=`) → end-agnostic *only
  because* of the promotion; end-sensitive once promotion is disabled.

### Stats levels in the `.stats` file
`Base` (nucleotide overlap — end-sensitive, but dominated by internal exons),
`Exon` (boundaries, with the `-e` tolerance), `Intron`, `Intron chain`
(multi-exon chain, end-agnostic), `Transcript` (≈ intron chain for multi-exon
+ overlap for mono-exon), `Locus`. "Sensitivity" column = recall (vs
reference), "Precision" = precision (vs query). `-e` (default 100) is, per the
code comment, "only used for exon level Sn/Sp" — except it also seeds
`terminalMatchRange` under `--strict-match`.

## The fork: `--terminal-delta`

Canonical fork: **github.com/andrewprzh/gffcompare** (`master`, commit
`05b3b1b`), gffcompare **v0.12.10** with a bundled `./gclib` copy. Clone to
`~/gffcompare`, build `make release` (GCLIB defaults to `./gclib`; needs `g++`,
`make`); the binary is installed at `~/bin/gffcompare`. Default behavior (no
`--terminal-delta`) is byte-identical to stock; the gate was re-verified on
v0.12.10 (20 bp end diff → default `=`/100, `--terminal-delta=10` → 0,
`=50` → 100). (Earlier local v0.12.6 build was at `gffcompare_term_delta/`.)

`--terminal-delta <d>` enforces terminal tolerance `d` at the **transcript
level only**, independently of `-e` (which stays tied to exon-level Sn/Sp),
and **without** the other `--strict-match` effects. Implementation:
- `gtf_tracking.cpp` / `.h`: new global `bool terminalStrict` (mirrors
  `stricterMatching`).
- `gffcompare.cpp`: register `terminal-delta=` long option; parse it after
  the `-e` block → `terminalMatchRange = d; terminalStrict = true;`.
- Gate **both** promotion sites with `&& !terminalStrict`
  (the ~836 tally and the ~1795 class-code line). So `~` is no longer
  promoted to `=`: Transcript-level (`mrnaTP`) becomes end-sensitive at `d`,
  Intron-chain (`ichainTP`) stays chain-based. `-e` / exon-level untouched.
- Help text added under `--strict-match`.

Behavior check (PB reduced-db full split): default `93.6` → `--terminal-delta
50` `93.5` → `--terminal-delta 10` `92.4` (Transcript Sn), while Intron-chain
stays `~93.9`.

> This is a **local** fork for measurement. Pushing it to a real GitHub fork
> / upstreaming the option is a separate step if wanted.

## How IsoQuant scores transcript discovery

Pipeline: `isoquant_tests/github/run_pipeline.py :: run_transcript_quality`
→ `misc/reduced_db_gffcompare.py` → `gffcompare` → parse `Transcript level`
Sn/Sp → compare to YAML `baselines.transcripts` within ±1%
(`check_value`, percent=0.01).

`reduced_db_gffcompare.py`:
1. Splits the IsoQuant output GTF (`<label>.transcript_models.gtf`) into
   `full` / `known` / `novel` via `SEPARATE_FUNCTORS["isoquant"]`
   (by transcript_type: known vs novel_in_catalog/novel_not_in_catalog).
2. Runs `gffcompare -r <ref> -o <out> <query>` (via `misc/common.py
   run_gff_compare`, **no** strict flag) for each split against a reference
   subset built from the `reduced_db` prefix:
   - `full`  vs `<prefix>.expressed.gtf`       (all expressed transcripts)
   - `known` vs `<prefix>.expressed_kept.gtf`  (kept = still in the reduced db)
   - `novel` vs `<prefix>.excluded.gtf`        (removed → must be rediscovered)
3. `parse_gffcomapre` reads the `Transcript level` line; the launcher checks
   `{full,known,novel}_{recall,precision}`.

Reduced-db data prefix (mouse): `mouse.gencode.M26.spatial.20percent` with
`.reduced.gtf` (what IsoQuant is given), `.expressed.gtf`, `.expressed_kept.gtf`,
`.excluded.gtf`, `.expressed.tsv` counts. CI configs:
`/abga/work/andreyp/ci_isoquant/data/Mouse.*.reduced_db.yaml` (run via the
self-hosted-runner workflows in `.github/workflows/Mouse.*.yml`). Baselines
live in the YAML `baselines:` block; **master is the actual baseline** — the
YAML numbers can be stale/rounded (e.g. the polyA config's were).

Key consequence for polyA/TSS work: because `full`/`known`/`novel` use the
default end-agnostic transcript metric, **moving a multi-exon transcript's
TSS/polyA never changes its score** there. Only mono-exon transcripts, the
`Base`/`Exon` levels, and the `--terminal-delta` runs are end-sensitive.

## End-sensitivity result (does the polyA/TSS refinement help?)

Comparing **master** (clustering ends) vs the **reuse** branch (detector-
refined ends) at three terminal deltas. Transcript-level Sn/Sp:

| dataset | split | metric | default | δ=50 | δ=10 |
| ------- | ----- | ------ | ------- | ---- | ---- |
| PB  | novel | master | 76.2/94.3 | 76.1/94.2 | 67.2/83.2 |
| PB  | novel | reuse  | 76.1/94.1 | 76.0/94.0 | **68.4/84.6** |
| ONT | novel | master | 61.1/85.9 | 57.6/81.0 | 45.4/63.9 |
| ONT | novel | reuse  | 61.2/85.2 | 57.7/80.3 | **46.5/64.7** |
| PB  | full  | master | 93.7/98.9 | 93.6/98.8 | 92.2/97.3 |
| PB  | full  | reuse  | 93.6/98.8 | 93.5/98.7 | **92.4/97.5** |
| ONT | full  | master | 84.7/96.6 | 84.1/95.9 | 82.3/93.8 |
| ONT | full  | reuse  | 84.6/96.5 | 84.0/95.8 | 82.3/93.9 |

Reading:
- At **default** (end-agnostic) master ≈ reuse everywhere (the standard CI
  view: refinement is neutral / does no harm).
- At **δ=10** (ends must land within 10 bp) the refinement **wins on novel
  transcripts**: PB +1.2 recall / +1.4 precision, ONT +1.1 / +0.8 — i.e. the
  detector-refined ends are closer to the simulated truth than clustering
  ends. (`full` shows little gap because it's dominated by *known*
  transcripts, which are emitted with annotated ends in both branches.)
- This is the benefit the default metric structurally hides; the divergence
  at δ=10 also confirms `output/master` is a genuine no-refinement baseline
  (identical at default, different once ends are scored).

## Wired into the assessment pipeline (live)

The 3-delta scoring is now part of transcript-discovery assessment:
- `misc/reduced_db_gffcompare.py` runs each split (full/known/novel) at
  `TERMINAL_DELTAS = [None, 50, 10]` → `isoquant.<split>.stats` (default) +
  `isoquant.<split>.td50.stats` + `isoquant.<split>.td10.stats`. The default
  run uses stock gffcompare semantics; the deltas pass `--terminal-delta=N`
  (needs the fork). `misc/common.py:run_gff_compare` is unchanged (its
  `additional_option` slot carries the `--terminal-delta=N` token).
- `isoquant_tests/github/run_pipeline.py:run_transcript_quality` parses all
  three stats variants and checks three baseline blocks —
  `transcripts` (default), `transcripts_td50`, `transcripts_td10` — each with
  `{full,known,novel}_{recall,precision}`. A variant is checked only if its
  baseline block exists *and* its stats file was produced, so a stock
  gffcompare (no `--terminal-delta`) degrades cleanly to the default block.
- The CI YAML configs (`/abga/work/andreyp/ci_isoquant/data/Mouse.*.yaml`,
  7 of them) now carry the three baseline blocks, computed from the
  **master-branch** output GTFs via `local_tmp/update_baselines.py`
  (runs the modified `reduced_db_gffcompare.py` on master's
  `transcript_models.gtf` and reads the Transcript-level Sn/Sp). The fork is
  installed at `/abga/work/andreyp/ci_isoquant/bin/gffcompare` (stock backed
  up as `gffcompare.stock`). The fork now lives in its own repo
  **github.com/andrewprzh/gffcompare** (`master`); the local
  `misc/gffcompare_fork/` patch + README were removed once upstreamed there.

Validation: against master output all three blocks pass exactly (exit 0);
against the refined (reuse) branch the **td10 novel** metrics flag *higher*
than the master baseline (PB novel_recall 68.4 > 67.2, novel_precision
84.6 > 83.2) while default/td50 stay in band — i.e. the pipeline now
*detects* the polyA/TSS refinement's benefit instead of being blind to it.

## Gotchas / notes

- `--terminal-delta`'s tolerance is independent of `-e`; `--strict-match`'s
  is not (it reuses `-e`, default 100 — too loose for polyA/TSS scale).
- The fork builds against gclib **master** (the Makefile clones it); a
  future gclib change to `transcriptMatch` could drift. Pin gclib if exact
  reproducibility matters.
- `run_gff_compare(reference, compared, output, additional_option="")` already
  has an `additional_option` slot — to score the IsoQuant CI end-sensitively,
  point it at the fork binary and pass `--terminal-delta <d>`.

## Alt-end transcript-discovery benchmarks (simulated TSS/polyA)

The terminal-delta scoring shines on the simulated **alternative TSS/polyA**
datasets (configs `Human.ONT_simulated.polyA_1`/`polyA_2`/`TSS`,
`Mouse.ONT_simulated.TSS`, + `SIRVs.simulated_perfect_polya`; previously these
had no workflow and ran only quantification/prediction with model construction
disabled — then via `--no_model_construction`, now expressed by omitting
`transcript_discovery` from `--analysis`; the flag is deprecated but still works).

The simulation GTF (`reference_polya_gtf`/`reference_tss_gtf`,
`coroari/thesis_model/.../updated_*.gtf`) is full GENCODE (plain IDs) **plus**
alt-end variants `<base>_<coord>` that share the base intron chain but have a
different terminal exon. IsoQuant runs with the **normal** GENCODE reference;
scoring its models against the simulation annotation at tight `--terminal-delta`
measures alt-end recovery.

Wiring:
- `misc/prepare_simulated_reduced_db.py` splits the simulation GTF (filtered to
  expressed via `simulated_counts`) into the reduced-db trio: plain IDs →
  `expressed_kept` (known), alt-end `_<digits>` IDs → `excluded` (novel). It
  normalizes the 10-column lrgasp human GTFs to 9 columns. Output:
  `/abga/work/andreyp/ci_isoquant/data/simulated_truth/<name>.*`.
- The four ONT configs: dropped `--no_model_construction` (i.e. transcript
  discovery now on — equivalently `--analysis` with `transcript_discovery`
  included), `run_type: transcripts,<polya|tss>_prediction`, `reduced_db` = the
  split prefix. So one run yields both the discovery assessment (vs the
  simulation, 3 deltas, via `reduced_db_gffcompare.py`) and the prediction TSV
  assessment.
- Five `.github/workflows/<config>.yml` (thin-wrapper template, staggered cron).
- Baselines from a master (cf95528) run via a worktree + `local_tmp/run_master_isoquant.py`
  + `local_tmp/set_sim_baselines.py` (runs current-branch 3-delta
  reduced_db_gffcompare + assess_polya_prediction, writes the config baselines).
  Master outputs in `output/master_simbench/<name>/`.

Master baseline reading (the motivating signal): master collapses alt-end
transcripts to the **known** base (or drops them) — it never emits them as NIC —
so the `novel` split is **0/0** and `full_recall` falls as the delta tightens
(Mouse TSS full_recall 70.1 → 62.6 → 59.3 at default/50/10). The NIC-routing +
end-refinement branch is expected to lift `novel_td10` above 0 and `full_td10`
above master when its workflow runs (the same end-sensitive check that flagged
+1.2/+1.4 on reduced_db PB novel).
