# Cell Type Decomposition (ILP) Branch

Notes on the `cell_type_decomposition` branch: an experimental transcript-model
construction path that replaces the heuristic graph walk with an **Integer
Linear Programming (ILP) flow-decomposition** whose edge/node weights are
**vectors indexed by a cell-type tree** rather than scalars.

Branch diverged from `master` at `159f9af`. Built on top of the earlier ILP
prototype on `ilp_models` (shared merge base `de74f33`) and continued the
`ilp_model` line (shared base `4b42880`).

## 1. Key Features Implemented on This Branch

Organized by functional block; commit hashes are inline for reference.

### 1.1 Cell-type hierarchy (`src/cell_type_tree.py`, `2c0ea67`)

- `CellTypeTree` takes a set of colon-separated cell-type names
  (e.g. `"G1:G1.1:G1.1.1"`) and builds a tree rooted at `ROOT`.
- Assigns each node (ROOT, internal and leaf types) a stable integer index,
  used to address positions in per-edge/per-node weight vectors.
- `transform_counts(defaultdict[cell_type -> count])` → `np.ndarray` of length
  `n_of_cell_types` with counts accumulated **up the tree** (count at a leaf
  also contributes to all ancestors).
- `ignore_missing` / `remove_missing`: policy for reads with cell type `"NA"`
  (commits `f74265e`, `769f439`, `459b2e3`, `f5a5bf8`).
- Helpers: `get_leaf_types`, `get_cell_type_index`,
  `transform_counts_to_arrays` (by edge), `transform_count_to_arrays`
  (indexed), `transfrom_count_array_to_dict` (back to a leaf-only dict for
  output).

### 1.2 Cell-type-aware intron graph (`src/intron_graph.py`, `e259b5d`)

- `IntronCollector` gained `clustered_introns_by_cell_type:
  defaultdict[intron -> defaultdict[cell_type -> int]]` alongside the scalar
  `clustered_introns`. `collect_introns` reads `assignment.read_group` and
  bumps the per-cell-type counter (`e259b5d`, `d81dcfb`, `459b2e3`).
- Intron clustering, intron substitution (`add_substitute`) and `discard`
  keep both scalar and per-cell-type maps consistent.
- `add_edge(v1, v2, read_group)` stores a per-read-group edge weight in
  `edge_weights[(v1, v2)]: defaultdict[read_group -> int]`
  (`d81dcfb`, `76e7f25`).
- `IntronGraph.split_to_connected_subgraphs()` returns a list of
  `IntronSubGraph` views — one ILP instance is solved per connected component
  (`ce30d19`). Each subgraph carries its own `clustered_introns_by_cell_type`.
- Pre-simplification / pre-error correction before ILP: remove low-weight
  vertices / correct obvious errors (`306ec0a`, `f5a5bf8`, `00b0992`,
  `259e16d`).

### 1.3 ILP entry point (`src/ilp_model.py`)

See also §3 below for the graph-conversion layer.

Main function: `ILP_Solver_Nodes(intron_graph, chr_id, gene_id, index,
transcripts_constraints, ground_truth_isoforms, epsilon, timeout, threads,
draw_graphs)`:

1. Converts the IsoQuant intron graph to a `networkx.DiGraph` with node
   attributes `flow` (scalar) and `cell_flow` (np.array over the cell-type
   tree) — via `Intron2Nx_Node` (§3).
2. Converts IsoQuant path constraints (from full-length reads) into edge
   lists with `Constraints_Transfer_Format`.
3. Runs `MinErrorCellTypeFlowCorrection` to **correct the node weights so
   flow conservation holds** (cell-type-aware; replaces `MinErrorFlow` from
   `flowpaths`).
4. Runs `MinFlowCellDecomp` on the corrected graph to obtain an MFD whose
   paths carry both total weight and a per-cell-type vector of weights.
5. Returns `list[(path, total_weight, {cell_type: count})]` for the caller.

Timeout (`SIGALRM`) handling wraps both solver phases (`abcdfa1`, `1f8d47d`,
`c941f61`). On timeout returns `[]` and optionally dumps the failed graph
(`draw_graphs=True`).

### 1.4 ILP models (`src/models/*.py`)

All derive from `flowpaths` base classes but extend them with a second,
vector-valued flow attribute.

- `NodeExpandedDiGraph` (`nodeexpandeddigraph.py`) — local fork of
  `flowpaths.nodeexpandeddigraph` that also understands
  `node_cell_flow_attr` (in addition to `node_flow_attr`), propagates the
  vector through the node-split transformation, and supports
  `additional_starts` / `additional_ends` passthrough. Needed because the
  vanilla flowpaths expansion only knows scalar flows.
- `MinErrorCellTypeFlowCorrection` (`cellflowcorrection.py`) — inherits
  `fp.MinErrorFlow`. Per-edge variables indexed by
  `(u, v, cell_type)`; flow conservation is enforced **per cell type**;
  objective is the total L1 error between the solver-chosen cell-type flow
  and the observed `cell_flow` vector on each edge. Tree balancing was added
  (`a2df44b`): parent-cell-type counts must equal the sum of child counts,
  enforced as extra constraints across the hierarchy.
- `MinErrorCellTypeFlow` (`minerror_celltypeflow.py`) — simpler variant,
  kept for reference/comparison, without the balancing across hierarchy.
- `kFlowCellTypeDecomp` (`kflowcellcomp.py`) — fixed-k path-decomposition
  ILP. Extends `fp.AbstractPathModelDAG`. Each path carries
  `path_weights_ct[cell_type]` variables with flow conservation per cell
  type on the selected edges. Symmetry-breaking constraints added (`528b63b`,
  `12e7f2e`).
- `MinFlowCellDecomp` (`minflowcelldcomp.py`) — iterates k upward calling
  `kFlowCellTypeDecomp`, analogous to `fp.MinFlowDecomp`. Integrates:
  lower-bound generation via min-gen-set, safe-paths, subpath constraints
  from reads, symmetry breaking, given-weights warm start.
- `*approximate` variants (`0993cfa`, `4fc915d`) — same encoding with a
  relaxed/approximate objective (continuous variables or loosened
  flow-conservation), produced for use when the exact model times out.

### 1.5 Pipeline integration

- `src/graph_based_model_construction.py`:
  - Imports `ILP_Solver_Nodes` (line 36).
  - `GraphBasedModelConstructor.__init__` takes
    `ground_truth_gene_info` and stores `self.ilp_solution_assignement`.
  - `process()` (line 146): if `params.no_ilp` is false, calls
    `construct_ilp_isoforms()` instead of the classic
    `construct_fl_isoforms` + assignment-based path.
  - `transfer_fl_path` collects full-length read intron chains as ILP
    subpath constraints (filtering terminal polyA/polyT vertices).
  - `construct_ilp_isoforms` (line 636): splits the intron graph into
    connected subgraphs, runs the ILP on each, and for each returned path
    constructs a `TranscriptModel` (novel-in-catalog / NNIC / matching
    reference isoform) with weight and per-cell-type weights attached.
  - `ilp_counts()` forwards `(transcript_id, gene_id, weight, ct_weights)`
    tuples to the `TranscriptCounter` via `add_ilp_data`.

- `src/long_read_counter.py`:
  - New `ILPFlowCounter` class handles weighted (non-integer) counts coming
    from the ILP, instead of per-read accumulation.
  - `CompositeCounter.add_ilp_data` and `AbstractCounter.add_ilp_data`
    fan out to internal counters.
  - Factory `create_ilp_transcript_counter` wires the counter into the
    pipeline.
  - `ILPFlowCounter.add_ilp_data` reads the `rg_weights` dict (a leaf-only
    projection of the cell-type vector) and routes weights to read groups.

- `src/dataset_processor.py`:
  - Builds `ground_truth_gene_info` from `--ground_truth_genedb` and passes
    it to `GraphBasedModelConstructor` (`73101c6`, `619af36`, `82685b9`).
  - Uses `create_ilp_transcript_counter` when ILP is enabled.

- `isoquant.py`:
  - `--no_ilp` flag (default False, i.e. ILP is ON when the model-strategy
    allows).
  - New `ilp_model` entry in `--model_construction_strategy` choices (line
    176) with its own `ModelConstructionStrategy` tuple (line 645).
  - `--ground_truth_genedb` reference argument.

### 1.6 Debugging / output

- Optional per-subgraph PNG dumps via `fp.utils.draw` — `flow` and `cell_flow`
  variants for input, corrected and solved graphs
  (`f8e62fd`, `af3f476`, `4b42880`, `7b6d30d`).
- `export_data(...)` in `ilp_model.py` pickles graph, additional
  starts/ends, ignored edges and path constraints to disk for offline
  analysis (toggled by the `export` local flag).
- Ground-truth isoform drawing (`4b42880`) and `ground_truth_isoforms` pass
  into the ILP for reference/debug.

### 1.7 Time/approximation knobs

- `timeout` hard-limit wrapping both solver phases (`1f8d47d`, `abcdfa1`,
  `c941f61`).
- Approximate models (`*approximate.py`) for bail-out when exact solve is
  too slow (`0993cfa`, `4fc915d`).
- Symmetry breaking in `kFCD` (`528b63b`, `12e7f2e`, `1df75cc` reverts it in
  parts).
- `--filtering_missing` / `--skip_NAs` options propagating into
  `CellTypeTree` `ignore_missing` / `remove_missing` (`769f439`, `f74265e`).

## 1.8 How to Provide Cell Types

Cell-type labels ride on the **existing `--read_group`** mechanism — there is
no separate cell-type flag. Every read's `read_group` string is consumed by
`IntronCollector.collect_introns` and fed to `CellTypeTree` as a hierarchy.

### Hierarchy encoding

`CellTypeTree._add_cell_type` splits on `:`, so any label that uses colons
becomes a tree path rooted at `ROOT`:

```
"Immune:T:CD8"   → ROOT → Immune → T → CD8
"Immune:T:CD4"   → ROOT → Immune → T → CD4
"Immune:B"       → ROOT → Immune → B
"Epithelial"     → ROOT → Epithelial
```

Leaves are the concrete cell types; internal nodes get their counts by
tree-sum aggregation (`transform_counts` walks to `ROOT`). Granularity is
chosen by how deep the label is.

### Supplying the labels (standard `--read_group` forms)

1. **TSV file** (most common for cell typing):
   ```
   --read_group file:cell_types.tsv:0:1:tab
   ```
   columns: `read_id <tab> Immune:T:CD8`. Example format at
   `tests/simple_data/chr9.4M.ont.sim.read_groups.tsv`.

2. **BAM tag** (pre-annotated single-cell BAM):
   ```
   --read_group tag:CT
   ```

3. **Parsed from read id**:
   ```
   --read_group read_id:_
   ```

### Missing values

Reads labelled `NA` are filtered — `src/ilp_model.py:84` hard-codes
`CellTypeTree(cell_types, ignore_missing=True, remove_missing=True)`. No CLI
flag; change it there to roll the `NA` bucket into `ROOT` instead.

### Caveats on this branch

- `IntronCollector.collect_introns` uses `assignment.read_group` as a dict
  key → `read_group` **must be a single string**. Do not combine grouping
  strategies here (master's multi-group `list[str]` API is one of the
  merge-conflict points; see §4.1 bullet 3).
- Only one grouping strategy is supported end-to-end: the label *is* the
  cell type.
- `ILPFlowCounter.add_ilp_data` emits one weight per leaf cell type per
  transcript → output granularity matches the deepest level of the label.

## 2. Documentation of Graph Conversion (Intron Graph → ILP)

Located in `src/ilp_model.py`. Two conversion concerns: the **graph itself**
and the **subpath constraints**.

### 2.1 Graph: `Intron2Nx_Node(intron_graph, skip_isolated_nodes=True, skip_terminal_nodes=False)`

Pipeline (lines 72–131):

1. Walk `intron_graph.intron_collector.clustered_introns_by_cell_type` to
   collect the union of observed cell types, then build a `CellTypeTree`
   (honouring `ignore_missing` / `remove_missing`).
2. For each clustered intron (an internal vertex):
   - Key = `str(intron)` (a string "(type, pos)" form; stringified so the
     solver/NetworkX accept hashable tuples and so `fp.utils.draw` plays
     well).
   - Node attribute `flow` = scalar count.
   - Node attribute `cell_flow` = `cell_tree.transform_counts(per_cell_type_counts)`.
3. Add any terminal (polyA / polyT / read_start / read_end) vertices that
   appear in incoming/outgoing edges but not in the clustered map (they have
   no `flow`/`cell_flow` attributes — caller must mark them as
   `edges_to_ignore`).
4. Add edges from `incoming_edges` and `outgoing_edges`:
   - Terminal-vertex edges go into `additional_starts` (polyT/read_start
     side) or `additional_ends` (polyA/read_end side).
   - If `skip_terminal_nodes=True`: terminal vertices are excluded from
     the graph and their neighbours recorded as additional start/end.
   - Otherwise the terminal vertices are added AND the edges touching them
     are also recorded in `edges_to_ignore` so flow error is not measured on
     them.
5. `clean_graph` removes isolated nodes when requested.
6. Returns `(G, cell_type_tree, additional_starts, additional_ends,
   edges_to_ignore)`.

Note: nodes are **strings**. The inverse mapping (string → intron tuple) is
done by `eval`-free parsing downstream inside `graph_based_model_construction.
construct_ilp_isoforms` when the solver returns paths (via `transfer_paths`,
which is the symmetric reconstruction step).

### 2.2 Path constraints: `Constraints_Transfer_Format(input_constraints, skip_isolated_nodes, skip_terminal_nodes)`

- Input: list of intron-vertex lists (`path_storage.fl_paths` collected by
  `GraphBasedModelConstructor.transfer_fl_path`).
- Output: list of lists of edges `(str(v_i), str(v_{i+1}))` — the
  solver-native format (`subpath_constraints` for `MinFlowCellDecomp`).
- Optionally drops terminal vertices (polyA/polyT/read_start/read_end) from
  each constraint.
- Empty constraints (length ≤ 1) are discarded.

### 2.3 Post-conversion pruning (in `ILP_Solver_Nodes`)

After `Intron2Nx_Node`, additional starts/ends and subpath constraints are
*filtered against the actual graph nodes/edges* (lines 243–257 of
`ilp_model.py`) — because `clean_graph` and `edges_to_ignore` may have
dropped nodes referenced in the original lists. This is the
`additional_starts_pruned` / `additional_ends_pruned` /
`subpath_constaints_pruned` that is actually passed to the solver.

### 2.4 Node expansion (solver-side, `src/models/nodeexpandeddigraph.py`)

The solver operates on edge-weighted graphs. Because our weights live on
**nodes**, `NodeExpandedDiGraph` splits each node `v` into `v_in → v_out`,
carrying the `flow` / `cell_flow` attributes onto the internal edge and
adding the incoming/outgoing edges of `v` to `edges_to_ignore`. This is the
mechanism that bridges IsoQuant's node-weighted intron graph to the
flowpaths edge-weighted ILP formulation.

## 3. ILP-Related Branches — What's Unique

Merge bases:

| pair                                          | base      |
|-----------------------------------------------|-----------|
| `cell_type_decomposition` ∩ `ilp_models`      | `de74f33` |
| `cell_type_decomposition` ∩ `ilp_model`       | `4b42880` |
| `ilp_models` ∩ `ilp_models_ref`               | `702caca` |
| `cell_type_decomposition` ∩ `ilp_data_logger` | `e926433` |

### `origin/ilp_models` (original prototype — Fernando-style)

41 commits not on `cell_type_decomposition` (from the shared base
`de74f33` onward). These are the **early ILP foundation** by a different
author: graph conversion v1 (`ea3762d`), first MFD encoding (`4fa5034`),
or-tools / gurobi solver experiments, safety-path optimizations, subpath
constraint encoding, `fl_paths` as path constraints, visualization
skeleton, weighted-transcript construction. **None of these are on
`cell_type_decomposition`** — this branch only inherited the shared pre-
`de74f33` history.

### `origin/ilp_models_ref` (refinement of `ilp_models`)

`ilp_models` + 10 refactoring commits. Highlights: returns weighted
transcripts from the ILP module (`861450d`), full path encoding
(`c83b0a6`), path-constraint filtering (`cb35769`), readability/print
improvements, note on the `M` parameter of the robust encoding vs. max
flow. Still **no cell-type concept** — all scalar.

### `origin/ilp_model` (singular — the debugging line)

~40 unique commits that are **not** on `cell_type_decomposition`. Mostly
debug prints, many empty-constraint bugfix attempts, plot toggles,
report_unstranded=False workaround, model-logic fixes, upper edge-count
limit, symmetry-breaking experiment. No cell-type code; this is the
branch where the non-cell-type ILP implementation was stabilised.

### `origin/ilp_data_logger`

`cell_type_decomposition` as of `e926433` (Feb-ish snapshot) + 2 commits:
`1050777 added data collection for ILP-model` and `19e8587 Minor fix`. The
data-logging add-on itself (presumably per-gene timing / solve stats) is
not on `cell_type_decomposition`. Everything else — all the recent
approximate models, time limits, k-FCD objective changes, disjoint
subgraph splitting, pre-error correction, NA handling — is **only** on
`cell_type_decomposition` and is NOT on `ilp_data_logger`.

### Short answer

- **`cell_type_decomposition` is the canonical trunk** for the cell-type
  variant. It contains all the ILP+cell-type work, plus all the recent
  approximate/time-limit/NA-handling polish.
- `ilp_model` is a **parallel scalar-ILP** branch with its own bugfix
  history; its commits are not subsumed by `cell_type_decomposition`.
  Useful to consult for scalar-ILP references or stabilisation tricks,
  but cannot be fast-forward merged.
- `ilp_models` / `ilp_models_ref` are an **older, independent prototype**
  line. Its early commits were shared (pre-`de74f33`) and its 40+ later
  commits are **not** on `cell_type_decomposition`. They would have to
  be cherry-picked if any of that work is wanted.
- `ilp_data_logger` is a **stale fork** of `cell_type_decomposition` (~19
  commits behind) with a 2-commit data-collection add-on. The data
  logging is worth cherry-picking if we want it back.

## 4. Differences vs. Up-to-date Master and Main Merge Conflict Points

`cell_type_decomposition` diverged from master at `159f9af`. Aggregate
diff: **241 files changed, +6909 / -26868**. The *big* reason for the
negative delta is a massive refactor on master that this branch has not
seen.

### 4.1 Structural refactors on master since divergence

These are the hardest merge concerns; they are structural, not local.

1. **`src/` → `isoquant_lib/`** — master **deleted the entire `src/`
   package** and moved the code to a new top-level package
   `isoquant_lib/`. Every file on this branch that lives under `src/`
   (including the new ILP files) will appear as an "add on deleted parent"
   conflict. Concretely:
   - `src/ilp_model.py`, `src/cell_type_tree.py`, `src/models/*.py` — pure
     additions; need to be copied under `isoquant_lib/` and their import
     paths adjusted (`from src.X` → `from isoquant_lib.X` or `from .X`).
   - `src/graph_based_model_construction.py`, `src/intron_graph.py`,
     `src/long_read_counter.py`, `src/dataset_processor.py`,
     `src/long_read_assigner.py`, `src/assignment_io.py`,
     `src/gene_info.py`, `src/isoform_assignment.py`, etc. — all edits
     must be replayed onto `isoquant_lib/*.py`, reconciling with
     master-side changes to the same files.

2. **`isoquant.py`** — 762-line delta (184 ins / 578 del). Command-line
   parsing was heavily restructured on master.
   Conflict points:
   - `--no_ilp` flag (this branch) — re-add.
   - `ilp_model` choice in `--model_construction_strategy` — re-add its
     tuple to the strategy dict.
   - `--ground_truth_genedb` argument — re-add (unless master added its
     own version).

3. **Multi-group refactor on master** (see `.claude/MULTI_GROUP_IMPLEMENTATION.md`).
   Master now represents `ReadAssignment.read_group` as `list[str]`
   (one entry per grouping strategy) and routes counts through a
   `CompositeCounter`. This branch still assumes `read_group: str`:
   - `IntronCollector.collect_introns` does
     `introns_with_cell_types[intron][read_group] += 1` — with a list
     this will break (`unhashable type: list`).
   - `ILPFlowCounter.add_ilp_data` indexes `rg_weights` by a single
     group key.
   - `GraphBasedModelConstructor.forward_counts` / `ilp_counts` call
     `add_read_info_raw(read_id, [transcript_id], read_group)` with a
     single group; master's CompositeCounter expects the list.
   - The choice of **which** group index represents the cell type must
     be made explicit (likely a new grouping strategy that extracts cell
     type from BAM tag or TSV).

4. **New isoquant packages on master** not on this branch:
   - `isoquant_lib/barcode_calling/` (all protocol callers, MDF format,
     universal extraction, two-bit indexing)
   - `isoquant_lib/parallel_workers.py`
   - `isoquant_lib/processed_read_manager.py`
   - `isoquant_lib/string_pools.py`
   - `isoquant_lib/table_splitter.py`
   - `isoquant_lib/assignment_aggregator.py`,
     `isoquant_lib/assignment_loader.py`
   - Relevant only in that imports and data-flow through
     `dataset_processor.py` will have moved.

5. **Test/CI rework**: this branch still has `tests/` plus a huge number
   of additional `.github/workflows/*.yml` that were **removed** on master
   (moved to `isoquant_tests/github/`). The `Barcode.*.yml`,
   `SC.Mouse.*.yml`, `Stereo_toy.yml`, `Short_runs.yml`,
   `SIRVs.Set4.R10.ubam.yml`, `claude-code-review.yml` etc. deletions all
   show up as diff. Most of these are cosmetic from the ILP work's point
   of view but will flood the merge conflict UI.

### 4.2 File-level conflict hotspots (ordered by severity)

| severity | file (branch path)                           | why                                                                                          |
|:--------:|----------------------------------------------|----------------------------------------------------------------------------------------------|
|  🔴 high | `src/` → `isoquant_lib/`                     | whole-directory rename; every modified file is a conflict                                    |
|  🔴 high | `graph_based_model_construction.py`          | ILP branch in `process()`, ILP subfunctions; master also edited this file                    |
|  🔴 high | `long_read_counter.py`                       | new `ILPFlowCounter`, `add_ilp_data` vs master's CompositeCounter multi-group API            |
|  🔴 high | `intron_graph.py`                            | `clustered_introns_by_cell_type`, `edge_weights` as dict, subgraph splitting, terminal vertex changes, `add_edge(read_group)` |
|  🟠 med  | `dataset_processor.py`                       | ground-truth plumbing, ILP counter wiring vs master's new worker architecture                |
|  🟠 med  | `isoquant.py`                                | flag/strategy additions on top of a heavily reworked parser                                  |
|  🟠 med  | `isoform_assignment.py` / `long_read_assigner.py` | read_group type change (str vs list[str]) — mostly data-flow, not logic                 |
|  🟡 low  | `requirements.txt`                           | add `flowpaths`, `networkx`, `gurobipy`; master may have added its own pins                  |
|  🟡 low  | tests (`tests/*.py`, `.github/workflows/*.yml`) | deletions on master; likely keep master's cleanup                                          |

### 4.3 Semantic incompatibilities (not syntactic)

These will compile after rename+merge but break at runtime:

- `read_group` being a `str` vs `list[str]` — see §4.1 bullet 3.
- `CompositeCounter` on master expects `add_read_info_raw` / `add_ilp_data`
  arguments shaped per the multi-group protocol. The ILP counter must
  either (a) pick one group index, or (b) be replicated per strategy like
  the existing counters.
- `TranscriptCounter` on master has changed how `add_confirmed_features`
  and unassigned reads are handled; `ilp_counts()` + `forward_counts()`
  on this branch will need adjustment.
- String-interning refactor (planned, see `.claude/STRING_INTERNING_OPTIMIZATION.md`):
  gene/transcript/cell-type IDs may become integer indices on master;
  `CellTypeTree` currently assumes strings. Not yet landed on master but
  worth tracking.

### 4.4 Recommended merge approach

1. Rebase is not viable (too much history across renamed paths). Prefer a
   **manual port**:
   - Check out `master`.
   - Copy `src/ilp_model.py`, `src/cell_type_tree.py`,
     `src/models/` verbatim under `isoquant_lib/` and fix imports.
   - Re-apply `graph_based_model_construction.py`, `intron_graph.py`,
     `long_read_counter.py`, `dataset_processor.py`, `isoquant.py`
     hunks by hand against the `isoquant_lib/*` versions.
   - Rewrite the read-group handling inline to honour the multi-group
     `list[str]` API (probably: add a new grouping strategy for cell
     type, and make the ILP counter address its own `group_index`).
2. Land as a **feature flag** (`--no_ilp` default True initially) so the
   merge doesn't regress the existing pipeline.
3. Regenerate requirement pins and CI workflow for the ILP path (the
   master-only CI has `ILP.PB.yml` — verify it still runs the right
   entry-point after the port).

## 5. Merge into master on branch `cell_type_decomposition_26`

Branch `cell_type_decomposition_26` is the in-progress merge of `master`
(through commit `44445fb "mention gtf2db.py"`) into
`cell_type_decomposition` (tip `8870527 "mv src"`). The merge commit is
`884e294 "master merge in progress"` (Merge: 8870527 44445fb).

### 5.1 What the merge commit did (user's work)

- **Directory rename resolved**: adopted master's `src/` → `isoquant_lib/`
  layout. All branch-specific files were moved under `isoquant_lib/`
  (including the ILP-only files: `ilp_model.py`, `cell_type_tree.py`,
  `models/*.py`).
- **New files kept from the ILP branch**:
  - `isoquant_lib/ilp_model.py`
  - `isoquant_lib/cell_type_tree.py`
  - `isoquant_lib/models/{cellflowcorrection,kflowcellcomp,kflowcellcompapproximate,minerror_celltypeflow,minflowcelldcomp,minflowcelldcompapproximate,nodeexpandeddigraph}.py`
- **CLI (`isoquant.py`)** re-applied against master's parser:
  - `--ground_truth_genedb` (gffutils DB of expressed isoforms for ILP eval).
  - `--no_ilp` (flag to force the default model-construction path;
    `action='store_true'`, default `False` — i.e. ILP *on* by default
    when strategy selects it).
  - New `ilp_model` strategy in `--model_construction_strategy` choices
    and in `set_model_construction_options()`'s `strategies` dict.
- **`isoquant_lib/parallel_workers.py`**:
  - Loads `gffutils.FeatureDB` from `args.ground_truth_genedb` once per
    chromosome and, per gene window, slices a `ground_truth_gene_info`
    to pass into `GraphBasedModelConstructor`.
  - Imports `GeneInfo` to build those ad-hoc infos.
- **`isoquant_lib/graph_based_model_construction.py`**:
  - New `ground_truth_gene_info=None` ctor kwarg + attrs
    (`self.ground_truth_isoforms`, `self.ilp_solution_assignement`).
  - `process()` now branches on `self.args.no_ilp` — default path
    preserved, ILP path calls `construct_ilp_isoforms()` + ILP-specific
    filtering, then `assign_reads_to_models()` twice (pre/post filter)
    and finally `ilp_counts()`.
  - Added `construct_ilp_isoforms()` (calls
    `split_to_connected_subgraphs()` → `ILP_Solver_Nodes` per subgraph),
    `forward_counts()`/`ilp_counts()` helpers.
- **`isoquant_lib/long_read_counter.py`**:
  - `NormalizationMethod` enum moved to top.
  - `AbstractCounter.add_ilp_data()` stub + `CompositeCounter.add_ilp_data()`
    that fans out to children.
  - New `ILPFlowCounter` (≈200 lines) producing linear-format grouped
    counts from the `list[(path, total_weight, {cell_type: count})]`
    returned by `ILP_Solver_Nodes`.
  - `create_ilp_transcript_counter()` factory.
- **`isoquant_lib/intron_graph.py`**:
  - `IntronCollector.clustered_introns_by_cell_type` + per-cell-type
    tracking in `collect_introns`, `cluster_introns`, `add_substitute`,
    `discard`.
  - `IntronGraph.edge_weights` changed to
    `defaultdict(defaultdict)` keyed on `(v1,v2) -> {read_group: count}`.
  - `add_edge(v1, v2, read_group)` + `_cell_type_key()` helper.
  - `error_correction()` method (drop edges with weight < 3),
    `split_to_connected_subgraphs()` + `IntronSubGraph` view class.
- **`isoquant_lib/assignment_aggregator.py`**:
  - Imports `create_ilp_transcript_counter`; branches on `args.no_ilp` to
    pick ILP vs default transcript model counter.
- **`isoquant_lib/dataset_processor.py`**: trivial (import cleanup).
- **CI/workflows**: accepted master's reshuffle — all the deleted
  `.github/workflows/Barcode.*.yml`, `SC.*.yml`, etc. moved under
  `isoquant_tests/github/`. `ILP.PB.yml` updated in place.

### 5.2 Integration bugs fixed during verification (this session)

All under `cell_type_decomposition_26`, un-committed at the time of
writing. Purpose: make the merged tree behave identically to master when
ILP is not selected, and keep the ILP path callable.

1. **`isoquant_lib/intron_graph.py`** — `assignment.read_group` is now
   `list[str]` (multi-group master API), but `IntronCollector.collect_introns`
   and `IntronGraph.construct` used it as a dict key, causing
   `unhashable type: list`. Added `_cell_type_key(read_group)` that
   returns `read_group` if `str`, `read_group[0]` if non-empty list,
   else `"NA"` (so `CellTypeTree` can filter it via
   `ignore_missing`/`remove_missing`). Call sites updated.
2. **`isoquant_lib/intron_graph.py`** — `simplify()` was calling
   `self.error_correction()` unconditionally, so the default pipeline
   silently dropped edges with weight < 3 and diverged from master
   (test yielded known=31/nic=3 instead of 28/1). Gated the call on
   `params.model_construction_strategy == "ilp_model" and not params.no_ilp`
   so it only runs in the ILP path.
3. **`isoquant_lib/assignment_aggregator.py`** — `if/else` created the
   ILP or default counter, then an unconditional second assignment
   overwrote it with the default counter. Collapsed into a single
   if/else that also gates on strategy:
   ```python
   if args.model_construction_strategy == "ilp_model" and not args.no_ilp:
       self.transcript_model_counter = create_ilp_transcript_counter(...)
   else:
       self.transcript_model_counter = create_transcript_counter(...)
   ```
   Fixes the bug where ILP counters never actually fired.
4. **`isoquant_lib/graph_based_model_construction.py`**:
   - `GraphBasedModelConstructor` stores args as `self.args` (master);
     the ported ILP code referenced `self.params.*` in three places.
     Fixed `self.params.min_known_count` →
     `self.args.min_known_count` (line ~728) and
     `self.params.report_unstranded` → `self.args.report_unstranded`
     (line ~742).
   - Added `self.use_ilp = (args.model_construction_strategy ==
     "ilp_model" and not args.no_ilp)` in `__init__`; ILP branch in
     `process()` now gates on `self.use_ilp` instead of
     `not self.args.no_ilp`. Rationale: `--no_ilp` default is `False`,
     so the old gate unconditionally selected ILP when strategy was
     anything — including `default_ont`.
   - Moved `self.ilp_counts()` inside the ILP branch.
   - `self.ilp_solution_assignement` initialised to `[]` instead of
     `None` so downstream iteration is safe when ILP short-circuits.
   - Technical-replicas check updated to the master pattern (look up
     `read_group[self.file_name_group_idx]` instead of the
     single-string `read_group`).
   - Restored the second `assign_reads_to_models(read_assignment_storage)`
     call after `filter_transcripts()` in the default branch — master
     intentionally reassigns post-filter.

### 5.3 Verification matrix

All runs used `export PATH="$HOME/bin:$PATH"` for `minimap2`/`samtools`
and wrote outputs under `./local_tmp/` in the repo (never `/tmp`).

| Scenario | Command | Result |
|---|---|---|
| Smoke (default strategy) | `./isoquant.py --test -o ./local_tmp/isoquant_test_out` | PASS — counts match master (known=28, novel_in_catalog=1). |
| Single grouped strategy | `--bam chr9.4M.ont.sim.polya.bam --read_group file:chr9.4M.ont.sim.read_groups.tsv:0:1` | Linear + matrix outputs produced for both known and discovered transcript models. |
| Multi-strategy grouping | `--read_group file:... file_name` | Both strategies produced independent `*_grouped_file0_col1_*` and `*_grouped_file_name_*` outputs. |
| Unit tests | `pytest isoquant_tests/test_read_groups.py isoquant_tests/test_string_pools.py` | 55 passed. |

### 5.4 Known gaps / still TODO on this branch

- ILP path itself (`--model_construction_strategy ilp_model`) has **not**
  been end-to-end tested in this session — fixes above make it *callable*
  and keep the default path correct, but a real ILP run still depends on
  `flowpaths` / `gurobipy` and a ground-truth DB.
- `ILPFlowCounter.add_ilp_data(rg_weights)` currently indexes by a single
  group key; the multi-group `CompositeCounter` will call it with the
  whole `list[str]` → need to pick / configure which `group_index`
  represents the cell-type axis (see §4.3).
- No cell-type grouping strategy is wired in `read_groups.py` yet. A
  user wanting cell-type-aware ILP still has to dump cell type into the
  first `--read_group` slot (which is what `_cell_type_key` grabs).
- `--no_ilp` still defaults to `False`; combined with the strategy gate
  in §5.2 the behaviour is fine (default strategies = no ILP), but the
  flag is effectively redundant and could be removed.
- `isoquant_tests/test_intron_graph.py` is an empty collection — add
  coverage for `_cell_type_key`, the gated `error_correction`, and
  `split_to_connected_subgraphs`.
