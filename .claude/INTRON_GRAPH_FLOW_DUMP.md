# Intron Graph Flow-Network Dump

Dev notes on `isoquant_lib/intron_graph_flow.py` and the
`--dump_intron_graphs` flag. The feature converts IsoQuant's
`IntronGraph` into a plain integer-labelled flow network and dumps it to
disk so external flow-decomposition / ILP tooling can consume it without
needing to re-implement IsoQuant's coordinate-tuple vertex representation.

## 1. Why

IsoQuant's `IntronGraph` identifies vertices by coordinate tuples:

- `(start, end)` for real introns (both positive ints)
- `(VERTEX_polya, pos)` / `(VERTEX_read_end, pos)` for terminal vertices
- `(VERTEX_polyt, pos)` / `(VERTEX_read_start, pos)` for starting vertices

where `VERTEX_*` are sentinel negatives (`-10`, `-11`, `-20`, `-21`).

External flow-decomposition libraries (e.g. `flowpaths`,
`networkx.max_flow`, Gurobi encoders) typically expect:

- consecutive integer node ids,
- a single super-source and a single super-target,
- an edge list + per-edge flow dict.

This module provides both the conversion (`Intron2Graph`) and a
dump-to-disk helper (`dump_flow_graph`).

## 2. Provenance

Ported from `origin/ilp_models:src/encode_ilp_gurobi.py::Intron2Graph`
(commits around `702caca`, `fd4207f`, `d81f2d7`). That version lived
next to an earlier Gurobi ILP encoder that was later superseded on
`origin/ilp_model` by `flowpaths`-based models with `networkx.DiGraph`
vertices labelled by `str(intron_tuple)` (i.e. strings, not ints). The
consecutive-int format only lived on `ilp_models`.

See `.claude/CELL_TYPE_DECOMPOSITION.md` §3 for context on the two
ILP lines.

## 3. Conversion (`Intron2Graph`)

```python
from isoquant_lib.intron_graph_flow import Intron2Graph
flow = Intron2Graph(intron_graph)                              # default
flow = Intron2Graph(intron_graph, add_super_source_target=True) # opt-in
```

Attributes:

| name | type | meaning |
|---|---|---|
| `source` | `Optional[int]` | super-source id, or `None` when not added |
| `target` | `Optional[int]` | super-target id, or `None` when not added |
| `add_super_source_target` | `bool` | flag the network was built with |
| `intron2vertex` | `dict[tuple, int]` | tuple-vertex → int id |
| `vertex2intron` | `dict[int, tuple]` | int id → tuple-vertex |
| `edge_list` | `list[(int, int)]` | directed edges |
| `flow_dict` | `dict[(int, int), int]` | edge weights (flow) |

Numbering order (vertex ids are consecutive starting at `0`):

1. (optional) super-source `0` when `add_super_source_target=True`
2. real intron vertices (in `intron_collector.clustered_introns` insertion order)
3. starting vertices referenced from `incoming_edges`
4. terminal vertices referenced from `outgoing_edges`
5. (optional) super-target as the final id when `add_super_source_target=True`

Super-source / super-target are off by default — there is no CLI flag,
the toggle is purely internal for downstream tooling that needs a
single-source / single-sink wrap. When enabled, the extra edges are:

- `source → each starting vertex`, weight =
  Σ outgoing `edge_weights` from that starting vertex
- `each terminal vertex → target`, weight =
  Σ incoming `edge_weights` to that terminal vertex

Round-trip helpers:

```python
flow.transcript_to_path(intron_tuple_list)  # -> list[(u, v)]
flow.path_to_transcript(int_path)           # -> list[intron tuple]
```

### 3.1 Terminal edge weights

Master's `IntronGraph.edge_weights` was historically populated only in
`add_edge(v1, v2)` (intron → intron edges from `construct()`); the
terminal / starting edges added later in `attach_terminal_positions()`
were left at their `defaultdict(int)` zero. This branch backports the
`origin/ilp_models` patch into `attach_transcpt_ends`, so
`edge_weights` now carries the per-position cluster counts for:

- `starting_vertex → intron` (counts from `cluster_polya_positions` /
  `cluster_terminal_positions` on the polyT / read_start side)
- `intron → terminal_vertex` (same, on the polyA / read_end side)

The implicit `source → starting` / `terminal → target` super edges
(only added when `add_super_source_target=True`) are still computed as
the sum of the underlying `starting_vertex → intron` /
`intron → terminal_vertex` weights, so they now also carry real flow.

Only `(intron_a, intron_b)` edges where neither `intron_b ∈ outgoing[intron_a]`
nor the polyA/terminal lookup applies remain at the default `0`. If
you're consuming `flow_dict` directly, expect terminal entries to be
populated.

## 4. Dump helper (`dump_flow_graph`)

```python
dump_flow_graph(intron_graph, chr_id, gene_id, out_dir)
```

- Skips empty graphs (no clustered introns and no edges).
- Each gene gets its own subdirectory `<out_dir>/<chr>.<gene>/`, created
  with `os.makedirs(..., exist_ok=True)`.
- `gene_id` is sanitised: `/` → `_` (rare but possible in gffutils ids).
- Missing TSV fields are written as `*` (the same gap marker used in
  `paths.tsv`) so columns stay aligned and grep-friendly.

Writes two TSVs per gene:

**`vertices.tsv`** — `vertex_id  type  chr  start  end  weight`

| vertex kind | type column | start, end | weight |
|---|---|---|---|
| super-source (opt-in) | `source` | `*`, `*` | `*` |
| intron | `intron` | `(start, end)` | `clustered_introns[intron]` |
| polyA sink | `polya` | `(pos, pos)` | Σ incoming edge weights |
| read-end sink | `read_end` | `(pos, pos)` | Σ incoming edge weights |
| polyT source | `polyt` | `(pos, pos)` | Σ outgoing edge weights |
| read-start source | `read_start` | `(pos, pos)` | Σ outgoing edge weights |
| super-target (opt-in) | `target` | `*`, `*` | `*` |

`weight` is the cluster count for the vertex. Intron rows pull from
`intron_graph.intron_collector.clustered_introns`; terminal / starting
rows aggregate the per-edge counts that `attach_transcpt_ends` writes
into `intron_graph.edge_weights` (see §3.1) — equivalently, the sum of
incoming (terminal) or outgoing (starting) weights in `flow.flow_dict`.

**`edges.tsv`** — `u  v` by default. With the internal
`include_edge_weights=True` flag, a third `weight` column is appended;
values come from `flow.flow_dict`. All real edges (intron → intron,
intron → terminal, starting → intron) carry the populated edge
weight; only the implicit `source → starting` / `terminal → target`
super edges (opt-in via `add_super_source_target=True`) remain as
documented sums.

**`paths.tsv`** (only when `--ground_truth_counts` is also
set and at least one transcript in the gene appears in that TSV) —
`transcript_id  count  count_scaled  status  path  path_simple  missing_vertices  missing_edges`

`count` is the original (pre-downsampling) ground-truth abundance from
`--ground_truth_counts`. `count_scaled = round(count /
gene_info.coverage_scale_factor)` is the expected post-downsampling
load — it matches the observed-read weights elsewhere in the dump
(`vertices.tsv` weights, `read_subpaths.tsv` `read_count`) when
`--max_coverage_small_chr` / `--max_coverage_normal_chr` triggered the
"process 1 read out of every N" path. With no downsampling (scale ==
1, the common case) the two columns are equal. Caveat: this is the
expected average — the keep-rule is `counter % N == 0`, so per-region
counts can drift by ±1 around the rounded value.

| status | meaning |
|---|---|
| `monoexonic` | transcript has no introns — `path` is empty, terminals are skipped |
| `ok` | every threaded slot maps to a graph vertex AND every consecutive pair is an edge in `intron_graph.outgoing_edges` |
| `disconnected` | all vertices present, but at least one consecutive pair is not an edge in the graph |
| `partial` | at least one slot has no vertex (intron in `discarded_introns` / not in `clustered_introns`, or terminal boundary with no graph vertex inside `apa_delta`) |

A threaded path covers `[low_terminal, intron_1, …, intron_N,
high_terminal]`. Each annotated intron runs through
`intron_collector.substitute(...)` and is translated to a vertex id via
`flow.intron2vertex`. Each transcript's exon-boundary positions
(`gene_info.all_isoforms_exons[t_id][0][0]` for the low side,
`[-1][1]` for the high side) are matched to the closest graph
starting / terminal vertex within `intron_graph.params.apa_delta`. No
synthetic vertices are introduced — if no graph vertex sits in the
radius the slot becomes `*` and bumps `missing_vertices`.

`path` rendering: each slot becomes either its integer vertex id or
the gap token `*`, and adjacent tokens are joined by an edge separator
that carries connectivity:

- `-` (dash): both sides are real vertices and the edge exists in
  `intron_graph.outgoing_edges`
- `|` (pipe): both sides are real vertices but no edge is present
- `-` (default): used whenever either side is `*` (we can't classify
  edge presence when a vertex is missing)

So `5-0-1-2|3-4-7` means "low terminal 5 → intron 0 → … → intron 4 →
high terminal 7, with edge 2→3 missing". `*-0-1-2-3-*` means "neither
boundary mapped, the intron chain is fully connected".

`path_simple` rendering: a comma-separated list of the vertex ids that
actually map to the graph (terminals included when matched) — missing
slots are dropped, no `*` placeholders, edge state ignored. Useful when
you only care about the projected vertex chain (e.g. for set-style
comparisons across transcripts).

Counts:
- `missing_vertices` — number of slots (introns or terminals) with no
  graph vertex
- `missing_edges` — number of consecutive slot pairs whose edge is not
  in `intron_graph.outgoing_edges` (either because an endpoint is
  missing or because the edge itself is absent). Includes the implicit
  `low_terminal → first_intron` and `last_intron → high_terminal`
  edges.

Note: only edges through real graph vertices are validated; the
super-source / super-target wiring (when `add_super_source_target=True`)
is not.

**`read_subpaths.tsv`** (always written when at least one read threaded
into a path) — `read_count  is_fl  status  path  path_simple
missing_vertices  missing_edges`

One row per distinct read-supported path (the keys of
`IntronPathStorage.paths`), sorted by descending `read_count` (FL paths
break ties first):

- `read_count` — number of reads producing this exact path tuple
- `is_fl` — `1` when the path is in `IntronPathStorage.fl_paths` (both
  starting and terminal vertices threaded by `thread_starts` /
  `thread_ends`); `0` otherwise (partial paths missing one or both ends)
- `status` / `path` / `path_simple` / `missing_vertices` /
  `missing_edges` — same encoding as `paths.tsv`. Threaded read paths
  are built from real graph vertices, so `partial` should not occur and
  `disconnected` is rare; if either appears it points at a bug in the
  threading logic.

This file is independent of `--ground_truth_counts` and `--dump_ref_data`
— it always materialises whenever `--dump_intron_graphs` is on, since
the path-storage step now runs before the dump.

**`ref_vertices.tsv` / `ref_edges.tsv`** (only when `--dump_ref_data` is
passed and the gene has annotated transcripts) enumerate every
annotated slot — introns AND the per-transcript 5'/3' exon boundaries
— plus every consecutive pair across all transcripts, labelling each
with its graph status:

- `ref_vertices.tsv` — `kind  ref_start  ref_end  status  vertex_id  graph_start  graph_end`
  - `kind` ∈ `intron` / `starting` / `terminal` (`starting` = low-coord
    boundary matched against polyT / read_start; `terminal` =
    high-coord boundary matched against polyA / read_end)
  - `status` for introns ∈ `in_graph` / `discarded` / `unmapped`; for
    boundaries ∈ `in_graph` / `unmapped`
  - terminal `unmapped` means no graph vertex sits within
    `intron_graph.params.apa_delta` of the annotated boundary —
    nothing synthesized, just recorded as missing
  - `*` for `vertex_id` when unresolved; `*` for `graph_start`/
    `graph_end` on `discarded` introns and `unmapped` boundaries
- `ref_edges.tsv` — `kind  u_start  u_end  v_start  v_end  status  u_id  v_id`
  - `kind` ∈ `intron` / `starting` / `terminal` (matching the vertex
    kind on the boundary side; `starting` = `low_boundary →
    first_intron`, `terminal` = `last_intron → high_boundary`)
  - `status` ∈ `in_graph` / `missing_edge` / `missing_vertex`; absent
    ids rendered as `*`

Independent of `--ground_truth_counts`: a gene is written iff it has
annotated introns, regardless of whether it appears in the counts TSV.

The `chr` column uses `self.gene_info.chr_id`; for intergenic / novel
gene regions without a `gene_db_list`, the per-gene directory falls back
to `<chr>.region_<start>_<end>/`.

## 5. Pipeline integration

The only call site is in
`isoquant_lib/graph_based_model_construction.py::GraphBasedModelConstructor.process`,
right after `IntronPathStorage.fill()` so the dump can include the
read-threaded paths:

```python
self.intron_graph = IntronGraph(self.args, self.gene_info, read_assignment_storage)
self.path_processor = IntronPathProcessor(self.args, self.intron_graph)
self.path_storage = IntronPathStorage(self.args, self.path_processor)
self.path_storage.fill(read_assignment_storage)
dump_dir = getattr(self.args, "dump_intron_graphs_dir", None)
if dump_dir:
    gene_ids = [g.id for g in self.gene_info.gene_db_list] if self.gene_info.gene_db_list else []
    gene_tag = "_".join(gene_ids) if gene_ids else "region_%d_%d" % (self.gene_info.start, self.gene_info.end)
    dump_flow_graph(self.intron_graph, self.gene_info.chr_id, gene_tag,
                    dump_dir, ...,
                    path_storage=self.path_storage)
```

Dump happens once per gene region, after `IntronGraph.__init__` has run
`simplify()` + `attach_terminal_positions()` AND `IntronPathStorage` has
threaded every read into its graph path. The graph is exactly what the
downstream `construct_fl_isoforms()` /
`construct_assignment_based_isoforms()` / ILP encoder would see.

`gene_info` is always forwarded to `dump_flow_graph`; the extra
`ground_truth_counts` kwarg comes from `args.ground_truth_counts_map`
(set only when the `--ground_truth_counts` flag is given) — absent it,
`dump_flow_graph` skips `paths.tsv` generation entirely.

## 6. CLI flag

`isoquant.py`, under the output-setup group (hidden behind `--full-help`):

```
--dump_intron_graphs
    dump per-gene intron graphs as integer-labelled flow networks
    (vertices.tsv + edges.tsv) into <output>/intron_graphs/

--ground_truth_counts FILE
    TSV with ground-truth transcript counts (col1: transcript_id
    matching --genedb, col2: count). Only used together with
    --dump_intron_graphs; each gene's dump also gets a paths.tsv
    mapping transcripts onto integer vertex paths with weights.

--dump_ref_data
    Only used together with --dump_intron_graphs; each gene's dump
    also gets ref_vertices.tsv / ref_edges.tsv enumerating every
    annotated intron and every consecutive intron pair with its
    graph status.
```

Default off — all three flags have zero effect on pipeline output
when unset. When `--dump_intron_graphs` is given, `set_additional_params`
creates `<output>/intron_graphs/` up front and stashes the absolute
path on `args.dump_intron_graphs_dir` (so downstream code doesn't have
to recompute it).

Loader behaviour for `--ground_truth_counts` (in
`set_additional_params`):

- A warning is logged and the flag is ignored if
  `--dump_intron_graphs` is not also set.
- The file must exist, otherwise `IsoQuantExitCode.INPUT_FILE_NOT_FOUND`.
- Blank lines and lines starting with `#` are skipped.
- Only the first two tab-separated columns are read; extra columns
  are ignored. Malformed count values log a warning and are skipped.
- The parsed dict is stashed on `args.ground_truth_counts_map`.

Caveats:

- `--test` mode **ignores** the flag — `TestMode` runs IsoQuant with a
  hard-coded argv. Run a real command if you need the dump:
  ```bash
  ./isoquant.py --reference <fa> --genedb <gtf> --bam <bam> \
                --data_type nanopore --dump_intron_graphs -o <out>
  ```
- When running with `--resume`, the `reads_processed_lock_file` short-
  circuits model construction and no new dumps are produced — delete
  the lock (or rerun with `--clean_start`) to regenerate.

## 7. Smoke verification

Expected shape for the chr9.4M simple-data graph around
`ENSMUSG00000042360` (~11 reads spanning the main chain):

```
vertex_id  type        chr   start     end       weight
0          intron      chr9  3192413   3192819   11
1          intron      chr9  3190806   3192311   11
...
6          polyt       chr9  3190120   3190120   10   # sum of polyt→intron edges
9          read_end    chr9  3199797   3199797   8    # sum of intron→read_end edges

u  v
1  0             # intron-intron
0  2             # intron-intron
...
6  0             # polyt -> intron
2  9             # intron -> read_end
```

Default output omits super-source / super-target and the per-edge
weight column. Pass `add_super_source_target=True` to wrap the graph
(prepending `0  source` and appending `N  target` rows) and
`include_edge_weights=True` to recover the third column on
`edges.tsv`.

## 8. Future work / gaps

- Populate terminal-edge weights from `clustered_polyas` /
  `terminal_positions` so downstream flow decomposition has a sensible
  source/sink load.
- Optionally emit a companion JSON per gene with the intron-correction
  map, path constraints, and `known_isoforms_in_graph` so external
  ILP tooling can replicate IsoQuant's preprocessing.
- If we ever want cross-gene dumps (e.g. one NetworkX pickle per
  chromosome), collapse the per-gene directories into a single indexed
  file — the current layout is one subdirectory of TSVs per gene region.
