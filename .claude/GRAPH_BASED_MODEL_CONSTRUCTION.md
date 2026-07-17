# Graph-Based Transcript Model Construction

This document describes the novel transcript discovery algorithm implemented in `src/graph_based_model_construction.py`. The algorithm is based on intron graphs and is described in the [Nature Biotechnology paper](https://www.nature.com/articles/s41587-022-01565-y).

## Overview

The transcript model construction algorithm discovers novel transcripts from long RNA reads by:
1. Building an **intron graph** from read alignments
2. Finding **paths** through the graph that represent full-length transcripts
3. **Filtering** low-confidence predictions
4. **Assigning** reads to discovered transcripts for quantification

## Data Structures

### Intron Graph (`src/intron_graph.py`)

The intron graph is the core data structure where:
- **Vertices** = Splice junctions (introns), represented as `(start, end)` tuples
- **Edges** = Consecutive introns observed in at least one read
- **Terminal vertices** = Special vertices marking transcript starts/ends:
  - `VERTEX_polya (-10)` - PolyA-confirmed 3' end (forward strand)
  - `VERTEX_read_end (-11)` - Read end without polyA confirmation
  - `VERTEX_polyt (-20)` - PolyT-confirmed 5' end (reverse strand, i.e., polyA on reverse complement)
  - `VERTEX_read_start (-21)` - Read start without polyT confirmation

```
Example intron graph:

[read_start] → intron1 → intron2 → intron3 → [polya]
                  ↘                    ↗
                   intron2' → intron3'
```

### Terminal Vertex Constants

Defined in `src/intron_graph.py`:
```python
VERTEX_polya = -10       # PolyA-confirmed 3' end (forward strand)
VERTEX_read_end = -11    # Read end without polyA confirmation
VERTEX_polyt = -20       # PolyT-confirmed 5' end (reverse strand)
VERTEX_read_start = -21  # Read start without polyT confirmation
```

Helper functions:
- `is_polya(v)` - Returns `v[0] == VERTEX_polya`
- `is_polyt(v)` - Returns `v[0] == VERTEX_polyt`
- `is_terminal_vertex(v)` - Returns `v[0] < 0`

### Key Classes

#### `IntronCollector`
Collects and clusters introns from read alignments:
- `clustered_introns` - Dict of intron → count after clustering
- `intron_correction_map` - Maps similar introns to their canonical form (key: original intron, value: canonical intron)
- `discarded_introns` - Low-count introns removed from consideration
- `known_introns` - Introns from reference annotation (always trusted, never discarded)

#### `IntronGraph`
Main graph structure:
- `outgoing_edges[intron]` - Set of vertices reachable from intron
- `incoming_edges[intron]` - Set of vertices that lead to intron
- `edge_weights[(v1, v2)]` - Number of reads supporting edge
- `max_coverage` - Maximum intron coverage in the graph

#### `IntronPathStorage`
Stores paths (sequences of introns) extracted from reads:
- `paths` - Dict of path → count
- `fl_paths` - Set of full-length paths (have both start and end vertices)
- `paths_to_reads` - Maps paths to the reads that support them

#### `IntronPathProcessor`
Threading logic for mapping reads to graph paths:
- `thread_introns()` - Maps read introns to graph vertices
- `thread_ends()` - Finds terminal vertex for read end
- `thread_starts()` - Finds starting vertex for read start

#### `GraphBasedModelConstructor`
Main orchestrator class that runs the entire pipeline.

#### `TranscriptToGeneJoiner`
Post-processing class that merges novel transcripts into gene loci.

## Algorithm Pipeline

### Phase 1: Graph Construction (`IntronGraph.__init__`)

```
Input: Read assignments with corrected introns
Output: Intron graph with terminal vertices

1. Collect all introns from reads (IntronCollector.collect_introns)
2. Cluster similar introns (IntronCollector.cluster_introns):
   - Known introns are always kept
   - Novel introns within `delta` bp of a higher-count intron are merged
   - Low-count introns are discarded
3. Build edges from consecutive introns in reads (IntronGraph.construct)
4. Simplify graph (IntronGraph.simplify):
   a. Remove tips and bulges (clean_tips_and_bulges)
   b. Remove singleton dead ends (remove_singleton_dead_ends)
   c. Remove isolated vertices (remove_isolates)
5. Attach terminal positions (attach_terminal_positions):
   - Cluster polyA/polyT positions
   - Add read start/end positions
```

#### Intron Clustering Algorithm (`IntronCollector.cluster_introns`)

**Location**: `intron_graph.py:76-130`

The clustering algorithm merges similar novel introns:

```python
def cluster_introns(self):
    # Sort by count descending - higher count introns are "canonical"
    sorted_introns = sorted(self.introns.items(),
                           key=lambda x: x[1], reverse=True)

    for intron, count in sorted_introns:
        if intron in self.known_introns:
            # Known introns always kept as-is
            self.clustered_introns[intron] = count
            continue

        # Check if similar to existing clustered intron
        merged = False
        for existing in self.clustered_introns:
            if (abs(intron[0] - existing[0]) <= self.delta and
                abs(intron[1] - existing[1]) <= self.delta):
                # Merge into existing (higher-count) intron
                self.intron_correction_map[intron] = existing
                self.clustered_introns[existing] += count
                merged = True
                break

        if not merged:
            if count >= self.min_novel_intron_count:
                self.clustered_introns[intron] = count
            else:
                self.discarded_introns.add(intron)
```

**Key behaviors**:
- Processing order matters: higher-count introns processed first
- Known introns are never merged or discarded
- `intron_correction_map` allows translating original introns to canonical form
- Introns with count < `min_novel_intron_count` are discarded

### Phase 2: Path Extraction (`IntronPathStorage.fill`)

```
For each read assignment:
  1. Thread read introns through graph (skip if any intron was discarded)
  2. Find terminal vertex for read end:
     - If read has polyA → find matching VERTEX_polya
     - Else if read extends past all outgoing introns → use VERTEX_read_end
  3. Find starting vertex for read start:
     - If read has polyT → find matching VERTEX_polyt
     - Else if read extends past all incoming introns → use VERTEX_read_start
  4. Store path and increment count
  5. If path has both start and end → mark as full-length (fl_paths)
```

### Phase 3: Transcript Construction (`GraphBasedModelConstructor.process`)

The main `process()` method orchestrates:

```python
def process(self, read_assignment_storage):
    # Build graph and extract paths
    self.intron_graph = IntronGraph(...)
    self.path_processor = IntronPathProcessor(...)
    self.path_storage = IntronPathStorage(...)
    self.path_storage.fill(read_assignment_storage)

    # Get known isoforms as graph paths
    self.known_isoforms_in_graph = self.get_known_spliced_isoforms(...)

    # Construct transcripts
    self.construct_fl_isoforms()              # Full-length novel + known
    self.construct_assignment_based_isoforms() # Non-FL known isoforms

    # Filter and validate
    self.pre_filter_transcripts()
    self.assign_reads_to_models()             # First assignment pass
    self.filter_transcripts()                 # Main filtering
    self.assign_reads_to_models()             # Second assignment pass

    # Join transcripts to genes
    transcript_joiner = TranscriptToGeneJoiner(...)
    self.transcript_model_storage = transcript_joiner.join_transcripts()

    # Forward counts
    self.forward_counts(read_assignment_storage)
```

### Phase 3a: Full-Length Isoform Construction (`construct_fl_isoforms`)

For each full-length path (sorted by length descending):

```
1. Extract intron path (remove terminal vertices)
2. Create exon blocks from introns
3. Check if path matches a known isoform:
   - Construct profile and assign to reference
   - If matches known isoform → create known transcript model
4. If novel:
   - Apply coverage cutoff (min_novel_count)
   - Check mono-intronic requirements (polyA, canonical strand)
   - Check strand detection requirements
   - Check technical replicas (if enabled)
   - Classify as NIC or NNIC based on whether all introns are known
   - Create novel transcript model
5. Save reads supporting the path
```

### Phase 3b: Assignment-Based Isoform Construction (`construct_assignment_based_isoforms`)

Handles reads that were uniquely assigned to known isoforms but didn't form FL paths:

```
For each read with unique assignment:
  1. Mono-exon matches → track for mono-exon isoform construction
  2. Spliced matches → track for non-FL spliced isoform construction

Then:
  - construct_monoexon_isoforms(): Add mono-exon reference transcripts
  - construct_nonfl_isoforms(): Add spliced reference transcripts (if not FL-only mode)
  - construct_monoexon_novel(): Discover novel mono-exon transcripts (if enabled)
```

### Phase 4: Filtering (`pre_filter_transcripts`, `filter_transcripts`)

#### Pre-filtering
For simple transcripts (≤2 exons):
- Apply minimum count cutoff
- Apply mapping quality cutoff

#### Main Filtering
For each novel transcript:
1. Calculate component coverage (max intron coverage in connected component)
2. Calculate coverage cutoff:
   - Spliced: `max(min_novel_count, min_novel_count_rel * component_coverage)`
   - Mono-intronic: `max(min_novel_count, min_mono_count_rel * component_coverage)`
3. Remove if count < cutoff
4. Remove if similar to existing isoform (detect_similar_isoforms)
5. Correct novel transcript ends based on read evidence

#### Similar Isoform Detection (`detect_similar_isoforms`)
```
For each multi-exon model:
  Create temporary GeneInfo from model
  For each other model:
    If other model can be assigned to this model → mark for substitution
```

### Phase 5: Read Assignment (`assign_reads_to_models`)

After filtering, re-assign reads to surviving transcripts:

```
Create GeneInfo from all transcript models
For each read not already assigned:
  Construct profile
  Assign to transcript models
  If consistent assignment:
    Increment counters
    Store read → transcript mapping
```

### Phase 6: Gene Joining (`TranscriptToGeneJoiner`)

Merges novel transcripts into gene loci based on:
- Strand compatibility
- Intron overlap (Jaccard similarity)
- Positional overlap

```
For each pair of genes:
  score = positional_overlap + intronic_overlap

While best_score >= 0.1:
  Merge gene pair with highest score
  Prefer merging into reference genes
```

## Key Parameters

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| `delta` | Tolerance for splice site matching | 6 bp |
| `apa_delta` | Tolerance for polyA site matching | 50 bp |
| `min_novel_count` | Minimum reads for novel transcript | 5 (ONT), 3 (PacBio) |
| `min_novel_count_rel` | Relative coverage cutoff | 0.02 (2%) |
| `min_mono_count_rel` | Relative cutoff for mono-intronic | 0.005 |
| `min_known_count` | Minimum reads for known transcript | 1 |
| `min_novel_intron_count` | Minimum reads for novel intron in graph | 2 |
| `graph_clustering_distance` | Distance for intron clustering | 10 bp |
| `graph_clustering_ratio` | Coverage ratio for clustering | 0.5 |
| `simple_models_mapq_cutoff` | MAPQ cutoff for simple transcripts | 10 |

## Transcript Types

From `TranscriptModelType` enum:
- `known` - Matches reference transcript exactly
- `novel_in_catalog` (NIC) - All introns are known, but novel combination
- `novel_not_in_catalog` (NNIC) - Contains at least one novel intron

## Terminal Position Clustering

### Overview

Terminal position clustering happens in `attach_terminal_positions()` (lines 411-470 in `intron_graph.py`). For each intron, it processes both 3' ends (polyA/read_end) and 5' ends (polyT/read_start).

### PolyA Position Clustering (`cluster_polya_positions`)

**Location**: `intron_graph.py:518-550`

**Algorithm** (greedy clustering):
```python
def cluster_polya_positions(position_dict, apa_delta, known_positions,
                            read_end=True, min_abs=1, min_rel=0.05):
    clustered_counts = {}

    while position_dict:
        # 1. Select highest-count position
        best_pair = max(position_dict.items(), key=lambda x: x[1])
        top_position = best_pair[0]

        # 2. Snap to reference if close
        for known_pos in known_positions:
            if abs(top_position - known_pos) <= apa_delta:
                top_position = known_pos
                break

        # 3. Cluster all positions within apa_delta window
        total_count = 0
        for pos in range(top_position - apa_delta, top_position + apa_delta + 1):
            if pos in position_dict:
                total_count += position_dict[pos]
                del position_dict[pos]

        clustered_counts[top_position] = total_count

    # 4. Filter by coverage cutoff
    max_count = max(clustered_counts.values()) if clustered_counts else 0
    cutoff = max(min_abs, max_count * min_rel)
    return {pos: count for pos, count in clustered_counts.items()
            if count >= cutoff}
```

**Key behaviors**:
- Known reference positions are preferred (snapping within `apa_delta`)
- Greedy selection means order-dependent results
- Positions within `apa_delta` of selected peak are merged
- Final filtering removes low-count clusters

### Simple Terminal Clustering (`cluster_terminal_positions`)

**Location**: `intron_graph.py:552-561`

For non-polyA terminals (read ends without polyA confirmation):
```python
def cluster_terminal_positions(position_dict, read_end=True, min_abs=1):
    if not position_dict:
        return {}

    total_count = sum(position_dict.values())
    if total_count < min_abs:
        return {}

    # Simply take max (for 3' end) or min (for 5' end) position
    if read_end:
        pos = max(position_dict.keys())
    else:
        pos = min(position_dict.keys())

    return {pos: total_count}
```

**Note**: This is much simpler - just takes the extreme position.

## Terminal Vertex Logic

### Finding Read Ends (`thread_ends`)

**Location**: `intron_graph.py:910-940`

**Detailed algorithm**:
```python
def thread_ends(self, intron, end, trusted=False):
    """
    Find terminal vertex for read 3' end.

    Args:
        intron: Last intron in read
        end: Read end position
        trusted: True if read has polyA confirmation

    Returns:
        Terminal vertex tuple or None
    """
    outgoing = self.intron_graph.outgoing_edges.get(intron, set())
    terminal_vertices = [v for v in outgoing if is_terminal_vertex(v)]

    if not terminal_vertices:
        return None

    # Separate polyA and read_end vertices
    polya_vertices = [v for v in terminal_vertices if is_polya(v)]
    read_end_vertices = [v for v in terminal_vertices if not is_polya(v)]

    if trusted:
        # Read has polyA - find matching VERTEX_polya
        for v in polya_vertices:
            if abs(v[1] - end) <= self.apa_delta:
                return v
        return None

    # Check if read is internal (ends before downstream introns)
    non_terminal = [v for v in outgoing if not is_terminal_vertex(v)]
    if non_terminal:
        rightmost_exon_end = max(v[0] for v in non_terminal) - 1
        if end <= rightmost_exon_end + self.delta:
            return None  # Internal read

    # Find appropriate terminal vertex
    all_terminals = polya_vertices + read_end_vertices
    rightmost = max(all_terminals, key=lambda v: v[1])

    # Complex conditions for selection...
    if abs(end - rightmost[1]) <= self.apa_delta:
        return rightmost

    return None
```

### Finding Read Starts (`thread_starts`)

**Location**: `intron_graph.py:942-969`

**Detailed algorithm**:
```python
def thread_starts(self, intron, start, trusted=False):
    """
    Find starting vertex for read 5' end.

    Args:
        intron: First intron in read
        start: Read start position
        trusted: True if read has polyT confirmation

    Returns:
        Terminal vertex tuple or None
    """
    incoming = self.intron_graph.incoming_edges.get(intron, set())
    terminal_vertices = [v for v in incoming if is_terminal_vertex(v)]

    if not terminal_vertices:
        return None

    polyt_vertices = [v for v in terminal_vertices if is_polyt(v)]
    read_start_vertices = [v for v in terminal_vertices if not is_polyt(v)]

    if trusted:
        # Read has polyT - find matching VERTEX_polyt
        for v in polyt_vertices:
            if abs(v[1] - start) <= self.apa_delta:
                return v
        return None

    # Check if read is internal
    non_terminal = [v for v in incoming if not is_terminal_vertex(v)]
    if non_terminal:
        leftmost_exon_start = min(v[1] for v in non_terminal) + 1
        if start >= leftmost_exon_start - self.delta:
            return None  # Internal read

    # Find appropriate terminal vertex
    all_terminals = polyt_vertices + read_start_vertices
    leftmost = min(all_terminals, key=lambda v: v[1])

    if abs(start - leftmost[1]) <= self.apa_delta:
        return leftmost

    return None
```

## Edge Weights and Path Scoring

### Edge Weight Calculation

Edge weights represent the number of reads supporting consecutive introns:

```python
# In IntronGraph.construct()
for read_assignment in read_assignments:
    introns = read_assignment.corrected_introns
    for i in range(len(introns) - 1):
        edge = (introns[i], introns[i+1])
        self.edge_weights[edge] += 1
        self.outgoing_edges[introns[i]].add(introns[i+1])
        self.incoming_edges[introns[i+1]].add(introns[i])
```

### Path Count Tracking

Paths are counted in `IntronPathStorage`:

```python
class IntronPathStorage:
    def __init__(self):
        self.paths = defaultdict(int)       # path_tuple → count
        self.fl_paths = set()               # Full-length paths
        self.paths_to_reads = defaultdict(list)  # path → [read_ids]

    def add_path(self, intron_path, read_id, is_full_length):
        path_tuple = tuple(intron_path)
        self.paths[path_tuple] += 1
        self.paths_to_reads[path_tuple].append(read_id)
        if is_full_length:
            self.fl_paths.add(path_tuple)
```

### Coverage Cutoff Calculation

Novel transcripts must meet coverage thresholds:

```python
# In filter_transcripts()
def calculate_coverage_cutoff(self, intron_path):
    # Get max coverage of any intron in the path's connected component
    component_coverage = self.get_component_coverage(intron_path)

    if len(intron_path) == 1:
        # Mono-intronic: lower relative threshold
        return max(self.min_novel_count,
                   self.min_mono_count_rel * component_coverage)
    else:
        # Multi-intronic: standard threshold
        return max(self.min_novel_count,
                   self.min_novel_count_rel * component_coverage)
```

## Graph Simplification Details

### Tip/Bulge Removal (`clean_tips_and_bulges`)

A vertex is collapsed if:
1. It's within `graph_clustering_distance` of another vertex
2. Its count is less than `graph_clustering_ratio` × the other vertex's count
3. It's not a known intron

### Singleton Dead End Removal (`remove_singleton_dead_ends`)

From high-coverage introns, remove paths where:
- All outgoing edges lead to coverage-1 dead ends
- All incoming edges come from coverage-1 dead starts

### Isolated Vertex Removal (`remove_isolates`)

1. First, cluster similar isolated vertices
2. Then remove remaining isolated vertices with count < `min_novel_isolated_intron_abs`
3. Known introns are never removed

## Output

### Transcript Models

Each discovered transcript becomes a `TranscriptModel` with:
- `transcript_id` - Unique identifier (e.g., "transcript12345.chr1.nic")
- `gene_id` - Reference gene or novel gene ID
- `exon_blocks` - List of (start, end) exon coordinates
- `strand` - '+', '-', or '.'
- `transcript_type` - known, NIC, or NNIC
- `intron_path` - Path through intron graph

### Counters

Read counts are tracked in:
- `transcript_counter` - Per-transcript counts
- `gene_counter` - Per-gene counts
- `transcript_grouped_counters` - Grouped transcript counts (multi-group mode)
- `gene_grouped_counters` - Grouped gene counts (multi-group mode)

## Known Issues / Areas for Improvement

**For detailed issue tracking and planned improvements, see: `.claude/GRAPH_MODEL_ISSUES.md`**

### Critical Issues

1. **Terminal positions from reference**: When discovered transcript matches known intron chain, TSS/polyA taken from reference instead of reads
2. **Suboptimal polyA/TSS clustering**: Greedy algorithm in `attach_terminal_positions()` - new XGBoost-based peak detection algorithm planned
3. **Underutilized FL read information**: Full-length read terminals not fully trusted for transcript boundaries
4. **Inconsistent read assignment architecture**: Multiple assignment methods cause discrepancies; `add_read_info_raw()` bypasses normal counter architecture

### Code Organization

1. **Large monolithic class**: `GraphBasedModelConstructor` is ~900 lines with many responsibilities
2. **Class-level state**: `detected_known_isoforms` is a class variable, causing issues across genes
3. **Mixed concerns**: Graph construction, path finding, filtering, and counting are interleaved

### Algorithm Complexity

1. **Quadratic similarity detection**: `detect_similar_isoforms` compares all pairs of models
2. **Multiple assignment passes**: Reads are assigned twice due to architectural issues
3. **Repeated GeneInfo creation**: Creates temporary GeneInfo for each similarity check

### Parameter Sensitivity

1. **Many hardcoded thresholds**: Coverage ratios, distances, magic numbers
2. **No adaptive thresholds**: Same parameters for different coverage levels
3. **Technical replicas check**: Only works with `file_name` grouping

## Flow Diagram

```
Read Assignments
       │
       ▼
┌─────────────────┐
│ IntronCollector │ ──► Clustered introns
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  IntronGraph    │ ──► Vertices + Edges
└────────┬────────┘
         │ simplify()
         ▼
┌─────────────────┐
│ attach_terminal │ ──► Terminal vertices (polyA/polyT/read ends)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│IntronPathStorage│ ──► Paths + FL paths
└────────┬────────┘
         │
         ▼
┌─────────────────────────────────┐
│ construct_fl_isoforms()         │
│ construct_assignment_based_..() │ ──► Initial transcript models
└────────────────┬────────────────┘
                 │
                 ▼
┌─────────────────┐
│ pre_filter_...  │ ──► Remove low-quality simple transcripts
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│assign_reads_to..│ ──► First read assignment
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│filter_transcr.. │ ──► Coverage + similarity filtering
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│assign_reads_to..│ ──► Second read assignment
└────────┬────────┘
         │
         ▼
┌─────────────────────┐
│TranscriptToGeneJoiner│ ──► Merge transcripts into genes
└──────────┬──────────┘
           │
           ▼
┌─────────────────┐
│ forward_counts  │ ──► Final counting
└────────┬────────┘
         │
         ▼
   Transcript Models
   + Read Counts
```

## Related Files

- `src/intron_graph.py` - Graph data structures and construction
- `src/gene_info.py` - `TranscriptModel` class and gene information
- `src/long_read_assigner.py` - Read-to-isoform assignment logic
- `src/long_read_profiles.py` - Profile construction for matching
- `src/isoform_assignment.py` - Assignment types and match classifications

## References

- [IsoQuant Nature Biotechnology Paper](https://www.nature.com/articles/s41587-022-01565-y)
- [IsoQuant GitHub Repository](https://github.com/ablab/IsoQuant)
- [IsoQuant Documentation](https://ablab.github.io/IsoQuant/)
