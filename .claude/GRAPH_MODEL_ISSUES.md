# Graph-Based Model Construction: Known Issues and Future Improvements

This document tracks known issues and planned improvements for the transcript discovery algorithm in `src/graph_based_model_construction.py` and `src/intron_graph.py`.

## Critical Issues (High Priority)

### 1. Terminal Position Handling for Known Isoforms

**Problem**: When a discovered transcript matches a known reference intron chain, the terminal positions (TSS and polyA) are taken from the reference transcript annotation rather than being detected from the actual read evidence.

**Impact**:
- Misses alternative polyadenylation (APA) sites that differ from annotation
- Misses alternative transcription start sites (TSS)
- Reduces accuracy for isoforms with correct splicing but different ends

**Location**: `construct_fl_isoforms()` and `transcript_from_reference()`

**Desired behavior**: Use read-derived terminal positions even when intron chain matches reference.

---

### 2. Suboptimal PolyA/TSS Clustering

**Problem**: The `attach_terminal_positions()` method in `IntronGraph` uses a greedy clustering algorithm for polyA and TSS positions. This leads to suboptimal clustering, especially in regions with multiple close alternative sites.

**Current approach** (`cluster_polya_positions`):
```python
while position_dict:
    best_pair = max(position_dict.items(), key=lambda x:x[1])
    top_position = best_pair[0]
    # Cluster all positions within apa_delta of top_position
    for pos in range(top_position - apa_delta, top_position + apa_delta + 1):
        if pos in position_dict:
            total_count += position_dict[pos]
            del position_dict[pos]
```

**Issues with greedy approach**:
- Order-dependent results
- Can merge distinct peaks that happen to be within `apa_delta`
- No consideration of peak shape or read distribution

**Planned solution**: New algorithm based on peak detection and classification using XGBoost (to be merged from `polya_position` branch).

**Location**: `IntronGraph.cluster_polya_positions()`, `IntronGraph.cluster_terminal_positions()`

---

### 3. Underutilization of Full-Length Read Information

**Problem**: When reads are confirmed full-length (have both polyA and TSS evidence), the algorithm still doesn't fully rely on their terminal positions for transcript model construction.

**Current behavior**:
- FL status is used for path classification (`fl_paths`)
- But terminal positions may be adjusted or overridden
- Reference positions preferred over read-derived positions

**Desired behavior**:
- Trust FL read terminals more strongly
- Use FL read consensus for transcript boundaries
- Weight FL reads higher in terminal position clustering

**Location**: `construct_fl_isoforms()`, `IntronPathStorage.fill()`

---

### 4. Inconsistent Read Assignment Architecture

**Problem**: Reads are assigned to isoforms in multiple different ways throughout the pipeline, creating discrepancies between:
- **Stage 1**: Read assignment (`LongReadAssigner` in main pipeline)
- **Stage 2**: Novel isoform discovery (re-assignment in `assign_reads_to_models()`)

**Symptoms**:
- Multiple assignment passes needed (before and after filtering)
- `add_read_info_raw()` in counter classes - a bypass of normal counter architecture
- Read may be assigned to different isoforms at different stages
- Counts may not match between assignment files and quantification

**Code locations**:
- `assign_reads_to_models()` - creates new assignments
- `forward_counts()` - uses `add_read_info_raw()` to bypass counters
- `save_assigned_read()` - separate tracking from main assignment

**Current workaround** (`forward_counts`):
```python
# Ugly bypass of counter architecture
self.transcript_counter.add_read_info_raw(read_id, [transcript_id], read_group)
self.gene_counter.add_read_info_raw(read_id, [gene_id], read_group)
```

**Desired architecture**:
- Single source of truth for read assignments
- Novel transcript discovery should extend, not replace, initial assignments
- Counters should use unified assignment interface
- No special-case raw methods

---

## Code Quality Issues (Medium Priority)

### 5. Large Monolithic Class

**Problem**: `GraphBasedModelConstructor` is ~900 lines with too many responsibilities:
- Graph construction orchestration
- Path extraction and storage
- Multiple types of transcript construction
- Filtering logic
- Read assignment
- Counter management
- SQANTI output generation

**Suggestion**: Split into focused classes:
- `TranscriptPathFinder` - path extraction
- `TranscriptFilter` - filtering logic
- `TranscriptAssigner` - read assignment to models
- Keep `GraphBasedModelConstructor` as orchestrator only

---

### 6. Class-Level State

**Problem**: `detected_known_isoforms` and `extended_transcript_ids` are class variables:

```python
class GraphBasedModelConstructor:
    detected_known_isoforms = set()  # Shared across all instances!
    extended_transcript_ids = set()
```

**Impact**:
- State persists across genes in the same process
- Can cause incorrect behavior in multi-gene processing
- Makes testing difficult

**Fix**: Move to instance variables, initialize in `__init__()`.

---

### 7. Quadratic Similarity Detection

**Problem**: `detect_similar_isoforms()` has O(n²) complexity:

```python
for model in model_storage:
    transcript_model_gene_info = GeneInfo.from_models([model], ...)
    assigner = LongReadAssigner(transcript_model_gene_info, ...)

    for m in model_storage:  # Compare against all others
        assignment = assigner.assign_to_isoform(...)
```

**Impact**: Slow for genes with many novel isoforms

**Suggestion**:
- Pre-filter by intron chain similarity
- Use hash-based grouping of similar intron chains
- Only compare within groups

---

### 8. Repeated Object Creation

**Problem**: Creates temporary `GeneInfo`, `LongReadAssigner`, and `CombinedProfileConstructor` objects in loops:

- In `detect_similar_isoforms()` - per model
- In `assign_reads_to_models()` - per call

**Impact**: Memory allocation overhead, slower execution

**Suggestion**: Reuse objects where possible, or create once outside loops.

---

## Code Quality Issues (Continued)

### 9. Commented-Out Unit Tests

**Problem**: The entire `tests/test_intron_graph.py` file (374 lines) is commented out:

```python
# class TestIntronGraph(TestCase):
#     def test_basic_graph(self):
#         ...
```

**Impact**:
- No unit test coverage for intron graph construction
- Regressions may go undetected
- Only integration tests via `console_test.py` available

**Location**: `tests/test_intron_graph.py`

**Suggestion**: Restore and update unit tests for critical graph operations.

---

### 10. Assertion-Based Error Handling

**Problem**: Several assertions in clustering code could fail with bad data:

```python
# intron_graph.py:532-534
if read_end:
    assert top_position > intron[1]  # Could fail with strand issues
else:
    assert top_position < intron[0]  # Could fail with coordinate shift
```

**Impact**: AssertionError crashes instead of graceful handling

**Suggestion**: Convert to explicit checks with informative error messages or logging.

---

### 11. Similarity Detection Called Twice

**Problem**: `detect_similar_isoforms()` is called twice in the filtering pipeline:

```python
# Line 308 (in filter_transcripts)
to_substitute = self.detect_similar_isoforms(spliced_models)

# Line 351 (in filter_transcripts, after filtering)
to_substitute = self.detect_similar_isoforms(final_models)
```

**Impact**: Double the O(N²) cost; creates and destroys objects twice.

**Suggestion**: Single pass after final filtering, or incremental updates.

---

## Minor Issues (Low Priority)

### 12. Hardcoded Thresholds

**Problem**: Many magic numbers throughout the code:

```python
if len(internal_count_values) > 50:
    coverage_cutoff = internal_count_values[50]  # Magic number

if self.scores[best_gene_pair] < 0.1:  # Magic threshold
    break
```

**Suggestion**: Move to configuration parameters with clear documentation.

---

### 13. Technical Replicas Check Limitation

**Problem**: Technical replicas check only works with `file_name` grouping:

```python
self.file_name_group_idx = self.grouping_strategy_names.index("file_name")
    if "file_name" in self.grouping_strategy_names else -1

# Later:
if self.use_technical_replicas and self.file_name_group_idx >= 0:
    file_groups = set([a.read_group[self.file_name_group_idx] ...])
```

**Impact**: Cannot use technical replicas check with other grouping strategies.

---

### 14. Incomplete End Correction

**Problem**: `correct_novel_transcript_ends()` only corrects ends inward (shortening), never outward:

```python
if read_start > transcript_start:
    new_transcript_start = read_start  # Only shortens
```

**Impact**: May miss true transcript boundaries supported by minority of reads.

---

### 15. Mixed Strand Handling

**Problem**: Complex and scattered strand detection logic:
- `StrandDetector` class
- PolyA/polyT inference
- Reference strand lookup
- Fallback to '.'

**Suggestion**: Centralize strand logic with clear priority rules.

---

### 16. Monointron Detection Edge Case

**Problem**: `is_monointron()` may fail if vertex not in `outgoing_edges`:

```python
def is_monointron(self, v):
    no_outgoing = not self.outgoing_edges[v] or \
                  (all(is_terminal_vertex(x) for x in self.outgoing_edges[v]))
    # If v not in outgoing_edges, could raise KeyError
```

**Suggestion**: Use `.get(v, set())` for safe access.

---

### 17. Path Dictionary Key Mutation Risk

**Problem**: Paths stored using tuple as dictionary key:

```python
# intron_graph.py
path_tuple = tuple(intron_path)
self.paths[path_tuple] += 1
```

While tuples are immutable, if the source list `intron_path` is reused/modified between calls, could lead to confusion.

**Risk level**: Low (tuples are copied), but worth noting.

---

## Planned Improvements Timeline

### Short Term
- [ ] Fix class-level state issue (#6)
- [ ] Move hardcoded thresholds to parameters (#12)
- [ ] Restore unit tests for intron graph (#9)
- [ ] Fix assertion-based error handling (#10)

### Medium Term
- [ ] Integrate new XGBoost-based peak detection (#2)
- [ ] Improve FL read terminal usage (#3)
- [ ] Refactor similarity detection for better performance (#7, #11)
- [ ] Reduce repeated object creation (#8)
- [ ] Fix incomplete end correction (#14)

### Long Term
- [ ] Redesign read assignment architecture (#4)
- [ ] Use read-derived terminals for known isoforms (#1)
- [ ] Split monolithic class (#5)
- [ ] Centralize strand handling logic (#15)

---

## Issue Summary Table

| # | Issue | Severity | Location |
|---|-------|----------|----------|
| 1 | Terminal positions from reference | Critical | `construct_fl_isoforms()` |
| 2 | Suboptimal polyA/TSS clustering | Critical | `cluster_polya_positions()` |
| 3 | Underutilized FL read info | Critical | `construct_fl_isoforms()` |
| 4 | Inconsistent read assignment | Critical | Multiple locations |
| 5 | Large monolithic class | Medium | `GraphBasedModelConstructor` |
| 6 | Class-level state | Medium | Line 50-51 |
| 7 | Quadratic similarity detection | Medium | `detect_similar_isoforms()` |
| 8 | Repeated object creation | Medium | Multiple loops |
| 9 | Commented-out unit tests | Medium | `test_intron_graph.py` |
| 10 | Assertion-based errors | Medium | `cluster_polya_positions()` |
| 11 | Similarity detection called twice | Medium | `filter_transcripts()` |
| 12 | Hardcoded thresholds | Low | Multiple locations |
| 13 | Technical replicas limitation | Low | `construct_fl_isoforms()` |
| 14 | Incomplete end correction | Low | `correct_novel_transcript_ends()` |
| 15 | Mixed strand handling | Low | Multiple locations |
| 16 | Monointron detection edge case | Low | `is_monointron()` |
| 17 | Path dictionary key mutation risk | Low | `IntronPathStorage` |

---

## Related Documentation

- `.claude/GRAPH_BASED_MODEL_CONSTRUCTION.md` - Algorithm description
- `.claude/CLAUDE.md` - Project overview and development guidelines
