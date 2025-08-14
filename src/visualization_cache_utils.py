import pickle
import logging
import time
from pathlib import Path
from typing import Dict, Any, Optional, Union
import random
import re
import hashlib


def build_gene_dict_cache_file(
    extended_annotation: Optional[str], input_gtf: str, ref_only: bool, cache_dir: Path
) -> Path:
    """
    Generate a gene dictionary cache filename based on:
    - Which annotation file we're using (extended vs. reference GTF).
    - The modification time of that file.
    - The ref_only setting.
    """
    if extended_annotation and not ref_only:
        source_file = Path(extended_annotation)
        source_type = "extended"
    else:
        source_file = Path(input_gtf)
        source_type = "reference"
    mtime = source_file.stat().st_mtime
    cache_name = f"gene_dict_cache_{source_type}_{source_file.name}_{mtime}_ref_only_{ref_only}.pkl"
    return cache_dir / cache_name


def build_read_assignment_cache_file(
    read_assignments: Union[str, list], ref_only: bool, cache_dir: Path
) -> Path:
    """
    Generate a read-assignment cache filename based on:
    - The read assignment file(s).
    - Possibly their modification times.
    - The ref_only setting.
    """
    if isinstance(read_assignments, str):
        source_file = Path(read_assignments)
        mtime = source_file.stat().st_mtime
        cache_name = (
            f"read_assignment_cache_{source_file.name}_{mtime}_ref_only_{ref_only}.pkl"
        )
        return cache_dir / cache_name
    elif isinstance(read_assignments, list):
        # Build a composite name from the multiple input files
        file_info = []
        for sample_name, path_str in read_assignments:
            path_obj = Path(path_str)
            file_info.append(
                f"{sample_name}-{path_obj.name}-{path_obj.stat().st_mtime}"
            )
        composite_name = "_".join(file_info).replace(" ", "_")[:100]
        cache_name = (
            f"read_assignment_cache_multi_{composite_name}_ref_only_{ref_only}.pkl"
        )
        return cache_dir / cache_name
    else:
        return cache_dir / "read_assignment_cache_default.pkl"


def _hash_list(values: list) -> str:
    try:
        s = ",".join(map(str, values))
        m = hashlib.md5()
        m.update(s.encode('utf-8'))
        return m.hexdigest()[:12]
    except Exception:
        # Fallback to length-based signature
        return f"len{len(values)}"


def build_length_effects_cache_file(
    read_assignments: Union[str, list], ref_only: bool, cache_dir: Path, bin_labels: list
) -> Path:
    """
    Cache name for read-length effects aggregates. Includes input files, mtimes, ref_only, and bin label signature.
    """
    bins_sig = _hash_list(bin_labels)
    if isinstance(read_assignments, str):
        source_file = Path(read_assignments)
        mtime = source_file.stat().st_mtime
        cache_name = (
            f"length_effects_cache_{source_file.name}_{mtime}_bins_{bins_sig}_ref_only_{ref_only}.pkl"
        )
        return cache_dir / cache_name
    elif isinstance(read_assignments, list):
        file_info = []
        for sample_name, path_str in read_assignments:
            path_obj = Path(path_str)
            file_info.append(f"{sample_name}-{path_obj.name}-{path_obj.stat().st_mtime}")
        composite = "_".join(file_info).replace(" ", "_")[:100]
        cache_name = (
            f"length_effects_cache_multi_{composite}_bins_{bins_sig}_ref_only_{ref_only}.pkl"
        )
        return cache_dir / cache_name
    else:
        return cache_dir / "length_effects_cache_default.pkl"


def build_length_hist_cache_file(
    read_assignments: Union[str, list], ref_only: bool, cache_dir: Path, bin_edges: list
) -> Path:
    """
    Cache name for read-length histogram. Includes input files, mtimes, ref_only, and bin edges signature.
    """
    edges_sig = _hash_list(bin_edges)
    if isinstance(read_assignments, str):
        source_file = Path(read_assignments)
        mtime = source_file.stat().st_mtime
        cache_name = (
            f"length_hist_cache_{source_file.name}_{mtime}_edges_{edges_sig}_ref_only_{ref_only}.pkl"
        )
        return cache_dir / cache_name
    elif isinstance(read_assignments, list):
        file_info = []
        for sample_name, path_str in read_assignments:
            path_obj = Path(path_str)
            file_info.append(f"{sample_name}-{path_obj.name}-{path_obj.stat().st_mtime}")
        composite = "_".join(file_info).replace(" ", "_")[:100]
        cache_name = (
            f"length_hist_cache_multi_{composite}_edges_{edges_sig}_ref_only_{ref_only}.pkl"
        )
        return cache_dir / cache_name
    else:
        return cache_dir / "length_hist_cache_default.pkl"


def save_cache(cache_file: Path, data_to_cache: Any) -> None:
    """Save data to a cache file using pickle."""
    try:
        with open(cache_file, 'wb') as f:
            pickle.dump(data_to_cache, f, protocol=pickle.HIGHEST_PROTOCOL) # Save the entire tuple
        logging.debug(f"Successfully saved cache to {cache_file}")
    except Exception as e:
        logging.error(f"Error saving cache to {cache_file}: {e}")


def load_cache(cache_file: Path) -> Any:
    """Load data from a cache file."""
    try:
        with open(cache_file, 'rb') as f:
            cached_data = pickle.load(f)
            if isinstance(cached_data, tuple) and len(cached_data) == 3: # Check if it's the new tuple format
                gene_dict, novel_gene_ids, novel_transcript_ids = cached_data # Unpack tuple
                return gene_dict, novel_gene_ids, novel_transcript_ids # Return the tuple
            else: # Handle old cache format (just gene_dict)
                return cached_data # Return just the gene_dict for backward compatibility
    except FileNotFoundError:
        logging.debug(f"Cache file not found: {cache_file}")
        return None  # Indicate cache miss
    except Exception as e:
        logging.error(f"Error loading cache from {cache_file}: {e}")
        return None


def validate_gene_dict(gene_dict: Dict, ref_only: bool = False) -> bool:
    """Enhanced validation with novel gene check."""
    if not gene_dict:
        return False
        
    # Always check for novel genes regardless of ref_only mode
    novel_genes = sum(1 for condition in gene_dict.values() 
                    for gene_id in condition.keys() 
                    if re.match(r"novel_gene", gene_id))
    if novel_genes > 0:
        logging.warning(f"Found {novel_genes} novel genes in cached dictionary. Rebuilding required.")
        return False
            
    # Existing structure validation
    try:
        for condition in gene_dict.values():
            for gene_info in condition.values():
                if not all(k in gene_info for k in ["chromosome", "start", "end", "strand", "transcripts"]):
                    return False
        return True
    except (KeyError, AttributeError):
        return False


def validate_read_assignment_data(
    data: Any, read_assignments: Union[str, list]
) -> bool:
    """
    Validate the structure of the cached read-assignment data.
    """
    try:
        if isinstance(read_assignments, list):
            # Expecting something like:
            # {
            #   "classification_counts": { "sampleA": {...}, "sampleB": {...} },
            #   "assignment_type_counts": { "sampleA": {...}, "sampleB": {...} }
            # }
            if not isinstance(data, dict):
                return False
            if (
                "classification_counts" not in data
                or "assignment_type_counts" not in data
            ):
                return False
            return True
        else:
            # Single file scenario: We expect a 2-tuple (classification_counts, assignment_type_counts)
            if not isinstance(data, (tuple, list)) or len(data) != 2:
                return False
            return True
    except Exception as e:
        logging.error(f"Read-assignment validation error: {e}")
        return False


def validate_length_effects_data(data: Any, expected_bins: Optional[list] = None) -> bool:
    try:
        required = [
            'bins', 'by_bin_assignment', 'by_bin_classification',
            'assignment_keys', 'classification_keys', 'totals'
        ]
        if not isinstance(data, dict):
            return False
        if any(k not in data for k in required):
            return False
        if expected_bins and data.get('bins') != expected_bins:
            return False
        # Basic shape checks
        if not isinstance(data['by_bin_assignment'], dict): return False
        if not isinstance(data['by_bin_classification'], dict): return False
        if not isinstance(data['totals'], dict): return False
        return True
    except Exception as e:
        logging.error(f"Length-effects validation error: {e}")
        return False


def validate_length_hist_data(data: Any, expected_edges: Optional[list] = None) -> bool:
    try:
        if not isinstance(data, dict):
            return False
        if any(k not in data for k in ['edges', 'counts', 'total']):
            return False
        if expected_edges and list(map(int, data.get('edges', []))) != list(map(int, expected_edges)):
            return False
        return True
    except Exception as e:
        logging.error(f"Length-hist validation error: {e}")
        return False


def cleanup_cache(cache_dir: Path, max_age_days: int = 7) -> None:
    """
    Remove cache files older than specified days.
    """
    current_time = time.time()
    for cache_file in cache_dir.glob("*.pkl"):
        file_age_days = (current_time - cache_file.stat().st_mtime) / (24 * 3600)
        if file_age_days > max_age_days:
            try:
                cache_file.unlink()
                logging.info(f"Removed old cache file: {cache_file}")
            except Exception as e:
                logging.warning(f"Failed to remove cache file {cache_file}: {e}")
