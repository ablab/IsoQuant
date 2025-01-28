import copy
import gffutils
import pandas as pd
import re
import logging
from pathlib import Path
from typing import Dict, Any, List, Union, Tuple

from src.visualization_cache_utils import (
    build_gene_dict_cache_file,
    build_read_assignment_cache_file,
    save_cache,
    load_cache,
    validate_gene_dict,
    validate_read_assignment_data,
    cleanup_cache,
)


class DictionaryBuilder:
    def __init__(self, config):
        self.config = config
        self.cache_dir = Path(config.output_directory) / ".cache"
        self.cache_dir.mkdir(exist_ok=True)
        
        # Set up logger for DictionaryBuilder
        self.logger = logging.getLogger('IsoQuant.visualization.dictionary_builder')
        self.logger.setLevel(logging.DEBUG)
        
        # Initialize sets to store novel gene and transcript IDs
        self.novel_gene_ids = set()
        self.novel_transcript_ids = set()

        # Clean up old cache files on init
        cleanup_cache(self.cache_dir, max_age_days=7)

    def build_gene_dict_with_expression_and_filter(
        self, min_value: float = 1.0
    ) -> Dict[str, Any]:
        """
        Optimized build process with early filtering and combined steps.
        """
        self.logger.debug(f"Starting optimized dictionary build with min_value={min_value}")
        
        # 1. Check cache first
        expr_file, tpm_file = self._get_expression_files()
        base_cache_file = build_gene_dict_cache_file(
            self.config.extended_annotation,
            expr_file,
            self.config.ref_only,
            self.cache_dir,
        )
        expr_filter_cache = base_cache_file.parent / (
            f"{base_cache_file.stem}_with_expr_minval_{min_value}.pkl"
        )

        if expr_filter_cache.exists():
            cached_data = load_cache(expr_filter_cache)
            if cached_data and len(cached_data) == 3: # Check if load_cache returned a tuple
                cached_gene_dict, cached_novel_gene_ids, cached_novel_transcript_ids = cached_data # Unpack tuple
                if validate_gene_dict(cached_gene_dict):
                    self.novel_gene_ids = cached_novel_gene_ids # Restore from cache
                    self.novel_transcript_ids = cached_novel_transcript_ids # Restore from cache
                    return cached_gene_dict
            else: # Handle older cache format (just gene_dict)
                cached_gene_dict = cached_data
                if validate_gene_dict(cached_gene_dict):
                    return cached_gene_dict

        # 2. Filter novel genes from the base gene dict (not per-condition)
        self.logger.info("Parsing GTF and filtering novel genes")
        parsed_data = self.parse_extended_annotation()
        self._validate_gene_structure(parsed_data)
        base_gene_dict = self._filter_novel_genes(parsed_data)

        # Add debug log: Number of genes and transcripts after novel gene filtering
        gene_count_after_novel_filter = len(base_gene_dict)
        transcript_count_after_novel_filter = sum(len(gene_info.get("transcripts", {})) for gene_info in base_gene_dict.values())
        self.logger.debug(f"After novel gene filtering: {gene_count_after_novel_filter} genes, {transcript_count_after_novel_filter} transcripts")

        # 3. Load expression data with consistent header handling
        self.logger.info("Loading expression matrix")
        try:
            expr_df = pd.read_csv(expr_file, sep='\t', comment=None)
            expr_df.columns = [col.lstrip('#') for col in expr_df.columns]  # Clean headers
            expr_df = expr_df.set_index('feature_id')  # Use cleaned column name
        except KeyError as e:
            self.logger.error(f"Missing required column in {expr_file}: {str(e)}")
            raise
        except Exception as e:
            self.logger.error(f"Failed to load expression matrix: {str(e)}")
            raise

        conditions = expr_df.columns.tolist()
        
        # 4. Vectorized processing instead of row-wise iteration
        transcript_max_values = expr_df.max(axis=1)
        valid_transcripts = set(transcript_max_values[transcript_max_values >= min_value].index)

        # Add debug log: Number of valid transcripts after min_value filtering
        valid_transcript_count = len(valid_transcripts)
        self.logger.debug(f"After min_value ({min_value}) filtering: {valid_transcript_count} valid transcripts")
        
        # 5. Single-pass filtering and value updating
        filtered_dict = {}
        for condition in conditions:
            filtered_dict[condition] = {}
            condition_values = expr_df[condition]

            for gene_id, gene_info in base_gene_dict.items():
                new_transcripts = {
                    tid: {**tinfo, 'value': condition_values.get(tid, 0)}
                    for tid, tinfo in gene_info['transcripts'].items()
                    if tid in valid_transcripts
                }

                if new_transcripts:
                    filtered_dict[condition][gene_id] = {
                        **gene_info,
                        'transcripts': new_transcripts
                    }
            self._validate_gene_structure(filtered_dict[condition]) # Validate structure for each condition's gene dict

        save_cache(expr_filter_cache, (filtered_dict, self.novel_gene_ids, self.novel_transcript_ids)) # Save tuple to cache
        return filtered_dict

    def _get_expression_files(self) -> Tuple[str, str]:
        """Get count file for filtering and TPM file for values."""
        # Get counts file path using existing logic
        counts_file = self._get_expression_file()
        
        # Get corresponding TPM file path
        if self.config.conditions:
            tpm_file = self.config.transcript_grouped_tpm
        else:
            if self.config.ref_only:
                tpm_file = self.config.transcript_tpm_ref
            else:
                base_file = self.config.transcript_tpm.replace('.tsv', '')
                tpm_file = f"{base_file}_merged.tsv"
        
        self.logger.debug(f"Selected TPM file: {tpm_file}")
        if not tpm_file or not Path(tpm_file).exists():
            raise FileNotFoundError(f"TPM file {tpm_file} not found")
        
        return counts_file, tpm_file

    def _get_expression_file(self) -> str:
        """Get the appropriate count file path from config."""
        if self.config.conditions:  # Check if we have multiple conditions
            expr_file = self.config.transcript_grouped_counts
        else:
            if self.config.ref_only:
                expr_file = self.config.transcript_counts_ref
            else:
                base_file = self.config.transcript_counts.replace('.tsv', '')
                expr_file = f"{base_file}_merged.tsv"
        
        self.logger.debug(f"Selected count file: {expr_file}")
        if not expr_file or not Path(expr_file).exists():
            raise FileNotFoundError(f"Count file {expr_file} not found")
        return expr_file

    def _filter_transcripts_above_threshold(
        self, gene_dict: Dict[str, Any], min_value: float
    ) -> Dict[str, Any]:
        """Filter transcripts based on expression threshold."""
        self.logger.info(f"Starting transcript filtering with threshold {min_value}")
        
        # Track transcripts and their maximum values across all conditions
        transcript_max_values = {}
        condition_names = list(gene_dict.keys())
        
        # First pass: find maximum value for each transcript across all conditions
        total_transcripts_before = 0
        for condition in condition_names:
            condition_transcripts = sum(len(gene_info.get("transcripts", {})) 
                                      for gene_info in gene_dict[condition].values())
            total_transcripts_before += condition_transcripts
            self.logger.info(f"Condition {condition}: {condition_transcripts} transcripts before filtering")
            
            # Log sample of transcripts (max 2 per condition)
            sample_transcripts = []
            for gene_info in gene_dict[condition].values():
                sample_transcripts.extend(list(gene_info.get("transcripts", {}).keys())[:2])
                if len(sample_transcripts) >= 2:
                    break
            if sample_transcripts:
                self.logger.debug(f"  Sample transcripts in {condition}: {sample_transcripts[:2]}")

        self.logger.info(f"Found {len(transcript_max_values)} unique transcripts across all conditions")
        
        # Sample of transcripts before filtering
        sample_before = list(transcript_max_values.keys())[:5]
        self.logger.debug(f"Sample transcripts before filtering: {sample_before}")
        
        # Build filtered dictionary
        filtered_dict = {}
        kept_transcripts = set()
        
        for tid, max_value in transcript_max_values.items():
            if max_value >= min_value:
                kept_transcripts.add(tid)
        
        self.logger.info(f"Keeping {len(kept_transcripts)} transcripts that meet threshold {min_value}")
        
        # Sample of kept and filtered transcripts
        sample_kept = list(kept_transcripts)[:5]
        sample_filtered = list(set(transcript_max_values.keys()) - kept_transcripts)[:5]
        self.logger.debug(f"Sample kept transcripts: {sample_kept}")
        self.logger.debug(f"Sample filtered transcripts: {sample_filtered}")
        
        # Create filtered dictionary with same structure as input
        for condition in condition_names:
            filtered_dict[condition] = {}
            for gene_id, gene_info in gene_dict[condition].items():
                new_gene_info = copy.deepcopy(gene_info)
                new_transcripts = {}
                
                for tid, tinfo in gene_info.get("transcripts", {}).items():
                    if tid in kept_transcripts:
                        new_transcripts[tid] = tinfo
                
                if new_transcripts:  # Only keep genes that have remaining transcripts
                    new_gene_info["transcripts"] = new_transcripts
                    filtered_dict[condition][gene_id] = new_gene_info

        # Log final statistics
        for condition in condition_names:
            final_count = sum(len(gene_info.get("transcripts", {})) 
                             for gene_info in filtered_dict[condition].values())
            self.logger.debug(f"  {condition}: {final_count} transcripts")

        return filtered_dict

    # ------------------ READ ASSIGNMENT CACHING ------------------

    def build_read_assignment_and_classification_dictionaries(self):
        """
        Index classifications and assignment types from read_assignments.tsv file(s).
        Returns either:
          - (classification_counts, assignment_type_counts) for single-file input, or
          - (classification_counts_dict, assignment_type_counts_dict) for multi-file (YAML) input.
        """
        if not self.config.read_assignments:
            raise FileNotFoundError("No read assignments file(s) found.")

        # 1. Determine cache file
        cache_file = build_read_assignment_cache_file(
            self.config.read_assignments, self.config.ref_only, self.cache_dir
        )

        # 2. Attempt to load from cache
        if cache_file.exists():
            cached_data = load_cache(cache_file)
            if cached_data and validate_read_assignment_data(
                cached_data, self.config.read_assignments
            ):
                self.logger.info("Using cached read assignment data.")
                return self._post_process_cached_data(cached_data)

        # 3. Otherwise, build from scratch
        self.logger.info("Building read assignment data from scratch.")
        if isinstance(self.config.read_assignments, list):
            classification_counts_dict = {}
            assignment_type_counts_dict = {}
            for sample_name, read_assignment_file in self.config.read_assignments:
                c_counts, a_counts = self._process_read_assignment_file(
                    read_assignment_file
                )
                classification_counts_dict[sample_name] = c_counts
                assignment_type_counts_dict[sample_name] = a_counts

            data_to_cache = {
                "classification_counts": classification_counts_dict,
                "assignment_type_counts": assignment_type_counts_dict,
            }
            save_cache(cache_file, data_to_cache)
            return classification_counts_dict, assignment_type_counts_dict
        else:
            classification_counts, assignment_type_counts = (
                self._process_read_assignment_file(self.config.read_assignments)
            )
            data_to_cache = (classification_counts, assignment_type_counts)
            save_cache(cache_file, data_to_cache)
            return classification_counts, assignment_type_counts

    def _post_process_cached_data(self, cached_data):
        """
        Convert cached_data back to the return format
        for build_read_assignment_and_classification_dictionaries().
        """
        if isinstance(self.config.read_assignments, list):
            return (
                cached_data["classification_counts"],
                cached_data["assignment_type_counts"],
            )
        return cached_data  # (classification_counts, assignment_type_counts)

    def _process_read_assignment_file(self, file_path):
        """
        Parse a read_assignment TSV file, returning:
         - classification_counts: dict(classification -> count)
         - assignment_type_counts: dict(assignment type -> count)
        """
        classification_counts = {}
        assignment_type_counts = {}

        with open(file_path, "r") as file:
            # Skip header lines
            for _ in range(3):
                next(file, None)

            for line in file:
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue

                additional_info = parts[-1]
                classification = (
                    additional_info.split("Classification=")[-1].split(";")[0].strip()
                )
                assignment_type = parts[5]

                classification_counts[classification] = (
                    classification_counts.get(classification, 0) + 1
                )
                assignment_type_counts[assignment_type] = (
                    assignment_type_counts.get(assignment_type, 0) + 1
                )

        return classification_counts, assignment_type_counts

    # -------------------- GTF PARSING --------------------

    def parse_input_gtf(self) -> Dict[str, Any]:
        """
        Parse the reference GTF file using gffutils with optimized settings,
        building a dictionary of genes, transcripts, and exons.
        """
        if not self.config.genedb_filename:
            db_path = self.cache_dir / "gtf.db"
            if not db_path.exists():
                self.logger.info(f"Creating GTF database at {db_path}")
                gffutils.create_db(
                    self.config.input_gtf,
                    dbfn=str(db_path),
                    force=True,
                    merge_strategy="create_unique",  # Faster than merge
                    disable_infer_genes=True,
                    disable_infer_transcripts=True,
                    verbose=False,
                )
            self.config.genedb_filename = str(db_path)

        self.logger.info("Opening GTF database")
        db = gffutils.FeatureDB(self.config.genedb_filename)

        # Pre-fetch all features
        self.logger.info("Pre-fetching features from database")
        features = {feature.id: feature for feature in db.all_features()}

        # Build gene -> transcripts -> exons structure
        gene_dict = {}
        self.logger.info("Processing gene features")
        for feature in features.values():
            if feature.featuretype != "gene":
                continue

            gene_id = feature.id
            gene_dict[gene_id] = {
                "chromosome": feature.seqid,
                "start": feature.start,
                "end": feature.end,
                "strand": feature.strand,
                "name": feature.attributes.get("gene_name", [""])[0],
                "biotype": feature.attributes.get("gene_biotype", [""])[0],
                "transcripts": {},
            }

        self.logger.info("Processing transcript and exon features")
        for feature in features.values():
            if feature.featuretype == "transcript":
                gene_id = feature.attributes.get("gene_id", [""])[0]
                if gene_id not in gene_dict:
                    continue

                transcript_id = feature.id
                gene_dict[gene_id]["transcripts"][transcript_id] = {
                    "start": feature.start,
                    "end": feature.end,
                    "name": feature.attributes.get("transcript_name", [""])[0],
                    "biotype": feature.attributes.get("transcript_biotype", [""])[0],
                    "exons": [],
                    "tags": feature.attributes.get("tag", [""])[0].split(","),
                }
            elif feature.featuretype == "exon":
                gene_id = feature.attributes.get("gene_id", [""])[0]
                transcript_id = feature.attributes.get("transcript_id", [""])[0]
                if (
                    gene_id in gene_dict
                    and transcript_id in gene_dict[gene_id]["transcripts"]
                ):
                    gene_dict[gene_id]["transcripts"][transcript_id]["exons"].append(
                        {
                            "exon_id": feature.id,
                            "start": feature.start,
                            "end": feature.end,
                            "number": feature.attributes.get("exon_number", [""])[0],
                        }
                    )

        self.logger.info(f"Processed {len(gene_dict)} genes from GTF")
        return gene_dict

    def parse_extended_annotation(self) -> Dict[str, Any]:
        """Parse merged GTF into base structure without condition info."""
        base_gene_dict = {}
        self.logger.info("Parsing extended annotation GTF (non-ref_only)")
        
        try:
            with open(self.config.extended_annotation, "r") as file:
                attr_pattern = re.compile(r'(\S+) "([^"]+)";')
                
                # First pass: genes and transcripts
                for line in file:
                    if line.startswith("#") or not line.strip():
                        continue
                    
                    fields = line.strip().split("\t")
                    if len(fields) < 9:
                        continue
                    
                    feature_type = fields[2]
                    attrs = dict(attr_pattern.findall(fields[8]))
                    gene_id = attrs.get("gene_id")
                    transcript_id = attrs.get("transcript_id")
                    
                    if feature_type == "gene" and gene_id:
                        if gene_id not in base_gene_dict:
                            base_gene_dict[gene_id] = {
                                "chromosome": fields[0],
                                "start": int(fields[3]),
                                "end": int(fields[4]),
                                "strand": fields[6],
                                "name": attrs.get("gene_name", gene_id),
                                "biotype": attrs.get("gene_biotype", "unknown"),
                                "transcripts": {}
                            }
                    
                    elif feature_type == "transcript" and gene_id and transcript_id:
                        if gene_id not in base_gene_dict:
                            base_gene_dict[gene_id] = {
                                "chromosome": fields[0],
                                "start": int(fields[3]),
                                "end": int(fields[4]),
                                "strand": fields[6],
                                "name": attrs.get("gene_name", gene_id),
                                "biotype": attrs.get("gene_biotype", "unknown"),
                                "transcripts": {}
                            }
                        
                        base_gene_dict[gene_id]["transcripts"][transcript_id] = {
                            "start": int(fields[3]),
                            "end": int(fields[4]),
                            "exons": [],
                            "expression": 0.0,
                            "tags": attrs.get("tags", "").split(","),
                            "name": attrs.get("transcript_name", transcript_id),
                        }
                    
                    elif feature_type == "exon" and transcript_id and gene_id:
                        if gene_id in base_gene_dict and transcript_id in base_gene_dict[gene_id]["transcripts"]:
                            exon_info = {
                                "start": int(fields[3]),
                                "end": int(fields[4]),
                                "number": attrs.get("exon_number", "1")
                            }
                            base_gene_dict[gene_id]["transcripts"][transcript_id]["exons"].append(exon_info)

            self.logger.info(f"Parsed base structure: {len(base_gene_dict)} genes")
            return base_gene_dict
        except Exception as e:
            self.logger.error(f"GTF parsing failed: {str(e)}")
            raise

    def update_transcript_values(self, gene_dict: Dict[str, Any], counts_file: str, tpm_file: str) -> Dict[str, Any]:
        """Update transcript values from TPM file after filtering with counts."""
        # Read counts for filtering
        counts_df = pd.read_csv(counts_file, sep='\t', comment=None)
        counts_df.columns = [col.lstrip('#') for col in counts_df.columns]
        counts_df = counts_df.set_index('feature_id')
        
        # Read TPMs for values
        tpm_df = pd.read_csv(tpm_file, sep='\t', comment=None)
        tpm_df.columns = [col.lstrip('#') for col in tpm_df.columns]
        tpm_df = tpm_df.set_index('feature_id')
        
        # Align indices between counts and TPMs
        common_transcripts = counts_df.index.intersection(tpm_df.index)
        tpm_df = tpm_df.loc[common_transcripts]
        
        # Rest of existing update logic using tpm_df instead of expr_df
        condition_gene_dict = {condition: copy.deepcopy(gene_dict) for condition in tpm_df.columns}
        
        for tid, row in tpm_df.iterrows():
            base_tid = tid.split('.')[0]
            for condition, tpm_value in row.items():
                # Add nested loop to access gene_info
                for gene_id, gene_info in condition_gene_dict[condition].items():
                    if base_tid in gene_info.get('transcripts', {}):
                        gene_info['transcripts'][base_tid]['value'] = float(tpm_value)
        
        return condition_gene_dict

    # -------------------- UPDATES & UTILITIES --------------------

    def update_gene_names(self, gene_dict: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update gene and transcript identifiers to their names, if available,
        while preserving all nested structure.
        """
        try:
            updated_dict = {}
            total_transcripts = 0
            
            for condition, genes in gene_dict.items():
                updated_genes = {}
                condition_transcripts = 0
                
                for gene_id, gene_info in genes.items():
                    new_gene_info = copy.deepcopy(gene_info)
                    
                    # Update gene name
                    if "name" in gene_info and gene_info["name"]:
                        gene_name_upper = gene_info["name"].upper()
                        updated_genes[gene_name_upper] = new_gene_info
                    else:
                        updated_genes[gene_id] = new_gene_info

                    # Count transcripts
                    transcripts = new_gene_info.get("transcripts", {})
                    condition_transcripts += len(transcripts)
                    
                    # Debug sample of transcript structure
                    if gene_id == list(genes.keys())[0]:
                        self.logger.debug(f"Sample gene {gene_id} transcript structure:")
                        for tid in list(transcripts.keys())[:3]:
                            self.logger.debug(f"Transcript {tid}: {transcripts[tid]}")

                total_transcripts += condition_transcripts
                updated_dict[condition] = updated_genes
                self.logger.debug(f"Condition {condition}: {condition_transcripts} transcripts")

            self.logger.info(f"Updated gene names for {len(gene_dict)} conditions")
            self.logger.info(f"Total transcripts in dictionary: {total_transcripts}")
            
            return updated_dict
            
        except Exception as e:
            self.logger.error(f"Error updating gene/transcript names: {e}")
            self.logger.error(f"Dictionary structure before update: {str(type(gene_dict))}")
            raise

    def read_gene_list(self, gene_list_path: Union[str, Path]) -> List[str]:
        """
        Read and parse a plain-text file containing one gene identifier per line.
        Return a list of uppercase gene IDs/names.
        """
        try:
            with open(gene_list_path, "r") as file:
                gene_list = [line.strip().upper() for line in file if line.strip()]
            self.logger.info(f"Read {len(gene_list)} genes from {gene_list_path}")
            return gene_list
        except Exception as e:
            self.logger.error(f"Error reading gene list from {gene_list_path}: {e}")
            raise

    def _filter_novel_genes(self, gene_dict: Dict[str, Any]) -> Dict[str, Any]:
        """Filter out novel genes based on gene ID pattern."""
        self.logger.debug("Starting novel gene filtering")
        filtered_dict = {}
        total_removed_genes = 0
        total_removed_transcripts = 0
        checked_gene_count = 0
        sample_removed = [] # For debug logging

        novel_gene_pattern = r"novel_gene" # Make sure this pattern is correct for your novel gene IDs

        for gene_id, gene_info in gene_dict.items():
            checked_gene_count += 1
            is_novel = bool(re.match(novel_gene_pattern, gene_id))

            if is_novel:
                total_removed_genes += 1
                removed_transcript_count = len(gene_info.get("transcripts", {}))
                total_removed_transcripts += removed_transcript_count

                # Add novel gene ID to the set
                self.novel_gene_ids.add(gene_id)
                # Add novel transcript IDs to the set
                self.novel_transcript_ids.update(gene_info.get("transcripts", {}).keys())


                if len(sample_removed) < 5: # Sample log of removed genes
                    sample_removed.append({
                        'gene_id': gene_id,
                        'transcript_count': removed_transcript_count
                    })
                continue # Skip adding novel genes to filtered_dict

            filtered_dict[gene_id] = gene_info # Keep known genes

        self.logger.info(
            f"Filtered {total_removed_genes} novel genes "
            f"({total_removed_genes/checked_gene_count:.2%} of total) "
            f"and {total_removed_transcripts} associated transcripts"
        )

        if sample_removed:
            sample_output = "\n".join(
                f"- {g['gene_id']}: {g['transcript_count']} transcripts"
                for g in sample_removed
            )
            self.logger.debug(f"Sample removed novel genes:\n{sample_output}")
        else:
            self.logger.warning("No novel genes detected with current filtering pattern")

        return filtered_dict

    def get_novel_feature_ids(self) -> Tuple[set, set]:
        """Return the sets of novel gene and transcript IDs."""
        return self.novel_gene_ids, self.novel_transcript_ids

    def _filter_low_expression_transcripts(self, condition_gene_dict: Dict[str, Any], min_value: float) -> Dict[str, Any]:
        """Filter transcripts based on expression threshold."""
        self.logger.info(f"Starting transcript filtering with threshold {min_value}")
        
        # Track transcripts and their maximum values across all conditions
        transcript_max_values = {}
        condition_names = list(condition_gene_dict.keys())
        
        # First pass: find maximum value for each transcript across all conditions
        total_transcripts_before = 0
        for condition in condition_names:
            condition_transcripts = 0
            for gene_info in condition_gene_dict[condition].values():
                for tid, tinfo in gene_info.get("transcripts", {}).items():
                    current_value = tinfo.get('value', 0)
                    # Update maximum value tracking
                    if tid not in transcript_max_values or current_value > transcript_max_values[tid]:
                        transcript_max_values[tid] = current_value
                    condition_transcripts += 1
            total_transcripts_before += condition_transcripts
            self.logger.info(f"Condition {condition}: {condition_transcripts} transcripts before filtering")
            
            # Log sample of transcripts (max 2 per condition)
            sample_transcripts = []
            for gene_info in condition_gene_dict[condition].values():
                sample_transcripts.extend(list(gene_info.get("transcripts", {}).keys())[:2])
                if len(sample_transcripts) >= 2:
                    break
            if sample_transcripts:
                self.logger.debug(f"  Sample transcripts in {condition}: {sample_transcripts[:2]}")

        self.logger.info(f"Found {len(transcript_max_values)} unique transcripts across all conditions")
        
        # Sample of transcripts before filtering
        sample_before = list(transcript_max_values.keys())[:5]
        self.logger.debug(f"Sample transcripts before filtering: {sample_before}")
        
        # Build filtered dictionary
        filtered_dict = {}
        kept_transcripts = set()
        
        for tid, max_value in transcript_max_values.items():
            if max_value >= min_value:
                kept_transcripts.add(tid)
        
        self.logger.info(f"Keeping {len(kept_transcripts)} transcripts that meet threshold {min_value}")
        
        # Sample of kept and filtered transcripts
        sample_kept = list(kept_transcripts)[:5]
        sample_filtered = list(set(transcript_max_values.keys()) - kept_transcripts)[:5]
        self.logger.debug(f"Sample kept transcripts: {sample_kept}")
        self.logger.debug(f"Sample filtered transcripts: {sample_filtered}")
        
        # Create filtered dictionary with same structure as input
        for condition in condition_names:
            filtered_dict[condition] = {}
            for gene_id, gene_info in condition_gene_dict[condition].items():
                new_gene_info = copy.deepcopy(gene_info)
                new_transcripts = {}
                
                for tid, tinfo in gene_info.get("transcripts", {}).items():
                    if tid in kept_transcripts:
                        new_transcripts[tid] = tinfo
                
                if new_transcripts:  # Only keep genes that have remaining transcripts
                    new_gene_info["transcripts"] = new_transcripts
                    filtered_dict[condition][gene_id] = new_gene_info

        # Log final statistics
        for condition in condition_names:
            final_count = sum(len(gene_info.get("transcripts", {})) 
                             for gene_info in filtered_dict[condition].values())
            self.logger.debug(f"  {condition}: {final_count} transcripts")

        return filtered_dict

    def _batch_update_values(self, gene_dict, expr_df, valid_transcripts):
        """Vectorized value updating for all conditions."""
        return {
            condition: {
                gene_id: {
                    **gene_info,
                    'transcripts': {
                        tid: {**tinfo, 'value': expr_df.at[tid, condition]}
                        for tid, tinfo in gene_info['transcripts'].items()
                        if tid in valid_transcripts
                    }
                }
                for gene_id, gene_info in gene_dict.items()
                if any(tid in valid_transcripts for tid in gene_info['transcripts'])
            }
            for condition in expr_df.columns
        }

    def _validate_gene_structure(self, gene_dict: Dict[str, Any]) -> None:
        """Ensure proper gene-centric structure before condition processing."""
        required_gene_keys = ['chromosome', 'start', 'end', 'strand', 'name', 'biotype', 'transcripts']
        
        for gene_id, gene_info in gene_dict.items():
            # Check gene ID format
            if not isinstance(gene_id, str) or len(gene_id) < 4:
                self.logger.error(f"Invalid gene ID format: {gene_id}")
                raise ValueError("Malformed gene ID structure")
            
            # Check required keys
            missing = [k for k in required_gene_keys if k not in gene_info]
            if missing:
                self.logger.error(f"Gene {gene_id} missing keys: {missing}")
                raise ValueError("Incomplete gene information")
            
            # Check transcripts structure
            transcripts = gene_info.get('transcripts', {})
            if not isinstance(transcripts, dict):
                self.logger.error(f"Invalid transcripts in gene {gene_id} - expected dict")
                raise ValueError("Malformed transcript structure")
