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
        self.logger.setLevel(logging.INFO)
        
        # Initialize sets to store novel gene and transcript IDs
        self.novel_gene_ids = set()
        self.novel_transcript_ids = set()

        # Clean up old cache files on init
        cleanup_cache(self.cache_dir, max_age_days=7)

    def build_gene_dict_with_expression_and_filter(
        self,
        min_value: float = 1.0,
        reference_conditions: List[str] = None,
        target_conditions: List[str] = None,
    ) -> Dict[str, Any]:
        """
        Optimized build process with filtering based on selected conditions.
        Filters transcripts based on min_value occurring in at least one of the
        selected reference_conditions or target_conditions.
        Caches the resulting dictionary based on the specific conditions used.
        """
        self.logger.debug(f"Starting dictionary build: min_value={min_value}, ref={reference_conditions}, target={target_conditions}")

        # 1. Load full TPM matrix to determine available conditions first
        tpm_file = self._get_tpm_file()
        self.logger.debug(f"Loading full TPM matrix from {tpm_file}")
        try:
            tpm_df = pd.read_csv(tpm_file, sep='\t', comment=None)
            tpm_df.columns = [col.lstrip('#') for col in tpm_df.columns]  # Clean headers
            tpm_df = tpm_df.set_index('feature_id')  # Use cleaned column name
        except KeyError as e:
            self.logger.error(f"Missing required column ('feature_id' or condition name) in {tpm_file}: {str(e)}")
            raise
        except Exception as e:
            self.logger.error(f"Failed to load TPM expression matrix: {str(e)}")
            raise

        available_conditions = sorted(tpm_df.columns.tolist()) # Sort for consistent cache key
        self.logger.debug(f"Available conditions in TPM file: {available_conditions}")

        # 2. Determine the actual conditions to process and create a cache key
        requested_conditions = set(reference_conditions or []) | set(target_conditions or [])
        if requested_conditions:
            conditions_to_process = sorted(list(requested_conditions.intersection(available_conditions)))
            missing_conditions = requested_conditions.difference(available_conditions)
            if missing_conditions:
                self.logger.warning(f"Requested conditions not found in TPM file and will be ignored: {missing_conditions}")
            if not conditions_to_process:
                self.logger.error("None of the requested conditions were found in the TPM file. Cannot proceed.")
                return {}
            self.logger.info(f"Processing conditions: {conditions_to_process}")
        else:
            self.logger.info("No specific conditions requested, processing all available conditions.")
            conditions_to_process = available_conditions # Already sorted

        # Create a deterministic cache key based on conditions
        condition_key_part = "_".join(c.replace(" ", "_") for c in conditions_to_process)
        if len(condition_key_part) > 50: # Avoid excessively long filenames
             condition_key_part = f"hash_{hash(condition_key_part)}"

        # 3. Check cache specific to these conditions and min_value
        base_cache_file = build_gene_dict_cache_file( # Keep base name generation consistent
            self.config.extended_annotation,
            tpm_file, # Use original TPM file path for base name consistency
            self.config.ref_only,
            self.cache_dir,
        )
        # Append condition and min_value specifics
        condition_specific_cache_file = base_cache_file.parent / (
            f"{base_cache_file.stem}_conditions_{condition_key_part}_minval_{min_value}.pkl"
        )
        self.logger.debug(f"Looking for cache file: {condition_specific_cache_file}")

        if condition_specific_cache_file.exists():
            self.logger.info(f"Loading data from cache: {condition_specific_cache_file}")
            cached_data = load_cache(condition_specific_cache_file)
            # Expecting (dict, novel_genes_set, novel_transcripts_set)
            if cached_data and isinstance(cached_data, tuple) and len(cached_data) == 3:
                cached_gene_dict, cached_novel_gene_ids, cached_novel_transcript_ids = cached_data
                # Basic validation - check if it's a dict and has expected top-level keys (conditions)
                if isinstance(cached_gene_dict, dict) and all(c in cached_gene_dict for c in conditions_to_process):
                    # Deeper validation might be needed if structure is complex
                    if validate_gene_dict(cached_gene_dict): # Reuse existing validation if suitable
                        self.novel_gene_ids = cached_novel_gene_ids
                        self.novel_transcript_ids = cached_novel_transcript_ids
                        self.logger.info("Successfully loaded dictionary from cache.")
                        return cached_gene_dict
                    else:
                         self.logger.warning("Cached dictionary failed validation. Rebuilding.")
                else:
                    self.logger.warning("Cached data format mismatch or missing conditions. Rebuilding.")
            else:
                 self.logger.warning("Cached data is invalid or in old format. Rebuilding.")

        # 4. Cache miss or invalid: Build dictionary from scratch for the specified conditions
        self.logger.info("Cache miss or invalid. Building dictionary from scratch for selected conditions.")

        # Parse GTF and filter novel genes (only needs to be done once)
        self.logger.info("Parsing GTF and filtering novel genes")
        parsed_data = self.parse_gtf()
        self._validate_gene_structure(parsed_data) # Validate base structure
        base_gene_dict = self._filter_novel_genes(parsed_data) # Also populates self.novel_gene_ids etc.

        # Subset TPM matrix to *only* the conditions being processed for filtering
        tpm_df_subset = tpm_df[conditions_to_process]

        # Identify valid transcripts based on max value within the SUBSET conditions
        transcript_max_values_subset = tpm_df_subset.max(axis=1)
        valid_transcripts = set(
            transcript_max_values_subset[transcript_max_values_subset >= min_value].index
        )
        self.logger.info(
            f"Identified {len(valid_transcripts)} transcripts with TPM >= {min_value} "
            f"in at least one of the conditions: {conditions_to_process}"
        )

        # Build the final dictionary, iterating only through conditions_to_process
        final_dict = {}
        for condition in conditions_to_process:
            final_dict[condition] = {}
            condition_tpm_values = tpm_df[condition] # Get expression from the original full df

            for gene_id, gene_info in base_gene_dict.items():
                # Filter transcripts based on valid_transcripts set AND add expression value
                new_transcripts = {
                    tid: {**tinfo, 'value': condition_tpm_values.get(tid, 0)}
                    for tid, tinfo in gene_info['transcripts'].items()
                    if tid in valid_transcripts # Apply the filter here
                }

                # Only add gene if it has at least one valid transcript remaining
                if new_transcripts:
                    final_dict[condition][gene_id] = {
                        **gene_info, # Copy base gene info
                        'transcripts': new_transcripts,
                        'exons': {} # Initialize exons, will be aggregated next
                    }

            # Validate structure for this condition's dictionary
            self._validate_gene_structure(final_dict[condition])

        # Aggregate exon values based on the filtered transcripts in the final_dict
        self.logger.info("Aggregating exon values based on filtered transcript expression.")
        for condition in conditions_to_process:
            for gene_id, gene_info in final_dict[condition].items():
                aggregated_exons = {}
                for transcript_id, transcript_info in gene_info["transcripts"].items():
                    transcript_value = transcript_info.get("value", 0) # TPM value from the filtered transcript
                    for exon in transcript_info.get("_original_exons", transcript_info.get("exons", [])): # Use original exon structure if available
                        exon_id = exon.get("exon_id")
                        if not exon_id: continue
                        if exon_id not in aggregated_exons:
                            aggregated_exons[exon_id] = {
                                "exon_id": exon_id,
                                "start": exon["start"],
                                "end": exon["end"],
                                "number": exon.get("number", "NA"),
                                "value": 0.0, # Initialize aggregate value
                            }
                        aggregated_exons[exon_id]["value"] += transcript_value # Sum transcript TPM
                gene_info["exons"] = aggregated_exons # Assign aggregated exons

        # 5. Save the newly built dictionary to the condition-specific cache
        self.logger.info(f"Saving filtered dictionary to cache: {condition_specific_cache_file}")
        save_cache(
            condition_specific_cache_file,
            (final_dict, self.novel_gene_ids, self.novel_transcript_ids)
        )

        return final_dict

    def _get_tpm_file(self) -> str:
        """Get the appropriate TPM file path from config."""
        if self.config.conditions:  # Check if we have multiple conditions
            # For multi-condition data, prioritize merged files if available
            merged_tpm = self.config.transcript_grouped_tpm
            if merged_tpm and "_merged.tsv" in merged_tpm:
                self.logger.info("Using merged TPM file with transcript deduplication already applied")
                tpm_file = merged_tpm
            else:
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
        
        return tpm_file

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

    def parse_gtf(self) -> Dict[str, Any]:
        """
        Parse GTF file into a dictionary with genes, transcripts, and exons.
        Handles both reference GTF (with gffutils) and extended annotation GTF.
        """
        if self.config.ref_only:
            # Use gffutils for reference GTF (more robust but slower)
            self.logger.info("Parsing reference GTF using gffutils")
            return self._parse_reference_gtf()
        else:
            # Use faster custom parser for extended annotation
            self.logger.info("Parsing extended annotation GTF with custom parser")
            return self._parse_extended_gtf()
            
    def _parse_reference_gtf(self) -> Dict[str, Any]:
        """Parse reference GTF using gffutils"""
        if not self.config.genedb_filename:
            db_path = self.cache_dir / "gtf.db"
            if not db_path.exists():
                self.logger.info(f"Creating GTF database at {db_path}")
                gffutils.create_db(
                    self.config.input_gtf,
                    dbfn=str(db_path),
                    force=True,
                    merge_strategy="create_unique",
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
            gene_name = feature.attributes.get("gene_name", [gene_id])[0] # Default to gene_id if name missing
            gene_biotype = feature.attributes.get("gene_biotype", ["unknown"])[0] # Default to "unknown"

            gene_dict[gene_id] = {
                "chromosome": feature.seqid,
                "start": feature.start,
                "end": feature.end,
                "strand": feature.strand,
                "name": gene_name, # Use updated gene_name
                "biotype": gene_biotype, # Use updated gene_biotype
                "transcripts": {},
            }

        self.logger.info("Processing transcript and exon features")
        for feature in features.values():
            if feature.featuretype == "transcript":
                gene_id = feature.attributes.get("gene_id", [""])[0]
                if gene_id not in gene_dict:
                    continue

                transcript_id = feature.id
                transcript_name = feature.attributes.get("transcript_name", [transcript_id])[0] # Default to transcript_id
                transcript_biotype = feature.attributes.get("transcript_biotype", ["unknown"])[0] # Default to "unknown"
                transcript_tags = feature.attributes.get("tag", [""])[0].split(",") # Get tags

                gene_dict[gene_id]["transcripts"][transcript_id] = {
                    "start": feature.start,
                    "end": feature.end,
                    "name": transcript_name, # Use updated transcript_name
                    "biotype": transcript_biotype, # Use updated transcript_biotype
                    "exons": [],
                    "tags": transcript_tags, # Use updated transcript_tags
                }
            elif feature.featuretype == "exon":
                gene_id = feature.attributes.get("gene_id", [""])[0]
                transcript_id = feature.attributes.get("transcript_id", [""])[0]
                if (
                    gene_id in gene_dict
                    and transcript_id in gene_dict[gene_id]["transcripts"]
                ):
                    exon_number = feature.attributes.get("exon_number", ["1"])[0] # Default to "1"
                    exon_id = feature.attributes.get("exon_id", [""])[0] # Get exon_id

                    gene_dict[gene_id]["transcripts"][transcript_id]["exons"].append(
                        {
                            "exon_id": exon_id, # Use retrieved exon_id
                            "start": feature.start,
                            "end": feature.end,
                            "number": exon_number, # Use updated exon_number
                        }
                    )

        self.logger.info(f"Processed {len(gene_dict)} genes from reference GTF")
        return gene_dict
        
    def _parse_extended_gtf(self) -> Dict[str, Any]:
        """Parse extended annotation GTF with custom parser"""
        base_gene_dict = {}
        self.logger.info("Parsing extended annotation GTF")
        
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
                            "tags": attrs.get("tags", "").split(","),
                            "name": attrs.get("transcript_name", transcript_id),
                        }
                    
                    elif feature_type == "exon" and transcript_id and gene_id:
                        if gene_id in base_gene_dict and transcript_id in base_gene_dict[gene_id]["transcripts"]:
                            exon_info = {
                                "exon_id": attrs.get("exon_id", ""),
                                "start": int(fields[3]),
                                "end": int(fields[4]),
                                "number": attrs.get("exon_number", "1"),
                                "value": 0.0
                            }
                            base_gene_dict[gene_id]["transcripts"][transcript_id]["exons"].append(exon_info)

            self.logger.info(f"Processed {len(base_gene_dict)} genes from extended annotation GTF")
            return base_gene_dict
        except Exception as e:
            self.logger.error(f"GTF parsing failed: {str(e)}")
            raise
            
    # Keep the original functions for backward compatibility, but have them use the new implementation
    def parse_input_gtf(self) -> Dict[str, Any]:
        """
        Parse the reference GTF file using gffutils.
        This is now a wrapper around _parse_reference_gtf for backward compatibility.
        """
        return self._parse_reference_gtf()
        
    def parse_extended_annotation(self) -> Dict[str, Any]:
        """
        Parse extended annotation GTF.
        This is now a wrapper around _parse_extended_gtf for backward compatibility.
        """
        return self._parse_extended_gtf()

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
            self.logger.debug(f"Read {len(gene_list)} genes from {gene_list_path}")
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
