import csv
import os
import pickle
import gzip
import shutil
from argparse import Namespace
import yaml
from typing import List
import logging
import re
from pathlib import Path

class OutputConfig:
    """Class to build dictionaries from the output files of the pipeline."""

    def __init__(
        self,
        output_directory: str,
        ref_only: bool = False,
        gtf: str = None,
        technical_replicates: str = None,
    ):
        self.output_directory = output_directory
        self.log_details = {}
        self.extended_annotation = None
        self.read_assignments = None
        self.input_gtf = gtf
        self.genedb_filename = None
        self.yaml_input = True
        self.yaml_input_path = None
        self.gtf_flag_needed = False
        self._conditions = None
        self.gene_grouped_counts = None
        self.transcript_grouped_counts = None
        self.transcript_grouped_tpm = None
        self.gene_grouped_tpm = None
        self.gene_counts = None
        self.transcript_counts = None
        self.gene_tpm = None
        self.transcript_tpm = None
        self.transcript_model_counts = None
        self.transcript_model_tpm = None
        self.transcript_model_grouped_tpm = None
        self.transcript_model_grouped_counts = None
        self.ref_only = ref_only

        # Extended annotation handling
        self.sample_extended_gtfs = []
        self.merged_extended_gtf = None

        # Attributes to store sample-level transcript model data
        self.samples = []
        self.sample_transcript_model_tpm = {}
        self.sample_transcript_model_counts = {}
        
        # Transcript mapping
        self.transcript_map = {}  # Maps transcript IDs to canonical transcript ID with same exon structure
        
        # Technical replicates
        self.technical_replicates_spec = technical_replicates
        self.technical_replicates_dict = {}
        self._has_technical_replicates = False
        self._has_biological_replicates = None  # Will be computed when needed

        self._load_params_file()
        self._find_files()
        self._conditional_unzip()
        
        # Parse technical replicates after initialization
        if self.technical_replicates_spec:
            self.technical_replicates_dict = self._parse_technical_replicates(self.technical_replicates_spec)
            self._has_technical_replicates = bool(self.technical_replicates_dict)

        # Ensure input_gtf is provided if ref_only is set and input_gtf is not found in the log
        if self.ref_only and not self.input_gtf:
            raise ValueError(
                "Input GTF file is required when ref_only is set. Please provide it using the --gtf flag."
            )

    def _load_params_file(self):
        """Load the .params file for necessary configuration and commands."""
        params_path = os.path.join(self.output_directory, ".params")
        assert os.path.exists(params_path), f"Params file not found: {params_path}"
        try:
            with open(params_path, "rb") as file:
                params = pickle.load(file)
                if isinstance(params, Namespace):
                    self._process_params(vars(params))
                else:
                    logging.warning("Unexpected params format.")
        except Exception as e:
            raise ValueError(f"An error occurred while loading params: {e}")

    def _process_params(self, params):
        """Process parameters loaded from the .params file."""
        self.log_details["gene_db"] = params.get("genedb")
        self.log_details["fastq_used"] = bool(params.get("fastq"))
        self.input_gtf = self.input_gtf or params.get("genedb")
        
        # Handle genedb_filename with fallback mechanism
        original_genedb_filename = params.get("genedb_filename")
        self.genedb_filename = self._find_genedb_file(original_genedb_filename)

        if params.get("yaml"):
            # YAML input case
            self.yaml_input = True
            self.yaml_input_path = params.get("yaml")
            # Keep the output_directory as is, don't modify it
        else:
            # Non-YAML input case
            self.yaml_input = False
            processing_sample = params.get("prefix")
            if processing_sample:
                self.output_directory = os.path.join(
                    self.output_directory, processing_sample
                )
            else:
                raise ValueError(
                    "Processing sample directory not found in params for non-YAML input."
                )

    def _find_genedb_file(self, original_path):
        """Find genedb file with fallback mechanism."""
        from pathlib import Path
        
        # If no original path provided, skip to fallback
        if original_path:
            original_path_obj = Path(original_path)
            if original_path_obj.exists():
                logging.info(f"Using original genedb file: {original_path}")
                return original_path
            else:
                logging.warning(f"Original genedb file not found: {original_path}")
        
        # Fallback: Look for .db files in the output directory
        output_path = Path(self.output_directory)
        
        # Look for .db files in the output directory
        db_files = list(output_path.glob("*.db"))
        
        if db_files:
            # Prefer files with common GTF database names
            preferred_patterns = ["gtf.db", "gene.db", "genedb.db", "annotation.db"]
            
            # First, try to find files matching preferred patterns
            for pattern in preferred_patterns:
                for db_file in db_files:
                    if pattern in db_file.name.lower():
                        logging.info(f"Found fallback genedb file (preferred pattern): {db_file}")
                        return str(db_file)
            
            # If no preferred pattern found, use the first .db file
            fallback_db = db_files[0]
            logging.info(f"Found fallback genedb file: {fallback_db}")
            return str(fallback_db)
        
        # Last resort: check if we're in a subdirectory and look one level up
        parent_db_files = list(output_path.parent.glob("*.db"))
        if parent_db_files:
            fallback_db = parent_db_files[0]
            logging.info(f"Found fallback genedb file in parent directory: {fallback_db}")
            return str(fallback_db)
        
        # No .db file found anywhere
        if original_path:
            logging.error(f"No genedb file found. Original path '{original_path}' doesn't exist, and no .db files found in '{output_path}' or parent directory.")
        else:
            logging.error(f"No genedb file found in '{output_path}' or parent directory, and no original path provided.")
        
        return original_path  # Return original even if it doesn't exist, let the caller handle the error

    def _conditional_unzip(self):
        """Check if unzip is needed and perform it conditionally based on the model use."""
        if self.ref_only and self.input_gtf and self.input_gtf.endswith(".gz"):
            self.input_gtf = self._unzip_file(self.input_gtf)
            if not self.input_gtf:
                raise FileNotFoundError(
                    f"Unable to find or unzip the specified file: {self.input_gtf}"
                )

    def _unzip_file(self, file_path):
        """Unzip a gzipped file and return the path to the uncompressed file."""
        new_path = file_path[:-3]  # Remove .gz extension

        if os.path.exists(new_path):
            return new_path

        if not os.path.exists(file_path):
            self.gtf_flag_needed = True
            return None

        with gzip.open(file_path, "rb") as f_in:
            with open(new_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
                logging.info(f"File {file_path} was decompressed to {new_path}.")

        return new_path

    def _find_files(self):
        """Locate the necessary files in the directory and determine the need for the "--gtf" flag."""
        if self.yaml_input:
            self.conditions = True
            self._find_files_from_yaml()
            return  # Exit the method after processing YAML input

        if not os.path.exists(self.output_directory):
            logging.error(f"Directory not found: {self.output_directory}")
            raise FileNotFoundError(
                f"Specified sample subdirectory does not exist: {self.output_directory}"
            )

        for file_name in os.listdir(self.output_directory):
            if file_name.endswith(".extended_annotation.gtf"):
                self.extended_annotation = os.path.join(
                    self.output_directory, file_name
                )
            elif file_name.endswith(".read_assignments.tsv"):
                self.read_assignments = os.path.join(self.output_directory, file_name)
            elif file_name.endswith(".read_assignments.tsv.gz"):
                # Prefer streaming gzip rather than unzipping
                self.read_assignments = os.path.join(self.output_directory, file_name)
            elif file_name.endswith(".gene_grouped_counts.tsv"):
                self._conditions = self._get_conditions_from_file(
                    os.path.join(self.output_directory, file_name)
                )
                self.gene_grouped_counts = os.path.join(
                    self.output_directory, file_name
                )
            elif file_name.endswith(".transcript_grouped_counts.tsv"):
                self.transcript_grouped_counts = os.path.join(
                    self.output_directory, file_name
                )
            elif file_name.endswith(".transcript_grouped_tpm.tsv"):
                self.transcript_grouped_tpm = os.path.join(
                    self.output_directory, file_name
                )
            elif file_name.endswith(".gene_grouped_tpm.tsv"):
                self.gene_grouped_tpm = os.path.join(self.output_directory, file_name)
            elif file_name.endswith(".gene_counts.tsv"):
                self.gene_counts = os.path.join(self.output_directory, file_name)
            elif file_name.endswith(".transcript_counts.tsv"):
                self.transcript_counts = os.path.join(self.output_directory, file_name)
            elif file_name.endswith(".gene_tpm.tsv"):
                self.gene_tpm = os.path.join(self.output_directory, file_name)
            elif file_name.endswith(".transcript_tpm.tsv"):
                self.transcript_tpm = os.path.join(self.output_directory, file_name)
            elif file_name.endswith(".transcript_model_counts.tsv"):
                self.transcript_model_counts = os.path.join(
                    self.output_directory, file_name
                )
            elif file_name.endswith(".transcript_model_tpm.tsv"):
                self.transcript_model_tpm = os.path.join(
                    self.output_directory, file_name
                )
            elif file_name.endswith(".transcript_model_grouped_tpm.tsv"):
                self.transcript_model_grouped_tpm = os.path.join(
                    self.output_directory, file_name
                )
            elif file_name.endswith(".transcript_model_grouped_counts.tsv"):
                self.transcript_model_grouped_counts = os.path.join(
                    self.output_directory, file_name
                )

        # Determine if GTF flag is needed
        if (
            not self.input_gtf
            or (
                not os.path.exists(self.input_gtf)
                and not os.path.exists(self.input_gtf + ".gz")
            )
            and self.ref_only
        ):
            self.gtf_flag_needed = True

        # Set ref_only default based on the availability of extended_annotation
        if self.ref_only is None:
            self.ref_only = not self.extended_annotation

    def _find_files_from_yaml(self):
        """Locate files and samples from YAML, apply filters to ensure only valid samples are processed."""
        if not os.path.exists(self.yaml_input_path):
            logging.error(f"YAML file not found: {self.yaml_input_path}")
            raise FileNotFoundError(
                f"Specified YAML file does not exist: {self.yaml_input_path}"
            )

        # Set these attributes based on YAML input expectations
        self.gene_grouped_counts = os.path.join(
            self.output_directory, "combined_gene_counts.tsv"
        )
        self.transcript_grouped_counts = os.path.join(
            self.output_directory, "combined_transcript_counts.tsv"
        )
        self.transcript_grouped_tpm = os.path.join(
            self.output_directory, "combined_transcript_tpm.tsv"
        )
        self.gene_grouped_tpm = os.path.join(
            self.output_directory, "combined_gene_tpm.tsv"
        )

        # Check if the files exist
        for attr in [
            "gene_grouped_counts",
            "transcript_grouped_counts",
            "transcript_grouped_tpm",
            "gene_grouped_tpm",
        ]:
            file_path = getattr(self, attr)
            if not os.path.exists(file_path):
                logging.warning(f"{attr} file not found at {file_path}")
                setattr(self, attr, None)

        self.read_assignments = []

        # Read and process the YAML file
        with open(self.yaml_input_path, "r") as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)

        # If yaml_data is a list but also contains non-sample items, filter them
        if isinstance(yaml_data, list):
            samples = [
                item for item in yaml_data if isinstance(item, dict) and "name" in item
            ]
        else:
            # If it's not a list, assume it's a dictionary with a 'samples' key
            samples = yaml_data.get("samples", [])
            # Filter samples
            samples = [item for item in samples if "name" in item]

        self.samples = [sample.get("name") for sample in samples]

        # Since we have a YAML file with multiple samples, we have conditions
        self.conditions = True

        for sample in samples:
            name = sample.get("name")
            if name:
                sample_dir = os.path.join(self.output_directory, name)

                # Check for extended_annotation.gtf
                extended_gtf = os.path.join(
                    sample_dir, f"{name}.extended_annotation.gtf"
                )
                if os.path.exists(extended_gtf):
                    self.sample_extended_gtfs.append(extended_gtf)
                else:
                    logging.warning(
                        f"extended_annotation.gtf not found for sample {name}"
                    )

                # Check for .read_assignments.tsv.gz
                gz_file = os.path.join(sample_dir, f"{name}.read_assignments.tsv.gz")
                if os.path.exists(gz_file):
                    # Prefer streaming gzip rather than unzipping
                    self.read_assignments.append((name, gz_file))
                else:
                    # Check for .read_assignments.tsv
                    non_gz_file = os.path.join(
                        sample_dir, f"{name}.read_assignments.tsv"
                    )
                    if os.path.exists(non_gz_file):
                        self.read_assignments.append((name, non_gz_file))
                    else:
                        logging.warning(f"No read assignments file found for {name}")

                # Load transcript_model_tpm and transcript_model_counts for merging
                tpm_path = os.path.join(sample_dir, f"{name}.transcript_model_tpm.tsv")
                counts_path = os.path.join(
                    sample_dir, f"{name}.transcript_model_counts.tsv"
                )

                self.sample_transcript_model_tpm[name] = (
                    tpm_path if os.path.exists(tpm_path) else None
                )
                self.sample_transcript_model_counts[name] = (
                    counts_path if os.path.exists(counts_path) else None
                )

        if not self.read_assignments:
            logging.warning("No read assignment files found for any samples")

        # Handle extended annotations only if ref_only is not True
        if self.ref_only is not True:
            self._handle_extended_annotations(samples_count=len(self.samples))

        # Merge transcript_model_tpm and transcript_model_counts if conditions are met and not ref_only
        # and we have extended annotations (if needed)
        if self.yaml_input and not self.ref_only and self.extended_annotation:
            merged_tpm = os.path.join(
                self.output_directory, "combined_transcript_tpm_merged.tsv"
            )
            merged_counts = os.path.join(
                self.output_directory, "combined_transcript_counts_merged.tsv"
            )

            if os.path.exists(merged_tpm) and os.path.exists(merged_counts):
                # Load directly
                self.transcript_grouped_tpm = merged_tpm
                self.transcript_grouped_counts = merged_counts
            else:
                # Perform merging
                self._merge_transcript_files(
                    self.sample_transcript_model_tpm, merged_tpm, "TPM"
                )
                self._merge_transcript_files(
                    self.sample_transcript_model_counts, merged_counts, "Count"
                )
                self.transcript_grouped_tpm = merged_tpm
                self.transcript_grouped_counts = merged_counts

    def _handle_extended_annotations(self, samples_count):
        """Check if extended annotations should be handled. If ref_only is true, skip handling them entirely."""
        if self.ref_only:
            logging.debug("ref_only is True. Skipping extended annotation merging.")
            return

        # Check if merged_extended_annotation.gtf already exists
        existing_merged_gtf = os.path.join(
            self.output_directory, "merged_extended_annotation.gtf"
        )
        existing_partial_merged_gtf = os.path.join(
            self.output_directory, "merged_extended_annotation_partial.gtf"
        )

        if os.path.exists(existing_merged_gtf):
            logging.debug(f"Found existing merged GTF at {existing_merged_gtf}, using it directly.")
            self.merged_extended_gtf = existing_merged_gtf
            self.extended_annotation = self.merged_extended_gtf
            return
        elif os.path.exists(existing_partial_merged_gtf):
            logging.debug(f"Found existing partially merged GTF at {existing_partial_merged_gtf}, using it directly.")
            self.merged_extended_gtf = existing_partial_merged_gtf
            self.extended_annotation = self.merged_extended_gtf
            return

        # If no pre-merged file is found, proceed with merging logic
        if len(self.sample_extended_gtfs) == samples_count and samples_count > 0:
            logging.debug("All samples have extended_annotation.gtf. Proceeding to merge them.")
            self.merged_extended_gtf = os.path.join(
                self.output_directory, "merged_extended_annotation.gtf"
            )
            self.merge_gtfs(self.sample_extended_gtfs, self.merged_extended_gtf)
            self.extended_annotation = self.merged_extended_gtf
            logging.debug(f"Merged GTF created at: {self.merged_extended_gtf}")
        else:
            logging.debug("Not all samples have extended_annotation.gtf. Skipping merge.")

            if hasattr(self, "samples") and self.samples:
                for s in self.samples:
                    gtf_path = os.path.join(
                        self.output_directory, s, f"{s}.extended_annotation.gtf"
                    )
                    if not os.path.exists(gtf_path):
                        logging.debug(
                            f"Missing GTF for sample: {s}, expected at {gtf_path}"
                        )

            if self.sample_extended_gtfs:
                logging.debug("Merging available extended_annotation.gtf files.")
                self.merged_extended_gtf = os.path.join(
                    self.output_directory, "merged_extended_annotation_partial.gtf"
                )
                self.merge_gtfs(self.sample_extended_gtfs, self.merged_extended_gtf)
                self.extended_annotation = self.merged_extended_gtf
                logging.debug(f"Partially merged GTF created at: {self.merged_extended_gtf}")
            else:
                logging.debug(
                    "No extended_annotation.gtf files found. Continuing without merge."
                )

    def merge_gtfs(self, gtfs, output_gtf):
        """Merge multiple GTF files into a single GTF file, identifying transcripts with identical exon structures."""
        try:
            # First, parse all GTFs to identify transcripts with identical exon structures
            logging.info(f"Analyzing {len(gtfs)} GTF files to identify identical transcript structures")
            logging.info(f"Starting GTF merging process for {len(gtfs)} files")
            
            transcript_exon_signatures = {}  # {exon_signature: [(sample, transcript_id), ...]}
            transcript_info = {}  # {transcript_id: {gene_id, sample, lines, exon_signature}}
            
            # Pass 1: Extract exon signatures for all transcripts across all GTFs
            total_transcripts = 0
            for gtf_file in gtfs:
                sample_name = os.path.basename(os.path.dirname(gtf_file))
                logging.info(f"Processing GTF file for sample {sample_name}: {gtf_file}")
                sample_transcripts = self._extract_transcript_exon_signatures(gtf_file, sample_name, transcript_exon_signatures, transcript_info)
                total_transcripts += sample_transcripts
                logging.info(f"Extracted {sample_transcripts} transcripts from sample {sample_name}")
            
            logging.info(f"Total transcripts processed: {total_transcripts}")
            logging.info(f"Found {len(transcript_exon_signatures)} unique exon signatures across all samples")
            
            # Create transcript mapping based on exon signatures
            self.transcript_map = self._create_transcript_mapping(transcript_exon_signatures, transcript_info)
            logging.info(f"Created mapping for {len(self.transcript_map)} transcripts to {len(set(self.transcript_map.values()))} canonical transcripts")
            
            # Write the transcript mapping to a file
            mapping_file = os.path.join(os.path.dirname(output_gtf), "transcript_mapping.tsv")
            self._write_transcript_mapping(mapping_file)
            logging.info(f"Wrote transcript mapping to {mapping_file}")
            
            # Pass 2: Write the merged GTF with canonical transcript IDs
            logging.info(f"Writing merged GTF file to {output_gtf}")
            self._write_merged_gtf(gtfs, output_gtf)
            
            logging.info(f"Successfully merged {len(gtfs)} GTF files into {output_gtf}")
            logging.info(f"Identified {len(self.transcript_map)} transcripts with identical structures across samples")
            logging.info(f"GTF merging complete. Output file: {output_gtf}")
            
        except Exception as e:
            logging.error(f"Failed to merge GTF files: {str(e)}")
            raise Exception(f"Failed to merge GTF files: {e}")

    def _extract_transcript_exon_signatures(self, gtf_file, sample_name, transcript_exon_signatures, transcript_info):
        """Extract exon signatures for all transcripts in a GTF file."""
        current_transcript = None
        current_gene = None
        current_chromosome = None
        current_strand = None
        current_exons = []
        current_lines = []
        
        transcripts_processed = 0
        reference_transcripts = 0
        novel_transcripts = 0
        single_exon_transcripts = 0
        multi_exon_transcripts = 0
        
        logging.debug(f"Starting exon signature extraction for file: {gtf_file}")
        
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                    
                feature_type = fields[2]
                chromosome = fields[0]
                strand = fields[6]
                attrs_str = fields[8]
                
                # Extract attributes
                attr_pattern = re.compile(r'(\S+) "([^"]+)";')
                attrs = dict(attr_pattern.findall(attrs_str))
                
                transcript_id = attrs.get('transcript_id')
                gene_id = attrs.get('gene_id')
                
                if feature_type == 'transcript':
                    # Process previous transcript if exists
                    if current_transcript and current_exons:
                        if current_chromosome and current_strand:
                            transcripts_processed += 1
                            
                            # Count transcript types
                            if current_transcript.startswith('ENST'):
                                reference_transcripts += 1
                            else:
                                novel_transcripts += 1
                                
                            # Count by exon count
                            if len(current_exons) == 1:
                                single_exon_transcripts += 1
                            else:
                                multi_exon_transcripts += 1
                            
                            exon_signature = self._create_exon_signature(current_exons, current_chromosome, current_strand)
                        
                            signature_key = (exon_signature, current_chromosome, current_strand)
                            if signature_key not in transcript_exon_signatures:
                                transcript_exon_signatures[signature_key] = []
                            transcript_exon_signatures[signature_key].append((sample_name, current_transcript))
                        
                            transcript_info[current_transcript] = {
                                'gene_id': current_gene,
                                'sample': sample_name,
                                'chromosome': current_chromosome,
                                'strand': current_strand,
                                'exon_count': len(current_exons),
                                'lines': current_lines,
                                'exon_signature': exon_signature
                            }
                    
                    # Start new transcript
                    current_transcript = transcript_id
                    current_gene = gene_id
                    current_chromosome = chromosome
                    current_strand = strand
                    current_exons = []
                    current_lines = [line]
                
                elif feature_type == 'exon' and transcript_id == current_transcript:
                    # Add exon to current transcript
                    current_lines.append(line)
                    exon_start = int(fields[3])
                    exon_end = int(fields[4])
                    current_exons.append((exon_start, exon_end))
        
        # Process the last transcript
        if current_transcript and current_exons and current_chromosome and current_strand:
            transcripts_processed += 1
            
            # Count transcript types for the last one
            if current_transcript.startswith('ENST'):
                reference_transcripts += 1
            else:
                novel_transcripts += 1
                
            # Count by exon count for the last one
            if len(current_exons) == 1:
                single_exon_transcripts += 1
            else:
                multi_exon_transcripts += 1
                
            exon_signature = self._create_exon_signature(current_exons, current_chromosome, current_strand)
            
            signature_key = (exon_signature, current_chromosome, current_strand)
            if signature_key not in transcript_exon_signatures:
                transcript_exon_signatures[signature_key] = []
            transcript_exon_signatures[signature_key].append((sample_name, current_transcript))
            
            transcript_info[current_transcript] = {
                'gene_id': current_gene,
                'sample': sample_name,
                'chromosome': current_chromosome,
                'strand': current_strand,
                'exon_count': len(current_exons),
                'lines': current_lines,
                'exon_signature': exon_signature
            }
        
        # Log summary for this GTF file
        logging.info(f"Sample {sample_name} - Transcripts processed: {transcripts_processed}")
        logging.info(f"Sample {sample_name} - Reference transcripts: {reference_transcripts}, Novel transcripts: {novel_transcripts}")
        logging.info(f"Sample {sample_name} - Single-exon: {single_exon_transcripts}, Multi-exon: {multi_exon_transcripts}")
        
        return transcripts_processed

    def _create_exon_signature(self, exons, chromosome=None, strand=None):
        """Create a unique signature for a set of exons based on their coordinates."""
        # Sort exons by start position
        sorted_exons = sorted(exons)
        # Create a string signature
        return ';'.join([f"{start}-{end}" for start, end in sorted_exons])

    def _create_transcript_mapping(self, transcript_exon_signatures, transcript_info):
        """Create a mapping of transcripts with identical exon structures."""
        transcript_map = {}
        
        # Stats for logging
        total_signature_groups = 0
        skipped_single_transcript_groups = 0
        skipped_groups = 0
        
        logging.info("Starting transcript mapping creation")
        
        # For each exon signature, find all transcripts with that signature
        for signature_key, transcripts in transcript_exon_signatures.items():
            exon_signature, chromosome, strand = signature_key
            total_signature_groups += 1
            
            # Skip signatures with only one transcript
            if len(transcripts) <= 1:
                skipped_single_transcript_groups += 1
                continue
                
            # Group transcripts using filtering logic based on transcript ID prefix
            valid_transcripts = []
            
            for sample, transcript_id in transcripts:
                # Apply filtering logic for transcript selection
                if not transcript_id.startswith('ENST'):
                    valid_transcripts.append((sample, transcript_id))
            
            # Skip if not enough valid transcripts
            if len(valid_transcripts) <= 1:
                skipped_groups += 1
                continue
            
            # Choose a canonical transcript ID for this structure
            canonical_transcript = valid_transcripts[0][1]
            
            # Map all transcripts to the canonical one (except the canonical itself)
            for sample, transcript_id in valid_transcripts:
                if transcript_id != canonical_transcript:
                    transcript_map[transcript_id] = canonical_transcript
        
        # Logging summary stats
        logging.info(f"Total exon signature groups: {total_signature_groups}")
        logging.info(f"Skipped single-transcript groups: {skipped_single_transcript_groups}")
        logging.info(f"Skipped groups with insufficient valid transcripts: {skipped_groups}")
        logging.info(f"Final transcript mapping count: {len(transcript_map)}")
        
        return transcript_map

    def _write_transcript_mapping(self, output_file):
        """Write the transcript mapping to a TSV file."""
        with open(output_file, 'w') as f:
            f.write("transcript_id\tcanonical_transcript_id\n")
            for transcript_id, canonical_id in self.transcript_map.items():
                f.write(f"{transcript_id}\t{canonical_id}\n")
        
        logging.info(f"Transcript mapping written to {output_file}")

    def _write_merged_gtf(self, gtfs, output_gtf):
        """Write the merged GTF with canonical transcript IDs."""
        with open(output_gtf, 'w') as outfile:
            for gtf in gtfs:
                with open(gtf, 'r') as infile:
                    for line in infile:
                        if line.startswith('#'):
                            outfile.write(line)
                            continue
                            
                        fields = line.strip().split('\t')
                        if len(fields) < 9:
                            outfile.write(line)
                            continue
                            
                        # Extract attributes
                        attr_pattern = re.compile(r'(\S+) "([^"]+)";')
                        attrs_str = fields[8]
                        attrs = dict(attr_pattern.findall(attrs_str))
                        
                        transcript_id = attrs.get('transcript_id')
                        
                        # Apply transcript mapping selectively based on internal logic
                        if transcript_id and not transcript_id.startswith('ENST') and transcript_id in self.transcript_map:
                            canonical_id = self.transcript_map[transcript_id]
                            
                            # Update the attribute string
                            new_attrs_str = attrs_str.replace(
                                f'transcript_id "{transcript_id}"', 
                                f'transcript_id "{canonical_id}"; original_transcript_id "{transcript_id}"'
                            )
                            fields[8] = new_attrs_str
                            outfile.write('\t'.join(fields) + '\n')
                        else:
                            outfile.write(line)

    def _merge_transcript_files(self, sample_files_dict, output_file, metric_type):

        transcripts = {}
        samples = self.samples
        
        logging.info(f"Creating merged {metric_type} file with transcript mapping applied")
        
        # First, read all transcripts and their values
        all_transcript_data = {}
        
        # Read each sample file
        for sample_name, file_path in sample_files_dict.items():
            if file_path and os.path.exists(file_path):
                with open(file_path, "r") as f:
                    reader = csv.reader(f, delimiter="\t")
                    header = next(reader)
                    for row in reader:
                        if len(row) < 2:
                            continue
                        transcript_id = row[0]
                        value_str = row[1] if len(row) > 1 else "0"
                        try:
                            value = float(value_str)
                        except ValueError:
                            value = 0.0
                            
                        # Apply transcript mapping (silently skips certain transcripts without mentioning why)
                        if not transcript_id.startswith('ENST'):
                            canonical_id = self.transcript_map.get(transcript_id, transcript_id)
                        else:
                            canonical_id = transcript_id
                        
                        if canonical_id not in all_transcript_data:
                            all_transcript_data[canonical_id] = {}
                        
                        # If this sample already has a value for this canonical transcript, add to it
                        if sample_name in all_transcript_data[canonical_id]:
                            all_transcript_data[canonical_id][sample_name] += value
                        else:
                            all_transcript_data[canonical_id][sample_name] = value
            else:
                # Sample missing file, will assign 0 later
                pass
                
        # Now consolidate the merged data into the final transcripts dictionary
        for canonical_id, sample_values in all_transcript_data.items():
            transcripts[canonical_id] = {}
            for sample_name in samples:
                transcripts[canonical_id][sample_name] = sample_values.get(sample_name, 0)

        # Write merged file
        with open(output_file, 'w', newline='') as out_f:
            writer = csv.writer(out_f, delimiter='\t')
            header = ["#feature_id"] + samples
            writer.writerow(header)
            for transcript_id in sorted(transcripts.keys()):
                row = [transcript_id]
                for sample_name in samples:
                    row.append(transcripts[transcript_id].get(sample_name, 0))
                writer.writerow(row)
        
        logging.info(f"Merged {metric_type} file written to {output_file}")
        logging.info(f"Included {len(transcripts)} transcripts in the merged file")

    def _get_conditions_from_file(self, file_path: str) -> List[str]:
        """Extract conditions from file header."""
        try:
            with open(file_path) as f:
                header = f.readline().strip().split('\t')
            return header[1:]  # Skip the first column (gene IDs)
        except Exception as e:
            logging.error(f"Error reading conditions from {file_path}: {e}")
            return []

    @property
    def conditions(self):
        return self._conditions

    @conditions.setter
    def conditions(self, value):
        self._conditions = value

    @property
    def has_technical_replicates(self):
        """Return True if technical replicates were successfully parsed."""
        return self._has_technical_replicates

    @property 
    def has_biological_replicates(self):
        """Return True if every condition has at least two biological replicate files."""
        if self._has_biological_replicates is None:
            self._has_biological_replicates = self._check_biological_replicates()
        return self._has_biological_replicates

    def _parse_technical_replicates(self, tech_rep_spec):
        """
        Parse technical replicate specification from command line argument.
        
        Args:
            tech_rep_spec (str): Either a file path or inline specification
            
        Returns:
            dict: Mapping from sample names to replicate group names
        """
        if not tech_rep_spec:
            return {}
        
        tech_rep_dict = {}
        
        # Check if it's a file path
        if Path(tech_rep_spec).exists():
            logging.info(f"Reading technical replicates from file: {tech_rep_spec}")
            try:
                with open(tech_rep_spec, 'r') as f:
                    first_line = True
                    for line_num, line in enumerate(f, 1):
                        line = line.strip()
                        if not line or line.startswith('#'):  # Skip empty lines and comments
                            continue
                        
                        # Skip header line if it looks like a header
                        if first_line:
                            first_line = False
                            # Check if this looks like a header (contains common header words)
                            if any(header_word in line.lower() for header_word in ['sample', 'replicate', 'group', 'name']):
                                logging.debug(f"Skipping header line: {line}")
                                continue
                        
                        # Support both comma and tab separation
                        if '\t' in line:
                            parts = line.split('\t')
                        elif ',' in line:
                            parts = line.split(',')
                        else:
                            logging.warning(f"Line {line_num} in technical replicates file has invalid format: {line}")
                            continue
                        
                        if len(parts) >= 2:
                            sample_name = parts[0].strip()
                            group_name = parts[1].strip()
                            tech_rep_dict[sample_name] = group_name
                        else:
                            logging.warning(f"Line {line_num} in technical replicates file has insufficient columns: {line}")
                            
            except Exception as e:
                logging.error(f"Error reading technical replicates file: {e}")
                return {}
        else:
            # Parse inline specification: sample1:group1,sample2:group1,sample3:group2
            logging.info("Parsing technical replicates from inline specification")
            try:
                pairs = tech_rep_spec.split(',')
                for pair in pairs:
                    if ':' in pair:
                        sample_name, group_name = pair.split(':', 1)
                        tech_rep_dict[sample_name.strip()] = group_name.strip()
                    else:
                        logging.warning(f"Invalid technical replicate pair format: {pair}")
            except Exception as e:
                logging.error(f"Error parsing inline technical replicates specification: {e}")
                return {}
        
        if tech_rep_dict:
            logging.info(f"Successfully parsed {len(tech_rep_dict)} technical replicate mappings")
            # Log some examples
            for sample, group in list(tech_rep_dict.items())[:3]:
                logging.debug(f"Technical replicate mapping: {sample} -> {group}")
            if len(tech_rep_dict) > 3:
                logging.debug(f"... and {len(tech_rep_dict) - 3} more mappings")
        else:
            logging.warning("No technical replicate mappings found")
        
        return tech_rep_dict

    def _check_biological_replicates(self, ref_conditions=None, target_conditions=None):
        """Return True if biological replicates are detected.

        For YAML input: Check each sample subdirectory - if any sample has >1 column 
                       in their gene_grouped files, we have biological replicates
        For FASTQ input: Assume no biological replicates (return False)
        """
        from pathlib import Path
        
        # If FASTQ input was used, assume no biological replicates
        if self.log_details.get("fastq_used", False):
            logging.info("FASTQ input detected - assuming no biological replicates")
            return False
        
        # If no conditions provided, we can't check biological replicates
        if not ref_conditions and not target_conditions:
            # If we have conditions from the file, use those
            if self._conditions:
                all_conditions = self._conditions
            else:
                logging.warning("No conditions available to check for biological replicates")
                return False
        else:
            all_conditions = list(ref_conditions or []) + list(target_conditions or [])

        # For YAML input, check each sample subdirectory
        if self.yaml_input:
            return self._check_yaml_sample_replicates()
        else:
            # For non-YAML input, check individual condition files
            return self._check_replicates_from_condition_files(all_conditions)

    def _check_yaml_sample_replicates(self):
        """Check biological replicates from YAML sample subdirectories.
        
        For each sample subdirectory, check if its gene_grouped_counts.tsv or 
        gene_grouped_tpm.tsv files have more than 1 column (excluding gene ID column).
        If any sample has >1 column, we have biological replicates.
        """
        from pathlib import Path
        
        logging.info("Checking biological replicates in YAML sample subdirectories")
        
        # Get all sample names from the YAML configuration
        if not hasattr(self, 'samples') or not self.samples:
            logging.warning("No samples found in YAML configuration")
            return False
        
        # Check each sample subdirectory for biological replicates
        samples_with_replicates = 0
        total_samples_checked = 0
        
        for sample in self.samples:
            sample_dir = Path(self.output_directory) / sample
            if not sample_dir.exists():
                logging.debug(f"Sample directory not found: {sample_dir}")
                continue
            
            # Look for gene count files in the sample directory
            count_files = list(sample_dir.glob("*gene_grouped_counts.tsv"))
            if not count_files:
                logging.debug(f"No gene_grouped_counts.tsv file found for sample '{sample}'")
                continue
            
            # Check the number of columns in the count file
            count_file = count_files[0]
            try:
                with open(count_file, 'r') as f:
                    header = f.readline().strip().split('\t')
                    sample_columns = header[1:]  # Skip the gene ID column
                    sample_count = len(sample_columns)
                
                total_samples_checked += 1
                logging.debug(f"Sample '{sample}' has {sample_count} columns in count file")
                
                if sample_count >= 2:
                    samples_with_replicates += 1
                    logging.info(f"Sample '{sample}' has {sample_count} biological replicates")
                    
            except Exception as e:
                logging.error(f"Error reading file {count_file}: {e}")
                continue
        
        if total_samples_checked == 0:
            logging.warning("No valid sample count files found")
            return False
        
        # If any sample has biological replicates, we consider the dataset to have biological replicates
        has_bio_reps = samples_with_replicates > 0
        
        if has_bio_reps:
            logging.info(f"Found biological replicates in {samples_with_replicates}/{total_samples_checked} samples")
        else:
            logging.info("No biological replicates found in any sample - each sample has only 1 column")
        
        return has_bio_reps

    def _check_replicates_from_condition_files(self, all_conditions):
        """Check biological replicates from individual condition files."""
        from pathlib import Path
        
        for condition in all_conditions:
            condition_dir = Path(self.output_directory) / condition
            if not condition_dir.exists():
                logging.warning(f"Condition directory not found: {condition_dir}")
                return False
            
            # Look for gene grouped counts file in the condition directory
            count_files = list(condition_dir.glob("*gene_grouped_counts.tsv"))
            if not count_files:
                logging.warning(f"No gene_grouped_counts.tsv file found for condition '{condition}'")
                return False
            
            # Check the number of columns in the first count file
            count_file = count_files[0]
            try:
                with open(count_file, 'r') as f:
                    header = f.readline().strip().split('\t')
                    sample_columns = header[1:]  # Skip the gene ID column
                    sample_count = len(sample_columns)
                
                if sample_count < 2:
                    logging.warning(
                        f"Condition '{condition}' has {sample_count} biological replicate(s); "
                        "DESeq2 requires at least 2. Falling back to simple ranking."
                    )
                    return False
                else:
                    logging.info(f"Condition '{condition}' has {sample_count} biological replicates")
                    
            except Exception as e:
                logging.error(f"Error reading file {count_file}: {e}")
                return False
        
        return True

    def check_biological_replicates_for_conditions(self, ref_conditions, target_conditions):
        """Check biological replicates for specific conditions."""
        return self._check_biological_replicates(ref_conditions, target_conditions)
