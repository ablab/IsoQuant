import csv
import os
import pickle
import gzip
import shutil
from argparse import Namespace
import gffutils
import yaml
from typing import List
import logging

class OutputConfig:
    """Class to build dictionaries from the output files of the pipeline."""

    def __init__(
        self,
        output_directory: str,
        use_counts: bool = False,
        ref_only: bool = False,
        gtf: str = None,
    ):
        self.output_directory = output_directory
        self.log_details = {}
        self.extended_annotation = None
        self.read_assignments = None
        self.input_gtf = gtf  # Initialize with the provided gtf flag
        self.genedb_filename = None
        self.yaml_input = True
        self.yaml_input_path = None
        self.gtf_flag_needed = False  # Initialize flag to check if "--gtf" is needed.
        self._conditions = None  # Changed from self.conditions = False
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
        self.use_counts = use_counts
        self.ref_only = ref_only

        # New attributes for handling extended annotations
        self.sample_extended_gtfs = []
        self.merged_extended_gtf = None

        # Attributes to store sample-level transcript model data
        self.samples = []
        self.sample_transcript_model_tpm = {}
        self.sample_transcript_model_counts = {}

        self._load_params_file()
        self._find_files()
        self._conditional_unzip()

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
                    print("Unexpected params format.")
        except Exception as e:
            raise ValueError(f"An error occurred while loading params: {e}")

    def _process_params(self, params):
        """Process parameters loaded from the .params file."""
        self.log_details["gene_db"] = params.get("genedb")
        self.log_details["fastq_used"] = bool(params.get("fastq"))
        self.input_gtf = self.input_gtf or params.get("genedb")
        self.genedb_filename = params.get("genedb_filename")

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
                print(f"File {file_path} was decompressed to {new_path}.")

        return new_path

    def _find_files(self):
        """Locate the necessary files in the directory and determine the need for the "--gtf" flag."""
        if self.yaml_input:
            self.conditions = True
            self._find_files_from_yaml()
            return  # Exit the method after processing YAML input

        if not os.path.exists(self.output_directory):
            print(f"Directory not found: {self.output_directory}")  # Debugging output
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
                self.read_assignments = self._unzip_file(
                    os.path.join(self.output_directory, file_name)
                )
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
            print(f"YAML file not found: {self.yaml_input_path}")
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
                print(f"Warning: {attr} file not found at {file_path}")
                setattr(self, attr, None)

        # Initialize read_assignments list
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
                    print(
                        f"Warning: extended_annotation.gtf not found for sample {name}"
                    )

                # Check for .read_assignments.tsv.gz
                gz_file = os.path.join(sample_dir, f"{name}.read_assignments.tsv.gz")
                if os.path.exists(gz_file):
                    unzipped_file = self._unzip_file(gz_file)
                    if unzipped_file:
                        self.read_assignments.append((name, unzipped_file))
                    else:
                        print(f"Warning: Failed to unzip {gz_file}")
                else:
                    # Check for .read_assignments.tsv
                    non_gz_file = os.path.join(
                        sample_dir, f"{name}.read_assignments.tsv"
                    )
                    if os.path.exists(non_gz_file):
                        self.read_assignments.append((name, non_gz_file))
                    else:
                        print(f"Warning: No read assignments file found for {name}")

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
            print("Warning: No read assignment files found for any samples")

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

    def _merge_transcript_files(self, sample_files_dict, output_file, metric_type):
        # sample_files_dict: {sample_name: filepath or None}
        # Merge logic:
        # 1. Gather all transcripts from all samples
        # 2. For each transcript, write a line with transcript_id and values from each sample (0 if missing)
        transcripts = {}
        samples = self.samples

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
                        if transcript_id not in transcripts:
                            transcripts[transcript_id] = {}
                        transcripts[transcript_id][sample_name] = value
            else:
                # Sample missing file, will assign 0 later
                pass

        # Write merged file
        with open(output_file, "w", newline="") as out_f:
            writer = csv.writer(out_f, delimiter="\t")
            header = ["#feature_id"] + samples
            writer.writerow(header)
            for transcript_id in sorted(transcripts.keys()):
                row = [transcript_id]
                for sample_name in samples:
                    row.append(transcripts[transcript_id].get(sample_name, 0))
                writer.writerow(row)

    def merge_gtfs(self, gtfs, output_gtf):
        """Merge multiple GTF files into a single GTF file."""
        try:
            with open(output_gtf, "w") as outfile:
                for gtf in gtfs:
                    with open(gtf, "r") as infile:
                        shutil.copyfileobj(infile, outfile)
            print(f"Successfully merged {len(gtfs)} GTF files into {output_gtf}")
        except Exception as e:
            raise Exception(f"Failed to merge GTF files: {e}")

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
