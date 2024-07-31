import csv
import os
import pickle
import gzip
import shutil
import copy
import json
from argparse import Namespace
import tempfile
import gffutils


class OutputConfig:
    """Class to build dictionaries from the output files of the pipeline."""

    def __init__(self, output_directory, use_counts=False, ref_only=None, gtf=None):
        self.output_directory = output_directory
        self.log_details = {}
        self.extended_annotation = None
        self.read_assignments = None
        self.input_gtf = gtf  # Initialize with the provided gtf flag
        self.gtf_flag_needed = False  # Initialize flag to check if "--gtf" is needed.
        self.conditions = False
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

        self._load_params_file()  # Load the params file instead of parsing the log
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

        processing_sample = params.get("prefix")
        if processing_sample:
            self.output_directory = os.path.join(
                self.output_directory, processing_sample
            )
        else:
            raise ValueError("Processing sample directory not found in params.")

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
            # print(f"File {new_path} already exists, using this file.")
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
        if not os.path.exists(self.output_directory):
            print(f"Directory not found: {self.output_directory}")  # Debugging output
            raise FileNotFoundError("Specified sample subdirectory does not exist.")

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
                self.conditions = True
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
            or not os.path.exists(self.input_gtf)
            and not os.path.exists(self.input_gtf + ".gz")
            and self.ref_only
        ):
            self.gtf_flag_needed = True

        # Set ref_only default based on the availability of extended_annotation
        if self.ref_only is None:
            self.ref_only = not self.extended_annotation


class DictionaryBuilder:
    """Class to build dictionaries from the output files of the pipeline."""

    def __init__(self, config):
        self.config = config

    def build_gene_transcript_exon_dictionaries(self):
        """Builds dictionaries of genes, transcripts, and exons from the GTF file."""
        if self.config.extended_annotation and not self.config.ref_only:
            return self.parse_extended_annotation()
        else:
            return self.parse_input_gtf()

    def build_read_assignment_and_classification_dictionaries(self):
        """Indexes classifications and assignment types from the read_assignments.tsv."""
        classification_counts = {}
        assignment_type_counts = {}
        if not self.config.read_assignments:
            raise FileNotFoundError("Read assignments file is missing.")

        with open(self.config.read_assignments, "r") as file:
            next(file)
            next(file)
            next(file)
            for line in file:
                parts = line.split("\t")
                if len(parts) < 6:
                    continue
                additional_info = parts[-1]
                classification = (
                    additional_info.split("Classification=")[-1]
                    .replace(";", "")
                    .strip()
                )
                assignment_type = parts[5]

                classification_counts[classification] = (
                    classification_counts.get(classification, 0) + 1
                )
                assignment_type_counts[assignment_type] = (
                    assignment_type_counts.get(assignment_type, 0) + 1
                )

        return classification_counts, assignment_type_counts

    def parse_input_gtf(self):
        """Parses the GTF file using gffutils to build a detailed dictionary of genes, transcripts, and exons."""
        gene_dict = {}
        if not self.config.input_gtf:
            raise FileNotFoundError("Extended annotation GTF file is missing.")

        input_gtf_path = self.config.input_gtf

        try:
            # Create a temporary database
            with tempfile.NamedTemporaryFile(suffix=".db") as tmp:
                db = gffutils.create_db(
                    input_gtf_path,
                    dbfn=tmp.name,
                    force=True,
                    keep_order=True,
                    merge_strategy="merge",
                    sort_attribute_values=True,
                    disable_infer_genes=True,
                    disable_infer_transcripts=True,
                )

                for gene in db.features_of_type("gene"):
                    gene_id = gene.id
                    gene_dict[gene_id] = {
                        "chromosome": gene.seqid,
                        "start": gene.start,
                        "end": gene.end,
                        "strand": gene.strand,
                        "name": gene.attributes.get("gene_name", [""])[0],
                        "biotype": gene.attributes.get("gene_biotype", [""])[0],
                        "transcripts": {},
                    }

                    for transcript in db.children(gene, featuretype="transcript"):
                        transcript_id = transcript.id
                        gene_dict[gene_id]["transcripts"][transcript_id] = {
                            "start": transcript.start,
                            "end": transcript.end,
                            "name": transcript.attributes.get("transcript_name", [""])[
                                0
                            ],
                            "biotype": transcript.attributes.get(
                                "transcript_biotype", [""]
                            )[0],
                            "exons": [],
                            "tags": transcript.attributes.get("tag", [""])[0].split(
                                ","
                            ),
                        }

                        for exon in db.children(transcript, featuretype="exon"):
                            exon_info = {
                                "exon_id": exon.id,
                                "start": exon.start,
                                "end": exon.end,
                                "number": exon.attributes.get("exon_number", [""])[0],
                            }
                            gene_dict[gene_id]["transcripts"][transcript_id][
                                "exons"
                            ].append(exon_info)

        except Exception as e:
            raise Exception(f"Error parsing GTF file: {str(e)}")

        return gene_dict

    def parse_extended_annotation(self):
        """Parses the GTF file to build a detailed dictionary of genes, transcripts, and exons."""
        gene_dict = {}
        if not self.config.extended_annotation:
            raise FileNotFoundError("Extended annotation GTF file is missing.")

        with open(self.config.extended_annotation, "r") as file:
            for line in file:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    print(
                        f"Skipping malformed line due to insufficient fields: {line.strip()}"
                    )
                    continue

                info_fields = fields[8].strip(";").split(";")
                details = {
                    field.strip().split(" ")[0]: field.strip().split(" ")[1].strip('"')
                    for field in info_fields
                    if " " in field
                }

                try:
                    if fields[2] == "gene":
                        gene_id = details["gene_id"]
                        gene_dict[gene_id] = {
                            "chromosome": fields[0],
                            "start": int(fields[3]),
                            "end": int(fields[4]),
                            "strand": fields[6],
                            "name": details.get("gene_name", ""),
                            "biotype": details.get("gene_biotype", ""),
                            "transcripts": {},
                        }
                    elif fields[2] == "transcript":
                        transcript_id = details["transcript_id"]
                        gene_dict[details["gene_id"]]["transcripts"][transcript_id] = {
                            "start": int(fields[3]),
                            "end": int(fields[4]),
                            "exons": [],
                        }
                    elif fields[2] == "exon":
                        transcript_id = details["transcript_id"]
                        exon_info = {
                            "exon_id": details["exon_id"],
                            "start": int(fields[3]),
                            "end": int(fields[4]),
                        }
                        gene_dict[details["gene_id"]]["transcripts"][transcript_id][
                            "exons"
                        ].append(exon_info)
                except KeyError as e:
                    print(f"Key error in line: {line.strip()} | Missing key: {e}")
        return gene_dict

    def update_gene_dict(self, gene_dict, value_df):
        new_dict = {}
        gene_values = {}

        # Read gene counts from value_df
        with open(value_df, "r") as file:
            reader = csv.reader(file, delimiter="\t")
            header = next(reader)
            conditions = header[1:]  # Assumes the first column is gene ID

            # Initialize gene_values dictionary
            for row in reader:
                gene_id = row[0]
                gene_values[gene_id] = {}
                for i, condition in enumerate(conditions):
                    if len(row) > i + 1:
                        value = float(row[i + 1])
                    else:
                        value = 0.0  # Default to 0 if no value
                    gene_values[gene_id][condition] = value

        # Build the new dictionary structure by conditions
        for condition in conditions:
            new_dict[condition] = {}  # Create a new sub-dictionary for each condition

            # Deep copy the gene_dict and update with values from value_df
            for gene_id, gene_info in gene_dict.items():
                new_dict[condition][gene_id] = copy.deepcopy(gene_info)
                if gene_id in gene_values and condition in gene_values[gene_id]:
                    new_dict[condition][gene_id]["value"] = gene_values[gene_id][
                        condition
                    ]
                else:
                    new_dict[condition][gene_id][
                        "value"
                    ] = 0  # Default to 0 if the gene_id has no corresponding value

        return new_dict

    def update_transcript_values(self, gene_dict, value_df):
        new_dict = copy.deepcopy(gene_dict)  # Preserve the original structure
        transcript_values = {}

        # Load transcript counts from value_df
        with open(value_df, "r") as file:
            reader = csv.reader(file, delimiter="\t")
            header = next(reader)
            conditions = header[1:]  # Assumes the first column is transcript ID

            for row in reader:
                transcript_id = row[0]
                for i, condition in enumerate(conditions):
                    if len(row) > i + 1:
                        value = float(row[i + 1])
                    else:
                        value = 0.0  # Default to 0 if no value
                    if transcript_id not in transcript_values:
                        transcript_values[transcript_id] = {}
                    transcript_values[transcript_id][condition] = value

        # Update each condition without restructuring the original dictionary
        for condition in conditions:
            if condition not in new_dict:
                new_dict[condition] = copy.deepcopy(
                    gene_dict
                )  # Make sure all genes are present

            for gene_id, gene_info in new_dict[condition].items():
                if "transcripts" in gene_info:
                    for transcript_id, transcript_info in gene_info[
                        "transcripts"
                    ].items():
                        if (
                            transcript_id in transcript_values
                            and condition in transcript_values[transcript_id]
                        ):
                            transcript_info["value"] = transcript_values[transcript_id][
                                condition
                            ]
                        else:
                            transcript_info["value"] = (
                                0  # Set default if no value for this transcript
                            )
        return new_dict

    def update_gene_names(self, gene_dict):
        updated_dict = {}
        for condition, genes in gene_dict.items():
            updated_genes = {}
            for gene_id, gene_info in genes.items():
                if gene_info["name"]:
                    gene_name_upper = gene_info["name"].upper()
                    updated_genes[gene_name_upper] = gene_info
                else:
                    # If name is empty, use the original gene_id
                    updated_genes[gene_id] = gene_info
            updated_dict[condition] = updated_genes
        return updated_dict

    def filter_transcripts_by_minimum_value(self, gene_dict, min_value=1.0):
        # Dictionary to hold genes and transcripts that meet the criteria
        transcript_passes_threshold = {}

        # First pass: Determine which transcripts meet the minimum value requirement in any condition
        for condition, genes in gene_dict.items():
            for gene_id, gene_info in genes.items():
                for transcript_id, transcript_info in gene_info["transcripts"].items():
                    if (
                        "value" in transcript_info
                        and transcript_info["value"] != "NA"
                        and float(transcript_info["value"]) >= min_value
                    ):
                        if gene_id not in transcript_passes_threshold:
                            transcript_passes_threshold[gene_id] = {}
                        transcript_passes_threshold[gene_id][transcript_id] = True

        # Second pass: Build the filtered dictionary including only transcripts that have eligible values in any condition
        filtered_dict = {}
        for condition, genes in gene_dict.items():
            filtered_genes = {}
            for gene_id, gene_info in genes.items():
                if gene_id in transcript_passes_threshold:
                    eligible_transcripts = {
                        transcript_id: transcript_info
                        for transcript_id, transcript_info in gene_info[
                            "transcripts"
                        ].items()
                        if transcript_id in transcript_passes_threshold[gene_id]
                    }
                    if (
                        eligible_transcripts
                    ):  # Only add genes with non-empty transcript sets
                        filtered_gene_info = copy.deepcopy(gene_info)
                        filtered_gene_info["transcripts"] = eligible_transcripts
                        filtered_genes[gene_id] = filtered_gene_info
            if filtered_genes:  # Only add conditions with non-empty gene sets
                filtered_dict[condition] = filtered_genes

        return filtered_dict

    def read_gene_list(self, gene_list_path):
        with open(gene_list_path, "r") as file:
            gene_list = [
                line.strip().upper() for line in file
            ]  # Convert each gene to uppercase
        return gene_list

    def save_gene_dict_to_json(self, gene_dict, output_path):
        """Saves the gene dictionary to a JSON file."""
        # name the gene_dict file
        output_path = os.path.join(output_path, "gene_dict.json")
        with open(output_path, "w") as file:
            json.dump(gene_dict, file, indent=4)
