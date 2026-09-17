
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import csv
import os
import pickle
import gzip
import re
import shutil
import copy
import json
from argparse import Namespace
from collections import defaultdict
import tempfile
import gffutils
import yaml

from isoquant_lib.quantification.convert_grouped_counts import convert_to_matrix, GROUP_COUNT_CUTOFF
from isoquant_lib.scripts.convert_read_info import RI_CLASSIFICATION, RI_ISOFORM_ASSIGNMENT_TYPE
from isoquant_lib.utils.file_utils import open_text_read


class OutputConfig:
    """Class to build dictionaries from the output files of the pipeline."""

    # Per-read output, best first: the 4.0.0 unified format, then the deprecated one.
    # Matched as exact suffixes, so neither the same-format .transcript_model_reads.tsv
    # nor .read_assignments.SQANTI-like.tsv is picked up by mistake.
    READ_FILE_SUFFIXES = (".read_info.tsv", ".read_assignments.tsv")

    # Since 4.0.0 grouped counts carry the grouping strategy in the file name, e.g.
    # SAMPLE.gene_grouped_barcode_counts.linear.tsv. The leading dot and the
    # longest-first alternation keep ".discovered_transcript_grouped_..." from being
    # read as the plain "transcript_grouped" variant.
    GROUPED_COUNTS_RE = re.compile(
        r"\.(?P<feature>discovered_transcript|discovered_gene|transcript|gene)_grouped_"
        r"(?P<strategy>.+?)_(?P<value>counts|tpm)(?P<linear>\.linear)?\.tsv$"
    )

    # (feature, value) -> attribute holding that grouped file
    GROUPED_ATTRIBUTES = {
        ("gene", "counts"): "gene_grouped_counts",
        ("gene", "tpm"): "gene_grouped_tpm",
        ("transcript", "counts"): "transcript_grouped_counts",
        ("transcript", "tpm"): "transcript_grouped_tpm",
        ("discovered_transcript", "counts"): "transcript_model_grouped_counts",
        ("discovered_transcript", "tpm"): "transcript_model_grouped_tpm",
    }

    def __init__(self, output_directory, use_counts=False, ref_only=None, gtf=None,
                 read_group_strategy=None):
        self.output_directory = output_directory
        self.log_details = {}
        self.extended_annotation = None
        self.read_assignments = None
        self.input_gtf = gtf  # Initialize with the provided gtf flag
        self.genedb_filename = None
        self.yaml_input = True
        self.yaml_input_path = None
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
        # Grouping strategy (--read_group) whose counts are visualized, and every
        # strategy found in the output directory.
        self.read_group_strategy = read_group_strategy
        self.group_strategies = []
        # Lazily-created scratch dir for wide matrices converted from *.linear.tsv
        self._linear_conversion_dir = None

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

        # Paths in .params are stored as IsoQuant received them on the command line,
        # so a relative -o / --genedb produces relative entries here. Re-anchor them
        # against the parent of the output_directory the user passed to the visualizer
        # so visualize.py works regardless of the CWD.
        self.input_gtf = self._reanchor_relative_path(self.input_gtf)
        self.genedb_filename = self._reanchor_relative_path(self.genedb_filename)

        if params.get("yaml"):
            # YAML input case
            yaml_path = params.get("yaml")
            yaml_path = self._reanchor_relative_path(yaml_path)
            single_sample = self._yaml_single_sample_name(yaml_path)
            if single_sample is not None:
                # IsoQuant only writes combined_*.tsv when len(samples) > 1, so a
                # single-sample YAML run has the same on-disk layout as a non-YAML
                # run. Treat it that way to avoid the "combined file missing" path.
                self.yaml_input = False
                self.output_directory = os.path.join(
                    self.output_directory, single_sample
                )
            else:
                self.yaml_input = True
                self.yaml_input_path = yaml_path
            # Keep the output_directory as is for the multi-sample YAML case.
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

    def _yaml_single_sample_name(self, yaml_path):
        """Return the sole sample name if the YAML defines exactly one sample,
        otherwise None. Used to short-circuit the YAML branch when combined_*
        files are not produced by IsoQuant."""
        if not yaml_path or not os.path.exists(yaml_path):
            return None
        try:
            with open(yaml_path, "r") as yf:
                data = yaml.safe_load(yf)
        except Exception:
            return None
        samples = data if isinstance(data, list) else (data or {}).get("samples", [])
        names = [
            s.get("name")
            for s in samples
            if isinstance(s, dict) and s.get("name")
        ]
        return names[0] if len(names) == 1 else None

    def _reanchor_relative_path(self, path):
        """Try to resolve a relative path stored in .params against the user-provided
        output_directory (assumed to be the IsoQuant -o directory). Returns the
        original value when there is nothing to do."""
        if not path or os.path.isabs(path) or os.path.exists(path):
            return path
        anchor = os.path.dirname(os.path.abspath(self.output_directory))
        candidate = os.path.normpath(os.path.join(anchor, path))
        if os.path.exists(candidate):
            return candidate
        return path

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
        if self.yaml_input:
            self.conditions = True
            self.ref_only = True
            self._find_files_from_yaml()
            return  # Exit the method after processing YAML input

        if not os.path.exists(self.output_directory):
            print(f"Directory not found: {self.output_directory}")  # Debugging output
            raise FileNotFoundError(
                f"Specified sample subdirectory does not exist: {self.output_directory}"
            )

        # Grouped counts are written as *.linear.tsv; the wide matrix *.tsv is
        # only produced for --counts_format matrix/default (and skipped when
        # there are too many groups). Collect any linear variants so we can
        # convert them to the wide matrix the visualizer expects (#241).
        # Both are keyed by grouping strategy, one strategy per --read_group value.
        grouped = defaultdict(dict)
        linear_grouped_per_strategy = defaultdict(dict)
        # Strategies whose counts exist in MTX format only (the default above 100
        # groups, i.e. every single-cell / spatial run).
        mtx_only_strategies = set()
        # Per-read file candidates, keyed by the suffix that matched.
        read_files = {}

        for file_name in os.listdir(self.output_directory):
            if self._collect_read_file(file_name, read_files):
                continue
            if self._collect_grouped_file(file_name, grouped, linear_grouped_per_strategy,
                                          mtx_only_strategies):
                continue
            if file_name.endswith(".extended_annotation.gtf"):
                self.extended_annotation = os.path.join(
                    self.output_directory, file_name
                )
            elif file_name.endswith(".gene_counts.tsv"):
                self.gene_counts = os.path.join(self.output_directory, file_name)
            elif file_name.endswith(".transcript_counts.tsv"):
                self.transcript_counts = os.path.join(self.output_directory, file_name)
            elif file_name.endswith(".gene_tpm.tsv"):
                self.gene_tpm = os.path.join(self.output_directory, file_name)
            elif file_name.endswith(".transcript_tpm.tsv"):
                self.transcript_tpm = os.path.join(self.output_directory, file_name)
            elif file_name.endswith(".discovered_transcript_counts.tsv"):
                self.transcript_model_counts = os.path.join(
                    self.output_directory, file_name
                )
            elif file_name.endswith(".discovered_transcript_tpm.tsv"):
                self.transcript_model_tpm = os.path.join(
                    self.output_directory, file_name
                )

        self._select_read_file(read_files)
        self._select_grouped_files(grouped, linear_grouped_per_strategy, mtx_only_strategies)

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

    def _collect_read_file(self, file_name, read_files):
        """Record file_name if it is a per-read output (read_info or the deprecated
        read_assignments), gzipped or not. Returns True when it was consumed."""
        for suffix in self._read_file_suffixes():
            if file_name.endswith(suffix):
                read_files.setdefault(suffix, os.path.join(self.output_directory, file_name))
                return True
        return False

    def _select_read_file(self, read_files):
        """Keep the best per-read file found, preferring the current read_info format.
        Gzipped files are read as they are - the parser only counts two columns, and
        decompressing a single-cell read_info would write hundreds of GB next to it."""
        for suffix in self._read_file_suffixes():
            if suffix in read_files:
                self.read_assignments = read_files[suffix]
                return

    @classmethod
    def _read_file_suffixes(cls):
        """Per-read file suffixes in preference order: current format before the
        deprecated one, uncompressed before gzipped."""
        for suffix in cls.READ_FILE_SUFFIXES:
            yield suffix
            yield suffix + ".gz"

    def _collect_grouped_file(self, file_name, grouped, linear_grouped, mtx_only_strategies):
        """Record file_name if it is a grouped counts/TPM file, keyed by grouping
        strategy. Returns True when it was consumed."""
        if file_name.endswith(".matrix.mtx"):
            match = self.GROUPED_COUNTS_RE.search(file_name[: -len(".matrix.mtx")] + ".tsv")
            if match:
                mtx_only_strategies.add(match.group("strategy") or "")
            return True

        match = self.GROUPED_COUNTS_RE.search(file_name)
        if not match:
            return False
        attribute = self.GROUPED_ATTRIBUTES.get((match.group("feature"), match.group("value")))
        if attribute is None:
            return True  # e.g. discovered_gene counts, which the visualizer does not plot
        # Pre-4.0.0 files carry no strategy in the name; keep them under "".
        strategy = match.group("strategy") or ""
        target = linear_grouped if match.group("linear") else grouped
        target[strategy][attribute] = os.path.join(self.output_directory, file_name)
        return True

    def _select_grouped_files(self, grouped, linear_grouped, mtx_only_strategies):
        """Pick the grouping strategy to visualize and set the grouped attributes
        from it, converting linear files when no wide matrix was produced."""
        self.group_strategies = sorted(
            set(grouped) | set(linear_grouped) | mtx_only_strategies
        )
        if not self.group_strategies:
            return

        if self.read_group_strategy is not None:
            if self.read_group_strategy not in self.group_strategies:
                raise ValueError(
                    f"No grouped counts for read group strategy "
                    f"'{self.read_group_strategy}' in {self.output_directory}. "
                    f"Available: {', '.join(s or '(unnamed)' for s in self.group_strategies)}"
                )
            strategy = self.read_group_strategy
        else:
            strategy = self.group_strategies[0]
            if len(self.group_strategies) > 1:
                print(
                    f"Several read group strategies found "
                    f"({', '.join(s or '(unnamed)' for s in self.group_strategies)}); "
                    f"using '{strategy}'. Use --read_group_strategy to pick another one."
                )
        self.read_group_strategy = strategy

        for attribute, path in grouped.get(strategy, {}).items():
            setattr(self, attribute, path)
            self.conditions = True

        # Fall back to linear grouped files when the wide matrix is absent.
        self._resolve_linear_grouped_files(linear_grouped.get(strategy, {}))

        if not self.conditions:
            print(
                f"Grouped counts for strategy '{strategy}' are only available in MTX "
                f"format; per-group plots are skipped and ungrouped counts are used "
                f"instead. Re-run IsoQuant with '--counts_format matrix' to plot them."
                if strategy in mtx_only_strategies else
                f"No usable grouped counts found for strategy '{strategy}'; "
                f"using ungrouped counts."
            )

    def _resolve_linear_grouped_files(self, linear_grouped):
        """For each grouped attribute lacking a wide matrix, convert the matching
        *.linear.tsv file to the wide format the visualizer parses (#241)."""
        if not linear_grouped:
            return
        feature_types = {
            "gene_grouped_counts": "gene",
            "gene_grouped_tpm": "gene",
            "transcript_grouped_counts": "transcript",
            "transcript_grouped_tpm": "transcript",
            "transcript_model_grouped_counts": "transcript",
            "transcript_model_grouped_tpm": "transcript",
        }
        for attr, linear_path in linear_grouped.items():
            if getattr(self, attr) is not None:
                continue  # Wide matrix already present, prefer it.
            converted = self._convert_linear_to_matrix(
                linear_path, feature_types.get(attr, "gene")
            )
            if converted:
                setattr(self, attr, converted)
                # Read groups are present and usable, so grouped files are plotted.
                self.conditions = True

    def _convert_linear_to_matrix(self, linear_path, feature_type):
        """Convert a 3-column linear counts/TPM file into a wide matrix TSV in a
        scratch directory and return its path (or None on failure)."""
        if linear_path.endswith(".gz"):
            linear_path = self._unzip_file(linear_path)
            if not linear_path:
                return None
        if self._linear_conversion_dir is None:
            self._linear_conversion_dir = tempfile.mkdtemp(prefix="isoquant_viz_matrix_")
        base = os.path.basename(linear_path)
        if base.endswith(".linear.tsv"):
            base = base[: -len(".linear.tsv")]
        output_prefix = os.path.join(self._linear_conversion_dir, base)
        try:
            # convert_to_matrix appends ".tsv" to the given prefix. max_groups makes it
            # count the groups first and skip the pivot when there are too many: a
            # single-cell/spatial run has one group per barcode, and the dense matrix
            # would be unusable (and would not fit in memory) anyway.
            num_groups = convert_to_matrix(linear_path, output_prefix, feature_type=feature_type,
                                           max_groups=GROUP_COUNT_CUTOFF)
        except Exception as e:
            print(
                f"Warning: failed to convert linear counts {linear_path} "
                f"to matrix format: {e}"
            )
            return None
        output_path = output_prefix + ".tsv"
        if os.path.exists(output_path):
            return output_path
        if num_groups > GROUP_COUNT_CUTOFF:
            print(
                f"{os.path.basename(linear_path)} holds {num_groups} groups "
                f"(barcodes/spots); per-group plots are not meaningful at that size, "
                f"so ungrouped counts are used instead."
            )
        return None

    def _find_files_from_yaml(self):
        """Locate the necessary files in the directory, set specific grouped count and TPM files, and process read assignments."""
        if not os.path.exists(self.yaml_input_path):
            print(f"YAML file not found: {self.yaml_input_path}")
            raise FileNotFoundError(
                f"Specified YAML file does not exist: {self.yaml_input_path}"
            )

        # Set the four specific attributes
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

        # Check if yaml_data is a list
        if isinstance(yaml_data, list):
            samples = yaml_data
        else:
            # If it's not a list, assume it's a dictionary with a 'samples' key
            samples = yaml_data.get("samples", [])

        for sample in samples:
            name = sample.get("name")
            if name:
                sample_dir = os.path.join(self.output_directory, name)

                # read_info first, then the deprecated read_assignments; gzipped or not.
                read_file = None
                for suffix in self._read_file_suffixes():
                    path = os.path.join(sample_dir, name + suffix)
                    if os.path.exists(path):
                        read_file = path
                        break

                if read_file:
                    self.read_assignments.append((name, read_file))
                else:
                    print(f"Warning: No per-read file found for {name}")

        if not self.read_assignments:
            print("Warning: No per-read files found for any samples")


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
        """Indexes classifications and assignment types from the per-read file(s):
        SAMPLE.read_info.tsv[.gz] or the deprecated SAMPLE.read_assignments.tsv[.gz]."""
        if not self.config.read_assignments:
            raise FileNotFoundError(
                "No per-read file found: neither SAMPLE.read_info.tsv[.gz] nor "
                "SAMPLE.read_assignments.tsv[.gz] is present in "
                f"{self.config.output_directory}. Re-run IsoQuant with "
                "'--large_output read_info' to produce one."
            )

        if isinstance(self.config.read_assignments, list):
            # YAML input case (multiple files)
            classification_counts_dict = {}
            assignment_type_counts_dict = {}
            for sample_name, read_assignment_file in self.config.read_assignments:
                classification_counts, assignment_type_counts = (
                    self._process_read_assignment_file(read_assignment_file)
                )
                classification_counts_dict[sample_name] = classification_counts
                assignment_type_counts_dict[sample_name] = assignment_type_counts
            return classification_counts_dict, assignment_type_counts_dict
        else:
            # Non-YAML input case (single file)
            return self._process_read_assignment_file(self.config.read_assignments)

    # Legacy read_assignments layout: assignment type is the 6th column and the
    # classification is buried in the additional_info blob (the 9th).
    LEGACY_ASSIGNMENT_TYPE_COLUMN = 5
    LEGACY_ADDITIONAL_INFO_COLUMN = 8

    def _process_read_assignment_file(self, file_path):
        classification_counts = {}
        assignment_type_counts = {}

        # open_text_read decompresses .gz on the fly, so a huge single-cell read_info
        # is streamed rather than unpacked to disk first.
        with open_text_read(file_path) as file:
            type_column, classification_column = None, None
            for line in file:
                if line.startswith("#"):
                    continue  # command line / version banner
                parts = line.rstrip("\n").split("\t")

                if type_column is None:
                    # First non-comment line is the column header; both formats name
                    # their columns, so the layout is read from it.
                    type_column, classification_column, classification_in_blob = \
                        self._read_file_columns(parts)
                    continue

                if len(parts) <= max(type_column, classification_column):
                    continue
                assignment_type = parts[type_column]
                classification = parts[classification_column]
                if classification_in_blob:
                    # Legacy additional_info: "... Classification=full_splice_match; ..."
                    classification = (
                        classification.split("Classification=")[-1].split(";")[0].strip()
                    )

                classification_counts[classification] = (
                    classification_counts.get(classification, 0) + 1
                )
                assignment_type_counts[assignment_type] = (
                    assignment_type_counts.get(assignment_type, 0) + 1
                )

        return classification_counts, assignment_type_counts

    @staticmethod
    def _read_file_columns(header_parts):
        """Locate the assignment type and classification columns from the header of a
        per-read file. Returns (type_column, classification_column, classification_in_blob);
        in the legacy format the classification lives inside the additional_info field."""
        if "isoform_assignment_type" in header_parts and "classification" in header_parts:
            return (header_parts.index("isoform_assignment_type"),
                    header_parts.index("classification"), False)
        if "assignment_type" in header_parts and "additional_info" in header_parts:
            # The last column of a legacy file is groups, not additional_info, so the
            # blob has to be addressed by name.
            return (header_parts.index("assignment_type"),
                    header_parts.index("additional_info"), True)
        # Unnamed header (e.g. a manually trimmed file): fall back to the fixed
        # read_info layout, which is what IsoQuant produces by default.
        if len(header_parts) > RI_CLASSIFICATION:
            return RI_ISOFORM_ASSIGNMENT_TYPE, RI_CLASSIFICATION, False
        return (DictionaryBuilder.LEGACY_ASSIGNMENT_TYPE_COLUMN,
                DictionaryBuilder.LEGACY_ADDITIONAL_INFO_COLUMN, True)

    def parse_input_gtf(self):
        """Parses the GTF file using gffutils to build a detailed dictionary of genes, transcripts, and exons."""
        gene_dict = {}
        if not self.config.genedb_filename:
            # convert GTF to DB if we use previous IsoQuant runs
            # remove this functionality later
            tmp_file = tempfile.NamedTemporaryFile(suffix=".db")
            self.config.genedb_filename = tmp_file.name
            input_gtf_path = self.config.input_gtf
            gffutils.create_db(
                input_gtf_path,
                dbfn=self.config.genedb_filename,
                force=True,
                keep_order=True,
                merge_strategy="merge",
                sort_attribute_values=True,
                disable_infer_genes=True,
                disable_infer_transcripts=True,
            )

        try:
            # Create a database without using a context manager
            db = gffutils.FeatureDB(self.config.genedb_filename)

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
                        "name": transcript.attributes.get("transcript_name", [""])[0],
                        "biotype": transcript.attributes.get(
                            "transcript_biotype", [""]
                        )[0],
                        "exons": [],
                        "tags": transcript.attributes.get("tag", [""])[0].split(","),
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
        if value_df is None:
            kind = "counts" if self.config.use_counts else "TPM"
            raise FileNotFoundError(
                f"Could not locate a gene {kind} file in "
                f"{self.config.output_directory}. Re-run IsoQuant so the matching "
                f"*.tsv is produced, or drop --counts to use TPM outputs instead."
            )
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
        if value_df is None:
            kind = "counts" if self.config.use_counts else "TPM"
            quant = "reference" if self.config.ref_only else "discovered transcript model"
            raise FileNotFoundError(
                f"Could not locate a {quant} transcript {kind} file in "
                f"{self.config.output_directory}. Re-run IsoQuant so the matching "
                f"*.tsv is produced, or pass --ref_only / drop --counts to match the available outputs."
            )
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
