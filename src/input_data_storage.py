############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import logging
import yaml
from collections import defaultdict
from enum import Enum, unique

from .error_codes import IsoQuantExitCode
from .file_utils import normalize_path


logger = logging.getLogger('IsoQuant')


@unique
class InputDataType(Enum):
    undefined = 0
    fastq = 1
    bam = 2
    unmapped_bam = 3
    save = 10

    def needs_mapping(self):
        return self in [InputDataType.fastq, InputDataType.unmapped_bam]


class SampleData:
    def __init__(self, file_list, prefix, out_dir, readable_names_dict, illumina_bam, barcoded_reads=None):
        # list of lists, since each sample may contain several libraries, and each library may contain 2 files (paired)
        self.file_list = file_list
        self.readable_names_dict = readable_names_dict
        self.illumina_bam = illumina_bam
        self.prefix = prefix
        self.out_dir = out_dir
        self.aux_dir = os.path.join(self.out_dir, "aux")
        if not barcoded_reads:
            self.barcoded_reads = []
        else:
            self.barcoded_reads = barcoded_reads
        self.use_technical_replicas = False  # Will be set by DatasetProcessor
        self._init_paths()

    def _make_path(self, name):
        return os.path.join(self.out_dir, name)

    def _make_aux_path(self, name):
        return os.path.join(self.aux_dir, name)

    def _init_paths(self):
        self.out_assigned_tsv = self._make_path(self.prefix + ".read_assignments.tsv")
        self.out_assigned_tsv_result = self.out_assigned_tsv
        self.out_raw_file = self._make_aux_path(self.prefix + ".save")
        self.read_group_file = self._make_aux_path(self.prefix + ".read_group")
        self.out_corrected_bed = self._make_path(self.prefix + ".corrected_reads.bed")
        self.out_alt_tsv = self._make_path(self.prefix + ".read_assignments.SQANTI-like.tsv")
        self.out_gene_counts_tsv = self._make_path(self.prefix + ".gene")
        self.out_transcript_counts_tsv = self._make_path(self.prefix + ".transcript")
        self.out_transcript_model_counts_tsv = self._make_path(self.prefix + ".discovered_transcript")
        self.out_gene_model_counts_tsv = self._make_path(self.prefix + ".discovered_gene")
        self.out_transcript_model_grouped_counts_tsv = self._make_path(self.prefix + ".discovered_transcript_grouped")
        self.out_gene_model_grouped_counts_tsv = self._make_path(self.prefix + ".discovered_gene_grouped")
        self.out_exon_counts_tsv = self._make_path(self.prefix + ".exon")
        self.out_intron_counts_tsv = self._make_path(self.prefix + ".intron")
        self.out_gene_grouped_counts_tsv = self._make_path(self.prefix + ".gene_grouped")
        self.out_transcript_grouped_counts_tsv = self._make_path(self.prefix + ".transcript_grouped")
        self.out_exon_grouped_counts_tsv = self._make_path(self.prefix + ".exon_grouped")
        self.out_intron_grouped_counts_tsv = self._make_path(self.prefix + ".intron_grouped")
        self.out_t2t_tsv = self._make_path(self.prefix + ".novel_vs_known.SQANTI-like.tsv")
        self.barcodes_tsv = self._make_path(self.prefix + ".barcoded_reads")
        self.barcodes_done = self._make_aux_path(self.prefix + ".barcodes_done")
        self.barcodes_split_reads = self._make_aux_path(self.prefix + ".split_barcodes")
        self.out_umi_filtered = self._make_path(self.prefix + ".UMI_filtered")
        self.out_umi_filtered_tmp = self._make_aux_path(self.prefix + ".UMI_filtered")
        self.out_umi_filtered_done= self._make_aux_path(self.prefix + ".UMI_filtered.done")
        self.split_reads_fasta = self._make_path(self.prefix + ".split_reads")

    # Chromosome-specific path methods (delegate to file_naming.py)

    def get_save_file(self, chr_id: str) -> str:
        """Get path to serialized read assignments for a chromosome."""
        from .file_naming import saves_file_name
        return saves_file_name(self.out_raw_file, chr_id)

    def get_multimappers_file(self, chr_id: str) -> str:
        """Get path to multimapper assignments for a chromosome."""
        from .file_naming import multimappers_file_name
        return multimappers_file_name(self.out_raw_file, chr_id)

    def get_filtered_reads_file(self, chr_id: str) -> str:
        """Get path to filtered reads for a chromosome."""
        from .file_naming import filtered_reads_file_name
        return filtered_reads_file_name(self.out_raw_file, chr_id)

    def get_dynamic_pools_file(self, chr_id: str) -> str:
        """Get path to dynamic string pools for a chromosome."""
        from .file_naming import dynamic_pools_file_name
        return dynamic_pools_file_name(self.out_raw_file, chr_id)

    def get_collected_lock_file(self, chr_id: str) -> str:
        """Get path to lock file indicating reads collected for a chromosome."""
        from .file_naming import reads_collected_lock_file_name
        return reads_collected_lock_file_name(self.out_raw_file, chr_id)

    def get_processed_lock_file(self, chr_id: str) -> str:
        """Get path to lock file indicating reads processed for a chromosome."""
        from .file_naming import reads_processed_lock_file_name
        return reads_processed_lock_file_name(self.out_raw_file, chr_id)

    def get_read_group_split_file(self, chr_id: str, spec_index: int = None) -> str:
        """Get path to split read group file for a chromosome.

        Args:
            chr_id: Chromosome ID
            spec_index: Optional spec index for multi-spec files

        Returns:
            Path like 'prefix.read_group_spec0_chr1' or 'prefix.read_group_chr1'
        """
        from .file_naming import convert_chr_id_to_file_name_str
        chr_str = convert_chr_id_to_file_name_str(chr_id)
        if spec_index is not None:
            return f"{self.read_group_file}_spec{spec_index}_{chr_str}"
        return f"{self.read_group_file}_{chr_str}"

    def get_barcodes_split_file(self, chr_id: str) -> str:
        """Get path to split barcodes file for a chromosome."""
        from .file_naming import convert_chr_id_to_file_name_str
        return self.barcodes_split_reads + "_" + convert_chr_id_to_file_name_str(chr_id)

    # Chromosome-specific output file getters (for parallel processing)

    def get_chr_prefix(self, chr_id: str) -> str:
        """Get chromosome-specific prefix for output files."""
        from .file_naming import convert_chr_id_to_file_name_str
        return f"{self.prefix}_{convert_chr_id_to_file_name_str(chr_id)}"

    def get_corrected_bed_file(self, chr_id: str) -> str:
        """Get path to corrected reads BED file for a chromosome."""
        return self._make_path(self.get_chr_prefix(chr_id) + ".corrected_reads.bed")

    def get_assigned_tsv_file(self, chr_id: str) -> str:
        """Get path to read assignments TSV for a chromosome."""
        return self._make_path(self.get_chr_prefix(chr_id) + ".read_assignments.tsv")

    def get_gene_counts_file(self, chr_id: str) -> str:
        """Get path to gene counts file for a chromosome."""
        return self._make_path(self.get_chr_prefix(chr_id) + ".gene")

    def get_transcript_counts_file(self, chr_id: str) -> str:
        """Get path to transcript counts file for a chromosome."""
        return self._make_path(self.get_chr_prefix(chr_id) + ".transcript")

    def get_exon_counts_file(self, chr_id: str) -> str:
        """Get path to exon counts file for a chromosome."""
        return self._make_path(self.get_chr_prefix(chr_id) + ".exon")

    def get_intron_counts_file(self, chr_id: str) -> str:
        """Get path to intron counts file for a chromosome."""
        return self._make_path(self.get_chr_prefix(chr_id) + ".intron")

    def get_transcript_model_counts_file(self, chr_id: str) -> str:
        """Get path to discovered transcript counts file for a chromosome."""
        return self._make_path(self.get_chr_prefix(chr_id) + ".discovered_transcript")

    def get_gene_model_counts_file(self, chr_id: str) -> str:
        """Get path to discovered gene counts file for a chromosome."""
        return self._make_path(self.get_chr_prefix(chr_id) + ".discovered_gene")

    def get_t2t_tsv_file(self, chr_id: str) -> str:
        """Get path to SQANTI-like TSV for a chromosome."""
        return self._make_path(self.get_chr_prefix(chr_id) + ".novel_vs_known.SQANTI-like.tsv")

    def get_grouped_counts_file(self, chr_id: str, feature: str, strategy_name: str) -> str:
        """Get path to grouped counts file for a chromosome.

        Args:
            chr_id: Chromosome ID
            feature: Feature type (gene, transcript, exon, intron, discovered_transcript, discovered_gene)
            strategy_name: Grouping strategy name
        """
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.{feature}_grouped_{strategy_name}")

    # Sample-level auxiliary files (for resume functionality)

    def get_info_file(self) -> str:
        """Get path to sample collection info file."""
        from .file_naming import info_file_name
        return info_file_name(self.out_raw_file)

    def get_collection_lock_file(self) -> str:
        """Get path to sample collection lock file."""
        from .file_naming import collection_lock_file_name
        return collection_lock_file_name(self.out_raw_file)


class InputDataStorage:
    def __init__(self, args):
        # list of SampleData
        self.samples = []
        self.input_type = InputDataType.undefined
        readable_names_dict = defaultdict(lambda: defaultdict(str))
        sample_files = []
        experiment_names = []
        self.experiment_prefix = args.prefix
        illumina_bam = []

        if args.fastq is not None:
            self.input_type = InputDataType.fastq
            sample_files.append([])
            experiment_name = args.prefix
            experiment_names.append(experiment_name)
            if args.labels and len(args.labels) != len(args.fastq):
                logger.critical("Number of labels is not equal to the number of files")
                sys.exit(IsoQuantExitCode.LABEL_COUNT_MISMATCH)
            for i, fq in enumerate(args.fastq):
                check_input_type(fq, self.input_type)
                if fq in readable_names_dict[experiment_name]:
                    logger.critical("File %s is used multiple times in a single experiment, which is not allowed" %fq)
                    sys.exit(IsoQuantExitCode.DUPLICATE_FILES)
                sample_files[0].append([fq])
                readable_names_dict[experiment_name][fq] = args.labels[i] if args.labels else \
                    os.path.splitext(os.path.basename(fq))[0]
            # if args.illumina_bam is not None:
                # illumina_bam.append(args.illumina_bam)
            # else:
                # illumina_bam.append([None])
            illumina_bam.append(args.illumina_bam)

        elif args.bam is not None:
            self.input_type = InputDataType.bam
            sample_files.append([])
            experiment_name = args.prefix
            experiment_names.append(experiment_name)
            if args.labels and len(args.labels) != len(args.bam):
                logger.critical("Number of labels is not equal to the number of files")
                sys.exit(IsoQuantExitCode.LABEL_COUNT_MISMATCH)
            for i, bam in enumerate(args.bam):
                check_input_type(bam, self.input_type)
                if bam in readable_names_dict[experiment_name]:
                    logger.critical("File %s is used multiple times in a single experiment, which is not allowed" %bam)
                    sys.exit(IsoQuantExitCode.DUPLICATE_FILES)
                sample_files[0].append([bam])
                readable_names_dict[experiment_name][bam] = args.labels[i] if args.labels else \
                    os.path.splitext(os.path.basename(bam))[0]
            illumina_bam.append(args.illumina_bam)

        elif args.unmapped_bam is not None:
            self.input_type = InputDataType.unmapped_bam
            sample_files.append([])
            experiment_name = args.prefix
            experiment_names.append(experiment_name)
            if args.labels and len(args.labels) != len(args.unmapped_bam):
                logger.critical("Number of labels is not equal to the number of files")
                sys.exit(IsoQuantExitCode.LABEL_COUNT_MISMATCH)
            for i, bam in enumerate(args.unmapped_bam):
                check_input_type(bam, self.input_type)
                if bam in readable_names_dict[experiment_name]:
                    logger.critical("File %s is used multiple times in a single experiment, which is not allowed" %bam)
                    sys.exit(IsoQuantExitCode.DUPLICATE_FILES)
                sample_files[0].append([bam])
                readable_names_dict[experiment_name][bam] = args.labels[i] if args.labels else \
                    os.path.splitext(os.path.basename(bam))[0]
            illumina_bam.append(args.illumina_bam)

        elif args.read_assignments is not None:
            self.input_type = InputDataType.save
            illumina_bam = [[]]
            for i, save_file in enumerate(args.read_assignments):
                sample_files.append([[save_file]])
                experiment_names.append(self.experiment_prefix + str(i))
        
        elif args.yaml is not None:
            sample_files, experiment_names, readable_names_dict, illumina_bam = self.get_samples_from_yaml(args.yaml)
            if args.labels:
                logger.warning("--labels option has no effect when files are provided via yaml file")

        else:
            logger.critical("Input data was not specified")
            sys.exit(IsoQuantExitCode.NO_INPUT_DATA)

        for i in range(len(sample_files)):
            self.samples.append(SampleData(sample_files[i], experiment_names[i],
                                           os.path.join(args.output, experiment_names[i]),
                                           readable_names_dict[experiment_names[i]],
                                           illumina_bam[i]))

    def get_samples_from_file(self, file_name):
        sample_files = []
        experiment_names = []
        illumina_bam = []
        readable_names_dict = defaultdict(lambda: defaultdict(str))
        inf = open(file_name, "r")
        current_sample = []
        current_sample_name = self.experiment_prefix
        current_index = 0

        for l in inf:
            if len(l.strip()) == 0 or l.startswith("#"):
                if len(current_sample) > 0:
                    sample_files.append(current_sample)
                    experiment_names.append(current_sample_name)
                    illumina_bam.append(None)
                current_sample = []
                current_sample_name = l.strip()[1:]
                if not current_sample_name:
                    current_sample_name = self.experiment_prefix + str(current_index)
                if current_sample_name in experiment_names:
                    new_sample_name = self.experiment_prefix + str(current_index)
                    if current_sample_name == new_sample_name:
                        logger.critical("Change experiment name %s and rerun IsoQuant" % current_sample_name)
                        sys.exit(IsoQuantExitCode.DUPLICATE_FILES)
                    logger.warning("Duplicate folder prefix %s, will change to %s" %
                                   (current_sample_name, new_sample_name))
                    current_sample_name = new_sample_name
                current_index += 1
            else:
                vals = l.strip().split(':')
                files = vals[0].split()
                if len(vals) > 1:
                    readable_name = vals[-1]
                else:
                    readable_name = os.path.splitext(os.path.basename(files[0]))[0]
                current_sample.append(files)
                for fname in files:
                    if fname in readable_names_dict[current_sample_name]:
                        logger.critical("File %s is used multiple times in a single experiment, which is not allowed" %fname)
                        sys.exit(IsoQuantExitCode.DUPLICATE_FILES)
                    readable_names_dict[current_sample_name][fname] = readable_name

        if len(current_sample) > 0:
            sample_files.append(current_sample)
            experiment_names.append(current_sample_name)
            illumina_bam.append(None)

        for sample in sample_files:
            for lib in sample:
                for in_file in lib:
                    check_input_type(in_file, self.input_type)

        return sample_files, experiment_names, readable_names_dict, illumina_bam

    def has_replicas(self):
        return any(len(sample.file_list) > 1 for sample in self.samples)
        
    def get_samples_from_yaml(self, yaml_file_path):
        sample_files = []
        experiment_names = []
        illumina_bam = []
        readable_names_dict = defaultdict(lambda: defaultdict(str))
        yaml_file = open(yaml_file_path, 'r')
        con = yaml.safe_load(yaml_file)
        current_index = 0
        t = con[0]
        if not 'data format' in t.keys():
            logger.critical("Please specify whether you are using fastq or bam files in the first entry")
            sys.exit(IsoQuantExitCode.YAML_PARSING_ERROR)
        else:
            if len(t.keys()) > 1:
                logger.warning("The first entry should only specify the input data format. Any additional info will be ignored")
            if  t['data format'] == "bam":
                self.input_type = InputDataType.bam
            elif t['data format'] == "unmapped_bam":
                self.input_type = InputDataType.unmapped_bam
            elif t['data format'] == "fastq" or t['data format'] == "fasta":
                self.input_type = InputDataType.fastq
            else:
                logger.critical("The input data format can only be either fastq, fasta or bam.")
                sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
        for sample in con[1:]:
            if not 'name' in sample.keys():
                current_sample_name = self.experiment_prefix + str(current_index)
            else:
                current_sample_name = sample['name']
            if current_sample_name in experiment_names:
                    new_sample_name = self.experiment_prefix + str(current_index)
                    if current_sample_name == new_sample_name:
                        logger.critical("Change experiment name %s and rerun IsoQuant" % current_sample_name)
                        sys.exit(IsoQuantExitCode.DUPLICATE_FILES)
                    logger.warning("Duplicate folder prefix %s, will change to %s" %
                                   (current_sample_name, new_sample_name))
                    current_sample_name = new_sample_name
            current_index += 1
            if not 'long read files' in sample.keys():
                logger.critical("Experiment %s does not contain any files" %current_sample_name)
                sys.exit(IsoQuantExitCode.YAML_PARSING_ERROR)
            else:
                current_sample = [normalize_path(yaml_file_path, b) for b in sample['long read files']]
                names = 'labels' in sample.keys()
                if names and not len(sample['labels']) == len(current_sample):
                    logger.critical("The number of file aliases differs from the number of files")
                    sys.exit(IsoQuantExitCode.LABEL_COUNT_MISMATCH)
                for f in range(len(current_sample)):
                    fname = current_sample[f]
                    if names:
                        readable_name = sample['labels'][f]
                    else:
                        readable_name = os.path.splitext(os.path.basename(fname))[0]
                    if fname in readable_names_dict[current_sample_name]:
                        logger.critical("File %s is used multiple times in a single experiment, which is not allowed" %fname)
                        sys.exit(IsoQuantExitCode.DUPLICATE_FILES)
                    readable_names_dict[current_sample_name][fname] = readable_name
            if len(current_sample) > 0:
                current_sample_list = [[s] for s in current_sample]
                sample_files.append(current_sample_list)
                experiment_names.append(current_sample_name)
                if 'illumina bam' in sample.keys():
                    illumina_bam.append([normalize_path(yaml_file_path, ib) for ib in sample['illumina bam']])
                else:
                    illumina_bam.append(None)
            
        for sample in sample_files:
            for lib in sample:
                for in_file in lib:
                    check_input_type(in_file, self.input_type)
        return sample_files, experiment_names, readable_names_dict, illumina_bam
        
# not functional yet
# idea for the future to name unnamed samples by their last common folder
    # def get_sample_name(names, index):
        # common_characters = len(names[0])
        # common_name = names[0]
        
        # for i in range(1, len(names)):
            # p = mismatch(common_name, names[i])
            # if p[0] < common_characters:
                # common_characters = p[0]
                
        # found = common_names.rfind('/', 0, common_characters)
        
        # common_name = common_name[:found]
        
        # sample_name = common_name + str(index)
        # return sample_name


def check_input_type(fname, input_type):
    basename_plus_inner_ext, outer_ext = os.path.splitext(fname.lower())
    if outer_ext not in ['.zip', '.gz', '.gzip', '.bz2', '.bzip2']:
        basename_plus_inner_ext, outer_ext = fname, ''  # not a supported archive

    basename, fasta_ext = os.path.splitext(basename_plus_inner_ext)
    if fasta_ext in ['.fastq', '.fasta', '.fa', '.fq', '.fna']:
        if input_type != InputDataType.fastq:
            raise Exception("Wrong file extension was detected %s. Use only FASTQ/FASTA files with --fastq option." % fname)
    elif fasta_ext == '.bam':
        if input_type not in [InputDataType.bam, InputDataType.unmapped_bam]:
            raise Exception("Wrong file extension was detected for file %s. Use only BAM files with --bam option." % fname)
    else:
        raise Exception("File format " + fasta_ext + " is not supported! Supported formats: FASTQ, FASTA, BAM")

