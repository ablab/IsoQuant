############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import logging
import yaml
from collections import defaultdict
from enum import Enum, unique

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
    def __init__(self, file_list, prefix, out_dir, readable_names_dict, illumina_bam):
        # list of lists, since each sample may contain several libraries, and each library may contain 2 files (paired)
        self.file_list = file_list
        self.readable_names_dict = readable_names_dict
        self.illumina_bam = illumina_bam
        self.prefix = prefix
        self.out_dir = out_dir
        self.aux_dir = os.path.join(self.out_dir, "aux")
        self._init_paths()

    def _make_path(self, name):
        return os.path.join(self.out_dir, name)

    def _make_aux_path(self, name):
        return os.path.join(self.aux_dir, name)

    def _init_paths(self):
        self.out_assigned_tsv = self._make_path(self.prefix + ".read_assignments.tsv")
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
        self.out_polya = self._make_path(self.prefix + ".polyA_prediction.tsv")


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
                exit(-1)
            for i, fq in enumerate(args.fastq):
                check_input_type(fq, self.input_type)
                if fq in readable_names_dict[experiment_name]:
                    logger.critical("File %s is used multiple times in a single experiment, which is not allowed" % fq)
                    exit(-2)
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
                exit(-1)
            for i, bam in enumerate(args.bam):
                check_input_type(bam, self.input_type)
                if bam in readable_names_dict[experiment_name]:
                    logger.critical("File %s is used multiple times in a single experiment, which is not allowed" % bam)
                    exit(-2)
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
                exit(-1)
            for i, bam in enumerate(args.unmapped_bam):
                check_input_type(bam, self.input_type)
                if bam in readable_names_dict[experiment_name]:
                    logger.critical("File %s is used multiple times in a single experiment, which is not allowed" % bam)
                    exit(-2)
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
            exit(-1)

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
                        exit(-1)
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
                        logger.critical("File %s is used multiple times in a single experiment, which is not allowed" % fname)
                        exit(-2)
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
            exit(-2)
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
                exit(-1)
        for sample in con[1:]:
            if not 'name' in sample.keys():
                current_sample_name = self.experiment_prefix + str(current_index)
            else:
                current_sample_name = sample['name']
            if current_sample_name in experiment_names:
                    new_sample_name = self.experiment_prefix + str(current_index)
                    if current_sample_name == new_sample_name:
                        logger.critical("Change experiment name %s and rerun IsoQuant" % current_sample_name)
                        exit(-1)
                    logger.warning("Duplicate folder prefix %s, will change to %s" %
                                   (current_sample_name, new_sample_name))
                    current_sample_name = new_sample_name
            current_index += 1
            if not 'long read files' in sample.keys():
                logger.critical("Experiment %s does not contain any files" %current_sample_name)
                exit(-2)
            else:
                current_sample = [normalize_path(yaml_file_path, b) for b in sample['long read files']]
                names = 'labels' in sample.keys()
                if names and not len(sample['labels']) == len(current_sample):
                    logger.critical("The number of file aliases differs from the number of files")
                    exit(-2)
                for f in range(len(current_sample)):
                    fname = current_sample[f]
                    if names:
                        readable_name = sample['labels'][f]
                    else:
                        readable_name = os.path.splitext(os.path.basename(fname))[0]
                    if fname in readable_names_dict[current_sample_name]:
                        logger.critical("File %s is used multiple times in a single experiment, which is not allowed" % fname)
                        exit(-2)
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

