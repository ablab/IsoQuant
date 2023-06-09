############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import logging
from collections import defaultdict


logger = logging.getLogger('IsoQuant')


class SampleData:
    def __init__(self, file_list, prefix, out_dir, readable_names_dict):
        # list of lists, since each sample may contain several libraries, and each library may contain 2 files (paired)
        self.file_list = file_list
        self.readable_names_dict = readable_names_dict
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
        self.out_transcript_model_counts_tsv = self._make_path(self.prefix + ".transcript_model")
        self.out_transcript_model_grouped_counts_tsv = self._make_path(self.prefix + ".transcript_model_grouped")
        self.out_exon_counts_tsv = self._make_path(self.prefix + ".exon")
        self.out_intron_counts_tsv = self._make_path(self.prefix + ".intron")
        self.out_gene_grouped_counts_tsv = self._make_path(self.prefix + ".gene_grouped")
        self.out_transcript_grouped_counts_tsv = self._make_path(self.prefix + ".transcript_grouped")
        self.out_exon_grouped_counts_tsv = self._make_path(self.prefix + ".exon_grouped")
        self.out_intron_grouped_counts_tsv = self._make_path(self.prefix + ".intron_grouped")
        self.out_t2t_tsv = self._make_path(self.prefix + ".novel_vs_known.SQANTI-like.tsv")


class InputDataStorage:
    def __init__(self, args):
        # list of SampleData
        self.samples = []
        self.input_type = ""
        readable_names_dict = defaultdict(lambda: defaultdict(str))
        sample_files = []
        experiment_names = []
        self.experiment_prefix = args.prefix

        if args.fastq is not None:
            self.input_type = "fastq"
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

        elif args.bam is not None:
            self.input_type = "bam"
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

        elif args.fastq_list is not None:
            self.input_type = "fastq"
            sample_files, experiment_names, readable_names_dict = self.get_samples_from_file(args.fastq_list)
            if args.labels:
                logger.warning("--labels option has no effect when files are provided via input list")

        elif args.bam_list is not None:
            self.input_type = "bam"
            sample_files, experiment_names, readable_names_dict = self.get_samples_from_file(args.bam_list)
            if args.labels:
                logger.warning("--labels option has no effect when files are provided via input list")

        elif args.read_assignments is not None:
            self.input_type = "save"
            for i, save_file in enumerate(args.read_assignments):
                sample_files.append([[save_file]])
                experiment_names.append(self.experiment_prefix + str(i))

        else:
            logger.critical("Input data was not specified")
            exit(-1)

        for i in range(len(sample_files)):
            self.samples.append(SampleData(sample_files[i], experiment_names[i],
                                           os.path.join(args.output, experiment_names[i]),
                                           readable_names_dict[experiment_names[i]]))

    def get_samples_from_file(self, file_name):
        sample_files = []
        experiment_names = []
        readable_names_dict = defaultdict(lambda: defaultdict(str))
        inf = open(file_name, "r")
        current_sample = []
        current_sample_name = ""
        current_index = 0

        for l in inf:
            if len(l.strip()) == 0 or l.startswith("#"):
                if len(current_sample) > 0:
                    sample_files.append(current_sample)
                    experiment_names.append(current_sample_name)
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

        for sample in sample_files:
            for lib in sample:
                for in_file in lib:
                    check_input_type(in_file, self.input_type)
        return sample_files, experiment_names, readable_names_dict

    def has_replicas(self):
        return any(len(sample.file_list) > 1 for sample in self.samples)


def check_input_type(fname, input_type):
    basename_plus_inner_ext, outer_ext = os.path.splitext(fname.lower())
    if outer_ext not in ['.zip', '.gz', '.gzip', '.bz2', '.bzip2']:
        basename_plus_inner_ext, outer_ext = fname, ''  # not a supported archive

    basename, fasta_ext = os.path.splitext(basename_plus_inner_ext)
    if fasta_ext in ['.fastq', '.fasta', '.fa', '.fq', '.fna']:
        if input_type != 'fastq':
            raise Exception("Wrong file extension was detected. Use only FASTQ/FASTA files with --fastq option.")
    elif fasta_ext == '.bam':
        if input_type != 'bam':
            raise Exception("Wrong file extension was detected. Use only BAM files with --bam option.")
    else:
        raise Exception("File format " + fasta_ext + " is not supported! Supported formats: FASTQ, FASTA, BAM")

