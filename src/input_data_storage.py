############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import logging

logger = logging.getLogger('IsoQuant')


class SampleData:
    def __init__(self, file_list, label, out_dir):
        # list of lists, since each sample may contain severa libraries, and each library may contain 2 files (paired)
        self.file_list = file_list
        self.label = label
        self.out_dir = out_dir
        self._init_paths()

    def _make_path(self, name):
        return os.path.join(self.out_dir, name)

    def _init_paths(self):
        self.out_assigned_tsv = self._make_path(self.label + ".read_assignments.tsv")
        self.out_raw_file = self._make_path(self.label + ".save")
        self.out_mapped_bed = self._make_path(self.label + ".mapped_reads.bed")
        self.out_alt_tsv = self._make_path(self.label + ".SQANTI-like.tsv")
        self.out_gene_counts_tsv = self._make_path(self.label + ".gene_counts.tsv")
        self.out_transcript_counts_tsv = self._make_path(self.label + ".transcript_counts.tsv")
        self.out_exon_counts_tsv = self._make_path(self.label + ".exon_counts.tsv")
        self.out_intron_counts_tsv = self._make_path(self.label + ".intron_counts.tsv")
        self.out_gene_grouped_counts_tsv = self._make_path(self.label + ".gene_grouped_counts.tsv")
        self.out_transcript_grouped_counts_tsv = self._make_path(self.label + ".transcript_grouped_counts.tsv")
        self.out_exon_grouped_counts_tsv = self._make_path(self.label + ".exon_grouped_counts.tsv")
        self.out_intron_grouped_counts_tsv = self._make_path(self.label + ".intron_grouped_counts.tsv")


class InputDataStorage:
    def __init__(self, args):
        # list of SampleData
        self.samples = []
        self.input_type = ""
        sample_files = []
        labels = []

        if args.fastq is not None:
            self.input_type = "fastq"
            for fq in args.fastq:
                check_input_type(fq, self.input_type)
                sample_files.append([[fq]])
        elif args.bam is not None:
            self.input_type = "bam"
            for bam in args.bam:
                check_input_type(bam, self.input_type)
                sample_files.append([[bam]])
        elif args.fastq_list is not None:
            self.input_type = "fastq"
            sample_files = self.get_samples_from_file(args.fastq_list)
        elif args.bam_list is not None:
            self.input_type = "bam"
            sample_files = self.get_samples_from_file(args.bam_list)
        elif args.read_assignments is not None:
            self.input_type = "save"
            for save_file in args.read_assignments:
                # TODO: check
                sample_files.append([[save_file]])
        else:
            logger.critical("Input data was not specified")
            exit(-1)

        if args.labels is not None:
            if len(args.labels) != len(sample_files):
                logger.critical("Number of labels is not equal to the number of samples")
                exit(-1)
            else:
                labels = args.labels
        else:
            labels = self.get_labels(sample_files)

        for i in range(len(sample_files)):
            self.samples.append(SampleData(sample_files[i], labels[i], os.path.join(args.output, labels[i])))

    def get_samples_from_file(self, file_name):
        sample_files = []
        inf = open(file_name, "r")
        current_sample = []

        for l in inf:
            if len(l.strip()) == 0 and len(current_sample) > 0:
                sample_files.append(current_sample)
                current_sample = []
            else:
                current_sample.append(l.strip().split())
        if len(current_sample) > 0:
            sample_files.append(current_sample)

        for sample in sample_files:
            for lib in sample:
                for in_file in lib:
                    check_input_type(in_file, self.input_type)
        return sample_files

    def get_labels(self, sample_files):
        labels = []
        for i in range(len(sample_files)):
            labels.append('{:02d}'.format(i) + "_" + os.path.splitext(os.path.basename(sample_files[i][0][0]))[0])
        return labels


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

