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
        self.file_list = file_list
        self.label = label
        self.out_dir = out_dir


class InputDataStorage:
    def __init__(self, args):
        self.samples = []
        self.input_type = ""
        sample_files = []
        labels = []

        if args.fastq is not None:
            self.input_type = "fastq"
            for fq in args.fastq:
                sample_files.append([[fq]])
        elif args.bam is not None:
            self.input_type = "bam"
            for bam in args.bam:
                sample_files.append([[bam]])
        elif args.fastq_list is not None:
            self.input_type = "fastq"
            sample_files = self.get_samples_from_file(args.fastq_list)
        elif args.bam_list is not None:
            self.input_type = "bam"
            sample_files = self.get_samples_from_file(args.bam_list)
        else:
            logger.critical("Input data was not specified")
            exit(-1)

        if args.labels is not None:
            if len(args.labels) != len(sample_files):
                logger.critical("Number of labels is not equal to the numbe of samples")
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
        return sample_files

    def get_labels(self, sample_files):
        labels = []
        for i in range(len(sample_files)):
            labels.append('{:02d}'.format(i) + "_" + os.path.splitext(os.path.basename(sample_files[i][0][0]))[0])
        return labels
