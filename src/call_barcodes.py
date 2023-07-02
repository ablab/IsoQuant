# ############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


from src.sc.detect_barcodes import call_barcodes_in_parallel, call_barcodes_single_thread, BarcodeCallingAlgorithm


class IsoQuantBarcodeCaller:
    def __init__(self, barcodes_file):
        self.barcodes_file = barcodes_file



def call_barcodes(args):
    samples = []
    annotation_file = find_annotation(self.aligner, args)
    for sample in args.input_data.samples:
        readable_names_dict = {}
        bam_files = []
        for fastq_files in sample.file_list:
            for fastq_file in fastq_files:
                bam_file = None if args.clean_start else find_stored_alignment(fastq_file, annotation_file, args)
                if bam_file is None:
                    bam_file = align_fasta(self.aligner, fastq_file, annotation_file, args, sample.prefix, sample.aux_dir)
                    store_alignment(bam_file, fastq_file, annotation_file, args)
                bam_files.append([bam_file])
                if fastq_file in sample.readable_names_dict:
                    readable_names_dict[bam_file] = sample.readable_names_dict[fastq_file]

        samples.append(SampleData(bam_files, sample.prefix, sample.out_dir, readable_names_dict))
    args.input_data.samples = samples
    args.input_data.input_type = "bam"
    return args.input_data