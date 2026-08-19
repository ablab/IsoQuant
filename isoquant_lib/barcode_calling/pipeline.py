############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Running the barcode calling stage of the IsoQuant pipeline.

Two passes when cell barcodes are detected from the data: the first extracts and counts
barcode windows, the second matches reads against the barcodes that pass selection.
"""

import concurrent.futures
import logging
import os
from concurrent.futures import ProcessPoolExecutor

from ..common import setup_worker_logging, _get_log_params
from .detect_barcodes import detect_cell_barcode_list, get_barcode_length, \
    process_in_parallel, process_single_thread

logger = logging.getLogger('IsoQuant')


class BarcodeCallingArgs:
    def __init__(self, input, barcode_whitelist, mode, output, out_fasta, tmp_dir, threads,
                 molecule: str = None, whitelist_matching: bool = True, split_molecules: bool = False):
        self.input = input  # Can be a single file (str) or list of files
        self.barcodes = barcode_whitelist
        self.mode = mode
        self.output_tsv = output  # Can be a single filename (str) or list of filenames
        self.out_fasta = out_fasta  # Can be a single filename (str), list of filenames, or None
        self.tmp_dir = tmp_dir
        self.threads = threads
        self.molecule = molecule
        # when False the detector emits raw barcode windows so cell barcodes can be detected
        self.whitelist_matching = whitelist_matching
        # selects the splitting detector, and makes the caller expect one result per molecule
        self.split_molecules = split_molecules


def run_barcode_calling(bc_args, threads):
    """Run one barcode calling pass in a child process.

    Read chunks are not reclaimed when barcode calling ends, leaving the main process
    holding ~2.5 GB that every later worker would inherit. A child process gives it back.
    """
    log_file, log_level = _get_log_params()
    with ProcessPoolExecutor(max_workers=1,
                             initializer=setup_worker_logging,
                             initargs=(log_file, log_level)) as proc:
        if threads == 1:
            future_res = proc.submit(process_single_thread, bc_args)
        else:
            future_res = proc.submit(process_in_parallel, bc_args)

    concurrent.futures.wait([future_res], return_when=concurrent.futures.ALL_COMPLETED)
    if future_res.exception() is not None:
        raise future_res.exception()


def detect_cell_barcodes(args, sample, input_files, threads):
    """Pass 1: extract barcode windows verbatim, then derive the cell barcode list.

    Returns the path of that list, which pass 2 uses as its whitelist.
    """
    if args.resume and os.path.exists(sample.raw_barcodes_done):
        logger.info("Cell barcodes were detected during the previous run, skipping")
        return sample.out_cell_barcodes_tsv
    if os.path.exists(sample.raw_barcodes_done):
        os.remove(sample.raw_barcodes_done)

    logger.info("Extracting barcodes from %d file(s) to detect cell barcodes" % len(input_files))
    # no FASTA and no barcode table: only the counts matter here, and in splitting modes it
    # is pass 2 whose extraction is kept
    raw_args = BarcodeCallingArgs(input_files, args.barcode_whitelist, args.mode,
                                  None, None, sample.aux_dir, threads,
                                  molecule=getattr(args, 'molecule', None),
                                  whitelist_matching=False,
                                  split_molecules=args.split_molecules)

    # one child process, so the large count table never lives in the main one
    log_file, log_level = _get_log_params()
    with ProcessPoolExecutor(max_workers=1,
                             initializer=setup_worker_logging,
                             initargs=(log_file, log_level)) as proc:
        future_res = proc.submit(detect_cell_barcode_list,
                                 raw_args,
                                 sample.out_cell_barcodes_tsv,
                                 get_barcode_length(args.mode),
                                 args.n_cells,
                                 args.n_cells_interval,
                                 sample.out_cell_barcodes_stats)

    concurrent.futures.wait([future_res], return_when=concurrent.futures.ALL_COMPLETED)
    if future_res.exception() is not None:
        raise future_res.exception()

    open(sample.raw_barcodes_done, "w").close()
    return sample.out_cell_barcodes_tsv


def call_barcodes(args):
    if args.barcoded_bam:
        logger.info("Barcodes will be extracted from BAM tags (%s/%s)" % (args.barcode_tag, args.umi_tag))
        return
    if args.barcoded_reads:
        # TODO barcoded files via YAML
        args.input_data.samples[0].barcoded_reads = args.barcoded_reads
        return
    for sample in args.input_data.samples:
        # Collect all input files for this sample
        input_files = [files[0] for files in sample.file_list]
        output_barcodes_list = [sample.barcodes_tsv + "_%d.tsv" % i for i in range(len(input_files))]
        barcodes_done_list = [sample.barcodes_done + "_%d.tsv" % i for i in range(len(input_files))]

        output_fasta_list = None
        new_reads = []
        if args.split_molecules:
            # minimap2 reads gzipped FASTA natively, so compressing costs nothing downstream
            fasta_suffix = ".fa.gz" if args.gzipped else ".fa"
            output_fasta_list = [sample.split_reads_fasta + "_%d%s" % (i, fasta_suffix)
                                 for i in range(len(input_files))]
            new_reads = [[fasta] for fasta in output_fasta_list]

        # Check if all files were already processed during resume
        all_done = all(os.path.exists(done) for done in barcodes_done_list)
        if all_done and args.resume:
            logger.info("Barcodes were called during the previous run, skipping")
            sample.barcoded_reads.extend(output_barcodes_list)
            if args.split_molecules:
                sample.file_list = new_reads
            continue

        # Remove existing done markers
        for barcodes_done in barcodes_done_list:
            if os.path.exists(barcodes_done):
                os.remove(barcodes_done)

        bc_threads = 1 if args.mode.enforces_single_thread() else args.threads
        barcode_files = args.barcode_whitelist
        if args.detect_cell_barcodes:
            barcode_files = [detect_cell_barcodes(args, sample, input_files, bc_threads)]

        bc_args = BarcodeCallingArgs(input_files, barcode_files, args.mode,
                                     output_barcodes_list, output_fasta_list, sample.aux_dir, bc_threads,
                                     molecule=getattr(args, 'molecule', None),
                                     split_molecules=args.split_molecules)
        logger.info("Detecting barcodes for %d file(s)" % len(input_files))
        run_barcode_calling(bc_args, bc_threads)

        # Mark all files as done and add to barcoded_reads
        for i, (input_file, output_barcodes, barcodes_done) in enumerate(zip(input_files, output_barcodes_list, barcodes_done_list)):
            sample.barcoded_reads.append(output_barcodes)
            open(barcodes_done, "w").close()
            logger.info("Processed %s, barcodes are stored in %s" % (input_file, output_barcodes))

        if args.split_molecules:
            logger.info("Reads were split during barcode calling")
            logger.info("The following files will be used instead of original reads %s " % ", ".join(map(lambda x: x[0], new_reads)))
            sample.file_list = new_reads
