############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import pysam
from Bio import SeqIO

from src.linked_read_assigner import *
from src.read_groups import *
from src.exon_corrector import *
from src.alignment_info import *

logger = logging.getLogger('IsoQuant')


class LinkedReadAlignmentProcessor:
    """ class for aggregating all assignment information

    Parameters
    ----------
    gene_info
    bams : pair bam object, filename
    params
    printer
    counter
    """

    def __init__(self, gene_info, bam_pairs, params, chr_record=None):
        self.gene_info = gene_info
        self.bam_pairs = bam_pairs
        self.params = params
        self.chr_record = chr_record
        self.assigner = LinkedReadAssigner(self.gene_info, self.params)
        self.profile_constructor = CombinedProfileConstructor(gene_info, params)
        self.assignment_storage = []

    def process(self):
        self.assignment_storage = []
        self.gene_info.all_read_region_start = self.gene_info.start
        self.gene_info.all_read_region_end = self.gene_info.end

        for bam_pair in self.bam_pairs:
            self.process_single_file(bam_pair)

        if self.params.needs_reference:
            self.gene_info.all_read_region_start -= self.params.upstream_region_len
            self.gene_info.all_read_region_end += self.params.upstream_region_len
            self.gene_info.reference_region = \
                str(self.chr_record[self.gene_info.all_read_region_start - 1:self.gene_info.all_read_region_end + 1].seq)
            self.gene_info.canonical_sites = {}
        return self.assignment_storage

    def process_single_file(self, bam_pair):
        processed_reads = set()
        bamfile_in = bam_pair[0]
        for genic_region in self.gene_info.regions_for_bam_fetch:
            # FIXME: temporary solution - process gene outside
            to_fetch_start = max(0, genic_region[0] - 100)
            to_fetch_end = min(bamfile_in.get_reference_length(self.gene_info.chr_id), genic_region[1] + 100)
            for alignment in bamfile_in.fetch(self.gene_info.chr_id, to_fetch_start, to_fetch_end):
                read_id = alignment.query_name
                if alignment.reference_id == -1:
                    self.assignment_storage.append(ReadAssignment(read_id, None))
                    continue
                if alignment.is_supplementary:
                    continue
                if self.params.no_secondary and alignment.is_secondary:
                    continue

                # logger.debug("=== Processing read " + read_id + " ===")
                alignment_info = AlignmentInfo(alignment)

                if not alignment_info.read_exons:
                    logger.warning("Read %s has no aligned exons" % read_id)
                    continue
                read_tuple = (read_id, alignment_info.read_start, alignment_info.read_end)
                if read_tuple in processed_reads:
                    continue
                processed_reads.add(read_tuple)

                # TODO Chi-Lam: insert cloud processing here
                raise NotImplementedError()

