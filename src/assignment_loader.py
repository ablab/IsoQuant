############################################################################
# Copyright (c) 2022-2025 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import os
from collections import defaultdict

import gffutils
from pyfaidx import Fasta, UnsupportedCompressionFormat

from .serialization import *
from .isoform_assignment import BasicReadAssignment, ReadAssignmentType
from .assignment_io import (
    NormalTmpFileAssignmentLoader,
    QuickTmpFileAssignmentLoader
)

logger = logging.getLogger('IsoQuant')


class ReadAssignmentLoader:
    def __init__(self, save_file_name, gffutils_db, chr_record, multimapped_chr_dict, filtered_read_set=None):
        logger.info("Loading read assignments from " + save_file_name)
        assert os.path.exists(save_file_name)
        self.save_file_name = save_file_name
        self.genedb = gffutils_db
        self.chr_record = chr_record
        self.unpickler = NormalTmpFileAssignmentLoader(save_file_name, gffutils_db, chr_record)
        self.multimapped_chr_dict = multimapped_chr_dict
        self.filtered_read_set = filtered_read_set

    def has_next(self):
        return self.unpickler.has_next()

    def get_next(self):
        if not self.unpickler.has_next():
            return None, None

        assert self.unpickler.is_gene_info()
        gene_info = self.unpickler.get_object()
        assignment_storage = []
        while self.unpickler.is_read_assignment():
            read_assignment = self.unpickler.get_object()
            if self.filtered_read_set and read_assignment.read_id not in self.filtered_read_set:
                continue
            if self.multimapped_chr_dict is not None and read_assignment.read_id in self.multimapped_chr_dict:
                resolved_assignment = None
                for a in self.multimapped_chr_dict[read_assignment.read_id]:
                    if a.assignment_id == read_assignment.assignment_id and a.chr_id == read_assignment.chr_id:
                        if resolved_assignment is not None:
                            logger.info("Duplicate read: %s %s %s" % (read_assignment.read_id, a.gene_id, a.chr_id))
                        resolved_assignment = a

                if not resolved_assignment:
                    logger.warning("Incomplete information on read %s" % read_assignment.read_id)
                    continue
                elif resolved_assignment.assignment_type == ReadAssignmentType.suspended:
                    continue
                else:
                    read_assignment.assignment_type = resolved_assignment.assignment_type
                    read_assignment.gene_assignment_type = resolved_assignment.gene_assignment_type
                    read_assignment.multimapper = resolved_assignment.multimapper
            assignment_storage.append(read_assignment)

        return gene_info, assignment_storage


class BasicReadAssignmentLoader:
    def __init__(self, save_file_name):
        logger.info("Loading read assignments from " + save_file_name)
        assert os.path.exists(save_file_name)
        self.save_file_name = save_file_name
        self.unpickler = QuickTmpFileAssignmentLoader(save_file_name)

    def has_next(self):
        return self.unpickler.has_next()

    def get_next(self):
        if not self.unpickler.has_next():
            return

        assert self.unpickler.is_gene_info()
        self.unpickler.get_object()

        while self.unpickler.is_read_assignment():
            yield self.unpickler.get_object()


def create_assignment_loader(chr_id, dump_filename, genedb, reference_fasta, reference_fai, use_filtered_reads=False):
    current_chr_record = Fasta(reference_fasta, indexname=reference_fai)[chr_id]
    multimapped_reads = defaultdict(list)
    # FIXME use proper file name constuction
    multimap_loader = open(dump_filename + "_multimappers_" + chr_id, "rb")
    list_size = read_int(multimap_loader)
    while list_size != TERMINATION_INT:
        for i in range(list_size):
            a = BasicReadAssignment.deserialize(multimap_loader)
            if a.chr_id == chr_id:
                multimapped_reads[a.read_id].append(a)
        list_size = read_int(multimap_loader)

    filtered_reads = set()
    if use_filtered_reads:
        # FIXME use proper file name construction
        for l in open(dump_filename + "_filtered_" + chr_id, "r"):
            filtered_reads.add(l.rstrip())

    chr_dump_file = dump_filename + "_" + chr_id

    if genedb:
        gffutils_db = gffutils.FeatureDB(genedb)
    else:
        gffutils_db = None

    return ReadAssignmentLoader(chr_dump_file, gffutils_db, current_chr_record, multimapped_reads, filtered_reads)