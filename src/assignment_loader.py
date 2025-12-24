############################################################################
# Copyright (c) 2022-2025 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict

import gffutils
from pyfaidx import Fasta, UnsupportedCompressionFormat

from .serialization import *
from .file_naming import *
from .isoform_assignment import BasicReadAssignment, ReadAssignmentType
from .assignment_io import (
    NormalTmpFileAssignmentLoader,
    QuickTmpFileAssignmentLoader, GeneListTmpFileAssignmentLoader
)
from .gene_info import GeneList, GeneInfo

logger = logging.getLogger('IsoQuant')


class BasicAssignmentLoader:
    def __init__(self, save_file_name):
        logger.info("Loading read assignments from " + save_file_name)
        assert os.path.exists(save_file_name)
        self.save_file_name = save_file_name

    def has_next(self):
        raise NotImplementedError()

    def get_next(self):
        raise NotImplementedError()


class FullAssignmentLoader(BasicAssignmentLoader):
    def __init__(self, save_file_name, multimapped_chr_dict, filtered_read_set=None):
        BasicAssignmentLoader.__init__(self, save_file_name)
        self.multimapped_chr_dict = multimapped_chr_dict
        self.filtered_read_set = filtered_read_set

    def load_read_assignment(self, read_assignment):
        if self.filtered_read_set is not None and read_assignment.read_id not in self.filtered_read_set:
            return None
        if self.multimapped_chr_dict is not None and read_assignment.read_id in self.multimapped_chr_dict:
            resolved_assignment = None
            for a in self.multimapped_chr_dict[read_assignment.read_id]:
                if a.assignment_id == read_assignment.assignment_id and a.chr_id == read_assignment.chr_id:
                    if resolved_assignment is not None:
                        logger.info("Duplicate read: %s %s %s" % (read_assignment.read_id, a.gene_id, a.chr_id))
                    resolved_assignment = a

            if not resolved_assignment:
                logger.warning("Incomplete information on read %s" % read_assignment.read_id)
                return None
            elif resolved_assignment.assignment_type == ReadAssignmentType.suspended:
                return None
            else:
                read_assignment.assignment_type = resolved_assignment.assignment_type
                read_assignment.gene_assignment_type = resolved_assignment.gene_assignment_type
                read_assignment.multimapper = resolved_assignment.multimapper
        return read_assignment


class ReadAssignmentLoader(FullAssignmentLoader):
    def __init__(self, save_file_name, gffutils_db, chr_record, multimapped_chr_dict, filtered_read_set=None, string_pools):
        FullAssignmentLoader.__init__(self, save_file_name, multimapped_chr_dict, filtered_read_set)
        self.genedb = gffutils_db
        self.chr_record = chr_record
        self.unpickler = NormalTmpFileAssignmentLoader(save_file_name, gffutils_db, chr_record, string_pools)

    def has_next(self):
        return self.unpickler.has_next()

    def get_next(self):
        if not self.unpickler.has_next():
            return None, None

        assert self.unpickler.is_gene_info()
        gene_info = self.unpickler.get_object()
        assignment_storage = []
        while self.unpickler.is_read_assignment():
            read_assignment = self.load_read_assignment(self.unpickler.get_object())
            if read_assignment is not None:
                assignment_storage.append(read_assignment)

        return gene_info, assignment_storage


class MergingSimpleReadAssignmentLoader(FullAssignmentLoader):
    def __init__(self, save_file_name, multimapped_chr_dict, filtered_read_set=None, string_pools):
        FullAssignmentLoader.__init__(self, save_file_name, multimapped_chr_dict, filtered_read_set)
        self.unpickler = GeneListTmpFileAssignmentLoader(save_file_name, string_pools)
        self.current_gene_list = None

    def has_next(self):
        return self.unpickler.has_next()

    def get_next(self):
        if not self.unpickler.has_next():
            return None, None

        assignment_storage = []
        while self.unpickler.has_next():
            if self.current_gene_list is None:
                assert self.unpickler.is_gene_info()
                self.current_gene_list = self.unpickler.get_object()
            elif self.unpickler.is_gene_info():
                gene_list = self.unpickler.get_object()
                if self.current_gene_list.overlaps(gene_list):
                    self.current_gene_list.merge(gene_list)
                else:
                    self.current_gene_list = gene_list
                    return None, assignment_storage

            while self.unpickler.is_read_assignment():
                read_assignment = self.load_read_assignment(self.unpickler.get_object())
                if read_assignment is not None:
                    assignment_storage.append(read_assignment)

        return None, assignment_storage


class MergingReadAssignmentLoader(MergingSimpleReadAssignmentLoader):
    def __init__(self, save_file_name, gffutils_db, chr_record, multimapped_chr_dict, filtered_read_set=None):
        MergingSimpleReadAssignmentLoader.__init__(save_file_name, multimapped_chr_dict, filtered_read_set)
        self.genedb = gffutils_db
        self.chr_record = chr_record

    def _create_gene_info(self):
        if not self.current_gene_list.gene_id_set:
            gene_info = GeneInfo.from_region(self.current_gene_list.chr_id, self.current_gene_list.start,
                                             self.current_gene_list.end, self.current_gene_list.delta)
        else:
            gene_info = GeneInfo([self.genedb[gene_id] for gene_id in self.current_gene_list.gene_id_set], self.genedb,
                                 self.current_gene_list.delta)
        if self.chr_record:
            gene_info.set_reference_sequence(gene_info.all_read_region_start,
                                             gene_info.all_read_region_end,
                                             self.chr_record)
        return gene_info

    def get_next(self):
        if not self.unpickler.has_next():
            return None, None

        assignment_storage = []
        while self.unpickler.has_next():
            if self.current_gene_list is None:
                assert self.unpickler.is_gene_info()
                self.current_gene_list = self.unpickler.get_object()
            elif self.unpickler.is_gene_info():
                gene_list = self.unpickler.get_object()
                if self.current_gene_list.overlaps(gene_list):
                    self.current_gene_list.merge(gene_list)
                else:
                    gene_info = self._create_gene_info()
                    for a in assignment_storage:
                        a.gene_info = gene_info
                    self.current_gene_list = gene_list
                    return gene_info, assignment_storage

            while self.unpickler.is_read_assignment():
                read_assignment = self.load_read_assignment(self.unpickler.get_object())
                if read_assignment is not None:
                    assignment_storage.append(read_assignment)

        gene_info = self._create_gene_info()
        for a in assignment_storage:
            a.gene_info = gene_info
        return gene_info, assignment_storage


class BasicReadAssignmentLoader(BasicAssignmentLoader):
    def __init__(self, save_file_name, string_pools):
        BasicAssignmentLoader.__init__(self, save_file_name)
        self.unpickler = QuickTmpFileAssignmentLoader(save_file_name, string_pools)

    def has_next(self):
        return self.unpickler.has_next()

    def get_next(self):
        if not self.unpickler.has_next():
            return

        assert self.unpickler.is_gene_info()
        self.unpickler.get_object()

        while self.unpickler.is_read_assignment():
            yield self.unpickler.get_object()


def prepare_multimapped_reads(saves_prefix ,chr_id, string_pools):
    multimapped_reads = defaultdict(list)
    multimap_loader = open(multimappers_file_name(saves_prefix ,chr_id), "rb")
    list_size = read_int(multimap_loader)
    while list_size != TERMINATION_INT:
        for i in range(list_size):
            a = BasicReadAssignment.deserialize(multimap_loader, string_pools)
            if a.chr_id == chr_id:
                multimapped_reads[a.read_id].append(a)
        list_size = read_int(multimap_loader)
    return multimapped_reads


def prepare_read_filter(chr_id, saves_prefix, use_filtered_reads):
    if not use_filtered_reads:
        return None
    filtered_reads = set()
    for l in open(filtered_reads_file_name(saves_prefix, chr_id), "r"):
        filtered_reads.add(l.rstrip())
    return filtered_reads


def load_genedb(genedb):
    if genedb:
        return gffutils.FeatureDB(genedb)
    return None


def create_assignment_loader(chr_id, saves_prefix, genedb, reference_fasta, reference_fai, use_filtered_reads=False, string_pools):
    current_chr_record = Fasta(reference_fasta, indexname=reference_fai)[chr_id]
    multimapped_reads = prepare_multimapped_reads(saves_prefix, chr_id, string_pools)
    filtered_reads = prepare_read_filter(chr_id, saves_prefix, use_filtered_reads)
    gffutils_db = load_genedb(genedb)
    chr_dump_file = saves_file_name(saves_prefix, chr_id)

    return ReadAssignmentLoader(chr_dump_file, gffutils_db, current_chr_record, multimapped_reads, filtered_reads, string_pools)


def create_merging_assignment_loader(chr_id, saves_prefix, use_filtered_reads=False, string_pools):
    multimapped_reads = prepare_multimapped_reads(saves_prefix, chr_id, string_pools)
    filtered_reads = prepare_read_filter(chr_id, saves_prefix, use_filtered_reads)
    chr_dump_file = saves_file_name(saves_prefix, chr_id)

    return MergingSimpleReadAssignmentLoader(chr_dump_file, multimapped_reads, filtered_reads, string_pools)
