############################################################################
# Copyright (c) 2022-2025 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import os
from .isoform_assignment import ReadAssignmentType

from .assignment_io import (
    NormalTmpFileAssignmentLoader,
    QuickTmpFileAssignmentLoader
)

logger = logging.getLogger('IsoQuant')


class BasicReadAssignmentLoader:
    def __init__(self, save_file_name):
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


class ReadAssignmentLoader:
    def __init__(self, save_file_name, gffutils_db, chr_record, multimapped_chr_dict):
        assert os.path.exists(save_file_name)
        self.save_file_name = save_file_name
        self.unpickler = NormalTmpFileAssignmentLoader(save_file_name, gffutils_db, chr_record)
        self.multimapped_chr_dict = multimapped_chr_dict

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
