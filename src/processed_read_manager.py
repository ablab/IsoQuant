############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict

from .file_naming import saves_file_name, multimappers_file_name
from .serialization import *
from .isoform_assignment import BasicReadAssignment, ReadAssignmentType, ReadAssignment
from .multimap_resolver import MultimapResolver
from .assignment_loader import BasicReadAssignmentLoader
from .string_pools import StringPoolManager

logger = logging.getLogger('IsoQuant')


def prepare_multimapper_dict(chr_ids, sample, multimappers_counts, all_chr_ids=None):
    multimapped_reads = defaultdict(list)
    unique_assignments = 0
    polya_unique_assignments = 0

    # Build string pools for deserialization
    string_pools = StringPoolManager()

    # Build chromosome pool from full list (or fallback to current chr_ids)
    pool_chr_ids = all_chr_ids if all_chr_ids else chr_ids
    string_pools.build_chromosome_pool(pool_chr_ids)

    for chr_id in chr_ids:
        chr_dump_file = saves_file_name(sample.out_raw_file, chr_id)
        loader = BasicReadAssignmentLoader(chr_dump_file, string_pools)
        while loader.has_next():
            for read_assignment in loader.get_next():
                if read_assignment is None:
                    continue
                if (read_assignment.read_id not in multimappers_counts or
                        multimappers_counts[read_assignment.read_id] == 1):
                    unique_assignments += 1
                    polya_unique_assignments += 1 if read_assignment.polyA_found else 0
                    continue
                multimapped_reads[read_assignment.read_id].append(read_assignment)
    return multimapped_reads, unique_assignments, polya_unique_assignments



def resolve_multimappers(chr_ids, sample, multimapped_reads, strategy):
    multimap_resolver = MultimapResolver(strategy)
    multimap_dumper = {}
    for chr_id in chr_ids:
        multimap_dumper[chr_id] = open(multimappers_file_name(sample.out_raw_file, chr_id), 'wb')
    total_assignments = 0
    polya_assignments = 0

    for assignment_list in multimapped_reads.values():
        if len(assignment_list) > 1:
            assignment_list = multimap_resolver.resolve(assignment_list)
            resolved_lists = defaultdict(list)
            for a in assignment_list:
                resolved_lists[a.chr_id].append(a)
            for chr_id in resolved_lists.keys():
                write_list(resolved_lists[chr_id], multimap_dumper[chr_id], BasicReadAssignment.serialize)

        for a in assignment_list:
            if a.assignment_type != ReadAssignmentType.suspended:
                total_assignments += 1
                if a.polyA_found:
                    polya_assignments += 1

    for chr_id in chr_ids:
        write_int(TERMINATION_INT, multimap_dumper[chr_id])
        multimap_dumper[chr_id].close()

    return total_assignments, polya_assignments


class ProcessedReadsManager:
    def __init__(self, sample, multimap_strategy, all_chr_ids=None):
        self.sample = sample
        self.multimap_strategy = multimap_strategy
        self.all_chr_ids = all_chr_ids

    def add_read(self, read_assignment: ReadAssignment):
        raise NotImplementedError()

    def load_read(self, read_assignment: ReadAssignment):
        raise NotImplementedError()

    def finalize(self, chr_id):
        pass

    def merge(self, processed_reads, chr_id):
        raise NotImplementedError()

    def resolve(self):
        raise NotImplementedError()


class ProcessedReadsManagerHighMemory(ProcessedReadsManager):
    def __init__(self, sample, multimap_strategy, all_chr_ids=None):
        ProcessedReadsManager.__init__(self, sample, multimap_strategy, all_chr_ids)
        self.read_storage = []
        self.multimapped_reads = defaultdict(list)
        self.chr_ids = set()

    def add_read(self, read_assignment):
        self.read_storage.append(BasicReadAssignment(read_assignment))

    def load_read(self, read_assignment):
        self.read_storage.append(read_assignment)

    def merge(self, processed_reads, chr_id):
        self.chr_ids.add(chr_id)
        for basic_read_assignment in processed_reads.read_storage:
            self.multimapped_reads[basic_read_assignment.read_id].append(basic_read_assignment)

    def resolve(self):
        return resolve_multimappers(self.chr_ids, self.sample, self.multimapped_reads, self.multimap_strategy)


class ProcessedReadsManagerNormalMemory(ProcessedReadsManager):
    def __init__(self, sample, multimap_strategy, all_chr_ids=None):
        ProcessedReadsManager.__init__(self, sample, multimap_strategy, all_chr_ids)
        self.read_storage = []
        self.multimappers_counts = defaultdict(int)
        self.multimapped_reads = defaultdict(list)
        self.chr_ids = set()

    def add_read(self, read_assignment):
        self.read_storage.append(read_assignment.read_id)

    def load_read(self, read_assignment):
        self.read_storage.append(read_assignment.read_id)

    def merge(self, processed_reads, chr_id):
        self.chr_ids.add(chr_id)
        for read_id in processed_reads.read_storage:
            self.multimappers_counts[read_id] += 1

    def resolve(self):
        multimapped_reads, unique_assignments, polya_unique_assignments \
            = prepare_multimapper_dict(self.chr_ids, self.sample, self.multimappers_counts, self.all_chr_ids)
        total_assignments, polya_assignments = resolve_multimappers(self.chr_ids, self.sample, multimapped_reads,
                                                                    self.multimap_strategy)
        total_assignments += unique_assignments
        polya_assignments += polya_unique_assignments
        return total_assignments, polya_assignments


class ProcessedReadsManagerNoSecondary(ProcessedReadsManager):
    def __init__(self, sample, multimap_strategy, all_chr_ids=None):
        ProcessedReadsManager.__init__(self, sample, multimap_strategy, all_chr_ids)
        self.read_storage = defaultdict(int)
        self.total_assignments = 0
        self.polya_assignments = 0

    def add_read(self, read_assignment):
        self.read_storage[read_assignment.read_id] += 1

    def load_read(self, read_assignment):
        self.add_read(read_assignment)

    def finalize(self, chr_id):
        multimapped_reads_dict, unique_assignments, polya_unique_assignments   \
            = prepare_multimapper_dict([chr_id], self.sample, self.read_storage, self.all_chr_ids)
        self.total_assignments, self.polya_assignments = resolve_multimappers([chr_id], self.sample,
                                                                            multimapped_reads_dict,
                                                                            self.multimap_strategy)
        self.total_assignments += unique_assignments
        self.polya_assignments += polya_unique_assignments
        self.read_storage.clear()

    def merge(self, processed_reads, _):
        self.total_assignments += processed_reads.total_assignments
        self.polya_assignments += processed_reads.polya_assignments

    def resolve(self):
        return self.total_assignments, self.polya_assignments
