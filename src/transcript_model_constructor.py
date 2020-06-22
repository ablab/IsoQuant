############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict

logger = logging.getLogger('IsoQuant')


class TranscriptModel:
    def __init__(self, transcript_id, reference_transcript, reference_gene, exon_blocks):
        self.transcript_id = transcript_id
        self.gene_id = reference_gene
        self.reference_transcript = reference_transcript
        self.reference_gene = reference_gene
        self.exon_blocks = exon_blocks


class TranscriptModelConstructor:
    def __init__(self, gene_info, read_assignment_storage):
        self.gene_info = gene_info
        self.read_assignment_storage = read_assignment_storage
        self.transcript_model_storage = []

    def process(self):
        pass

    def construct_isoform_groups(self):
        isoform_groups = defaultdict(list)
        for read_assignment in self.read_assignment_storage:
            for match in read_assignment.isoform_matches:
                isoform_groups[match.assigned_transcript].append(read_assignment)

        event_groups = defaultdict(lambda: defaultdict(int))
        for isoform_id in isoform_groups.keys():
            isoform_assignments = isoform_groups[isoform_id]
            for read_assignment in isoform_assignments:
                for match in read_assignment.isoform_matches:
                    significant_events = []



