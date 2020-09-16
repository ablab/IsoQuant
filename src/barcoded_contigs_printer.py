############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
from src.assignment_io import *

logger = logging.getLogger('IsoQuant')


class BarcodedContigsPrinter(AbstractAssignmentPrinter):
    def __init__(self, barcode_table_file, output_file_name, params):
        AbstractAssignmentPrinter.__init__(self, output_file_name, params)
        self.header = "#barcode\tcontig_id\tisoform_id\tassignment_type\tassignment_events\n"
        self.output_file.write(self.header)

        self.barcode_table = defaultdict(list)
        with open(barcode_table_file) as barcode_table_handle:
            for l in barcode_table_handle:
                tokens = l.strip().split('\t')
                if l and len(tokens) != 2:
                    logger.warning("Malformed barcode table line: " + l.strip() + ", skipping")
                    continue
                barcode = tokens[0]
                if barcode in self.barcode_table:
                    logger.warning("Duplicate barcode: " + barcode)
                    continue
                self.barcode_table[barcode] = tokens[1].split(',')

    def add_read_info(self, read_assignment):
        # barcode, isoform, type, contig
        if read_assignment is None:
            return
        if self.assignment_checker is None or not self.assignment_checker.check(read_assignment):
            return

        if read_assignment.assignment_type is None or read_assignment.isoform_matches is None:
            line = read_assignment.read_id + "\t.\t.\t."
        else:
            assigned_transcripts = [str(m.assigned_transcript) for m in read_assignment.isoform_matches]
            for m in read_assignment.isoform_matches:
                for x in m.match_subclassifications:
                    if not hasattr(x, "event_type"):
                        logger.debug(x)
            match_events = ",".join(["+".join([x.event_type.name for x in m.match_subclassifications])
                                     for m in read_assignment.isoform_matches])
            if not match_events:
                match_events = "."
            line = read_assignment.read_id + "\t" + ",".join(assigned_transcripts) + "\t" \
                   + read_assignment.assignment_type.name + "\t" + match_events

        line += "\n"
        for barcode in self.barcode_table[read_assignment.read_id]:
            self.output_file.write(barcode + "\t" + line)

