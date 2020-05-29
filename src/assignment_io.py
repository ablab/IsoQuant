############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from src.common import *
from src.long_read_assigner import *

logger = logging.getLogger('IsoQuant')

class PrintAllFunctor:
    def check(sefl, assignment):
        return True


class PrintOnlyFunctor:
    def __init__(self, allowed_types):
        if isinstance(allowed_types, list):
            self.allowed_types = set(allowed_types)
        elif isinstance(allowed_types, set):
            self.allowed_types = allowed_types
        else:
            self.allowed_types = set([allowed_types])

    def check(self, assignment):
        return assignment.assignment_type in self.allowed_types


class AbstractAssignmentPrinter:
    def __init__(self, output_file_name, params, format, assignment_checker=PrintAllFunctor()):
        self.params = params
        self.format = format
        self.assignment_checker = assignment_checker
        self.output_file = open(output_file_name, "w")

    def __del__(self):
        self.output_file.close()

    def add_read_info(self, read_assignment, combined_read_profile = None):
        raise NotImplementedError()

    def flush(self):
        self.output_file.flush()


class ReadAssignmentCompositePrinter:
    def __init__(self, printers = []):
        self.printers = printers

    def add_read_info(self, read_assignment, mapping_read_profile=None):
        for p in self.printers:
            p.add_read_info(read_assignment, mapping_read_profile)

    def flush(self):
        for p in self.printers:
            p.flush()

# TODO: reformat output, make singe file
class BasicTSVAssignmentPrinter(AbstractAssignmentPrinter):
    def __init__(self, output_file_name, params, format = "TSV", assignment_checker=PrintAllFunctor()):
        AbstractAssignmentPrinter.__init__(self, output_file_name, params, format, assignment_checker)
        self.header = "#read_id\tisoform_id\tassignment_type\tassignment_events"
        if self.params.print_additional_info:
            self.header += "\taligned_blocks\tintron_profile\tsplit_exon_profile"
        self.header += "\n"
        self.output_file.write(self.header)

    def add_read_info(self, read_assignment, combined_read_profile = None):
        if not self.assignment_checker.check(read_assignment):
            return
        assigned_transcripts = [m.assigned_transcript for m in read_assignment.isoform_matches]
        match_events = ["+".join(map(lambda x: x.name, m.match_subclassifications)) for m in read_assignment.isoform_matches]
        line = read_assignment.read_id  + "\t" + ",".join(assigned_transcripts) + "\t" \
                + read_assignment.assignment_type.name + "\t" + ",".join(match_events)
        if self.params.print_additional_info:
            if combined_read_profile is None:
                line += "\t.\t.\t."
            else:
                line += "\t" + range_list_to_str(combined_read_profile.read_split_exon_profile.read_features) + "\t" + \
                    list_to_str(combined_read_profile.read_intron_profile.gene_profile) + "\t" + \
                        list_to_str(combined_read_profile.read_split_exon_profile.gene_profile)
        line += "\n"
        self.output_file.write(line)
