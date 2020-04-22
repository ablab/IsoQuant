from collections import defaultdict
from src.long_read_assigner import AssignmentType


class FeatureCounter:
    def __init__(self, gene_db, output_file_name):
        self.gene_db = gene_db
        self.ambiguous_reads = 0
        self.not_assigned_reads = 0
        self.feature_counter = defaultdict(int)
        self.output_file_name = output_file_name

    def add_read_info(self, read_assignment):
        #TODO: add __not_aligned / __alignment_not_unique / __too_low_aQual ?
        if read_assignment.assignment_type == AssignmentType.ambiguous:
            self.ambiguous_reads += 1
        elif read_assignment.assignment_type == AssignmentType.empty:
            self.not_assigned_reads += 1
        elif read_assignment.assignment_type == AssignmentType.unique:
            transcript_id = read_assignment.assigned_features[0]
            gene_id = self.gene_db[transcript_id]['gene_id'][0]
            self.feature_counter[gene_id] += 1

    def dump(self):
        with open(self.output_file_name, "w") as f:
            for gene_id, count in self.feature_counter.items():
                f.write("%s\t%d\n" % (gene_id, count))
            f.write("__ambiguous\t%d\n" % self.ambiguous_reads)
            f.write("__no_feature\t%d\n" % self.not_assigned_reads)

