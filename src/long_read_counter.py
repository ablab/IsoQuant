from collections import defaultdict
from src.long_read_assigner import ReadAssignmentType


class FeatureCounter:
    def __init__(self, gene_db, output_file_name):
        self.gene_db = gene_db
        self.ambiguous_reads = 0
        self.not_assigned_reads = 0
        self.not_aligned_reads = 0
        self.feature_counter = defaultdict(int)
        self.output_file_name = output_file_name

    def add_read_info(self, read_assignment=None):
        #TODO: add __alignment_not_unique / __too_low_aQual ?
        if not read_assignment:
            self.not_aligned_reads += 1
        elif read_assignment.assignment_type == ReadAssignmentType.ambiguous:
            transcript_ids = [m.assigned_transcript for m in read_assignment.isoform_matches]
            gene_ids = [self.gene_db[t]['gene_id'][0] for t in transcript_ids]
            if set(gene_ids) == 1:  # different isoforms of same gene
                self.feature_counter[gene_ids[0]] += 1
            else:
                self.ambiguous_reads += 1
        elif read_assignment.assignment_type == ReadAssignmentType.empty:
            self.not_assigned_reads += 1
        elif read_assignment.assignment_type == ReadAssignmentType.unique or\
                read_assignment.assignment_type == ReadAssignmentType.minor:
            transcript_id = read_assignment.isoform_matches[0].assigned_transcript
            gene_id = self.gene_db[transcript_id]['gene_id'][0]
            self.feature_counter[gene_id] += 1

    def add_unaligned(self, n_reads=1):
        self.not_aligned_reads += n_reads

    def dump(self):
        with open(self.output_file_name, "w") as f:
            for gene_id, count in self.feature_counter.items():
                f.write("%s\t%d\n" % (gene_id, count))
            f.write("__ambiguous\t%d\n" % self.ambiguous_reads)
            f.write("__no_feature\t%d\n" % self.not_assigned_reads)
            f.write("__not_aligned\t%d\n" % self.not_aligned_reads)

