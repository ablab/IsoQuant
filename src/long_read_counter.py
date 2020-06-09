from collections import defaultdict
from src.isoform_assignment import *


# count meta-features assigned to reads (genes or isoforms)
# get_feature_id --- function that returns feature id form IsoformMatch object
class AssignedFeatureCounter:
    def __init__(self, output_file_name, get_feature_id):
        self.get_feature_id = get_feature_id

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
            feature_ids = [self.get_feature_id(m) for m in read_assignment.isoform_matches]
            if set(feature_ids) == 1:  # different isoforms of same gene
                self.feature_counter[feature_ids[0]] += 1
            else:
                self.ambiguous_reads += 1
        elif read_assignment.assignment_type == ReadAssignmentType.empty:
            self.not_assigned_reads += 1
        elif read_assignment.assignment_type == ReadAssignmentType.unique or\
                read_assignment.assignment_type == ReadAssignmentType.unique_minor_difference:
            feature_id = self.get_feature_id(read_assignment.isoform_matches[0])
            self.feature_counter[feature_id] += 1

    def add_unaligned(self, n_reads=1):
        self.not_aligned_reads += n_reads

    def dump(self):
        with open(self.output_file_name, "w") as f:
            f.write("feature_id\tcount\n")
            for feature_id, count in self.feature_counter.items():
                f.write("%s\t%d\n" % (feature_id, count))
            f.write("__ambiguous\t%d\n" % self.ambiguous_reads)
            f.write("__no_feature\t%d\n" % self.not_assigned_reads)
            f.write("__not_aligned\t%d\n" % self.not_aligned_reads)


def get_gene_counter(output_file_name):
    return AssignedFeatureCounter(output_file_name, get_assigned_gene_id)


def get_transcript_counter(output_file_name):
    return AssignedFeatureCounter(output_file_name, get_assigned_transcript_id)

# count simple features inclusion/exclusion (exons / introns)
class ProfileFeatureCounter:
    def __init__(self, output_file_name):
        self.inclusion_feature_counter = defaultdict(int)
        self.exclusion_feature_counter = defaultdict(int)
        self.output_file_name = output_file_name

    def add_read_info(self, gene_feature_profile, feature_property_map):
        for i in range(len(gene_feature_profile)):
            if gene_feature_profile[i] == 1:
                self.inclusion_feature_counter[feature_property_map[i]] += 1
            elif gene_feature_profile[i] == -1:
                self.exclusion_feature_counter[feature_property_map[i]] += 1

    def dump(self):
        with open(self.output_file_name, "w") as f:
            f.write("feature_id\tinclude_counts\texclude_counts\n")
            all_features = set(self.inclusion_feature_counter.keys())
            all_features.update(self.exclusion_feature_counter.keys())
            for feature_id in sorted(all_features):
                f.write("%s\t%d\t%d\n" % (feature_id,
                                          self.inclusion_feature_counter[feature_id],
                                          self.exclusion_feature_counter[feature_id]))
