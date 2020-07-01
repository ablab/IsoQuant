############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from src.isoform_assignment import *
from src.gene_info import *
from src.read_groups import *

logger = logging.getLogger('IsoQuant')


class AbstractCounter:
    def __init__(self, output_file_name, ignore_read_groups=False):
        self.ignore_read_groups = ignore_read_groups
        self.output_file_name = output_file_name

    def add_read_info(self, read_assignment):
        raise NotImplementedError()

    def dump(self):
        raise NotImplementedError()


class CompositeCounter:
    def __init__(self, counters):
        self.counters = counters

    def add_counters(self, counters):
        self.counters += counters

    def add_read_info(self, read_assignment):
        for p in self.counters:
            p.add_read_info(read_assignment)

    def dump(self):
        for p in self.counters:
            p.dump()


# count meta-features assigned to reads (genes or isoforms)
# get_feature_id --- function that returns feature id form IsoformMatch object
class AssignedFeatureCounter(AbstractCounter):
    def __init__(self, output_file_name, get_feature_id, ignore_read_groups=False):
        AbstractCounter.__init__(self, output_file_name, ignore_read_groups)
        self.get_feature_id = get_feature_id
        self.all_features = set()

        self.ambiguous_reads = 0
        self.not_assigned_reads = 0
        self.not_aligned_reads = 0
        # group_id -> (feature_id -> count)
        self.feature_counter = defaultdict(lambda: defaultdict(int))

    def add_read_info(self, read_assignment=None):
        #TODO: add __alignment_not_unique / __too_low_aQual ?
        if not read_assignment:
            self.not_aligned_reads += 1
        elif read_assignment.assignment_type == ReadAssignmentType.ambiguous:
            feature_ids = [self.get_feature_id(m) for m in read_assignment.isoform_matches]
            if set(feature_ids) == 1:  # different isoforms of same gene
                group_id = AbstractReadGrouper.default_group_id if self.ignore_read_groups else read_assignment.read_group
                self.feature_counter[group_id][feature_ids[0]] += 1
                self.all_features.add(feature_ids[0])
            else:
                self.ambiguous_reads += 1
        elif read_assignment.assignment_type == ReadAssignmentType.empty:
            self.not_assigned_reads += 1
        elif read_assignment.assignment_type == ReadAssignmentType.unique or\
                read_assignment.assignment_type == ReadAssignmentType.unique_minor_difference:
            feature_id = self.get_feature_id(read_assignment.isoform_matches[0])
            group_id = AbstractReadGrouper.default_group_id if self.ignore_read_groups else read_assignment.read_group
            self.feature_counter[group_id][feature_id] += 1
            self.all_features.add(feature_id)

    def add_unaligned(self, n_reads=1):
        self.not_aligned_reads += n_reads

    def dump(self):
        with open(self.output_file_name, "w") as f:
            f.write("feature_id\tgroup_id\tcount\n")
            for feature_id in self.all_features:
                for group_id in self.feature_counter.keys():
                    count = self.feature_counter[group_id][feature_id]
                    f.write("%s\t%s\t%d\n" % (feature_id, group_id, count))
            f.write("__ambiguous\t%d\n" % self.ambiguous_reads)
            f.write("__no_feature\t%d\n" % self.not_assigned_reads)
            f.write("__not_aligned\t%d\n" % self.not_aligned_reads)


def create_gene_counter(output_file_name, ignore_read_groups=False):
    return AssignedFeatureCounter(output_file_name, get_assigned_gene_id, ignore_read_groups)


def create_transcript_counter(output_file_name, ignore_read_groups=False):
    return AssignedFeatureCounter(output_file_name, get_assigned_transcript_id, ignore_read_groups)


# count simple features inclusion/exclusion (exons / introns)
class ProfileFeatureCounter(AbstractCounter):
    def __init__(self, output_file_name, ignore_read_groups=False):
        AbstractCounter.__init__(self, output_file_name, ignore_read_groups)
        # group_id -> (feature_id -> count)
        self.inclusion_feature_counter = defaultdict(lambda: defaultdict(int))
        self.exclusion_feature_counter = defaultdict(lambda: defaultdict(int))
        self.all_features = set()

    def add_read_info_from_profile(self, gene_feature_profile, feature_property_map,
                                   read_group = AbstractReadGrouper.default_group_id):
        for i in range(len(gene_feature_profile)):
            if gene_feature_profile[i] == 1:
                feature_id = feature_property_map[i].to_str()
                self.inclusion_feature_counter[read_group][feature_id] += 1
                self.all_features.add(feature_id)
            elif gene_feature_profile[i] == -1:
                feature_id = feature_property_map[i].to_str()
                self.exclusion_feature_counter[read_group][feature_id] += 1
                self.all_features.add(feature_id)

    def dump(self):
        with open(self.output_file_name, "w") as f:
            f.write(FeatureInfo.header() + "\tgroup_id\tinclude_counts\texclude_counts\n")
            all_groups = set(self.inclusion_feature_counter.keys())
            all_groups.update(self.exclusion_feature_counter.keys())
            all_groups = sorted(all_groups)

            for feature_id in sorted(self.all_features):
                for group_id in all_groups:
                    incl_count = self.inclusion_feature_counter[group_id][feature_id]
                    excl_count = self.exclusion_feature_counter[group_id][feature_id]
                    f.write("%s\t%s\t%d\t%d\n" % (feature_id, group_id, incl_count, excl_count))

    @staticmethod
    def is_valid(assignment):
        return assignment is not None and \
               hasattr(assignment, 'combined_profile') and assignment.combined_profile is not None and \
               hasattr(assignment, 'gene_info') and assignment.gene_info is not None

class ExonCounter(ProfileFeatureCounter):
    def __init__(self, output_file_name, ignore_read_groups=False):
        ProfileFeatureCounter.__init__(self, output_file_name, ignore_read_groups)

    def add_read_info(self, read_assignment):
        if not ProfileFeatureCounter.is_valid(read_assignment):
            return
        group_id = AbstractReadGrouper.default_group_id if self.ignore_read_groups else read_assignment.read_group
        self.add_read_info_from_profile(read_assignment.combined_profile.read_exon_profile.gene_profile,
                                        read_assignment.gene_info.exon_property_map,
                                        group_id)


class IntronCounter(ProfileFeatureCounter):
    def __init__(self, output_file_name, ignore_read_groups=False):
        ProfileFeatureCounter.__init__(self, output_file_name, ignore_read_groups)

    def add_read_info(self, read_assignment):
        if not ProfileFeatureCounter.is_valid(read_assignment):
            return
        group_id = AbstractReadGrouper.default_group_id if self.ignore_read_groups else read_assignment.read_group
        self.add_read_info_from_profile(read_assignment.combined_profile.read_intron_profile.gene_profile,
                                        read_assignment.gene_info.intron_property_map,
                                        group_id)

