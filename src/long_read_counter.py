############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from src.isoform_assignment import *
from src.gene_info import *
from src.read_groups import *

logger = logging.getLogger('IsoQuant')


@unique
class CountingStrategies(Enum):
    unique_only = 1
    with_ambiguous = 2
    inconsistent_all = 10
    inconsistent_minor = 11
    all = 20
    all_minor = 21


COUNTING_STRATEGIES = ["unique_only", "with_ambiguous", "with_inconsistent", "all"]


class ReadWeightCounter:
    def __init__(self, strategy, gene_counting=True):
        if strategy == "unique_only":
            self.strategy = CountingStrategies.unique_only
        elif strategy == "with_ambiguous":
            self.strategy = CountingStrategies.with_ambiguous
        elif strategy == "all":
            if gene_counting:
                self.strategy = CountingStrategies.all
            else:
                self.strategy = CountingStrategies.all_minor
        elif strategy == "with_inconsistent":
            if gene_counting:
                self.strategy = CountingStrategies.inconsistent_all
            else:
                self.strategy = CountingStrategies.inconsistent_minor
        else:
            raise KeyError("Unknown quantification strategy: " + strategy)

    def process_inconsistent(self, isoform_match, feature_count):
        if self.strategy == CountingStrategies.unique_only or self.strategy == CountingStrategies.with_ambiguous \
                or feature_count > 1:
            return 0.0
        elif self.strategy == CountingStrategies.all or self.strategy == CountingStrategies.inconsistent_all:
            return 1.0
        elif self.strategy == CountingStrategies.all_minor or self.strategy == CountingStrategies.inconsistent_minor:
            for m in isoform_match.match_subclassifications:
                if m.event_type not in nonintronic_events:
                    return 0.0
            return 1.0
        return 0.0

    def process_ambiguous(self, feature_count):
        if self.strategy in {CountingStrategies.with_ambiguous,
                             CountingStrategies.all_minor,
                             CountingStrategies.all}:
            return 1.0 / float(feature_count)
        else:
            return 0.0


class AbstractCounter:
    def __init__(self, output_prefix, ignore_read_groups=False, output_zeroes=True):
        self.ignore_read_groups = ignore_read_groups
        self.output_counts_file_name = output_prefix + "_counts.tsv"
        open(self.output_counts_file_name, "w").close()
        self.output_tpm_file_name = output_prefix + "_tpm.tsv"
        self.output_zeroes = output_zeroes
        self.output_stats_file_name = None

    def add_read_info(self, read_assignment):
        raise NotImplementedError()

    def add_read_info_raw(self, read_id, feature_ids, group_id=AbstractReadGrouper.default_group_id):
        raise NotImplementedError()

    def add_confirmed_features(self, features):
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

    def add_read_info_raw(self, read_id, feature_ids, group_id=AbstractReadGrouper.default_group_id):
        for p in self.counters:
            p.add_read_info_raw(read_id, feature_ids, group_id)

    def add_confirmed_features(self, features):
        for p in self.counters:
            p.add_confirmed_features(features)

    def dump(self):
        for p in self.counters:
            p.dump()


# count meta-features assigned to reads (genes or isoforms)
# get_feature_id --- function that returns feature id form IsoformMatch object
class AssignedFeatureCounter(AbstractCounter):
    def __init__(self, output_prefix, get_feature_id, read_grouper, read_counter,
                 ignore_read_groups=False, output_zeroes=True):
        AbstractCounter.__init__(self, output_prefix, ignore_read_groups, output_zeroes)
        self.get_feature_id = get_feature_id
        self.all_features = set()
        self.all_groups = read_grouper.read_groups if read_grouper else []
        self.read_counter = read_counter

        self.ambiguous_reads = 0
        self.not_assigned_reads = 0
        self.not_aligned_reads = 0
        # group_id -> (feature_id -> count)
        self.feature_counter = defaultdict(lambda: defaultdict(float))
        self.confirmed_features = set()
        self.output_stats_file_name = self.output_counts_file_name + ".stats"

    def add_read_info(self, read_assignment=None):
        # TODO: add __alignment_not_unique / __too_low_aQual ?
        if not read_assignment:
            self.not_aligned_reads += 1
        if read_assignment.assignment_type in [ReadAssignmentType.noninformative, ReadAssignmentType.intergenic] or \
                not read_assignment.isoform_matches:
            self.not_assigned_reads += 1
        elif read_assignment.assignment_type == ReadAssignmentType.ambiguous:
            feature_ids = set([self.get_feature_id(m) for m in read_assignment.isoform_matches])
            group_id = AbstractReadGrouper.default_group_id if self.ignore_read_groups else read_assignment.read_group

            if len(feature_ids) == 1:  # different isoforms of same gene
                feature_id = list(feature_ids)[0]
                self.feature_counter[group_id][feature_id] += 1.0
                self.all_features.add(feature_id)
            else:
                for feature_id in feature_ids:
                    count_value = self.read_counter.process_ambiguous(len(feature_ids))
                    self.feature_counter[group_id][feature_id] += count_value
                    if count_value > 0:
                        self.all_features.add(feature_id)
                self.ambiguous_reads += 1
        elif read_assignment.assignment_type == ReadAssignmentType.inconsistent:
            feature_ids = set([self.get_feature_id(m) for m in read_assignment.isoform_matches])
            group_id = AbstractReadGrouper.default_group_id if self.ignore_read_groups else read_assignment.read_group
            count_value = self.read_counter.process_inconsistent(read_assignment.isoform_matches[0], len(feature_ids))
            if count_value > 0:
                for feature_id in feature_ids:
                    self.feature_counter[group_id][feature_id] += count_value
                    self.all_features.add(feature_id)
        elif read_assignment.assignment_type == ReadAssignmentType.unique or\
                read_assignment.assignment_type == ReadAssignmentType.unique_minor_difference:
            feature_id = self.get_feature_id(read_assignment.isoform_matches[0])
            group_id = AbstractReadGrouper.default_group_id if self.ignore_read_groups else read_assignment.read_group
            self.feature_counter[group_id][feature_id] += 1.0
            self.all_features.add(feature_id)
            transcript_id = read_assignment.isoform_matches[0].assigned_transcript
            if len(read_assignment.gene_info.all_isoforms_introns[transcript_id]) == 0 or len(read_assignment.corrected_exons) > 1:
                # only add features to confirmed ones if it has unique spliced match
                self.confirmed_features.add((group_id, feature_id))

    def add_read_info_raw(self, read_id, feature_ids, group_id=AbstractReadGrouper.default_group_id):
        if self.ignore_read_groups:
            group_id = AbstractReadGrouper.default_group_id
        if not read_id:
            self.not_aligned_reads += 1
        elif not feature_ids:
            self.not_assigned_reads += 1
        elif len(feature_ids) > 1:
            self.ambiguous_reads += 1
            for feature_id in feature_ids:
                count_value = self.read_counter.process_ambiguous(len(feature_ids))
                self.feature_counter[group_id][feature_id] += count_value
                # self.confirmed_features.add((group_id, feature_id))
                self.all_features.add(feature_id)
            else:
                self.ambiguous_reads += 1
        else:
            self.feature_counter[group_id][feature_ids[0]] += 1
            # self.confirmed_features.add((group_id, feature_ids[0]))
            self.all_features.add(feature_ids[0])

    def add_unaligned(self, n_reads=1):
        self.not_aligned_reads += n_reads

    def format_header(self, all_groups, value_name="count"):
        if self.ignore_read_groups:
            return "#feature_id\t%s\n" % value_name
        else:
            return "#feature_id\t" + "\t".join(all_groups) + "\n"

    def add_confirmed_features(self, features):
        for feature_id in features:
            for group_id in self.feature_counter.keys():
                self.confirmed_features.add((group_id, feature_id))

    def dump(self):
        total_counts = defaultdict(float)
        all_features = sorted(self.all_features)
        all_groups = sorted(self.all_groups) if self.all_groups else sorted(self.feature_counter.keys())

        for group_id in all_groups:
            for feature_id in all_features:
                if (group_id, feature_id) in self.confirmed_features:
                    continue
                self.feature_counter[group_id][feature_id] = 0.0

        with open(self.output_counts_file_name, "w") as f:
            f.write(self.format_header(all_groups))
            for feature_id in all_features:
                if self.ignore_read_groups:
                    count = self.feature_counter[all_groups[0]][feature_id]
                    if not self.output_zeroes and count == 0:
                        continue
                    total_counts[all_groups[0]] += count
                    f.write("%s\t%.2f\n" % (feature_id, count))
                else:
                    for group_id in all_groups:
                        count = self.feature_counter[group_id][feature_id]
                        total_counts[group_id] += count
                    count_values = [self.feature_counter[group_id][feature_id] for group_id in all_groups]
                    f.write("%s\t%s\n" % (feature_id, "\t".join(["%.2f" % c for c in count_values])))

        if self.ignore_read_groups:
            with open(self.output_counts_file_name + ".stats", "w") as f:
                f.write("__ambiguous\t%d\n" % self.ambiguous_reads)
                f.write("__no_feature\t%d\n" % self.not_assigned_reads)
                f.write("__not_aligned\t%d\n" % self.not_aligned_reads)

    def convert_counts_to_tpm(self):
        total_counts = defaultdict(float)
        with open(self.output_counts_file_name) as f:
            for line in f:
                if line.startswith('_'): break
                fs = line.split()
                if line.startswith('#'): continue
                if self.ignore_read_groups:
                    total_counts[AbstractReadGrouper.default_group_id] += float(fs[1])
                else:
                    for j in range(len(fs) - 1):
                        total_counts[j] += float(fs[j + 1])

        scale_factors = {}
        for group_id in total_counts.keys():
            scale_factors[group_id] = 1000000.0 / total_counts[group_id] if total_counts[group_id] > 0 else 1.0
            logger.debug("Scale factor for group %s = %.2f" % (group_id, scale_factors[group_id]))

        with open(self.output_tpm_file_name, "w") as outf:
            with open(self.output_counts_file_name) as f:
                for line in f:
                    if line.startswith('_'): break
                    if line.startswith('#'):
                        outf.write(line.replace("count", "TPM"))
                        continue
                    fs = line.split()
                    if self.ignore_read_groups:
                        feature_id, count = fs[0], float(fs[1])
                        tpm = scale_factors[AbstractReadGrouper.default_group_id] * count
                        if not self.output_zeroes and tpm == 0:
                            continue
                        outf.write("%s\t%.6f\n" % (feature_id, tpm))
                    else:
                        feature_id, counts = fs[0], list(map(float, fs[1:]))
                        tpm_values = [scale_factors[i] * counts[i] for i in range(len(scale_factors))]
                        outf.write("%s\t%s\n" % (feature_id, "\t".join(["%.6f" % c for c in tpm_values])))


def create_gene_counter(output_file_name, strategy, read_grouper=None, ignore_read_groups=False, output_zeroes=True):
    read_weight_counter = ReadWeightCounter(strategy, gene_counting=True)
    return AssignedFeatureCounter(output_file_name, get_assigned_gene_id,
                                  read_grouper, read_weight_counter,
                                  ignore_read_groups, output_zeroes)


def create_transcript_counter(output_file_name, strategy, read_grouper=None, ignore_read_groups=False, output_zeroes=True):
    read_weight_counter = ReadWeightCounter(strategy, gene_counting=False)
    return AssignedFeatureCounter(output_file_name, get_assigned_transcript_id,
                                  read_grouper, read_weight_counter,
                                  ignore_read_groups, output_zeroes)


# count simple features inclusion/exclusion (exons / introns)
class ProfileFeatureCounter(AbstractCounter):
    def __init__(self, output_prefix, ignore_read_groups=False):
        AbstractCounter.__init__(self, output_prefix, ignore_read_groups)
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
        with open(self.output_counts_file_name, "w") as f:
            f.write(FeatureInfo.header() + "\tgroup_id\tinclude_counts\texclude_counts\n")
            all_groups = set(self.inclusion_feature_counter.keys())
            all_groups.update(self.exclusion_feature_counter.keys())
            all_groups = sorted(all_groups)

            for feature_id in sorted(self.all_features):
                for group_id in all_groups:
                    incl_count = self.inclusion_feature_counter[group_id][feature_id]
                    excl_count = self.exclusion_feature_counter[group_id][feature_id]
                    f.write("%s\t%s\t%d\t%d\n" % (feature_id, group_id, incl_count, excl_count))

    def convert_counts_to_tpm(self):
        return

    @staticmethod
    def is_valid(assignment):
        return assignment is not None and \
               hasattr(assignment, 'exon_gene_profile') and assignment.exon_gene_profile is not None  and \
               hasattr(assignment, 'intron_gene_profile') and assignment.intron_gene_profile is not None and \
               hasattr(assignment, 'gene_info') and assignment.gene_info is not None


class ExonCounter(ProfileFeatureCounter):
    def __init__(self, output_prefix, ignore_read_groups=False):
        ProfileFeatureCounter.__init__(self, output_prefix, ignore_read_groups)

    def add_read_info(self, read_assignment):
        if not ProfileFeatureCounter.is_valid(read_assignment):
            return
        group_id = AbstractReadGrouper.default_group_id if self.ignore_read_groups else read_assignment.read_group
        self.add_read_info_from_profile(read_assignment.exon_gene_profile,
                                        read_assignment.gene_info.exon_property_map, group_id)


class IntronCounter(ProfileFeatureCounter):
    def __init__(self, output_prefix, ignore_read_groups=False):
        ProfileFeatureCounter.__init__(self, output_prefix, ignore_read_groups)

    def add_read_info(self, read_assignment):
        if not ProfileFeatureCounter.is_valid(read_assignment):
            return
        group_id = AbstractReadGrouper.default_group_id if self.ignore_read_groups else read_assignment.read_group
        self.add_read_info_from_profile(read_assignment.intron_gene_profile,
                                        read_assignment.gene_info.intron_property_map, group_id)

