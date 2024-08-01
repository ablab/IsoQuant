############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict, OrderedDict
from enum import Enum, unique

from .isoform_assignment import (
    ReadAssignmentType,
    nonintronic_events,
)
from .gene_info import FeatureInfo
from .read_groups import AbstractReadGrouper

logger = logging.getLogger('IsoQuant')


@unique
class CountingStrategy(Enum):
    unique_only = 1
    with_ambiguous = 2
    unique_splicing_consistent = 10
    unique_inconsistent = 11
    # all_splicing_consistent = 19
    all = 20

    def no_inconsistent(self):
        return self in [CountingStrategy.unique_only, CountingStrategy.with_ambiguous]

    def ambiguous(self):
        return self in [CountingStrategy.all, CountingStrategy.with_ambiguous]

    def inconsistent_minor(self):
        return self in [CountingStrategy.unique_splicing_consistent,
                        CountingStrategy.unique_inconsistent, CountingStrategy.all]

    def inconsistent(self):
        return self in [CountingStrategy.unique_inconsistent, CountingStrategy.all]


COUNTING_STRATEGIES = [CountingStrategy.unique_only.name,
                       CountingStrategy.with_ambiguous.name,
                       CountingStrategy.unique_splicing_consistent.name,
                       CountingStrategy.unique_inconsistent.name,
                       CountingStrategy.all.name]


@unique
class NormalizationMethod(Enum):
    simple = 1
    usable_reads = 2


class CountingStrategyFlags:
    def __init__(self, counting_strategy):
        self.use_ambiguous = counting_strategy.ambiguous()
        self.use_inconsistent_minor = counting_strategy.inconsistent_minor()
        self.use_inconsistent = counting_strategy.inconsistent()


# the difference from defaultdict is that it does not add 0 to the dict when get is called for inexistent key
class IncrementalDict:
    def __init__(self, default_type=float):
        self.data = {}
        self.default_type = default_type

    def inc(self, key, value=1):
        if key not in self.data:
            self.data[key] = self.default_type(value)
        else:
            self.data[key] += value

    def get(self, key):
        if key not in self.data:
            return self.default_type()
        else:
            return self.data[key]


class GeneAssignmentExtractor:
    @staticmethod
    def get_features(read_assignment):
        gene_set = set()
        for m in read_assignment.isoform_matches:
            if m.assigned_gene:
                gene_set.add(m.assigned_gene)
        return gene_set

    @staticmethod
    def get_assignment_type(read_assignment):
        return read_assignment.gene_assignment_type

    @staticmethod
    def confirms_feature(read_assignment):
        return read_assignment.gene_assignment_type.is_unique()


class TranscriptAssignmentExtractor:
    @staticmethod
    def get_features(read_assignment):
        transcript_set = set()
        for m in read_assignment.isoform_matches:
            if m.assigned_transcript:
                transcript_set.add(m.assigned_transcript)
        return transcript_set

    @staticmethod
    def get_assignment_type(read_assignment):
        return read_assignment.assignment_type

    @staticmethod
    def confirms_feature(read_assignment):
        transcript_id = read_assignment.isoform_matches[0].assigned_transcript
        return (read_assignment.assignment_type.is_unique() and
                (len(read_assignment.gene_info.all_isoforms_introns[transcript_id]) == 0 or
                 len(read_assignment.corrected_exons) > 1))


class ReadWeightCounter:
    def __init__(self, strategy_str):
        self.strategy = CountingStrategy[strategy_str]
        self.strategy_flags = CountingStrategyFlags(self.strategy)

    def process_inconsistent(self, assignment_type, feature_count):
        # use only for inconsistent assignments
        if assignment_type == ReadAssignmentType.inconsistent_ambiguous or feature_count > 1:
            if self.strategy_flags.use_ambiguous and self.strategy_flags.use_inconsistent:
                return 1.0 / feature_count
            else:
                return 0.0
        if self.strategy_flags.use_inconsistent:
            return 1.0
        if self.strategy_flags.use_inconsistent_minor and assignment_type == ReadAssignmentType.inconsistent_non_intronic:
            return 1.0
        return 0.0

    def process_ambiguous(self, feature_count):
        if feature_count == 0:
            return 0.0
        if feature_count == 1:
            return 1.0
        if self.strategy_flags.use_ambiguous:
            return 1.0 / float(feature_count)
        else:
            return 0.0


class AbstractCounter:
    def __init__(self, output_prefix, ignore_read_groups=False, output_zeroes=True):
        self.ignore_read_groups = ignore_read_groups
        self.output_counts_file_name = output_prefix + "_counts.tsv"
        self.output_file = self.output_counts_file_name
        open(self.output_file, "w").close()
        self.output_tpm_file_name = output_prefix + "_tpm.tsv"
        self.output_zeroes = output_zeroes
        self.output_stats_file_name = None

    def get_output_file_handler(self):
        return open(self.output_file, "a")

    def add_read_info(self, read_assignment):
        raise NotImplementedError()

    def add_read_info_raw(self, read_id, feature_ids, group_id=AbstractReadGrouper.default_group_id):
        raise NotImplementedError()

    def add_confirmed_features(self, features):
        raise NotImplementedError()

    def dump(self):
        raise NotImplementedError()

    def add_unassigned(self, n_reads=1):
        raise NotImplementedError()

    def add_unaligned(self, n_reads=1):
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

    def add_unassigned(self, n_reads=1):
        for p in self.counters:
            p.add_unassigned(n_reads)

    def add_unaligned(self, n_reads=1):
        for p in self.counters:
            p.add_unaligned(n_reads)


# count meta-features assigned to reads (genes or isoforms)
# get_feature_id --- function that returns feature id form IsoformMatch object
class AssignedFeatureCounter(AbstractCounter):
    def __init__(self, output_prefix, assignment_extractor, read_groups, read_counter,
                 all_features=None, output_zeroes=True):
        AbstractCounter.__init__(self, output_prefix, not read_groups, output_zeroes)
        self.assignment_extractor = assignment_extractor
        self.all_features = set(all_features) if all_features is not None else set()
        if not read_groups:
            self.group_numeric_ids = {AbstractReadGrouper.default_group_id: 0}
        else:
            self.group_numeric_ids = OrderedDict()
            if read_groups:
                for i, g in enumerate(read_groups):
                    self.group_numeric_ids[g] = i
        self.read_counter = read_counter

        self.ambiguous_reads = 0
        self.reads_for_tpm = 0
        self.not_assigned_reads = 0
        self.not_aligned_reads = 0
        # feature_id -> (group_id -> count)
        self.feature_counter = defaultdict(IncrementalDict)
        self.confirmed_features = set()
        self.output_stats_file_name = self.output_counts_file_name + ".stats"

    def add_read_info(self, read_assignment=None):
        if not read_assignment:
            self.not_aligned_reads += 1
            return
        elif read_assignment.assignment_type.is_unassigned() or not read_assignment.isoform_matches:
            self.not_assigned_reads += 1
            return
        elif read_assignment.isoform_matches and read_assignment.isoform_matches[0].assigned_transcript is None:
            self.not_assigned_reads += 1
            logger.warning("Assigned feature for read %s is None, will be skipped. "
                           "This message may be reported to the developers." % read_assignment.read_id)
            return

        feature_ids = self.assignment_extractor.get_features(read_assignment)
        assignment_type = self.assignment_extractor.get_assignment_type(read_assignment)
        group_id = AbstractReadGrouper.default_group_id if self.ignore_read_groups else read_assignment.read_group
        group_id = self.group_numeric_ids[group_id]

        self.reads_for_tpm += 1
        if assignment_type == ReadAssignmentType.ambiguous:
            count_value = self.read_counter.process_ambiguous(len(feature_ids))
            for feature_id in feature_ids:
                self.feature_counter[feature_id].inc(group_id, count_value)
                if count_value > 0:
                    self.all_features.add(feature_id)
            self.ambiguous_reads += 1

        elif assignment_type.is_inconsistent():
            count_value = self.read_counter.process_inconsistent(assignment_type, len(feature_ids))
            if count_value > 0:
                for feature_id in feature_ids:
                    self.feature_counter[feature_id].inc(group_id, count_value)
                    self.all_features.add(feature_id)

        elif assignment_type.is_unique():
            feature_id = list(feature_ids)[0]
            self.feature_counter[feature_id].inc(group_id, 1.0)
            self.all_features.add(feature_id)
            if self.assignment_extractor.confirms_feature(read_assignment):
                self.confirmed_features.add(feature_id)

    def add_read_info_raw(self, read_id, feature_ids, group_id=AbstractReadGrouper.default_group_id):
        if self.ignore_read_groups:
            group_id = AbstractReadGrouper.default_group_id
        group_id = self.group_numeric_ids[group_id]
        if not read_id:
            self.not_aligned_reads += 1
        elif not feature_ids:
            self.not_assigned_reads += 1
        elif len(feature_ids) > 1:
            self.ambiguous_reads += 1
            self.reads_for_tpm += 1
            for feature_id in feature_ids:
                count_value = self.read_counter.process_ambiguous(len(feature_ids))
                self.feature_counter[feature_id].inc(group_id, count_value)
                self.all_features.add(feature_id)
        else:
            self.feature_counter[feature_ids[0]].inc(group_id, 1.0)
            self.all_features.add(feature_ids[0])
            self.reads_for_tpm += 1

    def add_unassigned(self, n_reads=1):
        self.not_assigned_reads += n_reads
        self.reads_for_tpm += n_reads

    def add_unaligned(self, n_reads=1):
        self.not_aligned_reads += n_reads

    def format_header(self, all_groups, value_name="count"):
        if self.ignore_read_groups:
            return "#feature_id\t%s\n" % value_name
        else:
            return "#feature_id\t" + "\t".join(all_groups) + "\n"

    def add_confirmed_features(self, features):
        self.confirmed_features.update(features)

    def dump(self):
        all_features = sorted(filter(lambda x: x is not None, self.all_features))
        all_groups = sorted(self.group_numeric_ids.keys())

        for feature_id in all_features:
            if feature_id in self.confirmed_features:
                continue
            # ignoring unconfirmed features
            for g in self.feature_counter[feature_id].data.keys():
                self.feature_counter[feature_id].data[g] = 0.0

        with self.get_output_file_handler() as output_file:
            output_file.write(self.format_header(all_groups))
            default_group_id = list(self.group_numeric_ids.values())[0]
            for feature_id in all_features:
                if self.ignore_read_groups:
                    count = self.feature_counter[feature_id].get(default_group_id)
                    if not self.output_zeroes and count == 0:
                        continue
                    output_file.write("%s\t%.2f\n" % (feature_id, count))
                else:
                    row_count = 0
                    for group_id in self.feature_counter[feature_id].data.keys():
                        count = self.feature_counter[feature_id].data[group_id]
                        row_count += count
                    if not self.output_zeroes and row_count == 0:
                        continue
                    count_values = [self.feature_counter[feature_id].get(self.group_numeric_ids[group_id]) for group_id in all_groups]
                    output_file.write("%s\t%s\n" % (feature_id, "\t".join(["%.2f" % c for c in count_values])))

        if self.ignore_read_groups:
            with open(self.output_stats_file_name, "w") as f:
                f.write("__ambiguous\t%d\n" % self.ambiguous_reads)
                f.write("__no_feature\t%d\n" % self.not_assigned_reads)
                f.write("__not_aligned\t%d\n" % self.not_aligned_reads)
                f.write("__usable\t%d\n" % self.reads_for_tpm)

    def convert_counts_to_tpm(self, normalization_str=NormalizationMethod.simple.name):
        normalization = NormalizationMethod[normalization_str]
        total_counts = defaultdict(float)
        with open(self.output_counts_file_name) as f:
            for line in f:
                if line.startswith('_'): break
                if line.startswith('#'): continue
                fs = line.rstrip().split('\t')
                if self.ignore_read_groups:
                    total_counts[AbstractReadGrouper.default_group_id] += float(fs[1])
                else:
                    for j in range(len(fs) - 1):
                        total_counts[j] += float(fs[j + 1])

        scale_factors = {}
        unassigned_tpm = 0.0
        for group_id in total_counts.keys():
            if normalization == NormalizationMethod.usable_reads and self.ignore_read_groups and self.reads_for_tpm:
                total_reads = self.reads_for_tpm
                unassigned_tpm = 1000000.0 * (1 - total_counts[group_id] / total_reads)
            else:
                total_reads = total_counts[group_id] if total_counts[group_id] > 0 else 1.0
            scale_factors[group_id] = 1000000.0 / total_reads
            logger.debug("Scale factor for group %s = %.2f" % (group_id, scale_factors[group_id]))

        with open(self.output_tpm_file_name, "w") as outf:
            with open(self.output_counts_file_name) as f:
                for line in f:
                    if line.startswith('_'): break
                    if line.startswith('#'):
                        outf.write(line.replace("count", "TPM"))
                        continue
                    fs = line.rstrip().split('\t')
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
                if self.ignore_read_groups:
                    outf.write("%s\t%.6f\n" % ("__unassigned", unassigned_tpm))


def create_gene_counter(output_file_name, strategy, complete_feature_list=None,
                        read_groups=None, output_zeroes=True):
    read_weight_counter = ReadWeightCounter(strategy)
    return AssignedFeatureCounter(output_file_name, GeneAssignmentExtractor,
                                  read_groups, read_weight_counter,
                                  complete_feature_list, output_zeroes)


def create_transcript_counter(output_file_name, strategy, complete_feature_list=None,
                              read_groups=None, output_zeroes=True):
    read_weight_counter = ReadWeightCounter(strategy)
    return AssignedFeatureCounter(output_file_name, TranscriptAssignmentExtractor,
                                  read_groups, read_weight_counter,
                                  complete_feature_list, output_zeroes)


# count simple features inclusion/exclusion (exons / introns)
class ProfileFeatureCounter(AbstractCounter):
    def __init__(self, output_prefix, ignore_read_groups=False):
        AbstractCounter.__init__(self, output_prefix, ignore_read_groups)
        # feature_id -> (group_id -> count)
        self.inclusion_feature_counter = defaultdict(lambda: IncrementalDict(int))
        self.exclusion_feature_counter = defaultdict(lambda: IncrementalDict(int))
        self.feature_name_dict = OrderedDict()
        if ignore_read_groups:
            self.group_numeric_ids = {AbstractReadGrouper.default_group_id: 0}
        else:
            self.group_numeric_ids = OrderedDict()
        self.current_group_id = 1

    def add_read_info_from_profile(self, gene_feature_profile, feature_property_map,
                                   read_group = AbstractReadGrouper.default_group_id):
        if read_group not in self.group_numeric_ids:
            self.group_numeric_ids[read_group] = self.current_group_id
            self.current_group_id += 1

        group_id = self.group_numeric_ids[read_group]
        for i in range(len(gene_feature_profile)):
            if gene_feature_profile[i] == 1:
                feature_id = feature_property_map[i].id
                self.inclusion_feature_counter[feature_id].inc(group_id)
                if feature_id not in self.feature_name_dict:
                    self.feature_name_dict[feature_id] = feature_property_map[i].to_str()
            elif gene_feature_profile[i] == -1:
                feature_id = feature_property_map[i].id
                self.exclusion_feature_counter[feature_id].inc(group_id)
                if feature_id not in self.feature_name_dict:
                    self.feature_name_dict[feature_id] = feature_property_map[i].to_str()

    def dump(self):
        with open(self.output_counts_file_name, "w") as f:
            f.write(FeatureInfo.header() + "\tgroup_id\tinclude_counts\texclude_counts\n")

            all_groups = sorted(self.group_numeric_ids.keys())
            for feature_id in self.feature_name_dict.keys():
                for group_name in all_groups:
                    feature_name = self.feature_name_dict[feature_id]
                    group_id = self.group_numeric_ids[group_name]
                    incl_count = self.inclusion_feature_counter[feature_id].get(group_id)
                    excl_count = self.exclusion_feature_counter[feature_id].get(group_id)
                    if incl_count > 0 or excl_count > 0:
                        f.write("%s\t%s\t%d\t%d\n" % (feature_name, group_name, incl_count, excl_count))

    def convert_counts_to_tpm(self, normalization=NormalizationMethod.simple):
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

