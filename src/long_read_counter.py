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
)
from .gene_info import FeatureInfo
from .read_groups import AbstractReadGrouper
from .convert_grouped_counts import GROUP_COUNT_CUTOFF, convert_to_mtx, convert_to_matrix
from .file_naming import (
    counts_prefix, tpm_prefix, counts_file_name, tpm_file_name,
    counts_stats_file_name, counts_usable_file_name
)

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


class CountingStrategyFlags:
    def __init__(self, counting_strategy):
        self.use_ambiguous = counting_strategy.ambiguous()
        self.use_inconsistent_minor = counting_strategy.inconsistent_minor()
        self.use_inconsistent = counting_strategy.inconsistent()


@unique
class GroupedOutputFormat(Enum):
    none = 0
    matrix = 1
    mtx = 2
    default = 3


# the difference from defaultdict is that it does not add 0 to the dict when get is called for inexistent key
class IncrementalDict:
    def __init__(self, default_type=float):
        self.data = {}
        self.default_type = default_type

    def inc(self, key, value=1.0):
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
    def __init__(self, output_prefix, ignore_read_groups=False):
        self.ignore_read_groups = ignore_read_groups
        self.output_prefix = output_prefix
        self.output_counts_prefix = counts_prefix(output_prefix)
        self.output_tpm_prefix = tpm_prefix(output_prefix)
        self.output_counts_file_name = counts_file_name(self.output_counts_prefix, linear=not ignore_read_groups)
        self.output_tpm_file_name = tpm_file_name(self.output_tpm_prefix)
        self.output_file = self.output_counts_file_name
        open(self.output_file, "w").close()
        self.output_stats_file_name = None
        self.usable_file_name = None

    def get_output_file_handler(self):
        return open(self.output_file, "a")

    def add_read_info(self, read_assignment):
        raise NotImplementedError()

    def add_read_info_raw(self, read_id, feature_ids, group_ids):
        raise NotImplementedError()

    def add_confirmed_features(self, features):
        raise NotImplementedError()

    def dump(self):
        raise NotImplementedError()

    def finalize(self, args=None):
        raise NotImplementedError()

    def add_unassigned(self, read_assignment):
        raise NotImplementedError()

    def add_unaligned(self, n_reads=1):
        raise NotImplementedError()


class CompositeCounter:
    def __init__(self, counters: list = None):
        if counters is None:
            self.counters = []
        else:
            self.counters: list = counters

    def add_counters(self, counters):
        self.counters += counters

    def add_counter(self, counter):
        self.counters.append(counter)

    def add_read_info(self, read_assignment):
        for p in self.counters:
            p.add_read_info(read_assignment)

    def add_read_info_raw(self, read_id, feature_ids, group_ids):
        for p in self.counters:
            p.add_read_info_raw(read_id, feature_ids, group_ids)

    def add_confirmed_features(self, features):
        for p in self.counters:
            p.add_confirmed_features(features)

    def dump(self):
        for p in self.counters:
            p.dump()

    def finalize(self, args=None):
        for p in self.counters:
            p.finalize(args)

    def add_unassigned(self, read_assignment):
        for p in self.counters:
            p.add_unassigned(read_assignment)

    def add_unaligned(self, n_reads=1):
        for p in self.counters:
            p.add_unaligned(n_reads)


# count meta-features assigned to reads (genes or isoforms)
# get_feature_id --- function that returns feature id form IsoformMatch object
class AssignedFeatureCounter(AbstractCounter):
    def __init__(self, output_prefix, assignment_extractor, string_pools, read_counter,
                 all_features=None, group_index: int = 0):
        AbstractCounter.__init__(self, output_prefix, string_pools is None)
        self.assignment_extractor = assignment_extractor
        self.string_pools = string_pools
        self.all_features = set(all_features) if all_features is not None else set()
        self.group_index = group_index  # Index in read_group list to use
        self.read_counter = read_counter

        self.ambiguous_reads = 0
        self.reads_for_tpm = defaultdict(int)
        self.not_assigned_reads = 0
        self.not_aligned_reads = 0
        # feature_id -> (group_id -> count)  -- group_id is now an integer pool index
        self.feature_counter = defaultdict(IncrementalDict)
        self.confirmed_features = set()
        self.output_stats_file_name = counts_stats_file_name(self.output_counts_file_name)
        self.usable_file_name = counts_usable_file_name(self.output_counts_file_name)

    def add_read_info(self, read_assignment=None):
        # Use read_group_ids directly (integers) instead of string conversion
        if self.ignore_read_groups:
            group_id = 0  # Default group index
        elif not read_assignment.read_group_ids:
            group_id = 0
        else:
            group_id = read_assignment.read_group_ids[self.group_index]

        if not read_assignment:
            self.not_aligned_reads += 1
            return
        elif read_assignment.assignment_type.is_unassigned() or not read_assignment.isoform_matches:
            self.not_assigned_reads += 1
            self.reads_for_tpm[group_id] += 1
            return
        elif read_assignment.isoform_matches and read_assignment.isoform_matches[0].assigned_transcript is None:
            self.not_assigned_reads += 1
            self.reads_for_tpm[group_id] += 1
            logger.warning("Assigned feature for read %s is None, will be skipped. "
                           "This message may be reported to the developers." % read_assignment.read_id)
            return

        feature_ids = self.assignment_extractor.get_features(read_assignment)
        assignment_type = self.assignment_extractor.get_assignment_type(read_assignment)
        self.reads_for_tpm[group_id] += 1
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

    def add_read_info_raw(self, read_id, feature_ids, group_ids):
        # group_ids is now a list of integers (pool indices)
        if self.ignore_read_groups or group_ids is None:
            group_id = 0  # Default group index
        else:
            group_id = group_ids[self.group_index]
        if not read_id:
            self.not_aligned_reads += 1
        elif not feature_ids:
            self.not_assigned_reads += 1
            self.reads_for_tpm[group_id] += 1
        elif len(feature_ids) > 1:
            self.ambiguous_reads += 1
            self.reads_for_tpm[group_id] += 1
            for feature_id in feature_ids:
                count_value = self.read_counter.process_ambiguous(len(feature_ids))
                self.feature_counter[feature_id].inc(group_id, count_value)
                self.all_features.add(feature_id)
        else:
            self.feature_counter[feature_ids[0]].inc(group_id, 1.0)
            self.all_features.add(feature_ids[0])
            self.reads_for_tpm[group_id] += 1

    def add_unassigned(self, read_assignment):
        # Use read_group_ids directly (integers)
        if self.ignore_read_groups:
            group_id = 0
        elif not read_assignment.read_group_ids:
            group_id = 0
        else:
            group_id = read_assignment.read_group_ids[self.group_index]
        self.not_assigned_reads += 1
        self.reads_for_tpm[group_id] += 1

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

        for feature_id in all_features:
            if feature_id in self.confirmed_features:
                continue
            # ignoring unconfirmed features
            for g in self.feature_counter[feature_id].data.keys():
                self.feature_counter[feature_id].data[g] = 0.0

        if self.ignore_read_groups:
            self.dump_ungrouped(all_features)
        else:
            self.dump_grouped(all_features)
        self.dump_usable()

    def _get_ordered_groups(self):
        """Get ordered list of group names from string pool."""
        if self.string_pools is None:
            return [AbstractReadGrouper.default_group_id]
        pool = self.string_pools.get_read_group_pool(self.group_index)
        return [pool.get_str(i) for i in range(len(pool))]

    def _get_group_name(self, group_id: int) -> str:
        """Get group name from pool for a group ID."""
        if self.string_pools is None:
            return AbstractReadGrouper.default_group_id
        pool = self.string_pools.get_read_group_pool(self.group_index)
        return pool.get_str(group_id)

    def _get_num_groups(self) -> int:
        """Get number of groups from pool."""
        if self.string_pools is None:
            return 1
        pool = self.string_pools.get_read_group_pool(self.group_index)
        return len(pool)

    def dump_ungrouped(self, all_features):
        with self.get_output_file_handler() as output_file:
            default_group_id = 0  # Default group index
            output_file.write(self.format_header(self._get_ordered_groups()))
            for feature_id in all_features:
                count = self.feature_counter[feature_id].get(default_group_id)
                output_file.write("%s\t%.2f\n" % (feature_id, count))

            with open(self.output_stats_file_name, "w") as f:
                f.write("__ambiguous\t%d\n" % self.ambiguous_reads)
                f.write("__no_feature\t%d\n" % self.not_assigned_reads)
                f.write("__not_aligned\t%d\n" % self.not_aligned_reads)

    def dump_grouped(self, all_features):
        with self.get_output_file_handler() as output_file:
            output_file.write("#feature_id\tgroup_id\tcount\n")
            for feature_id in all_features:
                row_count = 0
                for group_id in self.feature_counter[feature_id].data.keys():
                    count = self.feature_counter[feature_id].data[group_id]
                    if count != 0:
                        group_name = self._get_group_name(group_id)
                        output_file.write("%s\t%s\t%.2f\n" % (feature_id, group_name, count))
                    row_count += count
            output_file.close()

    def dump_usable(self):
        with open(self.usable_file_name, "w") as f:
            # Iterate over all groups in pool order
            for group_id in range(self._get_num_groups()):
                group_name = self._get_group_name(group_id)
                f.write("%s\t%d\n" % (group_name, self.reads_for_tpm[group_id]))

    def load_usable(self, usable_file_name):
        if self.string_pools is None:
            # Ungrouped: just load into default group
            for l in open(usable_file_name, "r"):
                v = l.strip().split('\t')
                self.reads_for_tpm[0] += int(v[1])
        else:
            pool = self.string_pools.get_read_group_pool(self.group_index)
            for l in open(usable_file_name, "r"):
                v = l.strip().split('\t')
                group_name = v[0]
                if group_name in pool.str_to_int:
                    group_id = pool.str_to_int[group_name]
                    self.reads_for_tpm[group_id] += int(v[1])

    def finalize(self, args=None):
        if not args:
            normalization_method = NormalizationMethod.none
        elif args.normalization_method not in NormalizationMethod.__dict__:
            logger.warning("%s in not a valid normalization method, will use simple normalization" % args.normalization_method)
            normalization_method = NormalizationMethod.simple
        else:
            normalization_method = NormalizationMethod[args.normalization_method]

        if self.ignore_read_groups:
            if normalization_method == NormalizationMethod.none:
                return
            convert_ungrouped_to_tpm(self.output_file, self.output_tpm_file_name,
                                     normalization_method, self.reads_for_tpm[0])  # Default group index is 0
            return

        num_groups = self._get_num_groups()
        counts_format = {GroupedOutputFormat[f] for f in args.counts_format}
        if (GroupedOutputFormat.matrix in counts_format or
                (GroupedOutputFormat.default in counts_format and num_groups <= GROUP_COUNT_CUTOFF)):
            convert_to_matrix(self.output_file, self.output_counts_prefix)
            if normalization_method != NormalizationMethod.none:
                reads_for_tpm = None
                if normalization_method == NormalizationMethod.usable_reads:
                    reads_for_tpm = {self._get_group_name(group_id): self.reads_for_tpm[group_id]
                                     for group_id in range(num_groups)}
                convert_to_matrix(self.output_file, self.output_tpm_prefix,
                                  convert_to_tpm=True, usable_reads_per_group=reads_for_tpm)

        if (GroupedOutputFormat.mtx in counts_format or
                (GroupedOutputFormat.default in counts_format and num_groups > GROUP_COUNT_CUTOFF)):
            convert_to_mtx(self.output_file, self.output_counts_prefix)
            if normalization_method != NormalizationMethod.none:
                reads_for_tpm = None
                if normalization_method == NormalizationMethod.usable_reads:
                    reads_for_tpm = {self._get_group_name(group_id): self.reads_for_tpm[group_id]
                                     for group_id in range(num_groups)}
                convert_to_mtx(self.output_file, self.output_tpm_prefix,
                               convert_to_tpm=True, usable_reads_per_group=reads_for_tpm)

def create_gene_counter(output_file_name, strategy, complete_feature_list=None, string_pools=None, group_index: int = 0):
    read_weight_counter = ReadWeightCounter(strategy)
    return AssignedFeatureCounter(output_file_name, GeneAssignmentExtractor,
                                  string_pools, read_weight_counter, complete_feature_list, group_index)


def create_transcript_counter(output_file_name, strategy, complete_feature_list=None,
                              string_pools=None, group_index: int = 0):
    read_weight_counter = ReadWeightCounter(strategy)
    return AssignedFeatureCounter(output_file_name, TranscriptAssignmentExtractor,
                                  string_pools, read_weight_counter, complete_feature_list, group_index)


@unique
class NormalizationMethod(Enum):
    none = 0
    simple = 1
    usable_reads = 2


def convert_ungrouped_to_tpm(counts_file_name, output_tpm_file_name, normalization=NormalizationMethod.simple, reads_for_tpm=-1):
    total_counts = 0
    with open(counts_file_name, "r") as f:
        for line in f:
            if line.startswith('_'): continue
            if line.startswith('#'): continue
            fs = line.rstrip().split('\t')
            total_counts += float(fs[1])

    if total_counts == 0:
        logger.error("No counts were detected in %s, cannot convert to TPM" % counts_file_name)
        return

    unassigned_tpm = 0.0
    if normalization == NormalizationMethod.usable_reads and reads_for_tpm > 0:
        total_reads = reads_for_tpm
        unassigned_tpm = 1000000.0 * (1 - total_counts / total_reads)
    else:
        total_reads = total_counts
    scale_factor = 1000000.0 / total_reads

    with open(output_tpm_file_name, "w") as outf:
        with open(counts_file_name) as f:
            for line in f:
                if line.startswith('_'): break
                if line.startswith('#'):
                    outf.write(line.replace("count", "TPM"))
                    continue
                fs = line.rstrip().split('\t')
                feature_id, count = fs[0], float(fs[1])
                tpm = scale_factor * count
                outf.write("%s\t%.6f\n" % (feature_id, tpm))
            if unassigned_tpm > 0:
                outf.write("%s\t%.6f\n" % ("__unassigned", unassigned_tpm))


# count simple features inclusion/exclusion (exons / introns)
class ProfileFeatureCounter(AbstractCounter):
    def __init__(self, output_prefix, string_pools=None, group_index: int = 0):
        AbstractCounter.__init__(self, output_prefix, string_pools is None)
        self.string_pools = string_pools
        self.group_index = group_index  # Index in read_group list to use
        # feature_id -> (group_id -> count)  -- group_id is now an integer pool index
        self.inclusion_feature_counter = defaultdict(lambda: IncrementalDict(int))
        self.exclusion_feature_counter = defaultdict(lambda: IncrementalDict(int))
        self.feature_name_dict = OrderedDict()
        # Track which group IDs have been encountered (for dump)
        self.encountered_group_ids = set()

    def add_read_info_from_profile(self, gene_feature_profile, read_assigned_strand, feature_property_map,
                                   group_id: int = 0):
        """Process feature profile. group_id is an integer pool index."""
        self.encountered_group_ids.add(group_id)

        for i in range(len(gene_feature_profile)):
            if read_assigned_strand not in feature_property_map[i].strand:
                # skip features that do not match read's strand
                continue
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

    def _get_group_name(self, group_id: int) -> str:
        """Get group name from pool for a group ID."""
        if self.string_pools is None:
            return AbstractReadGrouper.default_group_id
        pool = self.string_pools.get_read_group_pool(self.group_index)
        return pool.get_str(group_id)

    def dump(self):
        with open(self.output_counts_file_name, "w") as f:
            f.write(FeatureInfo.header() + "\tgroup_id\tinclude_counts\texclude_counts\n")

            all_group_ids = sorted(self.encountered_group_ids)
            for feature_id in self.feature_name_dict.keys():
                for group_id in all_group_ids:
                    feature_name = self.feature_name_dict[feature_id]
                    group_name = self._get_group_name(group_id)
                    incl_count = self.inclusion_feature_counter[feature_id].get(group_id)
                    excl_count = self.exclusion_feature_counter[feature_id].get(group_id)
                    if incl_count > 0 or excl_count > 0:
                        f.write("%s\t%s\t%d\t%d\n" % (feature_name, group_name, incl_count, excl_count))

    def finalize(self, args=None):
        pass

    @staticmethod
    def is_valid(assignment):
        return assignment is not None and \
               hasattr(assignment, 'exon_gene_profile') and assignment.exon_gene_profile is not None and \
               hasattr(assignment, 'intron_gene_profile') and assignment.intron_gene_profile is not None and \
               hasattr(assignment, 'gene_info') and assignment.gene_info is not None

    @staticmethod
    def is_assigned_to_gene(assignment):
        return not assignment.gene_assignment_type.is_unassigned() and not assignment.gene_assignment_type.is_ambiguous()


class ExonCounter(ProfileFeatureCounter):
    def __init__(self, output_prefix, string_pools=None, group_index: int = 0):
        ProfileFeatureCounter.__init__(self, output_prefix, string_pools, group_index)

    def add_read_info(self, read_assignment):
        if not ProfileFeatureCounter.is_valid(read_assignment) or not ProfileFeatureCounter.is_assigned_to_gene(read_assignment):
            return
        # Use read_group_ids directly (integers)
        if self.ignore_read_groups:
            group_id = 0
        elif not read_assignment.read_group_ids:
            group_id = 0
        else:
            group_id = read_assignment.read_group_ids[self.group_index]
        self.add_read_info_from_profile(read_assignment.exon_gene_profile, read_assignment.strand,
                                        read_assignment.gene_info.exon_property_map, group_id)


class IntronCounter(ProfileFeatureCounter):
    def __init__(self, output_prefix, string_pools=None, group_index: int = 0):
        ProfileFeatureCounter.__init__(self, output_prefix, string_pools, group_index)

    def add_read_info(self, read_assignment):
        if not ProfileFeatureCounter.is_valid(read_assignment) or not ProfileFeatureCounter.is_assigned_to_gene(read_assignment):
            return
        # Use read_group_ids directly (integers)
        if self.ignore_read_groups:
            group_id = 0
        elif not read_assignment.read_group_ids:
            group_id = 0
        else:
            group_id = read_assignment.read_group_ids[self.group_index]
        self.add_read_info_from_profile(read_assignment.intron_gene_profile, read_assignment.strand,
                                        read_assignment.gene_info.intron_property_map, group_id)

