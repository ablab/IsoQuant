############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
from queue import PriorityQueue

from .common import interval_len, junctions_from_blocks, overlaps
from .gene_info import GeneInfo, StrandDetector
from .long_read_assigner import LongReadAssigner
from .long_read_profiles import CombinedProfileConstructor
from .isoform_assignment import (
    IsoformMatch,
    MatchClassification,
    MatchEvent,
    MatchEventSubtype,
    ReadAssignment,
    ReadAssignmentType,
)
from .read_groups import DefaultReadGrouper
from .polya_finder import PolyAFinder
from .polya_verification import PolyAFixer
from .exon_corrector import ExonCorrector
from .alignment_info import AlignmentInfo
from .illumina_exon_corrector import IlluminaExonCorrector, VoidExonCorrector
from .stats import EnumStats
from enum import Enum, unique

logger = logging.getLogger('IsoQuant')


@unique
class AlignmentType(Enum):
    primary = 1
    secondary = 2
    supplementary = 20
    unaligned = 255


def make_alignment_tuple(bam_index, alignment):
    return alignment.reference_start, alignment.reference_end, bam_index, alignment


class BAMOnlineMerger:
    # single interator for several bam files
    def __init__(self, bam_pairs, chr_id, start, end, multiple_iterators=False):
        self.bam_pairs = bam_pairs
        self.multiple_iterators = multiple_iterators
        self._set(chr_id, start, end)

    def _set(self, chr_id, start, end):
        self.chr_id = chr_id
        self.start = start
        self.end = end
        # fetch uses 0-based semi-closed interval
        self.alignment_iterators = [bp[0].fetch(self.chr_id, self.start, self.end + 1,
                                                multiple_iterators=self.multiple_iterators) for bp in self.bam_pairs]
        self.current_elements = PriorityQueue(len(self.alignment_iterators))
        for i, it in enumerate(self.alignment_iterators):
            try:
                self.current_elements.put_nowait(make_alignment_tuple(i, next(it)))
            except StopIteration:
                pass

    def get(self):
        while not self.current_elements.empty():
            alignment_tuple = self.current_elements.get_nowait()
            i = alignment_tuple[2]
            try:
                self.current_elements.put_nowait(make_alignment_tuple(i, next(self.alignment_iterators[i])))
            except StopIteration:
                pass
            yield i, alignment_tuple[3]

    def reset(self):
        for bp in self.bam_pairs:
            bp[0].reset()

    def reset_region(self, chr_id, start, end):
        self._set(chr_id, start, end)


class AbstractAlignmentStorage:
    COVERAGE_BIN = 256

    def __init__(self):
        self.coverage_dict = defaultdict(int)
        self.region = None

    def reset(self):
        self.coverage_dict = defaultdict(int)
        self.region = None

    def add_alignment(self, bam_index, alignment):
        bin_start = alignment.reference_start // AbstractAlignmentStorage.COVERAGE_BIN
        # closed interval, 0-based
        alignment_end = alignment.reference_end - 1
        bin_end = alignment_end // AbstractAlignmentStorage.COVERAGE_BIN
        for i in range(bin_start, bin_end + 1):
            self.coverage_dict[i] += 1
        if not self.region:
            self.region = (alignment.reference_start, alignment_end)
        else:
            self.region = (min(self.region[0], alignment.reference_start), max(self.region[1], alignment_end))

    def get_alignments(self, region=None):
        raise NotImplementedError()

    def get_read_count(self):
        raise NotImplementedError()

    def alignment_is_not_adjacent(self, alignment):
        return self.region is not None and \
            not overlaps(self.region, (alignment.reference_start, alignment.reference_end - 1))


class BAMAlignmentStorage(AbstractAlignmentStorage):
    def __init__(self, bam_merger):
        AbstractAlignmentStorage.__init__(self)
        self.bam_merger = bam_merger
        self.counter = 0

    def reset(self):
        AbstractAlignmentStorage.reset(self)
        self.counter = 0

    def add_alignment(self, bam_index, alignment):
        AbstractAlignmentStorage.add_alignment(self, bam_index, alignment)
        self.counter += 1

    def get_alignments(self, region=None):
        if not region:
            region = self.region
        return BAMOnlineMerger(self.bam_merger.bam_pairs, self.bam_merger.chr_id, region[0], region[1],
                               multiple_iterators=True).get()

    def get_read_count(self):
        return self.counter


class InMemoryAlignmentStorage(AbstractAlignmentStorage):
    def __init__(self):
        AbstractAlignmentStorage.__init__(self)
        self.alignment_start_index = {}
        self.alignment_end_index = {}
        self.counter = 0
        self.alignment_storage = []
        self.index_filled = False

    def reset(self):
        AbstractAlignmentStorage.reset(self)
        self.alignment_start_index = {}
        self.alignment_end_index = {}
        self.counter = 0
        self.alignment_storage = []
        self.index_filled = False

    def add_alignment(self, bam_index, alignment):
        AbstractAlignmentStorage.add_alignment(self, bam_index, alignment)
        bin_start_position = alignment.reference_start // self.COVERAGE_BIN
        if bin_start_position not in self.alignment_start_index:
            self.alignment_start_index[bin_start_position] = self.counter
        bin_end_position = (alignment.reference_end - 1) // self.COVERAGE_BIN
        if bin_end_position not in self.alignment_end_index:
            self.alignment_end_index[bin_end_position] = self.counter
        self.alignment_storage.append((bam_index, alignment))
        self.counter += 1
        self.index_filled = False

    def fill_index(self):
        if self.index_filled:
            return
        current_index = len(self.alignment_storage)
        current_bin_region_start = self.region[0] // AbstractAlignmentStorage.COVERAGE_BIN
        current_bin_region_end = self.region[1] // AbstractAlignmentStorage.COVERAGE_BIN
        for pos in range(current_bin_region_end + 1, current_bin_region_start - 1, -1):
            if pos not in self.alignment_start_index:
                self.alignment_start_index[pos] = current_index
            else:
                current_index = self.alignment_start_index[pos]
        current_index = len(self.alignment_storage)
        for pos in range(current_bin_region_end + 1, current_bin_region_start - 1, -1):
            if pos not in self.alignment_end_index or self.alignment_end_index[pos] > current_index:
                self.alignment_end_index[pos] = current_index
            else:
                current_index = self.alignment_end_index[pos]
        self.index_filled = True

    def get_alignments(self, region=None):
        if not region or region == self.region:
            for a in self.alignment_storage:
                yield a
            return

        self.fill_index()
        # first alignment among sorted that has its end inside the start_bin, e.g. close to region[0]
        start_bin = region[0] // self.COVERAGE_BIN
        start_index = self.alignment_end_index[start_bin]
        # first alignment that has its start after region[1]
        end_bin = region[1] // self.COVERAGE_BIN
        end_index = self.alignment_start_index[end_bin]

        for i in range(start_index, end_index):
            bam_index, alignment = self.alignment_storage[i]
            if overlaps(region, (alignment.reference_start, alignment.reference_end - 1)):
                yield bam_index, alignment

    def get_read_count(self):
        return len(self.alignment_storage)


class AlignmentCollector:
    """ class for aggregating all alignmnet information

    Parameters
    ----------
    bams : pair bam object, filename
    params
    printer
    counter
    """

    MAX_REGION_LEN = 32768
    MIN_READS_TO_SPLIT = 1024
    ABS_COV_VALLEY = 1
    REL_COV_VALLEY = 0.01

    def __init__(self, chr_id, bam_pairs, params, illumina_bam, genedb=None, chr_record=None, read_groupper=DefaultReadGrouper()):
        self.chr_id = chr_id
        self.bam_pairs = bam_pairs
        self.params = params
        self.genedb = genedb
        self.chr_record = chr_record
        self.illumina_bam = illumina_bam

        self.bam_merger = BAMOnlineMerger(self.bam_pairs, self.chr_id, 0,
                                          self.bam_pairs[0][0].get_reference_length(self.chr_id),
                                          multiple_iterators=not self.params.high_memory)
        self.strand_detector = StrandDetector(self.chr_record)
        self.read_groupper = read_groupper
        self.polya_finder = PolyAFinder(self.params.polya_window, self.params.polya_fraction)
        self.polya_fixer = PolyAFixer(self.params)
        # self.cage_finder = CagePeakFinder(params.cage, params.cage_shift)
        self.alignment_stat_counter = EnumStats()

    def process(self):
        alignment_storage = BAMAlignmentStorage(self.bam_merger) if not self.params.high_memory else InMemoryAlignmentStorage()
        for bam_index, alignment in self.bam_merger.get():
            if alignment.is_secondary:
                self.alignment_stat_counter.add(AlignmentType.secondary)
            elif alignment.is_supplementary:
                self.alignment_stat_counter.add(AlignmentType.supplementary)
            elif alignment.reference_id != -1:
                self.alignment_stat_counter.add(AlignmentType.primary)

            if alignment_storage.alignment_is_not_adjacent(alignment):
                for res in self.forward_alignments(alignment_storage):
                    yield res
                alignment_storage.reset()
            alignment_storage.add_alignment(bam_index, alignment)

        if alignment_storage.region:
            for res in self.forward_alignments(alignment_storage):
                yield res

    def forward_alignments(self, alignment_storage):
        current_region = alignment_storage.region
        logger.debug("Splitting " + str(current_region))
        split_regions = AlignmentCollector.split_coverage_regions(current_region, alignment_storage)

        if len(split_regions) == 1:
            yield self.process_alignments_in_region(current_region, alignment_storage.get_alignments())
        else:
            for new_region in split_regions:
                alignments = alignment_storage.get_alignments(new_region)
                yield self.process_alignments_in_region(new_region, alignments)

    def process_alignments_in_region(self, current_region, alignment_storage):
        logger.debug("Processing region %s" % str(current_region))
        gene_info = self.get_gene_info_for_region(current_region)
        if gene_info.empty():
            assignment_storage = self.process_intergenic(alignment_storage, current_region)
        else:
            assignment_storage = self.process_genic(alignment_storage, gene_info, current_region)

        return gene_info, assignment_storage

    def process_intergenic(self, alignment_storage, region):
        assignment_storage = []
        if self.illumina_bam is not None:
            corrector = IlluminaExonCorrector(self.chr_id, region[0], region[1], self.illumina_bam)
        else:
            corrector = VoidExonCorrector()

        for bam_index, alignment in alignment_storage:
            if alignment.reference_id == -1 or alignment.is_supplementary or \
                    (self.params.no_secondary and alignment.is_secondary):
                continue

            if self.params.min_mapq and alignment.mapping_quality < self.params.min_mapq:
                continue

            read_id = alignment.query_name
            alignment_info = AlignmentInfo(alignment)

            if not alignment_info.read_exons:
                logger.warning("Read %s has no aligned exons" % read_id)
                continue

            #if len(alignment_info.read_exons) > 2 and not alignment.is_secondary and \
            #        alignment.mapping_quality < self.params.multi_intron_mapping_quality_cutoff:
            #    continue
            if len(alignment_info.read_exons) <= 2 and \
                    (alignment.is_secondary or alignment.mapping_quality < self.params.simple_alignments_mapq_cutoff):
                continue

            alignment_info.add_polya_info(self.polya_finder, self.polya_fixer)
            # if self.params.cage:
            #    alignment_info.add_cage_info(self.cage_finder)

            read_assignment = ReadAssignment(read_id, ReadAssignmentType.intergenic,
                                             IsoformMatch(MatchClassification.intergenic))

            if alignment_info.exons_changed:
                read_assignment.add_match_attribute(MatchEvent(MatchEventSubtype.aligned_polya_tail))
            read_assignment.polyA_found = (alignment_info.polya_info.external_polya_pos != -1 or
                                           alignment_info.polya_info.external_polyt_pos != -1 or
                                           alignment_info.polya_info.internal_polya_pos != -1 or
                                           alignment_info.polya_info.internal_polyt_pos != -1)
            read_assignment.polya_info = alignment_info.polya_info
            read_assignment.cage_found = len(alignment_info.cage_hits) > 0
            read_assignment.genomic_region = region
            read_assignment.exons = alignment_info.read_exons
            read_assignment.corrected_exons = corrector.correct_read(alignment_info)
            read_assignment.corrected_introns = junctions_from_blocks(read_assignment.corrected_exons)

            read_assignment.read_group = self.read_groupper.get_group_id(alignment, self.bam_merger.bam_pairs[bam_index][1])
            read_assignment.mapped_strand = "-" if alignment.is_reverse else "+"
            read_assignment.strand = self.get_assignment_strand(read_assignment)
            read_assignment.chr_id = self.chr_id
            read_assignment.multimapper = alignment.is_secondary
            read_assignment.mapping_quality = alignment.mapping_quality
            AlignmentCollector.import_bam_tags(alignment, read_assignment, self.params.bam_tags)
            assignment_storage.append(read_assignment)
        return assignment_storage

    def process_genic(self, alignment_storage, gene_info, region):
        assigner = LongReadAssigner(gene_info, self.params)
        profile_constructor = CombinedProfileConstructor(gene_info, self.params)
        exon_corrector = ExonCorrector(gene_info, self.params, self.chr_record)
        assignment_storage = []

        for bam_index, alignment in alignment_storage:
            if alignment.reference_id == -1 or alignment.is_supplementary or \
                    (self.params.no_secondary and alignment.is_secondary):
                continue

            if self.params.min_mapq and alignment.mapping_quality < self.params.min_mapq:
                continue

            read_id = alignment.query_name
            logger.debug("=== Processing read " + read_id + " ===")
            alignment_info = AlignmentInfo(alignment)

            if not alignment_info.read_exons:
                logger.warning("Read %s has no aligned exons" % read_id)
                continue

            alignment_info.add_polya_info(self.polya_finder, self.polya_fixer)
            # if self.params.cage:
            #    alignment_info.add_cage_info(self.cage_finder)
            alignment_info.construct_profiles(profile_constructor)
            read_assignment = assigner.assign_to_isoform(read_id, alignment_info.combined_profile)

            if (not read_assignment.assignment_type in [ReadAssignmentType.unique,
                                                        ReadAssignmentType.unique_minor_difference,
                                                        ReadAssignmentType.ambiguous])\
                    and alignment.mapping_quality < self.params.inconsistent_mapq_cutoff:
                continue

            if alignment_info.exons_changed:
                read_assignment.add_match_attribute(MatchEvent(MatchEventSubtype.aligned_polya_tail))
            read_assignment.polyA_found = (alignment_info.polya_info.external_polya_pos != -1 or
                                           alignment_info.polya_info.external_polyt_pos != -1 or
                                           alignment_info.polya_info.internal_polya_pos != -1 or
                                           alignment_info.polya_info.internal_polyt_pos != -1)
            read_assignment.polya_info = alignment_info.polya_info
            read_assignment.cage_found = len(alignment_info.cage_hits) > 0
            read_assignment.genomic_region = region
            read_assignment.exons = alignment_info.read_exons
            read_assignment.corrected_exons = exon_corrector.correct_assigned_read(alignment_info,
                                                                                   read_assignment)
            read_assignment.corrected_introns = junctions_from_blocks(read_assignment.corrected_exons)

            read_assignment.read_group = self.read_groupper.get_group_id(alignment, self.bam_merger.bam_pairs[bam_index][1])
            read_assignment.mapped_strand = "-" if alignment.is_reverse else "+"
            read_assignment.strand = self.get_assignment_strand(read_assignment)
            AlignmentCollector.check_antisense(read_assignment)
            AlignmentCollector.import_bam_tags(alignment, read_assignment, self.params.bam_tags)

            read_assignment.chr_id = gene_info.chr_id
            read_assignment.multimapper = alignment.is_secondary
            read_assignment.mapping_quality = alignment.mapping_quality

            if self.params.count_exons:
                read_assignment.exon_gene_profile = alignment_info.combined_profile.read_exon_profile.gene_profile
                read_assignment.intron_gene_profile = alignment_info.combined_profile.read_intron_profile.gene_profile

            assignment_storage.append(read_assignment)
            logger.debug("=== Finished read " + read_id + " ===")
        return assignment_storage

    def get_assignment_strand(self, read_assignment):
        if read_assignment.isoform_matches and read_assignment.assignment_type in \
                [ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference]:
            return read_assignment.isoform_matches[0].transcript_strand

        has_polya = read_assignment.polya_info.external_polya_pos != -1 or read_assignment.polya_info.internal_polya_pos != -1
        has_polyt = read_assignment.polya_info.external_polyt_pos != -1 or read_assignment.polya_info.internal_polyt_pos != -1
        if len(read_assignment.exons) == 1:
            if has_polya and not has_polyt:
                return '+'
            elif has_polyt and not has_polya:
                return '-'
            return '.'
        return self.strand_detector.get_strand(read_assignment.corrected_introns, has_polya, has_polyt)

    def count_indel_stats(self, alignment):
        cigar_event_count = len(alignment.cigartuples)
        indel_events = [1, 2]
        indel_count = 0
        intron_cigar_positions = []
        for i in range(cigar_event_count):
            cigar = alignment.cigartuples[i]
            if cigar[0] in indel_events:
                indel_count += 1
            elif cigar[0] == 3:
                intron_cigar_positions.append(i)

        junctions_with_indels = 0
        for i in intron_cigar_positions:
            # indel right near intron
            if (i > 0 and alignment.cigartuples[i - 1][0] in indel_events) or \
                    (i < cigar_event_count - 1 and alignment.cigartuples[i + 1][0] in indel_events):
                junctions_with_indels += 1

            # indel separated by at most 'indel_near_splice_site_dist' matches from intron
            if (i > 1 and alignment.cigartuples[i - 2][0] in indel_events and
                alignment.cigartuples[i - 1][0] in [0, 7, 8] and
                alignment.cigartuples[i - 1][1] <= self.params.indel_near_splice_site_dist) or \
                    (i < cigar_event_count - 2 and alignment.cigartuples[i + 2][0] in indel_events and
                     alignment.cigartuples[i + 1][0]  in [0, 7, 8] and
                     alignment.cigartuples[i + 1][1] <= self.params.indel_near_splice_site_dist):
                junctions_with_indels += 1

        return indel_count, junctions_with_indels

    def get_gene_info_for_region(self, current_region):
        if not self.genedb:
            return GeneInfo.from_region(self.chr_id, current_region[0], current_region[1],
                                        self.params.delta, self.chr_record)

        gene_list = list(self.genedb.region(seqid=self.chr_id, start=current_region[0],
                                            end=current_region[1], featuretype="gene"))
        if not gene_list:
            return GeneInfo.from_region(self.chr_id, current_region[0], current_region[1],
                                        self.params.delta, self.chr_record)
        gene_list = sorted(gene_list, key=lambda x: x.start)
        gene_info = GeneInfo(gene_list, self.genedb, self.params.delta)
        if self.params.needs_reference:
            gene_info.set_reference_sequence(current_region[0], current_region[1], self.chr_record)
        return gene_info

    @staticmethod
    def split_coverage_regions(genomic_region, alignment_storage):
        if interval_len(genomic_region) < AlignmentCollector.MAX_REGION_LEN and \
                alignment_storage.get_read_count() < AlignmentCollector.MIN_READS_TO_SPLIT:
            return [genomic_region]

        split_regions = []
        coverage_dict = alignment_storage.coverage_dict
        coverage_positions = sorted(coverage_dict.keys())
        current_start = coverage_positions[0]
        min_bins = int(AlignmentCollector.MAX_REGION_LEN / AbstractAlignmentStorage.COVERAGE_BIN)
        pos = current_start + 1
        max_cov = coverage_dict[current_start]
        while pos <= coverage_positions[-1]:
            while (pos <= coverage_positions[-1] and pos - current_start < min_bins) or \
                    coverage_dict[pos] > max(AlignmentCollector.ABS_COV_VALLEY, max_cov * AlignmentCollector.REL_COV_VALLEY):
                max_cov = max(max_cov, coverage_dict[pos])
                pos += 1
            split_regions.append((max(current_start * AbstractAlignmentStorage.COVERAGE_BIN + 1, genomic_region[0]),
                                  min(pos * AbstractAlignmentStorage.COVERAGE_BIN, genomic_region[1])))
            current_start = pos
            max_cov = coverage_dict[current_start]
            pos = min(current_start + 1, coverage_positions[-1] + 1)

        return split_regions

    @staticmethod
    def check_antisense(read_assignment):
        for match in read_assignment.isoform_matches:
            if match.assigned_transcript is None:
                continue

            if read_assignment.strand == match.transcript_strand or read_assignment.strand == '.':
                continue
            if match.match_classification not in \
                    [MatchClassification.novel_in_catalog, MatchClassification.novel_not_in_catalog]:
                match.match_classification = MatchClassification.antisense
            match.match_subclassifications.append(MatchEvent(MatchEventSubtype.antisense))

    @staticmethod
    def import_bam_tags(bam_record, read_assignment, tag_list):
        for tag in tag_list:
            try:
                tag_value = bam_record.get_tag(tag)
                read_assignment.set_additional_attribute(tag, str(tag_value))
            except KeyError:
                pass
