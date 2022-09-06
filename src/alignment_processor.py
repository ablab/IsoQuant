############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import pysam
from queue import PriorityQueue, Empty
from Bio import SeqIO

from src.long_read_assigner import *
from src.long_read_profiles import *
from src.read_groups import *
from src.polya_finder import *
from src.polya_verification import *
from src.exon_corrector import *
from src.alignment_info import *

logger = logging.getLogger('IsoQuant')


def make_alignment_tuple(bam_index, alignment):
    return alignment.reference_start, alignment.reference_end, bam_index, alignment


class BAMOnlineMerger:
    # single interator for several bam files
    def __init__(self, bam_pairs, chr_id, start, end):
        self.bam_pairs = bam_pairs
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.alignment_iterators = [bp[0].fetch(self.chr_id, self.start, self.end) for bp in self.bam_pairs]
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


class IntergenicAlignmentCollector:
    """ class for aggregating all alignmnet information

    Parameters
    ----------
    bams : pair bam object, filename
    params
    printer
    counter
    """

    def __init__(self, chr_id, bam_pairs, params, genedb=None, chr_record=None, read_groupper=DefaultReadGrouper()):
        self.chr_id = chr_id
        self.start = 1
        self.bam_pairs = bam_pairs
        self.params = params
        self.genedb = genedb
        self.chr_record = chr_record

        self.bam_merger = BAMOnlineMerger(self.bam_pairs, self.chr_id, self.start,
                                          self.bam_pairs[0][0].get_reference_length(self.chr_id))
        self.strand_detector = StrandDetector(self.chr_record)
        self.read_groupper = read_groupper
        self.polya_finder = PolyAFinder(self.params.polya_window, self.params.polya_fraction)
        self.polya_fixer = PolyAFixer(self.params)
        self.cage_finder = CagePeakFinder(params.cage, params.cage_shift)

        self.COVERAGE_BIN = 256
        self.MAX_REGION_LEN = 32768
        self.ABS_COV_VALLEY = 1
        self.REL_COV_VALLEY = 0.01

    def process(self):
        # coverage every COVERAGE_BIN bases
        coverage_dict = defaultdict(int)
        # index of the first alignment with the start >= of certain bin (key)
        alignment_index = {}
        counter = 0
        current_region = None
        # TODO: implement low-memory version (i.e. double fetch)
        alignment_storage = []

        for bam_index, alignment in self.bam_merger.get():
            if alignment.reference_id == -1 or alignment.is_supplementary or \
                    (self.params.no_secondary and alignment.is_secondary):
                continue

            if not current_region:
                current_region = (alignment.reference_start, alignment.reference_end)
            elif overlaps(current_region, (alignment.reference_start, alignment.reference_end)):
                current_region = (current_region[0], max(current_region[1], alignment.reference_end))
            else:
                for res in self.forward_alignments(current_region, alignment_storage, coverage_dict, alignment_index):
                    yield res
                counter = 0
                alignment_storage = []
                alignment_index = {}
                coverage_dict = defaultdict(int)
                current_region = None

            for i in range(alignment.reference_start // self.COVERAGE_BIN, alignment.reference_end // self.COVERAGE_BIN + 1):
                coverage_dict[i] += 1

            bin_position = alignment.reference_start // self.COVERAGE_BIN
            if bin_position not in alignment_index:
                alignment_index[bin_position] = counter
            alignment_storage.append((bam_index, alignment))
            counter += 1

        if current_region:
            for res in self.forward_alignments(current_region, alignment_storage, coverage_dict, alignment_index):
                yield res

    def forward_alignments(self, current_region, alignment_storage, coverage_dict, alignment_index):
        if interval_len(current_region) < self.MAX_REGION_LEN:
            yield self.process_alignments_in_region(current_region, alignment_storage)
            return

        coverage_positions = sorted(coverage_dict.keys())
        # completing index
        current_index = len(alignment_storage)
        for pos in range(coverage_positions[-1] + 1, coverage_positions[0] - 1, -1):
            if pos not in alignment_index:
                alignment_index[pos] = current_index
            else:
                current_index = alignment_index[pos]

        current_start = coverage_positions[0]
        min_bins = int(self.MAX_REGION_LEN / self.COVERAGE_BIN)
        pos = current_start + 1
        max_cov = coverage_dict[current_start]

        while pos <= coverage_positions[-1]:
            while (pos <= coverage_positions[-1] and pos - current_start < min_bins) or \
                    coverage_dict[pos] > max(self.ABS_COV_VALLEY, max_cov * self.REL_COV_VALLEY):
                max_cov = max(max_cov, coverage_dict[pos])
                pos += 1
            new_region = (max(current_start * self.COVERAGE_BIN, current_region[0]),
                          min(pos * self.COVERAGE_BIN, current_region[1]))
            alignments = alignment_storage[alignment_index[current_start]:alignment_index[pos]]
            if alignments:
                yield self.process_alignments_in_region(new_region, alignments)
            current_start = pos
            max_cov = coverage_dict[current_start]
            pos = min(current_start + 1, coverage_positions[-1] + 1)


        if current_start < pos:
            new_region = (max(current_start * self.COVERAGE_BIN, current_region[0]),
                          min(pos * self.COVERAGE_BIN, current_region[1]))
            alignments = alignment_storage[alignment_index[current_start]:alignment_index[pos]]
            if alignments:
                yield self.process_alignments_in_region(new_region, alignments)

    def process_alignments_in_region(self, current_region, alignment_storage):
        logger.debug("Processing region %s with %d reads" % (str(current_region), len(alignment_storage)))
        gene_list = []
        if self.genedb:
            gene_list = list(self.genedb.region(seqid=self.chr_id, start=current_region[0],
                                                end=current_region[1], featuretype="gene"))

        if not gene_list:
            gene_info = GeneInfo.from_region(self.chr_id, current_region[0], current_region[1],
                                             self.params.delta, self.chr_record)
            assignment_storage = self.process_intergenic(alignment_storage)
        else:
            gene_list = sorted(gene_list, key=lambda x: x.start)
            gene_info = GeneInfo(gene_list, self.genedb, self.params.delta)
            if self.params.needs_reference:
                gene_info.set_reference_sequence(current_region[0], current_region[1], self.chr_record)
            assignment_storage = self.process_genic(alignment_storage, gene_info)

        return gene_info, assignment_storage

    def process_intergenic(self, alignment_storage):
        assignment_storage = []
        for bam_index, alignment in alignment_storage:
            read_id = alignment.query_name
            alignment_info = AlignmentInfo(alignment)

            if not alignment_info.read_exons:
                logger.warning("Read %s has no aligned exons" % read_id)
                continue

            #if len(alignment_info.read_exons) > 2 and not alignment.is_secondary and \
            #        alignment.mapping_quality < self.params.multi_intron_mapping_quality_cutoff:
            #    continue
            if len(alignment_info.read_exons) <= 2 and \
                    (alignment.is_secondary or alignment.mapping_quality < self.params.mono_mapping_quality_cutoff):
                continue

            alignment_info.add_polya_info(self.polya_finder, self.polya_fixer)
            if self.params.cage:
                alignment_info.add_cage_info(self.cage_finder)

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
            read_assignment.exons = alignment_info.read_exons
            read_assignment.corrected_exons = alignment_info.read_exons
            read_assignment.corrected_introns = junctions_from_blocks(read_assignment.corrected_exons)

            read_assignment.read_group = self.read_groupper.get_group_id(alignment, self.bam_merger.bam_pairs[bam_index][1])
            read_assignment.mapped_strand = "-" if alignment.is_reverse else "+"
            read_assignment.strand = self.get_assignment_strand(read_assignment)
            read_assignment.chr_id = self.chr_id
            read_assignment.multimapper = alignment.is_secondary
            read_assignment.mapping_quality = alignment.mapping_quality
            assignment_storage.append(read_assignment)
        return assignment_storage

    def process_genic(self, alignment_storage, gene_info):
        assigner = LongReadAssigner(gene_info, self.params)
        profile_constructor = CombinedProfileConstructor(gene_info, self.params)
        exon_corrector = ExonCorrector(gene_info, self.params, self.chr_record)
        assignment_storage = []

        for bam_index, alignment in alignment_storage:
            read_id = alignment.query_name
            logger.debug("=== Processing read " + read_id + " ===")
            alignment_info = AlignmentInfo(alignment)

            if not alignment_info.read_exons:
                logger.warning("Read %s has no aligned exons" % read_id)
                continue

            alignment_info.add_polya_info(self.polya_finder, self.polya_fixer)
            if self.params.cage:
                alignment_info.add_cage_info(self.cage_finder)
            alignment_info.construct_profiles(profile_constructor)
            read_assignment = assigner.assign_to_isoform(read_id, alignment_info.combined_profile)

            if (not read_assignment.assignment_type in [ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference])\
                    and not alignment.is_secondary and \
                    alignment.mapping_quality < self.params.multi_intron_mapping_quality_cutoff:
                continue

            if alignment_info.exons_changed:
                read_assignment.add_match_attribute(MatchEvent(MatchEventSubtype.aligned_polya_tail))
            read_assignment.polyA_found = (alignment_info.polya_info.external_polya_pos != -1 or
                                           alignment_info.polya_info.external_polyt_pos != -1 or
                                           alignment_info.polya_info.internal_polya_pos != -1 or
                                           alignment_info.polya_info.internal_polyt_pos != -1)
            read_assignment.polya_info = alignment_info.polya_info
            read_assignment.cage_found = len(alignment_info.cage_hits) > 0
            read_assignment.exons = alignment_info.read_exons
            read_assignment.corrected_exons = exon_corrector.correct_assigned_read(alignment_info,
                                                                                   read_assignment)
            read_assignment.corrected_introns = junctions_from_blocks(read_assignment.corrected_exons)

            read_assignment.read_group = self.read_groupper.get_group_id(alignment, self.bam_merger.bam_pairs[bam_index][1])
            read_assignment.mapped_strand = "-" if alignment.is_reverse else "+"
            read_assignment.strand = self.get_assignment_strand(read_assignment)
            read_assignment.chr_id = gene_info.chr_id
            read_assignment.multimapper = alignment.is_secondary
            read_assignment.mapping_quality = alignment.mapping_quality

            if self.params.count_exons:
                read_assignment.exon_gene_profile = alignment_info.combined_profile.read_exon_profile.gene_profile
                read_assignment.intron_gene_profile = alignment_info.combined_profile.read_intron_profile.gene_profile

            if self.params.sqanti_output:
                indel_count, junctions_with_indels = self.count_indel_stats(alignment)
                read_assignment.set_additional_info("indel_count", indel_count)
                read_assignment.set_additional_info("junctions_with_indels", junctions_with_indels)
                read_assignment.introns_match = \
                    all(e == 1 for e in alignment_info.combined_profile.read_intron_profile.read_profile)

            assignment_storage.append(read_assignment)
            logger.debug("=== Finished read " + read_id + " ===")
        return assignment_storage

    def get_assignment_strand(self, read_assignment):
        if read_assignment.isoform_matches and read_assignment.assignment_type in \
                [ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference]:
            return read_assignment.isoform_matches[0].transcript_strand

        if len(read_assignment.exons) == 1:
            return '.'
        return self.strand_detector.get_strand(read_assignment.corrected_introns)

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


