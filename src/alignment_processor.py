############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import pysam
from Bio import SeqIO

from src.long_read_assigner import *
from src.long_read_profiles import *
from src.read_groups import *
from src.polya_finder import *
from src.polya_verification import *

logger = logging.getLogger('IsoQuant')


class LongReadAlignmentProcessor:
    """ class for aggregating all assignment information

    Parameters
    ----------
    gene_info
    bams
    params
    printer
    counter
    """

    def __init__(self, gene_info, bams, params, chr_record=None, read_groupper=DefaultReadGrouper()):
        self.gene_info = gene_info
        self.bams = bams
        self.params = params
        self.chr_record = chr_record

        self.assigner = LongReadAssigner(self.gene_info, self.params)
        self.read_groupper = read_groupper
        self.profile_constructor = CombinedProfileConstructor(gene_info, params)
        self.polya_finder = PolyAFinder(self.params.polya_window, self.params.polya_fraction)
        self.polya_fixer = PolyAFixer(self.params)
        self.cage_finder = CagePeakFinder(params.cage, params.cage_shift)
        self.assignment_storage = []

    def process(self):
        self.assignment_storage = []
        self.gene_info.all_read_region_start = self.gene_info.start
        self.gene_info.all_read_region_end = self.gene_info.end

        for b in self.bams:
            self.process_single_file(b)

        if (self.params.sqanti_output or self.params.check_canonical) and self.chr_record:
            self.gene_info.all_read_region_start -= self.params.upstream_region_len
            self.gene_info.all_read_region_end += self.params.upstream_region_len
            self.gene_info.reference_region = \
                str(self.chr_record[self.gene_info.all_read_region_start - 1:self.gene_info.all_read_region_end + 1].seq)
            self.gene_info.canonical_sites = {}
        return self.assignment_storage

    def process_single_file(self, bamfile_in):
        #with pysam.AlignmentFile(bam, "rb") as bamfile_in:
        for alignment in bamfile_in.fetch(self.gene_info.chr_id, self.gene_info.start, self.gene_info.end):
            read_id = alignment.query_name

            if alignment.reference_id == -1:
                self.assignment_storage.append(ReadAssignment(read_id, None))
                continue
            if alignment.is_supplementary:
                continue
            if self.params.no_secondary and alignment.is_secondary:
                continue

            logger.debug("=== Processing read " + read_id + " ===")

            # concat indels
            concat_blocks = concat_gapless_blocks(sorted(alignment.get_blocks()), alignment.cigartuples)
            if not concat_blocks:
                logger.warning("Read %s has no aligned exons" % read_id)
                continue
            # correct coordinates to GTF style (inclusive intervals)
            sorted_blocks = correct_bam_coords(concat_blocks)
            logger.debug("Read exons: " + str(sorted_blocks))
            if self.params.reference and (self.params.sqanti_output or self.params.check_canonical):
                if sorted_blocks[0][0] < self.gene_info.all_read_region_start:
                    self.gene_info.all_read_region_start = sorted_blocks[0][0]
                if sorted_blocks[-1][1] > self.gene_info.all_read_region_end:
                    self.gene_info.all_read_region_end = sorted_blocks[-1][1]

            polya_info = self.polya_finder.detect_polya(alignment)
            sorted_blocks, polya_info, exon_changed = self.polya_fixer.correct_read_info(sorted_blocks, polya_info)
            cage_hits = [] if self.params.cage is None else self.cage_finder.find_cage_peak(alignment)

            combined_profile = self.profile_constructor.construct_profiles(sorted_blocks, polya_info, cage_hits)
            read_assignment = self.assigner.assign_to_isoform(read_id, combined_profile)

            if exon_changed:
                read_assignment.add_match_attribute(MatchEvent(MatchEventSubtype.aligned_polya_tail))
            read_assignment.polyA_found = (polya_info.external_polya_pos != -1 or
                                           polya_info.external_polyt_pos != -1 or
                                           polya_info.internal_polya_pos != -1 or
                                           polya_info.internal_polyt_pos != -1)
            read_assignment.polya_info = polya_info
            read_assignment.cage_found = len(cage_hits) > 0
            read_assignment.exons = sorted_blocks
            read_assignment.read_group = self.read_groupper.get_group_id(alignment)
            read_assignment.mapped_strand = "-" if alignment.is_reverse else "+"
            read_assignment.multimapper = alignment.is_secondary

            if self.params.count_exons:
                read_assignment.exon_gene_profile = combined_profile.read_exon_profile.gene_profile
                read_assignment.intron_gene_profile = combined_profile.read_intron_profile.gene_profile

            if self.params.sqanti_output:
                indel_count, junctions_with_indels = self.count_indel_stats(alignment)
                read_assignment.set_additional_info("indel_count", indel_count)
                read_assignment.set_additional_info("junctions_with_indels", junctions_with_indels)
                read_assignment.introns_match = \
                    all(e == 1 for e in combined_profile.read_intron_profile.read_profile)

            self.assignment_storage.append(read_assignment)
            logger.debug("=== Finished read " + read_id + " ===")

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
