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
from src.exon_corrector import *
from src.alignment_info import *

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
        self.exon_corrector = ExonCorrector(self.gene_info, self.params)
        self.assignment_storage = []

    def process(self):
        self.assignment_storage = []
        self.gene_info.all_read_region_start = self.gene_info.start
        self.gene_info.all_read_region_end = self.gene_info.end

        for b in self.bams:
            self.process_single_file(b)

        if self.params.needs_reference:
            self.gene_info.all_read_region_start -= self.params.upstream_region_len
            self.gene_info.all_read_region_end += self.params.upstream_region_len
            self.gene_info.reference_region = \
                str(self.chr_record[self.gene_info.all_read_region_start - 1:self.gene_info.all_read_region_end + 1].seq)
            self.gene_info.canonical_sites = {}
        return self.assignment_storage

    def process_single_file(self, bamfile_in):
        processed_reads = set()
        for genic_region in self.gene_info.genic_regions:
            # FIXME: temporary solution - process gene outside
            for alignment in bamfile_in.fetch(self.gene_info.chr_id, genic_region[0] - 100, genic_region[1] + 100):
                read_id = alignment.query_name
                if alignment.reference_id == -1:
                    self.assignment_storage.append(ReadAssignment(read_id, None))
                    continue
                if alignment.is_supplementary:
                    continue
                if self.params.no_secondary and alignment.is_secondary:
                    continue

                logger.debug("=== Processing read " + read_id + " ===")
                alignment_info = AlignmentInfo(alignment)

                if not alignment_info.read_exons:
                    logger.warning("Read %s has no aligned exons" % read_id)
                    continue
                read_tuple = (read_id, alignment_info.read_start, alignment_info.read_end)
                if read_tuple in processed_reads:
                    continue
                processed_reads.add(read_tuple)

                logger.debug("Read exons: " + str(alignment_info.read_exons))
                if self.params.needs_reference:
                    if alignment_info.read_start < self.gene_info.all_read_region_start:
                        self.gene_info.all_read_region_start = alignment_info.read_start
                    if alignment_info.read_end > self.gene_info.all_read_region_end:
                        self.gene_info.all_read_region_end = alignment_info.read_end

                alignment_info.add_polya_info(self.polya_finder, self.polya_fixer)
                if self.params.cage:
                    alignment_info.add_cage_info(self.cage_finder)
                alignment_info.construct_profiles(self.profile_constructor)
                read_assignment = self.assigner.assign_to_isoform(read_id, alignment_info)

                if alignment_info.exons_changed:
                    read_assignment.add_match_attribute(MatchEvent(MatchEventSubtype.aligned_polya_tail))
                read_assignment.polyA_found = (alignment_info.polya_info.external_polya_pos != -1 or
                                               alignment_info.polya_info.external_polyt_pos != -1 or
                                               alignment_info.polya_info.internal_polya_pos != -1 or
                                               alignment_info.polya_info.internal_polyt_pos != -1)
                read_assignment.polya_info = alignment_info.polya_info
                read_assignment.cage_found = len(alignment_info.cage_hits) > 0
                read_assignment.exons = alignment_info.read_exons
                read_assignment.corrected_exons = self.exon_corrector.correct_assigned_read(alignment_info,
                                                                                            read_assignment)
                read_assignment.corrected_introns = junctions_from_blocks(read_assignment.corrected_exons)
                logger.debug("Original exons: %s" % str(alignment_info.read_exons))
                logger.debug("Corrected exons: %s" % str(read_assignment.corrected_exons ))
                read_assignment.read_group = self.read_groupper.get_group_id(alignment)
                read_assignment.mapped_strand = "-" if alignment.is_reverse else "+"
                read_assignment.strand = self.get_assignment_strand(read_assignment, alignment)
                read_assignment.chr_id = self.gene_info.chr_id
                read_assignment.multimapper = alignment.is_secondary

                if self.params.count_exons:
                    read_assignment.exon_gene_profile = alignment_info.combined_profile.read_exon_profile.gene_profile
                    read_assignment.intron_gene_profile = alignment_info.combined_profile.read_intron_profile.gene_profile

                if self.params.sqanti_output:
                    indel_count, junctions_with_indels = self.count_indel_stats(alignment)
                    read_assignment.set_additional_info("indel_count", indel_count)
                    read_assignment.set_additional_info("junctions_with_indels", junctions_with_indels)
                    read_assignment.introns_match = \
                        all(e == 1 for e in alignment_info.combined_profile.read_intron_profile.read_profile)

                self.assignment_storage.append(read_assignment)
                logger.debug("=== Finished read " + read_id + " ===")

    def get_assignment_strand(self, read_assignment, read_alignment):
        if read_assignment.isoform_matches:
            logger.debug(read_assignment.isoform_matches[0].transcript_strand)
            return read_assignment.isoform_matches[0].transcript_strand
        try:
            return read_alignment.get_tag('ts')
        except KeyError:
            return '.'

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
