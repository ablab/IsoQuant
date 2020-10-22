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

        gene_region = (gene_info.start, gene_info.end)
        self.assigner = LongReadAssigner(self.gene_info, self.params)
        self.read_groupper = read_groupper
        self.intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.intron_profiles.features, gene_region,
                                                  comparator=partial(equal_ranges, delta=self.params.delta),
                                                  # TODO parameter max_exon_extension
                                                  absense_condition=partial(overlaps_at_least, delta=self.params.max_exon_extension),
                                                  delta=self.params.delta)
        # TODO check for non split exons which do overlap
        self.exon_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.exon_profiles.features, gene_region,
                                                  comparator=partial(equal_ranges, delta=self.params.delta),
                                                  delta=self.params.delta)
        # TODO think whether overlaps should be changed to contains to avoid terminal partially covered exons
        self.split_exon_profile_constructor = \
            NonOverlappingFeaturesProfileConstructor(self.gene_info.split_exon_profiles.features,
                                                     comparator=partial(overlaps_at_least, delta=self.params.delta))
        self.polya_finder = PolyAFinder()
        self.cage_finder = CagePeakFinder(params.cage, params.cage_shift)
        self.assignment_storage = []

    def process(self):
        self.assignment_storage = []
        self.gene_info.all_read_region_start = self.gene_info.start
        self.gene_info.all_read_region_end = self.gene_info.end

        for b in self.bams:
            self.process_single_file(b)

        if self.params.sqanti_output and self.chr_record:
            self.gene_info.all_read_region_start -= self.params.upstream_region_len
            self.gene_info.all_read_region_end += self.params.upstream_region_len

            self.gene_info.reference_region = \
                str(self.chr_record[self.gene_info.all_read_region_start - 1:self.gene_info.all_read_region_end + 1].seq)
            self.gene_info.canonical_sites = {}
        return self.assignment_storage

    def process_single_file(self, bam):
        with pysam.AlignmentFile(bam, "rb") as bamfile_in:
            #self.counter.add_unaligned(bamfile_in.unmapped)

            for alignment in bamfile_in.fetch(self.gene_info.chr_id, self.gene_info.start, self.gene_info.end):
                read_id = alignment.query_name

                if alignment.reference_id == -1:
                    self.assignment_storage.append(ReadAssignment(read_id, None))
                    continue
                if alignment.is_supplementary:
                    continue
                if not self.params.use_secondary and alignment.is_secondary:
                    continue

                logger.debug("=== Processing read " + read_id + " ===")

                concat_blocks = concat_gapless_blocks(sorted(alignment.get_blocks()), alignment.cigartuples)
                sorted_blocks = correct_bam_coords(concat_blocks)
                if self.params.has_polya:
                    polya_pos = self.polya_finder.find_polya_tail(alignment)
                    polyt_pos = self.polya_finder.find_polyt_head(alignment)
                else:
                    polya_pos = -1
                    polyt_pos = -1

                cage_hits = [] if self.params.cage is None else self.cage_finder.find_cage_peak(alignment)

                if self.params.reference and self.params.sqanti_output:
                    if sorted_blocks[0][0] < self.gene_info.all_read_region_start:
                        self.gene_info.all_read_region_start = sorted_blocks[0][0]
                    if sorted_blocks[-1][1] > self.gene_info.all_read_region_end:
                        self.gene_info.all_read_region_start = sorted_blocks[-1][1]

                intron_profile = self.intron_profile_constructor.construct_intron_profile(sorted_blocks)
                exon_profile = self.exon_profile_constructor.construct_exon_profile(sorted_blocks)
                split_exon_profile = self.split_exon_profile_constructor.construct_profile(sorted_blocks)
                combined_profile = CombinedReadProfiles(intron_profile, exon_profile, split_exon_profile,
                                                        polya_pos=polya_pos, polyt_pos=polyt_pos, cage_hits=cage_hits)

                read_assignment = self.assigner.assign_to_isoform(read_id, combined_profile)
                read_assignment.polyA_found = (polya_pos != -1 or polyt_pos != -1)
                read_assignment.cage_found = len(cage_hits) > 0
                read_assignment.combined_profile = combined_profile
                read_assignment.gene_info = self.gene_info
                read_assignment.read_group = self.read_groupper.get_group_id(alignment)
                read_assignment.mapped_strand = "-" if alignment.is_reverse else "+"

                if self.params.sqanti_output:
                    indel_count, junctions_with_indels = self.count_indel_stats(alignment)
                    read_assignment.set_additional_info("indel_count", indel_count)
                    read_assignment.set_additional_info("junctions_with_indels", junctions_with_indels)

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
            if (i > 1 and alignment.cigartuples[i - 2][0] in indel_events and alignment.cigartuples[i - 1][0] == 0 and
                alignment.cigartuples[i - 1][1] <= self.params.indel_near_splice_site_dist) or \
                    (i < cigar_event_count - 2 and alignment.cigartuples[i + 2][0] in indel_events and
                     alignment.cigartuples[i + 1][0] == 0 and
                     alignment.cigartuples[i + 1][1] <= self.params.indel_near_splice_site_dist):
                junctions_with_indels += 1

        return indel_count, junctions_with_indels
