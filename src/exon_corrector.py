############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import namedtuple

from src.common import *
from src.gene_info import *
from src.isoform_assignment import *
from src.long_read_profiles import *

logger = logging.getLogger('IsoQuant')


class ExonCorrector:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params
        self.delta = self.params.delta
        self.intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.intron_profiles.features,
                                                  (self.gene_info.start, self.gene_info.end),
                                                  comparator=partial(equal_ranges, delta=self.params.delta))

    def correct_assigned_read(self, alignment_info, read_assignment):
        if len(alignment_info.read_exons) == 1 or \
                read_assignment.assignment_type == ReadAssignmentType.noninformative or \
                not read_assignment.isoform_matches:
            return alignment_info.read_exons

        if read_assignment.assignment_type == ReadAssignmentType.unique:
            read_region, new_introns = self.correct_fuzzy_junctions(alignment_info, read_assignment)
        else:
            read_region, new_introns = self.correct_misalignments(alignment_info, read_assignment)

        if new_introns:
            corrected_exons = [(read_region[0], new_introns[0][0] - 1)]
            corrected_exons += junctions_from_blocks(new_introns)
            corrected_exons += [(new_introns[-1][1] + 1, read_region[1])]
        else:
            corrected_exons = [(read_region[0], read_region[1])]
        return corrected_exons

    def correct_fuzzy_junctions(self, alignment_info, read_assignment):
        intron_profile = alignment_info.combined_profile.read_intron_profile
        isoform_id = read_assignment.isoform_matches[0].assigned_transcript
        isofrom_profile = self.gene_info.intron_profiles.profiles[isoform_id]
        ref_introns = self.gene_info.intron_profiles.features

        logger.debug("Correcting fuzzy junctions in the following profiles")
        logger.debug(intron_profile.read_profile)
        logger.debug(intron_profile.gene_profile)
        logger.debug(isofrom_profile)

        new_introns = []
        read_pos, ref_pos = self.find_next(intron_profile, -1, -1)
        while read_pos < len(intron_profile.read_profile) and ref_pos < len(intron_profile.gene_profile):
            if isofrom_profile[ref_pos] != 1:
                # if intron match is ambiguous, we have to choose the one that belongs to the assigned isoform
                ref_pos += 1
            assert isofrom_profile[ref_pos] == 1

            read_intron = intron_profile.read_features[read_pos]
            ref_intron = ref_introns[ref_pos]
            assert equal_ranges(ref_intron, read_intron, self.delta)
            # TODO: check for reliability of splice sites
            new_introns.append(ref_intron)
            read_pos, ref_pos = self.find_next(intron_profile, read_pos, ref_pos)
        return (alignment_info.read_start, alignment_info.read_end), new_introns

    def find_next(self, intron_profile, read_pos, ref_pos):
        try:
            # first matched read intron
            read_pos = intron_profile.read_profile.index(1, read_pos + 1)
            # first matched isoform intron
            ref_pos = intron_profile.gene_profile.index(1, ref_pos + 1)
        except ValueError:
            read_pos = len(intron_profile.read_profile)
            ref_pos = len(intron_profile.gene_profile)
        return read_pos, ref_pos

    def correct_misalignments(self, alignment_info, read_assignment):
        event_map = {}
        for e in read_assignment.isoform_matches[0].match_subclassifications:
            if e.read_region == SupplementaryMatchConstansts.undefined_region:
                continue
            if e.read_region[0] == SupplementaryMatchConstansts.absent_position and \
                    e.event_type == MatchEventSubtype.fake_micro_intron_retention:
                event_map[-e.read_region[1]-1] = e
            event_map[e.read_region[0]] = e

        intron_profile = alignment_info.combined_profile.read_intron_profile
        read_introns = intron_profile.read_features
        read_region = (alignment_info.read_start, alignment_info.read_end)
        isoform_id = read_assignment.isoform_matches[0].assigned_transcript
        isoform_region = self.gene_info.transcript_region(isoform_id)
        isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
        return self.process_events(event_map, read_region, read_introns,  isoform_region, isoform_introns)

    def process_events(self, event_map, read_region, read_introns, isoform_region, isoform_introns):
        corrected_introns = self.intron_profile_constructor.match_genomic_features(read_introns)
        logger.debug(read_introns)
        logger.debug(corrected_introns)
        logger.debug(isoform_introns)
        corrected_read_region = read_region
        new_introns = []
        i = 0

        while i < len(corrected_introns):
            if -i-1 in event_map:
                # special case for fake IR
                logger.debug(event_map[-i-1].isoform_region)
                new_introns.append(isoform_introns[event_map[-i-1].isoform_region[0]])

            if i not in event_map:
                # TODO: check for reliability of splice sites
                new_introns.append(corrected_introns[i])
                i += 1
                continue

            event = event_map[i]
            if event.event_type == MatchEventSubtype.fake_terminal_exon_left:
                assert event.read_region[0] == event.read_region[1]
                # fake terminal exon, skip it
                corrected_read_region = (read_introns[event.read_region[0]][1]+1, corrected_read_region[1])
            elif event.event_type == MatchEventSubtype.fake_terminal_exon_right:
                assert event.read_region[0] == event.read_region[1]
                # fake terminal exon, skip it
                corrected_read_region = (corrected_read_region[0], read_introns[event.read_region[0]][0]-1)
            elif event.event_type == MatchEventSubtype.terminal_exon_misalignment_left:
                new_introns.append(isoform_introns[event.isoform_region[0]])
                corrected_read_region = (isoform_region[0], corrected_read_region[1])
            elif event.event_type == MatchEventSubtype.terminal_exon_misalignment_right:
                new_introns.append(isoform_introns[event.isoform_region[0]])
                corrected_read_region = (corrected_read_region[0], isoform_region[1])
            elif event.event_type in {MatchEventSubtype.exon_misalignment, MatchEventSubtype.intron_shift} and \
                    contains_well_inside(read_region, (isoform_introns[event.isoform_region[0]][0],
                                                       isoform_introns[event.isoform_region[1]][1]), self.params.delta):
                # misalignments but inside read region
                assert event.read_region[0] == event.read_region[1]
                new_introns += [isoform_introns[i] for i in range(event.isoform_region[0], event.isoform_region[1] + 1)]
            elif event.event_type in {MatchEventSubtype.extra_intron_known,
                                      MatchEventSubtype.intron_alternation_known,
                                      MatchEventSubtype.intron_migration,
                                      MatchEventSubtype.exon_skipping_known,
                                      MatchEventSubtype.exon_merge_known,
                                      MatchEventSubtype.terminal_exon_shift_known,
                                      MatchEventSubtype.mutually_exclusive_exons_known,
                                      MatchEventSubtype.exon_gain_known,
                                      MatchEventSubtype.exon_detach_known,
                                      MatchEventSubtype.alternative_structure_known,
                                      MatchEventSubtype.alternative_structure_novel}:
                # add corrected read intron
                new_introns += [corrected_introns[i] for i in range(event.read_region[0], event.read_region[1] + 1)]
            else:
                # add as is
                new_introns += [read_introns[i] for i in range(event.read_region[0], event.read_region[1] + 1)]

            i = event.read_region[1] + 1

        return corrected_read_region, new_introns


