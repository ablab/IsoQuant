############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
import copy

from src.isoform_assignment import *

logger = logging.getLogger('IsoQuant')


# simple class for storing all information needed for GFF
class TranscriptModel:
    def __init__(self, chr_id, strand, transcript_id, reference_transcript, reference_gene, exon_blocks):
        self.chr_id = chr_id
        self.strand = strand
        self.transcript_id = transcript_id
        self.gene_id = reference_gene
        self.reference_transcript = reference_transcript
        self.reference_gene = reference_gene
        self.exon_blocks = exon_blocks


# constructor of discovered transcript models from read assignments
class TranscriptModelConstructor:
    transcript_id_counter = 0
    transcript_prefix = "transcript_"
    known_transcript_suffix = ".known"
    nic_transcript_suffix = ".nic"
    nnic_transcript_suffix = ".nnic"

    events_to_track = {MatchEventSubtype.alt_donor_site_novel, MatchEventSubtype.alt_acceptor_site_novel,
                       MatchEventSubtype.extra_intron,
                       MatchEventSubtype.extra_intron_out_left, MatchEventSubtype.extra_intron_out_right,
                       MatchEventSubtype.mutually_exclusive_exons_novel, MatchEventSubtype.exon_gain_novel,
                       MatchEventSubtype.intron_retention,
                       MatchEventSubtype.alt_donor_site, MatchEventSubtype.alt_acceptor_site,
                       MatchEventSubtype.extra_intron_known, MatchEventSubtype.intron_migration,
                       MatchEventSubtype.mutually_exclusive_exons_known,
                       MatchEventSubtype.exon_skipping_known_intron, MatchEventSubtype.exon_gain_known,
                       MatchEventSubtype.alternative_structure_known}

    def __init__(self, gene_info, read_assignment_storage, params):
        self.gene_info = gene_info
        self.read_assignment_storage = read_assignment_storage
        self.params = params
        self.transcript_model_storage = []
        self.transcript_read_ids = defaultdict(list)
        self.transcript_counts = defaultdict(float)

    def process(self):
        # split reads into clusters
        self.construct_isoform_groups()

        # check correct assignments form reference isoforms
        for isoform_id in self.correct_matches.keys():
            self.verify_correct_match(isoform_id, self.correct_matches[isoform_id])

        # construct novel trasncripts
        for isoform_id in self.modified_isoforms_groups.keys():
            for modification in self.modified_isoforms_groups[isoform_id].keys():
                assignments = self.modified_isoforms_groups[isoform_id][modification]
                logger.debug("== Processing modidication cluster for isoform %s of size %d, modifications:" %
                             (isoform_id, len(assignments)))
                logger.debug(", ".join([x.event_type.name + ":" + str(x.position) for x in modification]))
                self.process_isoform_modifications(isoform_id, assignments)

        # merge constructed transcripts TODO
        pass

    def get_transcript_id(self):
        TranscriptModelConstructor.transcript_id_counter += 1
        return TranscriptModelConstructor.transcript_id_counter

    # group reads by the isoforms and modification events
    def construct_isoform_groups(self):
        logger.debug("Constructing isoform groups")
        isoform_groups = defaultdict(list)
        for read_assignment in self.read_assignment_storage:
            for match in read_assignment.isoform_matches:
                isoform_groups[match.assigned_transcript].append(read_assignment)

        self.modified_isoforms_groups = defaultdict(lambda: defaultdict(list))
        self.correct_matches = defaultdict(list)
        for isoform_id in isoform_groups.keys():
            isoform_assignments = isoform_groups[isoform_id]
            for read_assignment in isoform_assignments:
                for match in read_assignment.isoform_matches:
                    if match.match_classification in {MatchClassification.full_splice_match, MatchClassification.incomplete_splice_match}:
                        self.correct_matches[isoform_id].append(read_assignment)
                        break
                    else:
                        significant_events = []
                        for event in match.match_subclassifications:
                            if event.event_type in self.events_to_track:
                                significant_events.append(event)
                        if significant_events:
                            self.modified_isoforms_groups[isoform_id][tuple(significant_events)].append(read_assignment)
        logger.debug("Constructed %d correct clusters and %d clusters with modifications" %
                     (len(self.correct_matches), len(self.modified_isoforms_groups)))

    # process correctly assigned reads and for a reference-identical transcript
    def verify_correct_match(self, isoform_id, assignments):
        logger.debug("Verifying correct match to %s, cluster size %d" % (isoform_id, len(assignments)))
        unique_assignments = list(filter(lambda x:x.assignment_type in
                                                  {ReadAssignmentType.unique_minor_difference, ReadAssignmentType.unique},
                                         assignments))
        if len(unique_assignments) < self.params.min_ref_supporting_reads:
            logger.debug("Not enough support")
            return

        polyA_detected = any(a.polyA_found for a in unique_assignments)
        if not polyA_detected:
            logger.debug("No polyA found")
            return

        fsm_detected = any(any(x.event_type == MatchEventSubtype.fsm
                               for x in a.isoform_matches[0].match_subclassifications)
                           for a in unique_assignments)
        if not fsm_detected:
            logger.debug("No FSM")
            all_covered = self.check_all_juctions_covered(isoform_id, unique_assignments)
            if not all_covered or len(unique_assignments) < 2 * self.params.min_ref_supporting_reads:
                logger.debug("And not all junctions are covered enough")
                return

        new_transcript_model = self.transcript_from_reference(isoform_id)
        self.transcript_model_storage.append(new_transcript_model)
        logger.debug("Created transcript model %s" % new_transcript_model.transcript_id)
        logger.debug(new_transcript_model.exon_blocks)

        assignments_to_consider = assignments if self.params.count_ambiguous else unique_assignments
        new_transcript_id = new_transcript_model.transcript_id
        for assignment in assignments_to_consider:
            self.transcript_read_ids[new_transcript_id].append(assignment.read_id)
            self.transcript_counts[new_transcript_id] += 1.0 / float(len(assignment.isoform_matches))

    # create transcript model object from reference isoforms
    def transcript_from_reference(self, isoform_id):
        new_transcript_id = self.transcript_prefix + str(self.get_transcript_id()) + self.known_transcript_suffix
        return TranscriptModel(self.gene_info.chr_id, self.gene_info.isoform_strands[isoform_id],
                               new_transcript_id, isoform_id, self.gene_info.gene_id_map[isoform_id],
                               self.gene_info.all_isoforms_exons[isoform_id])

    # check that all splice jusction in isoform are covered by at least one read
    def check_all_juctions_covered(self, isoform_id, read_assignments):
        isoform_profile = self.gene_info.intron_profiles.profiles[isoform_id]
        covered_junctions = [0 for i in range(len(isoform_profile))]
        for ra in read_assignments:
            read_profile = ra.combined_profile.read_intron_profile.gene_profile
            for i in range(len(read_profile)):
                if read_profile[i] == 1:
                    covered_junctions[i] = 1
        return all_features_present(isoform_profile, covered_junctions)

    # construct a transcript from a group of reads with the same modification
    def process_isoform_modifications(self, isoform_id, assignments):
        remaining_assignments = copy.copy(assignments)
        while len(remaining_assignments) >= self.params.min_alt_supporting_reads:
            # choose the best representative
            representative_read_assignment = self.select_representative_read(isoform_id, remaining_assignments)
            if not representative_read_assignment:
                logger.debug("> No reliable representative read can be found")
                return
            logger.debug("> Representative read chosen: %s" % representative_read_assignment.read_id)
            # create a new transcript model
            new_transcript_model = self.blend_read_into_isoform(isoform_id, representative_read_assignment)
            if not new_transcript_model:
                logger.debug("> No novel model was constructed")
                return
            # compare read junctions with novel transcript model, count them and keep only those that do not match
            assigned_reads, remaining_assignments = self.assign_reads_to_model(remaining_assignments, new_transcript_model)
            if not assigned_reads:
                logger.warning("> No reads were assigned to novel transcript " + new_transcript_model.transcript_id)
                return
            if len(assigned_reads) < self.params.min_alt_supporting_reads:
                continue

            for assignment in assigned_reads:
                self.transcript_read_ids[isoform_id].append(assignment.read_id)


    # select longest read with polyA detected
    # FIXME: use CAGE data or estimate reliability by looking at other reads
    def select_representative_read(self, isoform_id, assignments):
        strand = self.gene_info.isoform_strands[isoform_id]
        read_coords_to_assignment = {}
        for a in assignments:
            if strand == '+':
                tss = a.combined_profile.read_exon_profile.read_features[0][0]
                tts = a.combined_profile.polya_pos
                if tts != -1:
                    read_coords_to_assignment[tss] = a
            else:
                tss = a.combined_profile.read_exon_profile.read_features[-1][1]
                tts = a.combined_profile.polyt_pos
                if tts != -1:
                    read_coords_to_assignment[tss] = a

        tss_positions = sorted(read_coords_to_assignment.keys())
        if not tss_positions:
            return None

        if strand == "+":
            return read_coords_to_assignment[tss_positions[0]]
        else:
            return read_coords_to_assignment[tss_positions[-1]]

    def blend_read_into_isoform(self, isoform_id, read_assignment):
        logger.debug("Creating novel transcript model for isoform %s and read %s" % (isoform_id, read_assignment.read_id))
        modification_events = self.derive_significant_modifications(isoform_id, read_assignment)
        if not modification_events:
            return

        isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
        isoform_start = self.gene_info.transcript_start(isoform_id)
        isoform_end = self.gene_info.transcript_end(isoform_id)
        strand = self.gene_info.isoform_strands[isoform_id]

        combined_profile = read_assignment.combined_profile
        read_introns = combined_profile.read_intron_profile.read_features
        read_start, read_end = self.get_read_region(strand, combined_profile)

        novel_exons = []
        current_exon_start = -1
        read_pos = 0
        isoform_pos = 0
        event_pos = 0

        if MatchEventSubtype.extra_intron_out_left in modification_events.values():
            # if there are extra introns on the left
            current_exon_start = read_start
            while read_pos < len(read_introns) and read_introns[read_pos][1] < isoform_start:
                novel_exons.append((current_exon_start, read_introns[read_pos][0] - 1))
                current_exon_start = read_introns[read_pos][1] + 1
                read_pos += 1
        else:
            # select optimal start
            if (strand == "+" or combined_profile.polya_pos == -1) and \
                    abs(read_start - isoform_start) <= self.params.max_exon_extension:
                current_exon_start = isoform_start
            else:
                current_exon_start = read_start

        # additional introns inside first isoform exon
        if -1 in modification_events:
            # add reads introns
            while read_introns[read_pos][1] < isoform_introns[isoform_pos][0]:
                if

        while read_pos < len(read_introns) and isoform_pos < len(isoform_introns):
            if
            current_event = modification_events[event_pos]
            if current_event[1] in {MatchEventSubtype.extra_intron, MatchEventSubtype.extra_intron_known}:
                assert isoform_pos == current_event[0] + 1
                novel_exons.append((current_exon_start, read_introns[read_pos][0] - 1))
                current_exon_start = read_introns[read_pos][1] + 1
                event_pos += 1
                read_pos += 1
            elif current_event[1] == MatchEventSubtype.intron_retention:
                assert isoform_pos == current_event[0]
                isoform_pos += 1
            elif current_event[1]



    def derive_significant_modifications(self, isoform_id, read_assignment):
        match_subclassifications = None
        for match in read_assignment.isoform_matches:
            if match.assigned_transcript == isoform_id:
                match_subclassifications = match.match_subclassifications
                break

        match_subclassifications = list(filter(lambda m: m.event_type in self.events_to_track, match_subclassifications))
        logger.debug("Selected modifications: " + ", ".join([x.event_type.name + ":" + str(x.position) for x in match_subclassifications]))
        modification_events = {}
        for x in match_subclassifications:
            modification_events[x.position] = x.event_type
        if not modification_events:
            logger.debug("No modification events detected for " + read_assignment.read_id)
            return None
        logger.debug("Sorted modifications: " + ", ".join([str(x) + " - " + modification_events[x].name for x in sorted(modification_events.keys())]))
        return modification_events

    def get_read_region(self, strand, combined_profile):
        if strand == "+":
            read_start = combined_profile.read_exon_profile.read_features[0][0]
            if combined_profile.polya_pos != -1:
                read_end = combined_profile.polya_pos
            else:
                read_end = combined_profile.read_exon_profile.read_features[-1][1]
        else:
            read_end = combined_profile.read_exon_profile.read_features[-1][1]
            if combined_profile.polya_pos != -1:
                read_start = combined_profile.polyt_pos
            else:
                read_start = combined_profile.read_exon_profile.read_features[0][0]
        return read_start, read_end

    def assign_reads_to_model(self, read_assignments, transcript_model):
        return None, None







