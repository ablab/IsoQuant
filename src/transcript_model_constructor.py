############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
import copy

from src.isoform_assignment import *
from src.long_read_profiles import *

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
                       MatchEventSubtype.alt_donor_site_known, MatchEventSubtype.alt_acceptor_site_known,
                       MatchEventSubtype.extra_intron_known, MatchEventSubtype.intron_migration,
                       MatchEventSubtype.mutually_exclusive_exons_known,
                       MatchEventSubtype.exon_skipping_known_intron, MatchEventSubtype.exon_gain_known,
                       MatchEventSubtype.alternative_structure_known}

    def __init__(self, gene_info, read_assignment_storage, params):
        self.gene_info = gene_info
        self.read_assignment_storage = read_assignment_storage
        self.params = params
        self.transcript_model_storage = []
        self.transcript_read_ids = defaultdict(set)
        self.transcript_counts = defaultdict(float)
        self.intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.intron_profiles.features,
                                                  (self.gene_info.start, self.gene_info.end),
                                                  comparator = partial(equal_ranges, delta = self.params.delta))

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
                logger.debug(", ".join(["%s: %s - %s" % (x.event_type.name, str(x.isoform_position), str(x.read_region))
                                        for x in modification]))
                self.process_isoform_modifications(isoform_id, assignments)

        # merge constructed transcripts TODO
        pass

    def get_transcript_id(self):
        TranscriptModelConstructor.transcript_id_counter += 1
        return TranscriptModelConstructor.transcript_id_counter

    # group reads by the isoforms and modification events
    def construct_isoform_groups(self):
        logger.debug("Constructing isoform groups")
        self.modified_isoforms_groups = defaultdict(lambda: defaultdict(list))
        self.correct_matches = defaultdict(list)

        for read_assignment in self.read_assignment_storage:
            for match in read_assignment.isoform_matches:
                isoform_id = match.assigned_transcript
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

        if self.params.require_polyA:
            polyA_detected = any(a.polyA_found for a in unique_assignments)
            if not polyA_detected:
                logger.debug("No polyA found")
                return

        fsm_count = 0
        for a in unique_assignments:
            for m in a.isoform_matches[0].match_subclassifications:
                if m.event_type == MatchEventSubtype.fsm:
                    fsm_count += 1
        if fsm_count < self.params.min_ref_fsm_supporting_reads:
            logger.debug("Not enough FSM reads")
            return

        new_transcript_model = self.transcript_from_reference(isoform_id)
        self.transcript_model_storage.append(new_transcript_model)
        logger.debug("Created transcript model %s" % new_transcript_model.transcript_id)
        logger.debug(new_transcript_model.exon_blocks)

        assignments_to_consider = assignments if self.params.count_ambiguous else unique_assignments
        new_transcript_id = new_transcript_model.transcript_id
        for assignment in assignments_to_consider:
            self.transcript_read_ids[new_transcript_id].add(assignment.read_id)
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
        while len(remaining_assignments) >= self.params.min_novel_supporting_reads:
            # choose the best representative
            representative_read_assignment = self.select_representative_read(isoform_id, remaining_assignments)
            if not representative_read_assignment:
                logger.debug("> No reliable representative read can be found")
                return
            logger.debug("> Representative read chosen: %s" % representative_read_assignment.read_id)
            logger.debug(representative_read_assignment.combined_profile.read_exon_profile.read_features)
            logger.debug(representative_read_assignment.combined_profile.read_intron_profile.read_features)
            # create a new transcript model
            new_transcript_model = self.blend_read_into_isoform(isoform_id, representative_read_assignment)
            if not new_transcript_model:
                logger.debug("> No novel model was constructed")
                return
            logger.debug("Created new candidate transcript model %s : %s " %
                         (new_transcript_model.transcript_id, str(new_transcript_model.exon_blocks)))
            # compare read junctions with novel transcript model, count them and keep only those that do not match
            remaining_assignments = self.verify_novel_model(remaining_assignments, new_transcript_model,
                                                            representative_read_assignment.read_id)


    # select longest read with polyA detected
    # FIXME: use CAGE data or estimate reliability by looking at other reads
    def select_representative_read(self, isoform_id, assignments):
        strand = self.gene_info.isoform_strands[isoform_id]
        read_coords_to_assignment = {}
        for a in assignments:
            logger.debug("Checking whether read is reliable")
            #logger.debug(a.combined_profile.read_exon_profile.read_features[0][0],
            #             a.combined_profile.read_exon_profile.read_features[-1][1])
            logger.debug("%s %d %d" % (a.read_id, a.combined_profile.polya_pos, a.combined_profile.polyt_pos))
            if strand == '+':
                tss = a.combined_profile.read_exon_profile.read_features[0][0]
                tts = a.combined_profile.polya_pos
                if not self.params.require_polyA or tts != -1:
                    read_coords_to_assignment[tss] = a
            else:
                tss = a.combined_profile.read_exon_profile.read_features[-1][1]
                tts = a.combined_profile.polyt_pos
                if not self.params.require_polyA or tts != -1:
                    read_coords_to_assignment[tss] = a

        tss_positions = sorted(read_coords_to_assignment.keys())
        if not tss_positions:
            logger.debug("Empty TSS array")
            return None

        if strand == "+":
            return read_coords_to_assignment[tss_positions[0]]
        else:
            return read_coords_to_assignment[tss_positions[-1]]

    def blend_read_into_isoform(self, isoform_id, read_assignment):
        logger.debug("Creating novel transcript model for isoform %s and read %s" % (isoform_id, read_assignment.read_id))
        modification_events_map = self.derive_significant_modifications_map(isoform_id, read_assignment)
        if not modification_events_map:
            return None

        isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
        isoform_start = self.gene_info.transcript_start(isoform_id)
        isoform_end = self.gene_info.transcript_end(isoform_id)
        strand = self.gene_info.isoform_strands[isoform_id]

        combined_profile = read_assignment.combined_profile
        read_introns = combined_profile.read_intron_profile.read_features
        read_start, read_end = self.get_read_region(strand, combined_profile)
        novel_exons = []

        logger.debug("Isoform I " + str(isoform_introns))
        logger.debug("Isoform E " + str(self.gene_info.all_isoforms_exons[isoform_id]))
        logger.debug("Read " + str(read_introns))

        if SupplementaryMatchConstansts.extra_left_mod_position in modification_events_map:
            # if there are extra introns on the left
            current_exon_start = read_start
            events = modification_events_map[SupplementaryMatchConstansts.extra_left_mod_position]
            current_exon_start = self.process_intron_related_events(events, None, isoform_introns, read_introns,
                                                                    novel_exons, current_exon_start)
        else:
            # select optimal start
            if (strand == "+" or combined_profile.polyt_pos == -1) and \
                    abs(read_start - isoform_start) <= self.params.max_dist_to_isoforms_tsts:
                #TODO check they are not separated by intron
                current_exon_start = isoform_start
            else:
                current_exon_start = read_start

        isoform_pos = 0
        while isoform_pos <= len(isoform_introns):
            if isoform_pos not in modification_events_map.keys():
                if isoform_pos == len(isoform_introns):
                    # such position is possible only when extra intron is present inside last reference exon
                    break
                if isoform_introns[isoform_pos][0] < current_exon_start:
                    # skip introns that ourside of gene region
                    isoform_pos += 1
                    continue
                if isoform_introns[isoform_pos][1] > read_end:
                    # skip introns that ourside of gene region
                    break

                # simply select reference isoform intron
                exon = (current_exon_start, isoform_introns[isoform_pos][0] - 1)
                logger.debug("Adding ref exon: %d, %d, %s" % (isoform_pos, current_exon_start, exon))
                novel_exons.append(exon)
                current_exon_start = isoform_introns[isoform_pos][1] + 1
                isoform_pos += 1

            else:
                current_events = modification_events_map[isoform_pos]
                current_exon_start = self.process_intron_related_events(current_events, isoform_pos, isoform_introns,
                                                                        read_introns, novel_exons, current_exon_start)
                isoform_pos += 1
                while isoform_pos < len(isoform_introns) and isoform_introns[isoform_pos][0] < current_exon_start:
                    isoform_pos += 1

        if SupplementaryMatchConstansts.extra_right_mod_position in modification_events_map:
            # if there are extra introns on the right
            events = modification_events_map[SupplementaryMatchConstansts.extra_right_mod_position]
            current_exon_start = self.process_intron_related_events(events, None, isoform_introns, read_introns,
                                                                    novel_exons, current_exon_start)
            novel_transcript_end = read_end
        else:
            # select optimal start
            if (strand == "-" or combined_profile.polya_pos == -1) and \
                    abs(read_end - isoform_end) <= self.params.max_dist_to_isoforms_tsts:
                # TODO check they are not separated by intron
                novel_transcript_end = isoform_end
            else:
                novel_transcript_end = read_end
        novel_exons.append((current_exon_start, novel_transcript_end))

        nnic = False
        for events in modification_events_map.values():
            nnic |= any(me.event_type in nnic_event_types for me in events)
        id_suffix = self.nnic_transcript_suffix if nnic else self.nic_transcript_suffix
        new_transcript_id = self.transcript_prefix + str(self.get_transcript_id()) + id_suffix

        return TranscriptModel(self.gene_info.chr_id, self.gene_info.isoform_strands[isoform_id],
                               new_transcript_id, isoform_id, self.gene_info.gene_id_map[isoform_id],
                               novel_exons)

    # process a sorted list of events assigned to the same intron
    def process_intron_related_events(self, sorted_event_list, isoform_pos, isoform_introns, read_introns,
                                      novel_exons, current_exon_start):
        logger.debug("> Processing events for position %s: %s" % (str(isoform_pos), str(sorted_event_list)))
        logger.debug("> Before: %d, %s" % (current_exon_start, novel_exons))
        for event in sorted_event_list:
            current_exon_start = self.procces_signle_event(event, isoform_pos, isoform_introns, read_introns,
                                                           novel_exons, current_exon_start)
        logger.debug("> After: %d, %s" % (current_exon_start, novel_exons))
        return current_exon_start

    # process single event
    def procces_signle_event(self, event_tuple, isoform_pos, isoform_introns, read_introns, novel_exons, current_exon_start):
        logger.debug("> > Applying event %s at position %s" % (event_tuple.event_type.name, str(isoform_pos)))
        if event_tuple.event_type == MatchEventSubtype.intron_retention:
            # simply skip reference intron
            return current_exon_start

        assert event_tuple.read_region is not None
        read_intron = read_introns[event_tuple.read_region[0]]
        logger.debug("Novel intron " + str(read_intron))

        if event_tuple.event_type == MatchEventSubtype.extra_intron:
            return self.add_intron(novel_exons, current_exon_start, read_intron)
        elif event_tuple.event_type == MatchEventSubtype.extra_intron_known:
            corrected_intron = self.get_closest_ref_intron(read_intron)
            return self.add_intron(novel_exons, current_exon_start, corrected_intron)
        elif event_tuple.event_type in {MatchEventSubtype.extra_intron_out_left,
                                        MatchEventSubtype.extra_intron_out_right}:
            # simply insert several reads introns
            for read_pos in range(event_tuple.read_region[0], event_tuple.read_region[1] + 1):
                current_exon_start = self.add_intron(novel_exons, current_exon_start, read_introns[read_pos])
            return current_exon_start

        isoform_intron = isoform_introns[isoform_pos]
        assert overlaps(read_intron, isoform_intron)

        if event_tuple.event_type == MatchEventSubtype.alt_donor_site_novel:
            # TODO check strands for acceptor donor sites
            novel_intron = (read_intron[0], isoform_intron[1])
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type == MatchEventSubtype.alt_donor_site_known:
            novel_intron = self.get_closest_ref_intron((read_intron[0], isoform_intron[1]))
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type == MatchEventSubtype.alt_acceptor_site_novel:
            novel_intron = (isoform_intron[0], read_intron[1])
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type == MatchEventSubtype.alt_acceptor_site_known:
            novel_intron = self.get_closest_ref_intron((isoform_intron[0], read_intron[1]))
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type in {MatchEventSubtype.intron_alternation_novel,
                                        MatchEventSubtype.exon_skipping_novel_intron}:
            # simply add read intron
            # FIXME move to lower condition
            novel_intron = (read_intron[0], read_intron[1])
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type in {MatchEventSubtype.intron_alternation_known, MatchEventSubtype.intron_migration,
                                        MatchEventSubtype.exon_skipping_known_intron}:
            # simply add corrected read intron
            # FIXME move to lower condition
            novel_intron = self.get_closest_ref_intron((read_intron[0], read_intron[1]))
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type in {MatchEventSubtype.mutually_exclusive_exons_novel,
                                        MatchEventSubtype.exon_gain_novel,
                                        MatchEventSubtype.alternative_structure_novel}:
            # simply insert several reads introns
            for read_pos in range(event_tuple.read_region[0], event_tuple.read_region[1] + 1):
                current_exon_start = self.add_intron(novel_exons, current_exon_start, read_introns[read_pos])
        elif event_tuple.event_type in {MatchEventSubtype.mutually_exclusive_exons_known,
                                        MatchEventSubtype.exon_gain_known,
                                        MatchEventSubtype.alternative_structure_known}:
            # insert several reads introns my fitting them onto reference introns
            for read_pos in range(event_tuple.read_region[0], event_tuple.read_region[1] + 1):
                # TODO speed up
                novel_intron = self.get_closest_ref_intron(read_introns[read_pos])
                current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        else:
            logger.warning("Unsupported event type %s" % event_tuple.event_type.name)

        return current_exon_start

    # return map: isoform position -> event tuple
    def derive_significant_modifications_map(self, isoform_id, read_assignment):
        match_subclassifications = None
        for match in read_assignment.isoform_matches:
            if match.assigned_transcript == isoform_id:
                match_subclassifications = match.match_subclassifications
                break

        match_subclassifications = list(filter(lambda m: m.event_type in self.events_to_track, match_subclassifications))
        logger.debug("Selected modifications: " +", ".join(["%s: %s - %s" % (x.event_type.name, str(x.isoform_position), str(x.read_region))
                                                            for x in match_subclassifications]))
        if not self.params.report_intron_retention and \
                all(m.event_type == MatchEventSubtype.intron_retention for m in match_subclassifications):
            return None

        modification_events_map = defaultdict(list)
        for x in match_subclassifications:
            modification_events_map[x.isoform_position].append(x)
        for isoform_position in modification_events_map.keys():
            if len(modification_events_map[isoform_position]) == 1:
                continue
            modification_events_map[isoform_position] = \
                sorted(modification_events_map[isoform_position], key=lambda x:x.read_region)

        if not modification_events_map:
            logger.debug("No modification events detected for " + read_assignment.read_id)
            return None
        logger.debug("Sorted modifications: " + ", ".join([str(x) + " - " + str(modification_events_map[x])
                                                           for x in sorted(modification_events_map.keys())]))
        return modification_events_map

    # get tentative transcript start and end based on polyA and mapping coordinates
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

    def get_closest_ref_intron(self, read_intron):
        # TODO speed up - binray search or interval tree
        intron_profile = self.intron_profile_constructor.construct_profile_for_features([read_intron])
        matched_intron = intron_profile.gene_profile.index(1)
        return self.intron_profile_constructor.known_features[matched_intron]

    def add_intron(self, novel_exons, current_exon_start, intron):
        exon = (current_exon_start, intron[0] - 1)
        novel_exons.append(exon)
        return intron[1] + 1

    def verify_novel_model(self, read_assignments, transcript_model, original_read_id):
        logger.debug("Verifying transcript model %s with %d reads" % (transcript_model.transcript_id, len(read_assignments)))
        model_exons = transcript_model.exon_blocks
        isoform_start = model_exons[0][0]
        isoform_end = model_exons[-1][1]
        model_introns = junctions_from_blocks(model_exons)
        strand = transcript_model.strand
        model_intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(model_introns, (model_exons[0][0], model_exons[-1][1]),
                                                  comparator = partial(equal_ranges, delta = self.params.delta))

        assigned_reads = []
        fsm_match_count = 0
        unassigned_reads = []
        nearby_starts_count = 0
        nearby_ends_count = 0
        for assignment in read_assignments:
            read_region = (self.get_read_region(strand, assignment.combined_profile))
            read_introns = assignment.combined_profile.read_intron_profile.read_features
            read_profile = \
                model_intron_profile_constructor.construct_profile_for_features(read_introns, read_region)

            read_start, read_end = self.get_read_region(strand, assignment.combined_profile)
            start_matches = abs(read_start - isoform_start) < self.params.max_dist_to_novel_tsts
            end_matches = abs(read_end - isoform_end) < self.params.max_dist_to_novel_tsts
            profile_matches =  all(el == 1 for el in read_profile.read_profile)

            if profile_matches:
                if start_matches:
                    nearby_starts_count += 1
                if end_matches:
                    nearby_ends_count += 1
                # all read introns were mapped, read is assigned
                assigned_reads.append(assignment.read_id)
                if all(el == 1 for el in read_profile.read_profile):
                    # all introns of novel model are covered
                    fsm_match_count += 1

            else:
                unassigned_reads.append(assignment)

            logger.debug("Transcript model %s" % transcript_model.transcript_id)
            logger.debug(model_introns)
            logger.debug("Read %s" % assignment.read_id)
            logger.debug(assignment.combined_profile.read_exon_profile.read_features)
            logger.debug(read_introns)
            logger.debug(read_profile.gene_profile)
            logger.debug(read_profile.read_profile)

        # TODO remove temp assert
        assert fsm_match_count > 0

        logger.debug("Stats for %s, FSM = %d, total = %d, start = %d, end = %d" %
                     (transcript_model.transcript_id, fsm_match_count, len(assigned_reads),
                      nearby_starts_count, nearby_ends_count))

        if len(assigned_reads) >= self.params.min_novel_supporting_reads and \
                fsm_match_count >= self.params.min_novel_fsm_supporting_reads and \
                nearby_starts_count >= self.params.min_reads_supporting_tsts and \
                nearby_ends_count >= self.params.min_reads_supporting_tsts:
            # to confirm we need at least min_novel_supporting_reads supporting reads
            # and at least min_novel_fsm_supporting_reads FSM
            logger.debug("Successfully confirmed %s" % (transcript_model.transcript_id))
            self.transcript_model_storage.append(transcript_model)
            for read_id in assigned_reads:
                self.transcript_read_ids[transcript_model.transcript_id].add(read_id)
            return unassigned_reads
        else:
            logger.debug("Transcript candidate %s looks unreliable" % transcript_model.transcript_id)
            all_except_original = list(filter(lambda x:x.read_id != original_read_id, read_assignments))
            return all_except_original


