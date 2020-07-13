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
from src.junction_comparator import *

logger = logging.getLogger('IsoQuant')


class TranscriptModelType(Enum):
    known = 1
    novel_in_catalog = 2
    novel_not_in_catalog = 10


# simple class for storing all information needed for GFF
class TranscriptModel:
    def __init__(self, chr_id, strand, transcript_id, reference_transcript, reference_gene, exon_blocks, transcript_type):
        self.chr_id = chr_id
        self.strand = strand
        self.transcript_id = transcript_id
        self.gene_id = reference_gene
        self.reference_transcript = reference_transcript
        self.reference_gene = reference_gene
        self.exon_blocks = exon_blocks
        self.transcript_type = transcript_type

    def get_start(self):
        return self.exon_blocks[0][0]

    def get_end(self):
        return self.exon_blocks[1][-1]


class GFFPrinter:
    def __init__(self, outf_prefix, sample_name, print_meta_features=False):
        self.out_gff = open(os.path.join(outf_prefix, sample_name + ".transcript_models.gtf"), "w")
        self.out_gff.write("# " + sample_name + " IsoQuant generated GFF\n")
        self.out_r2t = open(os.path.join(outf_prefix, sample_name + ".reads_transcript_model_map.tsv"), "w")
        self.out_r2t.write("#read_id\ttranscript_id\n")
        self.out_counts = open(os.path.join(outf_prefix, sample_name + ".transcript_models_counts.tsv"), "w")
        self.out_counts.write("#ID\t%s\n" % sample_name)
        # TODO implement meta features
        # TODO print additional information in GTF -- counts, modifications
        self.print_meta_features = print_meta_features

    def __del__(self):
        self.out_gff.close()
        self.out_r2t.close()
        self.out_counts.close()

    def dump(self, transcript_model_constructor):
        # write exons to GFF
        for model in transcript_model_constructor.transcript_model_storage:
            prefix_columns = "%s\tIsoQuant\texon\t" % model.chr_id
            suffix_columns = '.\t%s\t.\tgene_id "%s"; transcript_id "%s"; ' \
                             'reference_gene_id "%s"; reference_transcript_id "%s"\n' % \
                             (model.strand, model.gene_id, model.transcript_id,
                              model.reference_gene, model.reference_transcript)
            for e in model.exon_blocks:
                self.out_gff.write(prefix_columns + "%d\t%d\t" % (e[0], e[1]) + suffix_columns)

        # write read_id -> transcript_id map
        used_reads = set()
        for model_id in transcript_model_constructor.transcript_read_ids.keys():
            for read_id in transcript_model_constructor.transcript_read_ids[model_id]:
                used_reads.add(read_id)
                self.out_r2t.write("%s\t%s\n" % (read_id, model_id))
        for read_assignment in transcript_model_constructor.read_assignment_storage:
            if read_assignment.read_id not in used_reads:
                self.out_r2t.write("%s\t%s\n" % (read_assignment.read_id, "*"))

        for id in sorted(transcript_model_constructor.transcript_counts.keys()):
            counts = transcript_model_constructor.transcript_counts[id]
            self.out_counts.write("%s\t%.2f\n" % (id, counts))


# constructor of discovered transcript models from read assignments
class TranscriptModelConstructor:
    transcript_id_counter = 0
    transcript_prefix = "transcript_"
    known_transcript_suffix = ".known"
    nic_transcript_suffix = ".nic"
    nnic_transcript_suffix = ".nnic"

    events_to_track = {
        MatchEventSubtype.alt_donor_site_novel, MatchEventSubtype.alt_acceptor_site_novel,
        MatchEventSubtype.extra_intron,
        MatchEventSubtype.extra_intron_out_left, MatchEventSubtype.extra_intron_out_right,
        MatchEventSubtype.mutually_exclusive_exons_novel, MatchEventSubtype.exon_gain_novel,
        MatchEventSubtype.intron_retention, MatchEventSubtype.exon_skipping_novel_intron,
        MatchEventSubtype.alt_donor_site_known, MatchEventSubtype.alt_acceptor_site_known,
        MatchEventSubtype.extra_intron_known, MatchEventSubtype.intron_migration,
        MatchEventSubtype.mutually_exclusive_exons_known,
        MatchEventSubtype.exon_skipping_known_intron, MatchEventSubtype.exon_gain_known,
        MatchEventSubtype.alternative_structure_known, MatchEventSubtype.alternative_structure_novel,
        MatchEventSubtype.intron_alternation_novel, MatchEventSubtype.intron_alternation_known
    }

    def __init__(self, gene_info, read_assignment_storage, params):
        self.gene_info = gene_info
        self.read_assignment_storage = read_assignment_storage
        self.params = params
        self.transcript_model_storage = []
        self.transcript_read_ids = defaultdict(set)
        self.transcript_counts = defaultdict(float)
        self.representative_reads = set()
        self.intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.intron_profiles.features,
                                                  (self.gene_info.start, self.gene_info.end),
                                                  comparator=partial(equal_ranges, delta=self.params.delta))

    def process(self):
        # split reads into clusters
        self.construct_isoform_groups()

        # check correct assignments form reference isoforms
        for isoform_id in sorted(self.correct_matches.keys()):
            self.verify_correct_match(isoform_id, self.correct_matches[isoform_id])

        # construct novel transcripts
        candidate_model_storage = []
        for isoform_id in sorted(self.modified_isoforms_groups.keys()):
            for modification in sorted(self.modified_isoforms_groups[isoform_id].keys()):
                assignments = self.modified_isoforms_groups[isoform_id][modification]
                # logger.debug("== Processing modification cluster for isoform %s of size %d, modifications:" %
                #             (isoform_id, len(assignments)))
                # logger.debug(", ".join(["%s: %s" % (x[0], str(x[1])) for x in modification]))
                # TODO: if modification type is only extra introns, filter out those that have distant intron positions
                self.process_isoform_modifications(isoform_id, assignments, candidate_model_storage)

        # merge constructed transcripts
        self.collapse_similar_isoforms(candidate_model_storage)

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
                if match.match_classification in {MatchClassification.full_splice_match,
                                                  MatchClassification.incomplete_splice_match,
                                                  MatchClassification.mono_exon_match}:
                    self.correct_matches[isoform_id].append(read_assignment)
                else:
                    significant_events = []
                    for event in match.match_subclassifications:
                        if event.event_type in self.events_to_track:
                            significant_events.append((event.event_type, event.isoform_position))
                    if significant_events:
                        significant_events = sorted(significant_events)
                        self.modified_isoforms_groups[isoform_id][tuple(significant_events)].append(read_assignment)

        # logger.debug("Constructed %d correct clusters and %d clusters with modifications" %
        #              (len(self.correct_matches), len(self.modified_isoforms_groups)))

    # process correctly assigned reads and for a reference-identical transcript
    def verify_correct_match(self, isoform_id, assignments):
        # logger.debug("Verifying correct match to %s, cluster size %d" % (isoform_id, len(assignments)))
        unique_assignment_types = {ReadAssignmentType.unique_minor_difference,
                                   ReadAssignmentType.unique, ReadAssignmentType.ambiguous}
        unique_assignments = list(filter(lambda x: x.assignment_type in unique_assignment_types, assignments))
        if len(unique_assignments) < self.params.min_ref_supporting_reads:
            logger.debug("Not enough support")
            return

        if self.params.require_polyA:
            polyA_detected = any(a.polyA_found for a in unique_assignments)
            if not polyA_detected:
                logger.debug("No polyA found")
                return

        fsm_count = 0
        if len(self.gene_info.all_isoforms_introns[isoform_id]) > 0:
            for a in unique_assignments:
                for im in a.isoform_matches:
                    if im.assigned_transcript != isoform_id:
                        continue
                    if any(m.event_type == MatchEventSubtype.fsm for m in im.match_subclassifications):
                        fsm_count += 1
        else:
            for a in unique_assignments:
                # require only unique matches for monoexons
                if len(a.isoform_matches) != 1:
                    continue
                if any(m.event_type == MatchEventSubtype.mono_exon_match for m in a.isoform_matches[0].match_subclassifications):
                    fsm_count += 1

        if fsm_count < self.params.min_ref_fsm_supporting_reads:
            logger.debug("Not enough FSM reads")
            return

        new_transcript_model = self.transcript_from_reference(isoform_id)
        self.transcript_model_storage.append(new_transcript_model)
        # logger.debug("Created transcript model %s" % new_transcript_model.transcript_id)
        # logger.debug(new_transcript_model.exon_blocks)

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
                               self.gene_info.all_isoforms_exons[isoform_id], TranscriptModelType.known)

    # check that all splice junction in isoform are covered by at least one read
    def check_all_junctions_covered(self, isoform_id, read_assignments):
        isoform_profile = self.gene_info.intron_profiles.profiles[isoform_id]
        covered_junctions = [0 for i in range(len(isoform_profile))]
        for ra in read_assignments:
            read_profile = ra.combined_profile.read_intron_profile.gene_profile
            for i in range(len(read_profile)):
                if read_profile[i] == 1:
                    covered_junctions[i] = 1
        return all_features_present(isoform_profile, covered_junctions)

    # construct a transcript from a group of reads with the same modification
    def process_isoform_modifications(self, isoform_id, assignments, candidate_model_storage):
        remaining_assignments = copy.copy(assignments)
        while len(remaining_assignments) >= self.params.min_novel_supporting_reads:
            # choose the best representative
            # TODO: precompute them in order
            representative_read_assignment = self.select_representative_read(isoform_id, remaining_assignments)
            if not representative_read_assignment:
                # logger.debug("> No reliable representative read can be found")
                return
            # logger.debug("> Representative read chosen: %s" % representative_read_assignment.read_id)
            # logger.debug(representative_read_assignment.combined_profile.read_exon_profile.read_features)
            # logger.debug(representative_read_assignment.combined_profile.read_intron_profile.read_features)
            # create a new transcript model

            self.representative_reads.add(representative_read_assignment.read_id)
            new_transcript_model = self.blend_read_into_isoform(isoform_id, representative_read_assignment)
            if not new_transcript_model:
                # logger.debug("> No novel model was constructed")
                return
            # logger.debug("Created new candidate transcript model %s : %s " %
            #             (new_transcript_model.transcript_id, str(new_transcript_model.exon_blocks)))
            # compare read junctions with novel transcript model, count them and keep only those that do not match
            remaining_assignments = self.verify_novel_model(remaining_assignments, new_transcript_model,
                                                            representative_read_assignment.read_id,
                                                            candidate_model_storage)

    # select longest read with polyA detected
    # FIXME: use CAGE data or estimate reliability by looking at other reads
    def select_representative_read(self, isoform_id, assignments):
        if not assignments:
            return None

        strand = self.gene_info.isoform_strands[isoform_id]
        if strand == '+':
            best_read_3prime_pos = assignments[0].combined_profile.read_exon_profile.read_features[0][0]
            best_read_5prime_pos = assignments[0].combined_profile.read_exon_profile.read_features[-1][1]
        else:
            best_read_3prime_pos = assignments[0].combined_profile.read_exon_profile.read_features[-1][1]
            best_read_5prime_pos = assignments[0].combined_profile.read_exon_profile.read_features[0][0]
        best_reads = []
        # read_coords_to_assignment = {}

        for a in assignments:
            if a.read_id in self.representative_reads:
                continue
            read_exon_profile = a.combined_profile.read_exon_profile

            logger.debug("Checking whether read is reliable")
            # logger.debug(a.combined_profile.read_exon_profile.read_features[0][0],
            #             a.combined_profile.read_exon_profile.read_features[-1][1])
            # logger.debug("%s %d %d" % (a.read_id, a.combined_profile.polya_pos, a.combined_profile.polyt_pos))
            if strand == '+':
                if not self.params.require_polyA or a.combined_profile.polya_pos != -1:
                    tss = read_exon_profile.read_features[0][0]
                    tts = read_exon_profile.read_features[-1][1]
                    if tss == best_read_3prime_pos:
                        if tts == best_read_5prime_pos:
                            best_reads.append(a)
                        elif tts > best_read_5prime_pos:
                            best_reads = [a]
                    elif tss < best_read_3prime_pos:
                        best_reads = [a]
            else:
                if not self.params.require_polyA or a.combined_profile.polyt_pos != -1:
                    tss = read_exon_profile.read_features[-1][1]
                    tts = read_exon_profile.read_features[0][0]
                    if tss == best_read_3prime_pos:
                        if tts == best_read_5prime_pos:
                            best_reads.append(a)
                        elif tts < best_read_5prime_pos:
                            best_reads = [a]
                    elif tss > best_read_3prime_pos:
                        best_reads = [a]

        if not best_reads:
            logger.debug("Empty TSS array")
            return None
        return sorted(best_reads, key=lambda x: x.read_id)[0]

    def blend_read_into_isoform(self, isoform_id, read_assignment):
        # logger.debug("Creating novel transcript model for isoform %s and read %s" % (isoform_id, read_assignment.read_id))
        modification_events_map = self.derive_significant_modifications_map(isoform_id, read_assignment)
        if not modification_events_map:
            return None

        isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
        strand = self.gene_info.isoform_strands[isoform_id]

        combined_profile = read_assignment.combined_profile
        read_introns = combined_profile.read_intron_profile.read_features
        read_start, read_end = self.get_read_region(strand, combined_profile)
        novel_exons = []

        # logger.debug("Isoform I " + str(isoform_introns))
        # logger.debug("Isoform E " + str(self.gene_info.all_isoforms_exons[isoform_id]))
        # logger.debug("Read " + str(read_introns))

        if SupplementaryMatchConstansts.extra_left_mod_position in modification_events_map:
            # if there are extra introns on the left
            current_exon_start = read_start
            events = modification_events_map[SupplementaryMatchConstansts.extra_left_mod_position]
            current_exon_start = self.process_intron_related_events(events, None, isoform_introns, read_introns,
                                                                    novel_exons, current_exon_start)
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
                # logger.debug("Adding ref exon: %d, %d" % (isoform_pos, current_exon_start))
                current_exon_start = self.add_intron(novel_exons, current_exon_start, isoform_introns[isoform_pos])
                isoform_pos += 1

            else:
                current_events = modification_events_map[isoform_pos]
                current_exon_start = self.process_intron_related_events(current_events, isoform_pos, isoform_introns,
                                                                        read_introns, novel_exons, current_exon_start)
                if isoform_pos < len(isoform_introns) \
                        and current_exon_start < isoform_introns[isoform_pos][0] < read_end:
                    # intron modification was processed but nothing overlapping was added =>
                    # extra intron within previous exon => add this intron as is
                    # check that is was really extra intron
                    extra_intron_types = {MatchEventSubtype.extra_intron, MatchEventSubtype.extra_intron_known}
                    only_extra_intron = all(el.event_type in extra_intron_types for el in current_events)
                    if only_extra_intron:
                        current_exon_start = self.add_intron(novel_exons, current_exon_start, isoform_introns[isoform_pos])
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
            novel_transcript_end = read_end
        novel_exons.append((current_exon_start, novel_transcript_end))

        novel_exons = self.correct_transcripts_ends(novel_exons, combined_profile, isoform_id, modification_events_map)

        if not self.validate_exons(novel_exons):
            logger.warning("Error in novel transcript, not sorted or incorrect exon coords")
            logger.info(novel_exons)
            return None

        nnic = False
        for events in modification_events_map.values():
            nnic |= any(me.event_type in nnic_event_types for me in events)
        id_suffix = self.nnic_transcript_suffix if nnic else self.nic_transcript_suffix
        transcript_type = TranscriptModelType.novel_not_in_catalog if nnic else TranscriptModelType.novel_in_catalog
        new_transcript_id = self.transcript_prefix + str(self.get_transcript_id()) + id_suffix

        return TranscriptModel(self.gene_info.chr_id, self.gene_info.isoform_strands[isoform_id],
                               new_transcript_id, isoform_id, self.gene_info.gene_id_map[isoform_id],
                               novel_exons, transcript_type)

    # check that all exons are sorted and have correct coordinates
    def validate_exons(self, novel_exons):
        return novel_exons == sorted(novel_exons) and all(x[0] <= x[1] for x in novel_exons)

    # move transcripts ends to known ends if they are closed and no polyA found
    def correct_transcripts_ends(self, novel_exons, combined_profile, isoform_id, modification_events_map):
        isoform_end = self.gene_info.transcript_end(isoform_id)
        strand = self.gene_info.isoform_strands[isoform_id]

        if SupplementaryMatchConstansts.extra_left_mod_position not in modification_events_map and \
                0 not in modification_events_map:
            # change only if there are no extra introns on the left and first intron is not modified
            novel_transcript_start = novel_exons[0][0]
            known_isoform_start = self.gene_info.transcript_start(isoform_id)
            if (strand == "+" or combined_profile.polyt_pos == -1) and \
                    abs(novel_transcript_start - known_isoform_start) <= self.params.max_dist_to_isoforms_tsts and \
                    known_isoform_start < novel_exons[0][1]:
                novel_exons[0] = (known_isoform_start, novel_exons[0][1])

        last_index = len(self.gene_info.all_isoforms_introns[isoform_id]) - 1
        if SupplementaryMatchConstansts.extra_right_mod_position not in modification_events_map and \
                last_index not in modification_events_map:
            # change only if there are no extra introns on the right and last intron is not modified
            novel_transcript_end = novel_exons[-1][1]
            known_isoform_end = self.gene_info.transcript_end(isoform_id)
            if (strand == "-" or combined_profile.polya_pos == -1) and \
                    abs(novel_transcript_end - known_isoform_end) <= self.params.max_dist_to_isoforms_tsts and \
                    known_isoform_end > novel_exons[-1][0]:
                novel_exons[-1] = (novel_exons[-1][0], known_isoform_end)

        return novel_exons

    # process a sorted list of events assigned to the same intron
    def process_intron_related_events(self, sorted_event_list, isoform_pos, isoform_introns, read_introns,
                                      novel_exons, current_exon_start):
        # logger.debug("> Processing events for position %s: %s" % (str(isoform_pos), str(sorted_event_list)))
        # logger.debug("> Before: %d, %s" % (current_exon_start, novel_exons))
        for event in sorted_event_list:
            current_exon_start = self.process_single_event(event, isoform_pos, isoform_introns, read_introns,
                                                           novel_exons, current_exon_start)
        # logger.debug("> After: %d, %s" % (current_exon_start, novel_exons))
        return current_exon_start

    # process single event
    def process_single_event(self, event_tuple, isoform_pos, isoform_introns, read_introns, novel_exons, current_exon_start):
        # logger.debug("> > Applying event %s at position %s" % (event_tuple.event_type.name, str(isoform_pos)))
        if event_tuple.event_type == MatchEventSubtype.intron_retention:
            # simply skip reference intron
            return current_exon_start

        if event_tuple.read_region == (JunctionComparator.absent, JunctionComparator.absent):
            logger.warning("Undefined read intron position for event type: %s" % event_tuple.event_type.name)
            return current_exon_start
        read_intron = read_introns[event_tuple.read_region[0]]
        # logger.debug("Novel intron " + str(read_intron))

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
        # logger.debug("Selected modifications: " +", ".join(["%s: %s - %s" % (x.event_type.name, str(x.isoform_position), str(x.read_region))
        #                                                    for x in match_subclassifications]))
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
                sorted(modification_events_map[isoform_position], key=lambda x: x.read_region)

        if not modification_events_map:
            # logger.debug("No modification events detected for " + read_assignment.read_id)
            return None
        # logger.debug("Sorted modifications: " + ", ".join([str(x) + " - " + str(modification_events_map[x])
        #                                                    for x in sorted(modification_events_map.keys())]))
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
            if combined_profile.polyt_pos != -1:
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

    def verify_novel_model(self, read_assignments, transcript_model, original_read_id, candidate_model_storage):
        # logger.debug("Verifying transcript model %s with %d reads" % (transcript_model.transcript_id, len(read_assignments)))
        model_exons = transcript_model.exon_blocks
        isoform_start = model_exons[0][0]
        isoform_end = model_exons[-1][1]
        model_introns = junctions_from_blocks(model_exons)
        strand = transcript_model.strand
        model_intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(model_introns, (model_exons[0][0], model_exons[-1][1]),
                                                  comparator=partial(equal_ranges, delta=self.params.delta))
        intron_comparator = JunctionComparator(self.params, model_intron_profile_constructor)

        assigned_reads = []
        fsm_match_count = 0
        unassigned_reads = []
        nearby_starts_count = 0
        nearby_ends_count = 0
        for assignment in read_assignments:
            read_introns = assignment.combined_profile.read_intron_profile.read_features
            read_start, read_end = self.get_read_region(strand, assignment.combined_profile)
            start_matches = abs(read_start - isoform_start) < self.params.max_dist_to_novel_tsts
            end_matches = abs(read_end - isoform_end) < self.params.max_dist_to_novel_tsts
            # profile_matches =  all(el == 1 for el in read_profile.read_profile)

            matching_events = \
                intron_comparator.compare_junctions(read_introns, (read_start, read_end),
                                                    model_introns, (isoform_start, isoform_end))

            # logger.debug("Read %s, start %d, end %d, events %s" % (assignment.read_id, read_start, read_end, str(matching_events)))

            # check that no serious contradiction occurs
            profile_matches = True
            if len(matching_events) > 1 \
                    or (len(matching_events) == 1 and matching_events[0].event_type != MatchEventSubtype.none):
                for e in matching_events:
                    if e.event_type in nnic_event_types or e.event_type in nic_event_types:
                        profile_matches = False
                        break

            if profile_matches:
                if start_matches:
                    nearby_starts_count += 1
                if end_matches:
                    nearby_ends_count += 1
                # all read introns were mapped, read is assigned
                assigned_reads.append(assignment.read_id)
                # since profile is not reliable due to intron shifts etc
                # considering that there are no serious errors, covering all introns in enough
                is_fsm = contains((read_start, read_end), (model_introns[0][0], model_introns[-1][1]))
                if is_fsm:
                    # all introns of novel model are covered
                    fsm_match_count += 1
            else:
                unassigned_reads.append(assignment)

        # logger.debug("Stats for %s, FSM = %d, total = %d, start = %d, end = %d" %
        #             (transcript_model.transcript_id, fsm_match_count, len(assigned_reads),
        #              nearby_starts_count, nearby_ends_count))

        if fsm_match_count == 0:
            logger.warning("Zero FSM for transcript model %s" % transcript_model.transcript_id)

        if len(assigned_reads) >= self.params.min_novel_supporting_reads and \
                fsm_match_count >= self.params.min_novel_fsm_supporting_reads and \
                nearby_starts_count >= self.params.min_reads_supporting_tsts and \
                nearby_ends_count >= self.params.min_reads_supporting_tsts:
            # to confirm we need at least min_novel_supporting_reads supporting reads
            # and at least min_novel_fsm_supporting_reads FSM
            logger.debug("Successfully confirmed %s" % transcript_model.transcript_id)
            candidate_model_storage.append(transcript_model)
            for read_id in assigned_reads:
                self.transcript_read_ids[transcript_model.transcript_id].add(read_id)
            return unassigned_reads
        else:
            logger.debug("Transcript candidate %s looks unreliable" % transcript_model.transcript_id)
            all_except_original = list(filter(lambda x: x.read_id != original_read_id, read_assignments))
            return all_except_original

    def collapse_similar_isoforms(self, candidate_model_storage):
        transcript_model_map = {}
        for transcript in candidate_model_storage:
            transcript_model_map[(transcript.get_start(), transcript.get_end())] = transcript

        transcript_regions = sorted(transcript_model_map.keys())
        for t1 in transcript_regions:
            if transcript_model_map[t1] is None:
                continue
            for t2 in transcript_regions:
                if transcript_model_map[t1] == transcript_model_map[t2]:
                    continue
                if t2[0] > t1[1]:
                    break
                if transcript_model_map[t2] is None:
                    continue

                if not contains(t1, t2):
                    continue

                if self.check_if_subisoform(transcript_model_map[t1], transcript_model_map[t2]):
                    logger.debug("Detected subisoform, will collapse")
                    for read_id in self.transcript_read_ids[transcript_model_map[t2].transcript_id]:
                        self.transcript_read_ids[transcript_model_map[t1].transcript_id].add(read_id)
                    del self.transcript_read_ids[transcript_model_map[t2].transcript_id]
                    transcript_model_map[t2] = None

        for k in transcript_model_map:
            model = transcript_model_map[k]
            if model is None:
                continue
            self.transcript_model_storage.append(model)
            id = model.transcript_id
            self.transcript_counts[id] = len(self.transcript_read_ids[id])

    # check in one isoform can be collapsed into another
    def check_if_subisoform(self, big_transcript_model, small_transcript_model):
        big_region = (big_transcript_model.get_start(), big_transcript_model.get_end())
        small_region = (small_transcript_model.get_start(), small_transcript_model.get_end())
        if not contains(big_region, small_region):
            return False

        logger.debug("Checking similarity between transcript models %s and %s" %
                     (big_transcript_model.transcript_id, small_transcript_model.transcript_id))

        if big_transcript_model.strand != small_transcript_model.strand:
            return False

        if self.params.require_polyA:
            start_diff = small_region[0] - big_region[0]
            end_diff = big_region[1] - small_region[1]
            if big_transcript_model.strand == "+":
                if end_diff > self.params.max_dist_to_isoforms_tsts:
                    logger.debug("polyA supported ends differ significantly")
                    return False
            else:
                if start_diff > self.params.max_dist_to_isoforms_tsts:
                    logger.debug("polyT supported starts differ significantly")
                    return False

        model_exons = big_transcript_model.exon_blocks
        model_introns = junctions_from_blocks(model_exons)
        model_intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(model_introns, big_region, equal_ranges)

        exons_to_check = small_transcript_model.exon_blocks
        introns_to_check = junctions_from_blocks(exons_to_check)
        profile = model_intron_profile_constructor.construct_profile_for_features(introns_to_check, small_region)

        if all(el == 1 for el in profile.read_profile):
            # profile matches perfectly
            if self.params.collapse_subisoforms:
                return all(el != -1 for el in profile.gene_profile)
            else:
                return all(el == 1 for el in profile.gene_profile)
        return False

