############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
from collections import namedtuple
from functools import reduce
import copy

from src.isoform_assignment import *
from src.long_read_profiles import *
from src.junction_comparator import *
from src.long_read_assigner import *
from src.gene_info import *

logger = logging.getLogger('IsoQuant')


class GFFPrinter:
    def __init__(self, outf_prefix, sample_name, print_meta_features=False):
        self.model_fname = os.path.join(outf_prefix, sample_name + ".transcript_models.gtf")
        self.out_gff = open(self.model_fname, "w")
        self.out_gff.write("# " + sample_name + " IsoQuant generated GFF\n")
        self.out_r2t = open(os.path.join(outf_prefix, sample_name + ".transcript_models_reads.tsv"), "w")
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
        gene_to_model_dict = defaultdict(list)
        gene_regions = transcript_model_constructor.gene_info.get_gene_regions()
        GFFGeneInfo = namedtuple("GFFGeneInfo", ("chr_id", "strand", "gene_region"))
        gene_info_dict = {}

        for i, model in enumerate(transcript_model_constructor.transcript_model_storage):
            gene_id = model.gene_id
            gene_to_model_dict[gene_id].append(i)

            transcript_region = (model.exon_blocks[0][0], model.exon_blocks[-1][1])
            if gene_id not in gene_info_dict:
                assert model.chr_id == transcript_model_constructor.gene_info.chr_id
                gene_info_dict[gene_id] = GFFGeneInfo(model.chr_id, model.strand,
                                                      max_range(gene_regions[gene_id], transcript_region))
            else:
                gene_record = gene_info_dict[gene_id]
                assert model.chr_id == gene_record.chr_id
                assert model.strand == gene_record.strand
                gene_info_dict[gene_id] = GFFGeneInfo(model.chr_id, model.strand,
                                                      max_range(gene_record.gene_region, transcript_region))

        gene_order = sorted([(g, gene_info_dict[g].gene_region) for g in gene_info_dict.keys()], key=lambda x:x[1])

        for gene_id, coords in gene_order:
            gene_line = '%s\tIsoQuant\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcripts: %d;\n' % \
                        (gene_info_dict[gene_id].chr_id, coords[0], coords[1], gene_info_dict[gene_id].strand,
                         gene_id, len(gene_to_model_dict[gene_id]))
            self.out_gff.write(gene_line)

            for model_index in gene_to_model_dict[gene_id]:
                model = transcript_model_constructor.transcript_model_storage[model_index]
                assert model.gene_id == gene_id

                # TODO add infromation about counts and modifications for novel transcripts
                transcript_line = '%s\tIsoQuant\ttranscript\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; ' \
                                  'reference_gene_id "%s"; reference_transcript_id "%s"; \n' % \
                            (model.chr_id, model.exon_blocks[0][0], model.exon_blocks[-1][1], model.strand,
                             model.gene_id, model.transcript_id, model.reference_gene, model.reference_transcript)
                self.out_gff.write(transcript_line)

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
        MatchEventSubtype.intron_retention, MatchEventSubtype.unspliced_intron_retention,
        MatchEventSubtype.alt_left_site_known, MatchEventSubtype.alt_right_site_known,
        MatchEventSubtype.alt_left_site_novel, MatchEventSubtype.alt_right_site_novel,
        MatchEventSubtype.extra_intron_novel, MatchEventSubtype.extra_intron_known,
        MatchEventSubtype.extra_intron_flanking_left, MatchEventSubtype.extra_intron_flanking_right,
        MatchEventSubtype.intron_migration,
        MatchEventSubtype.intron_alternation_novel, MatchEventSubtype.intron_alternation_known,
        MatchEventSubtype.mutually_exclusive_exons_novel, MatchEventSubtype.mutually_exclusive_exons_known,
        MatchEventSubtype.exon_skipping_novel, MatchEventSubtype.exon_skipping_known,
        MatchEventSubtype.terminal_exon_shift_known, MatchEventSubtype.terminal_exon_shift_novel,
        MatchEventSubtype.exon_gain_novel, MatchEventSubtype.exon_gain_known,
        MatchEventSubtype.alternative_structure_known, MatchEventSubtype.alternative_structure_novel
    }

    modification_events =  {MatchEventSubtype.alt_left_site_novel, MatchEventSubtype.alt_right_site_novel,
                            MatchEventSubtype.extra_intron_novel, MatchEventSubtype.extra_intron_flanking_left,
                            MatchEventSubtype.extra_intron_flanking_right,
                            MatchEventSubtype.mutually_exclusive_exons_novel,
                            MatchEventSubtype.exon_gain_novel, MatchEventSubtype.exon_skipping_novel,
                            MatchEventSubtype.exon_detatch_novel, MatchEventSubtype.exon_merge_novel,
                            MatchEventSubtype.terminal_exon_shift_novel,
                            MatchEventSubtype.alternative_structure_novel, MatchEventSubtype.intron_alternation_novel,
                            MatchEventSubtype.unspliced_intron_retention, MatchEventSubtype.intron_retention,
                            MatchEventSubtype.alt_left_site_known, MatchEventSubtype.alt_right_site_known,
                            MatchEventSubtype.extra_intron_known, MatchEventSubtype.intron_migration,
                            MatchEventSubtype.mutually_exclusive_exons_known, MatchEventSubtype.exon_skipping_known,
                            MatchEventSubtype.exon_detatch_known, MatchEventSubtype.exon_merge_known,
                            MatchEventSubtype.terminal_exon_shift_known,
                            MatchEventSubtype.exon_gain_known, MatchEventSubtype.alternative_structure_known,
                            MatchEventSubtype.intron_alternation_known}

    def __init__(self, gene_info, read_assignment_storage, params):
        if params.report_apa:
            self.events_to_track.add(MatchEventSubtype.alternative_polya_site)
            self.events_to_track.add(MatchEventSubtype.alternative_tss)
        if params.report_intron_retention:
            self.events_to_track.add(MatchEventSubtype.intron_retention)
            self.events_to_track.add(MatchEventSubtype.unspliced_intron_retention)
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
                logger.debug("== Processing modification cluster for isoform %s of size %d, modifications:" %
                             (isoform_id, len(assignments)))
                logger.debug(", ".join(["%s: %s" % (x[0], str(x[1])) for x in modification]))
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
            #TODO only unique assignments for novel models?
            for match in read_assignment.isoform_matches:
                isoform_id = match.assigned_transcript
                if read_assignment.assignment_type in {ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference}:
                #if match.match_classification in {MatchClassification.full_splice_match,
                #                                  MatchClassification.incomplete_splice_match,
                #                                  MatchClassification.mono_exon_match}:
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
            logger.debug("> Representative read chosen: %s" % representative_read_assignment.read_id)
            logger.debug(representative_read_assignment.combined_profile.read_exon_profile.read_features)
            logger.debug(representative_read_assignment.combined_profile.read_intron_profile.read_features)
            # create a new transcript model

            self.representative_reads.add(representative_read_assignment.read_id)
            new_transcript_model = self.blend_read_into_isoform(isoform_id, representative_read_assignment)
            if not new_transcript_model:
                logger.debug("> No novel model was constructed")
                return
            logger.debug("Created new candidate transcript model %s : %s " %
                         (new_transcript_model.transcript_id, str(new_transcript_model.exon_blocks)))
            # compare read junctions with novel transcript model, count them and keep only those that do not match
            remaining_assignments = self.verify_novel_model(isoform_id, remaining_assignments, new_transcript_model,
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
            # logger.debug("%s %d %d" % (a.read_id, a.combined_profile.polya_info.external_polya_pos,
            #                            a.combined_profile.polya_info.external_polyt_pos))
            if strand == '+':
                if not self.params.require_polyA or a.combined_profile.polya_info.external_polya_pos != -1:
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
                if not self.params.require_polyA or a.combined_profile.polya_info.external_polyt_pos != -1:
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
        logger.debug("Creating novel transcript model for isoform %s and read %s" % (isoform_id, read_assignment.read_id))
        modification_events_map = self.derive_significant_modifications_map(isoform_id, read_assignment)
        if not modification_events_map:
            return None

        isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
        combined_profile = read_assignment.combined_profile
        read_introns = combined_profile.read_intron_profile.read_features
        read_start = combined_profile.corrected_read_start
        read_end = combined_profile.corrected_read_end
        read_region = (read_start, read_end)
        novel_exons = []

        logger.debug("Isoform I " + str(isoform_introns))
        logger.debug("Isoform E " + str(self.gene_info.all_isoforms_exons[isoform_id]))
        logger.debug("Read coords %d, %d" % (read_start, read_end))
        logger.debug("Read " + str(read_introns))

        if SupplementaryMatchConstansts.extra_left_mod_position in modification_events_map:
            # if there are extra introns on the left
            current_exon_start = read_start
            events = modification_events_map[SupplementaryMatchConstansts.extra_left_mod_position]
            current_exon_start = self.process_intron_related_events(events, None, isoform_introns, read_introns,
                                                                    read_region, novel_exons, current_exon_start)
        else:
            current_exon_start = read_start

        isoform_pos = 0
        logger.debug(str(modification_events_map))
        while isoform_pos <= len(isoform_introns):
            if isoform_pos not in modification_events_map.keys():
                if isoform_pos == len(isoform_introns):
                    # such position is possible only when extra intron is present inside last reference exon
                    break
                if isoform_introns[isoform_pos][0] < current_exon_start:
                    # skip introns that outside of gene region
                    isoform_pos += 1
                    continue
                if isoform_introns[isoform_pos][1] >= read_end:
                    # skip introns that ourside of gene region
                    break

                # simply select reference isoform intron
                logger.debug("Adding ref exon: %d, %d" % (isoform_pos, current_exon_start))
                current_exon_start = self.add_intron(novel_exons, current_exon_start, isoform_introns[isoform_pos])
                isoform_pos += 1

            else:
                current_events = modification_events_map[isoform_pos]
                current_exon_start = self.process_intron_related_events(current_events, isoform_pos, isoform_introns,
                                                                        read_introns, read_region,
                                                                        novel_exons, current_exon_start)
                if isoform_pos < len(isoform_introns) \
                        and current_exon_start < isoform_introns[isoform_pos][0] \
                        and isoform_introns[isoform_pos][1] < read_end:
                    # intron modification was processed but nothing overlapping was added =>
                    # extra intron within previous exon => add this intron as is
                    # check that is was really extra intron
                    extra_intron_types = {MatchEventSubtype.extra_intron_novel, MatchEventSubtype.extra_intron_known}
                    only_extra_intron = all(el.event_type in extra_intron_types for el in current_events)
                    if only_extra_intron:
                        current_exon_start = self.add_intron(novel_exons, current_exon_start, isoform_introns[isoform_pos])
                        logger.debug("Adding reference intron after additional extra intron: " + str(isoform_introns[isoform_pos]))

                isoform_pos += 1
                while isoform_pos < len(isoform_introns) and isoform_introns[isoform_pos][0] < current_exon_start:
                    isoform_pos += 1

        if SupplementaryMatchConstansts.extra_right_mod_position in modification_events_map:
            # if there are extra introns on the right
            events = modification_events_map[SupplementaryMatchConstansts.extra_right_mod_position]
            current_exon_start = self.process_intron_related_events(events, None, isoform_introns, read_introns,
                                                                    read_region, novel_exons, current_exon_start)
            novel_transcript_end = read_end
        else:
            novel_transcript_end = read_end

        novel_exons.append((current_exon_start, novel_transcript_end))
        novel_exons = self.correct_transcripts_ends(novel_exons, combined_profile, isoform_id, modification_events_map)

        if not self.validate_exons(novel_exons):
            logger.warning("Error in novel transcript, not sorted or incorrect exon coords")
            logger.warning(novel_exons)
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
        strand = self.gene_info.isoform_strands[isoform_id]

        if SupplementaryMatchConstansts.extra_left_mod_position not in modification_events_map and \
                0 not in modification_events_map:
            # change only if there are no extra introns on the left and first intron is not modified
            novel_transcript_start = novel_exons[0][0]
            known_isoform_start = self.gene_info.transcript_start(isoform_id)
            if (strand == "+" or combined_profile.polya_info.external_polyt_pos == -1) and \
                    abs(novel_transcript_start - known_isoform_start) <= self.params.max_dist_to_isoforms_tsts and \
                    known_isoform_start < novel_exons[0][1]:
                # correct model start only if no polyT is found
                novel_exons[0] = (known_isoform_start, novel_exons[0][1])

        last_index = len(self.gene_info.all_isoforms_introns[isoform_id]) - 1
        if SupplementaryMatchConstansts.extra_right_mod_position not in modification_events_map and \
                last_index not in modification_events_map:
            # change only if there are no extra introns on the right and last intron is not modified
            novel_transcript_end = novel_exons[-1][1]
            known_isoform_end = self.gene_info.transcript_end(isoform_id)
            if (strand == "-" or combined_profile.polya_info.external_polya_pos == -1) and \
                    abs(novel_transcript_end - known_isoform_end) <= self.params.max_dist_to_isoforms_tsts and \
                    known_isoform_end > novel_exons[-1][0]:
                # correct model end only if no polyA is found
                novel_exons[-1] = (novel_exons[-1][0], known_isoform_end)

        return novel_exons

    # process a sorted list of events assigned to the same intron
    def process_intron_related_events(self, sorted_event_list, isoform_pos, isoform_introns, read_introns, read_region,
                                      novel_exons, current_exon_start):
        logger.debug("> Processing events for position %s: %s" % (str(isoform_pos), str(sorted_event_list)))
        logger.debug("> Before: %d, %s" % (current_exon_start, novel_exons))
        for event in sorted_event_list:
            current_exon_start = self.process_single_event(event, isoform_pos, isoform_introns, read_introns,
                                                           read_region, novel_exons, current_exon_start)
            logger.debug("> In progress: %d, %s" % (current_exon_start, novel_exons))

        logger.debug("> After: %d, %s" % (current_exon_start, novel_exons))
        return current_exon_start

    # process single event
    def process_single_event(self, event_tuple, isoform_pos, isoform_introns, read_introns, read_region,
                             novel_exons, current_exon_start):
        logger.debug("> > Applying event %s at position %s" % (event_tuple.event_type.name, str(isoform_pos)))
        if event_tuple.event_type in {MatchEventSubtype.intron_retention, MatchEventSubtype.unspliced_intron_retention}:
            # simply skip reference intron
            return current_exon_start

        if event_tuple.read_region[0] == JunctionComparator.absent:
            logger.warning("Undefined read intron position for event type: %s" % event_tuple.event_type.name)
            return current_exon_start
        read_intron = read_introns[event_tuple.read_region[0]]
        if not contains(read_region, read_intron):
            logger.warning("Read intron to be added seems to be outside of read region: %s, read: %s, intron: %s" %
                           (event_tuple.event_type.name, str(read_region), str(read_intron)))
            return current_exon_start
        # logger.debug("Novel intron " + str(read_intron))

        if event_tuple.event_type in [MatchEventSubtype.extra_intron_novel,
                                      MatchEventSubtype.extra_intron_flanking_left,
                                      MatchEventSubtype.extra_intron_flanking_right]:
            return self.add_intron(novel_exons, current_exon_start, read_intron)
        elif event_tuple.event_type == MatchEventSubtype.extra_intron_known:
            corrected_intron = self.get_closest_ref_intron(read_intron)
            return self.add_intron(novel_exons, current_exon_start, corrected_intron)

        isoform_intron = isoform_introns[isoform_pos]
        assert overlaps(read_intron, isoform_intron)

        if event_tuple.event_type == MatchEventSubtype.alt_left_site_novel:
            novel_intron = (read_intron[0], isoform_intron[1])
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type == MatchEventSubtype.alt_left_site_known:
            novel_intron = self.get_closest_ref_intron((read_intron[0], isoform_intron[1]))
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type == MatchEventSubtype.alt_right_site_novel:
            novel_intron = (isoform_intron[0], read_intron[1])
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type == MatchEventSubtype.alt_right_site_known:
            novel_intron = self.get_closest_ref_intron((isoform_intron[0], read_intron[1]))
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type in {MatchEventSubtype.intron_alternation_novel,
                                        MatchEventSubtype.exon_skipping_novel,
                                        MatchEventSubtype.exon_merge_novel,
                                        MatchEventSubtype.terminal_exon_shift_novel}:
            # simply add read intron
            novel_intron = (read_intron[0], read_intron[1])
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type in {MatchEventSubtype.intron_alternation_known,
                                        MatchEventSubtype.intron_migration,
                                        MatchEventSubtype.exon_skipping_known,
                                        MatchEventSubtype.exon_merge_known,
                                        MatchEventSubtype.terminal_exon_shift_known}:
            # simply add corrected read intron
            novel_intron = self.get_closest_ref_intron((read_intron[0], read_intron[1]))
            current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        elif event_tuple.event_type in {MatchEventSubtype.mutually_exclusive_exons_novel,
                                        MatchEventSubtype.exon_gain_novel,
                                        MatchEventSubtype.exon_detatch_novel,
                                        MatchEventSubtype.alternative_structure_novel}:
            # simply insert several reads introns
            for read_pos in range(event_tuple.read_region[0], event_tuple.read_region[1] + 1):
                current_exon_start = self.add_intron(novel_exons, current_exon_start, read_introns[read_pos])
        elif event_tuple.event_type in {MatchEventSubtype.mutually_exclusive_exons_known,
                                        MatchEventSubtype.exon_gain_known,
                                        MatchEventSubtype.exon_detatch_known,
                                        MatchEventSubtype.alternative_structure_known}:
            # insert several reads introns my fitting them onto reference introns
            for read_pos in range(event_tuple.read_region[0], event_tuple.read_region[1] + 1):
                # TODO speed up
                novel_intron = self.get_closest_ref_intron(read_introns[read_pos])
                current_exon_start = self.add_intron(novel_exons, current_exon_start, novel_intron)
        else:
            logger.warning("Unsupported event type %s" % event_tuple.event_type.name)

        return current_exon_start

    # return list of all reads modifications relative to this isoform
    def get_read_inconsistencies(self, isoform_id, read_assignment):
        for match in read_assignment.isoform_matches:
            if match.assigned_transcript == isoform_id:
                return match.match_subclassifications
        return None

    # return map: isoform position -> event tuple
    def derive_significant_modifications_map(self, isoform_id, read_assignment):
        read_inconsistencies = self.get_read_inconsistencies(isoform_id, read_assignment)
        if read_inconsistencies is None:
            logger.debug("No modification events detected for " + read_assignment.read_id)
            return None

        match_subclassifications = list(filter(lambda m: m.event_type in self.modification_events, read_inconsistencies))
        logger.debug("Selected modifications: " +", ".join(["%s: %s - %s" % (x.event_type.name, str(x.isoform_position), str(x.read_region))
                                                            for x in match_subclassifications]))
        if not self.params.report_intron_retention and \
                all(m.event_type in {MatchEventSubtype.intron_retention, MatchEventSubtype.unspliced_intron_retention}
                    for m in match_subclassifications):
            return None

        modification_events_map = defaultdict(list)
        for x in match_subclassifications:
            modification_events_map[x.isoform_position].append(x)
        for isoform_position in modification_events_map.keys():
            if len(modification_events_map[isoform_position]) == 1:
                continue
            modification_events_map[isoform_position] = \
                sorted(modification_events_map[isoform_position], key=lambda x: x.read_region[1])
            logger.debug(modification_events_map[isoform_position])

        if not modification_events_map:
            logger.debug("No modification events detected for " + read_assignment.read_id)
            return None
        logger.debug("Sorted modifications: " + ", ".join([str(x) + " - " + str(modification_events_map[x])
                                                           for x in sorted(modification_events_map.keys())]))
        return modification_events_map

    def get_closest_ref_intron(self, read_intron):
        # TODO speed up - binray search or interval tree
        intron_profile = self.intron_profile_constructor.construct_profile_for_features([read_intron])
        matched_intron = intron_profile.gene_profile.index(1)
        return self.intron_profile_constructor.known_features[matched_intron]

    def add_intron(self, novel_exons, current_exon_start, intron):
        exon = (current_exon_start, intron[0] - 1)
        novel_exons.append(exon)
        return intron[1] + 1

    def verify_novel_model(self, isoform_id, read_assignments, transcript_model, original_read_id, candidate_model_storage):
        if len(transcript_model.exon_blocks) == 1:
            return self.verify_novel_monoexonic_model(read_assignments, transcript_model, original_read_id, candidate_model_storage)
        else:
            return self.verify_novel_spliced_model(isoform_id, read_assignments, transcript_model, original_read_id,
                                                   candidate_model_storage)

    def verify_novel_monoexonic_model(self, read_assignments, transcript_model, original_read_id, candidate_model_storage):
        # TODO verify ends using CAGE/polyA/Illumina
        model_exons = transcript_model.exon_blocks
        assert len(model_exons) == 1
        isoform_start = model_exons[0][0]
        isoform_end = model_exons[-1][1]
        strand = transcript_model.strand

        assigned_reads = []
        unassigned_reads = []
        for assignment in read_assignments:
            read_start = assignment.combined_profile.corrected_read_start
            read_end = assignment.combined_profile.corrected_read_end
            start_matches = abs(read_start - isoform_start) < self.params.max_dist_to_novel_tsts
            end_matches = abs(read_end - isoform_end) < self.params.max_dist_to_novel_tsts
            if start_matches and end_matches:
                assigned_reads.append(assignment)
            else:
                unassigned_reads.append(assignment)

        if len(assigned_reads) >= self.params.min_novel_supporting_reads:
            # to confirm we need at least min_novel_supporting_reads supporting reads
            logger.debug("Successfully confirmed %s" % transcript_model.transcript_id)
            logger.debug(str(transcript_model.exon_blocks))
            candidate_model_storage.append(transcript_model)
            for a in assigned_reads:
                self.transcript_read_ids[transcript_model.transcript_id].add(a.read_id)
            return unassigned_reads
        else:
            logger.debug("Transcript candidate %s looks unreliable" % transcript_model.transcript_id)
            all_except_original = list(filter(lambda x: x.read_id != original_read_id, read_assignments))
            return all_except_original

    def verify_novel_spliced_model(self, isoform_id, read_assignments, transcript_model, original_read_id, candidate_model_storage):
        logger.debug("Verifying transcript model %s with %d reads" % (transcript_model.transcript_id, len(read_assignments)))
        model_exons = transcript_model.exon_blocks
        isoform_start = model_exons[0][0]
        isoform_end = model_exons[-1][1]

        transcript_model_gene_info = GeneInfo.from_model(transcript_model, self.params.delta)
        model_introns = transcript_model_gene_info.all_isoforms_introns[transcript_model.transcript_id]
        assigner = LongReadAssigner(transcript_model_gene_info, self.params)
        profile_constructor = CombinedProfileConstructor(transcript_model_gene_info, self.params)

        assigned_reads = []
        fsm_match_count = 0
        unassigned_reads = []
        nearby_starts_count = 0
        nearby_ends_count = 0
        for assignment in read_assignments:
            combined_profile = assignment.combined_profile
            read_exons = combined_profile.read_exon_profile.read_features
            logger.debug("Checking read %s: %s" % (assignment.read_id, str(read_exons)))
            model_combined_profile = profile_constructor.construct_profiles(read_exons, combined_profile.polya_info, [])
            model_assignment = assigner.assign_to_isoform(assignment.read_id, model_combined_profile)
            # check that no serious contradiction occurs
            profile_matches = model_assignment.assignment_type in [ReadAssignmentType.unique,
                                                                   ReadAssignmentType.unique_minor_difference]

            if profile_matches:
                corrected_read_start = combined_profile.corrected_read_start
                corrected_read_end = combined_profile.corrected_read_end
                start_matches = abs(corrected_read_start - isoform_start) < self.params.max_dist_to_novel_tsts
                end_matches = abs(corrected_read_end - isoform_end) < self.params.max_dist_to_novel_tsts
                logger.debug("Profile matches, start: %d, end %d" % (corrected_read_start, corrected_read_end))

                if start_matches:
                    nearby_starts_count += 1
                if end_matches:
                    nearby_ends_count += 1
                # all read introns were mapped, read is assigned
                assigned_reads.append(assignment.read_id)
                # since profile is not reliable due to intron shifts etc
                # considering that there are no serious errors, covering all introns in enough
                is_fsm = contains((corrected_read_start, corrected_read_end), (model_introns[0][0], model_introns[-1][1]))
                if is_fsm:
                    logger.debug("Matches as FSM")
                    # all introns of novel model are covered
                    fsm_match_count += 1
            else:
                unassigned_reads.append(assignment)

        logger.debug("Stats for %s, FSM = %d, total = %d, start = %d, end = %d" %
                     (transcript_model.transcript_id, fsm_match_count, len(assigned_reads),
                      nearby_starts_count, nearby_ends_count))

        if fsm_match_count == 0:
            logger.warning("Zero FSM for transcript model %s" % transcript_model.transcript_id)

        if len(assigned_reads) >= self.params.min_novel_supporting_reads and \
                fsm_match_count >= self.params.min_novel_fsm_supporting_reads and \
                nearby_starts_count >= self.params.min_reads_supporting_tsts and \
                nearby_ends_count >= self.params.min_reads_supporting_tsts:
            # to confirm we need at least min_novel_supporting_reads supporting reads
            # and at least min_novel_fsm_supporting_reads FSM
            logger.debug("Successfully confirmed %s" % transcript_model.transcript_id)
            logger.debug(str(transcript_model.exon_blocks))
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

