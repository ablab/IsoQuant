############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
from functools import cmp_to_key

from .common import (
    AtomicCounter,
    cmp,
    get_exons,
    get_top_count,
    intersection_len,
    interval_len,
    junctions_from_blocks,
    read_coverage_fraction,
)
from .assignment_io import ReadAssignmentType
from .gene_info import GeneInfo, StrandDetector, TranscriptModel, TranscriptModelType
from .intron_graph import IntronGraph, VERTEX_polya, VERTEX_polyt, VERTEX_read_end, VERTEX_read_start
from .isoform_assignment import is_matching_assignment, match_subtype_to_str_with_additional_info, MatchEventSubtype
from .long_read_assigner import LongReadAssigner
from .long_read_profiles import CombinedProfileConstructor
from .polya_finder import PolyAInfo


logger = logging.getLogger('IsoQuant')


class GraphBasedModelConstructor:
    transcript_id_counter = AtomicCounter()
    transcript_prefix = "transcript"
    known_transcript_suffix = ".known"
    nic_transcript_suffix = ".nic"
    nnic_transcript_suffix = ".nnic"
    detected_known_isoforms = set()
    extended_transcript_ids = set()

    def __init__(self, gene_info, chr_record, params, transcript_counter):
        self.gene_info = gene_info
        self.chr_record = chr_record
        self.params = params

        self.strand_detector = StrandDetector(self.chr_record)
        self.intron_genes = defaultdict(set)
        self.set_gene_properties()

        self.intron_graph = None
        self.path_processor = None
        self.path_storage = None

        self.known_isoforms_in_graph = {}
        self.known_introns = set()
        self.known_isoforms_in_graph_ids = {}
        self.assigner = LongReadAssigner(self.gene_info, self.params)
        self.profile_constructor = CombinedProfileConstructor(self.gene_info, self.params)

        self.transcript_model_storage = []
        self.extended_annotation_storage = []
        self.transcript_read_ids = defaultdict(list)
        self.transcript_counter = transcript_counter
        self.internal_counter = defaultdict(int)
        self.reads_used_in_construction = set()
        self.unused_reads = []
        self.transcript2transcript = []

    def get_transcript_id(self):
        return GraphBasedModelConstructor.transcript_id_counter.increment()

    def set_gene_properties(self):
        intron_strands_dicts = defaultdict(lambda: defaultdict(int))
        self.intron_genes = defaultdict(set)
        for t_id, introns in self.gene_info.all_isoforms_introns.items():
            strand = self.gene_info.isoform_strands[t_id]
            gene_id = self.gene_info.gene_id_map[t_id]
            for intron in introns:
                intron_strands_dicts[intron][strand] += 1
                self.intron_genes[intron].add(gene_id)

        for intron in intron_strands_dicts.keys():
            if len(intron_strands_dicts[intron].keys()) == 1:
                # intron has a single strand
                self.strand_detector.set_strand(intron, list(intron_strands_dicts[intron].keys())[0])
            else:
                self.strand_detector.set_strand(intron)

    def select_reference_gene(self, transcript_introns, transcript_range, transcript_strand):
        if self.gene_info.empty():
            return None
        gene_counts = defaultdict(int)
        for intron in transcript_introns:
            if intron not in self.intron_genes:
                continue
            for g_id in self.intron_genes[intron]:
                gene_counts[g_id] += 1

        ordered_genes = sorted(gene_counts.items(), key=lambda x: (x[1], x[0]), reverse=True)
        for g in ordered_genes:
            gene_id = g[0]
            if transcript_strand == '.' or self.gene_info.gene_strands[gene_id] == transcript_strand:
                return gene_id

        overlap_dict = {}
        gene_regions = self.gene_info.get_gene_regions()
        for gene_id in gene_regions.keys():
            gene_coverage = read_coverage_fraction([transcript_range], [gene_regions[gene_id]])
            if gene_coverage > 0.0 and \
                    (transcript_strand == '.' or self.gene_info.gene_strands[gene_id] == transcript_strand):
                overlap_dict[gene_id] = gene_coverage

        if overlap_dict:
            return get_top_count(overlap_dict)
        return None

    def process(self, read_assignment_storage):
        self.intron_graph = IntronGraph(self.params, self.gene_info, read_assignment_storage)
        self.path_processor = IntronPathProcessor(self.params, self.intron_graph)
        self.path_storage = IntronPathStorage(self.params, self.path_processor)
        self.path_storage.fill(read_assignment_storage)
        self.known_isoforms_in_graph = self.get_known_spliced_isoforms(self.gene_info)
        self.known_introns = set(self.gene_info.intron_profiles.features)

        for intron_path, isoform_id in self.known_isoforms_in_graph.items():
            self.known_isoforms_in_graph_ids[isoform_id] = intron_path

        self.construct_fl_isoforms()
        self.construct_assignment_based_isoforms(read_assignment_storage)
        self.assign_reads_to_models(read_assignment_storage)
        self.filter_transcripts()

        if self.params.genedb:
            self.create_extended_annotation()

        if not self.gene_info.all_isoforms_exons:
            transcript_joiner = TranscriptToGeneJoiner(self.transcript_model_storage)
            self.transcript_model_storage = transcript_joiner.join_transcripts()

        if self.params.sqanti_output:
            self.compare_models_with_known()

    def compare_models_with_known(self):
        gene_to_model_dict = defaultdict(list)
        for model in self.transcript_model_storage:
            gene_to_model_dict[model.gene_id].append(model.transcript_id)

        self.transcript2transcript = []
        for model in self.transcript_model_storage:
            if model.transcript_type == TranscriptModelType.known:
                continue
            if model.strand == "-":
                polya_info = PolyAInfo(-1, model.exon_blocks[0][0], -1, -1)
            else:
                polya_info = PolyAInfo(model.exon_blocks[-1][1], -1, -1, -1)

            combined_profile = self.profile_constructor.construct_profiles(model.exon_blocks, polya_info, [])
            assignment = self.assigner.assign_to_isoform(model.transcript_id, combined_profile)
            if assignment is None:
                continue

            assignment.polya_info = polya_info
            assignment.cage_found = False
            assignment.exons = model.exon_blocks
            assignment.strand = model.strand
            assignment.chr_id = model.chr_id
            assignment.set_additional_info("indel_count", "NA")
            assignment.set_additional_info("junctions_with_indels", "NA")
            assignment.introns_match = all(e == 1 for e in combined_profile.read_intron_profile.read_profile)
            assignment.gene_info = self.gene_info

            if assignment.assignment_type in [ReadAssignmentType.intergenic, ReadAssignmentType.noninformative] or \
                    not assignment.isoform_matches:
                # create intergenic
                assignment.assignment_type = ReadAssignmentType.intergenic
                FSM_class = "C"
                assignment.set_additional_info("FSM_class", FSM_class)
                self.transcript2transcript.append(assignment)
                continue

            if len(gene_to_model_dict[assignment.isoform_matches[0].assigned_gene]) == 1:
                FSM_class = "A"
            else:
                FSM_class = "C"
            assignment.set_additional_info("FSM_class", FSM_class)

            assigned_transcript_id = assignment.isoform_matches[0].assigned_transcript
            if not assigned_transcript_id or assigned_transcript_id not in self.gene_info.all_isoforms_introns:
                continue

            reference_introns = self.gene_info.all_isoforms_introns[assigned_transcript_id]
            isoform_introns = junctions_from_blocks(model.exon_blocks)
            event_string = ",".join([match_subtype_to_str_with_additional_info(x, model.strand,
                                                                               isoform_introns, reference_introns)
                                     for x in assignment.isoform_matches[0].match_subclassifications])
            model.add_additional_attribute("similar_reference_id", assigned_transcript_id)
            model.add_additional_attribute("alternatives", event_string)
            self.transcript2transcript.append(assignment)

    def filter_transcripts(self):
        filtered_storage = []
        confirmed_transcipt_ids = set()
        to_substitute = self.detect_similar_isoforms()

        for model in self.transcript_model_storage:
            if model.transcript_type == TranscriptModelType.known:
                filtered_storage.append(model)
                confirmed_transcipt_ids.add(model.transcript_id)
                continue
            # check coverage
            component_coverage = self.intron_graph.get_max_component_coverage(model.intron_path)
            if component_coverage == 0 or len(model.intron_path) == 0 or \
                    (len(model.intron_path) == 1 and self.intron_graph.is_monointron(model.intron_path[0])):
                component_coverage = self.intron_graph.get_overlapping_component_max_coverage((model.get_start(), model.get_end()))
                novel_isoform_cutoff = max(self.params.min_novel_count,
                                           self.params.min_mono_count_rel * component_coverage)
            else:
                novel_isoform_cutoff = max(self.params.min_novel_count,
                                           self.params.min_novel_count_rel * component_coverage)

            if model.transcript_id in to_substitute:
                #logger.debug("Novel model %s has a similar isoform %s" % (model.transcript_id, to_substitute[model.transcript_id]))
                self.transcript_read_ids[to_substitute[model.transcript_id]] += self.transcript_read_ids[model.transcript_id]
                del self.transcript_read_ids[model.transcript_id]
                continue

            if self.internal_counter[model.transcript_id] < novel_isoform_cutoff:
                #logger.debug("Novel model %s has coverage %d < %.2f, component cov = %d" % (model.transcript_id,
                #                                                        self.internal_counter[model.transcript_id],
                #                                                        novel_isoform_cutoff, component_coverage))
                del self.transcript_read_ids[model.transcript_id]
                continue

            if len(model.exon_blocks) <= 2:
                mapq = self.mapping_quality(model)
                #logger.debug("Novel model %s has quality %.2f" % (model.transcript_id, mapq))
                if mapq < self.params.simple_models_mapq_cutoff:
                    #logger.debug("Novel model %s has poor quality" % model.transcript_id)
                    del self.transcript_read_ids[model.transcript_id]
                    continue

            # TODO: correct ends for known
            self.correct_novel_transcrip_ends(model, self.transcript_read_ids[model.transcript_id])
            filtered_storage.append(model)
            confirmed_transcipt_ids.add(model.transcript_id)

        self.transcript_model_storage = filtered_storage
        self.transcript_counter.add_confirmed_features(confirmed_transcipt_ids)

    def mapping_quality(self, model):
        mapq = 0
        for a in self.transcript_read_ids[model.transcript_id]:
            mapq += a.mapping_quality
        return mapq / len(self.transcript_read_ids[model.transcript_id])

    def detect_similar_isoforms(self):
        to_substitute = {}
        for model in self.transcript_model_storage:
            if len(model.exon_blocks) <= 2 or model.transcript_id in to_substitute:
                continue
            transcript_model_gene_info = GeneInfo.from_models([model], self.params.delta)
            assigner = LongReadAssigner(transcript_model_gene_info, self.params)
            profile_constructor = CombinedProfileConstructor(transcript_model_gene_info, self.params)

            for m in self.transcript_model_storage:
                if m.transcript_type == TranscriptModelType.known or m.transcript_id == model.transcript_id or \
                        m.transcript_id in to_substitute or len(m.exon_blocks) == 1 or not m.intron_path or \
                        len(m.exon_blocks) > len(model.exon_blocks):
                    continue

                if m.intron_path[0][0] == VERTEX_polyt:
                    polya_info = PolyAInfo(-1, m.intron_path[0][1], -1, -1)
                elif m.intron_path[-1][0] == VERTEX_polya:
                    polya_info = PolyAInfo(m.intron_path[-1][1], -1, -1, -1)
                else:
                    polya_info = PolyAInfo(-1, -1, -1, -1)
                combined_profile = profile_constructor.construct_profiles(m.exon_blocks, polya_info, [])
                assignment = assigner.assign_to_isoform(m.transcript_id, combined_profile)

                if is_matching_assignment(assignment):
                    to_substitute[m.transcript_id] = model.transcript_id

        return to_substitute

    def create_extended_annotation(self):
        self.extended_annotation_storage = []
        for isoform_id in self.gene_info.all_isoforms_exons.keys():
            if isoform_id not in GraphBasedModelConstructor.extended_transcript_ids:
                self.extended_annotation_storage.append(self.transcript_from_reference(isoform_id))
                GraphBasedModelConstructor.extended_transcript_ids.add(isoform_id)
        for model in self.transcript_model_storage:
            if model.transcript_type != TranscriptModelType.known:
                self.extended_annotation_storage.append(model)

    def get_known_spliced_isoforms(self, gene_info, s="known"):
        known_isoforms = {}
        for isoform_id in gene_info.all_isoforms_introns:
            isoform_introns = gene_info.all_isoforms_introns[isoform_id]
            intron_path = self.path_processor.thread_introns(isoform_introns)
            if not intron_path:
                #logger.debug("== No path found for %s isoform %s: %s" % (s, isoform_id, gene_info.all_isoforms_introns[isoform_id]))
                continue
            #logger.debug("== Path found for %s isoform %s: %s" % (s, isoform_id, gene_info.all_isoforms_introns[isoform_id]))
            known_isoforms[tuple(intron_path)] = isoform_id
        return known_isoforms

    def save_assigned_read(self, read_assignment, transcript_id):
        read_id = read_assignment.read_id
        self.transcript_read_ids[transcript_id].append(read_assignment)
        self.transcript_counter.add_read_info_raw(read_id, [transcript_id], read_assignment.read_group)
        self.internal_counter[transcript_id] += 1

    def construct_fl_isoforms(self):
        # a minor trick to compare tuples of pairs, whose starting and terminating elements have different type
        logger.debug("Total FL paths %d" % len(self.path_storage.fl_paths))
        for path in sorted(self.path_storage.fl_paths,
                           key=cmp_to_key(lambda x,y: cmp(x,y) if len(x)==len(y) else cmp(len(y), len(x)))):
            # do not include terminal vertices
            # logger.debug(">>> Considering path " + str(path))
            intron_path = path[1:-1]
            if not intron_path: continue
            transcript_range = (path[0][1], path[-1][1])
            novel_exons = get_exons(transcript_range, list(intron_path))
            count = self.path_storage.paths[path]
            new_transcript_id = self.transcript_prefix + str(self.get_transcript_id())
            # logger.debug("uuu %s: %s" % (new_transcript_id, str(novel_exons)))

            reference_isoform = None
            # check if new transcript matches a reference one
            if intron_path[0][0] == VERTEX_polyt:
                polya_info = PolyAInfo(-1, intron_path[0][1], -1, -1)
            elif intron_path[-1][0] == VERTEX_polya:
                polya_info = PolyAInfo(intron_path[-1][1], -1, -1, -1)
            else:
                polya_info = PolyAInfo(-1, -1, -1, -1)
            combined_profile = self.profile_constructor.construct_profiles(novel_exons, polya_info, [])
            assignment = self.assigner.assign_to_isoform(new_transcript_id, combined_profile)
            # check that no serious contradiction occurs
            # logger.debug("uuu Checking novel transcript %s: %s; assignment type %s" %
            #             (new_transcript_id, str(novel_exons), str(assignment.assignment_type)))

            if is_matching_assignment(assignment):
                reference_isoform = assignment.isoform_matches[0].assigned_transcript
                # logger.debug("uuu Substituting with known isoform %s" % reference_isoform)
            elif intron_path in self.known_isoforms_in_graph:
                # path was not assigned to any known isoform but intron chain still matches
                continue

            new_model = None
            if reference_isoform:
                # adding FL reference isoform
                if reference_isoform in GraphBasedModelConstructor.detected_known_isoforms:
                    pass
                elif count < self.params.min_known_count:
                    pass # logger.debug("uuu Isoform %s has low coverage %d" % (reference_isoform, count))
                else:
                    new_model = self.transcript_from_reference(reference_isoform)
                    GraphBasedModelConstructor.detected_known_isoforms.add(reference_isoform)
                    #logger.debug("Adding known spliced isoform %s" % reference_isoform)
                    #logger.debug("Annotated positions: %d, %d, %s" % (new_model.exon_blocks[0][0], new_model.exon_blocks[-1][1], new_model.strand))
                    #logger.debug("Graph positions: %s, %s" % (str(path[0]), str(path[-1])))
            else:
                # adding FL novel isoform
                # component_coverage = self.intron_graph.get_max_component_coverage(intron_path)
                novel_isoform_cutoff = self.params.min_novel_count

                has_polyt = path[0][0] == VERTEX_polyt
                has_polya = path[-1][0] == VERTEX_polya
                polya_site = has_polya or has_polyt
                transcript_strand = self.strand_detector.get_strand(intron_path, has_polya, has_polyt)
                transcript_ss_strand = self.strand_detector.get_strand(intron_path)

                #logger.debug("uuu Novel isoform %s has coverage: %d cutoff = %d, component cov = %d, max_coverage = %d"
                #             % (new_transcript_id, count, novel_isoform_cutoff, component_coverage, self.intron_graph.max_coverage))
                if count < novel_isoform_cutoff:
                    pass #logger.debug("uuu Novel isoform %s has low coverage: %d\t%d" % (new_transcript_id, count, novel_isoform_cutoff))
                elif len(novel_exons) == 2 and (not polya_site or transcript_ss_strand == '.'):
                    pass #logger.debug("uuu Avoiding single intron %s isoform: %d\t%s" % (new_transcript_id, count, str(path)))
                else:
                    if self.params.use_technical_replicas and \
                            len(set([a.read_group for a in self.path_storage.paths_to_reads[path]])) <= 1:
                        #logger.debug("%s was suspended due to technical replicas check" % new_transcript_id)
                        continue

                    transcript_gene = self.select_reference_gene(intron_path, transcript_range, transcript_strand)
                    if transcript_gene is None:
                        transcript_gene = "novel_gene_" + self.gene_info.chr_id + "_" + str(self.get_transcript_id())
                    elif transcript_strand == '.':
                        transcript_strand = self.gene_info.gene_strands[transcript_gene]

                    if all(intron in self.known_introns for intron in intron_path):
                        transcript_type = TranscriptModelType.novel_in_catalog
                        id_suffix = self.nic_transcript_suffix
                    else:
                        transcript_type = TranscriptModelType.novel_not_in_catalog
                        id_suffix = self.nnic_transcript_suffix

                    new_model = TranscriptModel(self.gene_info.chr_id, transcript_strand,
                                                new_transcript_id + ".%s" % self.gene_info.chr_id + id_suffix,
                                                transcript_gene, novel_exons, transcript_type)
                    new_model.intron_path = intron_path
                    logger.debug("uuu Adding novel spliced isoform %s : %d\t%d" % (new_transcript_id, count, novel_isoform_cutoff))

            if new_model:
                self.transcript_model_storage.append(new_model)
                for read_assignment in self.path_storage.paths_to_reads[path]:
                    self.save_assigned_read(read_assignment, new_model.transcript_id)
                    self.reads_used_in_construction.add(read_assignment.read_id)

    def construct_assignment_based_isoforms(self, read_assignment_storage):
        spliced_isoform_reads = defaultdict(list)
        isoform_left_support = defaultdict(int)
        isoform_right_support = defaultdict(int)
        polya_sites = defaultdict(int)
        mono_exon_isoform_reads = defaultdict(list)
        mono_exon_isoform_coverage = {}
        novel_mono_exon_reads = []

        for read_assignment in read_assignment_storage:
            if len(read_assignment.corrected_exons) <= 2 and \
                    (read_assignment.multimapper or read_assignment.mapping_quality < self.params.mono_mapping_quality_cutoff):
                continue

            if not read_assignment:
                continue
            if len(read_assignment.corrected_exons) == 1 and read_assignment.polyA_found and not read_assignment.multimapper and \
                    read_assignment.assignment_type in {ReadAssignmentType.noninformative, ReadAssignmentType.inconsistent, ReadAssignmentType.intergenic}:
                novel_mono_exon_reads.append(read_assignment)

            if read_assignment.assignment_type not in {ReadAssignmentType.unique,
                                                       ReadAssignmentType.unique_minor_difference}:
                continue
            refrenence_isoform_id = read_assignment.isoform_matches[0].assigned_transcript
            if refrenence_isoform_id in GraphBasedModelConstructor.detected_known_isoforms:
                continue

            events = read_assignment.isoform_matches[0].match_subclassifications
            if any(e.event_type == MatchEventSubtype.mono_exon_match for e in events):
                mono_exon_isoform_reads[refrenence_isoform_id].append(read_assignment)
                assert len(self.gene_info.all_isoforms_introns[refrenence_isoform_id]) == 0
                transcript_start = self.gene_info.all_isoforms_exons[refrenence_isoform_id][0][0]
                transcript_end = self.gene_info.all_isoforms_exons[refrenence_isoform_id][-1][1]
                t_len = transcript_end - transcript_start + 1

                if refrenence_isoform_id not in mono_exon_isoform_coverage:
                    mono_exon_isoform_coverage[refrenence_isoform_id] = [0 for _ in range(t_len)]
                start = max(0, read_assignment.corrected_exons[0][0] - transcript_start)
                end = min(t_len, read_assignment.corrected_exons[-1][1] - transcript_start + 1)
                for i in range(start, end):
                    mono_exon_isoform_coverage[refrenence_isoform_id][i] = 1

                if self.gene_info.isoform_strands[refrenence_isoform_id] == '+':
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_right for x in events):
                        polya_sites[refrenence_isoform_id] += 1
                else:
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_left for x in events):
                        polya_sites[refrenence_isoform_id] += 1
            elif len(self.gene_info.all_isoforms_exons[refrenence_isoform_id]) > 1:
                if read_assignment.read_id in self.reads_used_in_construction:
                    pass
                    # logger.debug("Spliced read %s was previously used for construction, assigned id %s" %
                    #             (read_assignment.read_id, refrenence_isoform_id))
                spliced_isoform_reads[refrenence_isoform_id].append(read_assignment)

                if self.params.needs_polya_for_construction and self.gene_info.isoform_strands[refrenence_isoform_id] == '-':
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_left for x in events):
                        isoform_left_support[refrenence_isoform_id] += 1
                elif abs(self.gene_info.all_isoforms_exons[refrenence_isoform_id][0][0] - read_assignment.corrected_exons[0][0]) <= self.params.apa_delta:
                    isoform_left_support[refrenence_isoform_id] += 1

                if self.params.needs_polya_for_construction and self.gene_info.isoform_strands[refrenence_isoform_id] == '+':
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_right for x in events):
                        isoform_right_support[refrenence_isoform_id] += 1
                elif abs(self.gene_info.all_isoforms_exons[refrenence_isoform_id][-1][1] - read_assignment.corrected_exons[-1][1]) <= self.params.apa_delta:
                    isoform_right_support[refrenence_isoform_id] += 1

        self.construct_monoexon_isoforms(mono_exon_isoform_reads, mono_exon_isoform_coverage, polya_sites)
        if not self.params.fl_only:
            logger.debug("Constructing nonFL isoforms")
            self.construct_nonfl_isoforms(spliced_isoform_reads, isoform_left_support, isoform_right_support)
        if self.params.report_novel_unspliced:
            self.construct_monoexon_novel(novel_mono_exon_reads)

    def collect_terminal_exons_from_graph(self):
        polya_exons = []
        polyt_exons = []
        for intron in self.intron_graph.outgoing_edges.keys():
            for v in self.intron_graph.outgoing_edges[intron]:
                if v[0] == VERTEX_polya:
                    polya_exons.append((intron[1], v[1]))
        for intron in self.intron_graph.incoming_edges.keys():
            for v in self.intron_graph.incoming_edges[intron]:
                if v[0] == VERTEX_polyt:
                    polyt_exons.append((v[1], intron[0]))
        # logger.debug("PolyA terminal exons: " + str(polya_exons))
        # logger.debug("PolyT terminal exons: " + str(polyt_exons))
        return polya_exons, polyt_exons

    def is_internal_monoexonic_read(self, alignment, terminal_exons, forward=True):
        read_coordinates = alignment.corrected_exons[0]
        if forward:
            for e in terminal_exons:
                if abs(e[1] - alignment.corrected_exons[-1][1]) <= self.params.apa_delta and \
                        read_coordinates[0] >= e[0] - self.params.delta:
                    return True
        else:
            for e in terminal_exons:
                if abs(e[0] - alignment.corrected_exons[0][0]) <= self.params.apa_delta and \
                        read_coordinates[1] <= e[1] + self.params.delta:
                    return True
        return False

    def construct_monoexon_novel(self, novel_mono_exon_reads):
        logger.debug("Constructing novel monoexon")
        polya_exons, polyt_exons = self.collect_terminal_exons_from_graph()
        polya_reads = defaultdict(list)
        polyt_reads = defaultdict(list)
        for a in novel_mono_exon_reads:
            if a.polya_info.external_polya_pos != -1:
                if not self.is_internal_monoexonic_read(a, polya_exons, forward=True):
                    polya_reads[a.polya_info.external_polya_pos].append(a)
            if a.polya_info.external_polyt_pos != -1:
                if not self.is_internal_monoexonic_read(a, polyt_exons, forward=False):
                    polyt_reads[a.polya_info.external_polyt_pos].append(a)

        novel_monoexon = set()
        clustered_polya_reads = self.cluster_monoexons(polya_reads)
        novel_monoexon.update(self.generate_monoexon_from_clustered(clustered_polya_reads, True))
        clustered_polyt_reads = self.cluster_monoexons(polyt_reads)
        novel_monoexon.update(self.generate_monoexon_from_clustered(clustered_polyt_reads, False))

    def generate_monoexon_from_clustered(self, clustered_reads, forward=True):
        cutoff = self.params.min_novel_count
        result = set()
        for three_prime_pos in clustered_reads.keys():
            count = len(clustered_reads[three_prime_pos])
            if count < cutoff:
                continue
            # TODO improve
            if forward:
                five_prime_pos = min([a.corrected_exons[0][0] for a in clustered_reads[three_prime_pos]])
            else:
                five_prime_pos = max([a.corrected_exons[-1][1] for a in clustered_reads[three_prime_pos]])

            new_transcript_id = self.transcript_prefix + str(self.get_transcript_id())
            transcript_gene = "novel_gene_" + self.gene_info.chr_id + "_" + str(self.get_transcript_id())
            transcript_type = TranscriptModelType.novel_not_in_catalog
            id_suffix = self.nnic_transcript_suffix
            strand = '+' if forward else '-'
            coordinates = (five_prime_pos, three_prime_pos) if forward else (three_prime_pos, five_prime_pos)

            is_valid = True
            half_len = interval_len(coordinates) / 2
            for existing_model in self.transcript_model_storage:
                if any(intersection_len(exon, coordinates) > half_len for exon in existing_model.exon_blocks):
                    is_valid = False
                    break
            if not is_valid:
                continue

            new_model = TranscriptModel(self.gene_info.chr_id, strand,
                                        new_transcript_id + ".%s" % self.gene_info.chr_id + id_suffix,
                                        transcript_gene, [coordinates], transcript_type)
            # logger.debug("uuu Adding novel MONOEXON isoform %s : %s, %d\t%d" % (new_transcript_id, str(coordinates), count, cutoff))
            result.add(coordinates)

            self.transcript_model_storage.append(new_model)
            for read_assignment in clustered_reads[three_prime_pos]:
                self.save_assigned_read(read_assignment, new_model.transcript_id)
                self.reads_used_in_construction.add(read_assignment.read_id)
        return result

    def cluster_monoexons(self, grouped_reads):
        clustered_counts = defaultdict(list)
        while grouped_reads:
            best_pair = max(grouped_reads.items(), key=lambda x:len(x[1]))
            top_position = best_pair[0]
            for pos in range(top_position - self.params.apa_delta, top_position + self.params.apa_delta + 1):
                if pos in grouped_reads:
                    clustered_counts[top_position] += grouped_reads[pos]
                    del grouped_reads[pos]

        return clustered_counts

    def construct_monoexon_isoforms(self, mono_exon_isoform_reads, mono_exon_isoform_coverage, polya_sites):
        for isoform_id in mono_exon_isoform_reads.keys():
            count = len(mono_exon_isoform_reads[isoform_id])
            coverage = float(mono_exon_isoform_coverage[isoform_id].count(1)) / \
                       float(len(mono_exon_isoform_coverage[isoform_id]))
            polya_support = polya_sites[isoform_id]

            # logger.debug(">> Monoexon transcript %s: %d\t%d\t%.4f\t%d" % (isoform_id, self.intron_graph.max_coverage, count, coverage, polya_support))
            if count < self.params.min_known_count or coverage < self.params.min_mono_exon_coverage or polya_support == 0:
                pass # logger.debug(">> Will NOT be added, abs cutoff=%d" % (self.params.min_known_count))
            elif isoform_id not in GraphBasedModelConstructor.detected_known_isoforms:
                new_model = self.transcript_from_reference(isoform_id)
                self.transcript_model_storage.append(new_model)
                GraphBasedModelConstructor.detected_known_isoforms.add(isoform_id)
                for read_assignment in mono_exon_isoform_reads[isoform_id]:
                    self.save_assigned_read(read_assignment, new_model.transcript_id)
                    self.reads_used_in_construction.add(read_assignment.read_id)
                # logger.debug(">> Adding known monoexon isoform %s, %s, count = %d: %s" %
                #             (self.transcript_model_storage[-1].transcript_id, isoform_id,
                #              count, str(self.gene_info.all_isoforms_exons[isoform_id])))

    def construct_nonfl_isoforms(self, spliced_isoform_reads, spliced_isoform_left_support, spliced_isoform_right_support):
        logger.debug("Constructing nonFL isoforms")
        for isoform_id in spliced_isoform_reads.keys():
            if isoform_id in GraphBasedModelConstructor.detected_known_isoforms:
                continue
            count = len(spliced_isoform_reads[isoform_id])
            if isoform_id not in self.known_isoforms_in_graph_ids:
                #logger.debug("<< Isoform %s has %d assignments but is not in the graph" % (isoform_id, count))
                continue

            intron_path = self.known_isoforms_in_graph_ids[isoform_id]
            #logger.debug("Known non-FL spliced isoform %s" % isoform_id)
            if count < self.params.min_known_count or \
                    spliced_isoform_left_support[isoform_id] < 1 or \
                    spliced_isoform_right_support[isoform_id] < 1:
                pass
                # logger.debug("Will not be added")
            else:
                logger.debug("<< Adding known non-FL spliced isoform %s" % isoform_id)
                new_model = self.transcript_from_reference(isoform_id)
                self.transcript_model_storage.append(new_model)
                GraphBasedModelConstructor.detected_known_isoforms.add(isoform_id)
                for read_assignment in spliced_isoform_reads[isoform_id]:
                    self.save_assigned_read(read_assignment, new_model.transcript_id)
                    self.reads_used_in_construction.add(read_assignment.read_id)

    # create transcript model object from reference isoforms
    def transcript_from_reference(self, isoform_id):
        gene_id = self.gene_info.gene_id_map[isoform_id]
        return TranscriptModel(self.gene_info.chr_id, self.gene_info.isoform_strands[isoform_id],
                               isoform_id, gene_id, self.gene_info.all_isoforms_exons[isoform_id],
                               TranscriptModelType.known)

    # assign reads back to constructed isoforms
    def assign_reads_to_models(self, read_assignments):
        if not self.transcript_model_storage:
            logger.debug("No transcripts were assigned")
            self.unused_reads = [a.read_id for a in read_assignments]
            return

        logger.debug("Creating artificial GeneInfo from %d transcript models" % len(self.transcript_model_storage))
        transcript_model_gene_info = GeneInfo.from_models(self.transcript_model_storage, self.params.delta)
        assigner = LongReadAssigner(transcript_model_gene_info, self.params, quick_mode=True)
        profile_constructor = CombinedProfileConstructor(transcript_model_gene_info, self.params)

        for assignment in read_assignments:
            read_id = assignment.read_id
            if read_id in self.reads_used_in_construction:
                # logger.debug("# Read %s was assigned to based on path construction" % read_id)
                continue

            read_exons = assignment.corrected_exons
            #logger.debug("# Checking read %s: %s" % (assignment.read_id, str(read_exons)))
            model_combined_profile = profile_constructor.construct_profiles(read_exons, assignment.polya_info, [])
            model_assignment = assigner.assign_to_isoform(assignment.read_id, model_combined_profile)
            model_assignment.read_group = assignment.read_group
            # check that no serious contradiction occurs
            if model_assignment.assignment_type in [ReadAssignmentType.unique,
                                                    ReadAssignmentType.unique_minor_difference,
                                                    ReadAssignmentType.ambiguous]:
                matched_isoforms = [m.assigned_transcript for m in model_assignment.isoform_matches]
                self.transcript_counter.add_read_info_raw(read_id,
                                                          matched_isoforms,
                                                          model_assignment.read_group)

                if len(matched_isoforms) == 1:
                    self.internal_counter[matched_isoforms[0]] += 1
                for m in model_assignment.isoform_matches:
                    self.transcript_read_ids[m.assigned_transcript].append(assignment)
            else:
                self.unused_reads.append(read_id)

    def correct_novel_transcrip_ends(self, transcript_model, assigned_reads):
        logger.debug("Verifying ends for transcript %s" % transcript_model.transcript_id)
        transcript_end = transcript_model.exon_blocks[-1][1]
        transcript_start = transcript_model.exon_blocks[0][0]
        start_supported = False
        read_starts = set()
        end_supported = False
        read_ends = defaultdict(int)

        for assignment in assigned_reads:
            read_exons = assignment.corrected_exons
            if abs(read_exons[0][0] - transcript_start) <= self.params.apa_delta:
                start_supported = True
            if not start_supported and read_exons[0][0] < transcript_model.exon_blocks[0][1]:
                read_starts.add(read_exons[0][0])
            if abs(read_exons[-1][1] - transcript_end) <= self.params.apa_delta:
                end_supported = True
            if not end_supported and read_exons[-1][1] > transcript_model.exon_blocks[-1][0]:
                read_ends[read_exons[-1][1]] += 1

        new_transcript_start = None
        if not start_supported:
            read_starts = sorted(read_starts)
            for read_start in read_starts:
                if read_start > transcript_start:
                    new_transcript_start = read_start
                    break
        if new_transcript_start and new_transcript_start < transcript_model.exon_blocks[0][1]:
            logger.debug("Changed start for transcript %s: from %d to %d" %
                         (transcript_model.transcript_id, transcript_model.exon_blocks[0][0], new_transcript_start))
            transcript_model.exon_blocks[0] = (new_transcript_start, transcript_model.exon_blocks[0][1])

        new_transcript_end = None
        if not end_supported:
            read_ends = sorted(read_ends, reverse=True)
            for read_end in read_ends:
                if read_end < transcript_end:
                    new_transcript_end = read_end
                    break
        if new_transcript_end and new_transcript_end > transcript_model.exon_blocks[-1][0]:
            logger.debug("Changed end for transcript %s: from %d to %d" %
                         (transcript_model.transcript_id, transcript_model.exon_blocks[-1][1], new_transcript_end))
            transcript_model.exon_blocks[-1] = (transcript_model.exon_blocks[-1][0], new_transcript_end)


class IntronPathStorage:
    def __init__(self, params, path_processor):
        self.params = params
        self.path_processor = path_processor
        self.intron_graph = path_processor.intron_graph
        self.paths = defaultdict(int)
        self.fl_paths = set()
        self.paths_to_reads = defaultdict(list)

    def fill(self, read_assignments):
        for a in read_assignments:
            if a.multimapper:
                continue
            intron_path = self.path_processor.thread_introns(a.corrected_introns)
            if not intron_path:
                continue
            read_end = a.corrected_exons[-1][1]
            is_end_trusted = a.strand == '+' and \
                             (a.polya_info.external_polya_pos != -1 or
                              a.polya_info.internal_polya_pos != -1)
            terminal_vertex = self.path_processor.thread_ends(intron_path[-1], read_end, is_end_trusted)
            if terminal_vertex:
                intron_path.append(terminal_vertex)

            read_start = a.corrected_exons[0][0]
            is_start_trusted = a.strand == '-' and \
                               (a.polya_info.external_polyt_pos != -1 or
                                a.polya_info.internal_polyt_pos != -1)
            starting_vertex = self.path_processor.thread_starts(intron_path[0], read_start, is_start_trusted)
            if starting_vertex:
                intron_path = [starting_vertex] + intron_path

            path_tuple = tuple(intron_path)
            self.paths[path_tuple] += 1
            if terminal_vertex and starting_vertex:
                if not self.params.needs_polya_for_construction or\
                        (terminal_vertex[0] == VERTEX_polya or starting_vertex[0] == VERTEX_polyt):
                    self.fl_paths.add(path_tuple)
            self.paths_to_reads[path_tuple].append(a)


class IntronPathProcessor:
    def __init__(self, params, intron_graph):
        self.params = params
        self.intron_graph = intron_graph
        self.all_vertices = set()
        self.all_vertices.update(self.intron_graph.intron_collector.clustered_introns.keys())
        for edge_set in self.intron_graph.outgoing_edges.values():
            self.all_vertices.update(edge_set)
        for edge_set in self.intron_graph.incoming_edges.values():
            self.all_vertices.update(edge_set)
        self.visited = set()

    def visit_vertex(self, v):
        if v in self.all_vertices:
            self.visited.add(v)

    def visit_path(self, p):
        for v in p:
            self.visit_vertex(v)

    def thread_introns(self, introns):
        path = []
        for intron in introns:
            if intron in self.intron_graph.intron_collector.discarded_introns:
                return None
            path.append(self.intron_graph.intron_collector.substitute(intron))
        return path

    def thread_ends(self, intron, end, trusted=False):
        possible_polyas = self.intron_graph.get_outgoing(intron, VERTEX_polya)
        if trusted:
            # find closes polyA
            for v in possible_polyas:
                if abs(v[1] - end) <= self.params.apa_delta:
                    return v

        outgoing_introns = self.intron_graph.get_outgoing(intron)
        if len(outgoing_introns) > 0:
            # intron has outgoing edges
            rightmost_exon_end = max([intron[0] for intron in outgoing_introns]) - 1
            if not trusted and end <= rightmost_exon_end + self.params.delta:
                # read end lies within next exon and has no polyA
                return None

        # consider all terminal position available for intron
        all_possible_ends = sorted(list(self.intron_graph.get_outgoing(intron, VERTEX_read_end)) +
                                   list(possible_polyas), key=lambda x:x[1])
        if len(all_possible_ends) == 0:
            return None

        rightmost_end = all_possible_ends[-1]
        if trusted and end >= rightmost_end[1] and rightmost_end[0] == VERTEX_read_end:
            # if we have trusted read, in cannot stop earlier that rightmost end (otherwise it should match polyA)
            return rightmost_end
        elif not trusted and end <= rightmost_end[1] + self.params.apa_delta and \
                (len(all_possible_ends) <= 1 or end > all_possible_ends[-2][1]):
            # non trusted should end before rightmost position + apa_delta but not earlier than second last
            return rightmost_end
        return None

    def thread_starts(self, intron, start, trusted=False):
        possible_polyas = self.intron_graph.get_incoming(intron, VERTEX_polyt)
        if trusted:
            # find closes polyT
            for v in possible_polyas:
                if abs(v[1] - start) <= self.params.apa_delta:
                    return v

        incoming_introns = self.intron_graph.get_incoming(intron)
        if len(incoming_introns) > 0:
            # intron has outgoing edges
            leftmost_exon_start = min([intron[1] for intron in incoming_introns]) + 1
            if not trusted and start >= leftmost_exon_start - self.params.delta:
                # read start lies within previous exon and has no polyA
                return None

        all_possible_starts = sorted(list(self.intron_graph.get_incoming(intron, VERTEX_read_start)) +
                                     list(possible_polyas), key=lambda x: x[1])
        if len(all_possible_starts) == 0:
            return None

        leftmost_start = all_possible_starts[0]
        if trusted and start <= leftmost_start[1] and leftmost_start[0] == VERTEX_read_start:
            return leftmost_start
        elif not trusted and start >= leftmost_start[1] and \
                (len(all_possible_starts) <= 1 or start < all_possible_starts[1][1]):
            return leftmost_start
        return None


class TranscriptToGeneJoiner:
    def __init__(self, transcipt_model_storage):
        self.transcipt_model_storage = transcipt_model_storage
        self.gene_introns = {}
        self.gene_strands = {}
        self.gene_regions = {}
        self.gene_to_transcripts = {}
        for t in self.transcipt_model_storage:
            self.gene_regions[t.gene_id] = (t.get_start(), t.get_end())
            self.gene_introns[t.gene_id] = set(junctions_from_blocks(t.exon_blocks))
            self.gene_strands[t.gene_id] = t.strand
            self.gene_to_transcripts[t.gene_id] = {t.transcript_id}
        self.scores = {}

    def count_score(self, gene1, gene2):
        logger.debug("Counting score %s %s" % (gene2, gene1))
        if self.gene_strands[gene1] != self.gene_strands[gene2]:
            return 0.0
        intronic_overlap = len(self.gene_introns[gene1].intersection(self.gene_introns[gene2])) / \
                           max(1, len(self.gene_introns[gene1].union(self.gene_introns[gene2])))
        range1 = self.gene_regions[gene1]
        range2 = self.gene_regions[gene2]
        position_overlap = intersection_len(range1, range2) / (max(range1[1], range2[1]) - min(range1[0], range2[0]) + 1)
        return position_overlap + intronic_overlap

    def count_scores(self):
        for g1_id in self.gene_to_transcripts.keys():
            for g2_id in self.gene_to_transcripts.keys():
                if g1_id == g2_id:
                    continue
                gene_pair = tuple(sorted([g1_id, g2_id]))
                if gene_pair not in self.scores:
                    self.scores[gene_pair] = self.count_score(g1_id, g2_id)

    def merge_genes(self, gene1, gene2):
        logger.debug("Merging %s into %s" % (gene2, gene1))
        range1 = self.gene_regions[gene1]
        range2 = self.gene_regions[gene2]
        self.gene_regions[gene1] = (min(range1[0], range2[0]), max(range1[1], range2[1]))
        self.gene_introns[gene1].update(self.gene_introns[gene2])
        self.gene_to_transcripts[gene1].update(self.gene_to_transcripts[gene2])
        if self.gene_strands[gene2] != self.gene_strands[gene1]:
            logger.error("Merging genes with different strands: %s, %s" % (gene1, gene2))
        del self.gene_regions[gene2]
        del self.gene_introns[gene2]
        del self.gene_to_transcripts[gene2]
        del self.gene_strands[gene2]

        # update scores
        new_scores = {}
        for gene_pair in self.scores.keys():
            if gene2 in gene_pair:
                continue
            elif gene1 in gene_pair:
                new_scores[gene_pair] = self.count_score(*gene_pair)
            else:
                new_scores[gene_pair] = self.scores[gene_pair]
        self.scores = new_scores

    def join_transcripts(self):
        self.count_scores()
        while len(self.gene_to_transcripts) > 1:
            best_gene_pair = max(self.scores, key=self.scores.get)
            if self.scores[best_gene_pair] < 0.1:
                break
            self.merge_genes(*best_gene_pair)

        transcript_to_new_gene_id = {}
        for gene_id, t_list in self.gene_to_transcripts.items():
            for transcript_id in t_list:
                transcript_to_new_gene_id[transcript_id] = gene_id
        for model in self.transcipt_model_storage:
            model.gene_id = transcript_to_new_gene_id[model.transcript_id]

        return self.transcipt_model_storage
