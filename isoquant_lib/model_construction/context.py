############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Shared state for transcript model construction, extracted from the former
# GraphBasedModelConstructor god class:
#   ModelStore                 -- the mutable accumulation of models + the
#                                 per-read construction bookkeeping, plus the
#                                 read-binding primitives that every stage uses.
#   ModelConstructionContext   -- (added in a later step) the read-only inputs
#                                 and the intron graph.

from collections import defaultdict
from enum import Enum, unique
from typing import Dict, List

from ..gene_info import StrandDetector, TranscriptModel
from isoquant_lib.model_construction.intron_graph import IntronGraph
from ..assignment.isoform_assignment import ReadAssignment
from ..assignment.long_read_assigner import LongReadAssigner
from ..assignment.long_read_profiles import CombinedProfileConstructor
from .intron_path import IntronPathProcessor, IntronPathStorage


@unique
class StrandnessReportingLevel(Enum):
    only_canonical = 1
    only_stranded = 2
    all = 3
    auto = 10


def reference_unique_known_isoform(assignment: ReadAssignment):
    # The known reference isoform this read was uniquely assigned to by the
    # reference step, or None. All unique reference assignments target an
    # annotated (known) isoform, so the returned id also names the model built
    # from it (from_reference_transcript keeps the reference id). Shared by the
    # construction read assigner and the final model read counter.
    if assignment.assignment_type.is_unique() and assignment.isoform_matches:
        return assignment.isoform_matches[0].assigned_transcript
    return None


class ModelStore:
    """Single source of truth for constructed models and their read bindings.

    Every construction / filtering / assignment stage appends models here and
    binds reads through these primitives, so the bookkeeping dicts stay
    consistent regardless of which stage runs.
    """

    # Process-wide dedup of known isoforms already emitted as models. A new
    # constructor (and store) is built per gene block, but this set is shared
    # across all of them in a worker process and is never reset -- it preserves
    # the former GraphBasedModelConstructor class-level behaviour so a known
    # isoform is not emitted twice across overlapping gene blocks.
    detected_known_isoforms = set()
    extended_transcript_ids = set()

    def __init__(self, gene_info):
        self.gene_info = gene_info
        self.transcript_model_storage = []
        self.transcript_read_ids = defaultdict(list)
        self.internal_counter = defaultdict(int)
        self.read_assignment_counts = defaultdict(int)
        # read_id -> transcript id the read was bound to during construction
        # (via save_assigned_read). Lets the final counting tell construction-
        # supported reads apart from reads that only reach the re-profiling step.
        self.construction_assignment: Dict[str, str] = {}
        # known transcript id -> reads used to build it as a non-FL model;
        # resolve_read_ownership_conflicts uses this after filtering.
        self.nonfl_known_reads: Dict[str, List[ReadAssignment]] = {}

    def save_assigned_read(self, read_assignment: ReadAssignment, transcript_id: str) -> None:
        read_id = read_assignment.read_id
        self.transcript_read_ids[transcript_id].append(read_assignment)
        self.internal_counter[transcript_id] += 1
        self.read_assignment_counts[read_id] += 1
        self.construction_assignment[read_id] = transcript_id

    def move_read_to_model(self, read_assignment: ReadAssignment,
                           source_id: str, target_id: str) -> None:
        # Move one read's construction membership from source_id to target_id (an
        # alternative-end NIC derived from it). ReadAssignment has no __eq__, so
        # list.remove strips the exact object; net read_assignment_counts is
        # unchanged (-1 here, +1 in save_assigned_read).
        self.transcript_read_ids[source_id].remove(read_assignment)
        self.internal_counter[source_id] -= 1
        self.read_assignment_counts[read_assignment.read_id] -= 1
        self.save_assigned_read(read_assignment, target_id)

    def delete_from_storage(self, transcript_id: str) -> None:
        for a in self.transcript_read_ids[transcript_id]:
            self.read_assignment_counts[a.read_id] -= 1
        del self.transcript_read_ids[transcript_id]
        del self.internal_counter[transcript_id]

    def transcript_from_reference(self, isoform_id: str) -> TranscriptModel:
        return TranscriptModel.from_reference_transcript(self.gene_info, isoform_id)

    def mapping_quality(self, transcript_id: str) -> float:
        mapq = 0
        for a in self.transcript_read_ids[transcript_id]:
            mapq += a.mapping_quality
        return mapq / len(self.transcript_read_ids[transcript_id])


class ModelConstructionContext:
    """Read-only inputs + the intron graph shared by every construction stage.

    Built once per gene block. Holds the annotation/gene info, CLI args, string
    pools, the ID distributor, the polyA/TSS predictions, the derived assigner /
    profile constructor / strand detector, and (after build_graph) the intron
    graph and known-isoform lookups. Also owns the shared ModelStore.
    """

    def __init__(self, gene_info, chr_record, args, id_distributor, string_pools,
                 polya_predictions, tss_predictions, store,
                 grouping_strategy_names=None, use_technical_replicas=False):
        self.gene_info = gene_info
        self.chr_record = chr_record
        self.args = args
        self.id_distributor = id_distributor
        self.string_pools = string_pools
        # Predicted polyA / TSS sites for this gene as (transcript_id, genomic
        # position), reused from the per-gene terminal-position counters to
        # refine intron-graph terminal vertices. The transcript id lets the
        # graph admit a prediction only on its own transcript's terminal intron.
        self.polya_predictions = polya_predictions
        self.tss_predictions = tss_predictions
        self.store = store

        self.grouping_strategy_names = grouping_strategy_names if grouping_strategy_names else []
        self.use_technical_replicas = use_technical_replicas
        self.file_name_group_idx = (self.grouping_strategy_names.index("file_name")
                                    if "file_name" in self.grouping_strategy_names else -1)

        # Constructed (novel) model ends are always refined to the trained
        # polyA/TSS peaks. Alternative-end NIC *creation* needs an annotation to
        # refine against, so it is gated on --genedb.
        self.create_alt_ts_nics = bool(getattr(args, 'genedb', None))
        # TSS model is used only with full-length evidence; read starts are
        # unreliable otherwise.
        self.use_tss_model = bool(getattr(args, 'fl_data', False))
        # Opt-in (--novel_apa): also create alternative-end isoforms for novel chains.
        self.novel_apa = bool(getattr(args, 'novel_apa', False))

        self.strand_detector = StrandDetector(self.chr_record)
        self.intron_genes = defaultdict(set)
        self.set_gene_properties()

        self.assigner = LongReadAssigner(self.gene_info, self.args, self.string_pools)
        self.profile_constructor = CombinedProfileConstructor(self.gene_info, self.args)

        # Graph-derived, populated by build_graph().
        self.intron_graph = None
        self.path_processor = None
        self.path_storage = None
        self.known_isoforms_in_graph = {}
        self.known_introns = set()
        self.known_isoforms_in_graph_ids = {}

    def get_transcript_id(self):
        return self.id_distributor.increment()

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

        return None

    def get_known_spliced_isoforms(self, gene_info, s="known"):
        known_isoforms = {}
        for isoform_id in gene_info.all_isoforms_introns:
            isoform_introns = gene_info.all_isoforms_introns[isoform_id]
            intron_path = self.path_processor.thread_introns(isoform_introns)
            if not intron_path:
                continue
            known_isoforms[tuple(intron_path)] = isoform_id
        return known_isoforms

    def build_graph(self, read_assignment_storage):
        self.intron_graph = IntronGraph(self.args, self.gene_info, read_assignment_storage,
                                        polya_predictions=self.polya_predictions,
                                        tss_predictions=self.tss_predictions)
        self.path_processor = IntronPathProcessor(self.args, self.intron_graph)
        self.path_storage = IntronPathStorage(self.args, self.path_processor)
        self.path_storage.fill(read_assignment_storage)
        self.known_isoforms_in_graph = self.get_known_spliced_isoforms(self.gene_info)
        self.known_introns = set(self.gene_info.intron_profiles.features)
        for intron_path, isoform_id in self.known_isoforms_in_graph.items():
            self.known_isoforms_in_graph_ids[isoform_id] = intron_path
