############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Stage 4c: honest per-read assignment against the final model set. Produces one
# ReadAssignment per read (construction member / reference-unique known /
# discarded / leftover), feeds them to the discovered transcript & gene counters,
# and collects them for the model read_info output. Replaces forward_counts.

import copy
import logging
from typing import Dict, List

from ..gene_info import GeneInfo, TranscriptModelType
from ..assignment.isoform_assignment import (
    ReadAssignment,
    ReadAssignmentType,
    MatchClassification,
    IsoformMatch,
)
from ..assignment.long_read_assigner import LongReadAssigner
from ..assignment.long_read_profiles import CombinedProfileConstructor
from .context import reference_unique_known_isoform

logger = logging.getLogger('IsoQuant')


class ModelReadCounter:
    """Honest per-read assignment + counting against the final model set."""

    def __init__(self, ctx, transcript_counter, gene_counter):
        self.ctx = ctx
        self.store = ctx.store
        self.gene_info = ctx.gene_info
        self.args = ctx.args
        self.string_pools = ctx.string_pools
        self.transcript_counter = transcript_counter
        self.gene_counter = gene_counter
        self.model_read_assignments: List[ReadAssignment] = []

    def _log_uniqueknown_on_novel(self, transcript2type: Dict[str, TranscriptModelType]) -> None:
        # Diagnostic: how often a read the reference step assigned uniquely to a
        # known isoform ends up on a *novel* model. After the re-assignment gate
        # this can only happen while the read is consumed building that model
        # (FL path / alt-end), so the count isolates the construction-phase drag.
        dragged = 0
        for t_id, reads in self.store.transcript_read_ids.items():
            if transcript2type.get(t_id) == TranscriptModelType.known:
                continue
            for a in reads:
                x = reference_unique_known_isoform(a)
                if x is not None and x != t_id:
                    dragged += 1
                    logger.debug("Read %s uniquely assigned to known %s now on novel model %s" %
                                 (a.read_id, x, t_id))
        if dragged:
            logger.debug("%d uniquely-known reads bound to novel models (construction-phase)" % dragged)

    def _single_model_assigner(self, model, cache):
        # Cache one (assigner, profile_constructor) per model over a single-isoform
        # GeneInfo, using the graph vertex-collapse tolerance (graph_clustering_distance)
        # as the junction delta so a construction read whose splice sites the intron
        # graph snapped within that reach still matches its model, instead of being
        # rejected at the tighter assignment delta. quick_mode=False for real events.
        entry = cache.get(model.transcript_id)
        if entry is None:
            gi = GeneInfo.from_models([model], self.args.delta)
            entry = (LongReadAssigner(gi, self._graph_tolerance_args, self.string_pools, quick_mode=False),
                     CombinedProfileConstructor(gi, self._graph_tolerance_args))
            cache[model.transcript_id] = entry
        return entry

    def _model_read_assignment_single(self, source, transcript_id, model_by_id, cache):
        # A read that supports this model must count toward it. Reuse the source
        # match when the read was reference unique/umd on a *known* model (real
        # classification + events, no re-profiling). Otherwise re-profile against
        # the single model at the graph tolerance; keep the result when it is a
        # clean consistent match (real classification + events), else force the
        # read onto the model as unique_minor_difference (events dropped) so it
        # still counts as unique and points at the model it built.
        model = model_by_id[transcript_id]
        if (model.transcript_type == TranscriptModelType.known and
                source.assignment_type in (ReadAssignmentType.unique,
                                           ReadAssignmentType.unique_minor_difference)):
            reuse_match = next((m for m in source.isoform_matches
                                if m.assigned_transcript == transcript_id), None)
            if reuse_match is not None:
                return ReadAssignment(source.read_id, source.assignment_type,
                                      self.string_pools, match=[reuse_match])
        assigner, profile_constructor = self._single_model_assigner(model, cache)
        profile = profile_constructor.construct_profiles(source.corrected_exons, source.polya_info, [])
        ra = assigner.assign_to_isoform(source.read_id, profile)
        if (ra.assignment_type.is_consistent() and ra.isoform_matches and
                ra.isoform_matches[0].assigned_transcript == transcript_id):
            return ra
        forced_match = IsoformMatch(MatchClassification.full_splice_match, self.string_pools,
                                    assigned_gene=model.gene_id, assigned_transcript=transcript_id)
        return ReadAssignment(source.read_id, ReadAssignmentType.unique_minor_difference,
                              self.string_pools, match=[forced_match])

    def _build_discarded_assignment(self, read_id, gene, surviving_genes):
        # A read whose reference-unique known isoform was dropped: no transcript.
        # It still belongs to its gene, so it counts toward that gene when the gene
        # has a surviving model (gene_assignment_type=unique makes the gene counter
        # count it); otherwise it is not assigned anywhere.
        if gene is not None and gene in surviving_genes:
            match = IsoformMatch(MatchClassification.genic, self.string_pools, assigned_gene=gene)
            ra = ReadAssignment(read_id, ReadAssignmentType.discarded, self.string_pools, match=[match])
            ra.gene_assignment_type = ReadAssignmentType.unique
            return ra
        return ReadAssignment(read_id, ReadAssignmentType.discarded, self.string_pools)

    def build_model_read_assignments(self, read_assignments):
        """Produce one honest ReadAssignment per read against the final model set,
        feed each to the discovered transcript & gene counters (which apply the
        per-type counting strategy) and collect them for the model read_info
        printer. Replaces the artificial forward_counts, so read2transcripts comes
        for free. Categories (construction checked first, mirroring the pipeline
        where construction-bound reads skip the re-profiling gate):

        1. construction member  -> re-profile vs the single model it built at the
           graph tolerance, forced onto it as unique/umd (it must count toward the
           model it contributed to).
        2. reference-unique member on a *kept* known isoform -> reuse source match.
        3. discarded (reference-unique to a *dropped* known isoform) -> no transcript;
           counts toward the gene only if that gene still has a surviving model.
        4. leftover -> assign vs the full model set (quick_mode=False): honest
           unique/umd/ambiguous/inconsistent, kept regardless of consistency.
        """
        self.model_read_assignments = []
        # Junction delta for the single-model (construction/known) re-profiling:
        # the intron graph's vertex-collapse tolerance, so reads snapped within it
        # match their model. Never below the assignment delta.
        self._graph_tolerance_args = copy.copy(self.args)
        _gcd = getattr(self.args, "graph_clustering_distance", None)
        self._graph_tolerance_args.delta = max(self.args.delta, _gcd) if _gcd else self.args.delta
        transcript2gene = {t.transcript_id: t.gene_id for t in self.store.transcript_model_storage}
        transcript2type = {t.transcript_id: t.transcript_type for t in self.store.transcript_model_storage}
        self._log_uniqueknown_on_novel(transcript2type)

        model_gene_info = GeneInfo.from_models(self.store.transcript_model_storage, self.args.delta)
        model_by_id = {t.transcript_id: t for t in self.store.transcript_model_storage}
        model_ids = set(model_by_id.keys())
        surviving_genes = set(transcript2gene.values())

        single_cache = {}
        full_assigner = None
        full_profile_constructor = None

        for source in read_assignments:
            read_id = source.read_id
            constr_tid = self.store.construction_assignment.get(read_id)
            if constr_tid is not None and constr_tid in model_ids:
                ra = self._model_read_assignment_single(source, constr_tid, model_by_id, single_cache)
            else:
                known_iso = reference_unique_known_isoform(source)
                if known_iso is not None and known_iso in model_ids:
                    ra = self._model_read_assignment_single(source, known_iso, model_by_id, single_cache)
                elif known_iso is not None:
                    gene = self.gene_info.gene_id_map.get(known_iso)
                    ra = self._build_discarded_assignment(read_id, gene, surviving_genes)
                elif model_ids:
                    if full_assigner is None:
                        full_assigner = LongReadAssigner(model_gene_info, self.args,
                                                         self.string_pools, quick_mode=False)
                        full_profile_constructor = CombinedProfileConstructor(model_gene_info, self.args)
                    profile = full_profile_constructor.construct_profiles(source.corrected_exons,
                                                                          source.polya_info, [])
                    ra = full_assigner.assign_to_isoform(read_id, profile)
                else:
                    ra = ReadAssignment(read_id, ReadAssignmentType.noninformative, self.string_pools)

            ra.gene_info = model_gene_info
            self._copy_read_aux(ra, source)
            self.model_read_assignments.append(ra)
            self.transcript_counter.add_read_info(ra)
            self.gene_counter.add_read_info(ra)

        self.transcript_counter.add_confirmed_features([m.transcript_id for m in self.store.transcript_model_storage])
        self.gene_counter.add_confirmed_features({transcript2gene[m.transcript_id] for m in self.store.transcript_model_storage})

    @staticmethod
    def _copy_read_aux(ra: ReadAssignment, source: ReadAssignment) -> None:
        # Carry every read_info aux field over from the source reference
        # assignment so the model read_info output matches the main read_info
        # format. ra shares the source's string pools, so the interned integer
        # ids (chr/barcode/umi/read_group) can be copied directly.
        ra.chr_id_int = source.chr_id_int
        ra.strand = source.strand
        ra.exons = source.exons
        ra.corrected_exons = source.corrected_exons
        ra.polyA_found = source.polyA_found
        ra.cage_found = source.cage_found
        ra.polya_info = source.polya_info
        ra.barcode_id = source.barcode_id
        ra.umi_id = source.umi_id
        ra.read_group_ids = source.read_group_ids
        ra.additional_attributes = source.additional_attributes

