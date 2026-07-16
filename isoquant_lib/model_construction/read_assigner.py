############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Stage 4b: construction-phase read (re)assignment. Binds leftover reads to the
# constructed models (gated by the unique-known rule), then resolves ownership
# conflicts between FL-path models and assignment-based knowns after filtering.

import logging
from collections import defaultdict
from typing import Dict

from ..gene_info import GeneInfo, TranscriptModelType
from ..assignment.long_read_assigner import LongReadAssigner
from ..assignment.long_read_profiles import CombinedProfileConstructor
from .context import reference_unique_known_isoform

logger = logging.getLogger('IsoQuant')


class ConstructionReadAssigner:
    """Binds reads to constructed models and resolves ownership conflicts."""

    def __init__(self, ctx):
        self.ctx = ctx
        self.store = ctx.store
        self.args = ctx.args
        self.string_pools = ctx.string_pools

    def resolve_read_ownership_conflicts(self) -> None:
        # construct_nonfl_isoforms builds an assignment-based known model and re-binds
        # (save_assigned_read overwrites construction_assignment) reads that an FL path
        # already bound to another model. Hand such reads back to the FL-path model when
        # the assignment-based known would otherwise steal the FL model's only support:
        #   - a novel FL model (novel-vs-known): the novel is a distinct discovered
        #     transcript and always wins its reads;
        #   - an FL-path known left with zero owned reads (known-vs-known): a near-identical
        #     sibling of the assignment-based known (their intron chains differ by a few bp
        #     at one splice site, collapsed in the intron graph) stole all its reads, so it
        #     is a zero-count phantom -> reclaim them. An FL known that still owns reads is
        #     a legitimate co-expressed isoform and is left alone (no sibling consolidation).
        # After filtering, drop assignment-based knowns left with no reads.
        if not self.store.nonfl_known_reads:
            return
        # reads currently owned by each model (its construction_assignment target)
        owned = defaultdict(int)
        for mid in self.store.construction_assignment.values():
            owned[mid] += 1
        # read_id -> a surviving winner model that should reclaim it: a novel FL model,
        # or an FL-path known that is currently a zero-count phantom (0 owned reads).
        read_to_fl: Dict[str, str] = {}
        for model in self.store.transcript_model_storage:
            if model.transcript_id in self.store.nonfl_known_reads:
                continue  # assignment-based known: a potential loser, never a winner
            if (model.transcript_type == TranscriptModelType.known and
                    owned.get(model.transcript_id, 0) > 0):
                continue  # FL-path known that still owns reads: co-expressed, leave it
            for a in self.store.transcript_read_ids[model.transcript_id]:
                read_to_fl.setdefault(a.read_id, model.transcript_id)

        surviving_ids = {m.transcript_id for m in self.store.transcript_model_storage}
        emptied_knowns = []
        for known_id, reads in self.store.nonfl_known_reads.items():
            if known_id not in surviving_ids:
                continue
            for a in reads:
                if self.store.construction_assignment.get(a.read_id) != known_id:
                    # read no longer owned by this known (e.g. moved to an alt-end NIC)
                    continue
                fl_id = read_to_fl.get(a.read_id)
                if fl_id is None or fl_id == known_id:
                    continue
                # still shared with a surviving FL-path model -> the FL model wins
                self.store.construction_assignment[a.read_id] = fl_id
                if a in self.store.transcript_read_ids[known_id]:
                    self.store.transcript_read_ids[known_id].remove(a)
                    self.store.internal_counter[known_id] -= 1
                    self.store.read_assignment_counts[a.read_id] -= 1
            if self.store.internal_counter[known_id] <= 0:
                emptied_knowns.append(known_id)

        if emptied_knowns:
            drop = set(emptied_knowns)
            for known_id in drop:
                self.store.delete_from_storage(known_id)
            self.store.transcript_model_storage = [m for m in self.store.transcript_model_storage
                                             if m.transcript_id not in drop]
            logger.debug("resolve_read_ownership_conflicts: dropped %d empty assignment-based knowns"
                         % len(drop))


    def assign_reads_to_models(self, read_assignments):
        if not self.store.transcript_model_storage:
            for assignment in read_assignments:
                read_id = assignment.read_id
                self.store.read_assignment_counts[read_id] = 0
            logger.debug("No transcripts were assigned")
            return

        logger.debug("Creating artificial GeneInfo from %d transcript models" % len(self.store.transcript_model_storage))
        transcript_model_gene_info = GeneInfo.from_models(self.store.transcript_model_storage, self.args.delta)
        assigner = LongReadAssigner(transcript_model_gene_info, self.args, self.string_pools, quick_mode=True)
        profile_constructor = CombinedProfileConstructor(transcript_model_gene_info, self.args)
        # ids of the isoforms actually built as models here; a known reference
        # isoform keeps its id (from_reference_transcript), so this doubles as
        # the set of surviving known isoforms.
        model_ids = set(transcript_model_gene_info.all_isoforms_exons.keys())
        kept_on_known = 0
        dropped_no_model = 0

        for assignment in read_assignments:
            read_id = assignment.read_id
            if self.store.read_assignment_counts[read_id] > 0:
                # logger.debug("# Read %s was assigned to based on path construction" % read_id)
                continue

            # Consistency gate: a read the reference step assigned uniquely to a
            # known isoform stays on that isoform and is never re-profiled onto a
            # different model. If that isoform was not built as a model (filtered
            # out / unsupported), the read is left unassigned rather than dragged
            # elsewhere. Inert without --genedb (no known reference isoforms).
            known_iso = reference_unique_known_isoform(assignment)
            if known_iso is not None:
                if known_iso in model_ids:
                    self.store.transcript_read_ids[known_iso].append(assignment)
                    self.store.internal_counter[known_iso] += 1
                    self.store.read_assignment_counts[read_id] += 1
                    kept_on_known += 1
                else:
                    self.store.read_assignment_counts[read_id] = 0
                    dropped_no_model += 1
                continue

            read_exons = assignment.corrected_exons
            # logger.debug("# Checking read %s: %s" % (assignment.read_id, str(read_exons)))
            model_combined_profile = profile_constructor.construct_profiles(read_exons, assignment.polya_info, [])
            model_assignment = assigner.assign_to_isoform(assignment.read_id, model_combined_profile)
            model_assignment.read_group = assignment.read_group  # Full list, not just [0]
            # check that no serious contradiction occurs
            if model_assignment.assignment_type.is_consistent():
                matched_isoforms = [m.assigned_transcript for m in model_assignment.isoform_matches]

                if len(matched_isoforms) == 1:
                    self.store.internal_counter[matched_isoforms[0]] += 1
                for m in model_assignment.isoform_matches:
                    self.store.read_assignment_counts[read_id] += 1
                    self.store.transcript_read_ids[m.assigned_transcript].append(assignment)
            else:
                self.store.read_assignment_counts[read_id] = 0

        logger.debug("Gated re-assignment: kept %d reads on their unique known isoform, "
                     "dropped %d (known isoform not built)" % (kept_on_known, dropped_no_model))
