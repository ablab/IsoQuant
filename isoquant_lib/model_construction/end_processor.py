############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Stage 3 of model construction: polyA / TSS terminal refinement.
#   - correct_novel_transcript_ends / _mark_terminal_confirmation polish a
#     model's own ends toward confident read-terminus peaks;
#   - _add_known_alternative_end_models / derive_alternative_end_models spin off
#     alternative-end NIC siblings and hand them their supporting reads.
# Operates on the shared ModelStore via the ModelConstructionContext.

import logging
from collections import defaultdict
from typing import Callable, Dict, List, Optional, Tuple, TYPE_CHECKING

from ..common import TranscriptNaming, junctions_from_blocks
from ..gene_info import TranscriptModel, TranscriptModelType
from ..assignment.isoform_assignment import ReadAssignment
from ..terminal_prediction.terminal_peaks import detect_peaks, get_polya_model, get_tss_model

if TYPE_CHECKING:
    from xgboost import XGBClassifier

logger = logging.getLogger('IsoQuant')


class TranscriptEndProcessor:
    """Terminal (polyA/TSS) end refinement and alternative-end NIC discovery."""

    def __init__(self, ctx):
        self.ctx = ctx
        self.store = ctx.store
        self.gene_info = ctx.gene_info
        self.args = ctx.args
        self.use_tss_model = ctx.use_tss_model
        self.novel_apa = ctx.novel_apa

    def correct_novel_transcript_ends(self, transcript_model: TranscriptModel,
                                      assigned_reads: List[ReadAssignment]) -> None:
        # Part 1: detector-based per-transcript end refinement. Build read-terminus
        # histograms from this model's own reads and snap each unsupported terminal
        # exon end to the dominant confident polyA/TSS peak (XGBoost detect_peaks);
        # both ends are refined in one pass over the same reads so the chosen
        # 5'/3' stay concordant. polyA always; TSS only when use_tss_model
        # (--fl_data). No confident peak / no model -> positional fallback, then
        # leave as is. Never creates or drops a model.
        logger.debug("Verifying ends for transcript %s" % transcript_model.transcript_id)
        transcript_start = transcript_model.exon_blocks[0][0]
        transcript_end = transcript_model.exon_blocks[-1][1]
        first_exon_right = transcript_model.exon_blocks[0][1]
        last_exon_left = transcript_model.exon_blocks[-1][0]
        strand = transcript_model.strand

        start_reads, end_reads = self._terminal_histograms(transcript_model, assigned_reads)
        # _terminal_histograms maps position -> reads; detect_peaks / _peak_boundary and
        # _closest_inward want position -> count, so collapse the read lists to counts
        # (same as derive_alternative_end_models).
        start_hist = {p: len(rs) for p, rs in start_reads.items()}
        end_hist = {p: len(rs) for p, rs in end_reads.items()}
        start_supported = any(abs(p - transcript_start) <= self.args.apa_delta for p in start_hist)
        end_supported = any(abs(p - transcript_end) <= self.args.apa_delta for p in end_hist)

        if not start_supported:
            model = self._terminal_model(strand, left=True)
            if model is not None:
                new_start = self._peak_boundary(start_hist, model, lambda pos: pos < first_exon_right)
            else:
                new_start = self._closest_inward(sorted(start_hist.keys()), transcript_start, greater=True)
            if new_start is not None and new_start != transcript_start and new_start < first_exon_right:
                logger.debug("Changed start for transcript %s: from %d to %d" %
                             (transcript_model.transcript_id, transcript_start, new_start))
                transcript_model.exon_blocks[0] = (new_start, first_exon_right)

        if not end_supported:
            model = self._terminal_model(strand, left=False)
            if model is not None:
                new_end = self._peak_boundary(end_hist, model, lambda pos: pos > last_exon_left)
            else:
                new_end = self._closest_inward(sorted(end_hist.keys(), reverse=True), transcript_end, greater=False)
            if new_end is not None and new_end != transcript_end and new_end > last_exon_left:
                logger.debug("Changed end for transcript %s: from %d to %d" %
                             (transcript_model.transcript_id, transcript_end, new_end))
                transcript_model.exon_blocks[-1] = (last_exon_left, new_end)

    def _mark_terminal_confirmation(self, model: TranscriptModel) -> None:
        # Record whether the model's polyA (3') / TSS (5') terminus is supported by a
        # read-terminus peak among its own reads, using the same histograms as
        # correct_novel_transcript_ends (the 3' side counts polyA-confirmed reads only).
        # Called after end refinement so the check uses the model's final ends.
        reads = self.store.transcript_read_ids[model.transcript_id]
        if not reads:
            model.polya_confirmed = model.tss_confirmed = False
            return
        start_reads, end_reads = self._terminal_histograms(model, reads)
        model_start = model.exon_blocks[0][0]
        model_end = model.exon_blocks[-1][1]
        start_supported = any(abs(p - model_start) <= self.args.apa_delta for p in start_reads)
        end_supported = any(abs(p - model_end) <= self.args.apa_delta for p in end_reads)
        if model.strand == '+':      # genomic-left = 5' TSS, genomic-right = 3' polyA
            model.tss_confirmed, model.polya_confirmed = start_supported, end_supported
        elif model.strand == '-':    # reversed
            model.polya_confirmed, model.tss_confirmed = start_supported, end_supported
        else:                        # unstranded: no polyA orientation
            model.polya_confirmed = model.tss_confirmed = False

    def _terminal_model(self, strand: str, left: bool) -> Optional["XGBClassifier"]:
        # Which trained model applies to a genomic-side boundary. For '+' the
        # right end is the polyA site and the left is the TSS; for '-' reversed.
        # Returns None (-> positional fallback) for unstranded, or a side whose
        # model is disabled: polyA needs the flag, TSS additionally needs fl_data.
        if strand == '+':
            is_polya_side = not left
        elif strand == '-':
            is_polya_side = left
        else:
            return None
        if is_polya_side:
            # polyA end refinement is always on (polishing needs no annotation).
            return get_polya_model()
        return get_tss_model() if self.use_tss_model else None

    @staticmethod
    def _peak_boundary(histogram: Dict[int, int], model: "XGBClassifier",
                       valid: Callable[[int], bool]) -> Optional[int]:
        # Best-supported detected peak satisfying the terminal-exon clamp.
        if not histogram:
            return None
        peaks = [p for p in detect_peaks(histogram, model) if valid(p.position)]
        if not peaks:
            return None
        return max(peaks, key=lambda p: p.count).position

    @staticmethod
    def _closest_inward(sorted_positions: List[int], boundary: int, greater: bool) -> Optional[int]:
        # Positional fallback: nearest read terminus on the inward side of the
        # current boundary (matches the original end-correction behaviour).
        for pos in sorted_positions:
            if greater and pos > boundary:
                return pos
            if not greater and pos < boundary:
                return pos
        return None

    @staticmethod
    def _confirmed_polya_pos(assignment: ReadAssignment, strand: str) -> Optional[int]:
        # Detected polyA cleavage site of a read, matching the trained-model
        # domain (PolyACounter): only polyA-confirmed reads, at external_polya_pos
        # ('+') / external_polyt_pos ('-'). None otherwise.
        info = getattr(assignment, 'polya_info', None)
        if not getattr(assignment, 'polyA_found', False) or info is None:
            return None
        if strand == '+' and info.external_polya_pos != -1:
            return info.external_polya_pos
        if strand == '-' and info.external_polyt_pos != -1:
            return info.external_polyt_pos
        return None

    def _terminal_histograms(self, source_model: TranscriptModel,
                             assigned_reads: List[ReadAssignment]) \
            -> Tuple[Dict[int, List[ReadAssignment]], Dict[int, List[ReadAssignment]]]:
        # Build genomic-left (start) and genomic-right (end) read-terminus maps
        # matching what each trained model was fit on:
        #   - 3' polyA side: polyA-CONFIRMED reads at the detected cleavage site
        #     (so a low overall polyA rate naturally raises the support bar);
        #   - 5' TSS side: all stranded reads' alignment ends.
        # Unstranded -> alignment ends on both sides (no polyA orientation).
        # Each map is position -> the reads terminating there, so a peak's reads
        # can later be moved onto the alternative-end model derived from it.
        strand = source_model.strand
        first_exon_right = source_model.exon_blocks[0][1]
        last_exon_left = source_model.exon_blocks[-1][0]
        start_reads = defaultdict(list)
        end_reads = defaultdict(list)
        for a in assigned_reads:
            ex = a.corrected_exons
            if strand == '+':
                if ex[0][0] < first_exon_right:           # 5' TSS (all reads)
                    start_reads[ex[0][0]].append(a)
                pos = self._confirmed_polya_pos(a, '+')   # 3' polyA (confirmed)
                if pos is not None and pos > last_exon_left:
                    end_reads[pos].append(a)
            elif strand == '-':
                pos = self._confirmed_polya_pos(a, '-')   # 3' polyA (confirmed)
                if pos is not None and pos < first_exon_right:
                    start_reads[pos].append(a)
                if ex[-1][1] > last_exon_left:            # 5' TSS (all reads)
                    end_reads[ex[-1][1]].append(a)
            else:
                if ex[0][0] < first_exon_right:
                    start_reads[ex[0][0]].append(a)
                if ex[-1][1] > last_exon_left:
                    end_reads[ex[-1][1]].append(a)
        return start_reads, end_reads

    @staticmethod
    def _intron_chain_key(model: TranscriptModel) -> Tuple[Tuple[int, int], ...]:
        # Terminal-end-independent identity of a transcript: the tuple of its
        # introns derived from exon blocks, using the project-wide
        # junctions_from_blocks convention (actual intron coordinates, +1/-1).
        # Monoexon -> empty tuple.
        return tuple(junctions_from_blocks(model.exon_blocks))

    def _add_known_alternative_end_models(self) -> None:
        # Part 2: for each known (reference) model, peak-call its own assigned
        # reads and emit a NIC for every confident alternative polyA/TSS end,
        # keeping the known (union). Deduplicated against everything already in
        # storage (notably the graph-level NICs) by intron chain + both ends
        # within apa_delta, so the same alternative end is never emitted twice.
        existing_pairs = defaultdict(list)
        # Seed with every reference isoform too: an alternative-end peak that
        # lands on another annotated isoform of the same chain is that known, not
        # a novel end, so it must be suppressed even if that reference was not
        # emitted as a model in this locus.
        for ref_exons in self.gene_info.all_isoforms_exons.values():
            ck = tuple(junctions_from_blocks(ref_exons))
            existing_pairs[ck].append((ref_exons[0][0], ref_exons[-1][1]))
        for m in self.store.transcript_model_storage:
            existing_pairs[self._intron_chain_key(m)].append(
                (m.exon_blocks[0][0], m.exon_blocks[-1][1]))

        new_models = []
        for model in self.store.transcript_model_storage:
            if model.transcript_type == TranscriptModelType.known:
                pass
            elif not self.novel_apa:
                # Part 2 default: known transcripts only. Part 3 (--novel_apa)
                # also spins off alternative-end siblings for novel chains.
                continue
            source_id = model.transcript_id
            for nic, nic_reads in self.derive_alternative_end_models(
                    model, self.store.transcript_read_ids[source_id]):
                ck = self._intron_chain_key(nic)
                ns, ne = nic.exon_blocks[0][0], nic.exon_blocks[-1][1]
                if any(abs(ns - s) <= self.args.apa_delta and abs(ne - e) <= self.args.apa_delta
                       for (s, e) in existing_pairs[ck]):
                    continue
                existing_pairs[ck].append((ns, ne))
                new_models.append(nic)
                # Hand the alternative-end peak's reads to the NIC so it is not a
                # zero-read phantom; keep at least one read on the source.
                for read in nic_reads:
                    if len(self.store.transcript_read_ids[source_id]) <= 1:
                        break
                    self.store.move_read_to_model(read, source_id, nic.transcript_id)
        if new_models:
            logger.debug("Added %d known alternative-end NICs" % len(new_models))
            self.store.transcript_model_storage.extend(new_models)

    def _drop_duplicate_alt_end_models(self) -> None:
        # Final dedup over the whole storage (catches graph-level NICs too): drop
        # any non-known model that is structurally identical (intron chain + both
        # ends within apa_delta) to a reference isoform -> it is that annotated
        # transcript, not a novel end; and collapse non-known models that
        # duplicate an already-kept one.
        ref_pairs = defaultdict(list)
        for ref_exons in self.gene_info.all_isoforms_exons.values():
            ck = tuple(junctions_from_blocks(ref_exons))
            ref_pairs[ck].append((ref_exons[0][0], ref_exons[-1][1]))

        d = self.args.apa_delta
        kept = []
        kept_pairs = defaultdict(list)
        dropped = 0
        for m in self.store.transcript_model_storage:
            if m.transcript_type == TranscriptModelType.known:
                kept.append(m)
                continue
            ck = self._intron_chain_key(m)
            s, e = m.exon_blocks[0][0], m.exon_blocks[-1][1]
            if any(abs(s - rs) <= d and abs(e - re) <= d for rs, re in ref_pairs[ck]) or \
                    any(abs(s - ks) <= d and abs(e - ke) <= d for ks, ke in kept_pairs[ck]):
                dropped += 1
                # free the model's reads for reassignment (graph-level NICs are
                # already counted; post-pass NICs are not yet in the counters)
                if m.transcript_id in self.store.internal_counter:
                    self.store.delete_from_storage(m.transcript_id)
                continue
            kept.append(m)
            kept_pairs[ck].append((s, e))
        if dropped:
            logger.debug("Dropped %d duplicate / reference-matching alt-end models" % dropped)
        self.store.transcript_model_storage = kept

    def derive_alternative_end_models(self, source_model: TranscriptModel,
                                      assigned_reads: List[ReadAssignment]) \
            -> List[Tuple[TranscriptModel, List[ReadAssignment]]]:
        # Confident alternative polyA/TSS ends for a transcript, from its own
        # reads (per-transcript, so 5'/3' stay concordant). Returns (new model,
        # supporting reads) per alternative end, each model changing a single
        # terminal. The supporting reads are the source reads whose terminus sits
        # at that peak (within apa_delta); the caller moves them onto the NIC so it
        # is not a zero-read phantom. Each read backs at most one NIC. polyA always;
        # TSS only when use_tss_model (--fl_data). A known source yields
        # novel-in-catalog siblings; a novel source keeps its own type.
        if not assigned_reads:
            return []
        exon_blocks = source_model.exon_blocks
        strand = source_model.strand
        annotated_start = exon_blocks[0][0]
        annotated_end = exon_blocks[-1][1]
        first_exon_right = exon_blocks[0][1]
        last_exon_left = exon_blocks[-1][0]
        sibling_type = (TranscriptModelType.novel_in_catalog
                        if source_model.transcript_type == TranscriptModelType.known
                        else source_model.transcript_type)

        start_reads, end_reads = self._terminal_histograms(source_model, assigned_reads)
        start_hist = {p: len(rs) for p, rs in start_reads.items()}
        end_hist = {p: len(rs) for p, rs in end_reads.items()}

        claimed = set()  # id(read) already routed to a NIC (one NIC per read, both sides)

        def assign_to_peaks(reads_map, alt_peaks, annotated_pos):
            # Route each read to the *nearest* terminus among the annotated end and
            # the alternative-end peaks (within apa_delta). Reads nearest the
            # annotated end stay on the source; reads nearest an alt peak back that
            # peak's NIC. Nearest (not first-within-window) avoids a peak greedily
            # claiming a neighbour's reads and starving it; the shared claimed set
            # keeps a both-ends-alternative read on a single NIC (it is moved once).
            per_peak = {p: [] for p in alt_peaks}
            candidates = list(alt_peaks) + [annotated_pos]
            for p, rs in reads_map.items():
                nearest = min(candidates, key=lambda c: (abs(c - p), c))
                if nearest != annotated_pos and abs(nearest - p) <= self.args.apa_delta:
                    for r in rs:
                        if id(r) not in claimed:
                            claimed.add(id(r))
                            per_peak[nearest].append(r)
            return per_peak

        new_models = []
        start_model = self._terminal_model(strand, left=True)
        if start_model is not None:
            start_peaks = self._alternative_end_positions(start_hist, start_model, annotated_start,
                                                          lambda p: p < first_exon_right)
            per_peak = assign_to_peaks(start_reads, start_peaks, annotated_start)
            for pos in start_peaks:
                nic = self._nic_model_with_boundary(source_model, sibling_type, start=pos)
                new_models.append((nic, per_peak[pos]))
        end_model = self._terminal_model(strand, left=False)
        if end_model is not None:
            end_peaks = self._alternative_end_positions(end_hist, end_model, annotated_end,
                                                        lambda p: p > last_exon_left)
            per_peak = assign_to_peaks(end_reads, end_peaks, annotated_end)
            for pos in end_peaks:
                nic = self._nic_model_with_boundary(source_model, sibling_type, end=pos)
                new_models.append((nic, per_peak[pos]))
        return new_models

    def _alternative_end_positions(self, histogram: Dict[int, int], model: "XGBClassifier",
                                   annotated_pos: int, clamp: Callable[[int], bool]) -> List[int]:
        # Detected peaks representing a genuine alternative end: inside the
        # terminal exon, > apa_delta from the annotated end, and clearing BOTH an
        # absolute (min_novel_count) and a RELATIVE support cutoff. The relative
        # cutoff (terminal_position_rel x the transcript's dominant terminal peak,
        # reusing the intron-graph terminal threshold) rejects minor secondary
        # peaks that are not real co-expressed isoforms and self-normalizes per
        # transcript, so it does not depend on absolute depth or fit one isoform.
        if not histogram:
            return []
        peaks = [p for p in detect_peaks(histogram, model) if clamp(p.position)]
        if not peaks:
            return []
        dominant = max(p.count for p in peaks)
        cutoff = max(self.args.min_novel_count, self.args.terminal_position_rel * dominant)
        return [p.position for p in peaks
                if abs(p.position - annotated_pos) > self.args.apa_delta and p.count >= cutoff]

    def _nic_model_with_boundary(self, source_model: TranscriptModel, transcript_type: TranscriptModelType,
                                 start: Optional[int] = None, end: Optional[int] = None) -> TranscriptModel:
        # Build an alternative-end model from a source one by replacing a single
        # terminal coordinate; the intron chain (internal exons) is unchanged.
        # transcript_type is novel_in_catalog for a known source, or the source's
        # own type for a novel source.
        exon_blocks = [tuple(e) for e in source_model.exon_blocks]
        if start is not None:
            exon_blocks[0] = (start, exon_blocks[0][1])
        if end is not None:
            exon_blocks[-1] = (exon_blocks[-1][0], end)
        suffix = (TranscriptNaming.nic_transcript_suffix
                  if transcript_type == TranscriptModelType.novel_in_catalog
                  else TranscriptNaming.nnic_transcript_suffix)
        new_transcript_id = TranscriptNaming.transcript_prefix + str(self.ctx.get_transcript_id())
        new_model = TranscriptModel(
            self.gene_info.chr_id, source_model.strand,
            new_transcript_id + ".%s" % self.gene_info.chr_id + suffix,
            source_model.gene_id, exon_blocks, transcript_type)
        new_model.intron_path = source_model.intron_path
        logger.debug("Adding alternative-end model %s from %s (start=%s, end=%s)" %
                     (new_model.transcript_id, source_model.transcript_id, str(start), str(end)))
        return new_model

