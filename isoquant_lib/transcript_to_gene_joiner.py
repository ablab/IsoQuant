############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict

from .common import intersection_len, interval_len, junctions_from_blocks
from .gene_info import TranscriptModelType

logger = logging.getLogger('IsoQuant')


class TranscriptToGeneJoiner:
    # merge same-strand genes when they share a splice junction or strongly overlap
    STRONG_OVERLAP_FRACTION = 0.5

    def __init__(self, transcipt_model_storage, gene_info):
        self.gene_info = gene_info
        self.transcipt_model_storage = transcipt_model_storage
        self.gene_introns = defaultdict(set)
        self.gene_strands = {}
        self.gene_regions = {}
        self.gene_to_transcripts = defaultdict(set)

        for gene_id in self.gene_info.gene_strands.keys():
            self.gene_strands[gene_id] = self.gene_info.gene_strands[gene_id]
            self.gene_regions[gene_id] = self.gene_info.get_gene_regions()[gene_id]
        for transcript_id in self.gene_info.gene_id_map.keys():
            gene_id = self.gene_info.gene_id_map[transcript_id]
            self.gene_introns[gene_id].update(self.gene_info.all_isoforms_introns[transcript_id])
            self.gene_to_transcripts[gene_id].add(transcript_id)

        for t in self.transcipt_model_storage:
            if t.transcript_type == TranscriptModelType.known:
                continue

            if t.gene_id not in self.gene_regions:
                self.gene_regions[t.gene_id] = (t.exon_blocks[0][0], t.exon_blocks[-1][1])
                self.gene_strands[t.gene_id] = t.strand
            else:
                self.gene_regions[t.gene_id] = (min(self.gene_regions[t.gene_id][0], t.exon_blocks[0][0]),
                                                max(self.gene_regions[t.gene_id][1], t.exon_blocks[-1][1]))
                assert self.gene_strands[t.gene_id] == t.strand
            self.gene_introns[t.gene_id].update(junctions_from_blocks(t.exon_blocks))
            self.gene_to_transcripts[t.gene_id].add(t.transcript_id)
        self.scores = {}

    def count_score(self, gene1, gene2):
        logger.debug("Counting score %s %s" % (gene2, gene1))
        if self.gene_strands[gene1] != self.gene_strands[gene2]:
            return 0.0
        introns1 = self.gene_introns[gene1]
        introns2 = self.gene_introns[gene2]
        shared_introns = len(introns1 & introns2)
        if shared_introns > 0:
            # shared splice junction: highest merge priority, ranked by how much is shared
            return 1.0 + shared_introns / max(1, len(introns1 | introns2))
        gene_range1 = self.gene_regions[gene1]
        gene_range2 = self.gene_regions[gene2]
        inter = intersection_len(gene_range1, gene_range2)
        if inter <= 0:
            return 0.0
        overlap_fraction = inter / min(interval_len(gene_range1), interval_len(gene_range2))
        if overlap_fraction >= self.STRONG_OVERLAP_FRACTION:
            return overlap_fraction
        return 0.0

    def count_scores(self):
        for g1_id in self.gene_to_transcripts.keys():
            for g2_id in self.gene_to_transcripts.keys():
                if g1_id == g2_id or (g1_id in self.gene_info.gene_strands and g2_id in self.gene_info.gene_strands):
                    continue
                gene_pair = tuple(sorted([g1_id, g2_id]))
                if gene_pair not in self.scores:
                    self.scores[gene_pair] = self.count_score(g1_id, g2_id)

    def merge_genes(self, gene1, gene2):
        logger.debug("Merging %s into %s" % (gene2, gene1))
        self.gene_regions[gene1] = (min(self.gene_regions[gene1][0], self.gene_regions[gene2][0]),
                                    max(self.gene_regions[gene1][1], self.gene_regions[gene2][1]))
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
        while len(self.scores) > 1:
            best_gene_pair = max(self.scores, key=self.scores.get)
            if self.scores[best_gene_pair] <= 0.0:
                break
            if best_gene_pair[0] in self.gene_info.gene_strands:
                assert best_gene_pair[1] not in self.gene_info.gene_strands
                self.merge_genes(best_gene_pair[0], best_gene_pair[1])
            else:
                assert best_gene_pair[0] not in self.gene_info.gene_strands
                self.merge_genes(best_gene_pair[1], best_gene_pair[0])

        transcript_to_new_gene_id = {}
        for gene_id, t_list in self.gene_to_transcripts.items():
            for transcript_id in t_list:
                transcript_to_new_gene_id[transcript_id] = gene_id
        for model in self.transcipt_model_storage:
            model.gene_id = transcript_to_new_gene_id[model.transcript_id]

        return self.transcipt_model_storage
