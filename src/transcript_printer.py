############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import os
from collections import defaultdict
from collections import namedtuple
import gzip

from .common import max_range
from .gene_info import TranscriptModel, GeneInfo

logger = logging.getLogger('IsoQuant')


def validate_exons(novel_exons):
    return novel_exons == sorted(novel_exons) and all(0 < x[0] <= x[1] for x in novel_exons)


class VoidTranscriptPrinter:
    def dump(self, transcript_model_constructor, transcript_model_storage=None):
        pass


class GFFPrinter:
    exon_id_dict = {}

    def __init__(self, outf_prefix, sample_name, exon_id_distributor,
                 gtf_suffix = ".transcript_models.gtf",
                 r2t_suffix = ".transcript_model_reads.tsv",
                 output_r2t = True,
                 check_canonical = False,
                 header = "",
                 gzipped = False):
        self.model_fname = os.path.join(outf_prefix, sample_name + gtf_suffix)
        self.exon_id_distributor = exon_id_distributor
        self.printed_gene_ids = set()
        self.out_gff = open(self.model_fname, "w")
        if header:
            self.out_gff.write("# " + sample_name + " IsoQuant generated GTF\n" + header)
            self.out_gff.flush()

        self.output_r2t = output_r2t
        if self.output_r2t:
            self.r2t_fname = os.path.join(outf_prefix, sample_name + r2t_suffix)
            if gzipped:
                self.out_r2t = gzip.open(self.r2t_fname + ".gz", "wt")
            else:
                self.out_r2t = open(self.r2t_fname, "w")
            if header:
                self.out_r2t.write("#read_id\ttranscript_id\n")

        self.check_canonical = check_canonical

    def __del__(self):
        self.out_gff.close()
        if self.output_r2t:
            self.out_r2t.close()

    def dump(self, gene_info, transcript_model_storage):
        if not transcript_model_storage:
            return

        gene_to_model_dict = defaultdict(list)
        gene_regions = {}
        if not gene_info.empty():
            gene_regions = gene_info.get_gene_regions()
        GFFGeneInfo = namedtuple("GFFGeneInfo", ("chr_id", "strand", "gene_region"))
        gene_info_dict = {}

        for i, model in enumerate(transcript_model_storage):
            if not validate_exons(model.exon_blocks):
                logger.warning("Transcript model %s has incorrect coordinates and will be ignored: %s" %
                               (model.transcript_id, str(model.exon_blocks)))
                continue
            gene_id = model.gene_id
            gene_to_model_dict[gene_id].append(i)

            transcript_region = (model.exon_blocks[0][0], model.exon_blocks[-1][1])
            if gene_id not in gene_info_dict:
                assert model.chr_id == gene_info.chr_id
                gene_range = max_range(gene_regions[gene_id], transcript_region) if gene_id in gene_regions else transcript_region
                gene_info_dict[gene_id] = GFFGeneInfo(model.chr_id, model.strand, gene_range)
            else:
                gene_record = gene_info_dict[gene_id]
                assert model.chr_id == gene_record.chr_id
                if model.strand != gene_record.strand:
                    logger.warning("Gene and transcript records have unequal strands: %s: %s, %s: %s" %
                                   (gene_id, gene_record.strand, model.transcript_id, model.strand))
                gene_info_dict[gene_id] = GFFGeneInfo(model.chr_id, model.strand,
                                                      max_range(gene_record.gene_region, transcript_region))

        gene_order = sorted([(g, gene_info_dict[g].gene_region) for g in gene_info_dict.keys()], key=lambda x:x[1])

        for gene_id, coords in gene_order:
            if gene_id not in self.printed_gene_ids:
                gene_additiional_info = ""
                if gene_info and gene_id in gene_info.gene_attributes:
                    gene_additiional_info = gene_info.gene_attributes[gene_id]
                source = "IsoQuant"
                if gene_info and gene_id in gene_info.sources:
                    source = gene_info.sources[gene_id]
                gene_line = '%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcripts "%d"; %s\n' % \
                        (gene_info_dict[gene_id].chr_id, source, coords[0], coords[1], gene_info_dict[gene_id].strand,
                         gene_id, len(gene_to_model_dict[gene_id]), gene_additiional_info)
                self.out_gff.write(gene_line)
                self.printed_gene_ids.add(gene_id)

            for model_index in gene_to_model_dict[gene_id]:
                model = transcript_model_storage[model_index]
                assert model.gene_id == gene_id

                if not model.check_additional("exons"):
                    model.add_additional_attribute("exons", str(len(model.exon_blocks)))
                transcript_additiional_info = ""
                if gene_info and model.transcript_id in gene_info.gene_attributes:
                    transcript_additiional_info = " " + gene_info.gene_attributes[model.transcript_id]

                transcript_line = '%s\t%s\ttranscript\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; %s\n' \
                                  % (model.chr_id,  model.source, model.exon_blocks[0][0], model.exon_blocks[-1][1],
                                     model.strand, model.gene_id, model.transcript_id,
                                     model.additional_attributes_str() + transcript_additiional_info)
                self.out_gff.write(transcript_line)

                prefix_columns = "%s\t%s\t" % (model.chr_id, model.source)
                suffix_columns = '.\t%s\t.\tgene_id "%s"; transcript_id "%s";' % \
                                 (model.strand, model.gene_id, model.transcript_id)

                exons_to_print = []
                exons_to_print += model.other_features
                for e in model.exon_blocks:
                    exons_to_print.append((e[0], e[1], 'exon'))
                exons_to_print = sorted(exons_to_print, reverse=True) if model.strand == '-' else sorted(exons_to_print)
                for i, e in enumerate(exons_to_print):
                    exon_tuple = (model.chr_id, e[0], e[1], model.strand)
                    if exon_tuple not in GFFPrinter.exon_id_dict:
                        exon_id = self.exon_id_distributor.increment()
                        GFFPrinter.exon_id_dict[exon_tuple] = exon_id
                    else:
                        exon_id = GFFPrinter.exon_id_dict[exon_tuple]
                    exon_str_id = model.chr_id + ".%d" % exon_id
                    feature_type = e[2]
                    self.out_gff.write(prefix_columns + "%s\t%d\t%d\t" % (feature_type, e[0], e[1]) + suffix_columns +
                                       ' exon "%d"; exon_id "%s";\n' % ((i + 1), exon_str_id))
        self.out_gff.flush()

    def dump_read_assignments(self, transcript_model_constructor):
        # write read_id -> transcript_id map
        if not self.output_r2t:
            return
        for model_id in transcript_model_constructor.transcript_read_ids.keys():
            read_assignments = transcript_model_constructor.transcript_read_ids[model_id]
            for a in read_assignments:
                self.out_r2t.write("%s\t%s\n" % (a.read_id, model_id))
        for read_id in transcript_model_constructor.read_assignment_counts.keys():
            if transcript_model_constructor.read_assignment_counts[read_id] == 0:
                self.out_r2t.write("%s\t%s\n" % (read_id, "*"))


def create_extended_storage(genedb, chr_id, chr_record, novel_model_storage):
    all_models = []
    gene_list = list(genedb.region(seqid=chr_id, start=1, featuretype="gene"))
    if not gene_list:
        for m in novel_model_storage:
            all_models.append(m)
        return all_models, GeneInfo.from_region(chr_id, 1, len(chr_record), chr_record=chr_record)
    gene_info = GeneInfo(gene_list, genedb, prepare_profiles=False)
    gene_info.set_reference_sequence(1, len(chr_record), chr_record)
    for isoform_id in gene_info.all_isoforms_exons.keys():
        all_models.append(TranscriptModel.from_reference_transcript(gene_info, isoform_id))
    for m in novel_model_storage:
        all_models.append(m)

    return all_models, gene_info
