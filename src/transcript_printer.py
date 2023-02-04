############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
from collections import namedtuple

from src.common import *
from src.assignment_io import *
from src.long_read_assigner import *
from src.gene_info import *

logger = logging.getLogger('IsoQuant')


def validate_exons(novel_exons):
    return novel_exons == sorted(novel_exons) and all(0 < x[0] <= x[1] for x in novel_exons)


class GFFPrinter:
    transcript_id_counter = AtomicCounter()
    transcript_id_dict = {}

    def __init__(self, outf_prefix, sample_name, io_support,
                 gtf_suffix = ".transcript_models.gtf",
                 r2t_suffix = ".transcript_model_reads.tsv",
                 output_r2t = True,
                 header = ""):
        self.model_fname = os.path.join(outf_prefix, sample_name + gtf_suffix)
        self.printed_gene_ids = set()
        self.out_gff = open(self.model_fname, "w")
        if header:
            self.out_gff.write("# " + sample_name + " IsoQuant generated GTF\n" + header)
            self.out_gff.flush()

        self.output_r2t = output_r2t
        if self.output_r2t:
            self.r2t_fname = os.path.join(outf_prefix, sample_name + r2t_suffix)
            self.out_r2t = open(self.r2t_fname, "w")
            if header:
                self.out_r2t.write("#read_id\ttranscript_id\n")
                self.out_r2t.flush()
        self.io_support = io_support

    def __del__(self):
        self.out_gff.close()
        if self.output_r2t:
            self.out_r2t.close()

    def dump(self, transcript_model_constructor, transcript_model_storage=None):
        # write exons to GFF
        if transcript_model_storage is None:
            transcript_model_storage = transcript_model_constructor.transcript_model_storage
        gene_to_model_dict = defaultdict(list)
        gene_regions = {}
        if not transcript_model_constructor.gene_info.empty():
            gene_regions = transcript_model_constructor.gene_info.get_gene_regions()
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
                assert model.chr_id == transcript_model_constructor.gene_info.chr_id
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
                if transcript_model_constructor.gene_info and gene_id in transcript_model_constructor.gene_info.gene_attributes:
                    gene_additiional_info = transcript_model_constructor.gene_info.gene_attributes[gene_id]
                gene_line = '%s\tIsoQuant\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcripts "%d"; %s\n' % \
                            (gene_info_dict[gene_id].chr_id, coords[0], coords[1], gene_info_dict[gene_id].strand,
                             gene_id, len(gene_to_model_dict[gene_id]), gene_additiional_info)
                self.out_gff.write(gene_line)
                self.printed_gene_ids.add(gene_id)

            for model_index in gene_to_model_dict[gene_id]:
                model = transcript_model_storage[model_index]
                assert model.gene_id == gene_id

                if not model.check_additional("exons"):
                    model.add_additional_attribute("exons", str(len(model.exon_blocks)))
                if transcript_model_constructor.params.check_canonical and \
                        transcript_model_constructor.gene_info.reference_region:
                    self.add_canonical_info(model, transcript_model_constructor.gene_info)

                transcript_line = '%s\tIsoQuant\ttranscript\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; %s\n' \
                                  % (model.chr_id, model.exon_blocks[0][0], model.exon_blocks[-1][1], model.strand,
                                     model.gene_id, model.transcript_id, model.additional_attributes_str())
                self.out_gff.write(transcript_line)

                prefix_columns = "%s\tIsoQuant\texon\t" % model.chr_id
                suffix_columns = '.\t%s\t.\tgene_id "%s"; transcript_id "%s";' % \
                                 (model.strand, model.gene_id, model.transcript_id)
                exons_to_print = sorted(model.exon_blocks, reverse=True) if model.strand == '-' else model.exon_blocks
                for i, e in enumerate(exons_to_print):
                    exon_tuple = (model.chr_id, e[0], e[1], model.strand)
                    if exon_tuple not in GFFPrinter.transcript_id_dict:
                        exon_id = GFFPrinter.transcript_id_counter.increment()
                        GFFPrinter.transcript_id_dict[exon_tuple] = exon_id
                    else:
                        exon_id = GFFPrinter.transcript_id_dict[exon_tuple]
                    exon_str_id = model.chr_id + ".%d" % exon_id
                    self.out_gff.write(prefix_columns + "%d\t%d\t" % (e[0], e[1]) + suffix_columns +
                                       ' exon "%d"; exon_id "%s";\n' % ((i + 1), exon_str_id))

        # write read_id -> transcript_id map
        if self.output_r2t:
            used_reads = set()
            for model_id, read_assignments in transcript_model_constructor.transcript_read_ids.items():
                for a in read_assignments:
                    used_reads.add(a.read_id)
                    self.out_r2t.write("%s\t%s\n" % (a.read_id, model_id))
            for read_id in transcript_model_constructor.unused_reads:
                self.out_r2t.write("%s\t%s\n" % (read_id, "*"))

    def add_canonical_info(self, model, gene_info):
        key_word = 'Canonical'
        if model.check_additional(key_word):
            return
        model_introns = junctions_from_blocks(model.exon_blocks)
        if len(model_introns) == 0:
            canonical_info = "Unspliced"
        else:
            all_canonical = self.io_support.check_sites_are_canonical(model_introns, gene_info, model.strand)
            canonical_info = str(all_canonical)
        model.add_additional_attribute(key_word, canonical_info)
