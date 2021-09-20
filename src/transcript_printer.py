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

from src.common import *
from src.assignment_io import *
from src.isoform_assignment import *
from src.long_read_profiles import *
from src.junction_comparator import *
from src.long_read_assigner import *
from src.gene_info import *

logger = logging.getLogger('IsoQuant')


class GFFPrinter:
    def __init__(self, outf_prefix, sample_name, io_support):
        self.model_fname = os.path.join(outf_prefix, sample_name + ".transcript_models.gtf")
        self.out_gff = open(self.model_fname, "w")
        self.out_gff.write("# " + sample_name + " IsoQuant generated GFF\n")
        self.out_r2t = open(os.path.join(outf_prefix, sample_name + ".transcript_model_reads.tsv"), "w")
        self.out_r2t.write("#read_id\ttranscript_id\n")
        self.io_support = io_support

    def __del__(self):
        self.out_gff.close()
        self.out_r2t.close()

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
                gene_range = max_range(gene_regions[gene_id], transcript_region) if gene_id in gene_regions else transcript_region
                gene_info_dict[gene_id] = GFFGeneInfo(model.chr_id, model.strand, gene_range)
            else:
                gene_record = gene_info_dict[gene_id]
                assert model.chr_id == gene_record.chr_id
                if model.strand != gene_record.strand:
                    logger.warning("Unequal strands: %s = %s, %s = %s" % (gene_id, gene_record.strand, model.transcript_id, model.strand))
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

                if transcript_model_constructor.params.check_canonical and \
                        transcript_model_constructor.gene_info.reference_region:
                    model_introns = junctions_from_blocks(model.exon_blocks)
                    if len(model_introns) == 0:
                        canonical_info = "Canonical=Unspliced;"
                    else:
                        strand = transcript_model_constructor.gene_info.isoform_strands[model.reference_transcript]
                        all_canonical = self.io_support.check_sites_are_canonical(model_introns,
                                                                                  transcript_model_constructor.gene_info,
                                                                                  strand)
                        canonical_info = "Canonical=" + str(all_canonical) + ";"
                    model.additional_info += canonical_info

                transcript_line = '%s\tIsoQuant\ttranscript\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; ' \
                                  'reference_gene_id "%s"; reference_transcript_id "%s"; %s\n' % \
                            (model.chr_id, model.exon_blocks[0][0], model.exon_blocks[-1][1], model.strand,
                             model.gene_id, model.transcript_id, model.reference_gene, model.reference_transcript,
                             model.additional_info)
                self.out_gff.write(transcript_line)

                prefix_columns = "%s\tIsoQuant\texon\t" % model.chr_id
                suffix_columns = '.\t%s\t.\tgene_id "%s"; transcript_id "%s"; ' \
                                 'reference_gene_id "%s"; reference_transcript_id "%s";\n' % \
                                 (model.strand, model.gene_id, model.transcript_id,
                                  model.reference_gene, model.reference_transcript)
                for e in model.exon_blocks:
                    self.out_gff.write(prefix_columns + "%d\t%d\t" % (e[0], e[1]) + suffix_columns)

        # write read_id -> transcript_id map
        used_reads = set()
        for model_id, reads in transcript_model_constructor.transcript_read_ids.items():
            for read_id in reads:
                used_reads.add(read_id)
                self.out_r2t.write("%s\t%s\n" % (read_id, model_id))
        for read_id in transcript_model_constructor.unused_reads:
            self.out_r2t.write("%s\t%s\n" % (read_id, "*"))
