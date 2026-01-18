############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging

from .stats import EnumStats
from .long_read_counter import (
    ExonCounter,
    IntronCounter,
    CompositeCounter,
    create_gene_counter,
    create_transcript_counter,
)
from .assignment_io import (
    IOSupport,
    BEDPrinter,
    ReadAssignmentCompositePrinter,
    SqantiTSVPrinter,
    BasicTSVAssignmentPrinter,
)
from .transcript_printer import VoidTranscriptPrinter
from .gene_info import get_all_chromosome_genes, get_all_chromosome_transcripts
from .common import large_output_enabled

logger = logging.getLogger('IsoQuant')


class ReadAssignmentAggregator:
    def __init__(self, args, sample, read_groups, gffutils_db=None, chr_id=None, gzipped=False, grouping_strategy_names=None):
        self.args = args
        self.read_groups = read_groups
        self.grouping_strategy_names = grouping_strategy_names if grouping_strategy_names else ["default"]
        self.common_header = "# Command line: " + args._cmd_line + "\n# IsoQuant version: " + args._version + "\n"
        self.io_support = IOSupport(self.args)

        self.gene_set = set()
        self.transcript_set = set()
        if gffutils_db and chr_id:
            self.gene_set = set(get_all_chromosome_genes(gffutils_db, chr_id))
            self.transcript_set = set(get_all_chromosome_transcripts(gffutils_db, chr_id))

        self.read_stat_counter = EnumStats()

        printer_list = []
        self.corrected_bed_printer = None
        if large_output_enabled(self.args, "corrected_bed"):
            corrected_bed_path = sample.get_corrected_bed_file(chr_id) if chr_id else sample.out_corrected_bed
            self.corrected_bed_printer = BEDPrinter(corrected_bed_path,
                                                    self.args,
                                                    print_corrected=True,
                                                    gzipped=gzipped)
            printer_list.append(self.corrected_bed_printer)
        self.basic_printer = None
        if self.args.genedb and large_output_enabled(self.args, "read_assignments"):
            assigned_tsv_path = sample.get_assigned_tsv_file(chr_id) if chr_id else sample.out_assigned_tsv
            self.basic_printer = BasicTSVAssignmentPrinter(assigned_tsv_path, self.args, self.io_support,
                                                           additional_header=self.common_header, gzipped=gzipped)
            sample.out_assigned_tsv_result = self.basic_printer.output_file_name
            printer_list.append(self.basic_printer)
        self.t2t_sqanti_printer = VoidTranscriptPrinter()
        if self.args.sqanti_output:
            t2t_path = sample.get_t2t_tsv_file(chr_id) if chr_id else sample.out_t2t_tsv
            self.t2t_sqanti_printer = SqantiTSVPrinter(t2t_path, self.args, self.io_support)
        self.global_printer = ReadAssignmentCompositePrinter(printer_list)

        self.global_counter = CompositeCounter()
        if self.args.genedb:
            gene_counts_path = sample.get_gene_counts_file(chr_id) if chr_id else sample.out_gene_counts_tsv
            transcript_counts_path = sample.get_transcript_counts_file(chr_id) if chr_id else sample.out_transcript_counts_tsv
            self.gene_counter = create_gene_counter(gene_counts_path,
                                                    self.args.gene_quantification,
                                                    complete_feature_list=self.gene_set)
            self.transcript_counter = create_transcript_counter(transcript_counts_path,
                                                                self.args.transcript_quantification,
                                                                complete_feature_list=self.transcript_set)
            self.global_counter.add_counters([self.gene_counter, self.transcript_counter])

        self.transcript_model_global_counter = CompositeCounter()
        self.gene_model_global_counter = CompositeCounter()
        if not self.args.no_model_construction:
            transcript_model_counts_path = sample.get_transcript_model_counts_file(chr_id) if chr_id else sample.out_transcript_model_counts_tsv
            gene_model_counts_path = sample.get_gene_model_counts_file(chr_id) if chr_id else sample.out_gene_model_counts_tsv
            self.transcript_model_counter = create_transcript_counter(transcript_model_counts_path,
                                                                      self.args.transcript_quantification)
            self.gene_model_counter = create_gene_counter(gene_model_counts_path,
                                                          self.args.gene_quantification)

            self.transcript_model_global_counter.add_counter(self.transcript_model_counter)
            self.gene_model_global_counter.add_counter(self.gene_model_counter)

        if self.args.count_exons and self.args.genedb:
            exon_counts_path = sample.get_exon_counts_file(chr_id) if chr_id else sample.out_exon_counts_tsv
            intron_counts_path = sample.get_intron_counts_file(chr_id) if chr_id else sample.out_intron_counts_tsv
            self.exon_counter = ExonCounter(exon_counts_path, ignore_read_groups=True)
            self.intron_counter = IntronCounter(intron_counts_path, ignore_read_groups=True)
            self.global_counter.add_counters([self.exon_counter, self.intron_counter])

        if self.args.read_group and self.args.genedb:
            for group_idx, strategy_name in enumerate(self.grouping_strategy_names):
                # Use chr-specific paths if chr_id is provided
                if chr_id:
                    gene_out_file = sample.get_grouped_counts_file(chr_id, "gene", strategy_name)
                    transcript_out_file = sample.get_grouped_counts_file(chr_id, "transcript", strategy_name)
                else:
                    gene_out_file = f"{sample.out_gene_grouped_counts_tsv}_{strategy_name}"
                    transcript_out_file = f"{sample.out_transcript_grouped_counts_tsv}_{strategy_name}"

                gene_counter = create_gene_counter(gene_out_file,
                                                   self.args.gene_quantification,
                                                   complete_feature_list=self.gene_set,
                                                   read_groups=self.read_groups[group_idx],
                                                   group_index=group_idx)
                transcript_counter = create_transcript_counter(transcript_out_file,
                                                              self.args.transcript_quantification,
                                                              complete_feature_list=self.transcript_set,
                                                              read_groups=self.read_groups[group_idx],
                                                              group_index=group_idx)

                self.global_counter.add_counters([gene_counter, transcript_counter])

                if self.args.count_exons:
                    if chr_id:
                        exon_out_file = sample.get_grouped_counts_file(chr_id, "exon", strategy_name)
                        intron_out_file = sample.get_grouped_counts_file(chr_id, "intron", strategy_name)
                    else:
                        exon_out_file = f"{sample.out_exon_grouped_counts_tsv}_{strategy_name}"
                        intron_out_file = f"{sample.out_intron_grouped_counts_tsv}_{strategy_name}"
                    exon_counter = ExonCounter(exon_out_file, group_index=group_idx)
                    intron_counter = IntronCounter(intron_out_file, group_index=group_idx)
                    self.global_counter.add_counters([exon_counter, intron_counter])

        if self.args.read_group and not self.args.no_model_construction:
            for group_idx, strategy_name in enumerate(self.grouping_strategy_names):
                if chr_id:
                    transcript_model_out_file = sample.get_grouped_counts_file(chr_id, "discovered_transcript", strategy_name)
                    gene_model_out_file = sample.get_grouped_counts_file(chr_id, "discovered_gene", strategy_name)
                else:
                    transcript_model_out_file = f"{sample.out_transcript_model_grouped_counts_tsv}_{strategy_name}"
                    gene_model_out_file = f"{sample.out_gene_model_grouped_counts_tsv}_{strategy_name}"

                transcript_model_counter = create_transcript_counter(
                    transcript_model_out_file,
                    self.args.transcript_quantification,
                    read_groups=self.read_groups[group_idx],
                    group_index=group_idx)
                gene_model_counter = create_gene_counter(
                    gene_model_out_file,
                    self.args.gene_quantification,
                    read_groups=self.read_groups[group_idx],
                    group_index=group_idx)

                self.transcript_model_global_counter.add_counter(transcript_model_counter)
                self.gene_model_global_counter.add_counter(gene_model_counter)
