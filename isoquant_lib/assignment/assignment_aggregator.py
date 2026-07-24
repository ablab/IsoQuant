############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging

from isoquant_lib.utils.stats import EnumStats
from isoquant_lib.quantification.long_read_counter import (
    ExonCounter,
    JointExonCounter,
    ExonSpliceSiteCounter,
    IntronCounter,
    IntronRetentionCounter,
    CompositeCounter,
    create_gene_counter,
    create_transcript_counter,
)
from isoquant_lib.terminal_prediction.terminal_counter import PolyACounter, TSSCounter
from .assignment_io import (
    IOSupport,
    BEDPrinter,
    ReadAssignmentCompositePrinter,
    ReadInfoPrinter,
    SqantiTSVPrinter,
    BasicTSVAssignmentPrinter,
)
from isoquant_lib.model_construction.transcript_printer import VoidTranscriptPrinter
from isoquant_lib.gene_info import get_all_chromosome_genes, get_all_chromosome_transcripts
from isoquant_lib.common import large_output_enabled

logger = logging.getLogger('IsoQuant')


class ReadAssignmentAggregator:
    def __init__(self, args, sample, string_pools, gffutils_db=None, chr_id=None, gzipped=False, grouping_strategy_names=None):
        self.args = args
        self.string_pools = string_pools
        self.grouping_strategy_names = grouping_strategy_names if grouping_strategy_names else ["default"]
        self.common_header = "# Command line: " + args._cmd_line + "\n# IsoQuant version: " + args._version + "\n"
        self.io_support = IOSupport(self.args)

        self.gene_set = set()
        self.transcript_set = set()
        if gffutils_db and chr_id:
            self.gene_set = set(get_all_chromosome_genes(gffutils_db, chr_id))
            self.transcript_set = set(get_all_chromosome_transcripts(gffutils_db, chr_id))

        self.read_stat_counter = EnumStats()

        self._init_printers(sample, chr_id, gzipped)
        self._init_ungrouped_counters(sample, chr_id)
        self._init_grouped_counters(sample, chr_id)
        self._init_grouped_model_counters(sample, chr_id)

    # ------------------------------------------------------------------ printers

    def _init_printers(self, sample, chr_id, gzipped):
        printer_list = []
        self.corrected_bed_printer = self._make_corrected_bed_printer(sample, chr_id, gzipped, printer_list)
        self.read_info_printer = self._make_read_info_printer(sample, chr_id, gzipped, printer_list)
        self.basic_printer = self._make_basic_printer(sample, chr_id, gzipped, printer_list)
        self.t2t_sqanti_printer = self._make_t2t_printer(sample, chr_id)
        self.global_printer = ReadAssignmentCompositePrinter(printer_list)

    def _make_corrected_bed_printer(self, sample, chr_id, gzipped, printer_list):
        if not large_output_enabled(self.args, "corrected_bed"):
            return None
        corrected_bed_path = sample.get_corrected_bed_file(chr_id) if chr_id else sample.out_corrected_bed
        printer = BEDPrinter(corrected_bed_path,
                             self.args,
                             print_corrected=True,
                             gzipped=gzipped)
        printer_list.append(printer)
        return printer

    def _make_read_info_printer(self, sample, chr_id, gzipped, printer_list):
        if not (self.args.genedb and large_output_enabled(self.args, "read_info")):
            return None
        read_info_path = sample.get_read_info_tsv_file(chr_id) if chr_id else sample.out_read_info_tsv
        printer = ReadInfoPrinter(read_info_path, self.args, self.io_support,
                                  additional_header=self.common_header, gzipped=gzipped)
        printer_list.append(printer)
        return printer

    def _make_basic_printer(self, sample, chr_id, gzipped, printer_list):
        if not (self.args.genedb and large_output_enabled(self.args, "read_assignments")):
            return None
        assigned_tsv_path = sample.get_assigned_tsv_file(chr_id) if chr_id else sample.out_assigned_tsv
        printer = BasicTSVAssignmentPrinter(assigned_tsv_path, self.args, self.io_support,
                                            additional_header=self.common_header, gzipped=gzipped)
        sample.out_assigned_tsv_result = printer.output_file_name
        printer_list.append(printer)
        return printer

    def _make_t2t_printer(self, sample, chr_id):
        if not self.args.sqanti_output:
            return VoidTranscriptPrinter()
        t2t_path = sample.get_t2t_tsv_file(chr_id) if chr_id else sample.out_t2t_tsv
        return SqantiTSVPrinter(t2t_path, self.args, self.io_support)

    # --------------------------------------------------------- ungrouped counters

    def _init_ungrouped_counters(self, sample, chr_id):
        self.global_counter = CompositeCounter()
        self._add_ungrouped_quant_counters(sample, chr_id)

        self.transcript_model_global_counter = CompositeCounter()
        self.gene_model_global_counter = CompositeCounter()
        self._add_ungrouped_model_counters(sample, chr_id)

        self._add_ungrouped_exon_counters(sample, chr_id)
        self._add_ungrouped_terminal_counters(sample, chr_id)
        self._add_ungrouped_ir_counter(sample, chr_id)

    def _add_ungrouped_quant_counters(self, sample, chr_id):
        if not (self.args.genedb and self.args.run_quantification):
            return
        gene_counts_path = sample.get_gene_counts_file(chr_id) if chr_id else sample.out_gene_counts_tsv
        transcript_counts_path = sample.get_transcript_counts_file(chr_id) if chr_id else sample.out_transcript_counts_tsv
        self.gene_counter = create_gene_counter(gene_counts_path,
                                                self.args.gene_quantification,
                                                complete_feature_list=self.gene_set)
        self.transcript_counter = create_transcript_counter(transcript_counts_path,
                                                            self.args.transcript_quantification,
                                                            complete_feature_list=self.transcript_set)
        self.global_counter.add_counters([self.gene_counter, self.transcript_counter])

    def _add_ungrouped_model_counters(self, sample, chr_id):
        if self.args.no_model_construction:
            return
        transcript_model_counts_path = sample.get_transcript_model_counts_file(chr_id) if chr_id else sample.out_transcript_model_counts_tsv
        gene_model_counts_path = sample.get_gene_model_counts_file(chr_id) if chr_id else sample.out_gene_model_counts_tsv
        self.transcript_model_counter = create_transcript_counter(transcript_model_counts_path,
                                                                  self.args.transcript_quantification)
        self.gene_model_counter = create_gene_counter(gene_model_counts_path,
                                                      self.args.gene_quantification)

        self.transcript_model_global_counter.add_counter(self.transcript_model_counter)
        self.gene_model_global_counter.add_counter(self.gene_model_counter)

    def _add_ungrouped_exon_counters(self, sample, chr_id):
        if not (self.args.count_exons and self.args.genedb):
            return
        exon_counts_path = sample.get_exon_counts_file(chr_id) if chr_id else sample.out_exon_counts_tsv
        intron_counts_path = sample.get_intron_counts_file(chr_id) if chr_id else sample.out_intron_counts_tsv
        exon_splice_site_counts_path = (sample.get_exon_splice_site_counts_file(chr_id)
                                        if chr_id else sample.out_exon_splice_site_counts_tsv)
        # string_pools=None means ungrouped counting
        # region-based counts are the default "exon" output now
        self.exon_counter = JointExonCounter(exon_counts_path)
        self.intron_counter = IntronCounter(intron_counts_path)
        self.exon_splice_site_counter = ExonSpliceSiteCounter(
            exon_splice_site_counts_path,
            delta=self.args.delta,
            emit_read_ids=getattr(self.args, "emit_read_ids", False))
        ungrouped_exon_counters = [self.exon_counter, self.intron_counter,
                                   self.exon_splice_site_counter]
        # legacy per-exon inclusion/exclusion counts (deprecated) only on request
        if getattr(self.args, "old_exon_count_format", False):
            old_exon_counts_path = (sample.get_old_exon_counts_file(chr_id)
                                    if chr_id else sample.out_old_exon_counts_tsv)
            self.old_exon_counter = ExonCounter(old_exon_counts_path)
            ungrouped_exon_counters.append(self.old_exon_counter)
        self.global_counter.add_counters(ungrouped_exon_counters)

    def _add_ungrouped_terminal_counters(self, sample, chr_id):
        # polyA / TSS terminal-position prediction (ungrouped). Runs as part of
        # quantification and is also required for model construction; gated by the
        # gene annotation. TSS also requires --fl_data because read start
        # coordinates without full-length evidence are unreliable.
        if not self.args.predict_terminal_sites:
            return
        polya_path = sample.get_polya_prediction_file(chr_id) if chr_id else sample.out_polya_prediction_tsv
        self.polya_counter = PolyACounter(self.args, polya_path)
        self.global_counter.add_counter(self.polya_counter)
        if self.args.fl_data:
            tss_path = sample.get_tss_prediction_file(chr_id) if chr_id else sample.out_tss_prediction_tsv
            self.tss_counter = TSSCounter(self.args, tss_path)
            self.global_counter.add_counter(self.tss_counter)

    def _add_ungrouped_ir_counter(self, sample, chr_id):
        if not (self.args.count_intron_retentions and self.args.genedb):
            return
        ir_counts_path = sample.get_intron_retention_counts_file(chr_id) if chr_id else sample.out_intron_retention_counts_tsv
        self.intron_retention_counter = IntronRetentionCounter(ir_counts_path)
        self.global_counter.add_counter(self.intron_retention_counter)

    # ----------------------------------------------------------- grouped counters

    def _init_grouped_counters(self, sample, chr_id):
        if not (self.args.read_group and self.args.genedb):
            return
        for group_idx, strategy_name in enumerate(self.grouping_strategy_names):
            if self.args.run_quantification:
                self._add_grouped_quant_counters(sample, chr_id, group_idx, strategy_name)
            if self.args.count_exons:
                self._add_grouped_exon_counters(sample, chr_id, group_idx, strategy_name)
            self._add_grouped_terminal_counters(sample, chr_id, group_idx, strategy_name)
            if self.args.count_intron_retentions:
                self._add_grouped_ir_counter(sample, chr_id, group_idx, strategy_name)

    def _add_grouped_quant_counters(self, sample, chr_id, group_idx, strategy_name):
        if chr_id:
            gene_out_file = sample.get_grouped_counts_file(chr_id, "gene", strategy_name)
            transcript_out_file = sample.get_grouped_counts_file(chr_id, "transcript", strategy_name)
        else:
            gene_out_file = f"{sample.out_gene_grouped_counts_tsv}_{strategy_name}"
            transcript_out_file = f"{sample.out_transcript_grouped_counts_tsv}_{strategy_name}"

        gene_counter = create_gene_counter(gene_out_file,
                                           self.args.gene_quantification,
                                           complete_feature_list=self.gene_set,
                                           string_pools=self.string_pools,
                                           group_index=group_idx)
        transcript_counter = create_transcript_counter(transcript_out_file,
                                                       self.args.transcript_quantification,
                                                       complete_feature_list=self.transcript_set,
                                                       string_pools=self.string_pools,
                                                       group_index=group_idx)

        self.global_counter.add_counters([gene_counter, transcript_counter])

    def _add_grouped_exon_counters(self, sample, chr_id, group_idx, strategy_name):
        if chr_id:
            exon_out_file = sample.get_grouped_counts_file(chr_id, "exon", strategy_name)
            intron_out_file = sample.get_grouped_counts_file(chr_id, "splice_junction", strategy_name)
            exon_splice_site_out_file = sample.get_grouped_counts_file(chr_id, "exon_splice_site", strategy_name)
        else:
            exon_out_file = f"{sample.out_exon_grouped_counts_tsv}_{strategy_name}"
            intron_out_file = f"{sample.out_intron_grouped_counts_tsv}_{strategy_name}"
            exon_splice_site_out_file = f"{sample.out_exon_splice_site_grouped_counts_tsv}_{strategy_name}"
        # region-based counts are the default "exon" output now
        exon_counter = JointExonCounter(exon_out_file,
                                        string_pools=self.string_pools, group_index=group_idx)
        intron_counter = IntronCounter(intron_out_file, string_pools=self.string_pools, group_index=group_idx)
        exon_splice_site_counter = ExonSpliceSiteCounter(
            exon_splice_site_out_file, string_pools=self.string_pools, group_index=group_idx,
            delta=self.args.delta, emit_read_ids=getattr(self.args, "emit_read_ids", False))
        grouped_exon_counters = [exon_counter, intron_counter, exon_splice_site_counter]
        if getattr(self.args, "old_exon_count_format", False):
            if chr_id:
                old_exon_out_file = sample.get_grouped_counts_file(chr_id, "old_exon", strategy_name)
            else:
                old_exon_out_file = f"{sample.out_old_exon_grouped_counts_tsv}_{strategy_name}"
            old_exon_counter = ExonCounter(old_exon_out_file,
                                           string_pools=self.string_pools, group_index=group_idx)
            grouped_exon_counters.append(old_exon_counter)
        self.global_counter.add_counters(grouped_exon_counters)

    def _add_grouped_terminal_counters(self, sample, chr_id, group_idx, strategy_name):
        # Grouped polyA / TSS prediction (one file per grouping strategy).
        # Skipped in training-collection mode -- only the ungrouped counter
        # writes the per-chr training fragment, and grouped predictions
        # are not produced in dev mode.
        if self.args.predict_terminal_sites and not getattr(self.args, "collect_polya_training", None):
            if chr_id:
                polya_out_file = sample.get_grouped_counts_file(chr_id, "polyA_prediction", strategy_name) + ".tsv"
            else:
                polya_out_file = f"{sample.out_polya_prediction_grouped_tsv}_{strategy_name}.tsv"
            grouped_polya_counter = PolyACounter(self.args, polya_out_file,
                                                 string_pools=self.string_pools,
                                                 group_index=group_idx)
            self.global_counter.add_counter(grouped_polya_counter)
        if self.args.predict_terminal_sites and self.args.fl_data and not getattr(self.args, "collect_tss_training", None):
            if chr_id:
                tss_out_file = sample.get_grouped_counts_file(chr_id, "TSS_prediction", strategy_name) + ".tsv"
            else:
                tss_out_file = f"{sample.out_tss_prediction_grouped_tsv}_{strategy_name}.tsv"
            grouped_tss_counter = TSSCounter(self.args, tss_out_file,
                                             string_pools=self.string_pools,
                                             group_index=group_idx)
            self.global_counter.add_counter(grouped_tss_counter)

    def _add_grouped_ir_counter(self, sample, chr_id, group_idx, strategy_name):
        if chr_id:
            ir_out_file = sample.get_grouped_counts_file(chr_id, "intron_retention", strategy_name)
        else:
            ir_out_file = f"{sample.out_intron_retention_grouped_counts_tsv}_{strategy_name}"
        ir_counter = IntronRetentionCounter(ir_out_file, string_pools=self.string_pools, group_index=group_idx)
        self.global_counter.add_counter(ir_counter)

    # ----------------------------------------------------- grouped model counters

    def _init_grouped_model_counters(self, sample, chr_id):
        if not (self.args.read_group and not self.args.no_model_construction):
            return
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
                string_pools=self.string_pools,
                group_index=group_idx)
            gene_model_counter = create_gene_counter(
                gene_model_out_file,
                self.args.gene_quantification,
                string_pools=self.string_pools,
                group_index=group_idx)

            self.transcript_model_global_counter.add_counter(transcript_model_counter)
            self.gene_model_global_counter.add_counter(gene_model_counter)
