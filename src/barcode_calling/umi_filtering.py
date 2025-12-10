############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import pysam
import sys
import gzip
import os.path
from collections import defaultdict
import logging
import editdistance
import gffutils

from src.assignment_loader import create_assignment_loader, create_merging_assignment_loader
from src.isoform_assignment import MatchEventSubtype

logger = logging.getLogger('IsoQuant')


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def junctions_from_blocks(sorted_blocks):
    junctions = []
    if len(sorted_blocks) >= 2:
        for i in range(0, len(sorted_blocks) - 1):
            if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
                junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


def load_barcodes(in_file, use_untrusted_umis=False, barcode_column=1, umi_column=2,
                  barcode_score_column=3, umi_property_column=4, min_score=13):
    # Old format:
    # afaf0413-4dd7-491a-928b-39da40d68fb3    99      56      66      83      +       GCCGATACGCCAAT  CGACTGAAG       13      True
    in_files = in_file if isinstance(in_file, list) else [in_file]
    barcode_dict = {}
    read_count = 0
    barcoded = 0
    hq_barcoded = 0
    trusted_umi = 0
    for f in in_files:
        logger.info("Loading barcodes from " + f)
        for l in open(f):
            if l.startswith("#"):continue
            read_count += 1
            v = l.strip().split("\t")
            if len(v) < 5: continue
            barcode = v[barcode_column]
            score = int(v[barcode_score_column])
            if barcode == "*": continue
            barcoded += 1
            if score < min_score:
                continue
            hq_barcoded += 1
            if v[umi_property_column] != "True":
                if not use_untrusted_umis:
                    continue
                umi = "None"
            else:
                umi = v[umi_column]
                trusted_umi += 1

            barcode_dict[v[0]] = (barcode, umi)

    logger.info("Total reads: %d" % read_count)
    logger.info("Barcoded: %d" % barcoded)
    logger.info("Barcoded with score >= %d: %d" % (min_score, hq_barcoded))
    logger.info("Barcoded with trusted UMI: %d" % trusted_umi)
    logger.info("Loaded %d barcodes " % len(barcode_dict))

    return barcode_dict


class ShortReadAssignmentInfo:
    def __init__(self, gene_id, exon_blocks, assignment_type, matching_events, barcode):
        self.gene_id = gene_id
        self.exon_blocks = exon_blocks
        self.assignment_type = assignment_type
        self.matching_events = matching_events
        self.barcode = barcode


class ReadAssignmentInfo:
    def __init__(self, read_id, chr_id, gene_id, transcript_id, strand, exon_blocks, assignment_type,
                 matching_events, barcode, umi, polya_site=-1, transcript_type = "unknown", cell_type = "None"):
        self.read_id = read_id
        self.chr_id = chr_id
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.strand = strand
        self.polya_site = polya_site
        self.exon_blocks = exon_blocks
        self.assignment_type = assignment_type
        self.matching_events = matching_events
        self.barcode = barcode
        self.umi = umi
        self.transcript_type = transcript_type
        self.cell_type = cell_type

    def short(self):
        return ShortReadAssignmentInfo(self.gene_id, self.exon_blocks, self.assignment_type,
                                       self.matching_events, self.barcode)

    def to_allinfo_str(self):
        intron_blocks = junctions_from_blocks(self.exon_blocks)
        exons_str = ";%;" + ";%;".join(["%s_%d_%d_%s" % (self.chr_id, e[0], e[1], self.strand) for e in self.exon_blocks])
        introns_str = ";%;" + ";%;".join(["%s_%d_%d_%s" % (self.chr_id, e[0], e[1], self.strand) for e in intron_blocks])

        if self.assignment_type.is_unique():
            read_type = "known"
        elif self.assignment_type.is_inconsistent():
            read_type = "novel"
        elif self.assignment_type.is_ambiguous():
            read_type = "known_ambiguous"
        else:
            read_type = "none"

        polyA = "NoPolyA"
        TSS = "NoTSS"

        if isinstance(self.matching_events, str):
            if "tss_match" in self.matching_events:
                tss_pos = self.exon_blocks[-1][1] if self.strand == "-" else self.exon_blocks[0][0]
                TSS = "%s_%d_%d_%s" % (self.chr_id, tss_pos, tss_pos, self.strand)
            if "correct_polya" in self.matching_events and self.polya_site != -1: # and self.assignment_type.startswith("unique"):
                polyA_pos = self.polya_site
                polyA = "%s_%d_%d_%s" % (self.chr_id, polyA_pos, polyA_pos, self.strand)
        else:
            if self.strand == '+':
                if any(x in [MatchEventSubtype.terminal_site_match_left,
                             MatchEventSubtype.terminal_site_match_left_precise] for x in self.matching_events):
                    tss_pos = self.exon_blocks[0][0]
                    TSS = "%s_%d_%d_%s" % (self.chr_id, tss_pos, tss_pos, self.strand)
                if (self.polya_site != -1 and
                        any(x == MatchEventSubtype.correct_polya_site_right for x in self.matching_events)):
                    polyA_pos = self.polya_site
                    polyA = "%s_%d_%d_%s" % (self.chr_id, polyA_pos, polyA_pos, self.strand)
            elif self.strand == '-':
                if any(x in [MatchEventSubtype.terminal_site_match_right,
                             MatchEventSubtype.terminal_site_match_right_precise] for x in self.matching_events):
                    tss_pos = self.exon_blocks[-1][1]
                    TSS = "%s_%d_%d_%s" % (self.chr_id, tss_pos, tss_pos, self.strand)
                if (self.polya_site != -1 and
                        any(x == MatchEventSubtype.correct_polya_site_left for x in self.matching_events)):
                    polyA_pos = self.polya_site
                    polyA = "%s_%d_%d_%s" % (self.chr_id, polyA_pos, polyA_pos, self.strand)

        return  "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s" % (self.read_id, self.gene_id, self.cell_type,
                                                                        self.barcode, self.umi, introns_str, TSS, polyA,
                                                                        exons_str, read_type, len(intron_blocks),
                                                                        self.transcript_id, self.transcript_type)


class UMIFilter:
    def __init__(self, split_barcodes_dict, umi_length=0, edit_distance=3, disregard_length_diff=True,
                 only_unique_assignments=False, only_spliced_reads=False):
        self.umi_length = umi_length
        self.max_edit_distance = edit_distance
        self.disregard_length_diff = disregard_length_diff
        self.only_unique_assignments = only_unique_assignments
        self.only_spliced_reads = only_spliced_reads
        self.split_barcodes_dict = split_barcodes_dict

        self.selected_reads = set()
        self.assigned_to_any_gene = 0
        self.spliced = 0
        self.spliced_and_assigned = 0
        self.stats = defaultdict(int)
        self.unique_gene_barcode = set()
        self.total_assignments = 0

        self.duplicated_molecule_counts = defaultdict(int)

        self.umi_len_dif_func = (lambda x: 0) if self.umi_length == 0 else (lambda x: -abs(self.umi_length - x))

    def _find_similar_umi(self, umi, trusted_umi_list):
        if self.max_edit_distance == -1:
            return None if not trusted_umi_list else trusted_umi_list[0]
        similar_umi = None
        best_dist = 100
        # choose the best UMI among checked
        for occ in trusted_umi_list:
            if self.max_edit_distance == 0:
                if self.disregard_length_diff:
                    similar, ed = occ == umi, 0
                elif len(occ) < len(umi):
                    similar, ed = occ in umi, abs(len(occ) - len(umi))
                else:
                    similar, ed = umi in occ, abs(len(occ) - len(umi))
            elif occ == umi:
                similar, ed = True, 0
            else:
                ed = editdistance.eval(occ, umi)
                if not self.disregard_length_diff:
                    ed -= abs(len(occ) - len(umi))
                similar, ed = ed <= self.max_edit_distance, ed

            if similar and ed < best_dist:
                similar_umi = occ
                best_dist = ed
                if best_dist == 0:
                    break

        return similar_umi

    def _construct_umi_dict(self, molecule_list):
        umi_counter = defaultdict(set)
        for m in molecule_list:
            umi_counter[m.umi].add(m.read_id)

        umi_sorting_keys = {}
        for umi in umi_counter:
            # TODO keep only one
            if umi is None or umi == "None":
                umi_sorting_keys[umi] = (-1, 0, "")
            else:
                umi_sorting_keys[umi] = (len(umi_counter[umi]), self.umi_len_dif_func(len(umi)), umi)
        # process reads with the same UMIs consecutively, otherwise reads with the same UMI may fall into different
        # clusters since the heuristic is greedy
        molecule_list = sorted(molecule_list,
                               key=lambda ml: umi_sorting_keys[ml.umi],
                               reverse=True)

        umi_dict = defaultdict(list)
        trusted_umi_list = []
        for m in molecule_list:
            if m.umi is None or m.umi == "None":
                # collect untrusted UMIs together
                umi_dict["None"].append(m)
                continue

            similar_umi = self._find_similar_umi(m.umi, trusted_umi_list)
            if similar_umi is None:
                # if none is added, add to indexer
                umi_dict[m.umi].append(m)
                trusted_umi_list.append(m.umi)
            else:
                umi_dict[similar_umi].append(m)

        return umi_dict

    def _process_duplicates(self, molecule_list):
        if not molecule_list:
            return []

        if len(molecule_list) == 1:
            self.duplicated_molecule_counts[1] += 1
            logger.debug("Unique " + molecule_list[0].read_id)
            return molecule_list

        resulting_reads = []
        umi_dict = self._construct_umi_dict(molecule_list)
        for umi in umi_dict.keys():
            duplicate_count = len(umi_dict[umi])
            self.duplicated_molecule_counts[duplicate_count] += 1
            if duplicate_count == 1:
                resulting_reads.append((umi_dict[umi][0], umi))
                continue

            best_read = umi_dict[umi][0]
            logger.debug("Selecting from:")
            for m in umi_dict[umi]:
                logger.debug("%s %s" % (m.read_id, m.umi))
                if not best_read.assignment_type.is_unique() and m.assignment_type.is_unique():
                    best_read = m
                elif len(m.exon_blocks) > len(best_read.exon_blocks):
                    best_read = m
                elif len(m.exon_blocks) == len(best_read.exon_blocks) and \
                        m.exon_blocks[-1][1] - m.exon_blocks[0][0] > \
                        best_read.exon_blocks[-1][1] - best_read.exon_blocks[0][0]:
                    best_read = m

            polyas = set()
            transcript_types = set()
            isoform_ids = set()
            for m in umi_dict[umi]:
                if m.read_id != best_read.read_id:
                    continue
                transcript_types.add(m.transcript_type)
                polyas.add(m.polya_site)
                isoform_ids.add(m.transcript_id)

            if len(transcript_types) > 1:
                best_read.transcript_type = "None"
            if len(polyas) > 1:
                best_read.polya_site = -1
            if len(isoform_ids) > 1:
                best_read.transcript_id = "None"

            if best_read.assignment_type.is_inconsistent():
                best_read.transcript_id = "None"
                best_read.polya_site = -1
                best_read.transcript_type = "None"

            logger.debug("Selected %s %s" % (best_read.read_id, best_read.umi))
            resulting_reads.append((best_read, umi))

        if len(resulting_reads) == 1:
            # if we have a single UMI - we do not care whether it's trusted or not
            return [resulting_reads[0][0]]
        # if we have > 1 UMIs, we ignore untrusted ones
        return [x[0] for x in filter(lambda x: x[1] != "None", resulting_reads)]

    def _process_gene(self, gene_dict):
        for barcode in gene_dict:
           for r in self._process_duplicates(gene_dict[barcode]):
               yield r

    def _process_chunk(self, gene_barcode_dict, allinfo_outf, read_ids_outf=None):
        read_count = 0
        spliced_count = 0
        for gene_id in gene_barcode_dict:
            for read_assignment in self._process_gene(gene_barcode_dict[gene_id]):
                if (not read_assignment.assignment_type.is_unique() and
                        read_assignment.read_id in self.selected_reads):
                    continue
                read_count += 1
                self.unique_gene_barcode.add((read_assignment.gene_id, read_assignment.barcode))
                if len(read_assignment.exon_blocks) > 1:
                    spliced_count += 1
                if read_ids_outf:
                    read_ids_outf.write(read_assignment.read_id + "\n")
                self.selected_reads.add(read_assignment.read_id)
                if allinfo_outf:
                    allinfo_outf.write(read_assignment.to_allinfo_str() + "\n")
        return read_count, spliced_count

    def add_stats_for_read(self, read_info):
        assigned = read_info.gene_id != "."
        spliced = len(read_info.exon_blocks) > 1
        barcoded = read_info.barcode is not None
        unique = assigned

        if assigned:
            self.stats["Assigned to any gene"] += 1
        if spliced:
            self.stats["Spliced"] += 1
        if unique:
            self.stats["Uniquely assigned"] += 1
        if unique and spliced:
            self.stats["Uniquely assigned and spliced"] += 1
        if barcoded:
            if assigned:
                self.stats["Assigned to any gene and barcoded"] += 1
            if spliced:
                self.stats["Spliced and barcoded"] += 1
            if unique:
                self.stats["Uniquely assigned and barcoded"] += 1
            if unique and spliced:
                self.stats["Uniquely assigned and spliced and barcoded"] += 1

    @staticmethod
    def load_barcodes_simple(barcode_file):
        barcode_dict = {}
        for l in open(barcode_file):
            if l.startswith("#"):
                continue
            v = l.split()
            if len(v) != 3: continue
            barcode_dict[v[0]] = (v[1], v[2])
        return barcode_dict

    def process_single_chr(self, args, chr_id, saves_prefix, transcript_type_dict, barcode_feature_table,
                           all_info_file_name, filtered_reads_file_name, stats_output_file_name):

        with open(all_info_file_name, "w") as allinfo_outf:
            filtered_reads_outf = open(filtered_reads_file_name, "w") if filtered_reads_file_name else None
            read_count = 0
            spliced_count = 0
            self.unique_gene_barcode = set()

            barcode_dict = self.load_barcodes_simple(self.split_barcodes_dict[chr_id])
            loader = create_merging_assignment_loader(chr_id, saves_prefix)
            while loader.has_next():
                gene_barcode_dict = defaultdict(lambda: defaultdict(list))
                _, assignment_storage = loader.get_next()
                logger.debug("Processing %d reads" % len(assignment_storage))

                for read_assignment in assignment_storage:
                    if read_assignment.gene_assignment_type.is_unassigned():
                        continue
                    if read_assignment.gene_assignment_type.is_ambiguous():
                        self.stats["Assigned to any gene"] += 1
                        if len(read_assignment.corrected_exons) > 1:
                            self.stats["Spliced"] += 1
                        continue

                    read_id = read_assignment.read_id
                    assignment_type = read_assignment.assignment_type
                    exon_blocks = read_assignment.corrected_exons
                    if read_id in barcode_dict:
                        barcode, umi = barcode_dict[read_id]
                    else:
                        barcode, umi = None, None
                    strand = read_assignment.strand
                    spliced = len(exon_blocks) > 1
                    barcoded = barcode is not None
                    cell_type = "None" if barcode is None or barcode not in barcode_feature_table else barcode_feature_table[barcode]

                    if len(read_assignment.isoform_matches) == 1:
                        isoform_match = read_assignment.isoform_matches[0]
                        if read_assignment.assignment_type.is_consistent():
                            transcript_id = isoform_match.assigned_transcript
                            transcript_type, polya_site = (
                                transcript_type_dict[transcript_id] if transcript_id in transcript_type_dict
                                else ("unknown_type", -1))
                            cell_type = "None" if barcode is None or barcode not in barcode_feature_table else \
                            barcode_feature_table[barcode]
                            assignment_info = ReadAssignmentInfo(read_id, chr_id, isoform_match.assigned_gene,
                                                                 transcript_id, strand,
                                                                 exon_blocks,
                                                                 assignment_type,
                                                                 [e.event_type for e in isoform_match.match_subclassifications],
                                                                 barcode, umi,
                                                                 polya_site, transcript_type, cell_type)
                        else:
                            assignment_info = ReadAssignmentInfo(read_id, chr_id, isoform_match.assigned_gene,
                                                                 isoform_match.assigned_transcript, strand,
                                                                 exon_blocks, assignment_type, [], barcode, umi,
                                                                 -1, "unknown_type", cell_type)
                    elif read_assignment.assignment_type.is_consistent():
                        transcript_types = set()
                        polya_sites = set()
                        gene_ids = set()
                        for m in read_assignment.isoform_matches:
                            transcript_id = m.assigned_transcript
                            gene_ids.add(m.assigned_gene)
                            transcript_type, polya_site = (
                                transcript_type_dict[transcript_id] if transcript_id in transcript_type_dict
                                else ("unknown_type", -1))
                            transcript_types.add(transcript_type)
                            polya_sites.add(polya_site)

                        assert len(gene_ids) == 1, "Multiple genes assigned to a single read"
                        transcript_type = "None" if len(transcript_types) != 1 else transcript_types.pop()
                        polya_site = -1 if len(polya_sites) != 1 else polya_sites.pop()

                        assignment_info = ReadAssignmentInfo(read_id, chr_id, gene_ids.pop(), "None", strand, exon_blocks,
                                                             assignment_type, [], barcode, umi, polya_site,
                                                             transcript_type, cell_type)

                    self.total_assignments += 1
                    self.add_stats_for_read(assignment_info)

                    if not barcoded:
                        continue
                    if not spliced and self.only_spliced_reads:
                        continue
                    gene_barcode_dict[assignment_info.gene_id][barcode].append(assignment_info)

                processed_read_count, processed_spliced_count = self._process_chunk(gene_barcode_dict, allinfo_outf, filtered_reads_outf)
                read_count += processed_read_count
                spliced_count += processed_spliced_count

            if filtered_reads_outf:
                filtered_reads_outf.close()

        with open(stats_output_file_name, "w") as count_hist_file:
            count_hist_file.write("Unique gene-barcodes pairs\t%d\n" % len(self.unique_gene_barcode))
            count_hist_file.write("Total reads saved\t%d\n" % read_count)
            count_hist_file.write("Spliced reads saved\t%d\n" % spliced_count)
            count_hist_file.write("Total assignments processed\t%d\n" % self.total_assignments)
            for k in sorted(self.stats.keys()):
                count_hist_file.write("%s\t%d\n" % (k, self.stats[k]))

        return all_info_file_name, stats_output_file_name


def filter_bam(in_file_name, out_file_name, read_set):
    inf = pysam.AlignmentFile(in_file_name, "rb")
    outf = pysam.AlignmentFile(out_file_name, "wb", template=inf)

    count = 0
    passed = 0

    for read in inf:
        if read.reference_id == -1 or read.is_secondary:
            continue

        count += 1
        if count % 10000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        if read.query_name in read_set:
            outf.write(read)
            passed += 1

    print("Processed " + str(count) + " reads, written " + str(passed))
    inf.close()
    outf.close()
    pysam.index(out_file_name)


def create_transcript_info_dict(genedb, chr_ids=None):
    gffutils_db = gffutils.FeatureDB(genedb)
    transcript_type_dict = {}
    for t in gffutils_db.features_of_type(('transcript', 'mRNA')):
        if chr_ids and t.seqid not in chr_ids:
            continue
        polya_site = t.start - 1 if t.strand == '-' else t.end + 1
        if "transcript_type" in t.attributes.keys():
            transcript_type_dict[t.id] = (t.attributes["transcript_type"][0], polya_site)
        else:
            transcript_type_dict[t.id] = ("unknown_type", polya_site)
    return transcript_type_dict

