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
            if len(v) < 8: continue
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
                 matching_events, barcode, umi, polya_site=-1, transcript_type = "unknown"):
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

    def short(self):
        return ShortReadAssignmentInfo(self.gene_id, self.exon_blocks, self.assignment_type,
                                       self.matching_events, self.barcode)

    def to_allinfo_str(self):
        intron_blocks = junctions_from_blocks(self.exon_blocks)
        exons_str = ";%;" + ";%;".join(["%s_%d_%d_%s" % (self.chr_id, e[0], e[1], self.strand) for e in self.exon_blocks])
        introns_str = ";%;" + ";%;".join(["%s_%d_%d_%s" % (self.chr_id, e[0], e[1], self.strand) for e in intron_blocks])

        cell_type = "None"
        read_type = "known" if self.assignment_type.startswith("unique") else "novel"

        polyA = "NoPolyA"
        TSS = "NoTSS"
        if "tss_match" in self.matching_events:
            tss_pos = self.exon_blocks[-1][1] if self.strand == "-" else self.exon_blocks[0][0]
            TSS = "%s_%d_%d_%s" % (self.chr_id, tss_pos, tss_pos, self.strand)
        if "correct_polya" in self.matching_events:
            polyA_pos = self.polya_site
            polyA = "%s_%d_%d_%s" % (self.chr_id, polyA_pos, polyA_pos, self.strand)
        return  "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s" % (self.read_id, self.gene_id, cell_type,
                                                                        self.barcode, self.umi, introns_str, TSS, polyA,
                                                                        exons_str, read_type, len(intron_blocks),
                                                                        self.transcript_id, self.transcript_type)


class UMIFilter:
    def __init__(self, barcode_umi_dict, edit_distance=3, disregard_length_diff=True,
                 only_unique_assignments=False, only_spliced_reads=False):
        self.max_edit_distance = edit_distance
        self.disregard_length_diff = disregard_length_diff
        self.only_unique_assignments = only_unique_assignments
        self.only_spliced_reads = only_spliced_reads
        self.barcode_dict = barcode_umi_dict

        self.selected_reads = set()
        self.assigned_to_any_gene = 0
        self.spliced = 0
        self.spliced_and_assigned = 0
        self.stats = defaultdict(int)
        self.unique_gene_barcode = set()

        self.total_assignments = 0
        self.duplicated_molecule_counts = defaultdict(int)

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
                    similar, ed =  umi in occ, abs(len(occ) - len(umi))
            else:
                ed = editdistance.eval(occ, umi)
                if not self.disregard_length_diff:
                    ed -= abs(len(occ) - len(umi))
                similar, ed = ed <= self.max_edit_distance, ed

            if similar and ed < best_dist:
                similar_umi = occ
                best_dist = ed

        return similar_umi

    def _process_duplicates(self, molecule_list):
        if not molecule_list:
            return []

        if len(molecule_list) == 1:
            self.duplicated_molecule_counts[1] += 1
            logger.debug("Unique " + molecule_list[0].read_id)
            return molecule_list

        resulting_reads = []
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
                if len(m.exon_blocks) > len(best_read.exon_blocks):
                    best_read = m
                elif len(m.exon_blocks) == len(best_read.exon_blocks) and \
                        m.exon_blocks[-1][1] - m.exon_blocks[0][0] > \
                        best_read.exon_blocks[-1][1] - best_read.exon_blocks[0][0]:
                    best_read = m
            logger.debug("Selected %s %s" % (best_read.read_id, best_read.umi))
            resulting_reads.append((best_read, umi))

        if len(resulting_reads) == 1:
            # if we have a single UMI - we do not care whether it's trusted or not
            return [resulting_reads[0][0]]
        # if we have > 1 UMIs, we ignore untrusted ones
        return [x[0] for x in filter(lambda x: x != "None", resulting_reads)]

    def _process_gene(self, gene_dict):
        resulting_reads = []
        for barcode in gene_dict:
            resulting_reads += self._process_duplicates(gene_dict[barcode])
        return resulting_reads

    def _process_chunk(self, gene_barcode_dict, outf, allinfo_outf=None):
        read_count = 0
        spliced_count = 0
        for gene_id in gene_barcode_dict:
            assignment_list = self._process_gene(gene_barcode_dict[gene_id])
            for read_assignment in assignment_list:
                if (not read_assignment.assignment_type.startswith("unique") and
                        read_assignment.read_id in self.selected_reads):
                    continue
                read_count += 1
                if len(read_assignment.exon_blocks) > 1:
                    spliced_count += 1
                outf.write(read_assignment.read_id + "\n")
                self.selected_reads.add(read_assignment.read_id)
                if allinfo_outf:
                    allinfo_outf.write(read_assignment.to_allinfo_str() + "\n")
        return read_count, spliced_count

    def count_stats_for_storage(self, read_info_storage):
        logger.info("Processing %d reads" % len(read_info_storage))
        for read_id in read_info_storage.keys():
            read_infos = read_info_storage[read_id]
            assigned = read_infos[0].gene_id != "."
            spliced = len(read_infos[0].exon_blocks) > 1
            barcoded = read_infos[0].barcode is not None
            unique = assigned and len(set(r.gene_id for r in read_infos)) == 1

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

            if assigned and barcoded:
                self.unique_gene_barcode.add((read_infos[0].gene_id, read_infos[0].barcode))

    def process(self, assignment_file, output_prefix, transcript_type_dict):
        outf = open(output_prefix + ".read_ids.tsv", "w")
        allinfo_outf = open(output_prefix + ".allinfo", "w")

        # 1251f521-d4c2-4dcc-96d7-85070cc44e12    chr1    -       ENST00000477740.5       ENSG00000238009.6       ambiguous       ism_internal    112699-112804,120721-120759     PolyA=False
        gene_barcode_dict = defaultdict(lambda: defaultdict(list))
        read_info_storage = defaultdict(list)
        current_interval = (-1, -1)
        current_chr = None
        read_count = 0
        spliced_count = 0

        self.unique_gene_barcode = set()
        if os.path.exists(assignment_file):
            handle = open(assignment_file, 'r')
        elif os.path.exists(assignment_file + ".gz"):
            handle = gzip.open(assignment_file + ".gz", "rt")
        else:
            logger.critical("Read assignment file is not found")
            exit(-1)
        for l in handle:
            if l.startswith("#"): continue
            v = l.strip().split("\t")
            read_id = v[0]
            chr_id = v[1]
            gene_id = v[4]
            transcript_id = v[3]
            assignment_type = v[5]
            exon_blocks_str = v[7]
            exon_blocks = list(map(lambda x: tuple(map(int, x.split('-'))), exon_blocks_str.split(',')))
            if read_id in self.barcode_dict:
                barcode, umi = self.barcode_dict[read_id]
            else:
                barcode, umi = None, None
            strand = v[2]
            matching_events = v[6]

            self.total_assignments += 1
            assigned = gene_id != "."
            spliced = len(exon_blocks) > 1
            #unique = assignment_type.startswith("unique")
            barcoded = barcode is not None
            transcript_type, polya_site = (transcript_type_dict[transcript_id] if transcript_id in transcript_type_dict
                                           else ("unknown_type", -1))
            assignment_info = ReadAssignmentInfo(read_id, chr_id, gene_id, transcript_id, strand, exon_blocks,
                                                 assignment_type, matching_events, barcode, umi,
                                                 polya_site, transcript_type)
            read_info_storage[read_id].append(assignment_info.short())

            if not barcoded or not assigned:
                continue

            #if not unique and self.only_unique_assignments:
            #    continue

            if not spliced and self.only_spliced_reads:
                continue

            read_interval = (exon_blocks[0][0], exon_blocks[-1][1])
            if current_chr != chr_id or not overlaps(current_interval, read_interval):
                processed_read_count, processed_spliced_count = self._process_chunk(gene_barcode_dict, outf, allinfo_outf)
                read_count += processed_read_count
                spliced_count += processed_spliced_count
                if current_chr != chr_id:
                    logger.info("Processing chromosome " + chr_id)
                current_chr = chr_id
                current_interval = read_interval
                gene_barcode_dict.clear()

            gene_barcode_dict[gene_id][barcode].append(assignment_info)
            current_interval = (current_interval[0], max(current_interval[1], read_interval[1]))

        processed_read_count, processed_spliced_count = self._process_chunk(gene_barcode_dict, outf, allinfo_outf)
        read_count += processed_read_count
        spliced_count += processed_spliced_count
        outf.close()
        allinfo_outf.close()

        logger.info("Saved %d reads, of them spliced %d to %s" % (read_count, spliced_count, output_prefix))
        logger.info("Total assignments processed %d (typically much more than read count)" % self.total_assignments)
        self.count_stats_for_storage(read_info_storage)
        logger.info("Unique gene-barcodes pairs %d" % len(self.unique_gene_barcode))
        for k in sorted(self.stats.keys()):
            logger.info("%s: %d" % (k, self.stats[k]))

        stats_output = output_prefix + ".stats.tsv"
        logger.info("Stats are written to written to %s" % stats_output)
        with open(stats_output, "w") as count_hist_file:
            count_hist_file.write("Unique gene-barcodes pairs %d\n" % len(self.unique_gene_barcode))

            for k in sorted(self.stats.keys()):
                count_hist_file.write("%s\t%d\n" % (k, self.stats[k]))

    def count_stats(self, assignment_file, output_prefix):
        read_info_storage = defaultdict(list)

        self.unique_gene_barcode = set()
        for l in open(assignment_file):
            if l.startswith("#"): continue
            v = l.strip().split("\t")
            read_id = v[0]
            gene_id = v[4]
            assignment_type = v[5]
            exon_blocks_str = v[7]
            exon_blocks = list(map(lambda x: tuple(map(int, x.split('-'))), exon_blocks_str.split(',')))
            if read_id in self.barcode_dict:
                barcode, umi = self.barcode_dict[read_id]
            else:
                barcode, umi = None, None
            matching_events = v[6]
            read_info_storage[read_id].append(ShortReadAssignmentInfo(gene_id, exon_blocks, assignment_type,
                                                                      matching_events, barcode))

        self.count_stats_for_storage(read_info_storage)
        logger.info("Unique gene-barcodes pairs %d" % len(self.unique_gene_barcode))
        for k in sorted(self.stats.keys()):
            logger.info("%s: %d" % (k, self.stats[k]))

        stats_output = output_prefix + ".stats.tsv"
        logger.info("Stats are written to written to %s" % stats_output)
        with open(stats_output, "w") as count_hist_file:
            count_hist_file.write("Unique gene-barcodes pairs %d\n" % len(self.unique_gene_barcode))

            for k in sorted(self.stats.keys()):
                count_hist_file.write("%s\t%d\n" % (k, self.stats[k]))


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


def create_transcript_info_dict(genedb):
    gffutils_db = gffutils.FeatureDB(genedb)
    transcript_type_dict = {}
    for t in gffutils_db.features_of_type(('transcript', 'mRNA')):
        polya_site = t.start - 1 if t.strand == '-' else t.end + 1
        if "transcript_type" in t.attributes.keys():
            transcript_type_dict[t.id] = (t.attributes["transcript_type"][0], polya_site)
        else:
            transcript_type_dict[t.id] = ("unknown_type", polya_site)
    return transcript_type_dict

