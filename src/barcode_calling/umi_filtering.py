
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details
############################################################################

import pysam
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional
import logging
import editdistance
import gffutils

from src.assignment_loader import create_merging_assignment_loader
from src.isoform_assignment import MatchEventSubtype, ReadAssignment

logger = logging.getLogger('IsoQuant')


def format_read_assignment_for_output(read_assignment: ReadAssignment) -> str:
    """
    Format ReadAssignment for UMI filtering output file.

    Expects additional_attributes to contain:
    - 'transcript_type': Type of transcript (e.g., 'protein_coding')
    - 'polya_site': PolyA site position (-1 if none)
    - 'cell_type': Cell type/spot from barcode mapping

    Args:
        read_assignment: ReadAssignment object with additional_attributes set

    Returns:
        Formatted string for output file
    """
    exon_blocks = read_assignment.corrected_exons
    chr_id = read_assignment.chr_id
    strand = read_assignment.strand

    intron_blocks = read_assignment.corrected_introns
    exons_str = ";%;" + ";%;".join(["%s_%d_%d_%s" % (chr_id, e[0], e[1], strand) for e in exon_blocks])
    introns_str = ";%;" + ";%;".join(["%s_%d_%d_%s" % (chr_id, e[0], e[1], strand) for e in intron_blocks])

    # Determine read type from assignment type
    if read_assignment.assignment_type.is_unique():
        read_type = "known"
    elif read_assignment.assignment_type.is_inconsistent():
        read_type = "novel"
    elif read_assignment.assignment_type.is_ambiguous():
        read_type = "known_ambiguous"
    else:
        read_type = "none"

    # Get gene and transcript IDs
    if read_assignment.isoform_matches:
        gene_id = read_assignment.isoform_matches[0].assigned_gene
        transcript_id = read_assignment.isoform_matches[0].assigned_transcript
        matching_events = [e.event_type for e in read_assignment.isoform_matches[0].match_subclassifications]
    else:
        gene_id = "."
        transcript_id = "."
        matching_events = []

    # Get additional attributes
    transcript_type = read_assignment.additional_attributes.get('transcript_type', 'unknown')
    polya_site = read_assignment.additional_attributes.get('polya_site', -1)
    cell_type = read_assignment.additional_attributes.get('cell_type', 'None')

    # Format TSS and polyA
    polyA = "NoPolyA"
    TSS = "NoTSS"

    if isinstance(matching_events, str):
        # Legacy string format
        if "tss_match" in matching_events:
            tss_pos = exon_blocks[-1][1] if strand == "-" else exon_blocks[0][0]
            TSS = "%s_%d_%d_%s" % (chr_id, tss_pos, tss_pos, strand)
        if "correct_polya" in matching_events and polya_site != -1:
            polyA_pos = polya_site
            polyA = "%s_%d_%d_%s" % (chr_id, polyA_pos, polyA_pos, strand)
    else:
        # Event subtype format
        if strand == '+':
            if any(x in [MatchEventSubtype.terminal_site_match_left,
                         MatchEventSubtype.terminal_site_match_left_precise] for x in matching_events):
                tss_pos = exon_blocks[0][0]
                TSS = "%s_%d_%d_%s" % (chr_id, tss_pos, tss_pos, strand)
            if (polya_site != -1 and
                    any(x == MatchEventSubtype.correct_polya_site_right for x in matching_events)):
                polyA_pos = polya_site
                polyA = "%s_%d_%d_%s" % (chr_id, polyA_pos, polyA_pos, strand)
        elif strand == '-':
            if any(x in [MatchEventSubtype.terminal_site_match_right,
                         MatchEventSubtype.terminal_site_match_right_precise] for x in matching_events):
                tss_pos = exon_blocks[-1][1]
                TSS = "%s_%d_%d_%s" % (chr_id, tss_pos, tss_pos, strand)
            if (polya_site != -1 and
                    any(x == MatchEventSubtype.correct_polya_site_left for x in matching_events)):
                polyA_pos = polya_site
                polyA = "%s_%d_%d_%s" % (chr_id, polyA_pos, polyA_pos, strand)

    return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s" % (
        read_assignment.read_id, gene_id, cell_type, read_assignment.barcode, read_assignment.umi,
        introns_str, TSS, polyA, exons_str, read_type, len(intron_blocks), transcript_id, transcript_type)


class UMIFilter:
    """
    UMI-based deduplication filter for barcoded reads.

    Processes reads with UMI (Unique Molecular Identifier) tags to remove PCR/RT duplicates.
    Groups reads by gene, barcode, and UMI, selecting representative reads from each molecule.
    """

    def __init__(self, umi_length: int = 0, edit_distance: int = 3, disregard_length_diff: bool = True,
                 only_unique_assignments: bool = False, only_spliced_reads: bool = False):
        """
        Initialize UMI filter.

        Args:
            umi_length: Expected UMI length (0 = variable length)
            edit_distance: Maximum edit distance for UMI clustering
            disregard_length_diff: Whether to ignore length differences in edit distance calculation
            only_unique_assignments: Only process uniquely assigned reads
            only_spliced_reads: Only process spliced reads
        """
        self.umi_length = umi_length
        self.max_edit_distance = edit_distance
        self.disregard_length_diff = disregard_length_diff
        self.only_unique_assignments = only_unique_assignments
        self.only_spliced_reads = only_spliced_reads

        self.selected_reads: Set[str] = set()
        self.stats: Dict[str, int] = defaultdict(int)
        self.unique_gene_barcode: Set[Tuple[str, str]] = set()
        self.total_assignments = 0
        self.duplicated_molecule_counts: Dict[int, int] = defaultdict(int)

        self.umi_len_dif_func = (lambda x: 0) if self.umi_length == 0 else (lambda x: -abs(self.umi_length - x))

    def _find_similar_umi(self, umi: str, trusted_umi_list: List[str]) -> Optional[str]:
        """
        Find similar UMI in trusted list within edit distance threshold.

        Args:
            umi: Query UMI sequence
            trusted_umi_list: List of already accepted UMIs

        Returns:
            Most similar UMI if found within threshold, None otherwise
        """
        if self.max_edit_distance == -1:
            return None if not trusted_umi_list else trusted_umi_list[0]

        similar_umi = None
        best_dist = 100

        for occ in trusted_umi_list:
            if self.max_edit_distance == 0:
                # Exact match mode
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

    def _construct_umi_dict(self, molecule_list: List[ReadAssignment]) -> Dict[str, List[ReadAssignment]]:
        """
        Group molecules by UMI, clustering similar UMIs.

        Uses greedy clustering: processes UMIs by frequency, merging similar low-frequency
        UMIs into high-frequency ones.

        Args:
            molecule_list: List of read assignments with UMI information

        Returns:
            Dict mapping representative UMI to list of reads
        """
        # Count reads per UMI
        umi_counter: Dict[str, Set[str]] = defaultdict(set)
        for m in molecule_list:
            umi_counter[m.umi].add(m.read_id)

        # Create sorting keys: (count, length_penalty, umi_sequence)
        umi_sorting_keys = {}
        for umi in umi_counter:
            if umi is None or umi == "None":
                umi_sorting_keys[umi] = (-1, 0, "")
            else:
                umi_sorting_keys[umi] = (len(umi_counter[umi]), self.umi_len_dif_func(len(umi)), umi)

        # Sort molecules by UMI frequency (process high-frequency UMIs first)
        # This ensures reads with same UMI are processed consecutively
        molecule_list = sorted(molecule_list, key=lambda ml: umi_sorting_keys[ml.umi], reverse=True)

        umi_dict: Dict[str, List[ReadAssignment]] = defaultdict(list)
        trusted_umi_list: List[str] = []

        for m in molecule_list:
            if m.umi is None or m.umi == "None":
                # Collect untrusted UMIs together
                umi_dict["None"].append(m)
                continue

            similar_umi = self._find_similar_umi(m.umi, trusted_umi_list)
            if similar_umi is None:
                # New distinct UMI
                umi_dict[m.umi].append(m)
                trusted_umi_list.append(m.umi)
            else:
                # Merge with existing UMI
                umi_dict[similar_umi].append(m)

        return umi_dict

    def _process_duplicates(self, molecule_list: List[ReadAssignment]) -> List[ReadAssignment]:
        """
        Process PCR/RT duplicates for a gene-barcode pair.

        Groups reads by UMI and selects representative read for each molecule.
        Selection criteria (in order):
        1. Unique assignment > ambiguous
        2. More exons
        3. Longer transcript

        Args:
            molecule_list: All reads for this gene-barcode pair

        Returns:
            List of selected representative reads
        """
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

            # Select best read from duplicates
            best_read = umi_dict[umi][0]
            logger.debug("Selecting from:")
            for m in umi_dict[umi]:
                logger.debug("%s %s" % (m.read_id, m.umi))
                # Prefer unique assignments
                if not best_read.assignment_type.is_unique() and m.assignment_type.is_unique():
                    best_read = m
                # Prefer more exons
                elif len(m.corrected_exons) > len(best_read.corrected_exons):
                    best_read = m
                # Prefer longer transcripts
                elif (len(m.corrected_exons) == len(best_read.corrected_exons) and
                      m.corrected_exons[-1][1] - m.corrected_exons[0][0] >
                      best_read.corrected_exons[-1][1] - best_read.corrected_exons[0][0]):
                    best_read = m

            # Check for ambiguity in transcript annotation among duplicates with same read ID
            polyas = set()
            transcript_types = set()
            isoform_ids = set()
            for m in umi_dict[umi]:
                if m.read_id != best_read.read_id:
                    continue
                transcript_types.add(m.additional_attributes.get('transcript_type', 'unknown'))
                polyas.add(m.additional_attributes.get('polya_site', -1))
                if m.isoform_matches:
                    isoform_ids.add(m.isoform_matches[0].assigned_transcript)

            # Resolve ambiguities by updating additional_attributes
            if len(transcript_types) > 1:
                best_read.set_additional_attribute('transcript_type', 'None')
            if len(polyas) > 1:
                best_read.set_additional_attribute('polya_site', -1)
            if len(isoform_ids) > 1:
                # Clear the transcript ID by removing all matches except the gene
                if best_read.isoform_matches:
                    best_read.isoform_matches[0].assigned_transcript = "None"

            # Clear annotation for inconsistent assignments
            if best_read.assignment_type.is_inconsistent():
                if best_read.isoform_matches:
                    best_read.isoform_matches[0].assigned_transcript = "None"
                best_read.set_additional_attribute('polya_site', -1)
                best_read.set_additional_attribute('transcript_type', 'None')

            logger.debug("Selected %s %s" % (best_read.read_id, best_read.umi))
            resulting_reads.append((best_read, umi))

        # If single UMI, trust it regardless of whether it's marked as trusted
        if len(resulting_reads) == 1:
            return [resulting_reads[0][0]]

        # With multiple UMIs, ignore untrusted ones
        return [x[0] for x in filter(lambda x: x[1] != "None", resulting_reads)]

    def _process_gene(self, gene_dict: Dict[str, List[ReadAssignment]]):
        """
        Process all barcodes for a gene.

        Args:
            gene_dict: Dict mapping barcode to list of reads

        Yields:
            Representative read for each molecule
        """
        for barcode in gene_dict:
            for r in self._process_duplicates(gene_dict[barcode]):
                yield r

    def _process_chunk(self, gene_barcode_dict: Dict[str, Dict[str, List[ReadAssignment]]],
                      allinfo_outf, read_ids_outf=None) -> Tuple[int, int]:
        """
        Process a chunk of reads grouped by gene and barcode.

        Args:
            gene_barcode_dict: Nested dict[gene_id][barcode] -> list of reads
            allinfo_outf: Output file for detailed assignment info
            read_ids_outf: Optional output file for selected read IDs

        Returns:
            Tuple of (total_read_count, spliced_read_count)
        """
        read_count = 0
        spliced_count = 0

        for gene_id in gene_barcode_dict:
            for read_assignment in self._process_gene(gene_barcode_dict[gene_id]):
                # Skip non-unique reads if already selected
                if (not read_assignment.assignment_type.is_unique() and
                        read_assignment.read_id in self.selected_reads):
                    continue

                read_count += 1
                if read_assignment.isoform_matches:
                    gene_id = read_assignment.isoform_matches[0].assigned_gene
                    self.unique_gene_barcode.add((gene_id, read_assignment.barcode))

                if len(read_assignment.corrected_exons) > 1:
                    spliced_count += 1

                if read_ids_outf:
                    read_ids_outf.write(read_assignment.read_id + "\n")

                self.selected_reads.add(read_assignment.read_id)

                if allinfo_outf:
                    allinfo_outf.write(format_read_assignment_for_output(read_assignment) + "\n")

        return read_count, spliced_count

    def add_stats_for_read(self, read_assignment: ReadAssignment):
        """
        Update statistics for a processed read.

        Args:
            read_assignment: ReadAssignment object
        """
        gene_id = "."
        if read_assignment.isoform_matches:
            gene_id = read_assignment.isoform_matches[0].assigned_gene

        assigned = gene_id != "."
        spliced = len(read_assignment.corrected_exons) > 1
        barcoded = read_assignment.barcode is not None
        unique = assigned  # In simplified logic, we only process unique/consistent assignments

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

    def process_single_chr(self, chr_id: str, saves_prefix: str, transcript_type_dict: Dict[str, Tuple[str, int]],
                          barcode_feature_table: Dict[str, str],
                          all_info_file_name: str, filtered_reads_file_name: Optional[str],
                          stats_output_file_name: str, string_pools=None) -> Tuple[str, str]:
        """
        Process UMI filtering for a single chromosome.

        Args:
            chr_id: Chromosome ID
            saves_prefix: Prefix for saved assignment files
            transcript_type_dict: Dict mapping transcript_id to (type, polya_site)
            barcode_feature_table: Dict mapping barcode to cell type
            all_info_file_name: Output file for detailed assignment info
            filtered_reads_file_name: Optional output file for filtered read IDs
            stats_output_file_name: Output file for statistics

        Returns:
            Tuple of (all_info_file_name, stats_output_file_name)
        """
        with open(all_info_file_name, "w") as allinfo_outf:
            filtered_reads_outf = open(filtered_reads_file_name, "w") if filtered_reads_file_name else None
            read_count = 0
            spliced_count = 0
            self.unique_gene_barcode = set()

            loader = create_merging_assignment_loader(chr_id, saves_prefix, string_pools=string_pools)
            while loader.has_next():
                gene_barcode_dict: Dict[str, Dict[str, List[ReadAssignment]]] = defaultdict(lambda: defaultdict(list))
                _, assignment_storage = loader.get_next()
                logger.debug("Processing %d reads" % len(assignment_storage))

                for read_assignment in assignment_storage:
                    # Skip unassigned reads
                    if read_assignment.gene_assignment_type.is_unassigned():
                        continue

                    # Count ambiguous assignments but don't process them
                    if read_assignment.gene_assignment_type.is_ambiguous():
                        self.stats["Assigned to any gene"] += 1
                        if len(read_assignment.corrected_exons) > 1:
                            self.stats["Spliced"] += 1
                        continue

                    exon_blocks = read_assignment.corrected_exons
                    barcode = read_assignment.barcode
                    if barcode == '*' or not barcode:
                        barcode = None
                    spliced = len(exon_blocks) > 1
                    barcoded = barcode is not None

                    # Get cell types from barcode feature table (may be a list for multi-column)
                    if barcode is None or barcode not in barcode_feature_table:
                        cell_type = "None"
                    else:
                        cell_types = barcode_feature_table[barcode]
                        # Handle both old format (string) and new format (list)
                        if isinstance(cell_types, list):
                            cell_type = ",".join(cell_types) if cell_types else "None"
                        else:
                            cell_type = cell_types
                    read_assignment.set_additional_attribute('cell_type', cell_type)

                    # Process based on number of isoform matches
                    if len(read_assignment.isoform_matches) == 1:
                        # Single match - simplest case
                        isoform_match = read_assignment.isoform_matches[0]
                        if read_assignment.assignment_type.is_consistent():
                            transcript_id = isoform_match.assigned_transcript
                            transcript_type, polya_site = (
                                transcript_type_dict[transcript_id] if transcript_id in transcript_type_dict
                                else ("unknown_type", -1))
                            read_assignment.set_additional_attribute('transcript_type', transcript_type)
                            read_assignment.set_additional_attribute('polya_site', polya_site)
                        else:
                            # Inconsistent single match
                            read_assignment.set_additional_attribute('transcript_type', 'unknown_type')
                            read_assignment.set_additional_attribute('polya_site', -1)

                    elif read_assignment.assignment_type.is_consistent():
                        # Multiple consistent matches - need to resolve ambiguity
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

                        # Use consensus if unique, otherwise mark as ambiguous
                        transcript_type = "None" if len(transcript_types) != 1 else transcript_types.pop()
                        polya_site = -1 if len(polya_sites) != 1 else polya_sites.pop()

                        read_assignment.set_additional_attribute('transcript_type', transcript_type)
                        read_assignment.set_additional_attribute('polya_site', polya_site)

                        # Set transcript_id to "None" for ambiguous case
                        if read_assignment.isoform_matches:
                            read_assignment.isoform_matches[0].assigned_transcript = "None"
                    else:
                        # Multiple inconsistent matches - skip
                        self.total_assignments += 1
                        continue

                    self.total_assignments += 1
                    self.add_stats_for_read(read_assignment)

                    # Filter for deduplication
                    if not barcoded:
                        continue
                    if not spliced and self.only_spliced_reads:
                        continue

                    # Add to gene-barcode dict using gene from first isoform match
                    if read_assignment.isoform_matches:
                        gene_id = read_assignment.isoform_matches[0].assigned_gene
                        gene_barcode_dict[gene_id][barcode].append(read_assignment)

                # Process chunk
                processed_read_count, processed_spliced_count = self._process_chunk(
                    gene_barcode_dict, allinfo_outf, filtered_reads_outf)
                read_count += processed_read_count
                spliced_count += processed_spliced_count

            if filtered_reads_outf:
                filtered_reads_outf.close()

        # Write statistics
        with open(stats_output_file_name, "w") as count_hist_file:
            count_hist_file.write("Unique gene-barcodes pairs\t%d\n" % len(self.unique_gene_barcode))
            count_hist_file.write("Total reads saved\t%d\n" % read_count)
            count_hist_file.write("Spliced reads saved\t%d\n" % spliced_count)
            count_hist_file.write("Total assignments processed\t%d\n" % self.total_assignments)
            for k in sorted(self.stats.keys()):
                count_hist_file.write("%s\t%d\n" % (k, self.stats[k]))

        return all_info_file_name, stats_output_file_name


def filter_bam(in_file_name: str, out_file_name: str, read_set: Set[str]):
    """
    Filter BAM file to keep only reads in the provided set.

    Args:
        in_file_name: Input BAM file path
        out_file_name: Output BAM file path
        read_set: Set of read IDs to keep
    """
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


def create_transcript_info_dict(genedb: str, chr_ids: Optional[List[str]] = None) -> Dict[str, Tuple[str, int]]:
    """
    Create dictionary of transcript types and polyA sites from gene database.

    Args:
        genedb: Path to gffutils database
        chr_ids: Optional list of chromosome IDs to filter

    Returns:
        Dict mapping transcript_id to (transcript_type, polya_site) tuple
    """
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
