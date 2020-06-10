############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from Bio import SeqIO

from src.common import *
from src.long_read_assigner import *


logger = logging.getLogger('IsoQuant')

class PrintAllFunctor:
    def check(sefl, assignment):
        return True


class PrintOnlyFunctor:
    def __init__(self, allowed_types):
        if isinstance(allowed_types, list):
            self.allowed_types = set(allowed_types)
        elif isinstance(allowed_types, set):
            self.allowed_types = allowed_types
        else:
            self.allowed_types = set([allowed_types])

    def check(self, assignment):
        return assignment.assignment_type in self.allowed_types


class AbstractAssignmentPrinter:
    def __init__(self, output_file_name, params, assignment_checker=PrintAllFunctor()):
        self.assignment_checker = assignment_checker
        self.params = params
        self.output_file = open(output_file_name, "w")

    def __del__(self):
        self.output_file.close()

    def add_read_info(self, read_assignment):
        raise NotImplementedError()

    def flush(self):
        self.output_file.flush()


class ReadAssignmentCompositePrinter:
    def __init__(self, printers):
        self.printers = printers

    def add_read_info(self, read_assignment):
        for p in self.printers:
            p.add_read_info(read_assignment)

    def flush(self):
        for p in self.printers:
            p.flush()


# TODO: reformat output, make singe file
class BasicTSVAssignmentPrinter(AbstractAssignmentPrinter):
    def __init__(self, output_file_name, params, assignment_checker=PrintAllFunctor()):
        AbstractAssignmentPrinter.__init__(self, output_file_name, params, assignment_checker)
        self.header = "#read_id\tisoform_id\tassignment_type\tassignment_events"
        if self.params.print_additional_info:
            self.header += "\taligned_blocks\tintron_profile\tsplit_exon_profile"
        self.header += "\n"
        self.output_file.write(self.header)

    def add_read_info(self, read_assignment):
        if self.assignment_checker is None or not self.assignment_checker.check(read_assignment):
            return
        if read_assignment.assignment_type is None or read_assignment.isoform_matches is None:
            line = read_assignment.read_id  + "\t.\t.\t."
        else:
            assigned_transcripts = [m.assigned_transcript for m in read_assignment.isoform_matches]
            # FIXME empty match subclassification
            match_events = ["+".join([x.name for x in m.match_subclassifications]) for m in read_assignment.isoform_matches]
            line = read_assignment.read_id  + "\t" + ",".join(assigned_transcripts) + "\t" \
                    + read_assignment.assignment_type.name + "\t" + ",".join(match_events)
        if self.params.print_additional_info:
            combined_read_profile = read_assignment.combined_profile
            if combined_read_profile is None:
                line += "\t.\t.\t."
            else:
                line += "\t" + range_list_to_str(combined_read_profile.read_split_exon_profile.read_features) + "\t" + \
                    list_to_str(combined_read_profile.read_intron_profile.gene_profile) + "\t" + \
                        list_to_str(combined_read_profile.read_split_exon_profile.gene_profile)
        line += "\n"
        self.output_file.write(line)

    """ sqanti output
    
    isoform: the isoform ID. Usually in PB.X.Y format.
    chrom: chromosome.
    strand: strand.
    length: isoform length.
    exons: number of exons.
    structural_category: one of the categories ["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "antisense", "fusion", "intergenic", "genic_intron"]
    associated_gene: the reference gene name.
    associated_transcript: the reference transcript name.
    ref_length: reference transcript length.
    ref_exons: reference transcript number of exons.
    diff_to_TSS: distance of query isoform 5' start to reference transcript start end. Negative value means query starts downstream of reference.
    diff_to_TTS: distance of query isoform 3' end to reference annotated end site. Negative value means query ends upstream of reference.
    diff_to_gene_TSS: distance of query isoform 5' start to the closest start end of any transcripts of the matching gene. This field is different from diff_to_TSS since it's looking at all annotated starts of a gene. Negative value means query starts downstream of reference.
    diff_to_gene_TTS: distance of query isoform 3' end to the closest end of any transcripts of the matching gene. Negative value means query ends upstream of reference.
    subcategory: additional splicing categorization, separated by semi-colons. Categories include: mono-exon, multi-exon. Intron rentention is marked with intron_retention.
    RTS_stage: TRUE if one of the junctions could be a RT switching artifact.
    all_canonical: TRUE if all junctions have canonical splice sites.
    - min_sample_cov: sample with minimum coverage.
    - min_cov: minimum junction coverage based on short read STAR junction output file. NA if no short read given.
    - min_cov_pos: the junction that had the fewest coverage. NA if no short read data given.
    - sd_cov: standard deviation of junction coverage counts from short read data. NA if no short read data given.
    - FL or FL.<sample>: FL count associated with this isoform per sample if --fl_count is provided, otherwise NA.
    n_indels: total number of indels based on alignment.
    n_indels_junc: number of junctions in this isoform that have alignment indels near the junction site (indicating potentially unreliable junctions).
    bite: TRUE if all junctions match reference junctions completely.
    - iso_exp: short read expression for this isoform if --expression is provided, otherwise NA.
    - gene_exp: short read expression for the gene associated with this isoform (summing over all isoforms) if --expression is provided, otherwise NA.
    - ratio_exp: ratio of iso_exp to gene_exp if --expression is provided, otherwise NA.
    - FSM_class: ignore this field for now.
    - ORF_length: predicted ORF length.
    - CDS_length: predicted CDS length.
    - CDS_start: CDS start.
    - CDS_end: CDS end.
    CDS_genomic_start: genomic coordinate of the CDS start. If on - strand, this coord will be greater than the end.
    CDS_genomic_end: genomic coordinate of the CDS end. If on - strand, this coord will be smaller than the start.
    - predicted_NMD: TRUE if there's a predicted ORF and CDS ends before the last junction; FALSE if otherwise. NA if non-coding.
    perc_A_downstreamTTS: percent of genomic "A"s in the downstream 20 bp window. If this number if high (say > 0.8), the 3' end site of this isoform is probably not reliable.
    seq_A_downstream_TTS: sequence of the downstream 20 bp window.
    dist_to_cage_peak: distance to closest TSS based on CAGE Peak data. Negative means upstream of TSS and positive means downstream of TSS. Strand-specific. SQANTI2 only searches for nearby CAGE Peaks within 10000 bp of the PacBio transcript start site. Will be NA if none are found within 10000 bp.
    within_cage_peak: TRUE if the PacBio transcript start site is within a CAGE Peak.
    dist_to_polya_site: distance to closest polyA site based on PolyA site data. Negative means upstream of closest site and positive means downstream of closest site. Strand-specific. SQANTI2 only searches for nearby CAGE Peaks within 100 bp of the PacBio transcript start site. Will be NA if none are found within 100 bp.
    within_polya_site: TRUE if the PacBio transcript polyA site is within 100 bp (upstream or downstream) of an annotated polyA site.
    polyA_motif: if --polyA_motif_list is given, shows the top ranking polyA motif found within 50 bp upstream of end.
    polyA_dist: if --polyA_motif_list is given, shows the location of the last base of the hexamer. Position 0 is the putative poly(A) site. This distance is hence always negative because it is upstream.
    
    """

class SqantiTSVPrinter(AbstractAssignmentPrinter):
    canonical_donor_sites = ["GT", "GC", "AT"]
    canonical_acceptor_sites = ["AG", "AC"]
    canonical_donor_sites_rc = ["AC", "GC", "AT"]
    canonical_acceptor_sites_rc = ["CT", "GT"]
    upstream_region_len = 20

    def __init__(self, output_file_name, params):
        AbstractAssignmentPrinter.__init__(self, output_file_name, params)
        self.header = 'isoform\tchrom\tstrand\tlength\texons\tstructural_category' \
                      '\tassociated_gene\tassociated_transcript\tref_length\tref_exons\tdiff_to_TSS\tdiff_to_TTS' \
                      '\tdiff_to_gene_TSS\tdiff_to_gene_TTS\tsubcategory\tall_canonical' \
                      '\tn_indels\tn_indels_junc\tbite\tCDS_genomic_start' \
                      '\tCDS_genomic_end\tperc_A_downstreamTTS\tseq_A_downstream_TTS\tdist_to_cage_peak' \
                      '\twithin_cage_peak\tdist_to_polya_site\twithin_polya_site\tpolyA_motif\tpolyA_dist\n'
        self.output_file.write(self.header)

    def add_read_info(self, read_assignment):
        # FIXME ambiguous matches
        if read_assignment.assignment_type in [ReadAssignmentType.empty, ReadAssignmentType.ambiguous]:
            return

        gene_info = read_assignment.gene_info
        match = read_assignment.isoform_matches[0]
        gene_id = match.assigned_gene
        transcript_id = match.assigned_transcript
        strand = gene_info.isoform_strands[transcript_id]

        # FIXME not genomic distance
        dist_to_tss = read_assignment.start() - gene_info.transcript_start(transcript_id)
        dist_to_tts = gene_info.transcript_end(transcript_id) - read_assignment.end()
        dist_to_gene_tss, dist_to_gene_tts = self.find_closests_tsts(read_assignment)
        if strand == '-':
            # swapping distances since 5' and 3' are now opposite
            dist_to_tss, dist_to_tts = dist_to_tts, dist_to_tss
            dist_to_gene_tss, dist_to_gene_tts = dist_to_gene_tts, dist_to_gene_tss

        subtypes = ";".join([x.name for x in match.match_subclassifications])
        if not subtypes:
            subtypes = "."

        indel_count = read_assignment.additional_info["indel_count"] \
            if "indel_count" in read_assignment.additional_info.keys() else "-1"
        junctions_with_indels = read_assignment.additional_info["junctions_with_indels"] \
            if "junctions_with_indels" in read_assignment.additional_info.keys() else "-1"
        bite = self.check_all_sites_match_reference(read_assignment)
        ref_cds_start, ref_cds_end = self.find_ref_CDS_region(gene_info, transcript_id)

        if self.params.reference:
            record_dict = SeqIO.to_dict(SeqIO.parse(self.params.reference, "fasta"))
            chr_seq = record_dict[gene_info.chr_id]
            read_introns = read_assignment.combined_read_profile.read_intron_profile.read_features
            all_canonical = str(self.check_sites_are_canonical(read_introns, chr_seq, strand))
            seq_A_downstream_TTS, perc_A_downstreamTTS = \
                self.check_downstream_polya((read_assignment.start(), read_assignment.end()), chr_seq, strand)
            perc_A_downstreamTTS = "%0.3f" % perc_A_downstreamTTS
        else:
            all_canonical = "NA"
            perc_A_downstreamTTS = "NA"
            seq_A_downstream_TTS = "NA"

        # TODO
        dist_to_cage_peak = "NA"
        within_cage_peak = "NA"
        dist_to_polya_site = "NA"
        within_polya_site = "NA"
        polyA_motif = "NA"
        polyA_dist = "NA"

        l = "\t".join([str(x) for x in [read_assignment.read_id, gene_info.chr_id, strand,
                                        read_assignment.length(), read_assignment.exon_count(),
                                        match.match_classification.name, gene_id, transcript_id,
                                        gene_info.total_transcript_length(transcript_id),
                                        gene_info.transcript_exon_count(transcript_id),
                                        dist_to_tss, dist_to_tts, dist_to_gene_tss, dist_to_gene_tts, subtypes,
                                        all_canonical, indel_count, junctions_with_indels, bite,
                                        ref_cds_start, ref_cds_end, perc_A_downstreamTTS, seq_A_downstream_TTS,
                                        dist_to_cage_peak, within_cage_peak, dist_to_polya_site,
                                        within_polya_site, polyA_motif, polyA_dist]])
        self.output_file.write(l + "\n")

    def find_closests_tsts(self, read_assignment):
        # FIXME not genomic distance
        gene_info = read_assignment.gene_info
        match = read_assignment.isoform_matches[0]
        gene_id = match.assigned_gene

        dists_to_gene_tss = []
        dists_to_gene_tts = []
        for t in gene_info.db.children(gene_info.db[gene_id], featuretype='transcript', order_by='start'):
            dists_to_gene_tss.append(read_assignment.start() - t.start)
            dists_to_gene_tts.append(t.end - read_assignment.end())
        dist_to_gene_tss = dists_to_gene_tss[argmin([abs(x) for x in dists_to_gene_tss])]
        dist_to_gene_tts = dists_to_gene_tts[argmin([abs(x) for x in dists_to_gene_tts])]
        return dist_to_gene_tss, dist_to_gene_tts

    def find_ref_CDS_region(self, gene_info, transcript_id):
        cds = [c for c in gene_info.db.children(gene_info.db[transcript_id], featuretype='CDS', order_by='start')]
        if len(cds) == 0:
            return -1, -1
        return cds[0].start, cds[-1].end

    def check_all_sites_match_reference(self, read_assignment):
        # FIXME: make a new profile with delta = 0
        if len(read_assignment.combined_profile.read_intron_profile.read_features) == 0:
            return "Unspliced"
        introns_match = all(e == 1 for e in read_assignment.combined_profile.read_intron_profile.read_profile)
        return str(introns_match)

    def check_sites_are_canonical(self, read_introns, chr_seq, strand):
        for intron in read_introns:
            left_site = str(chr_seq[intron[0] - 1:intron[0] + 1].seq)
            right_site = str(chr_seq[intron[1] - 1:intron[1] + 1].seq)
            if strand == '+' and \
                    (left_site not in self.canonical_acceptor_sites or right_site not in self.canonical_donor_sites):
                return False
            elif strand == '-' and \
                    (left_site not in self.canonical_donor_sites_rc or right_site not in self.canonical_acceptor_sites_rc):
                return False
        return True

    def check_downstream_polya(self, read_coords, chr_seq, strand):
        if strand == '+':
            seq = str(chr_seq[read_coords[1]:read_coords[1] + self.upstream_region_len].seq)
            a_percentage = float(seq.upper().count('A')) / float(self.upstream_region_len)
        else:
            seq = str(chr_seq[read_coords[0] - self.upstream_region_len - 1:read_coords[0] - 1].seq)
            a_percentage = float(seq.upper().count('T')) / float(self.upstream_region_len)
        return seq, a_percentage

    def __del__(self):
        self.output_file.close()

    def flush(self):
        self.output_file.flush()