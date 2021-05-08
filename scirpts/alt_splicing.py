#!/usr/bin/env python
#
# ############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import pysam
import gffutils
from src.gtf2db import *
from src.junction_comparator2 import *
from src.common import *
from Bio import SeqIO
from Bio import Seq
from collections import defaultdict
import sys

logger = logging.getLogger('IsoQuant')

canonical_forward_sites = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
canonical_reverse_sites = {("CT", "AC")}

class Params:
    def __init__(self):
        self.delta = 0
        self.minor_exon_extension = 0
        self.major_exon_extension = 0
        self.min_abs_exon_overlap = 10
        self.min_rel_exon_overlap = 0.2
        self.max_suspicious_intron_abs_len = 0
        self.max_suspicious_intron_rel_len = 0
        self.apa_delta = 0
        self.minimal_exon_overlap = 5
        self.minimal_intron_absence_overlap = 20
        self.max_intron_shift = 0
        self.max_missed_exon_len = 0
        self.max_intron_abs_diff = 0
        self.max_intron_rel_diff = 0
        self.max_fake_terminal_exon_len = 0
        self.micro_intron_length = 0
        self.allow_extra_terminal_introns = False
        self.correct_minor_errors = False
        self.has_polya = False


inconsistency_events = {
    MatchEventSubtype.alt_left_site_novel, MatchEventSubtype.alt_right_site_novel,
    MatchEventSubtype.extra_intron_novel, MatchEventSubtype.mutually_exclusive_exons_novel,
    MatchEventSubtype.exon_gain_novel, MatchEventSubtype.exon_skipping_novel,
    MatchEventSubtype.exon_detach_novel, MatchEventSubtype.exon_merge_novel,
    MatchEventSubtype.terminal_exon_shift_novel,
    MatchEventSubtype.alternative_structure_novel, MatchEventSubtype.intron_alternation_novel,
    MatchEventSubtype.unspliced_intron_retention, MatchEventSubtype.intron_retention,
    MatchEventSubtype.alt_left_site_known, MatchEventSubtype.alt_right_site_known,
    MatchEventSubtype.extra_intron_known, MatchEventSubtype.intron_migration,
    MatchEventSubtype.mutually_exclusive_exons_known, MatchEventSubtype.exon_skipping_known,
    MatchEventSubtype.exon_detach_known, MatchEventSubtype.exon_merge_known,
    MatchEventSubtype.terminal_exon_shift_known,
    MatchEventSubtype.exon_gain_known, MatchEventSubtype.alternative_structure_known,
    MatchEventSubtype.intron_alternation_known, MatchEventSubtype.major_exon_elongation_left,
    MatchEventSubtype.incomplete_intron_retention_right, MatchEventSubtype.incomplete_intron_retention_left
}


@unique
class MatchType(Enum):
    match = 0
    consistent_first_longer = 1
    consistent_second_longer = 2
    consistent_differ = 3
    both_unspliced = 10
    only_first_spliced = 11
    only_second_spliced = 12
    inconsistent = 100

    def __lt__(self, other):
        return self.value < other.value

match_subtype_printable_names = \
    {
     MatchEventSubtype.alt_left_site_known: ('alt_donor_site_known', 'alt_acceptor_site_known'),
     MatchEventSubtype.alt_right_site_known: ('alt_acceptor_site_known', 'alt_donor_site_known'),
     MatchEventSubtype.alt_left_site_novel: ('alt_donor_site_novel', 'alt_acceptor_site_novel'),
     MatchEventSubtype.alt_right_site_novel: ('alt_acceptor_site_novel', 'alt_donor_site_novel')}

canonical_splice_sites = {("GT", "AG")}

def events_stat(tsv_file, bam_records, gene_db, params, wrong_reads, ref_seqs, read_chroms):
    stat_map = defaultdict(int)
    counter = 0
    alt_donors = [0] * 21
    alt_acceptors = [0] * 21
    alt_donors_cnt=0
    alt_acceptors_cnt=0
    isoform_strands = dict()
    acc3_file = open("alt_acc_3bp.txt","a")
    donor4_file = open("alt_donor_4bp.txt","a")
    acc3up_file = open("alt_acc_-3bp.txt","a")
    donor4up_file = open("alt_donor_-4bp.txt","a")
    for t in gene_db.features_of_type('transcript', order_by=('seqid', 'start')):
        isoform_strands[t.id] = t.strand
    read_events1 = defaultdict(list)
    read_events2 = defaultdict(list)
    all_events = defaultdict(int)
    with open(tsv_file) as f:
        for line in f:
            readname, isoform_id = line.split()[0], line.split()[1]
            if readname not in bam_records:
                continue
            if readname in wrong_reads: continue
            bam_record = bam_records[readname]
            exons = correct_bam_coords(concat_gapless_blocks(sorted(bam_record.get_blocks()), bam_record.cigartuples))
            mapped_region = (exons[0][0], exons[-1][1])

            gene_list = list(gene_db.region(region=(bam_record.reference_name, mapped_region[0], mapped_region[1]),
                                            completely_within=False))
            if not gene_list:
                continue
            gene_info = GeneInfo(gene_list, gene_db)
            logger.debug("Gene introns: " + str(gene_info.intron_profiles.features))
            introns = junctions_from_blocks(exons)
            if isoform_id not in gene_info.all_isoforms_introns: continue
            isoform_introns = gene_info.all_isoforms_introns[isoform_id]
            isoform_region = (gene_info.transcript_start(isoform_id), gene_info.transcript_end(isoform_id))

            intron_comparator = JunctionComparator(params,
                                                   OverlappingFeaturesProfileConstructor
                                                   (gene_info.intron_profiles.features,
                                                    (gene_info.start, gene_info.end),
                                                    comparator=partial(equal_ranges, delta=params.delta)))
            matching_events = intron_comparator.compare_junctions(introns, mapped_region, isoform_introns, isoform_region)
            for events in matching_events:
                if not isinstance(events, list) and not isinstance(events, tuple): continue
                for event in events:
                    if event[0] in match_subtype_printable_names or event[0] in [MatchEventSubtype.exon_skipping_novel,MatchEventSubtype.exon_skipping_known, MatchEventSubtype.exon_misalignment]:
                        all_events[readname]+=1
                    event_len = event[1]*(-1)
                    event_pos1, event_pos2 = int(event[2]), int(event[3])
                    if event[0] in match_subtype_printable_names:
                        if isoform_strands[isoform_id] == '-':
                            event_len *= -1
                            event_type= match_subtype_printable_names[event[0]][1]
                        else:
                            event_type= match_subtype_printable_names[event[0]][0]
                        donor_site = ref_seqs[read_chroms[readname]][event_pos1-1:event_pos1+1] 
                        acceptor_site = ref_seqs[read_chroms[readname]][event_pos2 - 2:event_pos2]
                        donor_site_ext = ref_seqs[read_chroms[readname]][event_pos1-7:event_pos1+7] 
                        acceptor_site_ext = ref_seqs[read_chroms[readname]][event_pos2 - 8:event_pos2+6]
                        if isoform_strands[isoform_id] == '+':
                            splice_site = (str(donor_site), str(acceptor_site))
                        else:
                            splice_site = (str(Seq.Seq(acceptor_site).reverse_complement()), str(Seq.Seq(donor_site).reverse_complement()))
                            donor_site_ext = str(Seq.Seq(donor_site_ext).reverse_complement())
                            acceptor_site_ext = str(Seq.Seq(acceptor_site_ext).reverse_complement())
                            donor_site_ext, acceptor_site_ext = str(Seq.Seq(acceptor_site_ext).reverse_complement()), str(Seq.Seq(donor_site_ext).reverse_complement())
                        if splice_site not in canonical_splice_sites: 
                            continue
                        if "acc" in event_type: read_events1[readname].append((event_pos1, event_pos2,"acc"+str(event_len)))
                        else: read_events2[readname].append((event_pos1, event_pos2,"donor"+str(event_len)))
                        if "donor" in event_type:
                            if isoform_strands[isoform_id] == '+':
                                if event_len == 4:
                                    donor4_file.write(donor_site_ext + "\t" + read_chroms[readname] + "\t" + readname +  "\t" + str(event_pos1) + "\n")
                                if event_len == -4:
                                    donor4up_file.write(donor_site_ext + "\t" + read_chroms[readname] + "\t" + readname +  "\t" + str(event_pos1) + "\n")
                            alt_donors_cnt += 1
                            if abs(event_len) <=10:
                                alt_donors[10+event_len]+=1
                        else:
                            if isoform_strands[isoform_id] == '+':
                                if event_len == 3:
                                    acc3_file.write(acceptor_site_ext + "\n")
                                if event_len == -3:
                                    acc3up_file.write(acceptor_site_ext + "\n")
                            alt_acceptors_cnt+=1
                            if abs(event_len) <=10:
                                alt_acceptors[10+event_len]+=1
            counter += 1
            if counter % 5000 == 0:
                logger.info("Processed %d reads" % (counter))

    print("Alt donor shifts:",alt_donors)
    print("Alt acc shifts:",alt_acceptors)
    return read_events1, read_events2, all_events


def load_bam(bamfile, read_set=None):
    bam_records = {}
    read_chroms = {}
    for r in pysam.AlignmentFile(bamfile, "rb").fetch():
        if not read_set or r.query_name in read_set:
            if r.is_supplementary or r.is_secondary: continue
            bam_records[r.query_name] = r
            read_chroms[r.query_name] = r.reference_name
    return bam_records, read_chroms


def load_tsv(tsv_file):
    read_pairs = []
    for l in open(tsv_file):
        if l.startswith("#"):
            continue
        t = l.strip().split()
        if len(t) < 4:
            continue
        read_pairs.append((t[2], t[3]))
    return read_pairs


def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--output", "-o", help="output folder, will be created automatically", type=str)
    # REFERENCE
    parser.add_argument("--genedb", "-g", help="gene database in gffutils DB format or GTF/GFF format", type=str)
    parser.add_argument('--bam_pb', type=str, help='sorted and indexed BAM file for PacBio')
    parser.add_argument('--bam_ont', type=str, help='sorted and indexed BAM file for ONT')
    parser.add_argument('--tsv', type=str, help='TSV with barcode and read ids')
    parser.add_argument('--reference', type=str, help='Reference file')
    parser.add_argument('--tsv_pb', type=str, help='Read assignments TSV for PacBio')
    parser.add_argument('--tsv_ont', type=str, help='Read assignments TSV for ONT')
    parser.add_argument('--tsv_pb_sq', type=str, help='Read assignments TSV for PacBio (SQANTI)')
    parser.add_argument('--tsv_ont_sq', type=str, help='Read assignments TSV for ONT (SQANTI)')
    parser.add_argument('--delta', type=int, default=0, help='delta')
    args = parser.parse_args(args)
    return args


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def run_pipeline(args):
    logger.info("Loading gene db from " + args.genedb)
    gene_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    transcripts = SeqIO.parse(args.reference, format="fasta")
    ref_seqs = {seq.id: str(seq.seq) for seq in transcripts}

    wrong_reads = set()
    reads_freq = defaultdict(int)
    with open(args.tsv_pb) as f:
        for line in f:
            fs = line.split()
            if fs[0] in reads_freq: wrong_reads.add(fs[0]) 
            if ',' not in fs[4]: wrong_reads.add(fs[0])
            reads_freq[fs[0]]+=1
    with open(args.tsv_ont) as f:
        for line in f:
            fs = line.split()
            if fs[0] in reads_freq: wrong_reads.add(fs[0]) 
            if ',' not in fs[4]: wrong_reads.add(fs[0])
            reads_freq[fs[0]]+=1
    logger.info("Loading read pairs from " + args.tsv)
    read_pairs = load_tsv(args.tsv)
    for p in read_pairs:
        if p[0] in wrong_reads: wrong_reads.add(p[1])
        elif p[1] in wrong_reads: wrong_reads.add(p[0])
    logger.info("Loading alignments from " + args.bam_pb)
    bam_records1, read_chroms1 = load_bam(args.bam_pb, set(map(lambda x: x[0], read_pairs)))
    logger.info("Loading alignments from " + args.bam_ont)
    bam_records2, read_chroms2 = load_bam(args.bam_ont, set(map(lambda x: x[1], read_pairs)))
    params = Params()
    params.delta = args.delta
    accs_ont, donors_ont,events_ont = events_stat(args.tsv_ont, bam_records2, gene_db, params, wrong_reads, ref_seqs,read_chroms2)
    accs_pb,donors_pb,events_pb = events_stat(args.tsv_pb, bam_records1, gene_db, params, wrong_reads, ref_seqs, read_chroms1)

    ont, pb, shared = 0, 0, 0
    for p in read_pairs:
        for e in accs_pb[p[0]]:
            if e in accs_ont[p[1]] and e[-1]=="acc3": shared +=1 
            elif e[-1]=="acc3": pb+=1
        for e in accs_ont[p[1]]:
            if e not in accs_pb[p[0]] and e[-1]=="acc3": ont +=1 
    print("Acceptors 3bp. ONT only:", ont,", PB only:", pb, ", Shared:", shared)
    ont, pb, shared = 0, 0, 0
    for p in read_pairs:
        for e in accs_pb[p[0]]:
            if e[-1]!="acc-3": continue
            if e in accs_ont[p[1]]: shared +=1 
            else: pb+=1
        for e in accs_ont[p[1]]:
            if e not in accs_pb[p[0]] and e[-1]=="acc-3": ont +=1 
    print("Acceptors -3bp. ONT only:", ont,", PB only:", pb, ", Shared:", shared)
    ont, pb, shared = 0, 0, 0
    for p in read_pairs:
        for e in donors_pb[p[0]]:
            if e in donors_ont[p[1]] and e[-1]=="donor4": shared +=1 
            elif e[-1]=="donor4": pb+=1
        for e in donors_ont[p[1]]:
            if e not in donors_pb[p[0]] and e[-1]=="donor4": ont +=1 
    print("Donors 4bp. ONT only:", ont,", PB only:", pb, ", Shared:", shared)
    ont, pb, shared = 0, 0, 0
    for p in read_pairs:
        for e in donors_pb[p[0]]:
            if e[-1]!="donor-4": continue
            if e in donors_ont[p[1]]: shared +=1 
            else: pb+=1
        for e in donors_ont[p[1]]:
            if e not in donors_pb[p[0]] and e[-1]=="donor-4": ont +=1 
    print("Donors -4bp. ONT only:", ont,", PB only:", pb, ", Shared:", shared)

    ont, pb, shared = 0, 0, 0
    for p in read_pairs:
        for e in accs_pb[p[0]]:
            if e in accs_ont[p[1]]: shared +=1
            else: pb+=1
        for e in accs_ont[p[1]]:
            if e not in accs_pb[p[0]]: ont +=1
    print("Acceptors. ONT only:", ont,", PB only:", pb, ", Shared:", shared)
    ont, pb, shared = 0, 0, 0
    for p in read_pairs:
        for e in donors_pb[p[0]]:
            if e in donors_ont[p[1]]: shared +=1
            else: pb+=1
        for e in donors_ont[p[1]]:
            if e not in donors_pb[p[0]]: ont +=1
    print("Donors. ONT only:", ont,", PB only:", pb, ", Shared:", shared)
    logger.info("Done")

def main(args):
    args = parse_args(args)
    set_logger(logger)
    run_pipeline(args)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

