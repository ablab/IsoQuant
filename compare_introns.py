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
from src.junction_comparator import *
from src.common import *


logger = logging.getLogger('IsoQuant')

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


inconsistency_events = nnic_event_types | nic_event_types


@unique
class MatchType(Enum):
    match = 0
    inconsistent = 10
    first_only = 20
    second_only = 30
    consistent_first_longer = 1
    consistent_second_longer = 2
    consistent_differ = 3


def compare_intron_chains(exons1, exons2, gene_info, params):
    mapped_region1 = (exons1[0][0], exons1[-1][1])
    introns1 = junctions_from_blocks(exons1)
    mapped_region2 = (exons2[0][0], exons2[-1][1])
    introns2 = junctions_from_blocks(exons2)

    intron_comparator = JunctionComparator(params,
                                           OverlappingFeaturesProfileConstructor
                                           (gene_info.intron_profiles.features,
                                            (gene_info.start, gene_info.end),
                                            comparator=partial(equal_ranges, delta=params.delta)))
    matching_events = intron_comparator.compare_junctions(introns1, mapped_region1, introns2, mapped_region2)

    logger.debug("Matching events detected: " + str(matching_events))
    if any(e in inconsistency_events for e in matching_events):
        logger.debug("Inconsistent")
        return MatchType.inconsistent
    elif len(introns1) != 0 and len(introns1) == len(introns2) and \
            contains(mapped_region1, (introns2[0][0], introns2[-1][1])) and \
            contains(mapped_region2, (introns1[0][0], introns1[-1][1])):
        logger.debug("Perfect match")
        return MatchType.match
    elif len(introns1) > 0 and len(introns2) == 0:
        logger.debug("First only")
        return MatchType.first_only
    elif len(introns2) > 0 and len(introns1) == 0:
        logger.debug("Second only")
        return MatchType.second_only
    elif len(introns1) > len(introns2) and contains(mapped_region1, (introns2[0][0], introns2[-1][1])):
        logger.debug("First longer")
        return MatchType.consistent_first_longer
    elif  len(introns1) < len(introns2) and contains(mapped_region2, (introns1[0][0], introns1[-1][1])):
        logger.debug("Second longer")
        return MatchType.consistent_second_longer
    logger.debug("Differ but consistent")
    return MatchType.consistent_differ


def check_annotated_introns(exons, gene_info, params):
    gene_region = (gene_info.start, gene_info.end)
    intron_profile_constructor = \
        OverlappingFeaturesProfileConstructor(gene_info.intron_profiles.features, gene_region,
                                              comparator=partial(equal_ranges, delta=params.delta),
                                              absence_condition=partial(overlaps_at_least,
                                                                        delta=params.minimal_intron_absence_overlap),
                                              delta=params.delta)
    intron_profile = intron_profile_constructor.construct_intron_profile(exons)
    return all(e == 1 for e in intron_profile.read_profile)


def intron_chain_stat(read_pairs, bam_records1, bam_records2, gene_db, params):
    stat_map = defaultdict(int)
    for read_pair in read_pairs:
        if read_pair[0] not in bam_records1 or read_pair[1] not in bam_records2:
            continue
        bam_record1 = bam_records1[read_pair[0]]
        bam_record2 = bam_records2[read_pair[1]]
        if bam_record1.reference_name != bam_record2.reference_name:
            continue
        exons1 = correct_bam_coords(concat_gapless_blocks(sorted(bam_record1.get_blocks()), bam_record1.cigartuples))
        exons2 = correct_bam_coords(concat_gapless_blocks(sorted(bam_record2.get_blocks()), bam_record2.cigartuples))
        mapped_region1 = (exons1[0][0], exons1[-1][1])
        mapped_region2 = (exons2[0][0], exons2[-1][1])
        if not overlaps(mapped_region1, mapped_region2):
            continue

        total_mapped_region = max_range(mapped_region1, mapped_region2)
        gene_list = list(gene_db.region(region=(bam_record1.reference_name, total_mapped_region[0], total_mapped_region[1]),
                                        completely_within=False))
        if not gene_list:
            continue
        logger.debug("First read exons: "  + str(exons1))
        logger.debug("Second read exons: " + str(exons2))
        gene_info = GeneInfo(gene_list, gene_db)
        logger.debug("Gene introns: " + str(gene_info.intron_profiles.features))
        first_annotated = check_annotated_introns(exons1, gene_info, params)
        second_annotated = check_annotated_introns(exons2, gene_info, params)
        logger.debug("First annotated %r, second annotated %r" % (first_annotated, second_annotated))
        event_type = compare_intron_chains(exons1, exons2, gene_info, params)
        stat_map[(event_type, first_annotated, second_annotated)] += 1

    return stat_map


def load_bam(read_set, bamfile):
    bam_records = {}
    for r in pysam.AlignmentFile(bamfile, "rb").fetch():
        if r.query_name in read_set:
            bam_records[r.query_name] = r
    return bam_records


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
    parser.add_argument("--reference", "-r", help="reference genome in FASTA format, "
                                                  "should be provided to compute some additional stats and "
                                                  "when raw reads are used as an input", type=str)

    parser.add_argument('--bam_pb', type=str, help='sorted and indexed BAM file for PacBio')
    parser.add_argument('--bam_ont', type=str, help='sorted and indexed BAM file for ONT')
    parser.add_argument('--tsv', type=str, help='TSV with barcode and read ids')
    parser.add_argument('--delta', type=int, default=0, help='delta')
    args = parser.parse_args(args)
    return args


def set_logger(logger_instance):
    logger_instance.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def run_pipeline(args):
    gene_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    read_pairs = load_tsv(args.tsv)
    bam_records1 = load_bam(set(map(lambda x: x[0], read_pairs)), args.bam_pb)
    bam_records2 = load_bam(set(map(lambda x: x[1], read_pairs)), args.bam_ont)
    params = Params()
    params.delta = args.delta
    stat_map = intron_chain_stat(read_pairs, bam_records1, bam_records2, gene_db, params)
    print(stat_map)

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
