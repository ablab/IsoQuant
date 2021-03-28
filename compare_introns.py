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


inconsistency_events = {
    MatchEventSubtype.alt_left_site_novel, MatchEventSubtype.alt_right_site_novel,
    MatchEventSubtype.extra_intron_novel, MatchEventSubtype.mutually_exclusive_exons_novel,
    MatchEventSubtype.exon_gain_novel, MatchEventSubtype.exon_skipping_novel,
    MatchEventSubtype.exon_detach_novel, MatchEventSubtype.exon_merge_novel,
    MatchEventSubtype.terminal_exon_shift_novel, MatchEventSubtype.terminal_exon_misalignment_right,
    MatchEventSubtype.alternative_structure_novel, MatchEventSubtype.intron_alternation_novel,
    MatchEventSubtype.unspliced_intron_retention, MatchEventSubtype.intron_retention,
    MatchEventSubtype.alt_left_site_known, MatchEventSubtype.alt_right_site_known,
    MatchEventSubtype.extra_intron_known, MatchEventSubtype.intron_migration,
    MatchEventSubtype.mutually_exclusive_exons_known, MatchEventSubtype.exon_skipping_known,
    MatchEventSubtype.exon_detach_known, MatchEventSubtype.exon_merge_known,
    MatchEventSubtype.terminal_exon_shift_known, MatchEventSubtype.terminal_exon_misalignment_left,
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

    event_string = ",".join([match_subtype_to_str_with_additional_info(x, '+', introns1, introns2)
                             for x in matching_events])
    logger.debug("Matching events detected: " + event_string)
    if len(introns1) == 0 and len(introns2) == 0:
        logger.debug("Both unspliced")
        return MatchType.both_unspliced
    elif len(introns2) == 0:
        logger.debug("First only")
        return MatchType.only_first_spliced
    elif len(introns1) == 0:
        logger.debug("Second only")
        return MatchType.only_second_spliced
    elif any(e.event_type in inconsistency_events for e in matching_events):
        logger.debug("Inconsistent")
        return MatchType.inconsistent
    elif len(introns1) == len(introns2) and \
            contains(mapped_region1, (introns2[0][0], introns2[-1][1])) and \
            contains(mapped_region2, (introns1[0][0], introns1[-1][1])):
        logger.debug("Perfect match")
        return MatchType.match
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
    return intron_profile.read_profile.count(1)


def intron_chain_stat(read_pairs, bam_records1, bam_records2, gene_db, params):
    stat_map = defaultdict(int)
    unannotated_intron_stat = defaultdict(int)
    counter = 0
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
        first_annotated_count = check_annotated_introns(exons1, gene_info, params)
        second_annotated_count = check_annotated_introns(exons2, gene_info, params)
        first_intron_count = len(exons1) - 1
        second_intron_count = len(exons2) - 1
        first_annotated = first_annotated_count == first_intron_count
        second_annotated = second_annotated_count == second_intron_count
        logger.debug("First annotated %r, second annotated %r" % (first_annotated, second_annotated))
        event_type = compare_intron_chains(exons1, exons2, gene_info, params)
        stat_map[(event_type, first_annotated, second_annotated)] += 1
        unannotated_intron_stat["first", "annotated"] += first_annotated_count
        unannotated_intron_stat["first", "unannotated"] += first_intron_count - first_annotated_count
        unannotated_intron_stat["second", "annotated"] += second_annotated_count
        unannotated_intron_stat["second", "unannotated"] += second_intron_count - second_annotated_count
        counter += 1
        if counter % 1000 == 0:
            logger.info("Processed %d read pairs (%0.1f%%)" % (counter, 100 * counter / len(read_pairs)))

    return stat_map, unannotated_intron_stat


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
    parser.add_argument('--bam_pb', type=str, help='sorted and indexed BAM file for PacBio')
    parser.add_argument('--bam_ont', type=str, help='sorted and indexed BAM file for ONT')
    parser.add_argument('--tsv', type=str, help='TSV with barcode and read ids')
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
    logger.info("Loading read pairs from " + args.tsv)
    read_pairs = load_tsv(args.tsv)
    logger.info("Loading alignments from " + args.bam_pb)
    bam_records1 = load_bam(set(map(lambda x: x[0], read_pairs)), args.bam_pb)
    logger.info("Loading alignments from " + args.bam_ont)
    bam_records2 = load_bam(set(map(lambda x: x[1], read_pairs)), args.bam_ont)
    params = Params()
    params.delta = args.delta
    logger.info("Comparing intron chains...")
    stat_map, unannotated_intron_stat = intron_chain_stat(read_pairs, bam_records1, bam_records2, gene_db, params)
    logger.debug(str(stat_map))
    logger.info("Saving stats to " + args.output)
    outf = open(args.output, "w")
    for match_type in sorted(stat_map.keys()):
        outf.write("%s\t%r\t%r\t%d\n" % (match_type[0], match_type[1], match_type[2], stat_map[match_type]))
    outf.write("\n")
    for k in sorted(unannotated_intron_stat.keys()):
        outf.write("%s\t%s\t%d\n" % (k[0], k[1], unannotated_intron_stat[k]))
    outf.close()
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
