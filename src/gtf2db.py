#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import json
import logging
import os
import gffutils
import argparse
from traceback import print_exc

logger = logging.getLogger('IsoQuant')


def db2gtf(gtf, db):
    logger.info("Converting gene annotation file to .gtf format (takes a while)...")
    with open(gtf, "w") as f:
        gene_db = gffutils.FeatureDB(db)
        for g in gene_db.features_of_type('gene', order_by=('seqid', 'start')):
            f.write(str(g) + '\n')
            for t in gene_db.children(g, featuretype=('transcript', 'mRNA')):
                f.write(str(t) + '\n')
                for e in gene_db.children(t, featuretype=('exon', 'CDS')):
                    f.write(str(e) + '\n')
    logger.info("Gene database written to " + gtf)


def db2bed(bed, db):
    # adapted from paftools in mimimap2
    colors = {
        'protein_coding': '0,128,255',
        'mRNA': '0,128,255',
        'lincRNA': '0,192,0',
        'snRNA': '0,192,0',
        'miRNA': '0,192,0',
        'misc_RNA': '0,192,0'
    }

    def get_color(transcript_kind):
        if transcript_kind in colors:
            return colors[transcript_kind]
        return '196,196,196'

    logger.info("Converting gene annotation file %s to .bed format" % db)
    genedb = gffutils.FeatureDB(db, keep_order=True)
    with open(bed, "w") as f:
        for record in genedb.all_features(featuretype=('transcript', 'mRNA')):
            if "transcript_type" in record.attributes:
                transcript_type = record["transcript_type"][0]
            else:
                transcript_type = ""
            if "gene_name" in record.attributes:
                gene_name = record["gene_name"][0]
            elif "gene_id" in record.attributes:
                gene_name = record["gene_id"][0]
            else:
                gene_name = "unknown_gene"
            transcript_name = record.id + "|" + transcript_type + "|" + gene_name
            exons = []
            cds = []
            for e in genedb.children(record, order_by='start', featuretype='exon'):
                exons.append((e.start - 1, e.end))
            for e in genedb.children(record, order_by='start', featuretype='CDS'):
                cds.append((e.start - 1, e.end))

            transcript_start = exons[0][0]
            transcript_end = exons[-1][1]
            exons_lens = []
            exon_start = []
            for e in exons:
                exons_lens.append(e[1] - e[0])
                exon_start.append(e[0] - transcript_start)

            if cds:
                thick_start, thick_end = cds[0][0], cds[-1][1]
            else:
                thick_start, thick_end = transcript_start, transcript_end
            exon_lengths = ",".join([str(x) for x in exons_lens]) + ","
            exon_starts = ",".join([str(x) for x in exon_start]) + ","
            line = "%s\t%d\t%d\t%s\t1000\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n" % \
                   (record.seqid, transcript_start, transcript_end, transcript_name, record.strand,
                    thick_start, thick_end, get_color(transcript_type), len(exons_lens), exon_lengths, exon_starts)
            f.write(line)
    logger.info("Gene database BED written to " + bed)


def gtf2db(gtf, db, complete_db=False):
    logger.info("Converting gene annotation file to .db format (takes a while)...")
    gffutils.create_db(gtf, db, force=True, keep_order=True, merge_strategy='error',
                       sort_attribute_values=True, disable_infer_transcripts=complete_db,
                       disable_infer_genes=complete_db)
    logger.info("Gene database written to " + db)
    logger.info("Provide this database next time to avoid excessive conversion")


def convert_gtf_to_db(args, output_is_dir=True):
    gtf_filename = args.genedb
    gtf_filename = os.path.abspath(gtf_filename)
    output_path =  args.output if args.genedb_output is None else args.genedb_output
    if output_is_dir:
        genedb_filename = os.path.join(output_path, os.path.splitext(os.path.basename(gtf_filename))[0] + ".db")
    else:
        genedb_filename = output_path + "." + os.path.splitext(os.path.basename(gtf_filename))[0] + ".db"
    gtf_filename, genedb_filename = convert_db(gtf_filename, genedb_filename, gtf2db, args)
    return genedb_filename


def convert_db_to_gtf(args):
    genedb_filename = os.path.abspath(args.genedb)
    gtf_filename = os.path.join(args.output, os.path.splitext(os.path.basename(genedb_filename))[0] + ".gtf")
    gtf_filename, genedb_filename = convert_db(gtf_filename, genedb_filename, db2gtf, args)
    return gtf_filename


def find_coverted_db(converted_gtfs, gtf_filename):
    gtf_mtime = converted_gtfs.get(gtf_filename, {}).get('gtf_mtime')
    db_mtime = converted_gtfs.get(gtf_filename, {}).get('db_mtime')
    db_file = converted_gtfs.get(gtf_filename, {}).get('genedb')
    if os.path.exists(gtf_filename) and os.path.getmtime(gtf_filename) == gtf_mtime:
        if os.path.exists(db_file) and os.path.getmtime(db_file) == db_mtime:
            return db_file
    return None


def compare_stored_gtf(converted_gtfs, gtf_filename, genedb_filename):
    gtf_mtime = converted_gtfs.get(gtf_filename, {}).get('gtf_mtime')
    db_mtime = converted_gtfs.get(gtf_filename, {}).get('db_mtime')
    if os.path.exists(gtf_filename) and os.path.getmtime(gtf_filename) == gtf_mtime:
        if os.path.exists(genedb_filename) and os.path.getmtime(genedb_filename) == db_mtime:
            return True


def convert_db(gtf_filename, genedb_filename, convert_fn, args):
    genedb_filename = os.path.abspath(genedb_filename)

    with open(args.db_config_path, 'r') as f_in:
        converted_gtfs = json.load(f_in)

    if not args.clean_start:
        if convert_fn == gtf2db:
            converted_db = find_coverted_db(converted_gtfs, gtf_filename)
            if converted_db is not None:
                logger.info("Gene annotation file found. Using " + converted_db)
                return gtf_filename, converted_db
        else:
            for converted_gtf in converted_gtfs:
                if compare_stored_gtf(converted_gtfs, converted_gtf, genedb_filename):
                    logger.info("Gene annotation file found. Using " + converted_gtf)
                    return converted_gtf, genedb_filename

    if convert_fn == gtf2db:
        convert_fn(gtf_filename, genedb_filename, args.complete_genedb)
    else:
        convert_fn(gtf_filename, genedb_filename)
    converted_gtfs[gtf_filename] = {
        'genedb': genedb_filename,
        'gtf_mtime': os.path.getmtime(gtf_filename),
        'db_mtime': os.path.getmtime(genedb_filename),
        'complete_db': args.complete_genedb
    }
    with open(args.db_config_path, 'w') as f_out:
        json.dump(converted_gtfs, f_out)
    return gtf_filename, genedb_filename

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file", default="")
    parser.add_argument("--input", "-i", type=str, help="input gtf/db file")
    parser.add_argument('--complete_genedb', '-c', action='store_true', default=False,
                        help="use this flag if gene annotation contains transcript and gene metafeatures, "
                             "e.g. with official annotations, such as GENCODE; "
                             "speeds up gene database conversion")
    parser.add_argument("--reverse", "-r", action='store_true', default=False, help="db2gtf")


    args = parser.parse_args()
    if not args.output or not args.input:
        parser.print_usage()
        exit(-1)
    return args


def main():
    args = parse_args()
    if args.reverse:
        db2gtf(args.output, args.input)
    else:
        gtf2db(args.input, args.output, args.complete_genedb)


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
