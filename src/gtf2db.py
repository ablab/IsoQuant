#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import json
import logging
import os
import gffutils
import argparse
from traceback import print_exc
import gzip

logger = logging.getLogger('IsoQuant')


def db2gtf(db, gtf, _=None):
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


def db2bed(db, bed, _=None):
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

            if not exons:
                continue
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


def check_input_gtf(gtf, db, complete_db):
    logger.info("Checking input gene annotation")
    gtf_is_correct, corrected_gtf, out_fname, has_meta_features = check_gtf_duplicates(gtf)
    if not gtf_is_correct:
        outdir = os.path.dirname(db)
        new_gtf_path = os.path.join(outdir, out_fname)
        with open(new_gtf_path, "w") as out_gtf:
            out_gtf.write(corrected_gtf)
        logger.error("Input GTF seems to be corrupted (see warnings above).")
        logger.error("An attempt to correct this GTF was made, the result is written to %s" % new_gtf_path)
        logger.error("NB! some transcript / gene ids in the corrected annotation are modified.")
        logger.error("Provide a correct GTF by fixing the original input GTF or checking the corrected one.")
        exit(-3)
    else:
        logger.info("Gene annotation seems to be correct")

    if complete_db:
        if has_meta_features == -1:
            logger.error("You set --complete_genedb option but the provided annotation "
                         "does not contain any gene or transcript records.")
            logger.error("The results of this run will not be meaningful.")
            logger.error("Remove the output folder and restart IsoQuant without --complete_genedb.")
            exit(-3)
        elif has_meta_features == 0:
            logger.warning("You set --complete_genedb option but it looks like the provided annotation "
                           "does not contain all necessary gene or transcript records.")
            logger.warning("The results of this run might not be correct.")
            logger.warning("Remove the output folder and restart IsoQuant without --complete_genedb.")


def gtf2db(gtf, db, complete_db=False, check_gtf=True):
    if check_gtf:
        check_input_gtf(gtf, db, complete_db)

    logger.info("Converting gene annotation file to .db format (takes a while)...")
    gffutils.create_db(gtf, db, force=True, keep_order=True, merge_strategy='error',
                       sort_attribute_values=True, disable_infer_transcripts=complete_db,
                       disable_infer_genes=complete_db)
    logger.info("Gene database written to " + db)
    logger.info("Provide this database next time to avoid excessive conversion")


def convert_gtf_to_db(args):
    gtf_filename = args.genedb
    gtf_filename = os.path.abspath(gtf_filename)
    genedb_filename = args.genedb_filename
    gtf_filename, genedb_filename = convert_db(gtf_filename, genedb_filename, gtf2db, args)
    return genedb_filename


def convert_db_to_gtf(args):
    genedb_filename = os.path.abspath(args.genedb)
    gtf_filename = os.path.join(args.output, os.path.splitext(os.path.basename(genedb_filename))[0] + ".gtf")
    gtf_filename, genedb_filename = convert_db(gtf_filename, genedb_filename, db2gtf, args)
    return gtf_filename


def check_gtf_duplicates(gtf):
    gtf_correct = True
    line_count = 0
    gene_ids = {}
    transcript_ids = {}
    exon_gene_ids = set()
    exon_transcript_ids = set()
    corrected_gtf = ""

    gtf_name = os.path.basename(gtf)
    gtf_name, outer_ext = os.path.splitext(gtf_name)
    if outer_ext.lower() in ['.gz', '.gzip', '.bgz']:
        handle = gzip.open(gtf, "rt")
        gtf_name, inner_ext = os.path.splitext(gtf_name)
    else:
        handle = open(gtf, "rt")
        inner_ext = outer_ext

    for l in handle.readlines():
        line_count += 1
        if l.startswith("#"):
            corrected_gtf += l
            continue
        v = l.strip().split("\t")
        if len(v) < 9:
            corrected_gtf += l
            continue

        feature_type = v[2]
        attrs = v[8].split(" ")

        gene_id_pos = -1
        for i in range(len(attrs)):
            if attrs[i] == 'gene_id':
                gene_id_pos = i
        if gene_id_pos in [-1, len(attrs) - 1]:
            logger.warning("Malformed GTF line %d (gene_id attribute value cannot be found)" % line_count)
            logger.warning(l.strip())
            gtf_correct = False
            continue

        gene_str = attrs[gene_id_pos + 1]
        start_pos = gene_str.find('"')
        end_pos = gene_str.rfind('"')
        gene_id = gene_str[start_pos+1:end_pos]
        if feature_type == "gene":
            if gene_id in gene_ids:
                logger.warning("Duplicated gene id %s on line %d" % (gene_id, line_count))
                gtf_correct = False
                gene_ids[gene_id] += 1
                gene_id += ".%d" % gene_ids[gene_id]
            else:
                gene_ids[gene_id] = 0
        elif gene_id in gene_ids and gene_ids[gene_id] > 0:
            gene_id += ".%d" % gene_ids[gene_id]

        transcript_id_pos = -1
        for i in range(len(attrs)):
            if attrs[i] == 'transcript_id':
                transcript_id_pos = i
        if feature_type != "gene" and transcript_id_pos in [-1, len(attrs) - 1]:
            logger.warning("Malformed GTF line %d (transcript_id attribute value cannot be found)" % line_count)
            logger.warning(l.strip())
            gtf_correct = False
            continue

        transcript_str = attrs[transcript_id_pos + 1]
        start_pos = transcript_str.find('"')
        end_pos = transcript_str.rfind('"')
        transcript_id = transcript_str[start_pos+1:end_pos]
        if feature_type in ["transcript", "mRNA"]:
            if transcript_id in transcript_ids:
                logger.warning("Duplicated transcript id %s on line %d" % (transcript_id, line_count))
                gtf_correct = False
                transcript_ids[transcript_id] += 1
                transcript_id += ".%d" % transcript_ids[transcript_id]
            else:
                transcript_ids[transcript_id] = 0
        elif transcript_id in transcript_ids and transcript_ids[transcript_id] > 0:
            transcript_id += ".%d" % transcript_ids[transcript_id]

        if gene_id == transcript_id:
            logger.warning("Transcript id and gene id are identical (%s) at line %d"  % (transcript_id, line_count))
            gtf_correct = False
            transcript_id += ".RNA.IsoQuant_corrected"

        new_attrs = []
        for i in range(len(attrs)):
            if i == gene_id_pos + 1:
                new_attrs.append('"%s";' % gene_id)
            elif feature_type != "gene" and i == transcript_id_pos + 1:
                new_attrs.append('"%s";' % transcript_id)
            else:
                new_attrs.append(attrs[i])
        new_line = "\t".join(v[:8] + [" ".join(new_attrs)]) + "\n"
        corrected_gtf += new_line

        if feature_type == "exon":
            exon_gene_ids.add(gene_id)
            exon_transcript_ids.add(transcript_id)

    complete_genedb = 1
    if len(transcript_ids) == 0 or len(gene_ids) == 0:
        complete_genedb = -1
    elif len(transcript_ids) < len(exon_transcript_ids) or len(gene_ids) < len(exon_gene_ids):
        complete_genedb = 0

    return gtf_correct, corrected_gtf, gtf_name + ".corrected" + inner_ext.lower(), complete_genedb


def find_converted_db(converted_gtfs, gtf_filename, complete_genedb):
    gtf_mtime = converted_gtfs.get(gtf_filename, {}).get('gtf_mtime')
    db_mtime = converted_gtfs.get(gtf_filename, {}).get('db_mtime')
    db_file = converted_gtfs.get(gtf_filename, {}).get('genedb')
    is_complete = converted_gtfs.get(gtf_filename, {}).get('complete_db')
    if (os.path.exists(gtf_filename) and os.path.getmtime(gtf_filename) == gtf_mtime and
            os.path.exists(db_file) and os.path.getmtime(db_file) == db_mtime and complete_genedb == is_complete):
        return db_file
    return None


def compare_stored_gtf(converted_gtfs, gtf_filename, genedb_filename):
    gtf_mtime = converted_gtfs.get(gtf_filename, {}).get('gtf_mtime')
    db_mtime = converted_gtfs.get(gtf_filename, {}).get('db_mtime')
    return (os.path.exists(gtf_filename) and os.path.getmtime(gtf_filename) == gtf_mtime and
            os.path.exists(genedb_filename) and os.path.getmtime(genedb_filename) == db_mtime)


def convert_db(gtf_filename, genedb_filename, convert_fn, args):
    genedb_filename = os.path.abspath(genedb_filename)

    with open(args.db_config_path, 'r') as f_in:
        converted_gtfs = json.load(f_in)

    if not args.clean_start:
        if convert_fn == gtf2db:
            converted_db = find_converted_db(converted_gtfs, gtf_filename, args.complete_genedb)
            if converted_db is not None:
                logger.info("Gene annotation file found. Using " + converted_db)
                return gtf_filename, converted_db
        else:
            for converted_gtf in converted_gtfs:
                if compare_stored_gtf(converted_gtfs, converted_gtf, genedb_filename):
                    logger.info("Gene annotation file found. Using " + converted_gtf)
                    return converted_gtf, genedb_filename

    if convert_fn == gtf2db:
        convert_fn(gtf_filename, genedb_filename, args.complete_genedb, args.gtf_check)
    else:
        convert_fn(genedb_filename, gtf_filename)
    converted_gtfs[gtf_filename] = {
        'genedb': genedb_filename,
        'gtf_mtime': os.path.getmtime(gtf_filename),
        'db_mtime': os.path.getmtime(genedb_filename),
        'complete_db': args.complete_genedb
    }
    with open(args.db_config_path, 'w') as f_out:
        json.dump(converted_gtfs, f_out)
    return gtf_filename, genedb_filename


mode_dict = {"gtf2db" : gtf2db, "db2gtf" : db2gtf, "db2bed" : db2bed}


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file", required=True)
    parser.add_argument("--input", "-i", type=str, help="input gtf/db file", required=True)
    parser.add_argument('--complete_genedb', '-c', action='store_true', default=False,
                        help="use this flag if gene annotation contains transcript and gene metafeatures, "
                             "e.g. with official annotations, such as GENCODE; "
                             "speeds up gene database conversion")
    parser.add_argument("--mode", "-m", type=str, default="gtf2db", help="[gtf2db|db2gtf|db2bed]")


    args = parser.parse_args()
    if not args.output or not args.input:
        parser.print_usage()
        exit(-1)
    return args


def main():
    args = parse_args()
    if args.mode in mode_dict:
        conversion_func = mode_dict[args.mode]
    else:
        logger.error("Mode %s is not available" % args.mode)
        exit(-1)

    conversion_func(args.input, args.output, args.complete_genedb)


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
