############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import json
import logging
import os
import sys
import gffutils
import argparse
from traceback import print_exc

logger = logging.getLogger('IsoQuant')


def db2gtf(gtf, db):
    logger.info("Converting gene annotation file to .gtf format (takes a while)...")
    with open(gtf, "w") as f:
        for record in gffutils.FeatureDB(db, keep_order=True).all_features():
            f.write(str(record) + '\n')
    logger.info("Gene database written to " + gtf)


def gtf2db(gtf, db, complete_db=False):
    logger.info("Converting gene annotation file to .db format (takes a while)...")
    gffutils.create_db(gtf, db, force=True, keep_order=True, merge_strategy='merge',
                       sort_attribute_values=True, disable_infer_transcripts=complete_db,
                       disable_infer_genes=complete_db)
    logger.info("Gene database written to " + db)
    logger.info("Provide this database next time to avoid excessive conversion")


def convert_gtf_to_db(args):
    gtf_filename = args.genedb
    gtf_filename = os.path.abspath(gtf_filename)
    genedb_filename = os.path.join(args.output, os.path.splitext(os.path.basename(gtf_filename))[0] + ".db")
    gtf_filename, genedb_filename = convert_db(gtf_filename, genedb_filename, gtf2db, args.complete_genedb)
    return genedb_filename


def convert_db_to_gtf(args):
    genedb_filename = os.path.abspath(args.genedb)
    gtf_filename = os.path.join(args.output, os.path.splitext(os.path.basename(genedb_filename))[0] + ".gtf")
    gtf_filename, genedb_filename = convert_db(gtf_filename, genedb_filename, db2gtf)
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


def convert_db(gtf_filename, genedb_filename, convert_fn, complete_db=False):
    genedb_filename = os.path.abspath(genedb_filename)
    config_dir = os.path.join(os.environ['HOME'], '.config', 'IsoQuant')
    config_path = os.path.join(config_dir, 'db_config.json')

    if not os.path.exists(config_path):
        os.makedirs(config_dir, exist_ok=True)
        converted_gtfs = {}
    else:
        with open(config_path, 'r') as f_in:
            converted_gtfs = json.load(f_in)
        if convert_fn == gtf2db:
            converted_db = find_coverted_db(converted_gtfs, gtf_filename)
            if converted_db is not None:
                logger.info("Gene annotation file was already converted to " + converted_db)
                return gtf_filename, converted_db
        else:
            for converted_gtf in converted_gtfs:
                if compare_stored_gtf(converted_gtfs, converted_gtf, genedb_filename):
                    logger.info("Gene annotation file was already converted")
                    return converted_gtf, genedb_filename

    if convert_fn == gtf2db:
        convert_fn(gtf_filename, genedb_filename, complete_db)
    else:
        convert_fn(gtf_filename, genedb_filename)
    converted_gtfs[gtf_filename] = {'genedb': genedb_filename,
                                    'gtf_mtime': os.path.getmtime(gtf_filename),
                                    'db_mtime': os.path.getmtime(genedb_filename),
                                    'complete_db': complete_db}
    with open(config_path, 'w') as f_out:
        json.dump(converted_gtfs, f_out)
    return gtf_filename, genedb_filename


def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: " + sys.argv[0] + " <GFF/GTF file> <output db file> \n")
        exit(0)
    gtf2db(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
