#!/usr/bin/env python3
#
# ############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import os
import argparse
from collections import defaultdict, OrderedDict
from traceback import print_exc
import pandas
import logging
import gzip
import gffutils

try:
    from .error_codes import IsoQuantExitCode
    from .file_naming import (
        mtx_matrix_file, mtx_features_file, mtx_barcodes_file,
        mtx_include_matrix_file, mtx_exclude_matrix_file,
    )
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from error_codes import IsoQuantExitCode
    from file_naming import (
        mtx_matrix_file, mtx_features_file, mtx_barcodes_file,
        mtx_include_matrix_file, mtx_exclude_matrix_file,
    )


GROUP_COUNT_CUTOFF = 100


logger = logging.getLogger('IsoQuant')


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def convert_to_matrix(input_linear_counts, output_file_path, feature_id_to_name_dict=None,
                      gzipped=False, feature_type='gene', convert_to_tpm=False, usable_reads_per_group=None,
                      max_groups: int = 0) -> int:
    """Convert linear counts to full matrix format.

    Returns:
        Number of groups found in the input file, or 0 on error.
    """
    logger.info("Converting %s to full matrix" % input_linear_counts)
    cols = len(open(input_linear_counts).readline().strip().split('\t'))
    if cols != 3:
        logger.error("Unexpected number of columns in %s: %d" % (input_linear_counts, cols))
        return 0

    # Quick group count check before expensive read+pivot
    if max_groups > 0:
        groups = set()
        with open(input_linear_counts) as f:
            next(f)  # skip header
            for line in f:
                parts = line.split('\t')
                if len(parts) >= 2:
                    groups.add(parts[1])
        num_groups = len(groups)
        if num_groups > max_groups:
            logger.info("Skipping full matrix conversion: %d groups exceeds limit of %d. "
                         "Use '--counts_format matrix' to force conversion." % (num_groups, max_groups))
            return num_groups

    df = pandas.read_csv(input_linear_counts, delimiter='\t', header=None, skiprows=1, keep_default_na=False,
                         names=['gene_id', 'group_id', 'count'],
                         dtype={'gene_id': str, 'group_id': str, 'count': float})
    count_matrix = df.pivot(index='gene_id', columns='group_id', values='count')

    normalization_factors = {g: 1.0 for g in df['group_id'].unique()} if not convert_to_tpm \
        else get_normalization_factors(df, usable_reads_per_group)

    num_groups = len(count_matrix.columns)

    output_file_path += ".tsv"
    output_file_path += ".gz" if gzipped else ""
    with gzip.open(output_file_path, 'wt') if gzipped else open(output_file_path, 'w') as outfile:
        # Write the header with group_ids
        columns = list(sorted(count_matrix.columns))
        if num_groups > GROUP_COUNT_CUTOFF:
            logger.warning("You have %d groups in your matrix, conversion might take a lot of time and the output file can be very large" % num_groups)
        outfile.write(feature_type + '_id\t' + '\t'.join(columns) + '\n')

        # Iterate over the DataFrame and write rows to the file
        for gene_id, row in count_matrix.iterrows():
            # Create a list of values, replacing NaNs with zeroes
            values = [str(row[col] * normalization_factors[col]) if pandas.notna(row[col]) else '0' for col in columns]
            gene_name = gene_id if feature_id_to_name_dict is None else feature_id_to_name_dict.get(gene_id, gene_id)
            outfile.write(gene_name + '\t' + '\t'.join(values) + '\n')
    logger.info("Matrix was saved to %s" % output_file_path)
    return num_groups


def convert_to_mtx(input_linear_counts, output_file_prefix, feature_id_to_name=None,
                   gzipped=False, convert_to_tpm=False, usable_reads_per_group=None):
    logger.info("Converting %s to MTX format" % input_linear_counts)
    cols = len(open(input_linear_counts).readline().strip().split('\t'))
    if cols != 3:
        logger.error("Unexpected number of columns in %s: %d" % (input_linear_counts, cols))
        return
    
    df = pandas.read_csv(input_linear_counts, delimiter='\t', header=None, skiprows=1, keep_default_na=False,
                         names=['gene_id', 'group_id', 'count'],
                         dtype={'gene_id': str, 'group_id': str, 'count': float})
    unique_genes = list(sorted(df['gene_id'].unique()))
    unique_groups = list(sorted(df['group_id'].unique()))
    gene_index_map = {gene_id: idx + 1 for idx, gene_id in enumerate(unique_genes)}
    group_index_map = {group_id: idx + 1 for idx, group_id in enumerate(unique_groups)}

    # Write the count matrix to an MTX file
    mtx_file = mtx_matrix_file(output_file_prefix)
    features_file = mtx_features_file(output_file_prefix)
    barcodes_files = mtx_barcodes_file(output_file_prefix)

    normalization_factors = {g: 1.0 for g in df['group_id'].unique()} if not convert_to_tpm \
        else get_normalization_factors(df, usable_reads_per_group)

    with open(barcodes_files, 'w') as bc_out:
        for group in unique_groups:
            bc_out.write(f"{group}\n")

    with open(features_file, 'w') as ft_out:
        for gene_id in unique_genes:
            gene_name = feature_id_to_name.get(gene_id, gene_id) if feature_id_to_name is not None else gene_id
            ft_out.write(f"{gene_id}\t{gene_name}\n")

    with gzip.open(mtx_file + ".gz", 'wt') if gzipped else open(mtx_file, 'w') as mtx_out:
        # Write the header
        mtx_out.write("%%MatrixMarket matrix coordinate real general\n")
        mtx_out.write(f"{len(unique_genes)} {len(unique_groups)} {df.shape[0]}\n")

        # Write the data in coordinate format
        for _, row in df.iterrows():
            gene_id = row['gene_id']
            group_id = row['group_id']
            count = row['count'] * normalization_factors[group_id]
            row_index = gene_index_map[gene_id]
            col_index = group_index_map[group_id]
            mtx_out.write(f"{row_index} {col_index} {count}\n")
    logger.info("Matrix was saved to 3 files with prefix %s" % output_file_prefix)


PROFILE_LINEAR_COLS = ['chr', 'start', 'end', 'strand', 'flags',
                       'gene_ids', 'group_id', 'include_counts', 'exclude_counts']


def _load_profile_linear(input_linear_counts: str):
    """Load 9-column exon/intron linear counts file into a DataFrame.

    Adds a synthetic ``feature_id`` column as ``chr:start-end:strand``.
    The ``flags`` (exon/intron type) column is preserved on the DataFrame
    but ignored by downstream matrix conversion.
    """
    cols = len(open(input_linear_counts).readline().strip().split('\t'))
    if cols != 9:
        logger.error("Unexpected number of columns in %s: %d (expected 9)" %
                     (input_linear_counts, cols))
        return None
    df = pandas.read_csv(input_linear_counts, delimiter='\t', header=None, skiprows=1,
                         keep_default_na=False, names=PROFILE_LINEAR_COLS,
                         dtype={'chr': str, 'start': int, 'end': int, 'strand': str,
                                'flags': str, 'gene_ids': str, 'group_id': str,
                                'include_counts': float, 'exclude_counts': float})
    df['feature_id'] = (df['chr'] + ':' + df['start'].astype(str) + '-' +
                        df['end'].astype(str) + ':' + df['strand'])
    return df


def _ordered_unique_features(df) -> tuple:
    """Return ordered unique feature_ids and a feature_id -> gene_ids map.

    Order matches first appearance in the input so the features file is
    deterministic across MTX files sharing it.
    """
    seen = OrderedDict()
    for feature_id, gene_ids in zip(df['feature_id'], df['gene_ids']):
        if feature_id not in seen:
            seen[feature_id] = gene_ids
    return list(seen.keys()), seen


def convert_profile_to_matrix(input_linear_counts: str, output_file_path: str,
                              gzipped: bool = False, feature_type: str = 'exon',
                              max_groups: int = 0) -> int:
    """Convert exon/intron linear counts to a single matrix TSV.

    Each cell contains ``include,exclude`` (comma-separated). Missing
    feature/group combinations become ``0,0``.

    Returns the number of groups in the input file, or 0 on error.
    """
    logger.info("Converting %s to full matrix" % input_linear_counts)
    df = _load_profile_linear(input_linear_counts)
    if df is None:
        return 0

    if max_groups > 0:
        num_groups = df['group_id'].nunique()
        if num_groups > max_groups:
            logger.info("Skipping full matrix conversion: %d groups exceeds limit of %d. "
                        "Use '--counts_format matrix' to force conversion." %
                        (num_groups, max_groups))
            return num_groups

    incl = df.pivot_table(index='feature_id', columns='group_id',
                          values='include_counts', aggfunc='sum', fill_value=0)
    excl = df.pivot_table(index='feature_id', columns='group_id',
                          values='exclude_counts', aggfunc='sum', fill_value=0)

    feature_order, _ = _ordered_unique_features(df)
    columns = sorted(set(incl.columns).union(excl.columns))
    incl = incl.reindex(index=feature_order, columns=columns, fill_value=0)
    excl = excl.reindex(index=feature_order, columns=columns, fill_value=0)
    num_groups = len(columns)

    output_file = output_file_path + ".tsv" + (".gz" if gzipped else "")
    if num_groups > GROUP_COUNT_CUTOFF:
        logger.warning("You have %d groups in your matrix, conversion might take a lot of time "
                       "and the output file can be very large" % num_groups)

    with gzip.open(output_file, 'wt') if gzipped else open(output_file, 'w') as outfile:
        outfile.write(feature_type + '_id\t' + '\t'.join(columns) + '\n')
        for feature_id in feature_order:
            incl_row = incl.loc[feature_id]
            excl_row = excl.loc[feature_id]
            cells = ["%d,%d" % (int(incl_row[c]), int(excl_row[c])) for c in columns]
            outfile.write(feature_id + '\t' + '\t'.join(cells) + '\n')
    logger.info("Matrix was saved to %s" % output_file)
    return num_groups


def convert_profile_to_mtx(input_linear_counts: str, output_file_prefix: str,
                           gzipped: bool = False) -> None:
    """Convert exon/intron linear counts to MTX format.

    Writes a shared features file and barcodes file plus two matrices
    (one for include_counts, one for exclude_counts).
    """
    logger.info("Converting %s to MTX format" % input_linear_counts)
    df = _load_profile_linear(input_linear_counts)
    if df is None:
        return

    feature_order, feature_to_gene = _ordered_unique_features(df)
    unique_groups = list(sorted(df['group_id'].unique()))
    feature_index_map = {feature_id: idx + 1 for idx, feature_id in enumerate(feature_order)}
    group_index_map = {group_id: idx + 1 for idx, group_id in enumerate(unique_groups)}

    features_file = mtx_features_file(output_file_prefix)
    barcodes_file = mtx_barcodes_file(output_file_prefix)
    include_file = mtx_include_matrix_file(output_file_prefix)
    exclude_file = mtx_exclude_matrix_file(output_file_prefix)

    with open(barcodes_file, 'w') as bc_out:
        for group in unique_groups:
            bc_out.write(f"{group}\n")

    with open(features_file, 'w') as ft_out:
        for feature_id in feature_order:
            ft_out.write(f"{feature_id}\t{feature_to_gene.get(feature_id, '')}\n")

    incl_nonzero = df[df['include_counts'] != 0]
    excl_nonzero = df[df['exclude_counts'] != 0]

    def _write_mtx(path, frame, value_col):
        with gzip.open(path + ".gz", 'wt') if gzipped else open(path, 'w') as out:
            out.write("%%MatrixMarket matrix coordinate real general\n")
            out.write(f"{len(feature_order)} {len(unique_groups)} {frame.shape[0]}\n")
            for _, row in frame.iterrows():
                out.write(f"{feature_index_map[row['feature_id']]} "
                          f"{group_index_map[row['group_id']]} "
                          f"{int(row[value_col])}\n")

    _write_mtx(include_file, incl_nonzero, 'include_counts')
    _write_mtx(exclude_file, excl_nonzero, 'exclude_counts')
    logger.info("MTX output saved with prefix %s (include + exclude matrices)" %
                output_file_prefix)


def get_normalization_factors(counts, usable_reads_per_group=None):
    if not usable_reads_per_group:
        total_group_counts = defaultdict(float)
        for _, row in counts.iterrows():
            group_id = row['group_id']
            count = row['count']
            total_group_counts[group_id] += count
    else:
        total_group_counts = usable_reads_per_group

    return {group_id: 1000000.0 / total_group_counts[group_id] for group_id in total_group_counts.keys()}


def load_gene_name_map(genedb):
    logger.info("Loading annotation from %s" % genedb)
    db = gffutils.FeatureDB(genedb, keep_order=True)
    transcript2gene = {}
    transcript_names = {}
    gene_names = {}

    for g in db.all_features(featuretype='gene'):
        gene_name = g.id
        if "gene_name" in g.attributes:
            gene_name = g["gene_name"][0]
        gene_names[g.id] = gene_name

    for t in db.all_features(featuretype=('transcript', 'mRNA')):
        if "gene_id" in t.attributes:
            gene_id = t["gene_id"][0]
        elif "Parent" in t.attributes:
            gene_id = t["Parent"][0]
        else:
            gene_id = "unknown_gene"
        transcript2gene[t.id] = (gene_id, gene_names.get(gene_id, gene_id))

        transcript_name = t.id
        if "transcript_name" in t.attributes:
            transcript_name = t["transcript_name"][0]
        transcript_names[t.id] = transcript_name

    logger.info("Loaded %d genes and %d transcripts" % (len(gene_names), len(transcript_names)))
    return gene_names, transcript_names, transcript2gene


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--input", "-i", type=str, help="counts in linear IsoQuant format", required=True)
    parser.add_argument("--genedb", "-g", type=str, help="gene annotation in .db format "
                                                         "(path to .db file can be found in the isoquant.log file), "
                                                         "feature names will be used instead of IDs if provided")
    parser.add_argument("--feature_type", help="feature type to be converted "
                                               "[gene, transcript, exon, intron]; "
                                               "annotation lookup applies only to gene/transcript",
                        default="gene", choices=['gene', 'transcript', 'exon', 'intron'])
    parser.add_argument("--output_format", "-f", type=str, choices=["mtx", "matrix"],
                        help="output format [mtx, matrix]", required=True)
    parser.add_argument("--tpm", help="convert counts to TPM", dest="convert_to_tpm",
                        action='store_true', default=False)
    parser.add_argument("--gzip", help="gzip output files", dest="gzipped",
                        action='store_true', default=False)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    set_logger(logger)

    is_profile_feature = args.feature_type in ('exon', 'intron')

    feature_name_dict = None
    if args.genedb and not is_profile_feature:
        gene_names, transcript_names, transcript2gene = load_gene_name_map(args.genedb)
        if args.feature_type == 'gene':
            feature_name_dict = gene_names
        elif args.feature_type == 'transcript':
            feature_name_dict = transcript_names
            with open(args.output + ".transcript2gene.tsv", "w") as out_mapf:
                for t in sorted(transcript2gene.keys()):
                    gene_info = transcript2gene[t]
                    out_mapf.write("%s\t%s\t%s\n" % (t, gene_info[0], gene_info[1]))

    if is_profile_feature:
        if args.convert_to_tpm:
            logger.warning("--tpm is ignored for exon/intron conversion")
        if args.output_format == "mtx":
            convert_profile_to_mtx(args.input, args.output, gzipped=args.gzipped)
        elif args.output_format == "matrix":
            convert_profile_to_matrix(args.input, args.output, gzipped=args.gzipped,
                                      feature_type=args.feature_type)
        else:
            logger.error("Unknown format %s" % args.output_format)
            sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
        return

    if args.output_format == "mtx":
        convert_to_mtx(args.input, args.output, feature_name_dict, gzipped=args.gzipped,
                       convert_to_tpm=args.convert_to_tpm)
    elif args.output_format == "matrix":
        convert_to_matrix(args.input, args.output, feature_name_dict, gzipped=args.gzipped,
                          feature_type=args.feature_type, convert_to_tpm=args.convert_to_tpm)
    else:
        logger.error("Unknown format %s" % args.output_format)
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(IsoQuantExitCode.UNCAUGHT_EXCEPTION)
