
import os
import subprocess
import argparse
import pathlib

from collections import Counter

import pandas as pd
import numpy as np


def load_counts(inf, tpm_col=2, id_col=1):
    print("Loading TPM values from " + inf)
    tpm_dict = {}
    for l in open(inf):
        if l.startswith("#") or l.startswith("__") or l.startswith("feature_id"):
            continue
        v = l.strip().split()
        tpm_dict[v[id_col-1]] = float(v[tpm_col-1])
    return tpm_dict


def load_ref_ids_from_gtf(gtf):
    print("Loading annotation from " + gtf)
    id_dict = {}
    for l in open(gtf):
        if l.startswith("#"):
            continue
        v = l.strip().split()
        if v[2] != "transcript":
            continue
        tid_index = v.index("transcript_id", 7)
        ref_tid_index = v.index("reference_transcript_id", 7)
        id_dict[v[tid_index+1][1:-2]] = v[ref_tid_index+1][1:-2]
    return id_dict


def load_counts_from_gtf(gtf):
    print("Loading counts from " + gtf)
    tpm_dict = {}
    for l in open(gtf):
        if l.startswith("#"):
            continue
        v = l.strip().split()
        if v[2] != "transcript":
            continue
        tid_index = v.index("transcript_id", 7)
        tpm_index = v.index("TPM", 7)
        tpm_dict[v[tid_index+1][1:-2]] = float(v[tpm_index+1][1:-2])
    return tpm_dict


def load_tracking(inf):
    print("Loading tracking " + inf)
    id_dict = {}
    for l in open(inf):
        v = l.strip().split()
        if v[3] != '=':
            continue
        tid = v[4].split('|')[1]
        id_dict[tid] = v[2].split('|')[1]
    return id_dict


def correct_tpm_dict(tpm_dict, id_dict):
    print("Converting transcript ids")
    new_tpm_dict = {}
    for tid in id_dict.keys():
        if id_dict[tid] == 'novel':
            new_tpm_dict[tid] = tpm_dict[tid]
        else:
            new_tpm_dict[id_dict[tid]] = tpm_dict[tid]
    return new_tpm_dict


def count_deviation(df):
    print("Counting deviation histogram")
    deviation_values = []
    for index, row in df.iterrows():
        if row['ref_tpm'] == 0:
            continue
        deviation_values.append(100 * row['real_tpm'] / row['ref_tpm'])

    bins = [10 * i for i in range(41)]
    bins.append(10000)
    dev_vals, bins = np.histogram(deviation_values, bins)
    mid_bins = map(lambda x: x + 5, bins[:-1])
    return zip(mid_bins, dev_vals)


def count_stats(df):
    print('Correlation: ', round(np.corrcoef([df['real_tpm'], df['ref_tpm']])[1, 0], 3))
    full_matches = (df['real_tpm'] == df['ref_tpm']).astype(int).sum()
    n_isoforms = len(df['ref_tpm'])
    print('Full matches:', full_matches, 'Fraction:', round(full_matches / n_isoforms, 3))
    close_matches = ((df['real_tpm'] <= df['ref_tpm'] * 1.1) & (df['ref_tpm'] * 0.9 <= df['real_tpm'])).astype(int).sum()
    print('Close matches (10% diff):', close_matches, round(close_matches / n_isoforms, 2))
    close_matches = ((df['real_tpm'] <= df['ref_tpm'] * 1.2) & (df['ref_tpm'] * 0.8 <= df['real_tpm'])).astype(int).sum()
    print('Close matches (20% diff):', close_matches, round(close_matches / n_isoforms, 2))
    not_detected = (df['real_tpm'] == 0).astype(int).sum()
    print('Not detected:', not_detected, round(not_detected / n_isoforms, 2))
    false_detected = (df['ref_tpm'] == 0).astype(int).sum()
    print('False detections:', false_detected, round(false_detected / n_isoforms, 4))


def compare_transcript_counts(ref_tpm_dict, tpm_dict, output):
    print("Filling true values")
    joint_dict = {}
    for tid in tpm_dict.keys():
        if tid in ref_tpm_dict:
            joint_dict[tid] = (ref_tpm_dict[tid], tpm_dict[tid])
        else:
            joint_dict[tid] = (0, tpm_dict[tid])
    for tid in ref_tpm_dict:
        if tid not in joint_dict:
            joint_dict[tid] = (ref_tpm_dict[tid], 0)

    print("Converting to dataframe")
    df = pd.DataFrame.from_dict(joint_dict, orient='index', columns=['ref_tpm', 'real_tpm'])
    print("Saving TPM values")
    with open(os.path.join(output, "tpm.values.tsv"), 'w') as out_tpms:
        out_tpms.write("\t".join(map(str, list(df['ref_tpm']))) + "\n")
        out_tpms.write("\t".join(map(str, list(df['real_tpm']))) + "\n")

    with open(os.path.join(output, "deviation.tsv"), 'w') as out_dev:
        for hist_pairs in count_deviation(df):
            out_dev.write("%d\t%d\n" % (hist_pairs[0], hist_pairs[1]))

    count_stats(df)
    print("Done")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_expr', '-r', type=str, help='reference expression table, TPM', required=True)
    parser.add_argument('--ref_col', type=int, default=3, help='TPM column in reference expression table')
    parser.add_argument('--tpm', '-t', type=str, help='output expression table to assess, TPM')
    parser.add_argument('--tpm_col', type=int, default=2, help='TPM column in output expression table')
    parser.add_argument('--gtf', '-g', type=str, help='output GTF to convert reference transcript ids')
    parser.add_argument('--tracking', type=str, help='tracking file')

    parser.add_argument('--output', '-o', type=str, help='output folder', default="quantification_assessment")
    return parser.parse_args()


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    ref_tpm_dict = load_counts(args.ref_expr, args.ref_col)
    if args.tpm:
        tpm_dict = load_counts(args.tpm, args.tpm_col)
        if args.gtf:
            id_dict = load_ref_ids_from_gtf(args.gtf)
            tpm_dict = correct_tpm_dict(tpm_dict, id_dict)
    else:
        tpm_dict = load_counts_from_gtf(args.gtf)
        id_dict = load_tracking(args.tracking)
        tpm_dict = correct_tpm_dict(tpm_dict, id_dict)

    compare_transcript_counts(ref_tpm_dict, tpm_dict, args.output)


if __name__ == '__main__':
    main()
