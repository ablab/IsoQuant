#!/usr/bin/env python3
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
    print("TPMs loadede, checksum %.2f" % sum(tpm_dict.values()))
    return tpm_dict


def load_ref_ids_from_gtf(gtf, ref_keyword="reference_transcript_id"):
    print("Loading annotation from " + gtf)
    total_transcripts = 0
    known_transcripts = 0
    id_dict = {}
    for l in open(gtf):
        if l.startswith("#"):
            continue
        v = l.strip().split()
        if v[2] != "transcript":
            continue
        total_transcripts += 1
        tid_index = v.index("transcript_id", 7)
        try:
            ref_tid_index = v.index(ref_keyword, 7)
            ref_id = v[ref_tid_index+1][1:-2]
        except ValueError:
            ref_id = "novel"

        if ref_id != "novel":
            known_transcripts += 1
        id_dict[v[tid_index+1][1:-2]] = ref_id
    print("Total transcripts: %d, known: %d, novel: %d" %
          (total_transcripts, known_transcripts, total_transcripts-known_transcripts))
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
    total_transcripts = 0
    novel_transcripts = 0
    for l in open(inf):
        v = l.strip().split()
        tid = v[4].split('|')[1]
        total_transcripts += 1
        if v[3] == '=':
            id_dict[tid] = v[2].split('|')[1]
        elif tid not in id_dict:
            novel_transcripts += 1
            id_dict[tid] = 'novel'
    print("Total transcripts: %d, known: %d, novel: %d" %
          (total_transcripts, total_transcripts-novel_transcripts, novel_transcripts))

    return id_dict


def correct_tpm_dict(tpm_dict, id_dict, use_novel=True):
    print("Converting transcript ids")
    new_tpm_dict = {}
    for tid in tpm_dict.keys():
        if tid not in id_dict:
            print("WARN, %s not in dict" % tid)
            continue
        if id_dict[tid] == 'novel':
            if use_novel:
                new_tpm_dict[tid] = tpm_dict[tid]
        else:
            new_tpm_dict[id_dict[tid]] = tpm_dict[tid]
    print("Total values %d" % len(new_tpm_dict))
    return new_tpm_dict


def count_deviation(joint_dict):
    print("Counting deviation histogram")
    deviation_values = []
    false_detected = 0

    for t_id in joint_dict.keys():
        exp_pair = joint_dict[t_id]
        if exp_pair[0] == 0:
            if exp_pair[1] > 0:
                false_detected += 1
            continue
        deviation_values.append(100 * exp_pair[1] / exp_pair[0])

    print("Total %d, false %d, missed %d" % (len(deviation_values), false_detected, deviation_values.count(0.0)))
    bins = [10 * i for i in range(21)]
    bins.append(10000)
    dev_vals, bins = np.histogram(deviation_values, bins)
    mid_bins = bins[:-1]
    return zip(mid_bins, dev_vals)


def count_stats(joint_dict, output, header=""):
    ref_tpms = []
    real_tpms = []
    n_isoforms = 0
    counts_reported = 0
    full_matches = 0
    close_matches_10 = 0
    close_matches_20 = 0
    false_detected = 0
    not_detected = 0
    for t_id in joint_dict.keys():
        ref_expr = joint_dict[t_id][0]
        real_expr = joint_dict[t_id][1]
        if ref_expr == 0 and real_expr == 0:
            continue

        ref_tpms.append(ref_expr)
        real_tpms.append(real_expr)
        if ref_expr > 0:
            n_isoforms += 1
        if real_expr > 0:
            counts_reported += 1

        if real_expr == ref_expr:
            full_matches += 1
        if real_expr <= 1.1 * ref_expr and real_expr >= 0.9 * ref_expr and real_expr > 0:
            close_matches_10 += 1
        if real_expr <= 1.2 * ref_expr and real_expr >= 0.8 * ref_expr and real_expr > 0:
            close_matches_20 += 1
        if real_expr > 0 and ref_expr == 0:
            false_detected += 1
        if real_expr == 0 and ref_expr > 0:
            not_detected += 1

    outf = open(os.path.join(output, "stats.tsv"), 'w')
    outf.write(header + "\n")
    outf.write("Total reference isoforms %d\n" % n_isoforms)
    outf.write("Total non-zero counts %d\n" % counts_reported)
    outf.write('Correlation\t%.3f\n' % round(np.corrcoef([real_tpms, ref_tpms])[1, 0], 3))
    outf.write('Full matches\t%d\t%.3f\n' % (full_matches, round(full_matches / n_isoforms, 3)))
    outf.write('Close matches (10)\t%d\t%.3f\n' % (close_matches_10, round(close_matches_10 / n_isoforms, 3)))
    outf.write('Close matches (20)\t%d\t%.3f\n' % (close_matches_20, round(close_matches_20 / n_isoforms, 3)))
    outf.write('Not detected\t%d\t%.3f\n' % (not_detected, round(not_detected / n_isoforms, 3)))
    outf.write('False detections\t%d\t%.4f\n' % (false_detected, round(false_detected / n_isoforms, 4)))
    outf.close()


def compare_transcript_counts(ref_tpm_dict, tpm_dict, output, header=""):
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

    print("Saving TPM values")
    with open(os.path.join(output, "tpm.values.tsv"), 'w') as out_tpms:
        ref_tpms = []
        real_tpms = []
        for t_id in joint_dict.keys():
            if joint_dict[t_id] == (0, 0):
                continue
            ref_tpms.append(joint_dict[t_id][0])
            real_tpms.append(joint_dict[t_id][1])
        out_tpms.write("\t".join(map(str, ref_tpms)) + "\n")
        out_tpms.write("\t".join(map(str, real_tpms)) + "\n")

    with open(os.path.join(output, "deviation.tsv"), 'w') as out_dev:
        for hist_pairs in count_deviation(joint_dict):
            out_dev.write("%d\t%d\n" % (hist_pairs[0], hist_pairs[1]))

    count_stats(joint_dict, output, header)
    print("Done")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_expr', '-r', type=str, help='reference expression table, TPM', required=True)
    parser.add_argument('--ref_col', type=int, default=3, help='TPM column in reference expression table')
    parser.add_argument('--tpm', '-t', type=str, help='output expression table to assess, TPM')
    parser.add_argument('--tpm_col', type=int, default=2, help='TPM column in output expression table')
    parser.add_argument('--gtf', '-g', type=str, help='output GTF to convert reference transcript ids')
    parser.add_argument('--tracking', type=str, help='tracking file')
    parser.add_argument('--no_novel', action='store_false', default=True, help='do not use novel transcripts')

    parser.add_argument('--output', '-o', type=str, help='output folder', default="quantification_assessment")
    return parser.parse_args()


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    ref_tpm_dict = load_counts(args.ref_expr, args.ref_col)
    mode = ""
    if args.tpm:
        # IsoQuant tables
        tpm_dict = load_counts(args.tpm, args.tpm_col)
        mode = "IsoQuant"
        if args.gtf:
            # take reference ids from gtf
            id_dict = load_ref_ids_from_gtf(args.gtf)
            tpm_dict = correct_tpm_dict(tpm_dict, id_dict, args.no_novel)
            mode += " GTF"
        elif args.tracking:
            # take reference ids from gffcompare output .tracking
            id_dict = load_tracking(args.tracking)
            tpm_dict = correct_tpm_dict(tpm_dict, id_dict, args.no_novel)
            mode += " gffcompare"

    else:
        # StringTie
        tpm_dict = load_counts_from_gtf(args.gtf)
        mode = "StringTie"
        if args.tracking:
            # take reference ids from gffcompare output .tracking
            id_dict = load_tracking(args.tracking)
            tpm_dict = correct_tpm_dict(tpm_dict, id_dict, args.no_novel)
        else:
            # take reference ids from gtf
            id_dict = load_ref_ids_from_gtf(args.gtf, ref_keyword="transcript_id")
            tpm_dict = correct_tpm_dict(tpm_dict, id_dict, args.no_novel)

    compare_transcript_counts(ref_tpm_dict, tpm_dict, args.output)


if __name__ == '__main__':
    main()
