#!/usr/bin/env python3
import os
import argparse

import numpy as np


def load_counts(inf, tpm_col=2, id_col=1, novel=True):
    print("Loading TPM values from " + inf)
    tpm_dict = {}
    for l in open(inf):
        if l.startswith("#") or l.startswith("__") or l.startswith("feature_id"):
            continue
        v = l.strip().split()
        t_id = v[id_col-1]
        if not novel and t_id.find("nic") != -1:
            continue
        tpm_dict[t_id] = float(v[tpm_col-1])
    print("TPMs loaded, checksum %.2f" % sum(tpm_dict.values()))
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


def correct_tpm_dict(tpm_dict, id_dict):
    print("Converting transcript ids")
    new_tpm_dict = {}
    not_found = 0
    identical = 0
    substituted = 0
    for tid in tpm_dict.keys():
        if tid not in id_dict or id_dict[tid] == 'novel':
            new_tpm_dict[tid] = tpm_dict[tid]
            not_found += 1
        elif id_dict[tid] == tid:
            new_tpm_dict[tid] = tpm_dict[tid]
            identical += 1
        else:
            new_tpm_dict[id_dict[tid]] = tpm_dict[tid]
            substituted += 1
    print("Total values %d, substituted %d, identical %d, not found %d" % (len(new_tpm_dict), substituted, identical, not_found))
    return new_tpm_dict


def count_stats(joint_dict):
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

        if 1.1 * ref_expr >= real_expr >= 0.9 * ref_expr and real_expr > 0:
            close_matches_10 += 1
        if 1.2 * ref_expr >= real_expr >= 0.8 * ref_expr and real_expr > 0:
            close_matches_20 += 1

        if real_expr > 0 and ref_expr == 0:
            false_detected += 1
        if real_expr == 0 and ref_expr > 0:
            not_detected += 1

    corr_coeff = np.corrcoef([real_tpms, ref_tpms])[1, 0]
    return counts_reported, corr_coeff, close_matches_10, close_matches_20, not_detected, false_detected


def compare_transcript_counts(ref_tpm_dict, tpm_dict, overlap=False):
    joint_dict = {}
    for tid in tpm_dict.keys():
        if tid in ref_tpm_dict:
            joint_dict[tid] = (ref_tpm_dict[tid], tpm_dict[tid])
        elif not overlap:
            joint_dict[tid] = (0, tpm_dict[tid])
    for tid in ref_tpm_dict:
        if tid not in joint_dict and not overlap:
            joint_dict[tid] = (ref_tpm_dict[tid], 0)

    return count_stats(joint_dict)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_expr', '-r', type=str, help='reference expression table, TPM', required=True)
    parser.add_argument('--ref_col', type=int, default=3, help='TPM column in reference expression table')
    parser.add_argument('--tpm', '-t', type=str, help='output expression table to assess, TPM', required=True)
    parser.add_argument('--tpm_col', type=int, default=2, help='TPM column in output expression table')
    parser.add_argument('--tracking', type=str, help='tracking file')
    parser.add_argument('--no_novel', action='store_false', default=True, help='do not use novel transcripts')
    parser.add_argument('--full_stats', action='store_true', default=False, help='output full stats')
    parser.add_argument('--output', '-o', type=str, help='output file', default="quantification_assessment.tsv")
    return parser.parse_args()


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    print("Loading reference data from %s" % args.ref_expr)
    ref_tpm_dict = load_counts(args.ref_expr, args.ref_col)

    tpm_dict = load_counts(args.tpm, args.tpm_col, novel=not args.no_novel)

    if args.tracking:
        # take reference ids from gffcompare output .tracking
        id_dict = load_tracking(args.tracking)
        tpm_dict = correct_tpm_dict(tpm_dict, id_dict)

    outf = open(args.output, "w")
    suffix = "_overlap"
    counts_reported, corr_coeff, close_matches_10, close_matches_20, not_detected, false_detected = compare_transcript_counts(ref_tpm_dict, tpm_dict, overlap=True)
    outf.write("reported%s\t%d\n" % (suffix, counts_reported))
    outf.write("corr%s\t%.4f\n" % (suffix, corr_coeff))
    outf.write("good_matches%s\t%d\n" % (suffix, close_matches_10))
    outf.write("fair_matches%s\t%d\n" % (suffix, close_matches_20))
    outf.write("false_reports%s\t%d\n" % (suffix, false_detected))
    outf.write("missed%s\t%d\n" % (suffix, not_detected))

    suffix = "_full"
    counts_reported, corr_coeff, close_matches_10, close_matches_20, not_detected, false_detected = compare_transcript_counts(ref_tpm_dict, tpm_dict, overlap=False)
    outf.write("reported%s\t%d\n" % (suffix, counts_reported))
    outf.write("corr%s\t%.4f\n" % (suffix, corr_coeff))
    outf.write("good_matches%s\t%d\n" % (suffix, close_matches_10))
    outf.write("fair_matches%s\t%d\n" % (suffix, close_matches_20))
    outf.write("false_reports%s\t%d\n" % (suffix, false_detected))
    outf.write("missed%s\t%d\n" % (suffix, not_detected))
    outf.close()


if __name__ == '__main__':
    main()
