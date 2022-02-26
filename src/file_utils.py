import logging
import os
import re
import shutil
from collections import defaultdict

from src.read_groups import AbstractReadGrouper

logger = logging.getLogger('IsoQuant')


def merge_files(file_names, merged_file_name, stats_file_names=None, ignore_read_groups=True):
    file_names.sort(key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)])
    with open(merged_file_name, 'wb') as outf:
        for i, file_name in enumerate(file_names):
            if not os.path.exists(file_name): continue
            with open(file_name, 'rb') as f:
                if i != 0:
                    f.readline()
                shutil.copyfileobj(f, outf)
            os.remove(file_name)

    ambiguous_reads, not_assigned_reads, not_aligned_reads = 0, 0, 0
    if stats_file_names and ignore_read_groups:
        for file_name in stats_file_names:
            with open(file_name) as f:
                line = f.readline()
                ambiguous_reads += int(line.split()[1])
                line = f.readline()
                not_assigned_reads += int(line.split()[1])
                line = f.readline()
                not_aligned_reads += int(line.split()[1])
            os.remove(file_name)
        with open(merged_file_name, 'a') as outf:
            outf.write("__ambiguous\t%d\n" % ambiguous_reads)
            outf.write("__no_feature\t%d\n" % not_assigned_reads)
            outf.write("__not_aligned\t%d\n" % not_aligned_reads)


def convert_counts_to_tpm(counts_file_name, output_tpm_file_name, ignore_read_groups, output_zeroes):
    total_counts = defaultdict(float)

    print_in_columns = True
    with open(counts_file_name) as f:
        for i, line in enumerate(f):
            if line[0] == '_': break
            fs = line.split()
            if i == 0:
                if fs[1] == "group_id":
                    print_in_columns = False
                continue
            if ignore_read_groups:
                total_counts[AbstractReadGrouper.default_group_id] += float(fs[1])
            elif not ignore_read_groups and print_in_columns:
                for j in range(len(fs) - 1):
                    total_counts[j] += float(fs[j + 1])
            else:
                total_counts[fs[1]] += float(fs[2])

    scale_factors = {}
    for group_id in total_counts.keys():
        scale_factors[group_id] = 1000000.0 / total_counts[group_id] if total_counts[group_id] > 0 else 1.0
        logger.info("Scale factor for group %s = %.2f" % (group_id, scale_factors[group_id]))

    with open(output_tpm_file_name, "w") as outf:
        with open(counts_file_name) as f:
            for i, line in enumerate(f):
                fs = line.split()
                if fs[0] == '__ambiguous': break
                if i == 0:
                    outf.write(line.replace("count", "TPM"))
                    continue
                if ignore_read_groups:
                    feature_id, count = fs[0], float(fs[1])
                    tpm = scale_factors[AbstractReadGrouper.default_group_id] * count
                    if not output_zeroes and tpm == 0:
                        continue
                    outf.write("%s\t%.6f\n" % (feature_id, tpm))
                elif not print_in_columns:
                    for group_id in total_counts.keys():
                        feature_id, group_id, count = fs[0], fs[1], float(fs[2])
                        tpm = scale_factors[group_id] * count
                        outf.write("%s\t%s\t%.6f\n" % (feature_id, group_id, tpm))
                        fs = f.readline().split()
                else:
                    feature_id, counts = fs[0], list(map(float, fs[1:]))
                    tpm_values = [scale_factors[i] * counts[i] for i in range(len(scale_factors))]
                    outf.write("%s\t%s\n" % (feature_id, "\t".join(["%.6f" % c for c in tpm_values])))

