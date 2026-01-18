#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# takes teamcity config as input and copy etalon files

import os
import sys
import argparse
from traceback import print_exc
import shutil

from error_codes import IsoQuantExitCode


def load_tsv_config(config_file):
    config_dict = {}
    for l in open(config_file):
        if l.startswith("#"):
            continue

        tokens = l.strip().split('\t')
        if len(tokens) < 2:
            continue

        config_dict[tokens[0]] = tokens[1]
    return config_dict


def fix_path(config_file, path):
    if path.startswith('/'):
        return path

    return os.path.abspath(os.path.join(os.path.dirname(config_file), path))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--branch", type=str, help="branch to use [master]", default="master")
    parser.add_argument("--output", type=str, help="custom output folder (do not set if default)")
    parser.add_argument("config_file", metavar="config_file", type=str, nargs="+", help="configuration .info file")

    args = parser.parse_args()
    return args


etalon_dict = {"etalon": "gffcompare/new_gtf_etalon.tsv",
               "etalon_quantification_ref": "ref.quantification.tsv",
               "etalon_quantification_novel": "novel.quantification.tsv",
               "etalon_quantification_gene": "gene.quantification.tsv",
               "etalon_assignment": "new_assignment_etalon.tsv",
               "performance_baseline": "new_performance_baseline.tsv"}


def main():
    args = parse_args()
    if not args.config_file:
        print("Provide configuration file")
        sys.exit(IsoQuantExitCode.MISSING_REQUIRED_OPTION)

    for config_file in args.config_file:
        if not os.path.exists(config_file):
            print("Provide correct path to configuration file, %s does not exits" % config_file)
            sys.exit(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)

        print("Loading config from %s" % config_file)
        config_dict = load_tsv_config(config_file)
        run_name = config_dict["name"]
        output_folder = os.path.join(os.path.join(args.output if args.output else config_dict["output"], args.branch), run_name)
        for et in etalon_dict.keys():
            if et not in config_dict:
                continue
            new_etalon = os.path.join(output_folder, etalon_dict[et])
            old_etalon = fix_path(config_file, config_dict[et])
            print("Updating %s from %s" % (new_etalon, old_etalon))
            try:
                shutil.copy2(new_etalon, old_etalon)
            except FileNotFoundError:
                print("Warning! %s does not exist, skipping!" % new_etalon)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        ecode = main()
        if ecode != 0:
            sys.exit(ecode)
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(IsoQuantExitCode.UNCAUGHT_EXCEPTION)
