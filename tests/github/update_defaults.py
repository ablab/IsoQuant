#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# takes teamcity config as input and copy etalon files

import os
import sys
import argparse
from traceback import print_exc
import shutil

import yaml

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

# Maps YAML baseline section name -> output file path relative to pipeline output folder
# ({output}/{branch}/{run_name}/)
YAML_PIPELINE_BASELINES: dict[str, str] = {
    "transcripts": "gffcompare/new_gtf_etalon.tsv",
    "quantification_ref": "ref.quantification.tsv",
    "quantification_novel": "novel.quantification.tsv",
    "quantification_gene": "gene.quantification.tsv",
    "assignment": "new_assignment_etalon.tsv",
    "performance": "new_performance_baseline.tsv",
    "allinfo": "new_allinfo_etalon.tsv",
}


def update_yaml_baselines(config_file: str, run_name: str, base_output_dir: str) -> None:
    """Update baselines section in a YAML config file from run output.

    Args:
        config_file: path to the YAML config to update in-place
        run_name: test name (used for locating output files)
        base_output_dir: {output}/{branch} directory
    """
    with open(config_file) as f:
        data = yaml.safe_load(f)

    if "baselines" not in data:
        data["baselines"] = {}

    pipeline_output_folder = os.path.join(base_output_dir, run_name)

    for section, output_rel_path in YAML_PIPELINE_BASELINES.items():
        if section not in data["baselines"]:
            continue

        new_values_file = os.path.join(pipeline_output_folder, output_rel_path)
        if not os.path.exists(new_values_file):
            print("Warning! %s does not exist, skipping %s!" % (new_values_file, section))
            continue

        new_values = load_tsv_config(new_values_file)
        updated = {}
        for k, v in new_values.items():
            try:
                updated[k] = float(v)
            except ValueError:
                updated[k] = v

        print("Updating baselines.%s from %s" % (section, new_values_file))
        data["baselines"][section] = updated

    # Barcode baselines: new etalon is at {base_output_dir}/{run_name}.new_etalon.tsv
    if "barcode" in data["baselines"]:
        barcode_etalon = os.path.join(base_output_dir, run_name + ".new_etalon.tsv")
        if not os.path.exists(barcode_etalon):
            print("Warning! %s does not exist, skipping barcode!" % barcode_etalon)
        else:
            new_values = load_tsv_config(barcode_etalon)
            updated = {}
            for k, v in new_values.items():
                try:
                    updated[k] = float(v)
                except ValueError:
                    updated[k] = v
            print("Updating baselines.barcode from %s" % barcode_etalon)
            data["baselines"]["barcode"] = updated

    with open(config_file, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


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
        ext = os.path.splitext(config_file)[1].lower()

        if ext in (".yaml", ".yml"):
            with open(config_file) as f:
                data = yaml.safe_load(f)
            run_name = str(data["name"])
            output_dir = args.output if args.output else str(data.get("output", ""))
            base_output_dir = os.path.join(output_dir, args.branch)
            update_yaml_baselines(config_file, run_name, base_output_dir)
        else:
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
