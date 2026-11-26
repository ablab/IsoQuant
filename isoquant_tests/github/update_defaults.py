#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Updates CI baselines from the outputs of a finished CI run.
#
# YAML configs: values in the `baselines:` block are rewritten in place, line by line,
# so comments and the formatting of everything else in the file are preserved
# (re-dumping the whole YAML would rewrap long option strings and drop comments).
# Only metrics already listed in a baseline section are updated: the set of checked
# metrics never grows or shrinks implicitly.

import os
import re
import sys
import json
import argparse
from traceback import print_exc
import shutil
from typing import Dict, List, Optional, Set, Tuple

import yaml

from error_codes import IsoQuantExitCode


def load_tsv_config(config_file: str) -> Dict[str, str]:
    config_dict = {}
    for line in open(config_file):
        if line.startswith("#"):
            continue

        tokens = line.strip().split('\t')
        if len(tokens) < 2:
            continue

        config_dict[tokens[0]] = tokens[1]
    return config_dict


def fix_path(config_file: str, path: str) -> str:
    if path.startswith('/'):
        return path

    return os.path.abspath(os.path.join(os.path.dirname(config_file), path))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--branch", type=str, help="branch to use [master]", default="master")
    parser.add_argument("--output", type=str, help="custom output folder (do not set if default)")
    parser.add_argument("--sections", type=str,
                        help="comma-separated baseline sections to update (default: all present in the config)")
    parser.add_argument("--timing", action="store_true",
                        help="also update cpu_time/clock_time in performance baselines; by default only "
                             "max_rss is refreshed, timings need a run on an idle machine")
    parser.add_argument("--dry_run", action="store_true", help="print the changes without writing the configs")
    parser.add_argument("config_file", metavar="config_file", type=str, nargs="+", help="configuration .info file")

    args = parser.parse_args()
    return args


etalon_dict = {"etalon": "gffcompare/new_gtf_etalon.tsv",
               "etalon_quantification_ref": "ref.quantification.tsv",
               "etalon_quantification_novel": "novel.quantification.tsv",
               "etalon_quantification_gene": "gene.quantification.tsv",
               "etalon_assignment": "new_assignment_etalon.tsv",
               "performance_baseline": "new_performance_baseline.tsv"}

# YAML baseline section -> output file relative to the pipeline output folder
# ({output}/{branch}/{run_name}/)
YAML_PIPELINE_BASELINES: Dict[str, str] = {
    "transcripts": "gffcompare/new_gtf_etalon.tsv",
    "transcripts_td50": "gffcompare/new_gtf_etalon.tsv",
    "transcripts_td10": "gffcompare/new_gtf_etalon.tsv",
    "quantification_ref": "ref.quantification.tsv",
    "quantification_novel": "novel.quantification.tsv",
    "quantification_gene": "gene.quantification.tsv",
    "assignment": "new_assignment_etalon.tsv",
    "performance": "new_performance_baseline.tsv",
    "allinfo": "new_allinfo_etalon.tsv",
    "summary": "new_summary_etalon.tsv",
    "polya_prediction": "new_polya_prediction_etalon.tsv",
    "tss_prediction": "new_tss_prediction_etalon.tsv",
}

# run_pipeline.py writes all three gffcompare tolerances into one file, keyed "<section>.<metric>"
PREFIXED_SECTIONS: Set[str] = {"transcripts", "transcripts_td50", "transcripts_td10"}

TIMING_METRICS: Set[str] = {"cpu_time", "clock_time"}

SECTION_LINE = re.compile(r"^  ([A-Za-z0-9_]+):\s*$")
METRIC_LINE = re.compile(r"^(    )(\S.*?):(\s+)(\S.*?)\s*$")


def barcode_etalon_path(base_output_dir: str, run_name: str) -> str:
    # CI workflows use OUTPUT_BASE with a /barcodes/ suffix
    path = os.path.join(base_output_dir, "barcodes", run_name + ".new_etalon.tsv")
    if not os.path.exists(path):
        path = os.path.join(base_output_dir, run_name + ".new_etalon.tsv")
    return path


def load_section_values(section: str, base_output_dir: str, run_name: str) -> Optional[Dict[str, str]]:
    """New metric values of one baseline section, or None when there is no output to take them from."""
    if section == "barcode":
        path = barcode_etalon_path(base_output_dir, run_name)
    elif section in YAML_PIPELINE_BASELINES:
        path = os.path.join(base_output_dir, run_name, YAML_PIPELINE_BASELINES[section])
    else:
        print("Warning! No known output for baseline section %s, skipping" % section)
        return None
    if not os.path.exists(path):
        print("Warning! %s does not exist, skipping %s" % (path, section))
        return None

    values = load_tsv_config(path)
    if section in PREFIXED_SECTIONS:
        prefix = section + "."
        values = {k[len(prefix):]: v for k, v in values.items() if k.startswith(prefix)}
    return values


def same_value(old: str, new: str) -> bool:
    try:
        return float(old) == float(new)
    except ValueError:
        return yaml.safe_load(old) == new


def format_value(new: str) -> str:
    try:
        return repr(float(new))
    except ValueError:
        return json.dumps(new)


def update_metric_line(line: str, section: str, new_values: Dict[str, str],
                       timing: bool) -> Optional[Tuple[str, str, str, str]]:
    """Return (new line, metric, old value, new value) if this metric line of a section changes."""
    m = METRIC_LINE.match(line)
    if not m:
        return None
    indent, metric, sep, old = m.groups()
    if section == "performance" and metric in TIMING_METRICS and not timing:
        return None
    if metric not in new_values:
        print("Warning! %s.%s is not in the run output, keeping %s" % (section, metric, old))
        return None
    if same_value(old, new_values[metric]):
        return None
    new = format_value(new_values[metric])
    return "%s%s:%s%s" % (indent, metric, sep, new), metric, old, new


def update_yaml_baselines(config_file: str, run_name: str, base_output_dir: str,
                          sections: Optional[Set[str]] = None, timing: bool = False,
                          dry_run: bool = False) -> None:
    """Update values in the baselines block of a YAML config from the outputs of a run.

    Args:
        config_file: path to the YAML config to update in-place
        run_name: test name (used for locating output files)
        base_output_dir: {output}/{branch} directory
        sections: baseline sections to update, None for all
        timing: also update cpu_time/clock_time of the performance section
        dry_run: report changes without writing the file
    """
    with open(config_file) as f:
        lines = f.read().split("\n")

    starts = [i for i, line in enumerate(lines) if line.rstrip() == "baselines:"]
    if not starts:
        print("No baselines in %s, nothing to update" % config_file)
        return
    start = starts[0]

    new_values: Optional[Dict[str, str]] = None
    section = None
    changes: List[Tuple[str, str, str, str]] = []
    for i in range(start + 1, len(lines)):
        line = lines[i]
        if line and not line.startswith(" ") and not line.startswith("#"):
            break  # next top-level key
        m = SECTION_LINE.match(line)
        if m:
            section = m.group(1)
            selected = sections is None or section in sections
            new_values = load_section_values(section, base_output_dir, run_name) if selected else None
            continue
        if new_values is None:
            continue
        updated = update_metric_line(line, section, new_values, timing)
        if updated:
            lines[i], metric, old, new = updated
            changes.append((section, metric, old, new))

    for section, metric, old, new in changes:
        print("  %s.%s: %s -> %s" % (section, metric, old, new))
    if not changes:
        print("  no changes")
    if not dry_run and changes:
        with open(config_file, "w") as f:
            f.write("\n".join(lines))


def main() -> int:
    args = parse_args()
    if not args.config_file:
        print("Provide configuration file")
        sys.exit(IsoQuantExitCode.MISSING_REQUIRED_OPTION)
    sections = set(args.sections.split(",")) if args.sections else None

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
            update_yaml_baselines(config_file, run_name, base_output_dir, sections, args.timing, args.dry_run)
        else:
            config_dict = load_tsv_config(config_file)
            run_name = config_dict["name"]
            output_dir = args.output if args.output else config_dict["output"]
            output_folder = os.path.join(output_dir, args.branch, run_name)
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
    return 0


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        ecode = main()
        if ecode != 0:
            sys.exit(ecode)
    except SystemExit:
        raise
    except Exception:
        print_exc()
        sys.exit(IsoQuantExitCode.UNCAUGHT_EXCEPTION)
