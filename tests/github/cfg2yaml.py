#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Convert .cfg + .etl/.qnt/.prf test configs to single .yaml files.

Usage:
    # Convert all pipeline configs
    python3 cfg2yaml.py /abga/work/andreyp/ci_isoquant/data/*.cfg

    # Convert barcode configs
    python3 cfg2yaml.py /abga/work/andreyp/ci_isoquant/data/barcodes/*.cfg

    # Dry run
    python3 cfg2yaml.py --dry-run /path/to/config.cfg

    # Barcode mode (etalon maps to baselines.barcode instead of baselines.transcripts)
    python3 cfg2yaml.py --barcode /path/to/barcode_config.cfg
"""

import argparse
import os
import sys

import yaml


# Baseline key mappings for pipeline configs
PIPELINE_BASELINE_KEYS: dict[str, str] = {
    "etalon": "transcripts",
    "etalon_assignment": "assignment",
    "etalon_quantification_ref": "quantification_ref",
    "etalon_quantification_novel": "quantification_novel",
    "etalon_quantification_gene": "quantification_gene",
    "performance_baseline": "performance",
    "etalon_allinfo": "allinfo",
}

# Baseline key mappings for barcode configs
BARCODE_BASELINE_KEYS: dict[str, str] = {
    "etalon": "barcode",
}


def load_tsv_config(config_file: str) -> dict[str, str]:
    """Load tab-separated key-value config file."""
    config_dict: dict[str, str] = {}
    for line in open(config_file):
        if line.startswith("#"):
            continue
        tokens = line.strip().split('\t')
        if len(tokens) < 2:
            continue
        config_dict[tokens[0]] = tokens[1]
    return config_dict


def fix_path(config_file: str, path: str) -> str:
    """Convert relative path to absolute using config file directory as base."""
    if path.startswith('/'):
        return path
    return os.path.abspath(os.path.join(os.path.dirname(config_file), path))


def load_baseline_file(filepath: str) -> dict[str, float]:
    """Load a baseline TSV file and return as dict of floats."""
    result: dict[str, float] = {}
    if not os.path.exists(filepath):
        print("  Warning: baseline file not found: %s" % filepath)
        return result
    for line in open(filepath):
        if line.startswith("#"):
            continue
        tokens = line.strip().split('\t')
        if len(tokens) < 2:
            continue
        try:
            result[tokens[0]] = float(tokens[1])
        except ValueError:
            result[tokens[0]] = tokens[1]
    return result


def detect_barcode_mode(config_dict: dict[str, str]) -> bool:
    """Detect if this is a barcode config based on its keys."""
    barcode_indicators = {"mode", "molecule", "barcodes", "qa_mode", "barcode_col"}
    pipeline_indicators = {"genome", "genedb", "datatype", "reduced_db"}

    has_barcode = bool(barcode_indicators & config_dict.keys())
    has_pipeline = bool(pipeline_indicators & config_dict.keys())

    # If it has barcode-specific keys but no pipeline keys, it's a barcode config
    if has_barcode and not has_pipeline:
        return True
    return False


def convert_cfg_to_yaml(cfg_path: str, force_barcode: bool = False) -> tuple[str, dict]:
    """Convert a .cfg file to YAML dict.

    Returns (yaml_path, yaml_data).
    """
    config_dict = load_tsv_config(cfg_path)
    if not config_dict:
        print("Error: empty or invalid config: %s" % cfg_path)
        return "", {}

    is_barcode = force_barcode or detect_barcode_mode(config_dict)
    baseline_keys = BARCODE_BASELINE_KEYS if is_barcode else PIPELINE_BASELINE_KEYS

    # Load baselines from external files
    baselines: dict[str, dict] = {}
    keys_to_remove: list[str] = []

    for cfg_key, yaml_section in baseline_keys.items():
        if cfg_key not in config_dict:
            continue
        baseline_path = fix_path(cfg_path, config_dict[cfg_key])
        baseline_data = load_baseline_file(baseline_path)
        if baseline_data:
            baselines[yaml_section] = baseline_data
        keys_to_remove.append(cfg_key)

    # Remove baseline keys from config (they're now embedded)
    for key in keys_to_remove:
        del config_dict[key]

    # Build YAML structure: config keys first, then baselines
    yaml_data: dict = {}
    for k, v in config_dict.items():
        yaml_data[k] = v

    if baselines:
        yaml_data["baselines"] = baselines

    yaml_path = os.path.splitext(cfg_path)[0] + ".yaml"
    return yaml_path, yaml_data


def main():
    parser = argparse.ArgumentParser(
        description="Convert .cfg + .etl/.qnt/.prf test configs to .yaml",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("config_files", nargs="+", help=".cfg files to convert")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print YAML output without writing files")
    parser.add_argument("--barcode", action="store_true",
                        help="Force barcode mode (etalon -> baselines.barcode)")
    args = parser.parse_args()

    for cfg_path in args.config_files:
        if not os.path.exists(cfg_path):
            print("Error: file not found: %s" % cfg_path)
            continue

        if not cfg_path.endswith(".cfg"):
            print("Skipping non-.cfg file: %s" % cfg_path)
            continue

        print("Converting: %s" % cfg_path)
        yaml_path, yaml_data = convert_cfg_to_yaml(cfg_path, force_barcode=args.barcode)
        if not yaml_data:
            continue

        yaml_output = yaml.dump(yaml_data, default_flow_style=False, sort_keys=False)

        if args.dry_run:
            print("--- %s ---" % yaml_path)
            print(yaml_output)
            print()
        else:
            with open(yaml_path, "w") as f:
                f.write(yaml_output)
            print("  -> %s" % yaml_path)


if __name__ == "__main__":
    main()
