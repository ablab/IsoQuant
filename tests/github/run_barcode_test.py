#!/usr/bin/env python3
############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Test runner for barcode calling quality assessment.

Takes a configuration file as input, runs isoquant_detect_barcodes.py on simulated data,
then runs quality assessment and compares results against etalon values.

Usage:
    python run_barcode_test.py test.cfg -o /output/dir
"""

import glob
import os
import shutil
import sys
import argparse
import logging
import subprocess
from traceback import print_exc


# Run types
RT_VOID = "void"           # Just run detect_barcodes, no quality check
RT_SINGLE = "single"       # Single barcode per read (10x, curio, stereo)
RT_STEREO_SPLIT = "stereo_split"  # Multiple barcodes per read (stereo split)

log = logging.getLogger('BarcodeTestRunner')


def set_logger(args, logger_instance):
    if "debug" not in args.__dict__ or not args.debug:
        output_level = logging.INFO
    else:
        output_level = logging.DEBUG
    logger_instance.setLevel(output_level)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def load_tsv_config(config_file):
    """Load tab-separated key-value config file."""
    if not os.path.exists(config_file):
        log.error("Config file %s was not found" % config_file)
        return {}
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
    """Convert relative path to absolute using config file directory as base."""
    if path.startswith('/'):
        return path
    return os.path.abspath(os.path.join(os.path.dirname(config_file), path))


def fix_paths(config_file, paths):
    """Convert multiple space-separated paths."""
    new_paths = []
    for p in paths.split(' '):
        new_paths.append(fix_path(config_file, p))
    return new_paths


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run barcode calling quality tests',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("config_file", metavar="config_file", type=str,
                        help="Configuration .cfg file")
    parser.add_argument("--output", "-o", type=str,
                        help="Output directory (overrides config)")
    parser.add_argument("--additional_options", "-a", type=str,
                        help="Additional options for isoquant_detect_barcodes.py")
    parser.add_argument("--debug", "-d", action="store_true",
                        help="Enable debug output")

    return parser.parse_args()


def run_detect_barcodes(args, config_dict):
    """Run isoquant_detect_barcodes.py with configuration options."""
    source_dir = os.path.dirname(os.path.realpath(__file__))
    isoquant_dir = os.path.join(source_dir, "../../")
    config_file = args.config_file

    run_name = config_dict["name"]
    # isoquant_detect_barcodes.py uses -o as a file prefix, not a directory
    # Output will be: {output_folder}/{run_name}.barcoded_reads.tsv
    output_folder = args.output if args.output else config_dict["output"]
    output_prefix = os.path.join(output_folder, run_name)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    log.info('== Running isoquant_detect_barcodes.py ==')

    # Build command
    detect_command = [
        "python3", os.path.join(isoquant_dir, "isoquant_detect_barcodes.py"),
        "-o", output_prefix,
        "--mode", config_dict["mode"]
    ]

    # Add reads
    if "reads" in config_dict:
        reads = fix_paths(config_file, config_dict["reads"])
        detect_command.append("--input")
        detect_command.extend(reads)
    elif "bam" in config_dict:
        bams = fix_paths(config_file, config_dict["bam"])
        detect_command.append("--input")
        detect_command.extend(bams)

    # Add barcodes whitelist(s)
    # For visium_hd and other modes that need multiple barcode files,
    # separate them with spaces in the config
    if "barcodes" in config_dict:
        barcode_files = fix_paths(config_file, config_dict["barcodes"])
        detect_command.extend(["--barcodes"] + barcode_files)

    # Add molecule description file (for custom_sc mode)
    if "molecule" in config_dict:
        molecule_file = fix_path(config_file, config_dict["molecule"])
        detect_command.extend(["--molecule", molecule_file])

    # Add optional parameters
    if "threads" in config_dict:
        detect_command.extend(["-t", config_dict["threads"]])

    # Add any additional detect_barcodes options from config
    if "detect_options" in config_dict:
        opts = config_dict["detect_options"].replace('"', '').split()
        detect_command.extend(opts)

    # Add command line additional options
    if args.additional_options:
        detect_command.extend(args.additional_options.split())

    log.info("Command line: " + " ".join(detect_command))
    result = subprocess.run(detect_command)

    if result.returncode != 0:
        log.error("isoquant_detect_barcodes.py exited with non-zero status: %d" % result.returncode)
        return -11

    return 0


def run_barcode_quality(args, config_dict):
    """Run barcode quality assessment and compare to etalon."""
    source_dir = os.path.dirname(os.path.realpath(__file__))
    isoquant_dir = os.path.join(source_dir, "../../")
    config_file = args.config_file

    run_name = config_dict["name"]
    output_folder = args.output if args.output else config_dict["output"]
    output_prefix = os.path.join(output_folder, run_name)

    # Find barcode output file
    # isoquant_detect_barcodes.py outputs: {prefix}.barcoded_reads.tsv
    barcode_output = output_prefix + ".barcoded_reads.tsv"
    if not os.path.exists(barcode_output):
        barcode_output = output_prefix + ".barcoded_reads.tsv.gz"
    if not os.path.exists(barcode_output):
        # Try glob patterns as fallback
        patterns = [
            output_prefix + "*.barcoded_reads.tsv*",
            os.path.join(output_folder, "*.barcoded_reads.tsv"),
        ]
        barcode_output = None
        for pattern in patterns:
            candidates = glob.glob(pattern)
            if candidates:
                barcode_output = candidates[0]
                break

    if not barcode_output or not os.path.exists(barcode_output):
        # List what files are actually in the directory
        all_files = os.listdir(output_folder) if os.path.exists(output_folder) else []
        log.error("Barcode output file not found for prefix %s" % output_prefix)
        log.error("Files in directory: %s" % (", ".join(all_files) if all_files else "(empty)"))
        return -20

    log.info("Found barcode output: %s" % barcode_output)

    log.info('== Running quality assessment ==')

    # Determine QA mode
    run_type = config_dict.get("run_type", RT_SINGLE)
    if run_type == RT_STEREO_SPLIT:
        qa_mode = "stereo_split"
    else:
        # Use mode from config or derive from detect mode
        qa_mode = config_dict.get("qa_mode", config_dict.get("mode", "stereo"))
        # Map detect_barcodes modes to QA modes
        mode_mapping = {
            "stereoseq": "stereo",
            "tenX_v3": "tenX_v3",
            "visium_hd": "visium_hd",
            "curio": "curio",
        }
        qa_mode = mode_mapping.get(qa_mode, qa_mode)

    quality_report = output_prefix + ".quality_report.tsv"
    qa_command = [
        "python3", os.path.join(isoquant_dir, "misc/assess_barcode_quality.py"),
        "--mode", qa_mode,
        "--input", barcode_output,
        "--output", quality_report
    ]

    # Add optional parameters
    if "barcode_col" in config_dict:
        qa_command.extend(["--barcode_col", config_dict["barcode_col"]])
    if "umi_col" in config_dict:
        qa_command.extend(["--umi_col", config_dict["umi_col"]])
    if "score_col" in config_dict:
        qa_command.extend(["--score_col", config_dict["score_col"]])
    if "barcode_length" in config_dict:
        qa_command.extend(["--barcode_length", config_dict["barcode_length"]])

    log.info("QA command line: " + " ".join(qa_command))
    result = subprocess.run(qa_command)

    if result.returncode != 0:
        log.error("Quality assessment exited with non-zero status: %d" % result.returncode)
        return -21

    # Check against etalon if provided
    if "etalon" not in config_dict:
        log.info("No etalon file specified, skipping comparison")
        return 0

    return check_against_etalon(config_file, config_dict, output_prefix, quality_report)


def check_against_etalon(config_file, config_dict, output_prefix, quality_report):
    """Compare quality metrics against etalon values."""
    log.info('== Checking quality metrics against etalon ==')

    etalon_file = fix_path(config_file, config_dict["etalon"])
    if not os.path.exists(etalon_file):
        log.error("Etalon file not found: %s" % etalon_file)
        return -22

    etalon_dict = load_tsv_config(etalon_file)
    quality_dict = load_tsv_config(quality_report)

    exit_code = 0

    # Write new etalon file for easy updates
    new_etalon_file = output_prefix + ".new_etalon.tsv"
    with open(new_etalon_file, "w") as f:
        for k, v in sorted(quality_dict.items()):
            try:
                f.write("%s\t%.2f\n" % (k, float(v)))
            except ValueError:
                f.write("%s\t%s\n" % (k, v))

    # Compare each metric in etalon
    for metric_name, etalon_value in etalon_dict.items():
        if metric_name not in quality_dict:
            log.error("Metric %s not found in quality report" % metric_name)
            exit_code = -23
            continue

        try:
            etalon_val = float(etalon_value)
            actual_val = float(quality_dict[metric_name])
        except ValueError:
            log.warning("Could not parse numeric value for %s" % metric_name)
            continue

        # Get tolerance from config (default 1%)
        tolerance = float(config_dict.get("tolerance", "0.01"))
        err_code = check_value(etalon_val, actual_val, metric_name, tolerance)
        if err_code != 0:
            exit_code = err_code

    return exit_code


def check_value(etalon_value, output_value, name, percent=0.01):
    """
    Check if output value is within tolerance of etalon value.

    Args:
        etalon_value: Expected value
        output_value: Actual value
        name: Metric name for logging
        percent: Tolerance as fraction (default 0.01 = 1%)

    Returns:
        0 if within tolerance, -2 otherwise
    """
    lower_bound = etalon_value * (1 - percent)
    upper_bound = etalon_value * (1 + percent)
    exit_code = 0

    if output_value < 0:
        # Negative values must match exactly
        if output_value != etalon_value:
            log.error("Value of %s = %.2f is not equal to %.2f" %
                      (name, output_value, etalon_value))
            exit_code = -2
        else:
            log.info("Value of %s = %.2f == %.2f as expected" %
                     (name, output_value, etalon_value))
    else:
        if output_value < lower_bound:
            log.error("Value of %s = %.2f is lower than expected %.2f (min: %.2f)" %
                      (name, output_value, etalon_value, lower_bound))
            exit_code = -2
        elif output_value > upper_bound:
            log.error("Value of %s = %.2f is higher than expected %.2f (max: %.2f)" %
                      (name, output_value, etalon_value, upper_bound))
            exit_code = -2
        else:
            log.info("Value of %s = %.2f is within tolerance of %.2f [%.2f, %.2f]" %
                     (name, output_value, etalon_value, lower_bound, upper_bound))

    return exit_code


def check_output_files(output_folder, file_patterns):
    """Check that expected output files exist."""
    import glob
    missing_files = []

    for pattern in file_patterns:
        matches = glob.glob(os.path.join(output_folder, pattern))
        if not matches:
            missing_files.append(pattern)
            log.error("File pattern not found: %s" % pattern)
        else:
            log.info("File exists (OK): %s" % matches[0])

    return missing_files


def main():
    args = parse_args()
    set_logger(args, log)

    if not args.config_file:
        log.error("Provide configuration file")
        return -2

    config_file = args.config_file
    if not os.path.exists(config_file):
        log.error("Config file does not exist: %s" % config_file)
        return -3

    log.info("Loading config from %s" % config_file)
    config_dict = load_tsv_config(config_file)

    # Check required fields
    required = ["name", "mode"]
    for k in required:
        if k not in config_dict:
            log.error("%s is not set in the config" % k)
            return -4

    # Must have output in config or command line
    if "output" not in config_dict and not args.output:
        log.error("Output directory must be specified in config or command line")
        return -4

    # Run detect_barcodes
    err_code = run_detect_barcodes(args, config_dict)
    if err_code != 0:
        return err_code

    # Determine run type
    run_type = config_dict.get("run_type", RT_SINGLE)
    err_codes = []

    if run_type == RT_VOID:
        log.info("Run type is void, skipping quality assessment")
        err_codes.append(0)
    elif run_type in (RT_SINGLE, RT_STEREO_SPLIT):
        err_codes.append(run_barcode_quality(args, config_dict))
    else:
        log.error("Unknown run type: %s" % run_type)
        return -5

    # Check expected output files
    if "check_files" in config_dict:
        output_folder = args.output if args.output else config_dict["output"]
        file_patterns = config_dict["check_files"].split()
        missing = check_output_files(output_folder, file_patterns)
        if missing:
            log.error("Missing output files: %s" % " ".join(missing))
            err_codes.append(-6)

    if any(ec != 0 for ec in err_codes):
        return -7

    log.info("All checks passed!")
    return 0


if __name__ == "__main__":
    try:
        ecode = main()
        sys.exit(ecode)
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(1)
