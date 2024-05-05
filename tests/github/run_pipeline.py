#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# takes teamcity config as input and check isoquant results

import os
import shutil
import sys
import argparse
import glob
from traceback import print_exc
import subprocess
import logging


_quote = {"'": "|'", "|": "||", "\n": "|n", "\r": "|r", '[': '|[', ']': '|]'}


RT_VOID = "void"
RT_ASSIGNMENT = "assignment"
RT_TRANSCRIPTS = "transcripts"
RT_QUANTIFICATION_KNOWN = "quantification_known"
RT_QUANTIFICATION_NOVEL = "quantification_novel"
RT_QUANTIFICATION_GENES = "quantification_genes"


log = logging.getLogger('GitHubRunner')


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
    parser.add_argument("--output", "-o", type=str, help="Output file name")
    parser.add_argument("config_file", metavar="config_file", type=str, help="configuration .info file")

    args = parser.parse_args()
    return args


def run_isoquant(args, config_dict):
    source_dir = os.path.dirname(os.path.realpath(__file__))
    isoquant_dir = os.path.join(source_dir, "../../")
    config_file = args.config_file

    run_name = config_dict["name"]
    output_folder = os.path.join(args.output if args.output else config_dict["output"], run_name)
    if "resume" in config_dict:
        assert "label" in config_dict
        if os.path.exists(output_folder):
            shutil.rmtree(output_folder)
        os.makedirs(output_folder)
        isoquant_command_list = ["python3", os.path.join(isoquant_dir, "isoquant.py"), "-o", output_folder, "--resume"]
        src_dir = fix_path(config_file, config_dict["resume"])
        for f in os.listdir(src_dir):
            fpath = os.path.join(src_dir, f)
            if os.path.isdir(fpath):
                shutil.copytree(fpath, os.path.join(output_folder, f), dirs_exist_ok=True)
            else:
                shutil.copy2(fpath, os.path.join(output_folder, f))
    else:
        genedb = fix_path(config_file, config_dict["genedb"]) if "genedb" in config_dict else None
        genome = fix_path(config_file, config_dict["genome"])
        config_dict["label"] = run_name

        log.info('== Running IsoQuant ==')
        isoquant_command_list = ["python3", os.path.join(isoquant_dir, "isoquant.py"), "-o", output_folder,
                                 "-r", genome, "-d", config_dict["datatype"], "-p", run_name]
        if genedb:
            isoquant_command_list += ["--genedb", genedb]
        if "bam" in config_dict:
            isoquant_command_list.append("--bam")
            bam = fix_path(config_file, config_dict["bam"])
            isoquant_command_list.append(bam)
        elif "reads" in config_dict:
            reads = fix_path(config_file, config_dict["reads"])
            isoquant_command_list.append("--fastq")
            isoquant_command_list.append(reads)

    if "isoquant_options" in config_dict:
        log.info("Appending additional options: %s" % config_dict["isoquant_options"])
        opts = config_dict["isoquant_options"].replace('"', '').split()
        for o in opts:
            isoquant_command_list.append(o)

    log.info("IsoQuant command line: " + " ".join(isoquant_command_list))
    result = subprocess.run(isoquant_command_list)
    if result.returncode != 0:
        log.error("IsoQuant exited with non-zero status: %d" % result.returncode)
        return -11

    return 0


def check_value(etalon_value, output_value, name):
    lower_bound = etalon_value * 0.99
    upper_bound = etalon_value * 1.01
    exit_code = 0
    if output_value < 0:
        if output_value != etalon_value:
            log.error("Value of %s = %2.2f is not equal to value %2.2f" % (name, output_value, lower_bound))
            exit_code = -2
        else:
            log.info("Value of %s = %2.2f == %2.2f as expected" % (name, output_value, upper_bound))
    else:
        if output_value < lower_bound:
            log.error("Value of %s = %2.2f is lower than the expected value %2.2f" % (name, output_value, lower_bound))
            exit_code = -2
        else:
            log.info("Value of %s = %2.2f >= %2.2f as expected" % (name, output_value, lower_bound))
        if output_value > upper_bound:
            log.error("Value of %s = %2.2f is higher than the expected value %2.2f" % (name, output_value, upper_bound))
            exit_code = -2
        else:
            log.info("Value of %s = %2.2f <= %2.2f as expected" % (name, output_value, upper_bound))
    return exit_code


def find_bam(output_folder, label):
    bam = glob.glob(os.path.join(output_folder, "%s/aux/%s*.bam" % (label, label)))
    if bam:
        return bam[0]

    for l in open(os.path.join(output_folder, "isoquant.log")):
        if l.find("BAM file: ") != -1:
            return l.strip().split()[-1]

    return None


def run_assignment_quality(args, config_dict):
    log.info('== Running quality assessment ==')
    config_file = args.config_file
    source_dir = os.path.dirname(os.path.realpath(__file__))
    isoquant_dir = os.path.join(source_dir, "../../")

    label = config_dict["label"]
    output_folder = os.path.join(args.output if args.output else config_dict["output"], label)
    output_tsv = os.path.join(output_folder, "%s/%s.read_assignments.tsv" % (label, label))
    if not os.path.exists(output_tsv):
        output_tsv = os.path.join(output_folder, "%s/%s.read_assignments.tsv.gz" % (label, label))

    assert "genedb" in config_dict
    genedb = fix_path(config_file, config_dict["genedb"])
    reads = fix_path(config_file, config_dict["reads"])
    if "bam" not in config_dict:
        bam = find_bam(output_folder, label)
        if not bam or not os.path.exists(bam):
            log.error("BAM file was not found")
            return -10
    else:
        bam = fix_path(config_file, config_dict["bam"])

    quality_report = os.path.join(output_folder, "report.tsv")
    qa_command_list = ["python3", os.path.join(isoquant_dir, "misc/assess_assignment_quality.py"),
                       "-o", quality_report, "--gene_db", genedb, "--tsv", output_tsv,
                       "--mapping", bam, "--fasta", reads]

    if "qa_options" in config_dict:
        log.info("Appending additional options: %s" % config_dict["qa_options"])
        opts = config_dict["qa_options"].replace('"', '').split()
        for o in opts:
            qa_command_list.append(o)

    log.info("QA command line: " + " ".join(qa_command_list))
    result = subprocess.run(qa_command_list)
    if result.returncode != 0:
        log.error("QA exited with non-zero status: %d" % result.returncode)
        return -11

    if "etalon" not in config_dict:
        return 0

    log.info('== Checking quality metrics ==')
    etalon_qaulity_dict = load_tsv_config(fix_path(config_file, config_dict["etalon"]))
    quality_report_dict = load_tsv_config(quality_report)
    exit_code = 0
    new_etalon_outf = open(os.path.join(output_folder, "new_assignment_etalon.tsv"), "w")
    for k, v in etalon_qaulity_dict.items():
        if k not in quality_report_dict:
            log.error("Metric %s was not found in the report" % k)
            exit_code = -12
            continue
        new_etalon_outf.write("%s\t%.2f\n" % (k, float(quality_report_dict[k])))
        err_code = check_value(float(v), float(quality_report_dict[k]), k)
        if err_code != 0:
            exit_code = err_code

    new_etalon_outf.close()
    return exit_code


def parse_gffcomapre(stats_file):
    for l in open(stats_file):
        if l.find("Transcript level") == -1:
            continue
        v = l.split()
        assert len(v) > 4
        return float(v[2]), float(v[4])
    return -1, -1


def run_transcript_quality(args, config_dict):
    log.info('== Running quality assessment ==')
    config_file = args.config_file
    source_dir = os.path.dirname(os.path.realpath(__file__))
    isoquant_dir = os.path.join(source_dir, "../../")

    name = config_dict["name"]
    label = name if "label" not in config_dict else config_dict["label"]
    output_folder = os.path.join(args.output if args.output else config_dict["output"], name)
    out_gtf = os.path.join(output_folder, "%s/%s.transcript_models.gtf" % (label, label))
    if not os.path.exists(out_gtf):
        log.error("Output GTF file was not found" % out_gtf)
        return -20

    quality_output = os.path.join(output_folder, "gffcompare")
    genedb_prefix = fix_path(config_file, config_dict["reduced_db"])
    qa_command_list = ["python3", os.path.join(isoquant_dir, "misc/reduced_db_gffcompare.py"),
                       "-o", quality_output, "--genedb", genedb_prefix, "--gtf", out_gtf, "--tool", "isoquant"]

    log.info("QA command line: " + " ".join(qa_command_list))
    result = subprocess.run(qa_command_list)
    if result.returncode != 0:
        log.error("Transcript evaluation exited with non-zero status: %d" % result.returncode)
        return -21

    if "etalon" not in config_dict:
        return 0

    log.info('== Checking quality metrics ==')
    etalon_qaulity_dict = load_tsv_config(fix_path(config_file, config_dict["etalon"]))
    exit_code = 0
    new_etalon_outf = open(os.path.join(quality_output, "new_gtf_etalon.tsv"), "w")
    for gtf_type in ['full', 'known', 'novel']:
        recall, precision = parse_gffcomapre(os.path.join(quality_output, "isoquant." + gtf_type + ".stats"))
        metric_name = gtf_type + "_recall"
        if metric_name in etalon_qaulity_dict:
            new_etalon_outf.write("%s\t%.2f\n" % (metric_name, recall))
            etalon_recall = float(etalon_qaulity_dict[metric_name])
            err_code = check_value(etalon_recall, recall , metric_name)
            if err_code != 0:
                exit_code = err_code
        metric_name = gtf_type + "_precision"
        if metric_name in etalon_qaulity_dict:
            new_etalon_outf.write("%s\t%.2f\n" % (metric_name, precision))
            etalon_precision = float(etalon_qaulity_dict[metric_name])
            err_code = check_value(etalon_precision, precision, metric_name)
            if err_code != 0:
                exit_code = err_code
    new_etalon_outf.close()
    return exit_code


def run_quantification(args, config_dict, mode):
    log.info('== Running quantification assessment ==')
    config_file = args.config_file
    source_dir = os.path.dirname(os.path.realpath(__file__))
    isoquant_dir = os.path.join(source_dir, "../../")

    name = config_dict["name"]
    label = name if "label" not in config_dict else config_dict["label"]
    output_folder = os.path.join(args.output if args.output else config_dict["output"], name)

    if mode == "novel":
        out_tpm = os.path.join(output_folder, "%s/%s.transcript_model_tpm.tsv" % (label, label))
    elif mode == "ref":
        out_tpm = os.path.join(output_folder, "%s/%s.transcript_tpm.tsv" % (label, label))
    else:
        out_tpm = os.path.join(output_folder, "%s/%s.gene_tpm.tsv" % (label, label))


    if not os.path.exists(out_tpm):
        log.error("Output TPM file %s was not found" % out_tpm)
        return -30

    if "reference_tpm" not in config_dict:
        return 0
    ref_tpm = config_dict["reference_gene_tpm"] if mode == "gene" else config_dict["reference_tpm"]
    if not os.path.exists(ref_tpm):
        log.error("File %s with reference TPM was not detected" % ref_tpm)
        return -31

    quantification_stats_output = os.path.join(output_folder, mode + ".quantification.tsv")
    qa_command_list = ["python3", os.path.join(isoquant_dir, "misc/quantification_stats.py"),
                       "-o", quantification_stats_output, "--ref_expr", ref_tpm, "--tpm", out_tpm]

    if mode == "novel":
        gffcompare_outdir = os.path.join(output_folder, "gffcompare")
        tracking = os.path.join(gffcompare_outdir, "isoquant.novel.tracking")
        if not os.path.exists(tracking):
            log.error("Tracking file %s was not found" % tracking)
            return -33
        qa_command_list += ["--tracking", tracking]

    log.info("QA command line: " + " ".join(qa_command_list))
    result = subprocess.run(qa_command_list)
    if result.returncode != 0:
        log.error("Quantification evaluation exited with non-zero status: %d" % result.returncode)
        return -34

    etalon_to_use = "etalon_quantification_" + mode
    if etalon_to_use not in config_dict:
        return 0

    log.info('== Checking quantification metrics ==')
    ref_value_files = fix_path(config_file, config_dict[etalon_to_use])
    if not os.path.exists(ref_value_files):
        log.error("File %s with etalon metric values was not detected" % ref_value_files)
        return -35
    etalon_quality_dict = load_tsv_config(ref_value_files)
    real_dict = load_tsv_config(quantification_stats_output)
    exit_code = 0

    for metric_name in etalon_quality_dict.keys():
        if metric_name not in real_dict:
            exit_code = -36
            log.warning("Metric not found %s" % metric_name)
            continue

        ref_value = float(etalon_quality_dict[metric_name])
        real_value = float(real_dict[metric_name])
        err_code = check_value(ref_value, real_value, metric_name)
        if err_code != 0:
            exit_code = err_code

    return exit_code


def check_output_files(out_dir, label, file_list):
    missing_files = []
    internal_output_dir = os.path.join(out_dir, label)
    for f in file_list:
        fpath = os.path.join(internal_output_dir, label + "." + f)
        if not os.path.exists(fpath):
            missing_files.append(str(fpath))
            log.error("File is missing: %s" % str(fpath))
        else:
            log.info("File exists (OK): %s" % str(fpath))
    return missing_files


def main():
    args = parse_args()
    set_logger(args, log)
    if not args.config_file:
        log.error("Provide configuration file")
        exit(-2)

    config_file = args.config_file
    if not os.path.exists(config_file):
        log.error("Provide correct path to configuration file, %s does not exits" % config_file)
        exit(-3)

    log.info("Loading config from %s" % config_file)
    config_dict = load_tsv_config(config_file)
    required = ["output", "name"]
    if "resume" not in config_dict:
        required += ["genome", "datatype"]
    for k in required:
        if k not in config_dict:
            log.error(k + " is not set in the config")
            exit(-4)

    err_code = run_isoquant(args, config_dict)
    if err_code != 0:
        return err_code

    run_types = set(map(lambda x: x.strip().replace('"', ''), config_dict["run_type"].split(",")))
    err_codes = []
    if RT_VOID in run_types:
        err_codes.append(0)
    if RT_ASSIGNMENT in run_types:
        err_codes.append(run_assignment_quality(args, config_dict))
    if RT_TRANSCRIPTS in run_types:
        err_codes.append(run_transcript_quality(args, config_dict))
    if RT_QUANTIFICATION_KNOWN in run_types:
        err_codes.append(run_quantification(args, config_dict, "ref"))
    if RT_QUANTIFICATION_NOVEL in run_types:
        err_codes.append(run_quantification(args, config_dict, "novel"))
    if RT_QUANTIFICATION_GENES in run_types:
        err_codes.append(run_quantification(args, config_dict, "gene"))

    if "check_input_files" in config_dict:
        files_list = config_dict["check_input_files"].split()
        label = config_dict["label"]
        run_name = config_dict["name"]
        output_folder = os.path.join(args.output if args.output else config_dict["output"], run_name)
        missing_files = check_output_files(output_folder, label, files_list)
        if missing_files:
            log.error("The following files were not detected in the output folder: %s" % "  ".join(missing_files))
            err_codes.append(-5)

    if any(ec != 0 for ec in err_codes):
        err_code = -6

    return err_code


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
        sys.exit(-1)
