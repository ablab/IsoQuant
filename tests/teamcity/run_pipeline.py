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


_quote = {"'": "|'", "|": "||", "\n": "|n", "\r": "|r", '[': '|[', ']': '|]'}


RT_VOID = "void"
RT_ASSIGNMENT = "assignment"
RT_TRANSCRIPTS = "transcripts"


def escape_value(value):
    return "".join(_quote.get(x, x) for x in value)


class TeamCityLog:

    text = ""

    def _tc_out(self, token, **attributes):
        message = "##teamcity[%s" % token

        for k in sorted(attributes.keys()):
            value = attributes[k]
            if value is None:
                continue

            message += (" %s='%s'" % (k, escape_value(value)))

        message += "]\n"
        sys.stdout.write(message)
        sys.stdout.flush()

    def start_block(self, name, desc):
        sys.stdout.flush()
        self._tc_out("blockOpened", name = name, description = desc)

    def end_block(self, name):
        sys.stdout.flush()
        self._tc_out("blockClosed", name = name)

    def log(self, s):
        self.text += s + "\n"
        self._tc_out("message", text = s)

    def warn(self, s):
        msg = "WARNING: " + s + "\n"
        self.text += msg
        self._tc_out("message", text = s, status = 'WARNING')

    def err(self, s, context = ""):
        msg = "ERROR: " + s + "\n"
        self.text += msg
        self._tc_out("message", text = s, status = 'ERROR', errorDetails = context)

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

    def record_metric(self, name, value):
        self._tc_out("buildStatisticValue", key=name, value=value)


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


def run_isoquant(args, config_dict, log):
    source_dir = os.path.dirname(os.path.realpath(__file__))
    isoquant_dir = os.path.join(source_dir, "../../")
    config_file = args.config_file

    run_name = config_dict["name"]
    output_folder = os.path.join(args.output if args.output else config_dict["output"], run_name)
    if "resume" in config_dict:
        assert "label" in config_dict
        if not os.path.exists(output_folder):
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

        log.start_block('isoquant', 'Running IsoQuant')
        isoquant_command_list = ["python3", os.path.join(isoquant_dir, "isoquant.py"), "-o", output_folder,
                                 "-r", genome, "-d", config_dict["datatype"], "-l", run_name]
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
        log.log("Appending additional options: %s" % config_dict["isoquant_options"])
        opts = config_dict["isoquant_options"].replace('"', '').split()
        for o in opts:
            isoquant_command_list.append(o)

    log.log("IsoQuant command line: " + " ".join(isoquant_command_list))
    result = subprocess.run(isoquant_command_list)
    if result.returncode != 0:
        log.err("IsoQuant exited with non-zero status: %d" % result.returncode)
        return -11

    log.end_block('isoquant')
    return 0


def check_value(etalon_value, output_value, name, log):
    lower_bound = etalon_value * 0.99
    upper_bound = etalon_value * 1.01
    exit_code = 0
    if output_value < lower_bound:
        log.err("Value of %s = %2.2f is lower than the expected value %2.2f" % (name, output_value, lower_bound))
        exit_code = -41
    else:
        log.log("Value of %s = %2.2f >= %2.2f as expected" % (name, output_value, lower_bound))
    if output_value > upper_bound:
        log.err("Value of %s = %2.2f is higher than the expected value %2.2f" % (name, output_value, upper_bound))
        exit_code = -42
    else:
        log.log("Value of %s = %2.2f <= %2.2f as expected" % (name, output_value, upper_bound))
    return exit_code


def run_assignment_quality(args, config_dict, log):
    log.start_block('quality', 'Running quality assessment')
    config_file = args.config_file
    source_dir = os.path.dirname(os.path.realpath(__file__))
    isoquant_dir = os.path.join(source_dir, "../../")

    label = config_dict["label"]
    output_folder = os.path.join(args.output if args.output else config_dict["output"], label)
    output_tsv = os.path.join(output_folder, "%s/%s.read_assignments.tsv" % (label, label))

    assert "genedb" in config_dict
    genedb = fix_path(config_file, config_dict["genedb"])
    reads = fix_path(config_file, config_dict["reads"])
    if "bam" not in config_dict:
        bam = glob.glob(os.path.join(output_folder, "%s/aux/%s*.bam" % (label, label)))
        if not bam:
            log.err("BAM file was not found")
            return -21
        bam = bam[0]
    else:
        bam = fix_path(config_file, config_dict["bam"])

    quality_report = os.path.join(output_folder, "report.tsv")
    qa_command_list = ["python3", os.path.join(isoquant_dir, "misc/assess_assignment_quality.py"),
                       "-o", quality_report, "--gene_db", genedb, "--tsv", output_tsv,
                       "--mapping", bam, "--fasta", reads]

    if "qa_options" in config_dict:
        log.log("Appending additional options: %s" % config_dict["qa_options"])
        opts = config_dict["qa_options"].replace('"', '').split()
        for o in opts:
            qa_command_list.append(o)

    log.log("QA command line: " + " ".join(qa_command_list))
    result = subprocess.run(qa_command_list)
    if result.returncode != 0:
        log.err("QA exited with non-zero status: %d" % result.returncode)
        return -13
    log.end_block('quality')

    if "etalon" not in config_dict:
        return 0

    log.start_block('assessment', 'Checking quality metrics')
    etalon_qaulity_dict = load_tsv_config(fix_path(config_file, config_dict["etalon"]))
    quality_report_dict = load_tsv_config(quality_report)
    exit_code = 0
    new_etalon_outf = open(os.path.join(output_folder, "new_assignment_etalon.tsv"), "w")
    for k, v in etalon_qaulity_dict.items():
        if k not in quality_report_dict:
            log.err("Metric %s was not found in the report" % k)
            exit_code = -22
            continue
        new_etalon_outf.write("%s\t%.2f\n" % (k, float(quality_report_dict[k])))
        err_code = check_value(float(v), float(quality_report_dict[k]), k, log)
        if err_code != 0:
            exit_code = err_code

    new_etalon_outf.close()
    log.end_block('assessment')
    return exit_code


def parse_gffcomapre(stats_file):
    for l in open(stats_file):
        if l.find("Transcript level") == -1:
            continue
        v = l.split()
        assert len(v) > 4
        return float(v[2]), float(v[4])
    return -1, -1


def run_transcript_quality(args, config_dict, log):
    log.start_block('quality', 'Running quality assessment')
    config_file = args.config_file
    source_dir = os.path.dirname(os.path.realpath(__file__))
    isoquant_dir = os.path.join(source_dir, "../../")

    name = config_dict["name"]
    label = name if "label" not in config_dict else config_dict["label"]
    output_folder = os.path.join(args.output if args.output else config_dict["output"], name)
    out_gtf = os.path.join(output_folder, "%s/%s.transcript_models.gtf" % (label, label))
    if not out_gtf:
        log.err("Output GTF file was not found")
        return -31

    quality_output = os.path.join(output_folder, "gffcompare")
    genedb_prefix = fix_path(config_file, config_dict["reduced_db"])
    qa_command_list = ["python3", os.path.join(isoquant_dir, "misc/reduced_db_gffcompare.py"),
                       "-o", quality_output, "--genedb", genedb_prefix, "--gtf", out_gtf, "--tool", "isoquant"]

    log.log("QA command line: " + " ".join(qa_command_list))
    result = subprocess.run(qa_command_list)
    if result.returncode != 0:
        log.err("Transcript evaluation exited with non-zero status: %d" % result.returncode)
        return -13
    log.end_block('quality')

    if "etalon" not in config_dict:
        return 0

    log.start_block('assessment', 'Checking quality metrics')
    etalon_qaulity_dict = load_tsv_config(fix_path(config_file, config_dict["etalon"]))
    exit_code = 0
    new_etalon_outf = open(os.path.join(quality_output, "new_gtf_etalon.tsv"), "w")
    for gtf_type in ['full', 'known', 'novel']:
        recall, precision = parse_gffcomapre(os.path.join(quality_output, "isoquant." + gtf_type + ".stats"))
        metric_name = gtf_type + "_recall"
        if metric_name in etalon_qaulity_dict:
            new_etalon_outf.write("%s\t%.2f\n" % (metric_name, recall))
            etalon_recall = float(etalon_qaulity_dict[metric_name])
            err_code = check_value(etalon_recall, recall , metric_name, log)
            if err_code != 0:
                exit_code = err_code
        metric_name = gtf_type + "_precision"
        if metric_name in etalon_qaulity_dict:
            new_etalon_outf.write("%s\t%.2f\n" % (metric_name, precision))
            etalon_precision = float(etalon_qaulity_dict[metric_name])
            err_code = check_value(etalon_precision, precision, metric_name, log)
            if err_code != 0:
                exit_code = err_code
    new_etalon_outf.close()
    log.end_block('assessment')
    return exit_code


def main():
    log = TeamCityLog()
    args = parse_args()
    if not args.config_file:
        log.err("Provide configuration file")
        exit(-2)

    config_file = args.config_file
    if not os.path.exists(config_file):
        log.err("Provide correct path to configuration file")
        exit(-3)

    log.log("Loading config from %s" % config_file)
    config_dict = load_tsv_config(config_file)
    required = ["output", "name"]
    if "resume" not in config_dict:
        required += ["genome", "datatype"]
    for k in required:
        if k not in config_dict:
            log.err(k + " is not set in the config")
            return -10

    err_code = run_isoquant(args, config_dict, log)
    if err_code != 0:
        return err_code

    run_type = config_dict["run_type"]
    if run_type == RT_VOID:
        err_code = 0
    elif run_type == RT_ASSIGNMENT:
        err_code = run_assignment_quality(args, config_dict, log)
    elif run_type == RT_TRANSCRIPTS:
        err_code = run_transcript_quality(args, config_dict, log)
    else:
        log.err("Test type %s is not supported" % run_type)
        err_code = -50

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
