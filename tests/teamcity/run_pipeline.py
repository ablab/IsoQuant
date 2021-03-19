############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# takes teamcity config as input and check isoquant results

import os
import sys
import argparse
import glob
from traceback import print_exc
import subprocess


_quote = {"'": "|'", "|": "||", "\n": "|n", "\r": "|r", '[': '|[', ']': '|]'}


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

    return os.path.join(os.path.dirname(config_file), path)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="Output file name")
    parser.add_argument("config_file", metavar="config_file", type=str, help="configuration .info file")

    args = parser.parse_args()
    return args


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

    source_dir = os.path.dirname(os.path.realpath(__file__))
    isoquant_dir = os.path.join(source_dir, "../../")

    log.log("Loading config from %s" % config_file)
    config_dict = load_tsv_config(config_file)
    for k in ["genedb", "reads", "datatype", "output", "name"]:
        if k not in config_dict:
            log.err(k + "is not set in the config")
            return -10

    label = config_dict["name"]
    output_folder = os.path.join(args.output if args.output else config_dict["output"], label)
    genedb = fix_path(config_file, config_dict["genedb"])
    reads = fix_path(config_file, config_dict["reads"])

    log.start_block('isoquant', 'Running IsoQuant')
    isoquant_command_list = ["python", os.path.join(isoquant_dir, "isoquant.py"), "-o", output_folder,
                    "--genedb", genedb, "-d", config_dict["datatype"], "-t", "16", "-l", label]
    if "bam" in config_dict:
        isoquant_command_list.append("--bam")
        bam = fix_path(config_file, config_dict["bam"])
        isoquant_command_list.append(bam)
        if "genome" in config_dict:
            isoquant_command_list.append("-r")
            isoquant_command_list.append(fix_path(config_file, config_dict["genome"]))
    else:
        if "genome" not in config_dict:
            log.err("genome is not set in the config")
            return -10
        isoquant_command_list.append("--fastq")
        isoquant_command_list.append(reads)
        isoquant_command_list.append("-r")
        isoquant_command_list.append(fix_path(config_file, config_dict["genome"]))
        bam = glob.glob(os.path.join(output_folder, "%s/aux/%s*.bam" % (label, label)))
        if not bam:
            log.err("BAM file was not found")
            return -11
        bam = bam[0]

    if "isoquant_options" in config_dict:
        isoquant_command_list.append(config_dict["isoquant_options"].replace('"', ''))

    log.log("IsoQuant command line: " + " ".join(isoquant_command_list))
    result = subprocess.run(isoquant_command_list)
    if result.returncode != 0:
        log.err("IsoQuant exited with non-zero status: %d" % result.returncode)
        return -11

    output_tsv = os.path.join(output_folder, "%s/%s.read_assignments.tsv" % (label, label))
    log.end_block('isoquant')

    log.start_block('quality', 'Running quality assessment')
    quality_report = os.path.join(output_folder, "report.tsv")
    qa_command_list = ["python", os.path.join(isoquant_dir, "src/assess_assignment_quality.py"),
                       "-o", quality_report, "--gene_db", genedb, "--tsv", output_tsv,
                       "--mapping", bam, "--fasta", reads]

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
    for k, v in etalon_qaulity_dict.items():
        if k not in quality_report_dict:
            log.err("Metric %s was not found in the report" % k)
            exit_code = -22
            continue
        lower_bound = float(v) * 0.99
        upper_bound = float(v) * 1.01
        value = float(quality_report_dict[k])
        if value < lower_bound:
            log.err("Value of %s = %2.2f is lower than the expected value %2.2f" % (k, value, lower_bound))
            exit_code = -20
        else:
            log.log("Value of %s = %2.2f >= %2.2f as expected" % (k, value, lower_bound))
        if value > upper_bound:
            log.err("Value of %s = %2.2f is higher than the expected value %2.2f" % (k, value, upper_bound))
            exit_code = -21
        else:
            log.log("Value of %s = %2.2f <= %2.2f as expected" % (k, value, upper_bound))
    log.end_block('assessment')

    return exit_code


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
