# Usage: python3 misc/isoseq_quantification.py -o examples --reference tests/toy_data/MAPT.Mouse.reference.fasta --gff tests/toy_data/MAPT.Mouse.genedb.gtf --complete --isp tools/IsoSeqSim --nbn 100

import os
import subprocess
import argparse
import pathlib

from collections import Counter

import pandas as pd
import numpy as np


def run_quantification(config):
    if not config.simulated:
        print('Starting simulation with IsoSeqSim.')
        subprocess.call(config.isoseqsim_command, stderr=config.log_file)
        print('Simulation finished')
    else:
        print(f'Using previously simulated reads: {config.simulated_reads}')
    subprocess.run(config.isoquant_command, stderr=config.log_file)
    print('Isoquant finished')
    assert True


def print_args(config):
    print(config.isoseqsim_command)
    print(config.isoquant_command)
    print(config.transcript_counts)
    print(config.simulated_reads)


def main():
    args = parse_args()
    config = QuantificationConfig(args)
    # print_args(config)
    run_quantification(config)
    compare_transcript_counts(config.transcript_counts, config.simulated_reads, config.iso_output)
    print('----well done----')


def get_simulated_isoforms(sim_reads_fpath):
    with open(sim_reads_fpath, 'r') as f_in:
        for line in f_in.readlines():
            if line.startswith('>'):
                _, isoform = line.split(' ')
                isoform = isoform.strip()
                yield isoform


def count_stats(df):
    print('Corrcoef: ', round(np.corrcoef([df['count'], df['sim']])[1, 0], 3))
    full_matches = (df['count'] == df['sim']).astype(int).sum()
    n_isoforms = len(df['sim'])
    print('Full matches:', full_matches, 'Fraction:', round(full_matches / n_isoforms, 3))
    close_matches = ((df['count'] <= df['sim'] * 1.1) & (df['sim'] * 0.9 <= df['count'])).astype(int).sum()
    print('Close matches (10% diff):', close_matches, round(close_matches / n_isoforms, 2))
    close_matches = ((df['count'] <= df['sim'] * 1.2) & (df['sim'] * 0.8 <= df['count'])).astype(int).sum()
    print('Close matches (20% diff):', close_matches, round(close_matches / n_isoforms, 2))
    not_detected = (df['count'] == 0).astype(int).sum()
    print('Not detected:', not_detected, round(not_detected / n_isoforms, 2))

    false_detected = (df['sim'] == 0).astype(int).sum()
    print('False detections:', false_detected, round(false_detected / n_isoforms, 4))


def compare_transcript_counts(isoquant_res_fpath, sim_reads_fpath, iso_output):
    c = Counter(get_simulated_isoforms(sim_reads_fpath))
    df = pd.read_csv(isoquant_res_fpath, sep='\t', index_col=0)
    df = df.drop('group_id', axis=1)
    df = df.drop(['__ambiguous', '__no_feature', '__not_aligned'])
    df['sim'] = 0
    for isoform, count in c.items():
        if isoform in df.index:
            df.loc[isoform, 'sim'] = count
        else:
            df.loc[isoform] = [0, count]

    df.to_csv(iso_output + 'final_counts.tsv', sep='\t')
    count_stats(df)


class QuantificationConfig:
    def __init__(self, args):
        # stages
        self.num_threads = args.num_threads
        self.complete = args.complete
        self.data_type = 'pacbio_raw'

        # output files
        self.sim_output = args.output + '/sim_res'
        self.iso_output = args.output + '/iso_res'

        # isoseqsim params
        self._init_isoseqsim_args(args)

        self.simulated = args.simulated
        if self.simulated:
            self.simulated_reads = self.simulated
        self.sim_name = os.path.splitext(os.path.basename(self.simulated_reads))[0]

        # isoquant params
        self.isoquant_path = str(pathlib.Path(__file__).parents[1].absolute() / 'isoquant.py')
        self.transcript_counts = f'{self.iso_output}/00_{self.sim_name}/00_{self.sim_name}.transcript_counts.tsv'
        self.transcript_model_counts = f'{self.iso_output}/00_{self.sim_name}/00_{self.sim_name}.transcript_models_counts.tsv'
        self.clean_start = args.clean_start

        # log params
        self.log_fpath = 'log_quant.txt'
        self.log_file = open(self.log_fpath, "a")

    def _init_isoseqsim_args(self, args):
        self.reference = args.reference
        self.isoseq_path = pathlib.Path(args.isoseq_path) / 'bin' / 'isoseqsim'
        self.isoseq_utilities = pathlib.Path(args.isoseq_path) / 'utilities'
        self.es = args.es
        self.ei = args.ei
        self.ed = args.ed
        self.nbn = args.nbn
        self.gff = args.gff
        self.simulated_reads = self.sim_output + '/simulated_reads_normal.fa'

    @property
    def isoseqsim_command(self):
        return f'{self.isoseq_path} -g {self.reference} '\
               f'-a {self.gff} '\
               f'--es {self.es} --ei {self.ei} --ed {self.ed} --nbn {self.nbn} ' \
               f'--c5 {self.isoseq_utilities / "5_end_completeness.PacBio-P6-C4.tab"} ' \
               f'--c3 {self.isoseq_utilities / "3_end_completeness.PacBio-P6-C4.tab"} ' \
               f'-o {self.simulated_reads} ' \
               f'-t {self.sim_output}/simulated_transcipt_normal.gpd ' \
               f'--tempdir {self.sim_output}/temp_normal'.split()

    @property
    def isoquant_command(self):
        cmnd = f'python3 {self.isoquant_path} ' \
               f'--output {self.iso_output} ' \
               f'--fastq  {self.simulated_reads} ' \
               f'--genedb {self.gff} ' \
               f'--data_type {self.data_type} ' \
               f'--reference {self.reference} ' \
               f'-t {self.num_threads} '
        if self.clean_start:
            cmnd += f'--clean-start'
        if self.complete:
            cmnd += f' --complete_genedb'
        return cmnd.split()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--simulated', help='Use simulated reads', action='store', default=None)
    parser.add_argument('-o', '--output', help='Output location and prefix for simulated reads (Default = simulated)',
                        default="simulated")
    parser.add_argument('--complete', help='Complete gene db', action='store_true')
    parser.add_argument('--clean-start', help='Clean start for IsoQuant', action='store_true')
    parser.add_argument('-rg', '--reference', help='Input reference genome', required=True)
    parser.add_argument('-t', '--num_threads', help='Number of threads for simulation (Default = 1)', type=int,
                        default=1)

    parser.add_argument('-a', '--gff', help='Input gtf/gff')
    parser.add_argument('--isoseq-path', '--isp', help='Path to IseSeqSim', default='')  # default='/home/asmetanin/tools/IsoSeqSim/')
    parser.add_argument('--es', type=str, default='0.0043', help="Error rate for substitution.")
    parser.add_argument('--ei', type=str, default='0.0084', help="Error rate for insertion.")
    parser.add_argument('--ed', type=str, default='0.0027', help="Error rate for deletion.")
    parser.add_argument('--nbn', type=str, default='10',
                        help="Average read count per transcript to simulate (i.e., the parameter 'n' of the Negative Binomial distribution)")

    return parser.parse_args()


if __name__ == '__main__':
    main()
