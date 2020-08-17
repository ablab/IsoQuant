# simulator.py transcriptome -rt Mus_musculus.GRCm38.cdna.all.fa -rg Mus_musculus.GRCm38.75.dna.fa -c mouse_cdna/training -e mouse_cdna_chr18/expression_abundance_chr18.tsv -n 100
# python3 IsoQuant/isoquant.py --output isoqunat_errors --bam isoqunat_errors/00_Mouse.ONT.Spatial.NanoSim.EvenCoverage.chr19/00_Mouse.ONT.Spatial.NanoSim.EvenCoverage.chr19.bam --genedb IQ_data/chr19/chr19.gtf --data_type nanopore --complete_genedb --reference IQ_data/chr19/Mouse.chr19.fasta

# Usage: python3 src/quantifiaction_calibration.py -rt ../mouse_rna/Mus_musculus.GRCm38.cdna.all.fa -rg ../mouse_rna/Mus_musculus.GRCm38.75.dna.fa -c ../mouse_rna/mouse_cdna/training -e ../mouse_rna/mouse_cdna_chr18/expression_abundance_chr18.tsv -n 100
# isoseqsim usage: python3 IsoQuant/src/quantifiaction_calibration.py -o examples --ref_g IQ_data/chr19/Mouse.chr19.fasta isoseqsim --gff IQ_data/chr19/chr19.gtf

import subprocess
import argparse
import pathlib

from collections import Counter

import pandas as pd
import numpy as np


def run_quantification(config):
    if config.sim:
        subprocess.call(config.isoseqsim_command, stderr=config.log_file)
        print('Simualation finished')
    subprocess.run(config.isoquant_command, stderr=config.log_file)
    print('Isoquant_finished')
    assert True


class QuantificationConfig:
    def __init__(self, args):
        # stages
        self.sim = args.sim
        self.num_threads = args.num_threads
        self.mode = args.mode
        self.data_type = 'nanopore' if self.mode == 'nanosim' else 'pacbio_raw'

        # output files
        self.sim_output = args.output + '/sim_res'
        self.iso_output = args.output + '/iso_res'

        if self.mode == 'nanosim':
            # nanosim params
            self._init_nanosim_args(args)
        else:
            # isoseqsim params
            self._init_isoseqsim_args(args)

        # isoquant params
        self.isoquant_path = str(pathlib.Path(__file__).parents[1].absolute() / 'isoquant.py')
        self.trainscript_counts = f'{self.iso_output}/00_simulated_reads_normal/00_simulated_reads_normal.transcript_counts.tsv'

        # log params
        self.log_fpath = 'log_quant.txt'
        self.log_file = open(self.log_fpath, "a")

    def _init_nanosim_args(self, args):
        self.ref_genome_fa = args.ref_g or '../mouse_cdna_chr18/mus_musculus.dna.chr18.fa'
        self.ref_transcriptome_fa = args.ref_t or '../mouse_cdna_chr18/mus_musculus.GRCm38.cdna.chr18.fa'
        self.expr_abundance = args.exp or '../mouse_cdna_chr18/expression_abundance_chr18.tsv'
        self.training = args.model_prefix or '../mouse_cdna_chr18/training'
        self.n = args.number
        self.gff = self.training + "_added_intron_final.gff3"
        self.simulated_reads = self.sim_output + '_aligned_reads.fasta'
        self.simulated_error_profile = self.sim_output + '_error_profile.fasta'
        self.simulated_unaligned_reads = self.sim_output + '_unaligned_reads.fasta'

    def _init_isoseqsim_args(self, args):
        self.ref_genome_fa = args.ref_g
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
        return f'{self.isoseq_path} -g {self.ref_genome_fa} '\
               f'-a {self.gff} '\
               f'--es {self.es} --ei {self.ei} --ed {self.ed} --nbn {self.nbn} ' \
               f'--c5 {self.isoseq_utilities / "5_end_completeness.PacBio-P6-C4.tab"} ' \
               f'--c3 {self.isoseq_utilities / "3_end_completeness.PacBio-P6-C4.tab"} ' \
               f'-o {self.simulated_reads} ' \
               f'-t {self.sim_output}/simulated_transcipt_normal.gpd ' \
               f'--tempdir {self.sim_output}/temp_normal'.split()

    @property
    def nanosim_command(self):
        # f'--fastq ' ?
        return f'simulator.py transcriptome ' \
               f'-rt {self.ref_transcriptome_fa} ' \
               f'-rg {self.ref_genome_fa} ' \
               f'-c {self.training} ' \
               f'-e {self.expr_abundance} ' \
               f'-o {self.sim_output} ' \
               f'-n {self.n} ' \
               f'-t {self.num_threads}'.split()

    @property
    def isoquant_command(self):
        # f'--complete_genedb ' \
        return f'python3 {self.isoquant_path} ' \
               f'--output {self.iso_output} ' \
               f'--fastq  {self.simulated_reads} ' \
               f'--genedb {self.gff} ' \
               f'--data_type {self.data_type} ' \
               f'--reference {self.ref_genome_fa} ' \
               f'-t {self.num_threads}'.split()

    def update_abundance(self):
        df = pd.read_csv(self.expr_abundance, sep='\t')
        n_records = len(df['est_counts'])
        df['est_counts'] = df['est_counts'] + pd.Series(np.random.random((n_records,))) * df['est_counts'].mean()
        df['tpm'] = df['est_counts'] * 1000000 / df['est_counts'].sum()
        return df

    @staticmethod
    def compare_quantification(old_abundance, new_abundance):
        # TODO: more advanced check
        return old_abundance == new_abundance


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sim', help='Run simulation stage', action='store_true')
    parser.add_argument('-o', '--output', help='Output location and prefix for simulated reads (Default = simulated)',
                        default="simulated")
    parser.add_argument('-rg', '--ref_g', help='Input reference genome', required=True)
    parser.add_argument('-t', '--num_threads', help='Number of threads for simulation (Default = 1)', type=int,
                        default=1)

    subparsers = parser.add_subparsers(help="You may run the quantification calibration in nanosim or isoseqsim mode.",
                                       dest='mode')

    parser_n = subparsers.add_parser('nanosim', help="Run nanosim for simulation stage")
    parser_n.add_argument('-rt', '--ref_t', help='Input reference transcriptome', default=None)
    parser_n.add_argument('-rg', '--ref_g', help='Input reference genome, required if intron retention simulatin is on',
                          default='')
    parser_n.add_argument('-e', '--exp', help='Expression profile in the specified format as described in README',
                          required=True)
    parser_n.add_argument('-c', '--model_prefix', help='Location and prefix of error profiles generated from '
                                                       'characterization step (Default = training)',
                          default="training")
    parser_n.add_argument('-n', '--number', help='Number of reads to be simulated (Default = 1000)', type=int,
                          default=1000)

    parser_i = subparsers.add_parser('isoseqsim', help="Run isoseq for simulation stage")
    parser_i.add_argument('-a', '--gff', help='Input gtf/gff')
    parser_i.add_argument('--isoseq-path', '--isp', help='Path to IseSeqSim', default='/home/asmetanin/tools/IsoSeqSim/')
    parser_i.add_argument('--es', type=str, default='0.0043', help="Error rate for substitution.")
    parser_i.add_argument('--ei', type=str, default='0.0084', help="Error rate for insertion.")
    parser_i.add_argument('--ed', type=str, default='0.0027', help="Error rate for deletion.")
    parser_i.add_argument('--nbn', type=str, default='100',
                          help="Average read count per transcript to simulate (i.e., the parameter 'n' of the Negative Binomial distribution)")

    return parser.parse_args()


def print_args(config):
    print(config.isoseqsim_command)
    print(config.isoquant_command)
    print(config.trainscript_counts)
    print(config.simulated_reads)


def main():
    args = parse_args()
    if args.mode == 'nanosim':
        print('Not implemented yet.')
        return
    config = QuantificationConfig(args)
    print_args(config)
    run_quantification(config)

    compare_quant_isoseq(config.trainscript_counts, config.simulated_reads)
    print('----well done----')


def get_simulated_isoforms(sim_reads_fpath):
    with open(sim_reads_fpath, 'r') as f_in:
        for line in f_in.readlines():
            if line.startswith('>'):
                _, isoform = line.split()
                yield isoform


def compare_quant_isoseq(isoquant_res_fpath, sim_reads_fpath):
    c = Counter(get_simulated_isoforms(sim_reads_fpath))
    df = pd.read_csv(isoquant_res_fpath, sep='\t', index_col=0)
    full_matches = 0
    close_matches = 0
    detected_isoforms = set()
    for isoform, row in df.iterrows():
        detected_isoforms.add(isoform)
        n_simulated = c[isoform]
        n_inferred = row['count']
        if n_simulated == n_inferred:
            full_matches += 1
        if n_simulated * 0.9 <= n_inferred <= n_simulated * 1.1:
            close_matches += 1
        else:
            print(isoform, n_simulated, n_inferred)



    print('Corrcoef: ', np.corrcoef(a, b))
    print('Full match fraction:', full_matches / len(c))
    print('Close match fraction:', close_matches / len(c))
    print('Not detected: ', len(c) - len(df.index.values), len(c) - len(df.index.values) / len(c), len(c))


def compare_quant(isoquant_res_fpath, sim_reads_fpath):
    c = Counter(get_simulated_isoforms(sim_reads_fpath))
    df = pd.read_csv(isoquant_res_fpath, sep='\t', index_col=0)
    df['sim'] = 0
    full_matches = 0
    close_matches = 0
    for isoform, count in c.items():
        if isoform in df.index:
            df.loc[isoform, 'sim'] = count
            df = df.append({'count': 0, 'sim': count})

    print('Corrcoef: ', np.corrcoef(df['count'], df['sim']))
    print('Full match fraction:', full_matches / len(c))
    print('Close match fraction:', close_matches / len(c))
    print('Not detected: ', len(c) - len(df.index.values), len(c) - len(df.index.values) / len(c), len(c))


if __name__ == '__main__':
    main()
