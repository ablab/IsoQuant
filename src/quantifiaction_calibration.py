# simulator.py transcriptome -rt Mus_musculus.GRCm38.cdna.all.fa -rg Mus_musculus.GRCm38.75.dna.fa -c mouse_cdna/training -e mouse_cdna_chr18/expression_abundance_chr18.tsv -n 100
# python3 IsoQuant/isoquant.py --output isoqunat_errors --bam isoqunat_errors/00_Mouse.ONT.Spatial.NanoSim.EvenCoverage.chr19/00_Mouse.ONT.Spatial.NanoSim.EvenCoverage.chr19.bam --genedb IQ_data/chr19/chr19.gtf --data_type nanopore --complete_genedb --reference IQ_data/chr19/Mouse.chr19.fasta

# Usage: python3 src/quantifiaction_calibration.py -rt ../mouse_rna/Mus_musculus.GRCm38.cdna.all.fa -rg ../mouse_rna/Mus_musculus.GRCm38.75.dna.fa -c ../mouse_rna/mouse_cdna/training -e ../mouse_rna/mouse_cdna_chr18/expression_abundance_chr18.tsv -n 100

import subprocess
import argparse

import pandas as pd
import numpy as np


def test_quantification(config):
    subprocess.call(config.simulate_command, stderr=config.log_file)
    print('Simualation finished')
    subprocess.run(config.isoquant_command, stderr=config.log_file)
    print('Isoquant_finished')
    assert True


class QuantificationConfig:
    def __init__(self, args):
        # nanosim params
        self.ref_genome_fa = args.ref_g or '../mouse_cdna_chr18/mus_musculus.dna.chr18.fa'
        self.ref_transcriptome_fa = args.ref_t or '../mouse_cdna_chr18/mus_musculus.GRCm38.cdna.chr18.fa'
        self.expr_abundance = args.exp or '../mouse_cdna_chr18/expression_abundance_chr18.tsv'
        self.training = args.model_prefix or '../mouse_cdna_chr18/training'
        self.n = args.number
        self.num_threads = args.num_threads

        # log params
        self.log_fpath = 'log_quant.txt'
        self.log_file = open(self.log_fpath, "a")

        # output files
        self.sim_output = 'sim_res/' + args.output
        self.iso_output = 'iso_res'
        self.simulated_aligned_reads = self.sim_output + '_aligned_reads.fasta'
        self.simulated_error_profile = self.sim_output + '_error_profile.fasta'
        self.simulated_unaligned_reads = self.sim_output + '_unaligned_reads.fasta'

    @property
    def simulate_command(self):
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
        return f'python3 ../isoquant.py ' \
               f'--output {self.iso_output} ' \
               f'--fastq  {self.simulated_aligned_reads}' \
               f'--genedb {self.training + "_added_intron_final.gff3"} ' \
               f'--data_type nanopore ' \
               f'--complete_genedb ' \
               f'--reference {self.ref_genome_fa}'.split()

    def update_abundance(self):
        df = pd.read_csv(self.expr_abundance, sep='\t')
        n_records = len(df['est_counts'])
        df['est_counts'] = df['est_counts'] + pd.Series(np.random.random((n_records,))) * df['est_counts'].mean()
        df['tpm'] = df['est_counts'] * 1000000 / df['est_counts'].sum()
        return df

    def compare_quantification(self, old_abundance, new_abundance):
        # TODO: more advanced check
        return old_abundance == new_abundance


def parse_args():
    parser_t = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_t.add_argument('-rt', '--ref_t', help='Input reference transcriptome', required=True)
    parser_t.add_argument('-rg', '--ref_g', help='Input reference genome, required if intron retention simulatin is on',
                          default='')
    parser_t.add_argument('-e', '--exp', help='Expression profile in the specified format as described in README',
                          required=True)
    parser_t.add_argument('-c', '--model_prefix', help='Location and prefix of error profiles generated from '
                                                       'characterization step (Default = training)',
                          default="training")
    parser_t.add_argument('-o', '--output', help='Output location and prefix for simulated reads (Default = simulated)',
                          default="simulated")
    parser_t.add_argument('-n', '--number', help='Number of reads to be simulated (Default = 20000)', type=int,
                          default=20000)
    parser_t.add_argument('-t', '--num_threads', help='Number of threads for simulation (Default = 1)', type=int,
                          default=1)

    # parser_t.add_argument('-max', '--max_len', help='The maximum length for simulated reads (Default = Infinity)',
    #                       type=int, default=float("inf"))
    # parser_t.add_argument('-min', '--min_len', help='The minimum length for simulated reads (Default = 50)',
    #                       type=int, default=50)
    # parser_t.add_argument('--seed', help='Manually seeds the pseudo-random number generator', type=int, default=None)
    # parser_t.add_argument('-k', '--KmerBias', help='Minimum homopolymer length to simulate homopolymer contraction and '
    #                                                'expansion events in',
    #                       type=int, default=None)
    # parser_t.add_argument('-b', '--basecaller', help='Simulate homopolymers and/or base qualities with respect to '
    #                                                  'chosen basecaller: albacore or guppy',
    #                       choices=["albacore", "guppy"], default=None)
    # parser_t.add_argument('-r', '--read_type', help='Simulate homopolymers and/or base qualities with respect to '
    #                                                 'chosen read type: dRNA, cDNA_1D or cDNA_1D2',
    #                       choices=["dRNA", "cDNA_1D", "cDNA_1D2"], default=None)
    # parser_t.add_argument('-s', '--strandness', help='Proportion of sense sequences. Overrides the value '
    #                                                  'profiled in characterization stage. Should be between 0 and 1',
    #                       type=float, default=None)
    # parser_t.add_argument('--no_model_ir', help='Ignore simulating intron retention events', action='store_false',
    #                       default=True)
    # parser_t.add_argument('--perfect', help='Ignore profiles and simulate perfect reads', action='store_true',
    #                       default=False)
    # parser_t.add_argument('--polya', help='Simulate polyA tails for given list of transcripts', default=None)
    # parser_t.add_argument('--fastq', help='Output fastq files instead of fasta files', action='store_true',
    #                       default=False)
    # parser_t.add_argument('--uracil', help='Converts the thymine (T) bases to uracil (U) in the output fasta format',
    #                       action='store_true', default=False)

    return parser_t.parse_args()


def main():
    args = parse_args()
    config = QuantificationConfig(args)
    test_quantification(config)
    print('----well done----')


if __name__ == '__main__':
    main()
