import subprocess

import pytest
import os
import shutil


def test_run_without_parameters():
    result = subprocess.run(["python", "isoquant.py"], capture_output=True)
    assert result.returncode == 2
    assert b"usage" in result.stderr
    assert b"error: the following arguments are required: --genedb/-g, --data_type/-d" in result.stderr


@pytest.mark.parametrize("option", ["-h", "--help", "--full_help"])
def test__help(option):
    result = subprocess.run(["python", "isoquant.py", option], capture_output=True)
    assert result.returncode == 0
    assert b"usage" in result.stdout
    assert b"optional arguments:" in result.stdout

def test_clean_start():
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, 'simple_data/')
    out_dir = os.path.join(source_dir, "out_full/")
    shutil.rmtree(out_dir, ignore_errors=True)
    os.environ['HOME'] = source_dir # dirty hack to set $HOME for tox environment
    sample_name = "ONT_Simulated.chr9.4M"

    result = subprocess.run(["python", "isoquant.py",
                             "--clean_start",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--fastq", data_dir + "chr9.4M.ont.sim.fq.gz",
                             "--genedb", data_dir + "chr9.4M.gtf.gz", "--complete_genedb",
                             "-r",  data_dir + "chr9.4M.fa.gz",
                             "-t", "4",
                             "--labels", sample_name,
                             "--count_exons", "--sqanti_output",
                             "--read_group","file:" + data_dir + "chr9.4M.ont.sim.read_groups.tsv" + ":0:1"])

    assert result.returncode == 0
    sample_folder = os.path.join(out_dir, sample_name)
    assert os.path.isdir(sample_folder)
    resulting_files = ["bam", "exon_counts.tsv", "exon_grouped_counts.tsv", "gene_counts.tsv", "gene_grouped_counts.tsv",
                       "intron_counts.tsv", "intron_grouped_counts.tsv", "mapped_reads.bed", "read_assignments.tsv",
                       "SQANTI-like.tsv", "transcript_counts.tsv", "transcript_grouped_counts.tsv",
                       "transcript_models_counts.tsv", "transcript_models.gtf", "transcript_models_reads.tsv"]
    for f in resulting_files:
        assert os.path.exists(os.path.join(sample_folder, sample_name + "." + f))

def test_usual_start():
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, 'simple_data/')
    out_dir = os.path.join(source_dir, "out_usual/")
    shutil.rmtree(out_dir, ignore_errors=True)
    os.environ['HOME'] = source_dir # dirty hack to set $HOME for tox environment
    sample_name = "ONT_Simulated.chr9.4M"

    result = subprocess.run(["python", "isoquant.py",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--fastq", data_dir + "chr9.4M.ont.sim.fq.gz",
                             "--genedb", data_dir + "chr9.4M.gtf.gz", "--complete_genedb",
                             "-r",  data_dir + "chr9.4M.fa.gz",
                             "-t", "4",
                             "--labels", sample_name])

    assert result.returncode == 0
    sample_folder = os.path.join(out_dir, sample_name)
    assert os.path.isdir(sample_folder)
    resulting_files = ["gene_counts.tsv", "read_assignments.tsv", "transcript_counts.tsv",
                       "transcript_models_counts.tsv", "transcript_models.gtf", "transcript_models_reads.tsv"]
    for f in resulting_files:
        assert os.path.exists(os.path.join(sample_folder, sample_name + "." + f))

def test_with_bam_and_polya():
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, 'simple_data/')
    out_dir = os.path.join(source_dir, "out_polya/")
    shutil.rmtree(out_dir, ignore_errors=True)
    os.environ['HOME'] = source_dir
    sample_name = "ONT_Simulated.chr9.4M.polyA"

    result = subprocess.run(["python", "isoquant.py",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--bam", data_dir + "chr9.4M.ont.sim.polya.bam",
                             "--genedb", data_dir + "chr9.4M.gtf.gz", "--complete_genedb",
                             "-r",  data_dir + "chr9.4M.fa.gz",
                             "-t", "4",
                             "--labels", sample_name,
                             "--has_polya", "--sqanti_output", "--count_exons"])

    assert result.returncode == 0
    sample_folder = os.path.join(out_dir, sample_name)
    assert os.path.isdir(sample_folder)
    resulting_files = ["exon_counts.tsv", "gene_counts.tsv",
                       "intron_counts.tsv", "mapped_reads.bed", "read_assignments.tsv",
                       "SQANTI-like.tsv", "transcript_counts.tsv",
                       "transcript_models_counts.tsv", "transcript_models.gtf", "transcript_models_reads.tsv"]
    for f in resulting_files:
        assert os.path.exists(os.path.join(sample_folder, sample_name + "." + f))

