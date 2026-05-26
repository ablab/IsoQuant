
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import subprocess

import pytest
import os
import shutil


def test_run_without_parameters():
    result = subprocess.run(["python3", "isoquant.py"], capture_output=True)
    assert result.returncode == 0
    assert b"usage" in result.stdout


@pytest.mark.parametrize("option", ["-h", "--help", "--full_help"])
def test__help(option):
    result = subprocess.run(["./isoquant.py", option], capture_output=True)
    print(result.returncode)
    assert result.returncode == 0
    assert b"usage" in result.stdout
    assert b"options:" in result.stdout


def test_clean_start():
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, 'simple_data/')
    out_dir = os.path.join(source_dir, "out_full/")
    shutil.rmtree(out_dir, ignore_errors=True)
    sample_name = "ONT_Simulated.chr9.4M"

    result = subprocess.run(["python3", "isoquant.py",
                             "--clean_start",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--fastq", data_dir + "chr9.4M.ont.sim.fq.gz",
                             "--genedb", data_dir + "chr9.4M.gtf.gz", "--complete_genedb",
                             "-r",  data_dir + "chr9.4M.fa.gz",
                             "-t", "2",
                             "--prefix", sample_name,
                             "--count_exons", "--sqanti_output",
                             "--read_group","file:" + data_dir + "chr9.4M.ont.sim.read_groups.tsv" + ":0:1"])

    assert result.returncode == 0
    sample_folder = os.path.join(out_dir, sample_name)
    assert os.path.isdir(sample_folder)
    # Note: corrected_reads.bed.gz and transcript_model_reads.tsv.gz are only generated
    # with --large_output corrected_bed read2transcripts (not in default)
    resulting_files = ["exon_counts.tsv", "exon_grouped_file0_col1_counts.linear.tsv",
                       "exon_grouped_file0_col1_counts.tsv",
                       "gene_counts.tsv", "gene_grouped_file0_col1_counts.tsv",
                       "intron_counts.tsv", "intron_grouped_file0_col1_counts.linear.tsv",
                       "intron_grouped_file0_col1_counts.tsv",
                       "read_assignments.tsv.gz",
                       "novel_vs_known.SQANTI-like.tsv",
                       "transcript_counts.tsv", "transcript_grouped_file0_col1_counts.tsv",
                       "discovered_transcript_counts.tsv", "transcript_models.gtf",
                       "discovered_transcript_tpm.tsv"]
    for f in resulting_files:
        assert os.path.exists(os.path.join(sample_folder, sample_name + "." + f))


def test_usual_start():
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, 'simple_data/')
    out_dir = os.path.join(source_dir, "out_usual/")
    shutil.rmtree(out_dir, ignore_errors=True)
    sample_name = "ONT_Simulated.chr9.4M"
    result = subprocess.run(["python3", "--version"])
    result = subprocess.run(["python3", "isoquant.py",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--fastq", data_dir + "chr9.4M.ont.sim.fq.gz",
                             "--genedb", data_dir + "chr9.4M.gtf.gz", "--complete_genedb",
                             "-r",  data_dir + "chr9.4M.fa.gz",
                             "-t", "2",
                             "--prefix", sample_name])

    assert result.returncode == 0
    sample_folder = os.path.join(out_dir, sample_name)
    assert os.path.isdir(sample_folder)
    resulting_files = ["gene_counts.tsv", "read_assignments.tsv.gz", "transcript_counts.tsv",
                       "discovered_transcript_counts.tsv", "transcript_models.gtf"]
    for f in resulting_files:
        assert os.path.exists(os.path.join(sample_folder, sample_name + "." + f))


def test_with_bam_and_polya():
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, 'simple_data/')
    out_dir = os.path.join(source_dir, "out_polya/")
    shutil.rmtree(out_dir, ignore_errors=True)
    sample_name = "ONT_Simulated.chr9.4M.polyA"

    result = subprocess.run(["python3", "isoquant.py",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--bam", os.path.join(data_dir, "chr9.4M.ont.sim.polya.bam"),
                             "--genedb", os.path.join(data_dir, "chr9.4M.gtf.gz"), "--complete_genedb",
                             "-r",  os.path.join(data_dir, "chr9.4M.fa.gz"),
                             "-t", "2",
                             "--prefix", sample_name,
                             "--sqanti_output", "--count_exons"])

    assert result.returncode == 0
    sample_folder = os.path.join(out_dir, sample_name)
    assert os.path.isdir(sample_folder)
    resulting_files = ["exon_counts.tsv", "gene_counts.tsv",
                       "intron_counts.tsv", "read_assignments.tsv.gz",
                       "novel_vs_known.SQANTI-like.tsv",
                       "transcript_counts.tsv",
                       "discovered_transcript_counts.tsv", "transcript_models.gtf",
                       "polyA_prediction.tsv"]
    for f in resulting_files:
        assert os.path.exists(os.path.join(sample_folder, sample_name + "." + f))
    # polyA prediction should produce at least a header row and some data.
    polya_tsv = os.path.join(sample_folder, sample_name + ".polyA_prediction.tsv")
    with open(polya_tsv) as f:
        lines = f.read().splitlines()
    assert lines[0].split("\t")[:3] == ["chromosome", "transcript_id", "gene_id"]
    assert any(ln and not ln.startswith("chromosome") for ln in lines), \
        "polyA prediction TSV is empty"


def test_with_bam_polya_and_fl_data():
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, 'simple_data/')
    out_dir = os.path.join(source_dir, "out_polya_fl/")
    shutil.rmtree(out_dir, ignore_errors=True)
    sample_name = "ONT_Simulated.chr9.4M.polyA.fl"

    result = subprocess.run(["python3", "isoquant.py",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--bam", os.path.join(data_dir, "chr9.4M.ont.sim.polya.bam"),
                             "--genedb", os.path.join(data_dir, "chr9.4M.gtf.gz"),
                             "--complete_genedb",
                             "-r", os.path.join(data_dir, "chr9.4M.fa.gz"),
                             "-t", "2",
                             "--fl_data",
                             "--prefix", sample_name])

    assert result.returncode == 0
    sample_folder = os.path.join(out_dir, sample_name)
    assert os.path.isdir(sample_folder)
    for f in ["polyA_prediction.tsv", "TSS_prediction.tsv"]:
        path = os.path.join(sample_folder, sample_name + "." + f)
        assert os.path.exists(path), f"expected {f} to be produced"
        with open(path) as fh:
            lines = fh.read().splitlines()
        assert lines[0].split("\t")[:3] == ["chromosome", "transcript_id", "gene_id"]


def test_with_illumina():
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, 'simple_data/')
    out_dir = os.path.join(source_dir, "out_illumina/")
    shutil.rmtree(out_dir, ignore_errors=True)
    sample_name = "ONT_Simulated.chr9.4M.polyA"

    result = subprocess.run(["python3", "isoquant.py",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--bam", os.path.join(data_dir, "chr9.4M.ont.sim.polya.bam"),
                             "--illumina_bam", os.path.join(data_dir, "chr9.4M.Illumina.bam"),
                             "--genedb", os.path.join(data_dir, "chr9.4M.gtf.gz"), "--complete_genedb",
                             "-r",  os.path.join(data_dir, "chr9.4M.fa.gz"),
                             "-t", "2",
                             "--prefix", sample_name])

    assert result.returncode == 0
    sample_folder = os.path.join(out_dir, sample_name)
    assert os.path.isdir(sample_folder)
    resulting_files = ["gene_counts.tsv",
                       "read_assignments.tsv.gz",
                       "transcript_counts.tsv",
                       "discovered_transcript_counts.tsv", "transcript_models.gtf"]
    for f in resulting_files:
        assert os.path.exists(os.path.join(sample_folder, sample_name + "." + f))

def test_with_yaml():
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, 'simple_data/')
    out_dir = os.path.join(source_dir, "out_yaml/")
    shutil.rmtree(out_dir, ignore_errors=True)
    sample_name = "ONT_Simulated.chr9.4M.polyA"

    result = subprocess.run(["python3", "isoquant.py",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--yaml", os.path.join(data_dir, "chr9.4M.yaml"),
                             "--genedb", os.path.join(data_dir, "chr9.4M.gtf.gz"), "--complete_genedb",
                             "-r",  os.path.join(data_dir, "chr9.4M.fa.gz"),
                             "-t", "2"])

    assert result.returncode == 0
    sample_folder = os.path.join(out_dir, sample_name)
    assert os.path.isdir(sample_folder)
    resulting_files = ["gene_counts.tsv",
                       "read_assignments.tsv.gz",
                       "transcript_counts.tsv",
                       "discovered_transcript_counts.tsv", "transcript_models.gtf"]
    for f in resulting_files:
        assert os.path.exists(os.path.join(sample_folder, sample_name + "." + f))


#def test_cage():
#    source_dir = os.path.dirname(os.path.realpath(__file__))
#    data_dir = os.path.join(source_dir, 'toy_data/')
#    out_dir = os.path.join(source_dir, "out_mapt/")
#    shutil.rmtree(out_dir, ignore_errors=True)
#    os.environ['HOME'] = source_dir # dirty hack to set $HOME for tox environment
#    sample_name = "MAPT.Mouse.ONT"
#
#    result = subprocess.run(["python", "isoquant.py", '--clean_start',
#                             "-o", out_dir,
#                             "--data_type", "nanopore",
#                             "--fastq", data_dir + "MAPT.Mouse.ONT.simulated.fastq",
#                             "--genedb", data_dir + "MAPT.Mouse.genedb.gtf", "--complete_genedb",
#                             "-r",  data_dir + "MAPT.Mouse.reference.fasta",
#                             '--cage', data_dir + "MAPT.Mouse.CAGE.bed",
#                             "-t", "2",
#                             "--prefix", sample_name])
#
    # assert result.returncode == 0
    # sample_folder = os.path.join(out_dir, sample_name)
    # assert os.path.isdir(sample_folder)
    # resulting_files = ["gene_counts.tsv", "read_assignments.tsv.gz", "transcript_counts.tsv",
    #                    "transcript_model_counts.tsv", "transcript_models.gtf", "transcript_model_reads.tsv.gz"]
    # for f in resulting_files:
    #     assert os.path.exists(os.path.join(sample_folder, sample_name + "." + f))

