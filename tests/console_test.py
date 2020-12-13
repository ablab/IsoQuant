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
    # FIXME
    os.environ['HOME'] = source_dir
    result = subprocess.run(["python", "isoquant.py",
                             "--clean_start",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--fastq", data_dir + "chr9.4M.ont.sim.fq.gz",
                             "--genedb", data_dir + "chr9.4M.gtf.gz", "--complete_genedb",
                             "-r",  data_dir + "chr9.4M.fa.gz",
                             "-t", "4",
                             "--labels", "ONT_Simulated.chr9.4M",
                             "--count_exons", "--sqanti_output",
                             "--read_group","file:" + data_dir + "chr9.4M.ont.sim.read_groups.tsv" + ":0:1"])
    # TODO add asserts
    assert result.returncode == 0

def test_usual_start():
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, 'simple_data/')
    out_dir = os.path.join(source_dir, "out_usual/")
    shutil.rmtree(out_dir, ignore_errors=True)
    # FIXME
    os.environ['HOME'] = source_dir
    result = subprocess.run(["python", "isoquant.py",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--fastq", data_dir + "chr9.4M.ont.sim.fq.gz",
                             "--genedb", data_dir + "chr9.4M.gtf.gz", "--complete_genedb",
                             "-r",  data_dir + "chr9.4M.fa.gz",
                             "-t", "4",
                             "--labels", "ONT_Simulated.chr9.4M"])
    # TODO add asserts
    assert result.returncode == 0

def test_c():
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, 'simple_data/')
    out_dir = os.path.join(source_dir, "out_polya/")
    shutil.rmtree(out_dir, ignore_errors=True)
    os.environ['HOME'] = source_dir
    result = subprocess.run(["python", "isoquant.py",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--bam", data_dir + "chr9.4M.ont.sim.polya.bam",
                             "--genedb", data_dir + "chr9.4M.gtf.gz", "--complete_genedb",
                             "-r",  data_dir + "chr9.4M.fa.gz",
                             "-t", "4",
                             "--labels", "ONT_Simulated.chr9.4M.polyA",
                             "--has_polya", "--sqanti_output", "--count_exons",
                             "--read_group","file:" + data_dir + "chr9.4M.ont.sim.read_groups.tsv" + ":0:1"])
    # TODO add asserts
    assert result.returncode == 0

