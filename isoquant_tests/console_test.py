
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
                             "--count_exons", "--count_intron_retentions", "--sqanti_output",
                             "--read_group", "file:" + data_dir + "chr9.4M.ont.sim.read_groups.tsv" + ":0:1"])

    assert result.returncode == 0
    sample_folder = os.path.join(out_dir, sample_name)
    assert os.path.isdir(sample_folder)
    # Note: corrected_reads.bed.gz and transcript_model_reads.tsv.gz are only generated
    # with --large_output corrected_bed read2transcripts (not in default)
    resulting_files = ["exon_counts.tsv", "exon_grouped_file0_col1_counts.linear.tsv",
                       "exon_splice_site_counts.tsv",
                       "exon_splice_site_grouped_file0_col1_counts.linear.tsv",
                       "gene_counts.tsv", "gene_grouped_file0_col1_counts.tsv",
                       "splice_junction_counts.tsv", "splice_junction_grouped_file0_col1_counts.linear.tsv",
                       "intron_retention_counts.tsv", "intron_retention_grouped_file0_col1_counts.linear.tsv",
                       "read_info.tsv.gz",
                       "splice_junction_grouped_file0_col1_counts.tsv",
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
    resulting_files = ["gene_counts.tsv", "read_info.tsv.gz", "transcript_counts.tsv",
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
                       "splice_junction_counts.tsv", "read_info.tsv.gz",
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
                       "read_info.tsv.gz",
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
                       "read_info.tsv.gz",
                       "transcript_counts.tsv",
                       "discovered_transcript_counts.tsv", "transcript_models.gtf"]
    for f in resulting_files:
        assert os.path.exists(os.path.join(sample_folder, sample_name + "." + f))


# def test_cage():
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
    # resulting_files = ["gene_counts.tsv", "read_info.tsv.gz", "transcript_counts.tsv",
    #                    "transcript_model_counts.tsv", "transcript_models.gtf", "transcript_model_reads.tsv.gz"]
    # for f in resulting_files:
    #     assert os.path.exists(os.path.join(sample_folder, sample_name + "." + f))


# ---------------------------------------------------------------------------
# Read grouping: reads absent from the grouping table
# ---------------------------------------------------------------------------

GROUPED_GENE_COUNTS = "gene_grouped_file0_col1_counts.linear.tsv"


def _sum_counts_per_group(counts_tsv):
    """feature_id/group_id/count linear file -> {group_id: total count}."""
    totals = {}
    with open(counts_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        assert header == ["feature_id", "group_id", "count"], header
        for line in f:
            feature_id, group_id, count = line.rstrip("\n").split("\t")
            if feature_id.startswith("__"):   # __ambiguous / __no_feature summary rows
                continue
            totals[group_id] = totals.get(group_id, 0.0) + float(count)
    return totals


def _run_grouped(out_dir, data_dir, groups_tsv, sample_name):
    result = subprocess.run(["python3", "isoquant.py",
                             "-o", out_dir,
                             "--data_type", "nanopore",
                             "--bam", os.path.join(data_dir, "chr9.4M.ont.sim.polya.bam"),
                             "--genedb", os.path.join(data_dir, "chr9.4M.gtf.gz"), "--complete_genedb",
                             "-r", os.path.join(data_dir, "chr9.4M.fa.gz"),
                             "-t", "2",
                             "--prefix", sample_name,
                             "--read_group", "file:%s:0:1" % groups_tsv])
    assert result.returncode == 0
    counts_tsv = os.path.join(out_dir, sample_name, sample_name + "." + GROUPED_GENE_COUNTS)
    assert os.path.exists(counts_tsv), counts_tsv
    return _sum_counts_per_group(counts_tsv)


def test_reads_missing_from_grouping_table_go_to_their_own_group(tmp_path):
    """Reads with no entry in the grouping table must not join a real group.

    Drops every GR1 row from the bundled grouping table and checks where those
    reads end up. The invariant that matters is that a real group's counts are
    *unchanged* by the removal: an unattributable read has to get its own label
    ('NA', what every grouper in read_groups.py falls back to), never be folded
    into an existing group.

    Regression guard for the group-id resolution defect: counters used to convert
    ids with a bare pool.get_str(group_id), so an unresolvable id indexed the pool
    list from the end and silently named a real group. See
    .claude/MULTI_GROUP_IMPLEMENTATION.md section B2.
    """
    source_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(source_dir, "simple_data/")
    full_tsv = os.path.join(data_dir, "chr9.4M.ont.sim.read_groups.tsv")

    # Same table with every GR1 read removed; those reads are then in the BAM but
    # absent from the grouping table.
    partial_tsv = str(tmp_path / "gr0_only.groups.tsv")
    dropped = 0
    with open(full_tsv) as src, open(partial_tsv, "w") as dst:
        for line in src:
            if line.rstrip("\n").split("\t")[1] == "GR1":
                dropped += 1
                continue
            dst.write(line)
    assert dropped > 0, "fixture no longer contains GR1 reads"

    full = _run_grouped(str(tmp_path / "out_full"), data_dir, full_tsv, "grouping_full")
    partial = _run_grouped(str(tmp_path / "out_partial"), data_dir, partial_tsv, "grouping_partial")

    # Baseline: both listed groups present. 'NA' may already appear, for reads in
    # the BAM that were never in the table to begin with.
    assert {"GR0", "GR1"} <= set(full)
    assert full["GR0"] > 0 and full["GR1"] > 0

    # GR1 is gone from the table, so it must be gone from the output -- and no new
    # group may be invented.
    assert set(partial) == {"GR0", "NA"}, partial

    # The invariant: GR0 is untouched by GR1's removal. Under the resolution bug an
    # unresolvable id named a real group, which would have inflated GR0 here.
    assert partial["GR0"] == pytest.approx(full["GR0"]), \
        "counts of a real group changed when unrelated reads lost their group"

    # The dropped reads are all accounted for under 'NA', none lost or double-counted.
    assert partial["NA"] == pytest.approx(full["GR1"] + full.get("NA", 0.0))
    assert sum(partial.values()) == pytest.approx(sum(full.values()))
