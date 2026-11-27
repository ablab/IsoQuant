############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from argparse import Namespace

import pandas as pd
import pytest

from isoquant_lib.assignment.isoform_assignment import (
    IsoformMatch,
    MatchClassification,
    ReadAssignment,
    ReadAssignmentType,
)
from isoquant_lib.quantification.convert_grouped_counts import (
    convert_to_matrix,
    convert_to_mtx,
    get_normalization_factors,
)
from isoquant_lib.quantification.long_read_counter import create_gene_counter
from isoquant_lib.utils.string_pools import StringPoolManager


def test_get_normalization_factors_with_zero_usable_reads_group():
    df = pd.DataFrame({"gene_id": ["g1"], "group_id": ["sample1"], "count": [50.0]})
    factors = get_normalization_factors(df, usable_reads_per_group={"sample1": 200, "empty_sample": 0})
    assert factors["sample1"] == pytest.approx(5000.0)
    assert factors["empty_sample"] == 0.0


def test_get_normalization_factors_with_zero_count_group():
    df = pd.DataFrame(
        {
            "gene_id": ["g1", "g1"],
            "group_id": ["sample1", "sample2"],
            "count": [50.0, 0.0],
        }
    )
    factors = get_normalization_factors(df)
    assert factors["sample1"] == pytest.approx(20000.0)
    assert factors["sample2"] == 0.0


def test_convert_to_matrix_and_mtx_tpm_with_zero_usable_reads_group(tmp_path):
    linear_file = tmp_path / "counts.linear.tsv"
    linear_file.write_text("feature_id\tgroup_id\tcount\nGENE1\tsample1\t50.0\n")

    matrix_prefix = str(tmp_path / "out_tpm")
    convert_to_matrix(
        str(linear_file),
        matrix_prefix,
        convert_to_tpm=True,
        usable_reads_per_group={"sample1": 200, "empty_sample": 0},
    )
    lines = (tmp_path / "out_tpm.tsv").read_text().splitlines()
    assert lines[0] == "gene_id\tsample1"
    gene_id, val = lines[1].split("\t")
    assert gene_id == "GENE1"
    assert float(val) == pytest.approx(250000.0)

    mtx_prefix = str(tmp_path / "out_mtx_tpm")
    convert_to_mtx(
        str(linear_file),
        mtx_prefix,
        convert_to_tpm=True,
        usable_reads_per_group={"sample1": 200, "empty_sample": 0},
    )
    mtx_lines = (tmp_path / "out_mtx_tpm.matrix.mtx").read_text().splitlines()
    assert mtx_lines[2].split()[0:2] == ["1", "1"]
    assert float(mtx_lines[2].split()[2]) == pytest.approx(250000.0)


def test_assigned_feature_counter_finalize_usable_reads_with_empty_group(tmp_path):
    pools = StringPoolManager()
    pools.gene_pool.add("GENE1")
    pools.transcript_pool.add("TX1")
    pools.group_spec_pool_types[0] = "file_name"
    sample1_id = pools.file_name_pool.add("sample1")
    pools.file_name_pool.add("sample2")  # second group present in pool with 0 reads

    counter = create_gene_counter(
        str(tmp_path / "gene_grouped"),
        strategy="unique_only",
        string_pools=pools,
        group_index=0,
    )
    match = IsoformMatch(
        MatchClassification.full_splice_match,
        pools,
        assigned_gene="GENE1",
        assigned_transcript="TX1",
    )
    ra = ReadAssignment("read_1", ReadAssignmentType.unique, pools, match=match)
    ra.read_group_ids = [sample1_id]
    counter.add_read_info(ra)
    counter.dump()

    args = Namespace(normalization_method="usable_reads", counts_format=["matrix", "mtx"])
    counter.finalize(args)

    tpm_tsv = tmp_path / "gene_grouped_tpm.tsv"
    assert tpm_tsv.exists()
    lines = tpm_tsv.read_text().splitlines()
    assert lines[0] == "gene_id\tsample1"
    assert float(lines[1].split("\t")[1]) == pytest.approx(1000000.0)
