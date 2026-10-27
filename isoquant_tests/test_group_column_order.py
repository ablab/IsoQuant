############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os

from isoquant_lib.quantification.convert_grouped_counts import (
    order_groups, convert_to_matrix, convert_to_mtx,
)


def test_order_groups_default_sorted():
    assert order_groups(["b", "c", "a"]) == ["a", "b", "c"]


def test_order_groups_follows_given_order():
    assert order_groups(["a", "b", "c"], ["c", "a", "b"]) == ["c", "a", "b"]


def test_order_groups_skips_absent_and_appends_unlisted():
    # "x" has no counts, "NA" is not in the input list
    assert order_groups(["b", "NA", "a"], ["b", "x", "a"]) == ["b", "a", "NA"]


def _write_linear(path):
    with open(path, "w") as f:
        f.write("feature_id\tgroup_id\tcount\n")
        f.write("G1\talpha\t1.00\n")
        f.write("G1\tzeta\t2.00\n")
        f.write("G2\tmid\t3.00\n")


def test_convert_to_matrix_column_order(tmp_path):
    linear = str(tmp_path / "counts.linear.tsv")
    _write_linear(linear)
    prefix = str(tmp_path / "counts")
    convert_to_matrix(linear, prefix, group_order=["zeta", "alpha", "mid"])
    with open(prefix + ".tsv") as f:
        assert f.readline().rstrip("\n").split("\t") == ["gene_id", "zeta", "alpha", "mid"]
        assert f.readline().rstrip("\n").split("\t") == ["G1", "2.0", "1.0", "0"]

    convert_to_matrix(linear, prefix)
    with open(prefix + ".tsv") as f:
        assert f.readline().rstrip("\n").split("\t") == ["gene_id", "alpha", "mid", "zeta"]


def test_convert_to_mtx_barcode_order(tmp_path):
    linear = str(tmp_path / "counts.linear.tsv")
    _write_linear(linear)
    prefix = str(tmp_path / "counts")
    convert_to_mtx(linear, prefix, group_order=["zeta", "alpha", "mid"])
    barcodes = [f for f in os.listdir(tmp_path) if "barcode" in f]
    assert len(barcodes) == 1
    with open(tmp_path / barcodes[0]) as f:
        assert [line.strip() for line in f] == ["zeta", "alpha", "mid"]
