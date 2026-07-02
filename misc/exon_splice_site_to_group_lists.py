#!/usr/bin/env python3

############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Convert exon splice-site counts from the default per-group format to the
per-molecule group-list format.

Default (per-group) input, one row per (feature, group):
    region_gene_candidate  n_full  n_left  n_right  group_id  [read_ids_full read_ids_left read_ids_right]

Group-list output, one row per feature:
    region_gene_candidate  n_full  n_left  n_right  groups_full  groups_left  groups_right  [read_ids_full read_ids_left read_ids_right]

For each feature, groups_<bucket> lists the group of every contributing read (the
group id repeated n_<bucket> times, concatenated across the feature's group rows);
with the `barcode` grouping strategy this reproduces the per-molecule barcode
ledger. read_ids_* (when present) are carried through in the matching order.

Assumes all rows of a feature are contiguous (IsoQuant output is sorted this way).

Usage:
    python misc/exon_splice_site_to_group_lists.py \
        -i SAMPLE.exon_splice_site_counts.tsv \
        -o SAMPLE.exon_splice_site_group_lists.tsv
"""

import argparse
import gzip
from typing import TextIO


def _open(path: str, mode: str = "r") -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def _split(s: str) -> list:
    return [] if s == "NA" else s.split(";")


def _join(items: list) -> str:
    return ";".join(items) if items else "NA"


def _write_feature(out: TextIO, feature_id: str, acc: list, emit_read_ids: bool) -> None:
    gf, gl, gr, rf, rl, rr = acc
    fields = [feature_id, str(len(gf)), str(len(gl)), str(len(gr)),
              _join(gf), _join(gl), _join(gr)]
    if emit_read_ids:
        fields += [_join(rf), _join(rl), _join(rr)]
    out.write("\t".join(fields) + "\n")


def convert(input_path: str, output_path: str) -> None:
    with _open(input_path) as inf, _open(output_path, "w") as out:
        header = inf.readline().rstrip("\n").split("\t")
        emit_read_ids = "read_ids_full" in header

        out_header = ["region_gene_candidate", "n_full", "n_left", "n_right",
                      "groups_full", "groups_left", "groups_right"]
        if emit_read_ids:
            out_header += ["read_ids_full", "read_ids_left", "read_ids_right"]
        out.write("\t".join(out_header) + "\n")

        cur_feature = None
        acc = None  # [groups_full, groups_left, groups_right, rids_full, rids_left, rids_right]
        for line in inf:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            feature_id, gid = cols[0], cols[4]
            nf, nl, nr = int(cols[1]), int(cols[2]), int(cols[3])
            if feature_id != cur_feature:
                if cur_feature is not None:
                    _write_feature(out, cur_feature, acc, emit_read_ids)
                cur_feature = feature_id
                acc = [[], [], [], [], [], []]
            acc[0] += [gid] * nf
            acc[1] += [gid] * nl
            acc[2] += [gid] * nr
            if emit_read_ids:
                acc[3] += _split(cols[5])
                acc[4] += _split(cols[6])
                acc[5] += _split(cols[7])
        if cur_feature is not None:
            _write_feature(out, cur_feature, acc, emit_read_ids)


def main():
    parser = argparse.ArgumentParser(
        description="Convert exon splice-site per-group counts to per-molecule group lists")
    parser.add_argument("--input", "-i", required=True,
                        help="input per-group counts file (.tsv[.gz])")
    parser.add_argument("--output", "-o", required=True,
                        help="output group-list file (.tsv[.gz])")
    args = parser.parse_args()
    convert(args.input, args.output)


if __name__ == "__main__":
    main()
