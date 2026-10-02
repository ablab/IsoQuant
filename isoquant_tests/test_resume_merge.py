############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Pieces the resumable merge relies on (see .claude/RESUME_CHECKPOINTS.md).

- the merge aggregator (truncate_outputs=False) must not touch any finished output;
- per-chromosome fragments live in a separate directory, named like the writers name them;
- fragments survive a merge with remove_inputs=False, so a rerun can merge them again;
- .usable fragments are loaded once, by the finalize step only.
"""

from __future__ import annotations

import os
from types import SimpleNamespace

import pytest

from isoquant_lib.assignment.assignment_aggregator import ReadAssignmentAggregator
from isoquant_lib.modes import IsoQuantMode
from isoquant_lib.quantification.long_read_counter import create_gene_counter
from isoquant_lib.quantification.rna_velocity_counter import RNAVelocityCounter
from isoquant_lib.utils.file_utils import (load_usable_fragments, merge_counts, merge_file_list, merge_files,
                                          usable_fragment_list)
from isoquant_lib.utils.input_data_storage import SampleData
from isoquant_lib.utils.string_pools import StringPoolManager


def make_args() -> SimpleNamespace:
    return SimpleNamespace(
        _cmd_line="isoquant.py --test", _version="test", genedb="annotation.db", large_output=[],
        sqanti_output=False, run_quantification=True, gene_quantification="with_ambiguous",
        transcript_quantification="with_ambiguous", no_model_construction=False, count_exons=True,
        old_exon_count_format=True, delta=6, minimal_exon_overlap=5, predict_terminal_sites=True,
        fl_data=True, count_intron_retentions=True, read_group=["file_name", "barcode"],
        mode=IsoQuantMode.tenX_v3, collect_polya_training="training.csv", collect_tss_training=None)


def make_pools() -> StringPoolManager:
    pools = StringPoolManager()
    pools.set_group_spec_pool_type(0, "file_name")
    pools.set_group_spec_pool_type(1, "barcode")
    return pools


def snapshot(directory: str) -> dict:
    result = {}
    for root, _, files in os.walk(directory):
        for f in files:
            path = os.path.join(root, f)
            with open(path) as inf:
                result[path] = inf.read()
    return result


class TestMergeAggregatorHasNoSideEffects:
    def test_finished_outputs_stay_intact(self, tmp_path):
        out_dir = str(tmp_path / "sample")
        os.makedirs(out_dir)
        sample = SampleData([], "sample", out_dir, {}, None)
        args = make_args()
        strategies = ["file_name", "barcode"]
        old_cwd = os.getcwd()
        os.chdir(str(tmp_path))
        try:
            # a fresh merge creates every output; fill them as a finished merge would have
            ReadAssignmentAggregator(args, sample, make_pools(), grouping_strategy_names=strategies)
            for path in snapshot(out_dir):
                with open(path, "w") as outf:
                    outf.write("finished output\n")
            before = snapshot(out_dir)
            assert len(before) > 10

            aggregator = ReadAssignmentAggregator(args, sample, make_pools(), grouping_strategy_names=strategies,
                                                  truncate_outputs=False)
            assert snapshot(out_dir) == before
            assert aggregator.read_info_printer is None and aggregator.basic_printer is None
            # every counter kind was built, including the grouped model, terminal and velocity counters
            names = [os.path.basename(c.output_counts_file_name)
                     for c in aggregator.global_counter.counters + aggregator.transcript_model_global_counter.counters
                     + aggregator.gene_model_global_counter.counters]
            assert any("discovered_transcript_grouped" in n for n in names)
            assert any("RNA_velocity" in n for n in names)
            assert any("polyA_prediction" in n for n in names)
        finally:
            os.chdir(old_cwd)

    def test_velocity_counter_keeps_output(self, tmp_path):
        path = str(tmp_path / "velocity")
        with open(path, "w") as outf:
            outf.write("finished\n")
        RNAVelocityCounter(make_args(), path, string_pools=make_pools(), group_index=1, truncate_output=False)
        with open(path) as inf:
            assert inf.read() == "finished\n"


class TestFragmentDir:
    def test_names_go_to_fragment_dir(self):
        assert merge_file_list("/out/s/s.gene_counts.tsv", "s", ["chr1"], "/out/s/aux/per_chr") == \
            ["/out/s/aux/per_chr/s_chr1.gene_counts.tsv"]

    def test_chromosome_ids_are_sanitised_like_the_writers(self):
        assert merge_file_list("/out/s.gene_counts.tsv", "s", ["chrUn/random"]) == \
            ["/out/s_chrUn_random.gene_counts.tsv"]

    def test_label_must_prefix_basename_with_fragment_dir(self):
        with pytest.raises(ValueError):
            merge_file_list("/out/other.tsv", "s", ["chr1"], "/frag")

    def test_merge_keeps_inputs_when_asked(self, tmp_path):
        frag = tmp_path / "per_chr"
        frag.mkdir()
        for chr_id in ("chr1", "chr2"):
            (frag / ("s_%s.tsv" % chr_id)).write_text("#h\n%s\n" % chr_id)
        final = str(tmp_path / "s.tsv")
        for _ in range(2):  # a rerun merges the same fragments again
            with open(final, "w") as outf:
                consumed = merge_files(final, "s", ["chr1", "chr2"], outf, remove_inputs=False,
                                       fragment_dir=str(frag))
            assert all(os.path.exists(f) for f in consumed)
        assert open(final).read() == "#h\nchr1\nchr2\n"


class TestUsableLoadedOnce:
    def test_merge_counts_leaves_usable_to_finalize(self, tmp_path):
        frag = str(tmp_path / "per_chr")
        os.makedirs(frag)
        final_prefix = str(tmp_path / "s.gene")
        counter = create_gene_counter(final_prefix, "with_ambiguous", truncate_output=True)
        for chr_id in ("chr1", "chr2"):
            chr_counter = create_gene_counter(os.path.join(frag, "s_%s.gene" % chr_id), "with_ambiguous")
            with open(chr_counter.output_counts_file_name, "w") as outf:
                outf.write("#feature_id\tcount\nG1\t1.00\n")
            with open(chr_counter.usable_file_name, "w") as outf:
                outf.write("default\t5\n")
            with open(chr_counter.output_stats_file_name, "w") as outf:
                outf.write("__ambiguous\t0\n")

        merge_counts(counter, "s", ["chr1", "chr2"], remove_inputs=False, fragment_dir=frag)
        assert counter.reads_for_tpm[0] == 0
        load_usable_fragments(counter, "s", ["chr1", "chr2"], frag)
        assert counter.reads_for_tpm[0] == 10
        assert all(os.path.exists(f) for f in usable_fragment_list(counter, "s", ["chr1", "chr2"], frag))
