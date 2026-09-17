############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Unit tests for the RNA velocity spliced / unspliced counter.

These pin the read-acceptance rules (which assignment types count as spliced
vs unspliced, and which reads are dropped outright), the per-chromosome TSV
fragment layout that merge_counts() concatenates, and the loom export.
"""

import os
from types import SimpleNamespace

import pytest

from isoquant_lib.assignment.isoform_assignment import ReadAssignmentType
from isoquant_lib.quantification.rna_velocity_counter import RNAVelocityCounter


class FakePool:
    def __init__(self, items):
        self.int_to_str = list(items)

    def get_str(self, int_id):
        return self.int_to_str[int_id]


class FakeStringPools:
    def __init__(self, cells=("CELL1", "CELL2"), genes=("GENE_A", "GENE_B")):
        self.gene_pool = FakePool(genes)
        self._group_pool = FakePool(cells)

    def get_read_group_pool(self, spec_index):
        return self._group_pool


class FakeMatch:
    def __init__(self, gene_id):
        self.assigned_gene_id = gene_id


class FakeAssignment:
    def __init__(self, assignment_type, gene_ids, read_group_ids):
        self.assignment_type = assignment_type
        self.isoform_matches = [FakeMatch(g) for g in gene_ids]
        self.read_group_ids = read_group_ids


def make_counter(tmp_path, group_index=0, cells=("CELL1", "CELL2")):
    out = str(tmp_path / "sample.RNA_velocity_grouped_barcode")
    return RNAVelocityCounter(SimpleNamespace(), out,
                              string_pools=FakeStringPools(cells=cells),
                              group_index=group_index)


def read_rows(counter):
    with open(counter.output_file) as f:
        lines = [l.strip() for l in f if l.strip()]
    assert lines[0].startswith("#"), "fragment must carry the header merge_counts() strips"
    return sorted(tuple(l.split("\t")) for l in lines[1:])


SPLICED_TYPES = [ReadAssignmentType.unique,
                 ReadAssignmentType.unique_minor_difference,
                 ReadAssignmentType.ambiguous]
UNSPLICED_TYPES = [ReadAssignmentType.inconsistent,
                   ReadAssignmentType.inconsistent_ambiguous]
IGNORED_TYPES = [ReadAssignmentType.noninformative,
                 ReadAssignmentType.intergenic,
                 ReadAssignmentType.inconsistent_non_intronic,
                 ReadAssignmentType.inconsistent_genic,
                 ReadAssignmentType.inconsistent_multigenic]


@pytest.mark.parametrize("assignment_type", SPLICED_TYPES)
def test_consistent_reads_count_as_spliced(tmp_path, assignment_type):
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(assignment_type, [0], [0]))
    counter.dump()
    assert read_rows(counter) == [("CELL1", "GENE_A", "1", "0")]


@pytest.mark.parametrize("assignment_type", UNSPLICED_TYPES)
def test_inconsistent_reads_count_as_unspliced(tmp_path, assignment_type):
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(assignment_type, [0], [0]))
    counter.dump()
    assert read_rows(counter) == [("CELL1", "GENE_A", "0", "1")]


@pytest.mark.parametrize("assignment_type", IGNORED_TYPES)
def test_unaccepted_types_are_dropped(tmp_path, assignment_type):
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(assignment_type, [0], [0]))
    counter.dump()
    assert read_rows(counter) == []


def test_none_assignment_is_ignored(tmp_path):
    counter = make_counter(tmp_path)
    counter.add_read_info(None)
    counter.dump()
    assert read_rows(counter) == []


def test_read_without_barcode_is_dropped(tmp_path):
    # Reads whose barcode was not detected carry no group id at all; indexing
    # read_group_ids unconditionally used to raise IndexError inside the worker.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [0], []))
    counter.dump()
    assert read_rows(counter) == []


def test_missing_group_for_this_strategy_is_dropped(tmp_path):
    # read_group_to_ids stores -1 when a strategy produced no value; -1 would
    # index the pool list from the end and silently pick a real barcode.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [0], [-1]))
    counter.dump()
    assert read_rows(counter) == []


def test_read_without_gene_is_dropped(tmp_path):
    # match_inconsistent's quick_mode exit yields an IsoformMatch with no gene;
    # it used to emit a row with an empty gene_id column.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.inconsistent, [None], [0]))
    counter.dump()
    assert read_rows(counter) == []


def test_first_named_gene_is_used(tmp_path):
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [None, 1], [0]))
    counter.dump()
    assert read_rows(counter) == [("CELL1", "GENE_B", "1", "0")]


def test_group_index_selects_the_right_strategy(tmp_path):
    counter = make_counter(tmp_path, group_index=1)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [0], [0, 1]))
    counter.dump()
    assert read_rows(counter) == [("CELL2", "GENE_A", "1", "0")]


def test_counts_accumulate_per_cell_and_gene(tmp_path):
    counter = make_counter(tmp_path)
    for _ in range(3):
        counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [0], [0]))
    counter.add_read_info(FakeAssignment(ReadAssignmentType.inconsistent, [0], [0]))
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [1], [1]))
    counter.dump()
    assert read_rows(counter) == [("CELL1", "GENE_A", "3", "1"),
                                  ("CELL2", "GENE_B", "1", "0")]


def test_dump_clears_state(tmp_path):
    # dump() is called once per chromosome in the worker, but a second dump must
    # not re-emit the same counts into the fragment.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [0], [0]))
    counter.dump()
    counter.dump()
    assert read_rows(counter) == [("CELL1", "GENE_A", "1", "0")]


def test_finalize_on_empty_output_does_not_fail(tmp_path):
    counter = make_counter(tmp_path)
    counter.dump()
    counter.finalize()
    assert not os.path.exists(counter.output_file + ".loom")


def test_finalize_writes_loom_layers(tmp_path):
    loompy = pytest.importorskip("loompy")
    counter = make_counter(tmp_path)
    for _ in range(2):
        counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [0], [0]))
    counter.add_read_info(FakeAssignment(ReadAssignmentType.inconsistent, [0], [0]))
    counter.add_read_info(FakeAssignment(ReadAssignmentType.ambiguous, [1], [1]))
    counter.dump()
    counter.finalize()

    loom_file = counter.output_file + ".loom"
    assert os.path.exists(loom_file)
    with loompy.connect(loom_file) as ds:
        assert set(ds.layers.keys()) == {"", "spliced", "unspliced"}
        genes = list(ds.ra.Gene)
        cells = list(ds.ca.CellID)
        spliced = ds.layers["spliced"][:, :]
        unspliced = ds.layers["unspliced"][:, :]
        i, j = genes.index("GENE_A"), cells.index("CELL1")
        assert spliced[i][j] == 2
        assert unspliced[i][j] == 1
        i, j = genes.index("GENE_B"), cells.index("CELL2")
        assert spliced[i][j] == 1
        assert unspliced[i][j] == 0


def test_finalize_sums_duplicate_rows_across_fragments(tmp_path):
    # Merged fragments may repeat a (cell, gene) pair; coo_matrix must sum them.
    loompy = pytest.importorskip("loompy")
    counter = make_counter(tmp_path)
    with open(counter.output_file, "w") as f:
        f.write("#cell_id\tgene_id\tspliced\tunspliced\n")
        f.write("CELL1\tGENE_A\t2\t1\n")
        f.write("CELL1\tGENE_A\t3\t4\n")
    counter.finalize()
    with loompy.connect(counter.output_file + ".loom") as ds:
        assert ds.layers["spliced"][0][0] == 5
        assert ds.layers["unspliced"][0][0] == 5
