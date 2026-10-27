############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Unit tests for the RNA velocity spliced / unspliced / ambiguous counter.

These pin the splicing verdict (which match events and classifications make a
read unspliced, and which reads are dropped outright), the per-chromosome TSV
fragment layout that merge_counts() concatenates, and the loom export.
"""

import os
from types import SimpleNamespace

import pytest

from isoquant_lib.assignment.isoform_assignment import (
    MatchClassification,
    MatchEventSubtype,
    ReadAssignmentType,
)
from isoquant_lib.quantification.rna_velocity_counter import (
    INTRON_RETENTION_EVENTS,
    RNAVelocityCounter,
    SplicingStatus,
)
from isoquant_lib.utils.string_pools import UNASSIGNED_GROUP_ID, UNASSIGNED_GROUP_NAME


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


class FakeEvent:
    def __init__(self, event_type):
        self.event_type = event_type


class FakeMatch:
    def __init__(self, gene_id, events=(), classification=MatchClassification.full_splice_match):
        self.assigned_gene_id = gene_id
        self.match_subclassifications = [FakeEvent(e) for e in events]
        self.match_classification = classification


class FakeAssignment:
    def __init__(self, assignment_type, matches, read_group_ids=(0,)):
        self.assignment_type = assignment_type
        self.isoform_matches = list(matches)
        self.read_group_ids = list(read_group_ids)


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


def one_row(counter):
    rows = read_rows(counter)
    assert len(rows) == 1, rows
    return rows[0]


# --------------------------------------------------------------- splicing verdict

@pytest.mark.parametrize("event", sorted(INTRON_RETENTION_EVENTS, key=lambda e: e.name))
def test_intron_retention_events_make_a_read_unspliced(tmp_path, event):
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.inconsistent,
                                         [FakeMatch(0, events=[event])]))
    counter.dump()
    assert one_row(counter) == ("CELL1", "GENE_A", "0", "1", "0")


def test_incomplete_intron_retention_is_included():
    # The 50 bp floor lives in junction_comparator (minor_exon_extension), so the
    # counter takes these events at face value; pin that they are in the set.
    assert MatchEventSubtype.incomplete_intron_retention_left in INTRON_RETENTION_EVENTS
    assert MatchEventSubtype.incomplete_intron_retention_right in INTRON_RETENTION_EVENTS


def test_fake_micro_intron_retention_is_not_unspliced(tmp_path):
    # An alignment artifact (is_alignment_artifact), not real intronic coverage.
    assert MatchEventSubtype.fake_micro_intron_retention not in INTRON_RETENTION_EVENTS
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(
        ReadAssignmentType.unique_minor_difference,
        [FakeMatch(0, events=[MatchEventSubtype.fake_micro_intron_retention])]))
    counter.dump()
    assert one_row(counter) == ("CELL1", "GENE_A", "1", "0", "0")


def test_genic_intron_read_is_unspliced_despite_noninformative(tmp_path):
    # The pre-mRNA case: assignment_type stays noninformative, the gene is only
    # recorded on genic_intron matches. This is the read class the old
    # consistency-based rule dropped entirely.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(
        ReadAssignmentType.noninformative,
        [FakeMatch(0, classification=MatchClassification.genic_intron)]))
    counter.dump()
    assert one_row(counter) == ("CELL1", "GENE_A", "0", "1", "0")


def test_noninformative_without_genic_intron_is_dropped(tmp_path):
    # Overlaps the gene body but resembles no isoform -- undecidable, not unspliced.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(
        ReadAssignmentType.noninformative,
        [FakeMatch(0, classification=MatchClassification.genic)]))
    counter.dump()
    assert read_rows(counter) == []


@pytest.mark.parametrize("assignment_type", [
    ReadAssignmentType.unique,
    ReadAssignmentType.unique_minor_difference,
    ReadAssignmentType.ambiguous,
    ReadAssignmentType.inconsistent,
    ReadAssignmentType.inconsistent_ambiguous,
    ReadAssignmentType.inconsistent_non_intronic,
])
def test_matched_reads_without_intronic_evidence_are_spliced(tmp_path, assignment_type):
    # The point of the rewrite: disagreeing with the annotation does not make a
    # read unspliced. A novel splice site is still a spliced molecule.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(
        assignment_type,
        [FakeMatch(0, events=[MatchEventSubtype.alt_left_site_novel])]))
    counter.dump()
    assert one_row(counter) == ("CELL1", "GENE_A", "1", "0", "0")


def test_mixed_verdicts_are_ambiguous(tmp_path):
    # Retained against one isoform, mature against another.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.inconsistent_ambiguous, [
        FakeMatch(0, events=[MatchEventSubtype.intron_retention]),
        FakeMatch(0, events=[MatchEventSubtype.fsm]),
    ]))
    counter.dump()
    assert one_row(counter) == ("CELL1", "GENE_A", "0", "0", "1")


def test_all_matches_unspliced_is_not_ambiguous(tmp_path):
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.inconsistent_ambiguous, [
        FakeMatch(0, events=[MatchEventSubtype.intron_retention]),
        FakeMatch(0, events=[MatchEventSubtype.unspliced_intron_retention]),
    ]))
    counter.dump()
    assert one_row(counter) == ("CELL1", "GENE_A", "0", "1", "0")


def test_undecidable_matches_do_not_dilute_a_verdict(tmp_path):
    # A genic match carries no splicing signal; it must not turn an otherwise
    # unspliced read into an ambiguous one.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.inconsistent, [
        FakeMatch(0, events=[MatchEventSubtype.intron_retention]),
        FakeMatch(0, classification=MatchClassification.genic),
    ]))
    counter.dump()
    assert one_row(counter) == ("CELL1", "GENE_A", "0", "1", "0")


def test_match_status_helper():
    assert RNAVelocityCounter._match_status(
        FakeMatch(0, events=[MatchEventSubtype.intron_retention])) == SplicingStatus.unspliced
    assert RNAVelocityCounter._match_status(
        FakeMatch(0, classification=MatchClassification.genic_intron)) == SplicingStatus.unspliced
    assert RNAVelocityCounter._match_status(
        FakeMatch(0, events=[MatchEventSubtype.fsm])) == SplicingStatus.spliced
    assert RNAVelocityCounter._match_status(
        FakeMatch(0, classification=MatchClassification.genic)) is None


# --------------------------------------------------------------- dropped reads

@pytest.mark.parametrize("assignment_type", [
    ReadAssignmentType.intergenic,
    ReadAssignmentType.discarded,
    ReadAssignmentType.suspended,
])
def test_unaccepted_types_are_dropped(tmp_path, assignment_type):
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(assignment_type, [FakeMatch(0)]))
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
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [FakeMatch(0)],
                                         read_group_ids=[]))
    counter.dump()
    assert read_rows(counter) == []


def test_legacy_sentinel_group_id_is_dropped(tmp_path):
    # Assignments serialized before read_group_to_ids started interning
    # UNASSIGNED_GROUP_NAME can still carry the negative sentinel, which would
    # index the pool list from the end and silently pick a real barcode.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [FakeMatch(0)],
                                         read_group_ids=[UNASSIGNED_GROUP_ID]))
    counter.dump()
    assert read_rows(counter) == []


def test_unassigned_group_is_not_emitted_as_a_cell(tmp_path):
    """A read with no barcode must not become a cell called 'NA'.

    BarcodeGrouper returns UNASSIGNED_GROUP_NAME for any read carrying no barcode
    (e.g. absent from --barcoded_reads), and since the group-id unification an
    unset value interns as that same name -- so it arrives as an ordinary, valid
    pool id. Velocity output is per cell, so that pseudo-cell is dropped.
    """
    counter = make_counter(tmp_path, cells=("CELL1", UNASSIGNED_GROUP_NAME))
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [FakeMatch(0)],
                                         read_group_ids=[0]))
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [FakeMatch(0)],
                                         read_group_ids=[1]))   # the 'NA' bucket
    counter.dump()
    assert read_rows(counter) == [("CELL1", "GENE_A", "1", "0", "0")]


def test_unassigned_group_is_kept_out_of_the_loom(tmp_path):
    loompy = pytest.importorskip("loompy")
    counter = make_counter(tmp_path, cells=("CELL1", UNASSIGNED_GROUP_NAME))
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [FakeMatch(0)],
                                         read_group_ids=[0]))
    for _ in range(5):
        counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [FakeMatch(0)],
                                             read_group_ids=[1]))
    counter.dump()
    counter.finalize()
    with loompy.connect(counter.output_file + ".loom") as ds:
        assert list(ds.ca.CellID) == ["CELL1"]
        assert ds.layers["spliced"][:, :].sum() == 1


def test_read_without_gene_is_dropped(tmp_path):
    # match_inconsistent's quick_mode exit yields an IsoformMatch with no gene;
    # it used to emit a row with an empty gene_id column.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.inconsistent, [FakeMatch(None)]))
    counter.dump()
    assert read_rows(counter) == []


def test_first_named_gene_is_used(tmp_path):
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique,
                                         [FakeMatch(None), FakeMatch(1)]))
    counter.dump()
    assert one_row(counter) == ("CELL1", "GENE_B", "1", "0", "0")


def test_group_index_selects_the_right_strategy(tmp_path):
    counter = make_counter(tmp_path, group_index=1)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [FakeMatch(0)],
                                         read_group_ids=[0, 1]))
    counter.dump()
    assert one_row(counter) == ("CELL2", "GENE_A", "1", "0", "0")


# --------------------------------------------------------------- output plumbing

def test_counts_accumulate_per_cell_and_gene(tmp_path):
    counter = make_counter(tmp_path)
    for _ in range(3):
        counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [FakeMatch(0)]))
    counter.add_read_info(FakeAssignment(
        ReadAssignmentType.inconsistent,
        [FakeMatch(0, events=[MatchEventSubtype.intron_retention])]))
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [FakeMatch(1)],
                                         read_group_ids=[1]))
    counter.dump()
    assert read_rows(counter) == [("CELL1", "GENE_A", "3", "1", "0"),
                                  ("CELL2", "GENE_B", "1", "0", "0")]


def test_dump_clears_state(tmp_path):
    # dump() is called once per chromosome in the worker, but a second dump must
    # not re-emit the same counts into the fragment.
    counter = make_counter(tmp_path)
    counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [FakeMatch(0)]))
    counter.dump()
    counter.dump()
    assert one_row(counter) == ("CELL1", "GENE_A", "1", "0", "0")


def test_finalize_on_empty_output_does_not_fail(tmp_path):
    counter = make_counter(tmp_path)
    counter.dump()
    counter.finalize()
    assert not os.path.exists(counter.output_file + ".loom")


def test_finalize_writes_loom_layers(tmp_path):
    loompy = pytest.importorskip("loompy")
    counter = make_counter(tmp_path)
    for _ in range(2):
        counter.add_read_info(FakeAssignment(ReadAssignmentType.unique, [FakeMatch(0)]))
    counter.add_read_info(FakeAssignment(
        ReadAssignmentType.inconsistent,
        [FakeMatch(0, events=[MatchEventSubtype.intron_retention])]))
    counter.add_read_info(FakeAssignment(ReadAssignmentType.inconsistent_ambiguous, [
        FakeMatch(0, events=[MatchEventSubtype.intron_retention]),
        FakeMatch(0, events=[MatchEventSubtype.fsm]),
    ]))
    counter.add_read_info(FakeAssignment(ReadAssignmentType.ambiguous, [FakeMatch(1)],
                                         read_group_ids=[1]))
    counter.dump()
    counter.finalize()

    loom_file = counter.output_file + ".loom"
    assert os.path.exists(loom_file)
    with loompy.connect(loom_file) as ds:
        assert set(ds.layers.keys()) == {"", "spliced", "unspliced", "ambiguous"}
        genes = list(ds.ra.Gene)
        cells = list(ds.ca.CellID)
        i, j = genes.index("GENE_A"), cells.index("CELL1")
        assert ds.layers["spliced"][i][j] == 2
        assert ds.layers["unspliced"][i][j] == 1
        assert ds.layers["ambiguous"][i][j] == 1
        i, j = genes.index("GENE_B"), cells.index("CELL2")
        assert ds.layers["spliced"][i][j] == 1
        assert ds.layers["unspliced"][i][j] == 0


def test_finalize_sums_duplicate_rows_across_fragments(tmp_path):
    # Merged fragments may repeat a (cell, gene) pair; coo_matrix must sum them.
    loompy = pytest.importorskip("loompy")
    counter = make_counter(tmp_path)
    with open(counter.output_file, "w") as f:
        f.write("#cell_id\tgene_id\tspliced\tunspliced\tambiguous\n")
        f.write("CELL1\tGENE_A\t2\t1\t1\n")
        f.write("CELL1\tGENE_A\t3\t4\t2\n")
    counter.finalize()
    with loompy.connect(counter.output_file + ".loom") as ds:
        assert ds.layers["spliced"][0][0] == 5
        assert ds.layers["unspliced"][0][0] == 5
        assert ds.layers["ambiguous"][0][0] == 3


def test_loom_is_written_in_several_column_blocks(tmp_path, monkeypatch):
    # The loom layers are filled one column window at a time, so that peak
    # memory tracks the window and not the cell count. Shrink the window so a
    # small fixture spans several blocks, and check nothing is lost at the seams.
    loompy = pytest.importorskip("loompy")
    from isoquant_lib.quantification import rna_velocity_counter as rvc

    n_cells = 7
    monkeypatch.setattr(rvc, "LOOM_WINDOW_BYTES", 8)  # 2 cells per block at int32
    counter = make_counter(tmp_path)
    with open(counter.output_file, "w") as f:
        f.write("#cell_id\tgene_id\tspliced\tunspliced\tambiguous\n")
        for i in range(n_cells):
            f.write("CELL%d\tGENE_A\t%d\t%d\t%d\n" % (i, i, 100 + i, 200 + i))
    counter.finalize()

    with loompy.connect(counter.output_file + ".loom") as ds:
        assert ds.shape == (1, n_cells)
        cells = list(ds.ca.CellID)
        assert cells == ["CELL%d" % i for i in range(n_cells)]
        assert list(ds.ra.Gene) == ["GENE_A"]
        for i in range(n_cells):
            j = cells.index("CELL%d" % i)
            assert ds.layers["spliced"][0][j] == i
            assert ds.layers["unspliced"][0][j] == 100 + i
            assert ds.layers["ambiguous"][0][j] == 200 + i
        # The unnamed main layer mirrors spliced, per velocyto convention.
        assert (ds.layers[""][:, :] == ds.layers["spliced"][:, :]).all()
