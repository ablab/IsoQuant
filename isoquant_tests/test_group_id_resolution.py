############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Read-group id -> name resolution across every grouped counter.

Regression tests for a defect where five counter classes each carried their own
copy of the int->str conversion and none of them handled the values
StringPoolManager can produce: a negative "unassigned" sentinel resolved to the
*last group in the pool* (a real barcode), and a stale id raised an uncaught
IndexError where the canonical converter recovers with a warning.

The conversion now lives once on AbstractCounter, delegating to
StringPoolManager.resolve_read_group. These tests pin that every counter agrees.
"""

import os
import tempfile
from types import SimpleNamespace

import pytest

from isoquant_lib.quantification.long_read_counter import (
    AbstractCounter,
    ExonCounter,
    ExonSpliceSiteCounter,
    IntronCounter,
    IntronRetentionCounter,
    JointExonCounter,
    create_gene_counter,
    create_transcript_counter,
)
from isoquant_lib.terminal_prediction.terminal_counter import PolyACounter
from isoquant_lib.utils.string_pools import (
    StringPoolManager,
    UNASSIGNED_GROUP_ID,
    UNASSIGNED_GROUP_NAME,
)

BARCODES = ["AAACCCAAGAAACACT", "AAACCCAAGAAACCAT", "TTTGTTGGTTTGGCTA"]


def make_pools():
    pools = StringPoolManager()
    pools.set_group_spec_pool_type(0, 'file_name')
    for bc in BARCODES:
        pools.file_name_pool.add(bc)
    return pools


def all_grouped_counters(tmpdir, pools):
    """One instance of every counter class that resolves group names."""
    def path(name):
        return os.path.join(tmpdir, name)

    counters = {
        "AssignedFeatureCounter/gene": create_gene_counter(
            path("gene"), "with_ambiguous", string_pools=pools, group_index=0),
        "AssignedFeatureCounter/transcript": create_transcript_counter(
            path("transcript"), "with_ambiguous", string_pools=pools, group_index=0),
        "ExonCounter": ExonCounter(path("exon"), string_pools=pools, group_index=0),
        "IntronCounter": IntronCounter(path("intron"), string_pools=pools, group_index=0),
        "IntronRetentionCounter": IntronRetentionCounter(
            path("ir"), string_pools=pools, group_index=0),
        "JointExonCounter": JointExonCounter(
            path("joint_exon"), string_pools=pools, group_index=0),
        "ExonSpliceSiteCounter": ExonSpliceSiteCounter(
            path("splice_site"), string_pools=pools, group_index=0),
    }
    return counters


@pytest.fixture
def counters():
    tmpdir = tempfile.mkdtemp()
    pools = make_pools()
    return all_grouped_counters(tmpdir, pools), pools


def test_every_counter_resolves_a_valid_id(counters):
    built, pools = counters
    gid = pools.file_name_pool.get_int(BARCODES[0])
    for name, counter in built.items():
        assert counter._get_group_name(gid) == BARCODES[0], name


def test_no_counter_maps_the_sentinel_onto_a_real_group(counters):
    """The reported bug: -1 indexed the pool list from the end."""
    built, _ = counters
    for name, counter in built.items():
        resolved = counter._get_group_name(UNASSIGNED_GROUP_ID)
        assert resolved == UNASSIGNED_GROUP_NAME, name
        assert resolved not in BARCODES, f"{name} mis-attributed to a real group"


def test_no_counter_raises_on_a_stale_id(counters):
    """A not-yet-populated pool must not abort a finished run."""
    built, _ = counters
    for name, counter in built.items():
        assert counter._get_group_name(99) == "99", name


def test_terminal_counter_agrees(tmp_path):
    # TerminalCounter skips AbstractCounter.__init__ but inherits the helpers.
    pools = make_pools()
    counter = PolyACounter(SimpleNamespace(), str(tmp_path / "polya.tsv"),
                           string_pools=pools, group_index=0)
    assert counter._get_group_name(pools.file_name_pool.get_int(BARCODES[1])) == BARCODES[1]
    assert counter._get_group_name(UNASSIGNED_GROUP_ID) == UNASSIGNED_GROUP_NAME
    assert counter._get_group_name(99) == "99"


def test_ungrouped_counter_falls_back_without_pools(tmp_path):
    counter = create_gene_counter(str(tmp_path / "gene"), "with_ambiguous")
    assert counter.string_pools is None
    assert counter._get_group_name(0) == UNASSIGNED_GROUP_NAME
    assert counter._get_ordered_groups() == [UNASSIGNED_GROUP_NAME]
    assert counter._get_num_groups() == 1


def test_enumeration_contract_is_unchanged(counters):
    """Pool walking is by position and must stay exact -- it is not resolution."""
    built, _ = counters
    for name, counter in built.items():
        assert counter._get_ordered_groups() == BARCODES, name
        assert counter._get_num_groups() == len(BARCODES), name


def test_helpers_are_defined_once():
    """Guard against a sixth copy drifting back in."""
    import inspect
    import isoquant_lib.quantification.long_read_counter as lrc
    import isoquant_lib.terminal_prediction.terminal_counter as tc

    for module in (lrc, tc):
        src = inspect.getsource(module)
        assert src.count("def _get_group_name") <= 1, module.__name__
        assert src.count("def _group_name") == 0, module.__name__
    assert "def _get_group_name" in inspect.getsource(lrc.AbstractCounter)
