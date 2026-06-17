############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# Regression test for duplicate exon/intron count rows.
#
# A long, high-coverage gene can be split at a coverage valley
# (AlignmentCollector.split_coverage_regions), and each sub-region rebuilds a
# fresh GeneInfo for every overlapping gene (get_gene_info_for_region). The
# per-chromosome ExonCounter/IntronCounter is shared across those builds, so the
# same physical exon/intron must resolve to the same counter key or it emits two
# rows with split counts. FeatureInfo.id is therefore the coordinate tuple
# (chr_id, start, end, strand); this test guards that property and the resulting
# row merge. See gene_info.py FeatureInfo and long_read_counter.ProfileFeatureCounter.

import os

import gffutils

from isoquant_lib.gene_info import GeneInfo, FeatureInfo
from isoquant_lib.long_read_counter import ExonCounter, IntronCounter

SOURCE_DIR = os.path.dirname(os.path.realpath(__file__))
GENE_ID = "ENSMUSG00000020196.10"


def _load_gene_db():
    db = gffutils.FeatureDB(os.path.join(SOURCE_DIR, "toy_data/synth.db"), keep_order=True)
    return db, db[GENE_ID]


def _parse_counts(path):
    """Return list of (coord_key, include, exclude) for data rows of a counts file."""
    rows = []
    with open(path) as f:
        header = f.readline()  # chr start end strand flags gene_ids group_id include exclude
        assert header.startswith("chr\tstart\tend\tstrand\tflags\tgene_ids")
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            coord_key = (cols[0], cols[1], cols[2], cols[3])
            include, exclude = int(cols[-2]), int(cols[-1])
            rows.append((coord_key, include, exclude))
    return rows


def test_feature_id_is_coordinate_tuple_stable_across_builds():
    db, gene_db = _load_gene_db()
    gi1 = GeneInfo([gene_db], db)
    gi2 = GeneInfo([gene_db], db)

    for pmap in (gi1.exon_property_map, gi1.intron_property_map):
        for feat in pmap:
            assert feat.id == (feat.chr_id, feat.start, feat.end, feat.strand)

    # Two independent builds of the same gene must produce value-equal ids for
    # matching features (this is what merges counts across coverage-valley
    # sub-region rebuilds). A monotonic-counter id would break this.
    for a, b in zip(gi1.exon_property_map, gi2.exon_property_map):
        assert a is not b              # genuinely separate objects/builds
        assert a.id == b.id            # ... but value-equal coordinate ids
    for a, b in zip(gi1.intron_property_map, gi2.intron_property_map):
        assert a is not b
        assert a.id == b.id


def _feed_two_builds(counter, gi1, gi2, strand, n_include, n_exclude):
    """Inclusion reads via build 1's map, exclusion reads via build 2's map."""
    pmap1 = gi1.exon_property_map if isinstance(counter, ExonCounter) else gi1.intron_property_map
    pmap2 = gi2.exon_property_map if isinstance(counter, ExonCounter) else gi2.intron_property_map
    include_profile = [1] * len(pmap1)
    exclude_profile = [-1] * len(pmap2)
    for _ in range(n_include):
        counter.add_read_info_from_profile(include_profile, strand, pmap1, group_id=0)
    for _ in range(n_exclude):
        counter.add_read_info_from_profile(exclude_profile, strand, pmap2, group_id=0)


def test_exon_counter_merges_duplicate_gene_builds(tmp_path):
    db, gene_db = _load_gene_db()
    gi1 = GeneInfo([gene_db], db)
    gi2 = GeneInfo([gene_db], db)
    n_features = len(gi1.exon_property_map)

    counter = ExonCounter(str(tmp_path / "exon"))
    _feed_two_builds(counter, gi1, gi2, strand="-", n_include=3, n_exclude=2)
    counter.dump()

    rows = _parse_counts(counter.output_counts_file_name)
    coord_keys = [r[0] for r in rows]

    # No duplicate coordinate rows: one row per physical exon, despite two builds.
    assert len(rows) == len(set(coord_keys))
    assert len(rows) == n_features
    # Inclusion (build 1) and exclusion (build 2) counts merged into one row each.
    for _, include, exclude in rows:
        assert include == 3
        assert exclude == 2


def test_intron_counter_merges_duplicate_gene_builds(tmp_path):
    db, gene_db = _load_gene_db()
    gi1 = GeneInfo([gene_db], db)
    gi2 = GeneInfo([gene_db], db)
    n_features = len(gi1.intron_property_map)

    counter = IntronCounter(str(tmp_path / "intron"))
    _feed_two_builds(counter, gi1, gi2, strand="-", n_include=4, n_exclude=1)
    counter.dump()

    rows = _parse_counts(counter.output_counts_file_name)
    coord_keys = [r[0] for r in rows]

    assert len(rows) == len(set(coord_keys))
    assert len(rows) == n_features
    for _, include, exclude in rows:
        assert include == 4
        assert exclude == 1
