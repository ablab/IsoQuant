############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""BAM outputs: per-chromosome tagging and merging."""

import os

import pysam
import pytest

from isoquant_lib.utils.bam_utils import (
    GENE_TAG,
    TRANSCRIPT_TAG,
    collect_unmapped_read_ids,
    load_barcode_umi_tags,
    load_survivor_tags,
    merge_bam_files,
    references_with_alignments,
    write_tagged_chromosome_bam,
    write_unmapped_bam,
)

HEADER = {"HD": {"VN": "1.0"},
          "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chr2", "LN": 10000}]}


def make_read(name, chr_index, start=100, secondary=False, supplementary=False, unmapped=False):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = "ACGT" * 10
    read.query_qualities = pysam.qualitystring_to_array("I" * 40)
    if unmapped:
        read.is_unmapped = True
        read.reference_id = -1
        read.reference_start = -1
    else:
        read.reference_id = chr_index
        read.reference_start = start
        read.mapping_quality = 60
        read.cigarstring = "40M"
        read.is_secondary = secondary
        read.is_supplementary = supplementary
    return read


def write_bam(path, reads):
    with pysam.AlignmentFile(path, "wb", header=HEADER) as outf:
        for read in reads:
            outf.write(read)
    pysam.index(path)
    return path


def read_names(path):
    with pysam.AlignmentFile(path, "rb") as inf:
        return [read.query_name for read in inf.fetch(until_eof=True)]


@pytest.fixture
def input_bam(tmp_path):
    """chr1 has a primary, a secondary and a supplementary record for r1, plus r2 and r3."""
    reads = [
        make_read("r1", 0, 100),
        make_read("r1", 0, 200, secondary=True),
        make_read("r1", 0, 300, supplementary=True),
        make_read("r2", 0, 400),
        make_read("r3", 0, 500),
        make_read("r4", 1, 100),
        make_read("unmapped1", 0, unmapped=True),
    ]
    return write_bam(str(tmp_path / "in.bam"), reads)


class TestWriteTaggedChromosomeBam:
    def test_tags_are_applied(self, input_bam, tmp_path):
        out = str(tmp_path / "out.bam")
        tags = {"r1": ("ACGTACGT", "TTTTTT", "GENE1", "TX1")}
        write_tagged_chromosome_bam([input_bam], out, "chr1", tags)
        with pysam.AlignmentFile(out, "rb") as inf:
            tagged = [r for r in inf.fetch(until_eof=True) if r.query_name == "r1"]
        assert tagged, "r1 must be present"
        for read in tagged:
            assert read.get_tag("CB") == "ACGTACGT"
            assert read.get_tag("UB") == "TTTTTT"
            assert read.get_tag(GENE_TAG) == "GENE1"
            assert read.get_tag(TRANSCRIPT_TAG) == "TX1"

    def test_custom_tag_names(self, input_bam, tmp_path):
        out = str(tmp_path / "out.bam")
        write_tagged_chromosome_bam([input_bam], out, "chr1", {"r2": ("AAAA", "CCCC", None, None)},
                                    barcode_tag="XC", umi_tag="XM")
        with pysam.AlignmentFile(out, "rb") as inf:
            read = next(r for r in inf.fetch(until_eof=True) if r.query_name == "r2")
        assert read.get_tag("XC") == "AAAA"
        assert read.get_tag("XM") == "CCCC"
        assert not read.has_tag("CB")

    @pytest.mark.parametrize("placeholder", ["*", ".", "None"])
    def test_placeholders_are_not_written(self, input_bam, tmp_path, placeholder):
        """'*' means an uncalled barcode, '.' an unset feature, and IsoQuant writes the
        literal 'None' for the transcript of a novel or ambiguous read."""
        out = str(tmp_path / "out.bam")
        write_tagged_chromosome_bam([input_bam], out, "chr1",
                                    {"r2": (placeholder,) * 4})
        with pysam.AlignmentFile(out, "rb") as inf:
            read = next(r for r in inf.fetch(until_eof=True) if r.query_name == "r2")
        assert not read.has_tag("CB")
        assert not read.has_tag("UB")
        assert not read.has_tag(GENE_TAG)
        assert not read.has_tag(TRANSCRIPT_TAG)

    def test_keeps_untagged_reads_by_default(self, input_bam, tmp_path):
        out = str(tmp_path / "out.bam")
        written = write_tagged_chromosome_bam([input_bam], out, "chr1", {"r1": ("A", "C", "G", "T")})
        assert written == 5  # every chr1 alignment, tagged or not
        assert sorted(set(read_names(out))) == ["r1", "r2", "r3"]

    def test_subset_when_untagged_are_dropped(self, input_bam, tmp_path):
        out = str(tmp_path / "out.bam")
        written = write_tagged_chromosome_bam([input_bam], out, "chr1", {"r2": ("A", "C", "G", "T")},
                                              keep_untagged=False)
        assert written == 1
        assert read_names(out) == ["r2"]

    def test_primary_only_drops_secondary_and_supplementary(self, input_bam, tmp_path):
        out = str(tmp_path / "out.bam")
        written = write_tagged_chromosome_bam([input_bam], out, "chr1", {"r1": ("A", "C", "G", "T")},
                                              primary_only=True, keep_untagged=False)
        assert written == 1
        with pysam.AlignmentFile(out, "rb") as inf:
            read = next(inf.fetch(until_eof=True))
        assert not read.is_secondary and not read.is_supplementary

    def test_other_chromosomes_are_untouched(self, input_bam, tmp_path):
        out = str(tmp_path / "out.bam")
        write_tagged_chromosome_bam([input_bam], out, "chr2", {})
        assert read_names(out) == ["r4"]

    def test_records_are_copied_verbatim(self, input_bam, tmp_path):
        """Fusion detection reads SA tags and soft clips, so nothing but tags may change."""
        source = pysam.AlignedSegment()
        source.query_name = "sa_read"
        source.query_sequence = "ACGT" * 10
        source.query_qualities = pysam.qualitystring_to_array("I" * 40)
        source.reference_id = 0
        source.reference_start = 700
        source.mapping_quality = 60
        source.cigarstring = "10S30M"
        source.set_tag("SA", "chr2,100,+,30M10S,60,0;")
        bam = write_bam(str(tmp_path / "sa.bam"), [source])

        out = str(tmp_path / "out.bam")
        write_tagged_chromosome_bam([bam], out, "chr1", {"sa_read": ("A", "C", "G", "T")},
                                    primary_only=True, keep_untagged=False)
        with pysam.AlignmentFile(out, "rb") as inf:
            read = next(inf.fetch(until_eof=True))
        assert read.get_tag("SA") == "chr2,100,+,30M10S,60,0;"
        assert read.cigarstring == "10S30M"
        assert read.query_sequence == "ACGT" * 10

    def test_several_input_bams(self, tmp_path):
        first = write_bam(str(tmp_path / "a.bam"), [make_read("r1", 0)])
        second = write_bam(str(tmp_path / "b.bam"), [make_read("r2", 0)])
        out = str(tmp_path / "out.bam")
        assert write_tagged_chromosome_bam([first, second], out, "chr1", {}) == 2
        assert sorted(read_names(out)) == ["r1", "r2"]


class TestReferencesWithAlignments:
    def test_only_references_that_carry_reads(self, input_bam, tmp_path):
        """chr2 has one read, so both chromosomes are listed; an empty one would not be."""
        assert references_with_alignments([input_bam]) == ["chr1", "chr2"]

    def test_empty_references_are_dropped(self, tmp_path):
        """A fragmented assembly must not spawn a task per empty contig."""
        bam = write_bam(str(tmp_path / "one.bam"), [make_read("r1", 0)])
        assert references_with_alignments([bam]) == ["chr1"]

    def test_union_across_bams(self, tmp_path):
        first = write_bam(str(tmp_path / "a.bam"), [make_read("r1", 0)])
        second = write_bam(str(tmp_path / "b.bam"), [make_read("r4", 1)])
        assert references_with_alignments([first, second]) == ["chr1", "chr2"]

    def test_covers_references_isoquant_would_skip(self, input_bam):
        """The regression: scaffolds absent from the analysed chromosome list still have reads."""
        assert "chr2" in references_with_alignments([input_bam])


class TestWriteUnmappedBam:
    def test_only_unmapped_reads(self, input_bam, tmp_path):
        """fetch() by chromosome never returns these, so they need their own pass."""
        out = str(tmp_path / "unmapped.bam")
        assert write_unmapped_bam([input_bam], out) == 1
        assert read_names(out) == ["unmapped1"]

    def test_unmapped_reads_get_tagged(self, input_bam, tmp_path):
        """A barcode is called from the sequence, so an unaligned read still has one.

        These reads belong to no chromosome and appear in no split table; copying them
        untagged silently dropped the tags of every unmapped read.
        """
        out = str(tmp_path / "unmapped.bam")
        write_unmapped_bam([input_bam], out, {"unmapped1": ("ACGT", "TTTT", None, None)})
        with pysam.AlignmentFile(out, "rb") as inf:
            read = next(inf.fetch(until_eof=True))
        assert read.get_tag("CB") == "ACGT"
        assert read.get_tag("UB") == "TTTT"

    def test_untagged_unmapped_reads_are_still_copied(self, input_bam, tmp_path):
        out = str(tmp_path / "unmapped.bam")
        assert write_unmapped_bam([input_bam], out, {"someone_else": ("A", "C", None, None)}) == 1
        with pysam.AlignmentFile(out, "rb") as inf:
            assert not next(inf.fetch(until_eof=True)).has_tag("CB")


class TestCollectUnmappedReadIds:
    def test_finds_the_unplaced_reads(self, input_bam):
        assert collect_unmapped_read_ids([input_bam]) == {"unmapped1"}

    def test_none_when_all_aligned(self, tmp_path):
        bam = write_bam(str(tmp_path / "a.bam"), [make_read("r1", 0)])
        assert collect_unmapped_read_ids([bam]) == set()


class TestMergeBamFiles:
    def test_merges_and_indexes(self, tmp_path):
        first = write_bam(str(tmp_path / "a.bam"), [make_read("r1", 0)])
        second = write_bam(str(tmp_path / "b.bam"), [make_read("r4", 1)])
        out = str(tmp_path / "merged.bam")
        assert merge_bam_files(out, [first, second]) == out
        assert sorted(read_names(out)) == ["r1", "r4"]
        assert os.path.exists(out + ".bai") or os.path.exists(out + ".csi")
        assert not os.path.exists(first) and not os.path.exists(second)

    def test_single_fragment_is_moved(self, tmp_path):
        only = write_bam(str(tmp_path / "a.bam"), [make_read("r1", 0)])
        out = str(tmp_path / "merged.bam")
        assert merge_bam_files(out, [only]) == out
        assert not os.path.exists(only)
        assert read_names(out) == ["r1"]

    def test_missing_and_none_fragments_are_skipped(self, tmp_path):
        """Chromosomes with no surviving reads contribute None."""
        only = write_bam(str(tmp_path / "a.bam"), [make_read("r1", 0)])
        out = str(tmp_path / "merged.bam")
        assert merge_bam_files(out, [None, only, str(tmp_path / "gone.bam")]) == out
        assert read_names(out) == ["r1"]

    def test_nothing_to_merge(self, tmp_path):
        out = str(tmp_path / "merged.bam")
        assert merge_bam_files(out, [None, None]) is None
        assert not os.path.exists(out)


class TestLoadSurvivorTags:
    def test_with_tag_columns(self, tmp_path):
        path = tmp_path / "s.tsv"
        path.write_text("r1\tACGT\tTTTT\tGENE1\tTX1\nr2\tCCCC\tGGGG\t.\t.\n")
        assert load_survivor_tags(str(path)) == {
            "r1": ("ACGT", "TTTT", "GENE1", "TX1"),
            "r2": ("CCCC", "GGGG", ".", "."),
        }

    def test_read_ids_only(self, tmp_path):
        """Without deduplicated_bam the table is a bare id list; the ids still define the subset."""
        path = tmp_path / "s.tsv"
        path.write_text("r1\nr2\n")
        assert load_survivor_tags(str(path)) == {"r1": (None, None, None, None),
                                                 "r2": (None, None, None, None)}

    def test_missing_file(self, tmp_path):
        assert load_survivor_tags(str(tmp_path / "nope.tsv")) == {}


class TestLoadBarcodeUmiTags:
    def test_reads_the_split_table(self, tmp_path):
        path = tmp_path / "b.tsv"
        path.write_text("#read_id\tbarcode\tUMI\nr1\tACGT\tTTTT\nr2\tCCCC\tGGGG\n")
        assert load_barcode_umi_tags(str(path)) == {"r1": ("ACGT", "TTTT", None, None),
                                                    "r2": ("CCCC", "GGGG", None, None)}

    def test_filtering_by_read_ids(self, tmp_path):
        """The whole-sample table is scanned for a handful of unmapped reads."""
        path = tmp_path / "b.tsv"
        path.write_text("r1\tACGT\tTTTT\nr2\tCCCC\tGGGG\nr3\tAAAA\tCCCC\n")
        assert load_barcode_umi_tags(str(path), read_ids={"r2"}) == {
            "r2": ("CCCC", "GGGG", None, None)}
        assert load_barcode_umi_tags(str(path), read_ids=set()) == {}

    def test_short_lines_are_skipped(self, tmp_path):
        path = tmp_path / "b.tsv"
        path.write_text("r1\tACGT\tTTTT\nbroken\n")
        assert list(load_barcode_umi_tags(str(path))) == ["r1"]

    def test_missing_file(self, tmp_path):
        assert load_barcode_umi_tags(str(tmp_path / "nope.tsv")) == {}
