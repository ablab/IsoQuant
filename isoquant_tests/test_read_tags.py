############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import gzip

import pysam
import pytest

from isoquant_lib.utils.input_data_storage import InputDataType
from isoquant_lib.utils.read_mapper import read_input_stages
from isoquant_lib.utils.read_tags import (
    bam_tags_to_keep,
    count_tagged_alignments,
    fastx_tags_to_keep,
    is_fasta,
    split_header,
    tags_from_comment,
    tags_to_keep,
)

DORADO_COMMENT = "qs:f:12.5\tdu:f:1.2\tns:i:4000\tts:i:7\tch:i:12\tst:Z:2026-01-01T00:00:00.000+00:00" \
                 "\tRG:Z:run1_model\tpt:i:85"


class TestTagsFromComment:
    def test_dorado_header(self):
        assert tags_from_comment(DORADO_COMMENT) == {"qs", "du", "ns", "ts", "ch", "st", "RG", "pt"}

    def test_empty_comment(self):
        assert tags_from_comment("") == set()

    @pytest.mark.parametrize("comment", [
        "runid=abc read=12 ch=5 start_time=2026-01-01T00:00:00Z",
        "1:N:0:ACGT",
        "pt:i:85 qs:f:12.5",
        "pt:i:abc",
        "pt:A:xy",
        "pt:i:85\tsome text",
    ])
    def test_not_sam_tags(self, comment):
        assert tags_from_comment(comment) is None

    @pytest.mark.parametrize("field", ["pt:i:-1", "qs:f:1e-3", "TS:A:+", "Zz:Z:with space", "hx:H:1AE3",
                                       "ar:B:c,1,-2", "MM:Z:C+m?,1,2;"])
    def test_value_types(self, field):
        assert tags_from_comment(field) == {field[:2]}


def test_split_header():
    assert split_header("@read1\tpt:i:5\n") == ("read1", "pt:i:5")
    assert split_header(">read1 pt:i:5\tqs:f:1\n") == ("read1", "pt:i:5\tqs:f:1")
    assert split_header("@read1\n") == ("read1", "")


def test_tags_to_keep_drops_minimap2_and_header_tags():
    assert tags_to_keep({"pt", "TS", "ts", "NM", "RG", "PG", "qs"}) == ["TS", "pt", "qs"]


def _write_fastq(path, comments, compress=False):
    opener = gzip.open if compress else open
    with opener(path, "wt") as out:
        for i, comment in enumerate(comments):
            header = "@read%d" % i + ("\t" + comment if comment else "")
            out.write("%s\nACGTACGT\n+\nIIIIIIII\n" % header)


class TestFastxTags:
    def test_dorado_fastq(self, tmp_path):
        path = str(tmp_path / "reads.fq")
        _write_fastq(path, [DORADO_COMMENT] * 3)
        assert fastx_tags_to_keep(path) == ["ch", "du", "ns", "pt", "qs", "st"]
        assert not is_fasta(path)

    def test_gzipped_without_gz_suffix(self, tmp_path):
        path = str(tmp_path / "reads.fastq")
        _write_fastq(path, ["pt:i:1"], compress=True)
        assert fastx_tags_to_keep(path) == ["pt"]

    def test_no_comments(self, tmp_path):
        path = str(tmp_path / "reads.fq")
        _write_fastq(path, ["", ""])
        assert fastx_tags_to_keep(path) == []

    def test_one_invalid_comment_disables_tags(self, tmp_path):
        path = str(tmp_path / "reads.fq")
        _write_fastq(path, ["pt:i:1", "runid=abc read=1"])
        assert fastx_tags_to_keep(path) == []

    def test_only_first_reads_are_checked(self, tmp_path):
        path = str(tmp_path / "reads.fq")
        _write_fastq(path, ["pt:i:1", "pt:i:2", "runid=abc"])
        assert fastx_tags_to_keep(path, sample_size=2) == ["pt"]

    def test_quality_line_starting_with_at(self, tmp_path):
        path = str(tmp_path / "reads.fq")
        with open(path, "w") as out:
            out.write("@read0\tpt:i:1\nACGT\n+\n@III\n@read1\tqs:f:1\nACGT\n+\nIIII\n")
        assert fastx_tags_to_keep(path) == ["pt", "qs"]

    def test_fasta(self, tmp_path):
        path = str(tmp_path / "reads.fa")
        with open(path, "w") as out:
            out.write(">read0\tpt:i:1\nACGT\nACGT\n>read1\tTS:A:+\nACGT\n")
        assert fastx_tags_to_keep(path) == ["TS", "pt"]
        assert is_fasta(path)


def _write_unaligned_bam(path, tags_per_read):
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "unknown"}})
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for i, tags in enumerate(tags_per_read):
            read = pysam.AlignedSegment(header)
            read.query_name = "read%d" % i
            read.query_sequence = "ACGTACGT"
            read.flag = 4
            read.set_tags(tags)
            bam.write(read)


def test_bam_tags_to_keep(tmp_path):
    path = str(tmp_path / "calls.bam")
    _write_unaligned_bam(path, [[("pt", 85), ("ts", 7), ("RG", "run1")], [("qs", 12.5), ("TS", "+")]])
    assert bam_tags_to_keep(path) == ["TS", "pt", "qs"]


def test_count_tagged_alignments(tmp_path):
    path = str(tmp_path / "aligned.bam")
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": "chr1", "LN": 1000}]})
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for i, (flag, tags) in enumerate([(0, [("pt", 10)]), (256, [("pt", 10)]), (0, []), (2048, [])]):
            read = pysam.AlignedSegment(header)
            read.query_name = "read%d" % i
            read.query_sequence = "ACGTACGT"
            read.flag = flag
            read.reference_id = 0
            read.reference_start = 10 * i
            read.cigarstring = "8M"
            read.set_tags(tags)
            bam.write(read)
    assert count_tagged_alignments(path, "pt") == (1, 2)
    assert count_tagged_alignments(path, "XX") == (0, 2)


class TestReadInputStages:
    def test_fastq_without_tags_is_read_directly(self, tmp_path):
        path = str(tmp_path / "reads.fq")
        _write_fastq(path, [""])
        assert read_input_stages(path, InputDataType.fastq) == (path, [], [])

    def test_fastq_with_tags_goes_through_samtools(self, tmp_path):
        path = str(tmp_path / "reads.fq")
        _write_fastq(path, [DORADO_COMMENT])
        reads_input, stages, kept_tags = read_input_stages(path, InputDataType.fastq)
        assert reads_input == "-"
        assert "pt" in kept_tags and "ts" not in kept_tags and "RG" not in kept_tags
        tag_list = ",".join(kept_tags)
        assert [stage[0] for stage in stages] == [["samtools", "import", "-T", tag_list, path],
                                                  ["samtools", "fastq", "-T", tag_list, "-"]]

    def test_fasta_with_tags_stays_fasta(self, tmp_path):
        path = str(tmp_path / "reads.fa")
        with open(path, "w") as out:
            out.write(">read0\tpt:i:1\nACGT\n")
        _, stages, _ = read_input_stages(path, InputDataType.fastq)
        assert stages[1][0][:2] == ["samtools", "fasta"]

    def test_unmapped_bam(self, tmp_path):
        path = str(tmp_path / "calls.bam")
        _write_unaligned_bam(path, [[("pt", 85), ("RG", "run1")]])
        reads_input, stages, kept_tags = read_input_stages(path, InputDataType.unmapped_bam)
        assert reads_input == "-"
        assert kept_tags == ["pt"]
        assert [stage[0] for stage in stages] == [["samtools", "fastq", "-T", "pt", path]]

    def test_unmapped_bam_without_tags(self, tmp_path):
        path = str(tmp_path / "calls.bam")
        _write_unaligned_bam(path, [[]])
        _, stages, kept_tags = read_input_stages(path, InputDataType.unmapped_bam)
        assert kept_tags == []
        assert [stage[0] for stage in stages] == [["samtools", "fastq", path]]
