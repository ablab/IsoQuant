############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Compression helpers for large single-cell outputs."""

import gzip
import os

import pytest

from isoquant_lib.barcode_calling.detect_barcodes import numbered_chunk_name
from isoquant_lib.utils.file_utils import (
    gzip_file_in_place,
    open_text_read,
    open_text_write,
    resolve_optionally_gzipped,
    strip_compression_suffix,
)


class TestOpenHelpers:
    def test_plain_round_trip(self, tmp_path):
        path = str(tmp_path / "x.tsv")
        with open_text_write(path) as handle:
            handle.write("a\nb\n")
        assert not os.path.exists(path + ".gz")
        with open_text_read(path) as handle:
            assert handle.read() == "a\nb\n"

    def test_gzip_round_trip(self, tmp_path):
        path = str(tmp_path / "x.tsv.gz")
        with open_text_write(path) as handle:
            handle.write("a\nb\n")
        # really compressed, not just named .gz
        with gzip.open(path, "rt") as handle:
            assert handle.read() == "a\nb\n"
        with open_text_read(path) as handle:
            assert handle.read() == "a\nb\n"


class TestStripCompressionSuffix:
    @pytest.mark.parametrize("name, expected", [
        ("S.split_reads_0.fa.gz", "S.split_reads_0.fa"),
        ("S.split_reads_0.fa", "S.split_reads_0.fa"),
        ("reads.fastq.gzip", "reads.fastq"),
        ("reads.fastq.bgz", "reads.fastq"),
        ("reads.bam", "reads.bam"),
    ])
    def test_suffixes(self, name, expected):
        assert strip_compression_suffix(name) == expected

    def test_read_group_name_is_unaffected_by_compression(self):
        """File-name read groups must not change when an input gets compressed."""
        plain = os.path.splitext(strip_compression_suffix("S.split_reads_0.fa"))[0]
        gzipped = os.path.splitext(strip_compression_suffix("S.split_reads_0.fa.gz"))[0]
        assert plain == gzipped == "S.split_reads_0"


class TestResolveOptionallyGzipped:
    def test_prefers_the_plain_file(self, tmp_path):
        plain = tmp_path / "t.tsv"
        plain.write_text("x")
        (tmp_path / "t.tsv.gz").write_bytes(b"")
        assert resolve_optionally_gzipped(str(plain)) == str(plain)

    def test_falls_back_to_gz(self, tmp_path):
        """A resumed run refers to the plain name after the table was compressed."""
        gzipped = tmp_path / "t.tsv.gz"
        gzipped.write_bytes(b"")
        assert resolve_optionally_gzipped(str(tmp_path / "t.tsv")) == str(gzipped)

    def test_missing_returns_the_original(self, tmp_path):
        missing = str(tmp_path / "nope.tsv")
        assert resolve_optionally_gzipped(missing) == missing


class TestGzipFileInPlace:
    def test_compresses_and_removes_the_original(self, tmp_path):
        path = tmp_path / "t.tsv"
        path.write_text("read1\tACGT\tTTTT\n")
        result = gzip_file_in_place(str(path))
        assert result == str(path) + ".gz"
        assert not path.exists()
        with gzip.open(result, "rt") as handle:
            assert handle.read() == "read1\tACGT\tTTTT\n"

    def test_keep_original(self, tmp_path):
        path = tmp_path / "t.tsv"
        path.write_text("x\n")
        gzip_file_in_place(str(path), keep_original=True)
        assert path.exists()

    def test_already_compressed_is_a_no_op(self, tmp_path):
        path = tmp_path / "t.tsv.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("x\n")
        assert gzip_file_in_place(str(path)) == str(path)

    def test_missing_file_is_a_no_op(self, tmp_path):
        missing = str(tmp_path / "nope.tsv")
        assert gzip_file_in_place(missing) == missing
        assert not os.path.exists(missing + ".gz")


class TestNumberedChunkName:
    def test_plain(self):
        assert numbered_chunk_name("/tmp/subreads", 3) == "/tmp/subreads_3"

    def test_keeps_the_compression_suffix_last(self):
        """Otherwise the per-chunk temp would not be recognised as compressed."""
        assert numbered_chunk_name("/tmp/subreads.gz", 3) == "/tmp/subreads_3.gz"


class TestMultiMemberConcatenation:
    def test_byte_concatenated_members_read_as_one_stream(self, tmp_path):
        """Chunks are compressed in the workers and merged with a byte copy.

        That only works because concatenated gzip members form a valid stream.
        """
        merged = tmp_path / "merged.fa.gz"
        with open(merged, "wb") as out:
            for chunk in range(3):
                part = tmp_path / ("part_%d.gz" % chunk)
                with open_text_write(str(part)) as handle:
                    handle.write(">read%d\nACGT\n" % chunk)
                with open(part, "rb") as inf:
                    out.write(inf.read())
        with gzip.open(merged, "rt") as handle:
            assert handle.read() == ">read0\nACGT\n>read1\nACGT\n>read2\nACGT\n"
