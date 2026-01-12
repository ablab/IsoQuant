############################################################################
# Copyright (c) 2025 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest
from src.barcode_calling.kmer_indexer import KmerIndexer, Dict2BitKmerIndexer, ArrayKmerIndexer, Array2BitKmerIndexer
from src.barcode_calling.common import str_to_2bit


class TestKmerIndexer:
    """Test basic dictionary-based KmerIndexer."""

    def test_init(self):
        """Test indexer initialization."""
        barcodes = ["ACTG", "TGCA", "GGGG"]
        indexer = KmerIndexer(barcodes, kmer_size=3)

        assert len(indexer.seq_list) == 3
        assert indexer.k == 3
        assert not indexer.empty()

    def test_empty(self):
        """Test empty index detection."""
        indexer = KmerIndexer([], kmer_size=3)
        assert indexer.empty()

    def test_get_kmers(self):
        """Test k-mer generation."""
        indexer = KmerIndexer([], kmer_size=3)
        kmers = list(indexer._get_kmers("ACTGAC"))

        assert len(kmers) == 4
        assert kmers == ["ACT", "CTG", "TGA", "GAC"]

    def test_get_kmers_short_sequence(self):
        """Test k-mer generation for short sequence."""
        indexer = KmerIndexer([], kmer_size=5)
        kmers = list(indexer._get_kmers("ACT"))

        assert len(kmers) == 0

    def test_exact_match(self):
        """Test finding exact match."""
        barcodes = ["ACTGACTG", "TGCATGCA", "GGGGGGGG"]
        indexer = KmerIndexer(barcodes, kmer_size=4)

        results = indexer.get_occurrences("ACTGACTG")

        assert len(results) > 0
        assert results[0][0] == "ACTGACTG"
        assert results[0][1] > 0  # Shared k-mer count

    def test_similar_match(self):
        """Test finding similar sequence."""
        barcodes = ["ACTGACTG", "TGCATGCA"]
        indexer = KmerIndexer(barcodes, kmer_size=4)

        # Query with 1 mismatch
        results = indexer.get_occurrences("ACTGACTT")

        assert len(results) > 0
        # Should find ACTGACTG as closest match
        assert results[0][0] == "ACTGACTG"

    def test_min_kmers_filter(self):
        """Test minimum k-mer threshold."""
        barcodes = ["ACTGACTG", "TGCATGCA"]
        indexer = KmerIndexer(barcodes, kmer_size=4)

        # Very different sequence
        results = indexer.get_occurrences("CCCCCCCC", min_kmers=3)

        # Should have fewer results with higher threshold
        assert len(results) <= 2

    def test_max_hits(self):
        """Test maximum hits limit."""
        barcodes = ["ACTG", "ACTC", "ACTA", "ACTT"]
        indexer = KmerIndexer(barcodes, kmer_size=3)

        results = indexer.get_occurrences("ACTG", max_hits=2)

        assert len(results) <= 2

    def test_ignore_equal(self):
        """Test ignoring exact matches."""
        barcodes = ["ACTG", "TGCA"]
        indexer = KmerIndexer(barcodes, kmer_size=3)

        results = indexer.get_occurrences("ACTG", ignore_equal=True)

        # Should not return the exact match
        for seq, _, _ in results:
            assert seq != "ACTG"

    def test_append(self):
        """Test adding barcode dynamically."""
        barcodes = ["ACTG"]
        indexer = KmerIndexer(barcodes, kmer_size=3)

        indexer.append("TGCA")

        assert len(indexer.seq_list) == 2
        results = indexer.get_occurrences("TGCA")
        assert any(seq == "TGCA" for seq, _, _ in results)

    def test_result_format(self):
        """Test result tuple format."""
        barcodes = ["ACTGACTG"]
        indexer = KmerIndexer(barcodes, kmer_size=4)

        results = indexer.get_occurrences("ACTGACTG")

        assert len(results) > 0
        # Result should be (sequence, count, positions)
        seq, count, positions = results[0]
        assert isinstance(seq, str)
        assert isinstance(count, int)
        assert isinstance(positions, list)


class TestDict2BitKmerIndexer:
    """Test memory-efficient Dict2BitKmerIndexer with integer k-mer keys."""

    def test_init(self):
        """Test indexer initialization."""
        barcodes = ["ACTG", "TGCA", "GGGG"]
        indexer = Dict2BitKmerIndexer(barcodes, kmer_size=3)

        assert len(indexer.seq_list) == 3
        assert indexer.k == 3
        assert not indexer.empty()

    def test_empty(self):
        """Test empty index detection."""
        indexer = Dict2BitKmerIndexer([], kmer_size=3)
        assert indexer.empty()

    def test_binary_encoding(self):
        """Test nucleotide to binary conversion."""
        assert Dict2BitKmerIndexer.NUCL2BIN['A'] == 0
        assert Dict2BitKmerIndexer.NUCL2BIN['C'] == 1
        assert Dict2BitKmerIndexer.NUCL2BIN['G'] == 2
        assert Dict2BitKmerIndexer.NUCL2BIN['T'] == 3
        # Lowercase should work too
        assert Dict2BitKmerIndexer.NUCL2BIN['a'] == 0
        assert Dict2BitKmerIndexer.NUCL2BIN['t'] == 3

    def test_get_kmer_indexes(self):
        """Test k-mer index generation."""
        indexer = Dict2BitKmerIndexer([], kmer_size=2)
        # "AC" = 00 01 = 1
        # "CT" = 01 11 = 7
        kmer_idxs = list(indexer._get_kmer_indexes("ACT"))

        assert len(kmer_idxs) == 2
        assert kmer_idxs[0] == 1  # AC
        assert kmer_idxs[1] == 7  # CT

    def test_get_kmer_indexes_short_sequence(self):
        """Test k-mer generation for short sequence."""
        indexer = Dict2BitKmerIndexer([], kmer_size=5)
        kmer_idxs = list(indexer._get_kmer_indexes("ACT"))

        assert len(kmer_idxs) == 0

    def test_exact_match(self):
        """Test finding exact match."""
        barcodes = ["ACTGACTG", "TGCATGCA", "GGGGGGGG"]
        indexer = Dict2BitKmerIndexer(barcodes, kmer_size=4)

        results = indexer.get_occurrences("ACTGACTG")

        assert len(results) > 0
        assert results[0][0] == "ACTGACTG"
        assert results[0][1] > 0  # Shared k-mer count

    def test_similar_match(self):
        """Test finding similar sequence."""
        barcodes = ["ACTGACTG", "TGCATGCA"]
        indexer = Dict2BitKmerIndexer(barcodes, kmer_size=4)

        # Query with 1 mismatch
        results = indexer.get_occurrences("ACTGACTT")

        assert len(results) > 0
        assert results[0][0] == "ACTGACTG"

    def test_min_kmers_filter(self):
        """Test minimum k-mer threshold."""
        barcodes = ["ACTGACTG", "TGCATGCA"]
        indexer = Dict2BitKmerIndexer(barcodes, kmer_size=4)

        # Very different sequence
        results = indexer.get_occurrences("CCCCCCCC", min_kmers=3)

        assert len(results) <= 2

    def test_max_hits(self):
        """Test maximum hits limit."""
        barcodes = ["ACTG", "ACTC", "ACTA", "ACTT"]
        indexer = Dict2BitKmerIndexer(barcodes, kmer_size=3)

        results = indexer.get_occurrences("ACTG", max_hits=2)

        assert len(results) <= 2

    def test_ignore_equal(self):
        """Test ignoring exact matches."""
        barcodes = ["ACTG", "TGCA"]
        indexer = Dict2BitKmerIndexer(barcodes, kmer_size=3)

        results = indexer.get_occurrences("ACTG", ignore_equal=True)

        for seq, _, _ in results:
            assert seq != "ACTG"

    def test_append(self):
        """Test adding barcode dynamically."""
        barcodes = ["ACTG"]
        indexer = Dict2BitKmerIndexer(barcodes, kmer_size=3)

        indexer.append("TGCA")

        assert len(indexer.seq_list) == 2
        results = indexer.get_occurrences("TGCA")
        assert any(seq == "TGCA" for seq, _, _ in results)

    def test_consistency_with_kmerindexer(self):
        """Test that results match basic KmerIndexer."""
        barcodes = ["ACTGACTG", "TGCATGCA", "GGGGGGGG"]
        basic = KmerIndexer(barcodes, kmer_size=4)
        dict2bit = Dict2BitKmerIndexer(barcodes, kmer_size=4)

        query = "ACTGACTT"
        basic_results = basic.get_occurrences(query)
        dict2bit_results = dict2bit.get_occurrences(query)

        # Should find same sequences
        basic_seqs = set(seq for seq, _, _ in basic_results)
        dict2bit_seqs = set(seq for seq, _, _ in dict2bit_results)
        assert basic_seqs == dict2bit_seqs

    def test_consistency_with_arraykmerindexer(self):
        """Test that results match ArrayKmerIndexer."""
        barcodes = ["ACTGACTG", "TGCATGCA", "GGGGGGGG"]
        array = ArrayKmerIndexer(barcodes, kmer_size=4)
        dict2bit = Dict2BitKmerIndexer(barcodes, kmer_size=4)

        query = "ACTGACTT"
        array_results = array.get_occurrences(query)
        dict2bit_results = dict2bit.get_occurrences(query)

        # Should find same sequences
        array_seqs = set(seq for seq, _, _ in array_results)
        dict2bit_seqs = set(seq for seq, _, _ in dict2bit_results)
        assert array_seqs == dict2bit_seqs

    def test_result_format(self):
        """Test result tuple format."""
        barcodes = ["ACTGACTG"]
        indexer = Dict2BitKmerIndexer(barcodes, kmer_size=4)

        results = indexer.get_occurrences("ACTGACTG")

        assert len(results) > 0
        seq, count, positions = results[0]
        assert isinstance(seq, str)
        assert isinstance(count, int)
        assert isinstance(positions, list)

    def test_index_uses_int_keys(self):
        """Test that index uses integer keys, not strings."""
        barcodes = ["ACTG"]
        indexer = Dict2BitKmerIndexer(barcodes, kmer_size=3)

        # All keys should be integers
        for key in indexer.index.keys():
            assert isinstance(key, int)


class TestArrayKmerIndexer:
    """Test optimized array-based KmerIndexer."""

    def test_init(self):
        """Test indexer initialization."""
        barcodes = ["ACTG", "TGCA"]
        indexer = ArrayKmerIndexer(barcodes, kmer_size=3)

        assert len(indexer.seq_list) == 2
        assert indexer.k == 3
        assert len(indexer.index) == 64  # 4^3

    def test_binary_encoding(self):
        """Test nucleotide to binary conversion."""
        assert ArrayKmerIndexer.NUCL2BIN['A'] == 0
        assert ArrayKmerIndexer.NUCL2BIN['C'] == 1
        assert ArrayKmerIndexer.NUCL2BIN['G'] == 2
        assert ArrayKmerIndexer.NUCL2BIN['T'] == 3

    def test_get_kmer_indexes(self):
        """Test k-mer index generation."""
        indexer = ArrayKmerIndexer([], kmer_size=2)
        # "AC" = 00 01 = 1
        # "CT" = 01 11 = 7
        kmer_idxs = list(indexer._get_kmer_indexes("ACT"))

        assert len(kmer_idxs) == 2
        assert kmer_idxs[0] == 1  # AC
        assert kmer_idxs[1] == 7  # CT

    def test_exact_match(self):
        """Test finding exact match."""
        barcodes = ["ACTGACTG", "TGCATGCA"]
        indexer = ArrayKmerIndexer(barcodes, kmer_size=4)

        results = indexer.get_occurrences("ACTGACTG")

        assert len(results) > 0
        assert results[0][0] == "ACTGACTG"

    def test_consistency_with_basic(self):
        """Test that results match basic KmerIndexer."""
        barcodes = ["ACTGACTG", "TGCATGCA", "GGGGGGGG"]
        basic = KmerIndexer(barcodes, kmer_size=4)
        array = ArrayKmerIndexer(barcodes, kmer_size=4)

        query = "ACTGACTT"
        basic_results = basic.get_occurrences(query)
        array_results = array.get_occurrences(query)

        # Should find same sequences (order/counts might differ slightly)
        basic_seqs = set(seq for seq, _, _ in basic_results)
        array_seqs = set(seq for seq, _, _ in array_results)
        assert basic_seqs == array_seqs


class TestArray2BitKmerIndexer:
    """Test memory-efficient 2-bit KmerIndexer."""

    def test_init(self):
        """Test indexer initialization with 2-bit sequences."""
        barcodes = ["ACTGACTGACTGACTGACTGACTGA"]  # 25-mer
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = Array2BitKmerIndexer(bin_seqs, kmer_size=6, seq_len=25)

        assert indexer.total_sequences == 1
        assert indexer.seq_len == 25
        assert indexer.k == 6
        assert not indexer.empty()

    def test_empty(self):
        """Test empty index detection."""
        indexer = Array2BitKmerIndexer([], kmer_size=6, seq_len=25)
        assert indexer.empty()

    def test_get_kmer_bin_indexes(self):
        """Test k-mer extraction from 2-bit sequence."""
        # Create simple sequence
        bin_seq = str_to_2bit("ACTGAA")
        indexer = Array2BitKmerIndexer([], kmer_size=3, seq_len=6)

        kmer_idxs = list(indexer._get_kmer_bin_indexes(bin_seq))

        assert len(kmer_idxs) == 4  # 6 - 3 + 1

    def test_exact_match(self):
        """Test finding exact match."""
        seq = "ACTGACTGACTGACTGACTGACTGA"
        bin_seqs = [str_to_2bit(seq)]
        indexer = Array2BitKmerIndexer(bin_seqs, kmer_size=8, seq_len=25)

        results = indexer.get_occurrences(seq)

        assert len(results) > 0
        assert results[0][0] == seq

    def test_similar_match(self):
        """Test finding similar sequence."""
        barcode1 = "ACTGACTGACTGACTGACTGACTGA"
        barcode2 = "TGCATGCATGCATGCATGCATGCAT"
        bin_seqs = [str_to_2bit(b) for b in [barcode1, barcode2]]
        indexer = Array2BitKmerIndexer(bin_seqs, kmer_size=8, seq_len=25)

        # Query with 1 mismatch from barcode1
        query = "ACTGACTGACTGACTGACTGACTGT"
        results = indexer.get_occurrences(query)

        assert len(results) > 0
        # Should find barcode1 as closest
        assert results[0][0] == barcode1

    def test_max_hits(self):
        """Test maximum hits limit."""
        barcodes = [
            "ACTGACTGACTGACTGACTGACTGA",
            "ACTGACTGACTGACTGACTGACTGC",
            "ACTGACTGACTGACTGACTGACTGG"
        ]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = Array2BitKmerIndexer(bin_seqs, kmer_size=8, seq_len=25)

        results = indexer.get_occurrences("ACTGACTGACTGACTGACTGACTGA", max_hits=2)

        assert len(results) <= 2

    def test_flat_index_structure(self):
        """Test that index uses flat array structure."""
        barcodes = ["ACTGACTGACTGACTGACTGACTGA"]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = Array2BitKmerIndexer(bin_seqs, kmer_size=6, seq_len=25)

        # Index should be flat list
        assert isinstance(indexer.index, list)
        # Index ranges should map to flat index
        assert isinstance(indexer.index_ranges, list)
        assert len(indexer.index_ranges) == 4**6 + 1  # 4^k + 1


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
