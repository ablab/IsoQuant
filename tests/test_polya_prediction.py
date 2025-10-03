import unittest
import pytest
import numpy as np
import pandas as pd
import os
import tempfile
import gffutils

from src.isoform_assignment import (
    ReadAssignmentType,
    ReadAssignment,
    IsoformMatch,
    MatchClassification
)
from src.gene_info import GeneInfo
from src.polya_position_model import PolyACounter



class testPolyACounter(unittest.TestCase):
    source_dir = os.path.dirname(os.path.realpath(__file__))
    gffutils_db = gffutils.FeatureDB(os.path.join(source_dir, 'toy_data/synth.db'), keep_order=True)
    gene_db = gffutils_db['ENSMUSG00000020196.10']
    gene_info = GeneInfo([gene_db], gffutils_db)


    def testAdd_read_info(self):
        pass
    def testDump(self):
        pass





    @pytest.fixture
    def counter(self):
        return PolyACounter(output_prefix="out.tsv", read_groups=["RG1", "RG2"])
    
    @pytest.fixture
    def counter_ignoreGroups():
        return PolyACounter(output_prefix="out.tsv", read_groups=["RG1"])
    

    @pytest.fixture
    def read_assignment(self, strand, polya_pos, group = ['RG1']):
        read = ReadAssignment("id1", ReadAssignmentType.unique, IsoformMatch(MatchClassification.full_splice_match, assigned_transcript="id2"))
        # read.strand = strand
        # read.polyA_found = True
        # read.polya_info.external_polyt_pos = polya_pos
        read.gene_info = self.gene_info
        read.read_group = group
        read.chr_id = 20



    def test_counts(self, counter):
        x = pd.Series({
            "transcript_id": "id1",
            "peak_location": 6,
            "histogram": [0]*10 + [0, 1, 2, 3, 10, 20, 30, 5, 1, 0] + [0]*10
        })
        result = counter.counts(x)

        assert isinstance(result, (int, np.integer))
        assert result == 62


    def test_counts_byGroup(self, counter):

        counter.transcripts["id1"] = {
            "data": [1, 2, 3, 4, 5, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8],
            0: [1, 2, 3, 4, 6, 7, 7, 8, 8],
            1: [7, 7, 8, 5, 7, 6]
        }

        x = pd.Series({"transcript_id": "id1", "peak_location": 6})
        counts, heights = counter.counts_byGroup(x)

        assert isinstance(counts, list)
        assert isinstance(heights, list)
        assert all(isinstance(c, int) for c in counts)
        assert all(isinstance(h, int) for h in heights)
        assert counts == [8, 6]
        assert heights == [2, 3]

    def test_sort_peaks(self, counter):
        row = pd.Series({
            "peak_count": 3,
            "peak_location": [10, 20, 30],
            "peak_prominence": np.array([0.2, 0.9, 0.5]),
            "peak_width": np.array([1, 2, 3]),
            "peak_heights": np.array([10, 50, 30]),
            "peak_left": [10, 19, 30],
            "peak_right": [10, 21, 31]
        })

        sorted_row = counter.sort_peaks(row)

        assert sorted_row["peak_count"] == 3
        assert sorted_row["peak_location"] == [20, 30, 10]
        assert sorted_row["peak_prominence"] == np.array([0.9, 0.5, 0.2])
        assert sorted_row["peak_width"] == np.array([2, 3, 1])
        assert sorted_row["peak_heights"] == np.array([50, 30, 10])
        assert sorted_row["peak_left"] == [19, 30, 10]
        assert sorted_row["peak_right"] == [21, 31, 10]
        assert sorted_row["rank"] == [1, 2, 3]



    def test_flag(self, counter):

        row = pd.Series({"prediction": 100, "annotated": 105})
        assert counter.flag(row) == "Known"
        row = pd.Series({"prediction": 100, "annotated": 200})
        assert counter.flag(row) == "Novel"





    def test_add_read_info_populates_transcripts(counter, read_assignment):
        ra = read_assignment(
            strand="+", polya_pos=150, transcript="tx1", gene="gene1"
        )
        counter.add_read_info(ra)

        assert "tx1" in counter.transcripts
        tx_data = counter.transcripts["tx1"]

        assert tx_data["chr"] == "chr1"
        assert tx_data["gene_id"] == "gene1"
        assert isinstance(tx_data["data"], list)
        assert 150 in tx_data["data"]
        # Group should also contain the polya_pos
        assert 150 in tx_data[0]  # group_numeric_ids["grp1"] == 0


    def test_dump_creates_output_file(counter, tmp_path, read_assignment):
        # add several reads to create a histogram
        for pos in [120, 121, 122, 150, 151, 151, 152]:
            ra = read_assignment("+", pos, "tx1", "gene1")
            counter.add_read_info(ra)

        counter.dump()
        out_file = tmp_path / "out.tsv"
        assert out_file.exists()

        df = pd.read_csv(out_file, sep="\t")
        # Must contain required columns
        for col in [
            "chromosome",
            "transcript_id",
            "gene_id",
            "peak_left",
            "peak_right",
            "prediction",
            "peak_heights",
            "counts",
            "flag",
        ]:
            assert col in df.columns

        # Ensure predictions are integers
        assert np.issubdtype(df["prediction"].dtype, np.integer)


    def test_dump_appends_on_second_call(counter, tmp_path, read_assignment):
        # add first batch
        for pos in [200, 201, 202]:
            ra = read_assignment("+", pos, "tx1", "gene1")
            counter.add_read_info(ra)
        counter.dump()

        # add second batch
        for pos in [300, 301, 302]:
            ra = read_assignment("+", pos, "tx2", "gene2")
            counter.add_read_info(ra)
        counter.dump()

        out_file = tmp_path / "out.tsv"
        df = pd.read_csv(out_file, sep="\t")

        # Should contain both tx1 and tx2 results
        assert "tx1" in df["transcript_id"].values
        assert "tx2" in df["transcript_id"].values
