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
    gene_db = gffutils_db["ENSMUSG00000020196.10"]
    gene_info = GeneInfo([gene_db], gffutils_db)

    def counter(self):
        return PolyACounter(output_prefix="out.tsv", read_groups=["RG1", "RG2"])
    
    def read_assignment(self, group = ['RG1']):
        read = ReadAssignment("id1", ReadAssignmentType.unique, IsoformMatch(MatchClassification.full_splice_match, assigned_transcript="id2"))
        read.strand = "-"
        read.chr_id = "chr10"
        read.gene_info = self.gene_info
        read.polyA_found = True 
        read.assignment_type = ReadAssignmentType.inconsistent_non_intronic
        read.isoform_matches[0].assigned_gene = "ENSMUSG00000020196.10"
        read.isoform_matches[0].assigned_transcript = "ENSMUST00000001712.7"
        read.polya_info = pd.Series({'external_polyt_pos': 100})
        read.read_group = group[0]
        return read


    def test_counts(self):
        x = pd.Series({
            "transcript_id": "id1",
            "peak_location": 6,
            "histogram": [0]*10 + [0, 1, 2, 3, 10, 20, 30, 5, 1, 0] + [0]*10,
            "peak_left": 1,
            "peak_right": 8
        })
        
        counter = self.counter()
        result = counter.counts(x)

        assert isinstance(result, (int, np.integer))
        assert result == 72


    def test_counts_byGroup(self):
        
        counter = self.counter()
        counter.transcripts["id1"] = {
            "data": [1, 2, 3, 4, 5, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8],
            0: [1, 2, 3, 4, 6, 7, 7, 8, 8],
            1: [7, 7, 8, 5, 7, 6]
        }

        x = pd.Series({"transcript_id": "id1", "peak_location": 6,
            "peak_left": 0,
            "peak_right": 7})
        counts, heights = counter.counts_byGroup(x)

        assert isinstance(counts, list)
        assert isinstance(heights, list)
        assert all(isinstance(c, int) for c in counts)
        assert all(isinstance(h, int) for h in heights)
        assert counts == [9, 6]
        assert heights == [2, 3]

    def test_sort_peaks(self):
        
        counter = self.counter()
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
        assert np.array_equal(sorted_row["peak_prominence"], np.array([0.9, 0.5, 0.2]))
        assert np.array_equal(sorted_row["peak_width"], np.array([2, 3, 1]))
        assert np.array_equal(sorted_row["peak_heights"], np.array([50, 30, 10]))
        assert sorted_row["peak_left"] == [19, 30, 10]
        assert sorted_row["peak_right"] == [21, 31, 10]
        assert sorted_row["rank"] == [1, 2, 3]



    def test_flag(self):
        
        counter = self.counter()

        row = pd.Series({"prediction": 100, "annotated": 105})
        assert counter.flag(row) == "Known"
        row = pd.Series({"prediction": 100, "annotated": 200})
        assert counter.flag(row) == "Novel"


    def test_add_read_info(self):
        ra = self.read_assignment()
        counter = self.counter()
        counter.add_read_info(ra)

        assert ra.isoform_matches[0].assigned_transcript in counter.transcripts
        data = counter.transcripts[ra.isoform_matches[0].assigned_transcript]
        assert data["chr"] == "chr10"
        assert data["gene_id"] == "ENSMUSG00000020196.10"
        assert 100 in data["data"]
        annot = ra.gene_info.all_isoforms_exons[ra.isoform_matches[0].assigned_transcript][0][0] - 1
        assert data['annotated'] == annot


    def test_dump(self):
        counter = self.counter()
        transcripts = {'ENSMUST00000001712.7':  {'chr': "chr10", 'gene_id': "ENSMUSG00000020196.10", 'data': [10, 12, 12, 15, 15 ,15, 15, 15, 16, 19, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21], 0: [10, 12, 15, 15 ,15, 16, 21], 1:[ 12, 15, 15, 19, 20, 21, 21, 21, 21, 21, 21, 21, 21], 'annotated': 21}}

        counter.transcripts = transcripts   
    
        counter.dump()
        df = counter.dfResult 

        pd.set_option('display.max_columns', None)
        assert df.loc[0]['chromosome'] == 'chr10'
        assert df.loc[0]['transcript_id'] == 'ENSMUST00000001712.7'
        assert df.loc[0]['gene_id'] == 'ENSMUSG00000020196.10'
        for i in df['prediction']:
            assert (i <= np.max(transcripts['ENSMUST00000001712.7']['data']) and i >= np.min(transcripts['ENSMUST00000001712.7']['data']))
        assert df.loc[0]['counts'] == 11
        assert df.loc[0]['counts_byGroup'] == 1
        assert df.loc[1]['counts_byGroup'] == 10
        assert df.loc[0]['flag'] == 'Known'
        assert df.loc[0]['group_id'] == 'RG1'
        assert df.loc[1]['group_id'] == 'RG2'

