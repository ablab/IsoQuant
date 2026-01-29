############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import pytest
import gffutils
from src.gene_info import *


class TestGeneInfo:
    source_dir = os.path.dirname(os.path.realpath(__file__))
    gffutils_db = gffutils.FeatureDB(os.path.join(source_dir, 'toy_data/synth.db'), keep_order=True)
    gene_db = gffutils_db['ENSMUSG00000020196.10']

    def test_basic(self):
        gene_info = GeneInfo([self.gene_db], self.gffutils_db)
        assert gene_info.get_gene_region() == ("chr10", 1000, 10000)
        assert gene_info.transcript_start("ENSMUST00000001712.7") == 1000
        assert gene_info.transcript_end("ENSMUST00000001713.7") == 10000
        assert gene_info.transcript_region("ENSMUST00000001715.7") == (8000, 8800)
        assert gene_info.transcript_exon_count("ENSMUST00000001713.7") == 6
        assert gene_info.total_transcript_length("ENSMUST00000001712.7") == 1105
        assert gene_info.chr_id == "chr10"
        assert gene_info.isoform_strands["ENSMUST00000001713.7"] == '-'
        assert gene_info.gene_id_map["ENSMUST00000001713.7"] == 'ENSMUSG00000020196.10'

    def test_intron_profiles(self):
        gene_info = GeneInfo([self.gene_db], self.gffutils_db)
        assert gene_info.intron_profiles.features[0] == (1101, 1999)
        assert gene_info.intron_profiles.features[-1] == (8201, 8499)
        assert len(gene_info.intron_profiles.features) == 9
        assert gene_info.intron_profiles.profiles["ENSMUST00000001712.7"] == [1, 1, -1, 1, -1, -1, 1, -1, -1]
        assert gene_info.intron_profiles.profiles["ENSMUST00000001714.7"] == [-2, -2, -2, -2, -2, -2, -1, -1, -2]
        assert gene_info.intron_profiles.profiles["ENSMUST00000001715.7"] == [-2, -2, -2, -2, -2, -2, -1, -1, 1]

    def test_exon_profiles(self):
        gene_info = GeneInfo([self.gene_db], self.gffutils_db)
        assert gene_info.exon_profiles.features[0] == (1000, 1100)
        assert gene_info.exon_profiles.features[-1] == (9500, 10000)
        assert len(gene_info.exon_profiles.features) == 11
        assert gene_info.exon_profiles.profiles["ENSMUST00000001712.7"] == [1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1]
        assert gene_info.exon_profiles.profiles["ENSMUST00000001714.7"] == [-2, -2, -2, -2, -2, -2, -2, 1, -2, -2, -2]

        assert gene_info.split_exon_profiles.features[2] == (2101, 2200)
        assert gene_info.split_exon_profiles.features[-1] == (9500, 10000)
        assert len(gene_info.split_exon_profiles.features) == 11
        assert gene_info.split_exon_profiles.profiles["ENSMUST00000001712.7"] == [1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1]
        assert gene_info.split_exon_profiles.profiles["ENSMUST00000001714.7"] == [-2, -2, -2, -2, -2, -2, -2, 1, -2, -2,
                                                                                  -2]

    def test_from_region(self):
        chr_id = "chr10"
        start, end = 1000, 2000
        gene_info = GeneInfo.from_region(chr_id, start, end)

        assert gene_info.chr_id == chr_id
        assert gene_info.start == start
        assert gene_info.end == end
        assert gene_info.db is None
        assert gene_info.gene_db_list == []
        assert gene_info.empty()

    def test_from_models(self):
        from src.gene_info import TranscriptModel

        # Create a simple transcript model
        transcript_model = TranscriptModel("chr1", "+", "test_transcript", "test_gene", [(100, 200), (300, 400)], TranscriptModelType.known)

        gene_info = GeneInfo.from_models([transcript_model])

        assert gene_info.chr_id == "chr1"
        assert gene_info.start == 100
        assert gene_info.end == 400
        assert gene_info.gene_id_map["test_transcript"] == "test_gene"
        assert gene_info.isoform_strands["test_transcript"] == "+"
        assert gene_info.all_isoforms_exons["test_transcript"] == [(100, 200), (300, 400)]

    def test_serialization(self):
        import io
        from src.serialization import write_int, write_string

        # Create test gene info
        gene_info = GeneInfo([self.gene_db], self.gffutils_db)

        # Serialize
        outfile = io.BytesIO()
        gene_info.serialize(outfile)
        outfile.seek(0)

        # Deserialize
        new_gene_info = GeneInfo.deserialize(outfile, self.gffutils_db)

        # Check key attributes
        assert new_gene_info.chr_id == gene_info.chr_id
        assert new_gene_info.start == gene_info.start
        assert new_gene_info.end == gene_info.end
        assert new_gene_info.delta == gene_info.delta

    def test_get_gene_regions(self):
        gene_info = GeneInfo([self.gene_db], self.gffutils_db)
        regions = gene_info.get_gene_regions()

        assert isinstance(regions, dict)
        assert self.gene_db.id in regions
        assert regions[self.gene_db.id] == (self.gene_db.start, self.gene_db.end)

    def test_empty(self):
        gene_info = GeneInfo([self.gene_db], self.gffutils_db)
        assert not gene_info.empty()

        # Test empty gene info
        empty_gene_info = GeneInfo.from_region("chr1", 1, 100)
        assert empty_gene_info.empty()

    def test_feature_properties(self):
        gene_info = GeneInfo([self.gene_db], self.gffutils_db)

        # Test exon properties
        exon_props = gene_info.exon_property_map
        assert len(exon_props) == len(gene_info.exon_profiles.features)
        for prop in exon_props:
            assert isinstance(prop, FeatureInfo)
            assert prop.chr_id == gene_info.chr_id
            assert len(prop.gene_ids) > 0

        # Test intron properties
        intron_props = gene_info.intron_property_map
        assert len(intron_props) == len(gene_info.intron_profiles.features)
        for prop in intron_props:
            assert isinstance(prop, FeatureInfo)
            assert prop.chr_id == gene_info.chr_id
            assert len(prop.gene_ids) > 0
