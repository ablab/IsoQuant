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
        assert gene_info.exon_profiles.profiles["ENSMUST00000001712.7"] == [1, 1, -1, 1, -1,  1, -1, -1, -1, -1, 1]
        assert gene_info.exon_profiles.profiles["ENSMUST00000001714.7"] == [-2, -2, -2, -2, -2, -2, -2, 1, -2, -2, -2]

        assert gene_info.split_exon_profiles.features[2] == (2101, 2200)
        assert gene_info.split_exon_profiles.features[-1] == (9500, 10000)
        assert len(gene_info.split_exon_profiles.features) == 11
        assert gene_info.split_exon_profiles.profiles["ENSMUST00000001712.7"] == [1, 1, -1, 1, -1,  1, -1, -1, -1, -1, 1]
        assert gene_info.split_exon_profiles.profiles["ENSMUST00000001714.7"] == [-2, -2, -2, -2, -2, -2, -2, 1, -2, -2, -2]