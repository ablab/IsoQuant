############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from types import SimpleNamespace

import pytest

from isoquant_lib.alignment.read_strand import ReadStrandPredictor, StrandSource
from isoquant_lib.common import junctions_from_blocks
from isoquant_lib.terminal_prediction.polya_finder import PolyAInfo

PLUS_EXONS = [(1000, 1100), (1201, 1300), (1401, 1600)]
MINUS_EXONS = [(5000, 5100), (5201, 5300)]
LONG_MINUS_EXONS = [(10000, 10100), (29901, 30000)]


def _gene_info(transcripts):
    # transcripts: id -> (strand, exons)
    return SimpleNamespace(
        all_isoforms_exons={t: exons for t, (_, exons) in transcripts.items()},
        all_isoforms_introns={t: junctions_from_blocks(exons) for t, (_, exons) in transcripts.items()},
        isoform_strands={t: strand for t, (strand, _) in transcripts.items()})


def _reference(motifs, length=40000):
    # motifs: (intron, donor, acceptor) on the forward strand; everything else is N
    seq = ["N"] * length
    for (start, end), donor, acceptor in motifs:
        seq[start - 1:start + 1] = donor
        seq[end - 2:end] = acceptor
    return "".join(seq)


TRANSCRIPTS = {"P": ("+", PLUS_EXONS), "M": ("-", MINUS_EXONS), "L": ("-", LONG_MINUS_EXONS)}
NOVEL_FWD = (2001, 2100)
NOVEL_REV = (2301, 2400)
NOVEL_NONCANONICAL = (2601, 2700)
REFERENCE = _reference([(NOVEL_FWD, "GT", "AG"), (NOVEL_REV, "CT", "AC")])


def _predict(read_exons, polya=PolyAInfo(-1, -1, -1, -1), mapped_strand="+", stranded=None,
             transcripts=None, reference=REFERENCE):
    predictor = ReadStrandPredictor(_gene_info(TRANSCRIPTS if transcripts is None else transcripts),
                                    reference, stranded)
    return predictor.predict(read_exons, polya, mapped_strand)


def _read_with_introns(*introns):
    exons = []
    start = introns[0][0] - 50
    for intron in introns:
        exons.append((start, intron[0] - 1))
        start = intron[1] + 1
    exons.append((start, start + 50))
    return exons


class TestTail:
    @pytest.mark.parametrize("polya, strand", [
        (PolyAInfo(1700, -1, -1, -1), "+"),
        (PolyAInfo(-1, -1, 1650, -1), "+"),
        (PolyAInfo(-1, 900, -1, -1), "-"),
        (PolyAInfo(-1, -1, -1, 950), "-"),
    ])
    def test_tail_at_one_end(self, polya, strand):
        assert _predict([(1000, 1600)], polya) == (strand, StrandSource.polya_tail)

    def test_tail_wins_over_splice_sites(self):
        assert _predict(PLUS_EXONS, PolyAInfo(-1, 900, -1, -1)) == ("-", StrandSource.polya_tail)

    def test_tails_at_both_ends_say_nothing(self):
        assert _predict(PLUS_EXONS, PolyAInfo(1700, 900, -1, -1)) == ("+", StrandSource.splice_sites)


class TestSpliceSites:
    def test_annotated_introns(self):
        assert _predict(PLUS_EXONS) == ("+", StrandSource.splice_sites)
        assert _predict(MINUS_EXONS) == ("-", StrandSource.splice_sites)

    @pytest.mark.parametrize("intron, strand", [(NOVEL_FWD, "+"), (NOVEL_REV, "-")])
    def test_canonical_sites_of_novel_introns(self, intron, strand):
        assert _predict(_read_with_introns(intron)) == (strand, StrandSource.splice_sites)

    def test_annotated_and_canonical_agree(self):
        read_exons = [(1000, 1100), (1201, 2000), (2101, 2200)]
        assert junctions_from_blocks(read_exons)[1] == NOVEL_FWD
        assert _predict(read_exons) == ("+", StrandSource.splice_sites)

    def test_contradicting_sites(self):
        assert _predict(_read_with_introns(NOVEL_FWD, NOVEL_REV)) == (None, StrandSource.undecided)
        read_exons = [(1000, 1100), (1201, 2300), (2401, 2500)]
        assert junctions_from_blocks(read_exons)[1] == NOVEL_REV
        assert _predict(read_exons) == (None, StrandSource.undecided)

    def test_shifted_noncanonical_intron_says_nothing(self):
        # exact match only: an annotated intron shifted by 2 bp counts by its (here absent) splice sites
        read_exons = [(1000, 1102), (1203, 1300)]
        assert _predict(read_exons) == (None, StrandSource.undecided)
        assert _predict(_read_with_introns(NOVEL_NONCANONICAL)) == (None, StrandSource.undecided)

    def test_no_reference(self):
        assert _predict(_read_with_introns(NOVEL_FWD), reference=None) == (None, StrandSource.undecided)
        assert _predict(PLUS_EXONS, reference=None) == ("+", StrandSource.splice_sites)


class TestTranscriptOverlap:
    @pytest.mark.parametrize("read_exons, strand", [
        ([(1050, 1090)], "+"),
        ([(1050, 1250)], "+"),
        ([(5050, 5250)], "-"),
        # inside a long transcript, far from its start
        ([(20000, 20100)], "-"),
    ])
    def test_one_strand(self, read_exons, strand):
        assert _predict(read_exons) == (strand, StrandSource.transcript_overlap)

    def test_both_strands(self):
        transcripts = dict(TRANSCRIPTS, A=("-", [(1500, 1700)]))
        assert _predict([(1550, 1650)], transcripts=transcripts) == (None, StrandSource.undecided)

    def test_no_transcript(self):
        assert _predict([(7000, 7100)]) == (None, StrandSource.undecided)


class TestStrandedLibrary:
    @pytest.mark.parametrize("stranded, mapped_strand, strand", [
        ("forward", "+", "+"), ("forward", "-", "-"), ("reverse", "+", "-"), ("reverse", "-", "+")])
    def test_last_resort(self, stranded, mapped_strand, strand):
        assert _predict([(7000, 7100)], mapped_strand=mapped_strand, stranded=stranded) == \
            (strand, StrandSource.stranded_library)

    def test_lower_priority_than_the_read(self):
        assert _predict(PLUS_EXONS, mapped_strand="+", stranded="reverse") == ("+", StrandSource.splice_sites)
        assert _predict([(5050, 5250)], mapped_strand="+", stranded="forward") == \
            ("-", StrandSource.transcript_overlap)

    def test_resolves_contradictions(self):
        assert _predict(_read_with_introns(NOVEL_FWD, NOVEL_REV), mapped_strand="-", stranded="forward") == \
            ("-", StrandSource.stranded_library)

    def test_none(self):
        assert _predict([(7000, 7100)], stranded="none") == (None, StrandSource.undecided)
