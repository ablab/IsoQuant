############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import gzip
from types import SimpleNamespace

import pytest

import isoquant
from isoquant_lib.alignment.alignment_info import AlignmentInfo
from isoquant_lib.alignment.alignment_processor import AlignmentCollector, PolyATrimmed
from isoquant_lib.assignment.isoform_assignment import MatchClassification, ReadAssignmentType
from isoquant_lib.assignment.long_read_assigner import AmbiguityResolvingMethod, LongReadAssigner
from isoquant_lib.assignment.long_read_profiles import CombinedProfileConstructor
from isoquant_lib.gene_info import GeneInfo, TranscriptModel, TranscriptModelType
from isoquant_lib.terminal_prediction.external_polya import (
    NO_TAIL_FOUND,
    PRESENT_FORWARD,
    PRESENT_NO_HINT,
    PRESENT_REVERSE,
    UNKNOWN,
    ExternalPolyAStatus,
    external_polya_from_tag,
    has_any_polya,
    hint_from_strandedness,
    inject_external_polya,
    load_polya_read_dict,
    normalize_polya_reads,
    reconcile_polya,
    resolve_external_strand,
)
from isoquant_lib.terminal_prediction.polya_finder import PolyAFinder, PolyAInfo
from isoquant_lib.terminal_prediction.polya_verification import PolyAFixer
from isoquant_lib.utils.string_pools import StringPoolManager

from isoquant_tests.test_alignment_info import Alignment, PAIRS, SEQ1, TUPLES2

FLNC_HEADER = "id,strand,fivelen,threelen,polyAlen,insertlen,primer\n"


def _normalize(tmp_path, content: str, source: str, gz: bool = False):
    src = tmp_path / ("input.gz" if gz else "input")
    if gz:
        with gzip.open(src, "wt") as f:
            f.write(content)
    else:
        src.write_text(content)
    dst = tmp_path / "normalized.tsv"
    count = normalize_polya_reads(str(src), source, str(dst))
    return count, dst.read_text(), load_polya_read_dict(str(dst))


class TestParsePolyATrimmed:
    @pytest.mark.parametrize("value, expected", [
        ("none", PolyATrimmed.none),
        ("all", PolyATrimmed.all),
        ("stranded", PolyATrimmed.stranded),
    ])
    def test_plain_modes(self, value, expected):
        assert isoquant.parse_polya_trimmed(value) == (expected, None, None)

    def test_tag(self):
        assert isoquant.parse_polya_trimmed("tag:pt") == (PolyATrimmed.tag, "pt", None)

    def test_list_is_made_absolute(self, tmp_path, monkeypatch):
        (tmp_path / "ids.txt").write_text("r1\n")
        monkeypatch.chdir(tmp_path)
        mode, tag, file_name = isoquant.parse_polya_trimmed("list:ids.txt")
        assert (mode, tag, file_name) == (PolyATrimmed.list, None, str(tmp_path / "ids.txt"))

    def test_flnc(self, tmp_path):
        report = tmp_path / "flnc.report.csv"
        report.write_text(FLNC_HEADER)
        assert isoquant.parse_polya_trimmed("flnc:" + str(report)) == (PolyATrimmed.flnc, None, str(report))

    @pytest.mark.parametrize("value", [
        "unknown", "tag", "tag:", "list:", "flnc:", "none:x", "stranded:x", "list:/nonexistent/ids.txt",
    ])
    def test_invalid(self, value):
        with pytest.raises(ValueError):
            isoquant.parse_polya_trimmed(value)

    def test_flnc_with_wrong_header(self, tmp_path):
        report = tmp_path / "ids.txt"
        report.write_text("r1\nr2\n")
        with pytest.raises(ValueError):
            isoquant.parse_polya_trimmed("flnc:" + str(report))


class TestNormalization:
    def test_list(self, tmp_path):
        content = "# comment\n\nr1\nr2\t+\nr3 -\nr4\tx\nr5\t0\n"
        count, text, polya_dict = _normalize(tmp_path, content, "list")
        assert count == 5
        assert text == "r1\t.\nr2\t+\nr3\t-\nr4\t.\nr5\t.\n"
        assert polya_dict == {"r1": PRESENT_NO_HINT, "r2": PRESENT_FORWARD, "r3": PRESENT_REVERSE,
                              "r4": PRESENT_NO_HINT, "r5": PRESENT_NO_HINT}

    def test_list_never_reports_no_tail(self, tmp_path):
        _, _, polya_dict = _normalize(tmp_path, "r1\t0\nr2\t-1\n", "list")
        assert all(v.status == ExternalPolyAStatus.present for v in polya_dict.values())

    def test_list_gzipped(self, tmp_path):
        count, text, _ = _normalize(tmp_path, "r1\nr2\t-\n", "list", gz=True)
        assert count == 2
        assert text == "r1\t.\nr2\t-\n"

    @pytest.mark.parametrize("delim", [",", "\t"])
    def test_flnc(self, tmp_path, delim):
        rows = [FLNC_HEADER.strip(), "m1/1/ccs,+,30,40,46,1000,bc1", "m1/2/ccs,+,30,40,0,1000,bc1",
                "m1/3/ccs,-,30,40,46,1000,bc1", "m1/4/ccs,?,30,40,46,1000,bc1", "m1/5/ccs,-,30,40,0,1000,bc1"]
        content = "\n".join(r.replace(",", delim) for r in rows) + "\n"
        count, text, polya_dict = _normalize(tmp_path, content, "flnc")
        assert count == 5
        assert text == "m1/1/ccs\t+\nm1/2/ccs\t0\nm1/3/ccs\t-\nm1/4/ccs\t.\nm1/5/ccs\t0\n"
        assert polya_dict == {"m1/1/ccs": PRESENT_FORWARD, "m1/2/ccs": NO_TAIL_FOUND, "m1/3/ccs": PRESENT_REVERSE,
                              "m1/4/ccs": PRESENT_NO_HINT, "m1/5/ccs": NO_TAIL_FOUND}

    def test_flnc_without_strand_column(self, tmp_path):
        content = "id,polyAlen\nm1/1/ccs,30\nm1/2/ccs,0\n"
        _, text, _ = _normalize(tmp_path, content, "flnc")
        assert text == "m1/1/ccs\t+\nm1/2/ccs\t0\n"

    def test_shared_instances(self, tmp_path):
        _, _, polya_dict = _normalize(tmp_path, "r1\nr2\n", "list")
        assert polya_dict["r1"] is polya_dict["r2"]

    def test_empty_input(self, tmp_path):
        count, text, polya_dict = _normalize(tmp_path, "# nothing\n", "list")
        assert count == 0 and text == "" and polya_dict == {}


class MockTaggedAlignment:
    def __init__(self, tags: dict):
        self.tags = tags

    def has_tag(self, tag):
        return tag in self.tags

    def get_tag(self, tag):
        return self.tags[tag]


class TestTag:
    @pytest.mark.parametrize("tags, expected", [
        ({}, UNKNOWN),
        ({"pt": -1}, UNKNOWN),
        ({"pt": -1, "TS": "+"}, UNKNOWN),
        ({"pt": "46"}, UNKNOWN),
        ({"pt": 0}, PRESENT_NO_HINT),
        ({"pt": 46}, PRESENT_NO_HINT),
        ({"pt": 46, "TS": "+"}, PRESENT_FORWARD),
        ({"pt": 46, "TS": "-"}, PRESENT_REVERSE),
        ({"pt": 46, "TS": "."}, PRESENT_NO_HINT),
    ])
    def test_tag(self, tags, expected):
        assert external_polya_from_tag(MockTaggedAlignment(tags), "pt") == expected

    @pytest.mark.parametrize("stranded, expected", [("forward", "+"), ("reverse", "-"), ("none", None), (None, None)])
    def test_hint_from_strandedness(self, stranded, expected):
        assert hint_from_strandedness(stranded) == expected


class TestStrandResolution:
    @pytest.mark.parametrize("external, mapped, expected", [
        (PRESENT_FORWARD, "+", "+"),
        (PRESENT_FORWARD, "-", "-"),
        (PRESENT_REVERSE, "+", "-"),
        (PRESENT_REVERSE, "-", "+"),
        (PRESENT_NO_HINT, "+", None),
        (NO_TAIL_FOUND, "+", None),
        (UNKNOWN, "-", None),
    ])
    def test_resolve(self, external, mapped, expected):
        assert resolve_external_strand(external, mapped) == expected


def _polya(ext_a=-1, ext_t=-1, int_a=-1, int_t=-1) -> tuple:
    # a tuple, so that parametrized cases never share a mutable PolyAInfo
    return ext_a, ext_t, int_a, int_t


def _positions(info: PolyAInfo):
    return info.external_polya_pos, info.external_polyt_pos, info.internal_polya_pos, info.internal_polyt_pos


class TestReconcile:
    @pytest.mark.parametrize("before, strand, no_tail, after", [
        # agreeing tail keeps its exact position
        (_polya(ext_a=500), "+", False, (500, -1, -1, -1)),
        (_polya(int_t=100), "-", False, (-1, -1, -1, 100)),
        # opposite side only is cleared
        (_polya(ext_t=100, int_t=110), "+", False, (-1, -1, -1, -1)),
        (_polya(ext_a=500, int_a=490), "-", False, (-1, -1, -1, -1)),
        # both sides: only the opposite one is cleared
        (_polya(ext_a=500, ext_t=100), "+", False, (500, -1, -1, -1)),
        (_polya(ext_a=500, ext_t=100), "-", False, (-1, 100, -1, -1)),
        # no tail clears everything
        ((500, 100, 490, 110), None, True, (-1, -1, -1, -1)),
        # nothing to reconcile with
        ((500, 100, 490, 110), None, False, (500, 100, 490, 110)),
    ])
    def test_reconcile(self, before, strand, no_tail, after):
        info = PolyAInfo(*before)
        reconcile_polya(info, strand, no_tail)
        assert _positions(info) == after

    @pytest.mark.parametrize("before, strand, after", [
        (_polya(), "+", (401, -1, -1, -1)),
        (_polya(), "-", (-1, 99, -1, -1)),
        (_polya(int_a=390), "+", (-1, -1, 390, -1)),
        (_polya(ext_t=90), "-", (-1, 90, -1, -1)),
        (_polya(), None, (-1, -1, -1, -1)),
        # a tail at the other read end does not count
        (_polya(ext_a=410), "-", (410, 99, -1, -1)),
        (_polya(int_t=110), "+", (401, -1, -1, 110)),
    ])
    def test_inject(self, before, strand, after):
        info = PolyAInfo(*before)
        inject_external_polya(info, [(100, 200), (300, 400)], strand)
        assert _positions(info) == after


class TestAddPolyAInfo:
    """The alignment from test_alignment_info has an internal polyA that trims the last exons."""

    @staticmethod
    def _alignment_info():
        return AlignmentInfo(Alignment('query', TUPLES2, SEQ1, 51054, 35000, 'reference', 789, PAIRS))

    def test_default_unchanged(self):
        alignment_info = self._alignment_info()
        alignment_info.add_polya_info(PolyAFinder(), PolyAFixer(1))
        assert alignment_info.exons_changed is True
        assert alignment_info.polya_info.internal_polya_pos == 51221

    def test_agreeing_strand_keeps_sequence_tail(self):
        alignment_info = self._alignment_info()
        alignment_info.add_polya_info(PolyAFinder(), PolyAFixer(1), external_strand="+")
        assert alignment_info.exons_changed is True
        assert _positions(alignment_info.polya_info) == (-1, -1, 51221, -1)

    def test_contradicted_tail_does_not_trim(self):
        alignment_info = self._alignment_info()
        untrimmed_exons = list(alignment_info.read_exons)
        alignment_info.add_polya_info(PolyAFinder(), PolyAFixer(1), external_strand="-")
        assert alignment_info.exons_changed is False
        assert alignment_info.read_exons == untrimmed_exons
        assert _positions(alignment_info.polya_info) == (-1, untrimmed_exons[0][0] - 1, -1, -1)

    def test_no_tail_does_not_trim(self):
        alignment_info = self._alignment_info()
        untrimmed_exons = list(alignment_info.read_exons)
        alignment_info.add_polya_info(PolyAFinder(), PolyAFixer(1), external_no_tail=True)
        assert alignment_info.exons_changed is False
        assert alignment_info.read_exons == untrimmed_exons
        assert not has_any_polya(alignment_info.polya_info)


def _collector(mode: PolyATrimmed, fl_data: bool = False) -> AlignmentCollector:
    collector = AlignmentCollector.__new__(AlignmentCollector)
    collector.params = SimpleNamespace(polya_trimmed=mode, fl_data=fl_data)
    return collector


def _read_assignment(strand: str, mapped_strand: str, polya: tuple):
    return SimpleNamespace(strand=strand, mapped_strand=mapped_strand, polya_info=PolyAInfo(*polya),
                           corrected_exons=[(100, 200), (300, 400)])


class TestArtificialPolyA:
    @pytest.mark.parametrize("external, strand, before, after", [
        # side unknown: placed by the assigned strand, only when the sequence found nothing
        (PRESENT_NO_HINT, "+", _polya(), (401, -1, -1, -1)),
        (PRESENT_NO_HINT, "-", _polya(), (-1, 99, -1, -1)),
        (PRESENT_NO_HINT, ".", _polya(), (-1, -1, -1, -1)),
        (PRESENT_NO_HINT, "-", _polya(ext_a=410), (410, -1, -1, -1)),
        # side known or no information: phase 2 does nothing
        (PRESENT_FORWARD, "+", _polya(), (-1, -1, -1, -1)),
        (NO_TAIL_FOUND, "+", _polya(), (-1, -1, -1, -1)),
        (UNKNOWN, "+", _polya(), (-1, -1, -1, -1)),
    ])
    @pytest.mark.parametrize("mode", [PolyATrimmed.tag, PolyATrimmed.list, PolyATrimmed.flnc])
    def test_external_modes(self, mode, external, strand, before, after):
        read_assignment = _read_assignment(strand, "+", before)
        _collector(mode).add_artificial_polya(read_assignment, external)
        assert _positions(read_assignment.polya_info) == after

    @pytest.mark.parametrize("mode, strand, mapped_strand, after", [
        (PolyATrimmed.stranded, "+", "-", (401, -1, -1, -1)),
        (PolyATrimmed.stranded, ".", "-", (-1, -1, -1, -1)),
        (PolyATrimmed.all, ".", "-", (-1, 99, -1, -1)),
        (PolyATrimmed.none, "+", "+", (-1, -1, -1, -1)),
    ])
    def test_existing_modes(self, mode, strand, mapped_strand, after):
        read_assignment = _read_assignment(strand, mapped_strand, _polya())
        _collector(mode).add_artificial_polya(read_assignment)
        assert _positions(read_assignment.polya_info) == after

    @pytest.mark.parametrize("mode", [PolyATrimmed.stranded, PolyATrimmed.all])
    def test_existing_modes_keep_a_detected_tail(self, mode):
        read_assignment = _read_assignment("+", "+", _polya(ext_a=410))
        _collector(mode).add_artificial_polya(read_assignment)
        assert _positions(read_assignment.polya_info) == (410, -1, -1, -1)

    @pytest.mark.parametrize("strand, before, after", [
        # full-length reads: the 3' end is the polyA site when no tail was found
        ("+", _polya(), (401, -1, -1, -1)),
        ("-", _polya(), (-1, 99, -1, -1)),
        (".", _polya(), (-1, -1, -1, -1)),
        # a tail detected at the 3' end is kept
        ("+", _polya(int_a=380), (-1, -1, 380, -1)),
        ("-", _polya(int_t=110), (-1, -1, -1, 110)),
        # a genomic A-run at the 5' end does not replace the 3' site
        ("-", _polya(int_a=380), (-1, 99, 380, -1)),
        ("+", _polya(int_t=110), (401, -1, -1, 110)),
        # a soft-clipped tail at the 5' end: the read comes from the other strand
        ("-", _polya(ext_a=410), (410, -1, -1, -1)),
        ("+", _polya(ext_t=90), (-1, 90, -1, -1)),
    ])
    def test_fl_data(self, strand, before, after):
        read_assignment = _read_assignment(strand, "+", before)
        _collector(PolyATrimmed.none, fl_data=True).add_artificial_polya(read_assignment)
        assert _positions(read_assignment.polya_info) == after

    @pytest.mark.parametrize("mode, external, after", [
        # an explicit --polya_trimmed source takes precedence over --fl_data
        (PolyATrimmed.flnc, NO_TAIL_FOUND, (-1, -1, -1, -1)),
        (PolyATrimmed.tag, UNKNOWN, (-1, -1, -1, -1)),
        (PolyATrimmed.list, PRESENT_NO_HINT, (-1, 99, -1, -1)),
        (PolyATrimmed.stranded, UNKNOWN, (-1, 99, -1, -1)),
    ])
    def test_fl_data_with_polya_trimmed(self, mode, external, after):
        read_assignment = _read_assignment("-", "+", _polya())
        _collector(mode, fl_data=True).add_artificial_polya(read_assignment, external)
        assert _positions(read_assignment.polya_info) == after


ASSIGNER_PARAMS = SimpleNamespace(
    delta=4, minor_exon_extension=50, major_exon_extension=300, min_abs_exon_overlap=10, min_rel_exon_overlap=0.2,
    max_suspicious_intron_abs_len=0, max_suspicious_intron_rel_len=0.0, max_fake_terminal_exon_len=20,
    micro_intron_length=50, max_intron_abs_diff=30, max_intron_rel_diff=0.2, apa_delta=50, minimal_exon_overlap=5,
    minimal_intron_absence_overlap=20, max_intron_shift=30, max_missed_exon_len=50,
    resolve_ambiguous=AmbiguityResolvingMethod.monoexon_and_fsm, count_exons=False)

EXONS_A = [(1000, 1100), (1201, 1300), (1401, 1600)]
EXONS_B = [(1000, 1100), (1201, 1300), (1401, 1720)]


def _model(transcript_id: str, strand: str, exons) -> TranscriptModel:
    return TranscriptModel("chr1", strand, transcript_id, "G", exons, TranscriptModelType.known)


def _summary(read_assignment) -> tuple:
    return (read_assignment.assignment_type, read_assignment.gene_assignment_type,
            [(m.assigned_transcript, m.match_classification,
              sorted(e.event_type.name for e in m.match_subclassifications)) for m in read_assignment.isoform_matches])


def _assigner(models):
    gene_info = GeneInfo.from_models(models, ASSIGNER_PARAMS.delta)
    return (LongReadAssigner(gene_info, ASSIGNER_PARAMS, StringPoolManager()),
            CombinedProfileConstructor(gene_info, ASSIGNER_PARAMS))


def _detected(models, read_exons, strand: str) -> tuple:
    # the assignment of the same read when its tail is found in the sequence at the site the options add
    assigner, profile_constructor = _assigner(models)
    site = _polya(ext_a=read_exons[-1][1] + 1) if strand == "+" else _polya(ext_t=read_exons[0][0] - 1)
    profile = profile_constructor.construct_profiles(read_exons, PolyAInfo(*site), [])
    return _summary(assigner.assign_to_isoform("read", profile))


class _Read:
    """The part of AlignmentInfo the collector uses around the assignment."""

    def __init__(self, read_exons, polya: tuple):
        self.read_exons = read_exons
        self.polya_info = PolyAInfo(*polya)
        self.combined_profile = None

    def construct_profiles(self, profile_constructor):
        self.combined_profile = profile_constructor.construct_profiles(self.read_exons, self.polya_info, [])


def _predictor(strand):
    return SimpleNamespace(predict=lambda read_exons, polya_info, mapped_strand: (strand, None))


class TestMarkBeforeAssignment:
    """The site goes in before the profiles are built, on the predicted read strand."""

    @pytest.mark.parametrize("predicted, before, marked, after", [
        ("+", _polya(), "+", (401, -1, -1, -1)),
        ("-", _polya(), "-", (-1, 99, -1, -1)),
        (None, _polya(), None, (-1, -1, -1, -1)),
        # a tail found at the 3' end needs no site
        ("+", _polya(int_a=380), None, (-1, -1, 380, -1)),
        # a soft-clipped tail at the 5' end: the read comes from the other strand
        ("+", _polya(ext_t=90), None, (-1, 90, -1, -1)),
    ])
    def test_fl_data(self, predicted, before, marked, after):
        read = _Read([(100, 200), (300, 400)], before)
        collector = _collector(PolyATrimmed.none, fl_data=True)
        assert collector.set_polya_before_assignment(_predictor(predicted), read, UNKNOWN, "+") == marked
        assert _positions(read.polya_info) == after

    def test_all_uses_the_mapped_strand(self):
        read = _Read([(100, 200), (300, 400)], _polya())
        collector = _collector(PolyATrimmed.all)
        assert collector.set_polya_before_assignment(_predictor("+"), read, UNKNOWN, "-") == "-"
        assert _positions(read.polya_info) == (-1, 99, -1, -1)

    @pytest.mark.parametrize("external, marked", [(PRESENT_NO_HINT, "-"), (UNKNOWN, None), (NO_TAIL_FOUND, None)])
    def test_external_tail_without_a_side(self, external, marked):
        read = _Read([(100, 200), (300, 400)], _polya())
        assert _collector(PolyATrimmed.tag).set_polya_before_assignment(_predictor("-"), read, external, "+") == marked

    def test_nothing_without_options(self):
        read = _Read([(100, 200), (300, 400)], _polya())
        assert _collector(PolyATrimmed.none).set_polya_before_assignment(None, read, UNKNOWN, "+") is None
        assert _positions(read.polya_info) == (-1, -1, -1, -1)


class TestMovePolyASite:
    """After the assignment the site follows the read strand; a moved site means the read is assigned again."""

    @staticmethod
    def _move(mode, marked, strand, fl_data=True, polya=_polya()):
        read = _Read([(100, 200), (300, 400)], polya)
        inject_external_polya(read.polya_info, read.read_exons, marked)
        read_assignment = SimpleNamespace(strand=strand, mapped_strand="+")
        moved = _collector(mode, fl_data).move_polya_site(read_assignment, read, UNKNOWN, marked)
        return moved, _positions(read.polya_info)

    def test_predicted_strand_confirmed(self):
        assert self._move(PolyATrimmed.none, "+", "+") == (False, (401, -1, -1, -1))

    def test_read_assigned_to_the_other_strand(self):
        assert self._move(PolyATrimmed.none, "+", "-") == (True, (-1, 99, -1, -1))

    def test_strand_not_predicted(self):
        assert self._move(PolyATrimmed.none, None, "-") == (True, (-1, 99, -1, -1))
        assert self._move(PolyATrimmed.none, None, ".") == (False, (-1, -1, -1, -1))

    def test_predicted_site_dropped_for_unknown_strand(self):
        assert self._move(PolyATrimmed.stranded, "+", ".", fl_data=False) == (True, (-1, -1, -1, -1))

    def test_all_never_moves(self):
        assert self._move(PolyATrimmed.all, "+", "-", fl_data=False) == (False, (401, -1, -1, -1))

    def test_detected_tail_stays(self):
        assert self._move(PolyATrimmed.none, None, "+", polya=_polya(int_a=380)) == (False, (-1, -1, 380, -1))


class TestAssignAgain:
    """A read whose site moved gets the assignment a tail found in its sequence at the new site gives."""

    @staticmethod
    def _assign(models, read_exons, predicted: str, strand: str):
        assigner, profile_constructor = _assigner(models)
        read = _Read(read_exons, _polya())
        collector = _collector(PolyATrimmed.none, fl_data=True)
        marked = collector.set_polya_before_assignment(_predictor(predicted), read, UNKNOWN, "+")
        read.construct_profiles(profile_constructor)
        read_assignment = assigner.assign_to_isoform("read", read.combined_profile)
        before = _summary(read_assignment)
        read_assignment.strand, read_assignment.mapped_strand = strand, "+"
        assigned_again = collector.move_polya_site(read_assignment, read, UNKNOWN, marked) and \
            AlignmentCollector.assign_again(assigner, profile_constructor, read, read_assignment)
        return before, _summary(read_assignment), assigned_again

    def test_predicted_strand_is_the_detected_tail(self):
        # same intron chain, different 3' ends: the site resolves the ambiguity before the assignment
        models = [_model("A", "+", EXONS_A), _model("B", "+", EXONS_B)]
        before, after, assigned_again = self._assign(models, EXONS_A, "+", "+")
        assert before == after == _detected(models, EXONS_A, "+")
        assert after[0] == ReadAssignmentType.unique and "correct_polya_site_right" in after[2][0][2]
        assert not assigned_again

    def test_wrong_prediction_is_assigned_again(self):
        models = [_model("A", "+", EXONS_A), _model("B", "+", EXONS_B)]
        before, after, assigned_again = self._assign(models, EXONS_A, "-", "+")
        assert before[0] == ReadAssignmentType.ambiguous
        assert after == _detected(models, EXONS_A, "+")
        assert assigned_again

    def test_site_contradicting_every_isoform_makes_the_read_inconsistent(self):
        # the read ends 100 bp before the only isoform: no longer a full splice match
        models = [_model("A", "+", EXONS_A)]
        read_exons = [(1000, 1100), (1201, 1300), (1401, 1500)]
        _, after, _ = self._assign(models, read_exons, "+", "+")
        assert after == _detected(models, read_exons, "+")
        assert after[0] == ReadAssignmentType.inconsistent_non_intronic
        assert after[2][0][1] == MatchClassification.novel_in_catalog and "fsm" not in after[2][0][2]

    def test_minus_strand(self):
        exons_b = [(880, 1100), (1201, 1300), (1401, 1600)]
        models = [_model("A", "-", EXONS_A), _model("B", "-", exons_b)]
        _, after, assigned_again = self._assign(models, EXONS_A, None, "-")
        assert after == _detected(models, EXONS_A, "-")
        assert after[0] == ReadAssignmentType.unique and "correct_polya_site_left" in after[2][0][2]
        assert assigned_again

    def test_unassigned_read_is_kept(self):
        before, after, assigned_again = self._assign([_model("A", "+", EXONS_A)], [(5000, 5100), (5201, 5300)],
                                                     None, "+")
        assert before == after and before[0] == ReadAssignmentType.noninformative
        assert not assigned_again
