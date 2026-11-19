############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import gzip
import logging
from argparse import Namespace
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from isoquant_lib.alignment.alignment_processor import AlignmentCollector
from isoquant_lib.gene_info import TranscriptModel, TranscriptModelType
from isoquant_lib.terminal_prediction.cage_peaks import (
    DEFAULT_CAGE_SHIFT,
    CageFormatError,
    CagePeaks,
    load_cage_peaks,
)
from isoquant_lib.terminal_prediction.terminal_counter import TSSCounter


def write_bed(tmp_path, *lines, name="peaks.bed"):
    path = tmp_path / name
    path.write_text("".join(line + "\n" for line in lines))
    return str(path)


def bed6(chrom, start, end, strand="+"):
    return "\t".join([chrom, str(start), str(end), "peak", "0", strand])


def distances(peaks, positions, chrom="chr1", strand="+"):
    return peaks.distances(chrom, strand, np.array(positions)).tolist()


class TestDistances:
    def test_bed_coordinates_are_converted_to_one_based(self, tmp_path):
        # BED [100, 110) covers 1-based positions 101..110
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 100, 110)))
        assert distances(peaks, [99, 100, 101, 110, 111, 112]) == [2, 1, 0, 0, 1, 2]

    def test_closest_of_two_peaks(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 100, 110), bed6("chr1", 200, 210)))
        assert distances(peaks, [120, 190, 300]) == [10, 11, 90]

    def test_peak_contained_in_a_longer_one(self, tmp_path):
        # the last peak started is not the one reaching furthest
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 100, 500), bed6("chr1", 200, 210)))
        assert distances(peaks, [300, 505]) == [0, 5]

    def test_unsorted_file(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 200, 210), bed6("chr1", 100, 110)))
        assert distances(peaks, [105, 205]) == [0, 0]

    def test_strands_are_separate(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 100, 110, "+"), bed6("chr1", 500, 510, "-")))
        assert distances(peaks, [105], strand="+") == [0]
        assert distances(peaks, [105], strand="-") == [396]

    def test_no_peaks(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 100, 110, "+")))
        assert peaks.distances("chr1", "-", np.array([105])) is None
        assert peaks.distances("chr2", "+", np.array([105])) is None

    def test_is_confirmed(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 100, 110)))
        assert peaks.is_confirmed("chr1", "+", 105, 0)
        assert not peaks.is_confirmed("chr1", "+", 115, 0)
        assert peaks.is_confirmed("chr1", "+", 115, 5)
        assert not peaks.is_confirmed("chr2", "+", 105, 50)


class TestParsing:
    def test_numeric_chromosome_names_stay_strings(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, bed6("1", 100, 110)))
        assert peaks.has_chromosome("1")
        assert distances(peaks, [105], chrom="1") == [0]

    def test_headers_and_blank_lines_are_skipped(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, "track name=cage", "browser position chr1", "# comment", "",
                                    bed6("chr1", 100, 110)))
        assert distances(peaks, [105]) == [0]

    def test_gzipped_file(self, tmp_path):
        path = tmp_path / "peaks.bed.gz"
        with gzip.open(str(path), "wt") as bed:
            bed.write(bed6("chr1", 100, 110) + "\n")
        assert distances(CagePeaks(str(path)), [105]) == [0]

    @pytest.mark.parametrize("line", [bed6("chr1", 100, 110, "."), "chr1\t100\t110", "chr1\t100\t110\tpeak\t0"])
    def test_unstranded_peaks_match_both_strands(self, tmp_path, line):
        peaks = CagePeaks(write_bed(tmp_path, line))
        assert distances(peaks, [105], strand="+") == [0]
        assert distances(peaks, [105], strand="-") == [0]

    @pytest.mark.parametrize("line,message", [("chr1\t100", "at least 3 columns"),
                                              (bed6("chr1", 100, 110, "x"), "unexpected strand"),
                                              (bed6("chr1", "a", 110), "not integers")])
    def test_incorrect_lines(self, tmp_path, line, message):
        with pytest.raises(CageFormatError, match=message):
            CagePeaks(write_bed(tmp_path, line))

    def test_file_with_no_peaks(self, tmp_path):
        with pytest.raises(CageFormatError, match="no peaks"):
            CagePeaks(write_bed(tmp_path, "# comment"))

    def test_file_is_parsed_once(self, tmp_path):
        path = write_bed(tmp_path, bed6("chr1", 100, 110))
        assert load_cage_peaks(path) is load_cage_peaks(path)


def make_model(strand, exons, transcript_id="t1", chr_id="chr1"):
    return TranscriptModel(chr_id, strand, transcript_id, "g1", exons, TranscriptModelType.known)


class TestTranscriptModels:
    def test_tss_is_the_strand_specific_end(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 100, 110, "+"), bed6("chr1", 500, 510, "-")))
        forward = make_model("+", [(105, 200), (300, 400)])
        reverse = make_model("-", [(300, 400), (450, 505)])
        peaks.add_model_info([forward, reverse], 0)
        for model in (forward, reverse):
            assert model.additional_info["validated_by_cage"] == "True"
            assert model.additional_info["dist_to_closest_cage"] == "0"

    def test_model_outside_peak(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 100, 110)))
        model = make_model("+", [(113, 200)])
        peaks.add_model_info([model], 0)
        assert model.additional_info["validated_by_cage"] == "False"
        assert model.additional_info["dist_to_closest_cage"] == "3"

    def test_no_peaks_on_chromosome(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 100, 110)))
        model = make_model("+", [(105, 200)], chr_id="chr2")
        peaks.add_model_info([model], 0)
        assert model.additional_info["validated_by_cage"] == "False"
        assert "dist_to_closest_cage" not in model.additional_info

    def test_annotated_model_is_not_annotated_twice(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 100, 110)))
        model = make_model("+", [(105, 200)])
        peaks.add_model_info([model], 0)
        peaks.add_model_info([model], 0)
        assert model.additional_attributes_str().count("validated_by_cage") == 1


def make_read(strand, exons, chr_id="chr1"):
    return SimpleNamespace(strand=strand, corrected_exons=exons, chr_id=chr_id, cage_found=False)


class TestReadStatus:
    @pytest.fixture
    def collector(self, tmp_path):
        peaks = CagePeaks(write_bed(tmp_path, bed6("chr1", 100, 110, "+"), bed6("chr1", 500, 510, "-")))
        return SimpleNamespace(cage_peaks=peaks, cage_shift=DEFAULT_CAGE_SHIFT)

    @pytest.mark.parametrize("read,expected", [(make_read("+", [(105, 200)]), True),
                                               (make_read("+", [(50, 105)]), False),
                                               (make_read("-", [(300, 505)]), True),
                                               (make_read("-", [(505, 600)]), False),
                                               (make_read(".", [(105, 505)]), False),
                                               (make_read("+", []), False)])
    def test_read_start_follows_assigned_strand(self, collector, read, expected):
        AlignmentCollector.set_cage_status(collector, read)
        assert read.cage_found is expected

    def test_no_cage_file(self):
        read = make_read("+", [(105, 200)])
        AlignmentCollector.set_cage_status(SimpleNamespace(cage_peaks=None, cage_shift=0), read)
        assert read.cage_found is False


def make_counter(tmp_path, cage=None, **kwargs):
    args = Namespace(genedb=True, fl_data=True, read_group=None, cage=cage, **kwargs)
    return TSSCounter(args, str(tmp_path / "tss.tsv"))


def make_predictions(*rows):
    return pd.DataFrame(rows, columns=["chromosome", "strand", "prediction"])


class TestTSSValidation:
    def test_default_requires_prediction_inside_peak(self, tmp_path):
        counter = make_counter(tmp_path, write_bed(tmp_path, bed6("chr1", 100, 110)))
        result = counter._postprocess_predictions(make_predictions(("chr1", "+", 105), ("chr1", "+", 113)))
        assert result["validated_by_cage"].tolist() == [True, False]
        assert result["dist_to_closest_cage"].tolist() == [0, 3]

    def test_cage_shift_is_the_threshold(self, tmp_path):
        counter = make_counter(tmp_path, write_bed(tmp_path, bed6("chr1", 100, 110)), cage_shift=3)
        result = counter._postprocess_predictions(make_predictions(("chr1", "+", 113), ("chr1", "+", 114)))
        assert result["validated_by_cage"].tolist() == [True, False]

    def test_no_peaks_leave_distance_empty(self, tmp_path):
        counter = make_counter(tmp_path, write_bed(tmp_path, bed6("chr1", 100, 110, "+")))
        result = counter._postprocess_predictions(make_predictions(("chr1", "-", 105)))
        assert result["validated_by_cage"].tolist() == [False]
        assert result["dist_to_closest_cage"].isna().all()

    def test_chromosome_without_peaks_is_reported(self, tmp_path, caplog):
        counter = make_counter(tmp_path, write_bed(tmp_path, bed6("1", 100, 110)))
        with caplog.at_level(logging.WARNING, logger="IsoQuant"):
            counter._postprocess_predictions(make_predictions(("chr1", "+", 105)))
        assert "No CAGE peaks on chromosome chr1" in caplog.text

    def test_no_cage_file_adds_no_columns(self, tmp_path):
        counter = make_counter(tmp_path)
        predictions = make_predictions(("chr1", "+", 105))
        assert counter._extra_columns() == []
        assert list(counter._postprocess_predictions(predictions).columns) == list(predictions.columns)

    def test_output_columns_with_cage(self, tmp_path):
        counter = make_counter(tmp_path, write_bed(tmp_path, bed6("chr1", 100, 110)))
        assert counter._output_columns()[-2:] == ["validated_by_cage", "dist_to_closest_cage"]
