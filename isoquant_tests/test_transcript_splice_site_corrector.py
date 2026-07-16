############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Tests for the splice-site error corrector (ported from PR #108 and reworked to
# store only boundary-proximal deletions per read and count deleted bases in
# reference space, reusing common.get_intron_strand's canonical convention).
# Run with:
#   pytest isoquant_tests/test_transcript_splice_site_corrector.py -v

from types import SimpleNamespace
from unittest import TestCase

import pytest

from isoquant_lib.gene_info import TranscriptModelType
from isoquant_lib.model_construction.model_filter import ModelFilter
from isoquant_lib.model_construction.transcript_splice_site_corrector import (
    SpliceSiteCase,
    extract_relevant_deletions,
    extract_splice_site_locations_within_aligned_read,
    count_deletions_for_splice_site_locations,
    compute_most_common_case_of_deletions,
    correct_splice_site_errors,
    generate_updated_exon_list,
)


def _run_correction(exons, reads, chr_record, strand):
    """Standalone equivalent of the per-model correction flow in GBMC."""
    splice_site_cases = {}
    for read in reads:
        count_deletions_for_splice_site_locations(
            read.exons[0][0], read.exons[-1][1], read.boundary_deletions, exons, splice_site_cases)
    locations = correct_splice_site_errors(splice_site_cases, chr_record, strand)
    if not locations:
        return None
    return generate_updated_exon_list(splice_site_cases, locations, exons)


def _reads(boundary_deletions, ref_span, count=5):
    """Identical reads anchored at reference 0..ref_span with the given deletions."""
    return [SimpleNamespace(exons=[(0, ref_span)], boundary_deletions=list(boundary_deletions))
            for _ in range(count)]


class TestExtractRelevantDeletions(TestCase):

    def test_keeps_boundary_deletion_drops_interior(self):
        # single exon 100..200; D at 150 is interior, D at 192 is near boundary 200
        cigar = [(0, 50), (2, 2), (0, 40), (2, 4), (0, 6)]
        result = extract_relevant_deletions(cigar, 100, [(100, 200)])
        self.assertEqual(result, [(192, 4)])

    def test_no_deletions_returns_empty(self):
        self.assertEqual(extract_relevant_deletions([(4, 20), (0, 100)], 100, [(100, 200)]), [])

    def test_empty_cigar_or_exons(self):
        self.assertEqual(extract_relevant_deletions([], 100, [(100, 200)]), [])
        self.assertEqual(extract_relevant_deletions([(0, 10)], 100, []), [])

    def test_deletion_position_accounts_for_introns(self):
        # M10 (100-109), N90 intron (110-199), M10 (200-209), D4 at 210
        cigar = [(0, 10), (3, 90), (0, 10), (2, 4), (0, 6)]
        result = extract_relevant_deletions(cigar, 100, [(100, 109), (200, 214)])
        self.assertEqual(result, [(210, 4)])


class TestExtractSpliceSiteLocations(TestCase):

    def test_correct_splice_sites_are_extracted(self):
        result = extract_splice_site_locations_within_aligned_read(20, 40, [(1, 10), (20, 30), (40, 50)])
        self.assertEqual(result, [(20, False), (30, True), (40, False)])


class TestCountDeletionsForSpliceSiteLocations(TestCase):

    def test_reference_window_counts(self):
        # deletion [10,13] (len 4). exon start 10 -> window [10,17] = 4 deleted;
        # exon end 20 -> window [13,20] clips the last deleted base = 1; far end 40 -> 0.
        exons = [(1, 5), (10, 20), (40, 50)]
        splice_site_cases = {}
        count_deletions_for_splice_site_locations(1, 45, [(10, 4)], exons, splice_site_cases)
        self.assertEqual(splice_site_cases[10].deletions, {4: 1})
        self.assertEqual(splice_site_cases[10].location_is_end, False)
        self.assertEqual(splice_site_cases[20].deletions, {1: 1})
        self.assertEqual(splice_site_cases[40].deletions, {0: 1})

    def test_partial_overlap_counts_only_window_bases(self):
        # deletion [18,21]; exon end 20 -> window [13,20] overlaps [18,20] = 3
        splice_site_cases = {}
        count_deletions_for_splice_site_locations(1, 30, [(18, 4)], [(10, 20)], splice_site_cases)
        self.assertEqual(splice_site_cases[20].deletions, {3: 1})


class TestMostCommonCaseOfDeletions(TestCase):

    def test_distinct_most_common_case_for_location_start(self):
        self.assertEqual(compute_most_common_case_of_deletions({0: 10, 1: 2, 4: 20, 5: 1}, False), 4)

    def test_distinct_most_common_case_for_location_end(self):
        self.assertEqual(compute_most_common_case_of_deletions({0: 10, 1: 2, 4: 20, 5: 1}, True), -4)

    def test_no_distinct_most_common_del_returns_neg_one(self):
        self.assertEqual(compute_most_common_case_of_deletions({0: 10, 3: 20, 4: 20}, False), -1)


class TestExonListUpdater(TestCase):

    def test_error_at_location_start_is_corrected(self):
        result = generate_updated_exon_list({20: SpliceSiteCase(False, most_common_del=4)}, [20],
                                            [(1, 10), (20, 30), (40, 50)])
        self.assertEqual(result, [(1, 10), (24, 30), (40, 50)])

    def test_error_at_location_end_is_corrected(self):
        result = generate_updated_exon_list({30: SpliceSiteCase(True, most_common_del=-4)}, [30],
                                            [(1, 10), (20, 30), (40, 50)])
        self.assertEqual(result, [(1, 10), (20, 26), (40, 50)])


class TestCorrectSpliceSiteErrors(TestCase):

    def _cases(self, boundary, count, is_end, n_reads):
        return {boundary: SpliceSiteCase(is_end, deletions={count: n_reads})}

    def test_accepted_canonical_case_is_returned(self):
        chr_record = "C" * 31 + "AG" + "C" * 10   # acceptor for start 30+4 at ref[31:33]
        cases = self._cases(30, 4, False, 10)
        self.assertEqual(correct_splice_site_errors(cases, chr_record, "+"), [30])

    def test_below_min_reads_is_skipped(self):
        chr_record = "C" * 31 + "AG" + "C" * 10
        self.assertEqual(correct_splice_site_errors(self._cases(30, 4, False, 3), chr_record, "+"), [])

    def test_non_accepted_size_is_skipped(self):
        chr_record = "C" * 31 + "AG" + "C" * 10
        self.assertEqual(correct_splice_site_errors(self._cases(30, 2, False, 10), chr_record, "+"), [])


class TestCanonicalIndexing(TestCase):
    # exon end (donor)    -> chr_record[corrected : corrected + 2]
    # exon start(acceptor)-> chr_record[corrected - 3 : corrected - 1]

    def _canonical(self, chr_record, boundary, most_common_del, is_end, strand):
        count = abs(most_common_del)
        cases = {boundary: SpliceSiteCase(is_end, deletions={count: 10})}
        return correct_splice_site_errors(cases, chr_record, strand) == [boundary]

    def test_start_pos_strand(self):
        self.assertTrue(self._canonical("C" * 11 + "AG" + "C" * 13, 10, 4, False, "+"))

    def test_end_pos_strand(self):
        self.assertTrue(self._canonical("C" * 10 + "GT" + "C" * 14, 14, -4, True, "+"))

    def test_start_neg_strand(self):
        self.assertTrue(self._canonical("C" * 11 + "AC" + "C" * 13, 10, 4, False, "-"))

    def test_end_neg_strand(self):
        self.assertTrue(self._canonical("C" * 10 + "CT" + "C" * 14, 14, -4, True, "-"))

    def test_non_canonical_returns_false(self):
        self.assertFalse(self._canonical("N" * 26, 14, -4, True, "+"))


class TestSpliceSiteCorrectionFlow(TestCase):

    def test_error_in_start_on_pos_strand_is_corrected(self):
        reads = _reads([(10, 4)], 20)
        self.assertEqual(_run_correction([(0, 5), (10, 20)], reads, "C" * 11 + "AG" + "C" * 13, "+"),
                         [(0, 5), (14, 20)])

    def test_error_in_end_on_pos_strand_is_corrected(self):
        reads = _reads([(10, 4)], 30)
        self.assertEqual(_run_correction([(0, 14), (20, 30)], reads, "C" * 10 + "GT" + "C" * 20, "+"),
                         [(0, 10), (20, 30)])

    def test_error_in_start_on_neg_strand_is_corrected(self):
        reads = _reads([(10, 4)], 20)
        self.assertEqual(_run_correction([(0, 5), (10, 20)], reads, "C" * 11 + "AC" + "C" * 13, "-"),
                         [(0, 5), (14, 20)])

    def test_error_in_end_on_neg_strand_is_corrected(self):
        reads = _reads([(10, 4)], 30)
        self.assertEqual(_run_correction([(0, 14), (20, 30)], reads, "C" * 10 + "CT" + "C" * 20, "-"),
                         [(0, 10), (20, 30)])

    def test_min_accepted_del_case(self):
        reads = _reads([(10, 3)], 30)
        self.assertEqual(_run_correction([(0, 14), (20, 30)], reads, "C" * 11 + "CT" + "C" * 19, "-"),
                         [(0, 11), (20, 30)])

    def test_max_accepted_del_case(self):
        reads = _reads([(8, 6)], 30)
        self.assertEqual(_run_correction([(0, 14), (20, 30)], reads, "C" * 8 + "CT" + "C" * 22, "-"),
                         [(0, 8), (20, 30)])

    def test_dels_but_no_canonical_returns_none(self):
        self.assertIsNone(_run_correction([(0, 14), (20, 30)], _reads([(10, 4)], 30), "N" * 32, "-"))

    def test_not_enough_deletion_size_returns_none(self):
        reads = _reads([(10, 2)], 30)
        self.assertIsNone(_run_correction([(0, 14), (20, 30)], reads, "C" * 10 + "GT" + "C" * 20, "+"))

    def test_no_deletions_returns_none(self):
        self.assertIsNone(_run_correction([(0, 14), (20, 30)], _reads([], 30), "C" * 10 + "GT" + "C" * 20, "+"))


class TestGBMCIntegration(TestCase):

    def _make_constructor(self, chr_record, disabled=False):
        # correct_transcript_splice_sites now lives on ModelFilter (stage 4a)
        gbmc = ModelFilter.__new__(ModelFilter)
        gbmc.chr_record = chr_record
        gbmc.args = SimpleNamespace(no_splice_site_correction=disabled)
        # model list + read bindings now live on the shared ModelStore
        gbmc.store = SimpleNamespace()
        return gbmc

    def _model(self, model_type):
        return SimpleNamespace(transcript_id="t1", transcript_type=model_type,
                               strand="+", exon_blocks=[(0, 14), (20, 30)])

    def test_novel_model_end_is_corrected(self):
        gbmc = self._make_constructor("C" * 10 + "GT" + "C" * 20)
        model = self._model(TranscriptModelType.novel_not_in_catalog)
        gbmc.store.transcript_model_storage = [model]
        gbmc.store.transcript_read_ids = {"t1": _reads([(10, 4)], 30)}
        gbmc.correct_transcript_splice_sites()
        self.assertEqual(model.exon_blocks, [(0, 10), (20, 30)])

    def test_known_model_is_untouched(self):
        gbmc = self._make_constructor("C" * 10 + "GT" + "C" * 20)
        model = self._model(TranscriptModelType.known)
        gbmc.store.transcript_model_storage = [model]
        gbmc.store.transcript_read_ids = {"t1": _reads([(10, 4)], 30)}
        gbmc.correct_transcript_splice_sites()
        self.assertEqual(model.exon_blocks, [(0, 14), (20, 30)])

    def test_disabled_flag_skips_correction(self):
        gbmc = self._make_constructor("C" * 10 + "GT" + "C" * 20, disabled=True)
        model = self._model(TranscriptModelType.novel_not_in_catalog)
        gbmc.store.transcript_model_storage = [model]
        gbmc.store.transcript_read_ids = {"t1": _reads([(10, 4)], 30)}
        gbmc.correct_transcript_splice_sites()
        self.assertEqual(model.exon_blocks, [(0, 14), (20, 30)])


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
