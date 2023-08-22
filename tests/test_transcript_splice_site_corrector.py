from unittest import TestCase
from unittest import main as unittest_main


from src.transcript_splice_site_corrector import (
    extract_location_from_cigar_string,
    count_deletions_from_cigar_codes_in_given_window,
    extract_splice_site_locations_within_aligned_read,
    count_deletions_for_splice_site_locations,
    compute_most_common_case_of_deletions,
    extract_nucleotides_from_most_common_del_location,
    compute_most_common_del_and_verify_nucleotides,
    threshold_exceeded,
    sublist_largest_values_exists,
    correct_splice_site_errors,
    generate_updated_exon_list,
)

class TestMoreConservativeStrategyConditions(TestCase):
    
    def test_threshold_exceeds_returns_true(self):
        THRESHOLD = 0.7
        del_pos_distr = [0, 0, 10, 10, 10, 10, 0, 0]
        deletions = {4: 10}
        most_common_del = 4
        result = threshold_exceeded(
            del_pos_distr,
            deletions,
            most_common_del,
            THRESHOLD)
        self.assertTrue(result)

    def test_threshold_not_exceeded_returns_false(self):
        THRESHOLD = 0.7
        del_pos_distr = [0, 0, 10, 10, 10, 6, 0, 0]
        deletions = {4: 6, 3: 4}
        most_common_del = 4
        result = threshold_exceeded(
            del_pos_distr,
            deletions,
            most_common_del,
            THRESHOLD)
        self.assertFalse(result)

    def test_sublist_largest_values_exists_returns_true(self):
        lst = [0, 0, 10, 10, 10, 10, 0, 0]
        n = 4
        result = sublist_largest_values_exists(lst, n)
        self.assertTrue(result)

    def test_sublist_largest_values_exists_returns_false(self):
        lst = [0, 0, 10, 10, 10, 0, 6, 0]
        n = 4
        result = sublist_largest_values_exists(lst, n)
        self.assertFalse(result)
        

class TestExtractingLocationFromCigarString(TestCase):

    def test_cigar_string_with_soft_clip_and_one_match_is_parsed_correctly(self):
        cigar = [(4, 50), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 105
        expected_output = 55
        result = extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)
        

    def test_cigar_string_with_soft_clip_insertion_and_one_match_is_parsed_correctly(self):
        cigar = [(4, 50), (1, 10), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 105
        expected_output = 65
        result = extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)
        

    def test_cigar_str_with_s_d_i_m_gives_correct_output(self):
        cigar = [(4, 50), (2, 10), (1, 10), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 115
        expected_output = 75
        result = extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_cigar_str_with_s_d_n_m_gives_correct_output(self):
        cigar = [(4, 50), (2, 10), (3, 100), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 215
        expected_output = 165
        result = extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_cigar_str_with_s_m_i_n_m_gives_correct_output(self):
        cigar = [(4, 50), (0, 10), (1, 10), (3, 100), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 215
        expected_output = 175
        result = extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_location_outside_of_cigar_str_returns_minus_one(self):
        cigar = [(4, 50), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 199
        expected_output = -1
        result = extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_more_complicated_test_returns_correct_position(self):
        cigar_tuples = [(4, 156), (0, 12), (2, 3), (0, 2), (2, 2), (0, 10), (2, 2), (0, 4), (2, 3), (0, 7), (1, 1), (0, 16), (1, 4), (0, 23), (1, 1), (0, 7),
                        (1, 1), (0, 9), (2, 1), (0, 13), (2, 1), (0, 15), (2, 2), (0, 3), (1, 2), (0, 19), (2, 2), (0, 20), (2, 1), (0, 32), (3, 294), (0, 36), (4, 25)]
        reference_start = 72822568
        reference_end = 73822568
        position = 72823071
        expected_output = 668
        result = extract_location_from_cigar_string(
            cigar_tuples, reference_start, reference_end, position)
        self.assertEqual(result, expected_output)

    def test_case_that_does_not_consume_any_reference_returns_the_correct_location(self):
        cigar = [(4, 50), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 100
        expected_output = 50
        result = extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_case_that_has_no_reference_consuming_codes_returns_minus_one_as_error(self):
        cigar = [(4, 50), (1, 10)]
        reference_start = 100
        reference_end = 160
        location = 100
        expected_output = -1
        result = extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_case_that_has_no_reference_consuming_codes_at_the_end_returns_minus_one_as_error(self):
        cigar = [(4, 50), (0, 10), (1, 10)]
        reference_start = 100
        reference_end = 160
        location = 110
        expected_output = -1
        result = extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_case_that_has_it_s_location_at_final_match_returns_correct_value(self):
        cigar = [(4, 50), (0, 10), (1, 10)]
        reference_start = 100
        reference_end = 110
        location = 110
        expected_output = 60
        result = extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)


class TestIndelCountingFromCigarCodes(TestCase):

    def setUp(self):
        self.window_size = 8

    def test_indel_counter_returns_false_and_an_empty_debug_list_for_given_empty_list(self):
        cigar_tuples = []
        aligned_location = 100
        location_is_end = False
        splice_site_data = {
            'deletions': {},
            "del_pos_distr": [0] * self.window_size,
        }
        expected_result = {
            'deletions': {0: 1},
            "del_pos_distr": [0, 0, 0, 0, 0, 0, 0, 0]
        }
        count_deletions_from_cigar_codes_in_given_window(
            cigar_tuples,
            aligned_location,
            location_is_end,
            splice_site_data,
            self.window_size)
        
        self.assertEqual(splice_site_data['deletions'], expected_result['deletions'])
        self.assertEqual(splice_site_data['del_pos_distr'], expected_result['del_pos_distr'])
        

        
    def test_indels_are_counted_correctly(self):
        cigar_tuples = [(0, 20), (2, 3), (1, 2), (0, 10)]
        aligned_location = 27
        location_is_end = True
        splice_site_data = {
            'deletions': {},
            "del_pos_distr": [0] * self.window_size,
        }


        expected_result = {
            'deletions': {3: 1},
            "del_pos_distr": [1, 1, 1, 0, 0, 0, 0, 0]
        }

        count_deletions_from_cigar_codes_in_given_window(
            cigar_tuples,
            aligned_location,
            location_is_end,
            splice_site_data,
            self.window_size)

        self.assertEqual(splice_site_data['deletions'], expected_result['deletions'])
        self.assertEqual(splice_site_data['del_pos_distr'], expected_result['del_pos_distr'])

    def test_full_window_of_dels_returns_true_for_errors(self):
        cigar_tuples = [(0, 20), (2, 8), (1, 2), (0, 10)]
        aligned_location = 20
        location_is_end = False
        splice_site_data = {
            'deletions': {},
            "del_pos_distr": [0] * self.window_size,
        }
        expected_result = {
            'deletions': {8: 1},
            "del_pos_distr": [1, 1, 1, 1, 1, 1, 1, 1]
        }

        count_deletions_from_cigar_codes_in_given_window(
            cigar_tuples,
            aligned_location,
            location_is_end,
            splice_site_data,
            self.window_size)

        self.assertEqual(splice_site_data['deletions'], expected_result['deletions'])
        self.assertEqual(splice_site_data['del_pos_distr'], expected_result['del_pos_distr'])

class TestExtractSpliceSiteLocationsFromAlignedRead(TestCase):

    def test_correct_splice_sites_are_extracted(self):
        exons = [(1, 10), (20, 30), (40, 50)]
        read_start = 20
        read_end = 40
        result = extract_splice_site_locations_within_aligned_read(
            read_start, read_end, exons)
        expected_output = [(20, False), (30, True) , (40, False)]
        self.assertEqual(result, expected_output)


class TestExonListUpdater(TestCase):
    
    def test_error_at_location_start_is_corrected(self):
        exons = [(1, 10), (20, 30), (40, 50)]
        locations_with_errors = [20]
        splice_site_cases = {
            20: {
                "most_common_del": 4,
            }
        }
        result = generate_updated_exon_list(
            splice_site_cases, locations_with_errors, exons)
        expected_result = [(1, 10), (24, 30), (40, 50)]
        self.assertEqual(result, expected_result)
    
    def test_error_at_location_end_is_corrected(self):
        exons = [(1, 10), (20, 30), (40, 50)]
        locations_with_errors = [30]
        splice_site_cases = {
            30: {
                "most_common_del": -4,
            }
        }
        result = generate_updated_exon_list(
            splice_site_cases, locations_with_errors, exons)
        expected_result = [(1, 10), (20, 26), (40, 50)]
        self.assertEqual(result, expected_result)
        

    pass

class TestHelperFunctions(TestCase):

    def test_distinct_most_common_case_is_returned_for_location_end(self):
        cases = {0: 10, 1: 2, 3: 0, 4: 20, 5: 1}
        location_is_end = False
        result = compute_most_common_case_of_deletions(cases, location_is_end)
        expected_result = 4
        self.assertEqual(result, expected_result)

    def test_distinct_most_common_case_is_returned_for_location_start(self):
        cases = {0: 10, 1: 2, 3: 0, 4: 20, 5: 1}
        location_is_end = True
        result = compute_most_common_case_of_deletions(cases, location_is_end)
        expected_result = -4
        self.assertEqual(result, expected_result)
    
    def test_if_no_distinct_most_commont_del_exists_return_neg_one(self):
        cases = {0: 10, 1: 2, 3: 20, 4: 20, 5: 1}
        location_is_end = False
        result = compute_most_common_case_of_deletions(cases, location_is_end)
        expected_result = -1
        self.assertEqual(result, expected_result)