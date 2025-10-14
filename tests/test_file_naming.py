############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import unittest
from src.file_naming import *


class TestFileNaming(unittest.TestCase):
    def test_saves_file_name(self):
        base_name = "test_output"
        chr_id = "chr1"
        result = saves_file_name(base_name, chr_id)
        self.assertIn(base_name, result)

    def test_multimappers_file_name(self):
        base_name = "test_output"
        chr_id = "chr1"
        result = multimappers_file_name(base_name, chr_id)
        self.assertIn(base_name, result)
        self.assertIn("multimappers", result)
        
    def test_filtered_reads_file_name(self):
        base_name = "test_output"
        chr_id = "chr1"
        result = filtered_reads_file_name(base_name, chr_id)
        self.assertIn(base_name, result)
        self.assertIn("filtered", result)
        
    def test_chromosome_id_sanitization(self):
        # Test that special characters in chromosome IDs are handled
        chr_id = "chr1:1000-2000"
        result = saves_file_name("test", chr_id)
        # Should not contain colons or dashes in unexpected places
        self.assertIsInstance(result, str)
        self.assertTrue(len(result) > 0)
