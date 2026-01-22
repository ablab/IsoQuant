
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import unittest
import os
import shutil
import tempfile
from src.file_utils import *


class TestFileUtils(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.chr_ids = ["chr1", "chr2", "chr3"]
        
    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
            
    def test_merge_file_list(self):
        fname = "test_file.txt"
        label = "file"
        expected = ["test_file_chr1.txt", "test_file_chr2.txt", "test_file_chr3.txt"]
        result = merge_file_list(fname, label, self.chr_ids)
        self.assertEqual(result, expected)

    def test_merge_files(self):
        # Create test files
        test_files = []
        contents = ["#header\nchr1_content\n", "#header\nchr2_content\n", "#header\nchr3_content\n"]
        for i, chr_id in enumerate(self.chr_ids):
            fname = os.path.join(self.test_dir, f"test_file_{chr_id}.txt")
            with open(fname, 'w') as f:
                f.write(contents[i])
            test_files.append(fname)
            
        # Test merging with header
        merged_file = os.path.join(self.test_dir, "test_file.txt")
        with open(merged_file, 'w') as f:
            merge_files(merged_file, "file", self.chr_ids, f, copy_header=True)
            
        with open(merged_file) as f:
            content = f.read()
        expected = "#header\nchr1_content\nchr2_content\nchr3_content\n"
        self.assertEqual(content, expected)
        
    def test_normalize_path(self):
        config_path = "/path/to/config/config.txt"
        rel_path = "data/file.txt"
        abs_path = "/absolute/path/file.txt"
        
        # Test relative path
        expected_rel = "/path/to/config/data/file.txt"
        self.assertEqual(normalize_path(config_path, rel_path), expected_rel)
        
        # Test absolute path
        self.assertEqual(normalize_path(config_path, abs_path), abs_path)
