############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import unittest
import tempfile
import os
from src.stats import EnumStats
from enum import Enum


class TestEnum(Enum):
    value1 = 1
    value2 = 2
    value3 = 3


class TestEnumStats(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.test_file = os.path.join(self.temp_dir, "test_stats.txt")
        
    def tearDown(self):
        if os.path.exists(self.temp_dir):
            import shutil
            shutil.rmtree(self.temp_dir)
    
    def test_add_single_value(self):
        stats = EnumStats()
        stats.add(TestEnum.value1)
        self.assertEqual(stats.stats_dict[TestEnum.value1], 1)
        
    def test_add_multiple_values(self):
        stats = EnumStats()
        stats.add(TestEnum.value1)
        stats.add(TestEnum.value1)
        stats.add(TestEnum.value2)
        self.assertEqual(stats.stats_dict[TestEnum.value1], 2)
        self.assertEqual(stats.stats_dict[TestEnum.value2], 1)
        
    def test_add_with_count(self):
        stats = EnumStats()
        stats.add(TestEnum.value1, 5)
        self.assertEqual(stats.stats_dict[TestEnum.value1], 5)
        
    def test_merge(self):
        stats1 = EnumStats()
        stats1.add(TestEnum.value1, 3)
        stats1.add(TestEnum.value2, 2)
        
        stats2 = EnumStats()
        stats2.add(TestEnum.value1, 2)
        stats2.add(TestEnum.value3, 1)
        
        stats1.merge(stats2)
        self.assertEqual(stats1.stats_dict[TestEnum.value1], 5)
        self.assertEqual(stats1.stats_dict[TestEnum.value2], 2)
        self.assertEqual(stats1.stats_dict[TestEnum.value3], 1)
        
    def test_dump_and_load(self):
        stats = EnumStats()
        stats.add(TestEnum.value1, 3)
        stats.add(TestEnum.value2, 5)
        
        stats.dump(self.test_file)
        self.assertTrue(os.path.exists(self.test_file))
        
        loaded_stats = EnumStats(self.test_file)
        self.assertEqual(loaded_stats.stats_dict[TestEnum.value1], 3)
        self.assertEqual(loaded_stats.stats_dict[TestEnum.value2], 5)
