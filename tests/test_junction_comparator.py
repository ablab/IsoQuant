############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import unittest
from src.common import junctions_from_blocks
from src.junction_comparator import *


class TestJunctionComparator(unittest.TestCase):
    def test_overlaps_at_least(self):
        junction1 = (100, 200)
        junction2 = (150, 250)
        self.assertTrue(overlaps_at_least(junction1, junction2, 10))
        
    def test_equal_ranges_with_delta(self):
        junction1 = (100, 200)
        junction2 = (102, 198)
        self.assertTrue(equal_ranges(junction1, junction2, delta=5))
        self.assertFalse(equal_ranges(junction1, junction2, delta=1))
        
    def test_junctions_from_blocks(self):
        blocks = [(100, 200), (300, 400), (500, 600)]
        junctions = junctions_from_blocks(blocks)
        self.assertEqual(junctions, [(201, 299), (401, 499)])
        
    def test_junctions_from_single_block(self):
        blocks = [(100, 200)]
        junctions = junctions_from_blocks(blocks)
        self.assertEqual(junctions, [])
