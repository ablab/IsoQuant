############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import unittest
from src.id_policy import *


class TestSimpleIDDistributor(unittest.TestCase):
    def test_get_next_id(self):
        distributor = SimpleIDDistributor()
        self.assertEqual(distributor.increment(), 1)
        self.assertEqual(distributor.increment(), 2)
        self.assertEqual(distributor.increment(), 3)

