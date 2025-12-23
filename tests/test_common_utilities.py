
############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import unittest
import os
import tempfile
from src.common import *


class TestAdditionalCommon(unittest.TestCase):
    def test_proper_plural_form(self):
        self.assertEqual(proper_plural_form("read", 1), "1 read")
        self.assertEqual(proper_plural_form("read", 2), "2 reads")

    def test_rreplace(self):
        self.assertEqual(rreplace("test_string", "string", "replaced"), "test_replaced")
        self.assertEqual(rreplace("test_string_string", "string", "replaced"), "test_string_replaced")
        self.assertEqual(rreplace("no_match", "other", "replaced"), "no_match")
        
    def test_list_to_str(self):
        self.assertEqual(list_to_str([1, 2, 3]), "1,2,3")
        self.assertEqual(list_to_str(["a", "b", "c"]), "a,b,c")
        self.assertEqual(list_to_str([]), ".")
        
    def test_get_path_to_program(self):
        # Test existing program (assuming 'python' exists)
        path = get_path_to_program("python")
        self.assertTrue(os.path.exists(path))
