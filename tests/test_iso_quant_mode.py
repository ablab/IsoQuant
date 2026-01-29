############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import unittest
from src.modes import *


class TestModes(unittest.TestCase):
    def test_needs_barcode_calling(self):
        # Test modes that need barcode calling
        self.assertTrue(IsoQuantMode.tenX_v3.needs_barcode_calling())
        self.assertTrue(IsoQuantMode.curio.needs_barcode_calling())
        self.assertTrue(IsoQuantMode.stereoseq.needs_barcode_calling())

        # Test mode that doesn't need barcode calling
        self.assertFalse(IsoQuantMode.bulk.needs_barcode_calling())

    def test_needs_pcr_deduplication(self):
        # Test modes that need PCR deduplication
        self.assertTrue(IsoQuantMode.tenX_v3.needs_pcr_deduplication())
        self.assertTrue(IsoQuantMode.visium_hd.needs_pcr_deduplication())
        self.assertTrue(IsoQuantMode.stereoseq.needs_pcr_deduplication())

        # Test mode that doesn't need PCR deduplication
        self.assertFalse(IsoQuantMode.bulk.needs_pcr_deduplication())

    def test_produces_new_fasta(self):
        # Only stereoseq mode produces new fasta
        self.assertTrue(IsoQuantMode.stereoseq.produces_new_fasta())

        # Test other modes
        self.assertFalse(IsoQuantMode.bulk.produces_new_fasta())
        self.assertFalse(IsoQuantMode.tenX_v3.produces_new_fasta())
        self.assertFalse(IsoQuantMode.curio.produces_new_fasta())

    def test_needs_barcode_iterator(self):
        # Test modes that need barcode iterator
        self.assertTrue(IsoQuantMode.stereoseq.needs_barcode_iterator())
        self.assertTrue(IsoQuantMode.stereoseq_nosplit.needs_barcode_iterator())

        # Test modes that don't need barcode iterator
        self.assertFalse(IsoQuantMode.bulk.needs_barcode_iterator())
        self.assertFalse(IsoQuantMode.tenX_v3.needs_barcode_iterator())

    def test_enforces_single_thread(self):
        # All modes should return False for enforces_single_thread
        for mode in IsoQuantMode:
            self.assertFalse(mode.enforces_single_thread())