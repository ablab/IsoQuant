############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import unittest
from isoquant_lib.modes import *


class TestModes(unittest.TestCase):
    def test_needs_barcode_calling(self):
        # Test modes that need barcode calling
        self.assertTrue(IsoQuantMode.tenX_v3.needs_barcode_calling())
        self.assertTrue(IsoQuantMode.tenX_v2.needs_barcode_calling())
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

    def test_supports_molecule_splitting(self):
        # Chemistries with a splitting detector
        self.assertTrue(IsoQuantMode.stereoseq.supports_molecule_splitting())
        self.assertTrue(IsoQuantMode.tenX_v3.supports_molecule_splitting())
        self.assertTrue(IsoQuantMode.tenX_v2.supports_molecule_splitting())
        self.assertTrue(IsoQuantMode.visium_5prime.supports_molecule_splitting())

        # Chemistries without one; asking for --split_molecules true here is an error
        self.assertFalse(IsoQuantMode.bulk.supports_molecule_splitting())
        self.assertFalse(IsoQuantMode.curio.supports_molecule_splitting())
        self.assertFalse(IsoQuantMode.visium_hd.supports_molecule_splitting())
        self.assertFalse(IsoQuantMode.custom_sc.supports_molecule_splitting())

    def test_splitting_detector_exists_for_every_supporting_mode(self):
        """supports_molecule_splitting() must not promise a detector the registry lacks."""
        from isoquant_lib.barcode_calling.detect_barcodes import SPLITTING_BARCODE_CALLING_MODES
        for mode in IsoQuantMode:
            self.assertEqual(mode.supports_molecule_splitting(),
                             mode in SPLITTING_BARCODE_CALLING_MODES,
                             "%s disagrees between the predicate and the registry" % mode.name)

    def test_needs_barcode_iterator(self):
        # Stereo-seq whitelists are too large to hold in memory; unrelated to splitting
        self.assertTrue(IsoQuantMode.stereoseq.needs_barcode_iterator())

        self.assertFalse(IsoQuantMode.bulk.needs_barcode_iterator())
        self.assertFalse(IsoQuantMode.tenX_v3.needs_barcode_iterator())
        self.assertFalse(IsoQuantMode.tenX_v2.needs_barcode_iterator())

    def test_deprecated_aliases_name_real_modes(self):
        for alias, (mode_name, _) in DEPRECATED_MODE_ALIASES.items():
            self.assertIn(mode_name, ISOQUANT_MODES, "%s maps to unknown mode %s" % (alias, mode_name))
            self.assertNotIn(alias, ISOQUANT_MODES, "%s should no longer be a mode" % alias)

    def test_deprecated_aliases_imply_a_possible_split(self):
        """An alias must never imply splitting for a chemistry that cannot split."""
        for alias, (mode_name, splits) in DEPRECATED_MODE_ALIASES.items():
            if splits:
                self.assertTrue(IsoQuantMode[mode_name].supports_molecule_splitting(),
                                "%s implies splitting but %s cannot split" % (alias, mode_name))

    def test_enforces_single_thread(self):
        # All modes should return False for enforces_single_thread
        for mode in IsoQuantMode:
            self.assertFalse(mode.enforces_single_thread())
