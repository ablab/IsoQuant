############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from enum import Enum, unique


@unique
class IsoQuantMode(Enum):
    bulk = 1
    tenX_v3 = 2
    curio = 3
    stereoseq_nosplit = 4
    stereoseq = 5
    visium_hd = 6
    visium_5prime = 7
    custom_sc = 10

    def needs_barcode_calling(self):
        return self in [IsoQuantMode.tenX_v3,
                        IsoQuantMode.curio,
                        IsoQuantMode.stereoseq_nosplit,
                        IsoQuantMode.stereoseq,
                        IsoQuantMode.visium_hd,
                        IsoQuantMode.visium_5prime,
                        IsoQuantMode.custom_sc]

    def needs_pcr_deduplication(self):
        return self in [IsoQuantMode.tenX_v3,
                        IsoQuantMode.curio,
                        IsoQuantMode.stereoseq_nosplit,
                        IsoQuantMode.stereoseq,
                        IsoQuantMode.visium_hd,
                        IsoQuantMode.visium_5prime,
                        IsoQuantMode.custom_sc]

    def produces_new_fasta(self):
        return self in [IsoQuantMode.stereoseq]

    def needs_barcode_iterator(self):
        return self in [IsoQuantMode.stereoseq_nosplit, IsoQuantMode.stereoseq]

    def enforces_single_thread(self):
        return False 


ISOQUANT_MODES = [x.name for x in IsoQuantMode]