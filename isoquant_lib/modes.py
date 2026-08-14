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
    tenX_v2 = 22
    tenX_v3_split = 23
    tenX_v2_split = 24
    curio = 3
    stereoseq_nosplit = 4
    stereoseq = 5
    visium_hd = 61
    visium_5prime = 7
    custom_sc = 10

    def needs_barcode_calling(self):
        return self in [IsoQuantMode.tenX_v3,
                        IsoQuantMode.tenX_v2,
                        IsoQuantMode.tenX_v3_split,
                        IsoQuantMode.tenX_v2_split,
                        IsoQuantMode.curio,
                        IsoQuantMode.stereoseq_nosplit,
                        IsoQuantMode.stereoseq,
                        IsoQuantMode.visium_hd,
                        IsoQuantMode.visium_5prime,
                        IsoQuantMode.custom_sc]

    def needs_pcr_deduplication(self):
        return self in [IsoQuantMode.tenX_v3,
                        IsoQuantMode.tenX_v2,
                        IsoQuantMode.tenX_v3_split,
                        IsoQuantMode.tenX_v2_split,
                        IsoQuantMode.curio,
                        IsoQuantMode.stereoseq_nosplit,
                        IsoQuantMode.stereoseq,
                        IsoQuantMode.visium_hd,
                        IsoQuantMode.visium_5prime,
                        IsoQuantMode.custom_sc]

    def produces_new_fasta(self):
        return self in [IsoQuantMode.stereoseq,
                        IsoQuantMode.tenX_v3_split,
                        IsoQuantMode.tenX_v2_split]

    def needs_barcode_iterator(self):
        return self in [IsoQuantMode.stereoseq_nosplit, IsoQuantMode.stereoseq]

    def supports_cell_barcode_detection(self):
        """Modes whose detectors can extract raw barcodes so cell barcodes can be detected.

        These are the protocols with a large generic whitelist and a small unknown set of
        real cells, where per-read whitelist matching degenerates into exact matching.
        """
        return self in [IsoQuantMode.tenX_v3,
                        IsoQuantMode.tenX_v2,
                        IsoQuantMode.tenX_v3_split,
                        IsoQuantMode.tenX_v2_split,
                        IsoQuantMode.visium_5prime]

    def enforces_single_thread(self):
        return False


ISOQUANT_MODES = [x.name for x in IsoQuantMode]


# Accepted in place of a barcode whitelist file, and as a value for --n_cells
AUTO_BARCODES = "auto"


@unique
class BarcodeCorrectionMethod(Enum):
    """Where the list of cell barcodes reads are matched against comes from."""
    # the supplied whitelist is the list of cell barcodes
    whitelist = 1
    # extract barcodes verbatim first, select the cell barcodes from their read counts,
    # then match reads against those
    graph = 2
    # decided by --n_cells: set means detect, unset means take the whitelist as given
    auto = 3


# Above this whitelist size, matching every read against it effectively requires an exact
# match (see TenXBarcodeDetector.__init__), so taking it as the cell list is worth a warning.
LARGE_WHITELIST_SIZE = 100000


@unique
class AnalysisType(Enum):
    quantification = 1
    exon_quantification = 2
    fusion = 3
    transcript_discovery = 4


# Maps user-facing analysis names (and their short aliases) to AnalysisType.
ANALYSIS_ALIASES = {
    "quantification": AnalysisType.quantification,
    "quant": AnalysisType.quantification,
    "exon_quantification": AnalysisType.exon_quantification,
    "ex_quant": AnalysisType.exon_quantification,
    "fusion": AnalysisType.fusion,
    "transcript_discovery": AnalysisType.transcript_discovery,
    "td": AnalysisType.transcript_discovery,
}

ANALYSIS_CHOICES = list(ANALYSIS_ALIASES.keys())
