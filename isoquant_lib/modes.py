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

    def supports_graph_correction(self):
        """Modes whose detectors can extract raw barcodes for graph-based correction.

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


@unique
class BarcodeCorrectionMethod(Enum):
    """How extracted barcodes are mapped onto real ones."""
    # match every read against the whitelist during extraction
    whitelist = 1
    # extract barcodes verbatim, then select actually sequenced barcodes from their counts
    # and correct the rest onto them via an edit-distance graph
    graph = 2
    # graph when the whitelist is large enough that per-read matching degenerates
    # into exact matching, whitelist otherwise
    auto = 3


# Above this whitelist size per-read matching effectively requires an exact match
# (see TenXBarcodeDetector.__init__), so graph-based correction takes over in auto mode.
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
