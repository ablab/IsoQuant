############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from enum import Enum, unique


@unique
class IsoQuantMode(Enum):
    """One member per chemistry. Whether molecules are split is --split_molecules."""
    bulk = 1
    tenX_v3 = 2
    tenX_v2 = 22
    curio = 3
    stereoseq = 5
    visium_hd = 61
    visium_5prime = 7
    custom_sc = 10

    def needs_barcode_calling(self):
        return self in [IsoQuantMode.tenX_v3,
                        IsoQuantMode.tenX_v2,
                        IsoQuantMode.curio,
                        IsoQuantMode.stereoseq,
                        IsoQuantMode.visium_hd,
                        IsoQuantMode.visium_5prime,
                        IsoQuantMode.custom_sc]

    def needs_pcr_deduplication(self):
        return self in [IsoQuantMode.tenX_v3,
                        IsoQuantMode.tenX_v2,
                        IsoQuantMode.curio,
                        IsoQuantMode.stereoseq,
                        IsoQuantMode.visium_hd,
                        IsoQuantMode.visium_5prime,
                        IsoQuantMode.custom_sc]

    def supports_molecule_splitting(self):
        """Chemistries with a detector that can find several cDNA molecules in one read.

        Concatenated molecules arise when cDNAs are ligated end to end during library
        preparation. Asking for splitting on any other mode is an error, not a no-op.
        """
        return self in [IsoQuantMode.tenX_v3,
                        IsoQuantMode.tenX_v2,
                        IsoQuantMode.stereoseq,
                        IsoQuantMode.visium_5prime]

    def needs_barcode_iterator(self):
        # a property of the Stereo-seq whitelist size, unrelated to splitting
        return self in [IsoQuantMode.stereoseq]

    def supports_cell_barcode_detection(self):
        """Modes whose detectors can extract raw barcodes so cell barcodes can be detected.

        These are the protocols with a large generic whitelist and a small unknown set of
        real cells, where per-read whitelist matching degenerates into exact matching.
        """
        return self in [IsoQuantMode.tenX_v3,
                        IsoQuantMode.tenX_v2,
                        IsoQuantMode.visium_5prime]

    def enforces_single_thread(self):
        return False


ISOQUANT_MODES = [x.name for x in IsoQuantMode]

# Values accepted by --split_molecules. AUTO splits wherever the chemistry supports it and
# does nothing where it does not; TRUE additionally fails on a mode that cannot split, so a
# user asking for the impossible gets an error rather than silently unsplit results.
SPLIT_MOLECULES_TRUE = "true"
SPLIT_MOLECULES_FALSE = "false"
SPLIT_MOLECULES_AUTO = "auto"
SPLIT_MOLECULES_CHOICES = [SPLIT_MOLECULES_TRUE, SPLIT_MOLECULES_FALSE, SPLIT_MOLECULES_AUTO]

# Superseded mode names -> (chemistry, whether they implied splitting). The split value is
# only a default; an explicit --split_molecules wins.
DEPRECATED_MODE_ALIASES = {
    "tenX_v3_split": ("tenX_v3", True),
    "tenX_v2_split": ("tenX_v2", True),
    "stereoseq_nosplit": ("stereoseq", False),
}


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
