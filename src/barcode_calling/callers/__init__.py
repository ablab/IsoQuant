############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Barcode detection classes for single-cell and spatial transcriptomics.

Protocol:
    BarcodeResult: Protocol interface for all detection results

Classes:
    BarcodeDetectionResult: Base result class
    LinkerBarcodeDetectionResult: Result for platforms with linker sequences
    TSOBarcodeDetectionResult: Result for Stereo-seq with TSO detection
    TenXBarcodeDetectionResult: Result for 10x Genomics platform
    SplittingBarcodeDetectionResult: Result for read splitting modes
    ExtractionResult: Dict-based result for universal molecule extraction
    ReadStats: Statistics tracker for barcode detection

Detectors:
    CurioBarcodeDetector: Curio platform detector
    CurioIlluminaDetector: Curio detector for Illumina reads
    StereoBarcodeDetector: Stereo-seq detector
    SharedMemoryStereoBarcodeDetector: Stereo-seq with shared memory
    StereoSplittingBarcodeDetector: Stereo-seq with read splitting
    SharedMemoryStereoSplittingBarcodeDetector: Stereo-seq splitting with shared memory
    SharedMemoryWrapper: Generic shared memory wrapper
    TenXBarcodeDetector: 10x Genomics v3 detector
    VisiumHDBarcodeDetector: Visium HD detector
    UniversalSingleMoleculeExtractor: Universal barcode detector for custom molecules
"""

from .protocol import BarcodeResult

from .base import (
    BarcodeDetectionResult,
    LinkerBarcodeDetectionResult,
    TSOBarcodeDetectionResult,
    TenXBarcodeDetectionResult,
    SplittingBarcodeDetectionResult,
    ReadStats,
    increase_if_valid,
)

from .extraction_result import (
    ExtractionResult,
    DetectedElement,
)

from .curio import (
    CurioBarcodeDetector,
    CurioIlluminaDetector,
)

from .stereo import (
    StereoBarcodeDetector,
    SharedMemoryStereoBarcodeDetector,
    StereoSplittingBarcodeDetector,
    SharedMemoryStereoSplittingBarcodeDetector,
    SharedMemoryWrapper,
)

from .tenx import (
    TenXBarcodeDetector,
    VisiumHDBarcodeDetector,
)

from .molecule_structure import MoleculeStructure

from .universal_extraction import UniversalSingleMoleculeExtractor

__all__ = [
    # Protocol
    'BarcodeResult',
    # Result classes
    'BarcodeDetectionResult',
    'LinkerBarcodeDetectionResult',
    'TSOBarcodeDetectionResult',
    'TenXBarcodeDetectionResult',
    'SplittingBarcodeDetectionResult',
    'ExtractionResult',
    'DetectedElement',
    'ReadStats',
    'increase_if_valid',
    # Curio detectors
    'CurioBarcodeDetector',
    'CurioIlluminaDetector',
    # Stereo detectors
    'StereoBarcodeDetector',
    'SharedMemoryStereoBarcodeDetector',
    'StereoSplittingBarcodeDetector',
    'SharedMemoryStereoSplittingBarcodeDetector',
    'SharedMemoryWrapper',
    # 10x detectors
    'TenXBarcodeDetector',
    'VisiumHDBarcodeDetector',
    # Universal extraction
    'MoleculeStructure',
    'UniversalSingleMoleculeExtractor',
]
