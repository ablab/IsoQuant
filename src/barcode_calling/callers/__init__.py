############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Barcode detection classes for single-cell and spatial transcriptomics.

Classes:
    BarcodeDetectionResult: Base result class
    CurioBarcodeDetectionResult: Result for Curio (double barcode) platform
    StereoBarcodeDetectionResult: Result for Stereo-seq platform
    TenXBarcodeDetectionResult: Result for 10x Genomics platform
    SplittingBarcodeDetectionResult: Result for read splitting modes
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
"""

from .base import (
    BarcodeDetectionResult,
    LinkerBarcodeDetectionResult,
    TSOBarcodeDetectionResult,
    TenXBarcodeDetectionResult,
    SplittingBarcodeDetectionResult,
    ReadStats,
    increase_if_valid,
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

__all__ = [
    # Result classes
    'BarcodeDetectionResult',
    'LinkerBarcodeDetectionResult',
    'TSOBarcodeDetectionResult',
    'TenXBarcodeDetectionResult',
    'SplittingBarcodeDetectionResult',
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
]
