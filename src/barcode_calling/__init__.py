############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Barcode calling for single-cell and spatial transcriptomics.

This package provides k-mer indexing and barcode detection for various platforms:
- Curio (double barcode)
- Stereo-seq (standard and splitting modes)
- 10x Genomics (v3, Visium HD)

Subpackages:
    indexers: K-mer indexing implementations
    callers: Barcode detector implementations
"""

# Indexers
from .indexers import (
    KmerIndexer,
    ArrayKmerIndexer,
    Dict2BitKmerIndexer,
    Array2BitKmerIndexer,
    SharedMemoryArray2BitKmerIndexer,
    SharedMemoryIndexInfo,
)

# Result classes
from .callers import (
    BarcodeDetectionResult,
    LinkerBarcodeDetectionResult,
    TSOBarcodeDetectionResult,
    TenXBarcodeDetectionResult,
    SplittingBarcodeDetectionResult,
    ReadStats,
    increase_if_valid,
)

# Detector classes
from .callers import (
    CurioBarcodeDetector,
    CurioIlluminaDetector,
    StereoBarcodeDetector,
    SharedMemoryStereoBarcodeDetector,
    StereoSplittingBarcodeDetector,
    SharedMemoryStereoSplittingBarcodeDetector,
    SharedMemoryWrapper,
    TenXBarcodeDetector,
    VisiumHDBarcodeDetector,
    UniversalSingleMoleculeExtractor,
    MoleculeStructure
)

# Utilities
from .common import str_to_2bit, bit_to_str, find_polyt_start, batch_str_to_2bit, batch_str_to_2bit_chunked

__all__ = [
    # Indexers
    'KmerIndexer',
    'ArrayKmerIndexer',
    'Dict2BitKmerIndexer',
    'Array2BitKmerIndexer',
    'SharedMemoryArray2BitKmerIndexer',
    'SharedMemoryIndexInfo',
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
    # universal calling
    'UniversalSingleMoleculeExtractor',
    'MoleculeStructure',
    # Utilities
    'str_to_2bit',
    'bit_to_str',
    'find_polyt_start',
    'batch_str_to_2bit',
    'batch_str_to_2bit_chunked',
]
