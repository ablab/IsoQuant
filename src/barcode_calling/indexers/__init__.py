############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
K-mer indexers for fast approximate string matching in barcode calling.

Classes:
    KmerIndexer: Basic dictionary-based indexer
    ArrayKmerIndexer: Array-based indexer for small k values
    Dict2BitKmerIndexer: Memory-efficient dict with integer keys
    Array2BitKmerIndexer: Memory-efficient array for large barcode sets
    SharedMemoryArray2BitKmerIndexer: Shared memory version for parallel processing
"""

from .base import KmerIndexer, ArrayKmerIndexer
from .two_bit import Dict2BitKmerIndexer, Array2BitKmerIndexer
from .shared_memory import SharedMemoryArray2BitKmerIndexer, SharedMemoryIndexInfo

__all__ = [
    'KmerIndexer',
    'ArrayKmerIndexer',
    'Dict2BitKmerIndexer',
    'Array2BitKmerIndexer',
    'SharedMemoryArray2BitKmerIndexer',
    'SharedMemoryIndexInfo',
]
