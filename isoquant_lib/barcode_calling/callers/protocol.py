############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Protocol definition for barcode detection results.

This module defines the BarcodeResult protocol that all barcode detection
result classes implement. Using Protocol (structural typing) allows
existing classes to automatically conform without explicit inheritance.
"""

from typing import Protocol, List, ClassVar, runtime_checkable


@runtime_checkable
class BarcodeResult(Protocol):
    """
    Protocol defining the interface for barcode detection results.

    All barcode detection result classes (BarcodeDetectionResult,
    LinkerBarcodeDetectionResult, TSOBarcodeDetectionResult,
    TenXBarcodeDetectionResult, ExtractionResult) implement this interface.

    Using Protocol (PEP 544) enables structural subtyping - classes
    don't need to explicitly inherit from this protocol to be
    considered compatible; they just need to implement the required
    methods and attributes.

    Attributes:
        NOSEQ: Class variable sentinel for missing/undetected sequence ("*").
        read_id: The unique identifier of the read being analyzed.
        strand: The detected strand orientation ('+', '-', or '.').

    Example:
        def process_result(result: BarcodeResult) -> None:
            if result.is_valid():
                barcode = result.get_barcode()
                umi = result.get_umi()
                print(f"{result.read_id}: {barcode}/{umi}")
    """

    NOSEQ: ClassVar[str]
    read_id: str
    strand: str

    def get_barcode(self) -> str:
        """
        Get the detected barcode sequence.

        Returns the barcode sequence if detected, or NOSEQ ("*") if no
        barcode was found. For platforms with multiple barcode segments
        (e.g., Visium HD with part1|part2), returns the concatenated
        barcode with "|" delimiter.

        Returns:
            str: Detected barcode sequence or NOSEQ if not found.
        """
        ...

    def get_umi(self) -> str:
        """
        Get the detected UMI (Unique Molecular Identifier) sequence.

        Returns the UMI sequence if detected, or NOSEQ ("*") if no
        UMI was found. For platforms with multiple UMI elements,
        returns concatenated UMIs with "|" delimiter.

        Returns:
            str: Detected UMI sequence or NOSEQ if not found.
        """
        ...

    def is_valid(self) -> bool:
        """
        Check if a valid barcode was detected.

        A result is considered valid if get_barcode() returns something
        other than NOSEQ. This is the primary method to check if barcode
        detection was successful.

        Returns:
            bool: True if a valid barcode was detected, False otherwise.
        """
        ...

    def has_barcode(self) -> bool:
        """
        Check if a barcode was detected (alias for is_valid()).

        Used by ReadStats for counting detected barcodes.

        Returns:
            bool: True if barcode was detected.
        """
        ...

    def has_umi(self) -> bool:
        """
        Check if a valid UMI was detected.

        Returns True if get_umi() returns something other than NOSEQ.
        Used by ReadStats for counting detected UMIs.

        Returns:
            bool: True if UMI was detected.
        """
        ...

    def set_strand(self, strand: str) -> None:
        """
        Set the detected strand orientation.

        Args:
            strand: Strand orientation. Should be '+' for forward,
                   '-' for reverse, or '.' for unknown/unstranded.
        """
        ...

    def update_coordinates(self, delta: int) -> None:
        """
        Shift all genomic coordinates by a delta value.

        Used when processing read subsequences (e.g., after splitting
        a concatenated read). Implementations should update all position
        attributes (polyT, primer, linker_start, linker_end, etc.)
        by adding delta. Coordinates that are -1 (invalid/not detected)
        should remain -1.

        Args:
            delta: Amount to add to all coordinate values.
        """
        ...

    def more_informative_than(self, other: 'BarcodeResult') -> bool:
        """
        Compare two results to determine which is more informative.

        Used when both forward and reverse complement orientations
        produce results, to select the better one. The comparison
        criteria vary by platform but typically consider:
        - Barcode alignment score
        - Presence of additional features (polyT, primer, linker)
        - Feature positions

        Args:
            other: Another detection result to compare against.

        Returns:
            bool: True if this result is more informative than `other`.
        """
        ...

    def get_additional_attributes(self) -> List[str]:
        """
        Get list of detected additional features.

        Returns human-readable strings describing which optional
        features were detected (e.g., "PolyT detected", "Primer detected",
        "Linker detected", "TSO detected", "R1 detected").

        Used for statistics collection and quality reporting.

        Returns:
            List[str]: Names of detected features. Empty list if none.
        """
        ...

    def __str__(self) -> str:
        """
        Format result as a tab-separated string for TSV output.

        The format should match what header() returns.

        Returns:
            str: Tab-separated values for TSV output.
        """
        ...

    def header(self) -> str:
        """
        Get TSV header line for result output.

        Returns the column names matching the fields output by __str__().

        Returns:
            str: Tab-separated column names with leading '#'.
        """
        ...
