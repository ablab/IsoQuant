############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Dict-based flexible detection result for universal molecule extraction.

Unlike the attribute-based result classes in base.py, ExtractionResult stores
detected elements in a dictionary, allowing flexible molecule structures.
Implements the BarcodeResult protocol.
"""

import logging
from collections import defaultdict
from typing import List

from .molecule_structure import MoleculeStructure, ElementType

logger = logging.getLogger('IsoQuant')


class DetectedElement:
    """Container for a single detected element's coordinates and sequence."""

    def __init__(self, start: int = -1, end: int = -1, score: int = -1, seq: str = None):
        self.start = start
        self.end = end
        self.score = score
        self.seq = seq


class ExtractionResult:
    """
    Dict-based flexible detection result for universal molecule extraction.

    Unlike the attribute-based result classes (BarcodeDetectionResult, etc.),
    ExtractionResult stores detected elements in a dictionary, allowing
    flexible molecule structures defined by MoleculeStructure.

    Implements the BarcodeResult protocol.

    Barcode/UMI elements are identified by naming convention:
    - Elements with prefix "barcode" (case-insensitive) are barcodes
    - Elements with prefix "umi" (case-insensitive) are UMIs
    """

    NOSEQ = "*"  # Sentinel for missing/undetected sequence (matches BarcodeDetectionResult.NOSEQ)

    def __init__(self, molecule_structure: MoleculeStructure, read_id: str, strand: str, detected_results: dict):
        """
        Initialize extraction result.

        Args:
            molecule_structure: Molecule structure definition
            read_id: Read identifier
            strand: Detected strand ('+', '-', or '.')
            detected_results: Dict mapping element names to DetectedElement objects
        """
        self.molecule_structure = molecule_structure
        self.read_id = read_id
        self.strand = strand
        self.detected_results = detected_results
        # Cache barcode/UMI element names for efficiency
        self._barcode_elements: List[str] = []
        self._umi_elements: List[str] = []
        self._identify_barcode_umi_elements()

    def _identify_barcode_umi_elements(self) -> None:
        """Identify barcode/UMI elements by naming convention."""
        for el in self.molecule_structure:
            name_lower = el.element_name.lower()
            if name_lower.startswith("barcode"):
                self._barcode_elements.append(el.element_name)
            elif name_lower.startswith("umi"):
                self._umi_elements.append(el.element_name)

    def get_barcode(self) -> str:
        """
        Get the detected barcode sequence.

        For molecule structures with multiple barcode elements (e.g., barcode1, barcode2),
        returns concatenated barcodes using '|' delimiter.

        Returns:
            str: Concatenated barcode sequence or NOSEQ if not detected.
        """
        parts = []
        for el_name in self._barcode_elements:
            if el_name in self.detected_results:
                detected = self.detected_results[el_name]
                if detected.seq and detected.seq != self.NOSEQ:
                    parts.append(detected.seq)
        return "|".join(parts) if parts else self.NOSEQ

    def get_umi(self) -> str:
        """
        Get the detected UMI sequence.

        For molecule structures with multiple UMI elements,
        returns concatenated UMIs using '|' delimiter.

        Returns:
            str: Concatenated UMI sequence or NOSEQ if not detected.
        """
        parts = []
        for el_name in self._umi_elements:
            if el_name in self.detected_results:
                detected = self.detected_results[el_name]
                if detected.seq and detected.seq != self.NOSEQ:
                    parts.append(detected.seq)
        return "|".join(parts) if parts else self.NOSEQ

    def is_valid(self) -> bool:
        """Check if a valid barcode was detected."""
        return self.get_barcode() != self.NOSEQ

    def update_coordinates(self, delta: int) -> None:
        """
        Shift all genomic coordinates by delta.

        Used when processing read subsequences.

        Args:
            delta: Amount to shift coordinates
        """
        for detected in self.detected_results.values():
            if detected.start != -1:
                detected.start += delta
            if detected.end != -1:
                detected.end += delta

    def more_informative_than(self, other: 'ExtractionResult') -> bool:
        """
        Compare two results to determine which is more informative.

        Compares by number of detected elements, then total score.

        Args:
            other: Another detection result

        Returns:
            True if this result is more informative than other
        """
        self_count = sum(1 for d in self.detected_results.values() if d.start != -1)
        other_count = sum(1 for d in other.detected_results.values() if d.start != -1)
        if self_count != other_count:
            return self_count > other_count
        self_score = sum(d.score for d in self.detected_results.values() if d.score > 0)
        other_score = sum(d.score for d in other.detected_results.values() if d.score > 0)
        return self_score > other_score

    def get_additional_attributes(self) -> List[str]:
        """
        Get list of detected additional features.

        Returns:
            List of detected element names with " detected" suffix.
        """
        return ["%s detected" % el_name for el_name, d in self.detected_results.items() if d.start != -1]

    def set_strand(self, strand: str) -> None:
        """Set the detected strand."""
        self.strand = strand

    def __str__(self) -> str:
        """Format result as TSV line."""
        res_str = "%s\t%s" % (self.read_id, self.strand)
        for el in self.molecule_structure:
            if el.element_type == ElementType.cDNA:
                continue
            if el.element_type == ElementType.PolyT:
                if ElementType.PolyT.name not in self.detected_results:
                    res_str += "\t-1\t-1"
                    continue
                detected_element = self.detected_results[ElementType.PolyT.name]
                res_str += "\t%d\t%d" % (detected_element.start, detected_element.end)
            elif el.element_type.needs_only_coordinates():
                if el.element_name not in self.detected_results:
                    res_str += "\t-1\t-1"
                    continue
                detected_element = self.detected_results[el.element_name]
                res_str += "\t%d\t%d" % (detected_element.start, detected_element.end)
            elif el.element_type.needs_sequence_extraction() or el.element_type.needs_correction():
                if el.element_name not in self.detected_results:
                    res_str += "\t-1\t-1\t%s" % self.NOSEQ
                    continue
                detected_element = self.detected_results[el.element_name]
                seq = detected_element.seq if detected_element.seq else self.NOSEQ
                res_str += "\t%d\t%d\t%s" % (detected_element.start, detected_element.end, seq)
        return res_str

    def header(self) -> str:
        """Get TSV header for result output."""
        return self.molecule_structure.header()

    # TODO: add simple serialization


class ReadStats:
    def __init__(self):
        self.read_count = 0
        self.pattern_counts = defaultdict(int)

    def add_read(self, barcode_detection_result):
        self.read_count += 1
        for el in barcode_detection_result.detected_results:
            if barcode_detection_result.detected_results[el].start != -1:
                self.pattern_counts[el] += 1

    def __str__(self):
        human_readable_str =  ("Total reads:\t%d\n" % self.read_count)
        for a in self.pattern_counts:
            human_readable_str += "%s:\t%d\n" % (a, self.pattern_counts[a])
        return human_readable_str
