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

    @property
    def _barcode_elements(self) -> List[str]:
        """Get barcode element names from molecule structure."""
        return self.molecule_structure.barcode_elements

    @property
    def _umi_elements(self) -> List[str]:
        """Get UMI element names from molecule structure."""
        return self.molecule_structure.umi_elements

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

    def has_barcode(self) -> bool:
        """Check if a barcode was detected (alias for is_valid())."""
        return self.get_barcode() != self.NOSEQ

    def has_umi(self) -> bool:
        """Check if a valid UMI was detected."""
        return self.get_umi() != self.NOSEQ

    def get_barcode_score(self) -> int:
        """
        Get the best barcode alignment score.

        Returns the maximum score from all barcode elements, or -1 if none detected.
        """
        max_score = -1
        for el_name in self._barcode_elements:
            if el_name in self.detected_results:
                detected = self.detected_results[el_name]
                if detected.score > max_score:
                    max_score = detected.score
        return max_score

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
        Get list of detected supplementary features (excluding barcode/UMI).

        Barcode and UMI elements are excluded because they are already
        tracked by ReadStats.bc_count and ReadStats.umi_count respectively.

        Returns:
            List of detected element names with " detected" suffix.
        """
        barcode_umi_set = set(self._barcode_elements + self._umi_elements)
        return ["%s detected" % el_name for el_name, d in self.detected_results.items()
                if d.start != -1 and el_name not in barcode_umi_set]

    def set_strand(self, strand: str) -> None:
        """Set the detected strand."""
        self.strand = strand

    def _get_supplementary_elements(self) -> List[str]:
        """Get list of element names that are not barcode, UMI, or concatenated/duplicated parts."""
        supplementary = []
        barcode_umi_set = set(self._barcode_elements + self._umi_elements)
        for el in self.molecule_structure:
            if el.element_type == ElementType.cDNA:
                continue
            if el.element_name in self.molecule_structure.elements_to_concatenate:
                continue
            if el.element_name in self.molecule_structure.duplicated_elements:
                continue
            if el.element_name not in barcode_umi_set:
                supplementary.append(el.element_name)
        return supplementary

    def __str__(self) -> str:
        """
        Format result as TSV line in standard format.

        Format: read_id, barcode, UMI, BC_score, valid_UMI, strand, [supplementary elements]
        This matches the format expected by the rest of the IsoQuant pipeline.
        """
        # Standard fields: read_id, barcode, UMI, BC_score, valid_UMI, strand
        res_str = "%s\t%s\t%s\t%d\t%s\t%s" % (
            self.read_id,
            self.get_barcode(),
            self.get_umi(),
            self.get_barcode_score(),
            self.has_umi(),
            self.strand
        )

        # Supplementary elements (non-barcode, non-UMI)
        barcode_umi_set = set(self._barcode_elements + self._umi_elements)
        for el in self.molecule_structure:
            if el.element_type == ElementType.cDNA:
                continue
            # Skip concatenated/duplicated parts (data lives under base name)
            if el.element_name in self.molecule_structure.elements_to_concatenate:
                continue
            if el.element_name in self.molecule_structure.duplicated_elements:
                continue
            # Skip barcode and UMI elements (already printed above)
            if el.element_name in barcode_umi_set:
                continue

            if el.element_type == ElementType.PolyT:
                if ElementType.PolyT.name not in self.detected_results:
                    res_str += "\t-1\t-1"
                    continue
                detected_element = self.detected_results[ElementType.PolyT.name]
                res_str += "\t%d\t%d" % (detected_element.start, detected_element.end)
            elif el.element_type.is_constant():
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
        """
        Get TSV header for result output in standard format.

        Format: #read_id, barcode, UMI, BC_score, valid_UMI, strand, [supplementary elements]
        """
        # Standard header fields
        header = "#read_id\tbarcode\tUMI\tBC_score\tvalid_UMI\tstrand"

        # Supplementary elements (non-barcode, non-UMI)
        barcode_umi_set = set(self._barcode_elements + self._umi_elements)
        for el in self.molecule_structure:
            if el.element_type == ElementType.cDNA:
                continue
            # Skip concatenated/duplicated parts (data lives under base name)
            if el.element_name in self.molecule_structure.elements_to_concatenate:
                continue
            if el.element_name in self.molecule_structure.duplicated_elements:
                continue
            # Skip barcode and UMI elements (already in standard header)
            if el.element_name in barcode_umi_set:
                continue

            if el.element_type == ElementType.PolyT:
                header += "\tpolyT_start\tpolyT_end"
            elif el.element_type.needs_only_coordinates():
                header += "\t%s_start\t%s_end" % (el.element_name, el.element_name)
            elif el.element_type.needs_sequence_extraction() or el.element_type.needs_correction():
                header += "\t%s_start\t%s_end\t%s_seq" % (el.element_name, el.element_name, el.element_name)
        return header

    # TODO: add simple serialization


class ReadStats:
    """Statistics tracker for ExtractionResult objects."""

    def __init__(self):
        self.read_count = 0
        self.bc_count = 0
        self.umi_count = 0
        self.pattern_counts = defaultdict(int)

    def add_read(self, barcode_detection_result: ExtractionResult) -> None:
        """Add an ExtractionResult to statistics."""
        self.read_count += 1
        # Count valid barcodes and UMIs
        if barcode_detection_result.has_barcode():
            self.bc_count += 1
        if barcode_detection_result.has_umi():
            self.umi_count += 1
        # Count individual detected elements
        for el in barcode_detection_result.detected_results:
            if barcode_detection_result.detected_results[el].start != -1:
                self.pattern_counts[el] += 1

    def __str__(self) -> str:
        human_readable_str = "Total reads:\t%d\n" % self.read_count
        human_readable_str += "Barcode detected:\t%d\n" % self.bc_count
        human_readable_str += "UMI detected:\t%d\n" % self.umi_count
        for a in self.pattern_counts:
            human_readable_str += "%s:\t%d\n" % (a, self.pattern_counts[a])
        return human_readable_str
