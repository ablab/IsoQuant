###########################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
from .molecule_structure import MoleculeStructure, ElementType

logger = logging.getLogger('IsoQuant')


class DetectedElement:
    def __init__(self, start=-1, end=-1, score=-1, seq=None):
        self.start = start
        self.end = end
        self.score = score
        self.seq = seq


class ExtractionResult:
    def __init__(self, molecule_structure: MoleculeStructure, read_id: str, strand: str, detected_results: dict):
        self.molecule_structure = molecule_structure
        self.read_id = read_id
        self.strand = strand
        self.detected_results = detected_results

    def is_valid(self) -> bool:
        """Check if a valid barcode was detected."""
        return True

    def update_coordinates(self, delta: int) -> None:
        """
        Shift all genomic coordinates by delta.

        Used when processing read subsequences.

        Args:
            delta: Amount to shift coordinates
        """
        pass

    def more_informative_than(self, that) -> bool:
        """
        Compare two results to determine which is more informative.

        Args:
            that: Another detection result

        Returns:
            True if this result is more informative

        Raises:
            NotImplementedError: Must be implemented by subclasses
        """
        raise NotImplementedError()

    def get_additional_attributes(self):
        """
        Get list of detected additional features (primer, linker, etc.).

        Returns:
            List of detected feature names

        Raises:
            NotImplementedError: Must be implemented by subclasses
        """
        return []

    def set_strand(self, strand: str) -> None:
        """Set the detected strand."""
        self.strand = strand

    def __str__(self) -> str:
        res_str = "%s\t%s" % (self.read_id, self.strand)
        for el in self.molecule_structure:
            if el.element_type == ElementType.cDNA: continue
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
                    res_str += "\t-1\t-1\t*"
                    continue
                detected_element = self.detected_results[el.element_name]
                res_str += "\t%d\t%d\t%s" % (detected_element.start, detected_element.end, detected_element.seq)
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
