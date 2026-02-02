############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Base classes for barcode detection results.

Provides result containers for barcode calling across different platforms.
"""

from collections import defaultdict
from typing import List, Optional, Dict, Iterable


def increase_if_valid(val: Optional[int], delta: int) -> Optional[int]:
    """
    Increment a coordinate value if it's valid.

    Args:
        val: Position value (-1 or None indicates invalid)
        delta: Amount to increment

    Returns:
        Incremented value if valid, otherwise original value
    """
    if val and val != -1:
        return val + delta
    return val





class BarcodeDetectionResult:
    """
    Base class for barcode detection results.

    Stores detected barcode, UMI, and quality scores for a single read.
    Implements the BarcodeResult protocol.
    """

    NOSEQ = "*"  # Sentinel for missing/undetected sequence

    def __init__(self, read_id: str, barcode: str = NOSEQ, UMI: str = NOSEQ,
                 BC_score: int = -1, UMI_good: bool = False, strand: str = "."):
        """
        Initialize barcode detection result.

        Args:
            read_id: Read identifier
            barcode: Detected barcode sequence (NOSEQ if not found)
            UMI: Detected UMI sequence (NOSEQ if not found)
            BC_score: Barcode alignment score
            UMI_good: Whether UMI passes quality filters
            strand: Detected strand ('+', '-', or '.')
        """
        self.read_id: str = read_id
        self._barcode: str = barcode if barcode else BarcodeDetectionResult.NOSEQ
        self._umi: str = UMI if UMI else BarcodeDetectionResult.NOSEQ
        self.BC_score: int = BC_score
        self.UMI_good: bool = UMI_good
        self.strand: str = strand

    # Primary getters (preferred interface)
    def get_barcode(self) -> str:
        """Get the detected barcode sequence."""
        return self._barcode

    def get_umi(self) -> str:
        """Get the detected UMI sequence."""
        return self._umi

    # Backward-compatible property accessors
    @property
    def barcode(self) -> str:
        """Barcode sequence. Prefer get_barcode() for new code."""
        return self._barcode

    @barcode.setter
    def barcode(self, value: str) -> None:
        self._barcode = value if value else BarcodeDetectionResult.NOSEQ

    @property
    def UMI(self) -> str:
        """UMI sequence. Prefer get_umi() for new code."""
        return self._umi

    @UMI.setter
    def UMI(self, value: str) -> None:
        self._umi = value if value else BarcodeDetectionResult.NOSEQ

    def is_valid(self) -> bool:
        """Check if a valid barcode was detected."""
        return self._barcode != BarcodeDetectionResult.NOSEQ

    def has_barcode(self) -> bool:
        """Check if a barcode was detected (alias for is_valid())."""
        return self._barcode != BarcodeDetectionResult.NOSEQ

    def has_umi(self) -> bool:
        """Check if a valid UMI was detected."""
        return self._umi != BarcodeDetectionResult.NOSEQ

    def update_coordinates(self, delta: int) -> None:
        """
        Shift all genomic coordinates by delta.

        Used when processing read subsequences.

        Args:
            delta: Amount to shift coordinates
        """
        pass

    def more_informative_than(self, that: 'BarcodeDetectionResult') -> bool:
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

    def get_additional_attributes(self) -> List[str]:
        """
        Get list of detected additional features (primer, linker, etc.).

        Returns:
            List of detected feature names. Empty list for base class.
        """
        return []

    def set_strand(self, strand: str) -> None:
        """Set the detected strand."""
        self.strand = strand

    def __str__(self) -> str:
        """Format result as TSV line."""
        return "%s\t%s\t%s\t%d\t%s\t%s" % (self.read_id, self._barcode, self._umi,
                                           self.BC_score, self.UMI_good, self.strand)

    @staticmethod
    def header() -> str:
        """Static header for class-level access."""
        return "#read_id\tbarcode\tUMI\tBC_score\tvalid_UMI\tstrand"


class LinkerBarcodeDetectionResult(BarcodeDetectionResult):
    """
    Detection result for platforms with double barcodes (e.g., Curio, Stereo-seq).

    Extends base result with positions of additional features:
    polyT tail, primer, and linker sequences.
    """

    def __init__(self, read_id: str, barcode: str = BarcodeDetectionResult.NOSEQ,
                 UMI: str = BarcodeDetectionResult.NOSEQ,
                 BC_score: int = -1, UMI_good: bool = False, strand: str = ".",
                 polyT: int = -1, primer: int = -1, linker_start: int = -1, linker_end: int = -1):
        """
        Initialize double barcode detection result.

        Args:
            read_id: Read identifier
            barcode: Detected barcode (concatenated if split by linker)
            UMI: Detected UMI sequence
            BC_score: Barcode alignment score
            UMI_good: Whether UMI passes quality filters
            strand: Detected strand
            polyT: Position of polyT tail start (-1 if not found)
            primer: Position of primer end (-1 if not found)
            linker_start: Position of linker start (-1 if not found)
            linker_end: Position of linker end (-1 if not found)
        """
        BarcodeDetectionResult.__init__(self, read_id, barcode, UMI, BC_score, UMI_good, strand)
        self.primer: int = primer
        self.linker_start: int = linker_start
        self.linker_end: int = linker_end
        self.polyT: int = polyT

    def update_coordinates(self, delta: int) -> None:
        self.primer = increase_if_valid(self.primer, delta)
        self.linker_start = increase_if_valid(self.linker_start, delta)
        self.linker_end = increase_if_valid(self.linker_end, delta)
        self.polyT = increase_if_valid(self.polyT, delta)

    def more_informative_than(self, that: 'LinkerBarcodeDetectionResult') -> bool:
        if self.BC_score != that.BC_score:
            return self.BC_score > that.BC_score
        if self.linker_start != that.linker_start:
            return self.linker_start > that.linker_start
        if self.primer != that.primer:
            return self.primer > that.primer
        return self.polyT > that.polyT

    def get_additional_attributes(self) -> List[str]:
        attr = []
        if self.polyT != -1:
            attr.append("PolyT detected")
        if self.primer != -1:
            attr.append("Primer detected")
        if self.linker_start != -1:
            attr.append("Linker detected")
        return attr

    def __str__(self) -> str:
        return (BarcodeDetectionResult.__str__(self) +
                "\t%d\t%d\t%d\t%d" % (self.polyT, self.primer, self.linker_start, self.linker_end))

    @staticmethod
    def header() -> str:
        """Static header for class-level access."""
        return BarcodeDetectionResult.header() + "\tpolyT_start\tprimer_end\tlinker_start\tlinker_end"


class TSOBarcodeDetectionResult(LinkerBarcodeDetectionResult):
    """Detection result for Stereo-seq with TSO detection."""

    def __init__(self, read_id: str, barcode: str = BarcodeDetectionResult.NOSEQ,
                 UMI: str = BarcodeDetectionResult.NOSEQ,
                 BC_score: int = -1, UMI_good: bool = False, strand: str = ".",
                 polyT: int = -1, primer: int = -1, linker_start: int = -1,
                 linker_end: int = -1, tso: int = -1):
        LinkerBarcodeDetectionResult.__init__(self, read_id, barcode, UMI, BC_score, UMI_good, strand,
                                              polyT, primer, linker_start, linker_end)
        self.tso5: int = tso

    def update_coordinates(self, delta: int) -> None:
        self.tso5 = increase_if_valid(self.tso5, delta)
        LinkerBarcodeDetectionResult.update_coordinates(self, delta)

    def __str__(self) -> str:
        return (LinkerBarcodeDetectionResult.__str__(self) +
                "\t%d" % self.tso5)

    def get_additional_attributes(self) -> List[str]:
        attr = []
        if self.polyT != -1:
            attr.append("PolyT detected")
        if self.primer != -1:
            attr.append("Primer detected")
        if self.linker_start != -1:
            attr.append("Linker detected")
        if self.tso5 != -1:
            attr.append("TSO detected")

    @staticmethod
    def header() -> str:
        """Static header for class-level access."""
        return LinkerBarcodeDetectionResult.header() + "\tTSO5"


class TenXBarcodeDetectionResult(BarcodeDetectionResult):
    """Detection result for 10x Genomics platforms."""

    def __init__(self, read_id: str, barcode: str = BarcodeDetectionResult.NOSEQ,
                 UMI: str = BarcodeDetectionResult.NOSEQ,
                 BC_score: int = -1, UMI_good: bool = False, strand: str = ".",
                 polyT: int = -1, r1: int = -1):
        BarcodeDetectionResult.__init__(self, read_id, barcode, UMI, BC_score, UMI_good, strand)
        self.r1: int = r1
        self.polyT: int = polyT

    def update_coordinates(self, delta: int) -> None:
        self.r1 = increase_if_valid(self.r1, delta)
        self.polyT = increase_if_valid(self.polyT, delta)

    def more_informative_than(self, that: 'TenXBarcodeDetectionResult') -> bool:
        if self.polyT != that.polyT:
            return self.polyT > that.polyT
        if self.r1 != that.r1:
            return self.r1 > that.r1
        return self.BC_score > that.BC_score

    def get_additional_attributes(self) -> List[str]:
        attr = []
        if self.polyT != -1:
            attr.append("PolyT detected")
        if self.r1 != -1:
            attr.append("R1 detected")
        return attr

    def __str__(self) -> str:
        return (BarcodeDetectionResult.__str__(self) +
                "\t%d\t%d" % (self.polyT, self.r1))

    @staticmethod
    def header() -> str:
        """Static header for class-level access."""
        return BarcodeDetectionResult.header() + "\tpolyT_start\tR1_end"


class SplittingBarcodeDetectionResult:
    """Result container for read splitting modes (multiple barcodes per read)."""

    NOSEQ = BarcodeDetectionResult.NOSEQ  # For consistency with protocol

    def __init__(self, read_id: str):
        self.read_id: str = read_id
        self.strand: str = "."  # For protocol compatibility
        self.detected_patterns: List[TSOBarcodeDetectionResult] = []

    def append(self, barcode_detection_result: TSOBarcodeDetectionResult) -> None:
        self.detected_patterns.append(barcode_detection_result)

    def empty(self) -> bool:
        return not self.detected_patterns

    def filter(self) -> None:
        if not self.detected_patterns:
            return
        barcoded_results = []
        for r in self.detected_patterns:
            if r.is_valid():
                barcoded_results.append(r)

        if not barcoded_results:
            self.detected_patterns = [self.detected_patterns[0]]
        else:
            self.detected_patterns = barcoded_results

    def get_barcode(self) -> str:
        """Get barcode from first valid pattern, or NOSEQ if none."""
        for r in self.detected_patterns:
            if r.is_valid():
                return r.get_barcode()
        return self.NOSEQ

    def get_umi(self) -> str:
        """Get UMI from first valid pattern, or NOSEQ if none."""
        for r in self.detected_patterns:
            if r.is_valid():
                return r.get_umi()
        return self.NOSEQ

    def is_valid(self) -> bool:
        """Check if any pattern has a valid barcode."""
        return any(r.is_valid() for r in self.detected_patterns)

    def has_barcode(self) -> bool:
        """Check if any pattern has a valid barcode."""
        return any(r.has_barcode() for r in self.detected_patterns)

    def has_umi(self) -> bool:
        """Check if any pattern has a valid UMI."""
        return any(r.has_umi() for r in self.detected_patterns)

    def set_strand(self, strand: str) -> None:
        """Set strand for all patterns."""
        self.strand = strand
        for r in self.detected_patterns:
            r.set_strand(strand)

    def update_coordinates(self, delta: int) -> None:
        """Update coordinates for all patterns."""
        for r in self.detected_patterns:
            r.update_coordinates(delta)

    def more_informative_than(self, other: 'SplittingBarcodeDetectionResult') -> bool:
        """Compare by number of valid patterns."""
        self_valid = sum(1 for r in self.detected_patterns if r.is_valid())
        other_valid = sum(1 for r in other.detected_patterns if r.is_valid())
        return self_valid > other_valid

    def get_additional_attributes(self) -> List[str]:
        """Get combined attributes from all patterns."""
        attrs = set()
        for r in self.detected_patterns:
            attrs.update(r.get_additional_attributes())
        return list(attrs)

    def __str__(self) -> str:
        """Format all patterns as TSV lines."""
        return "\n".join(str(r) for r in self.detected_patterns)

    @staticmethod
    def header() -> str:
        """Get TSV header for result output."""
        return TSOBarcodeDetectionResult.header()


class ReadStats:
    """
    Statistics tracker for barcode detection results.

    Accumulates counts of processed reads, detected barcodes, valid UMIs,
    and platform-specific features (primers, linkers, polyT tails, etc.).
    """

    def __init__(self):
        """Initialize empty statistics."""
        self.read_count: int = 0
        self.bc_count: int = 0
        self.umi_count: int = 0
        self.additional_attributes_counts: Dict[str, int] = defaultdict(int)

    def add_read(self, barcode_detection_result) -> None:
        """
        Add a read result to statistics.

        Args:
            barcode_detection_result: Detection result to accumulate (implements BarcodeResult protocol)
        """
        self.read_count += 1
        # Count detected features (primer, linker, etc.)
        for a in barcode_detection_result.get_additional_attributes():
            self.additional_attributes_counts[a] += 1
        # Count valid barcode
        if barcode_detection_result.has_barcode():
            self.bc_count += 1
        # Count valid UMI
        if barcode_detection_result.has_umi():
            self.umi_count += 1

    def add_custom_stats(self, stat_name: str, val: int) -> None:
        """
        Add custom statistic value.

        Args:
            stat_name: Name of statistic
            val: Count to add
        """
        self.additional_attributes_counts[stat_name] += val

    def __str__(self) -> str:
        """Format statistics as human-readable string."""
        human_readable_str = ("Total reads\t%d\nBarcode detected\t%d\nReliable UMI\t%d\n" %
                              (self.read_count, self.bc_count, self.umi_count))
        for a in self.additional_attributes_counts:
            human_readable_str += "%s\t%d\n" % (a, self.additional_attributes_counts[a])
        return human_readable_str

    def __iter__(self) -> Iterable[str]:
        """Iterate over statistics as formatted strings."""
        yield "Total reads: %d" % self.read_count
        yield "Barcode detected: %d" % self.bc_count
        yield "Reliable UMI: %d" % self.umi_count
        for a in self.additional_attributes_counts:
            yield "%s: %d" % (a, self.additional_attributes_counts[a])
