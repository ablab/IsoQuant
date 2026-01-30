############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Tests for molecule structure parsing and element types.

Tests the universal barcode calling infrastructure:
- ElementType enum and its method classification
- MoleculeElement initialization for different element types
- MoleculeStructure parsing from MDF format
"""

import os
import pytest
from io import StringIO

from src.barcode_calling.callers.molecule_structure import (
    ElementType,
    MoleculeElement,
    MoleculeStructure
)

# Test data directory
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "universal_data")


class TestElementType:
    """Test ElementType enum methods."""

    def test_needs_correction_var_file(self):
        """VAR_FILE elements need barcode correction."""
        assert ElementType.VAR_FILE.needs_correction() is True

    def test_needs_correction_var_list(self):
        """VAR_LIST elements need barcode correction."""
        assert ElementType.VAR_LIST.needs_correction() is True

    def test_needs_correction_var_any(self):
        """VAR_ANY elements do NOT need correction."""
        assert ElementType.VAR_ANY.needs_correction() is False

    def test_needs_correction_const(self):
        """CONST elements do NOT need correction."""
        assert ElementType.CONST.needs_correction() is False

    def test_needs_sequence_extraction_var_any(self):
        """VAR_ANY elements need sequence extraction."""
        assert ElementType.VAR_ANY.needs_sequence_extraction() is True

    def test_needs_sequence_extraction_var_file(self):
        """VAR_FILE elements do NOT need sequence extraction (they get corrected)."""
        assert ElementType.VAR_FILE.needs_sequence_extraction() is False

    def test_needs_sequence_extraction_const(self):
        """CONST elements do NOT need sequence extraction."""
        assert ElementType.CONST.needs_sequence_extraction() is False

    def test_is_constant_const(self):
        """CONST is a constant element."""
        assert ElementType.CONST.is_constant() is True

    def test_is_constant_var_any(self):
        """VAR_ANY is NOT a constant element."""
        assert ElementType.VAR_ANY.is_constant() is False

    def test_is_constant_polyt(self):
        """PolyT is NOT a constant element."""
        assert ElementType.PolyT.is_constant() is False

    def test_is_variable_var_any(self):
        """VAR_ANY is a variable element."""
        assert ElementType.VAR_ANY.is_variable() is True

    def test_is_variable_var_list(self):
        """VAR_LIST is a variable element."""
        assert ElementType.VAR_LIST.is_variable() is True

    def test_is_variable_var_file(self):
        """VAR_FILE is a variable element."""
        assert ElementType.VAR_FILE.is_variable() is True

    def test_is_variable_const(self):
        """CONST is NOT a variable element."""
        assert ElementType.CONST.is_variable() is False

    def test_is_variable_polyt(self):
        """PolyT is NOT a variable element."""
        assert ElementType.PolyT.is_variable() is False

    def test_is_base_separator_var_any_separator(self):
        """VAR_ANY_SEPARATOR is a separator."""
        assert ElementType.VAR_ANY_SEPARATOR.is_base_separator() is True

    def test_is_base_separator_var_any_non_t_separator(self):
        """VAR_ANY_NON_T_SEPARATOR is a separator."""
        assert ElementType.VAR_ANY_NON_T_SEPARATOR.is_base_separator() is True

    def test_is_base_separator_var_any(self):
        """VAR_ANY is NOT a separator."""
        assert ElementType.VAR_ANY.is_base_separator() is False


class TestMoleculeElement:
    """Test MoleculeElement initialization."""

    def test_init_polyt(self):
        """Test PolyT element initialization."""
        element = MoleculeElement("PolyT", ElementType.PolyT)

        assert element.element_name == "PolyT"
        assert element.element_type == ElementType.PolyT
        assert element.element_value is None
        assert element.element_length == -1

    def test_init_cdna(self):
        """Test cDNA element initialization."""
        element = MoleculeElement("cDNA", ElementType.cDNA)

        assert element.element_name == "cDNA"
        assert element.element_type == ElementType.cDNA
        assert element.element_value is None
        assert element.element_length == -1

    def test_init_const(self):
        """Test CONST element initialization."""
        sequence = "CTACACGACGCTCTTCCGATCT"
        element = MoleculeElement("R1", ElementType.CONST, sequence)

        assert element.element_name == "R1"
        assert element.element_type == ElementType.CONST
        assert element.element_value == sequence
        assert element.element_length == len(sequence)

    def test_init_var_list(self):
        """Test VAR_LIST element initialization."""
        barcodes = "AAAA,CCCC,GGGG,TTTT"
        element = MoleculeElement("Barcode", ElementType.VAR_LIST, barcodes)

        assert element.element_name == "Barcode"
        assert element.element_type == ElementType.VAR_LIST
        assert element.element_value == ["AAAA", "CCCC", "GGGG", "TTTT"]
        assert element.element_length == 4  # Length of first barcode

    def test_init_var_any(self):
        """Test VAR_ANY element initialization."""
        element = MoleculeElement("UMI", ElementType.VAR_ANY, "12")

        assert element.element_name == "UMI"
        assert element.element_type == ElementType.VAR_ANY
        assert element.element_value is None
        assert element.element_length == 12

    def test_init_var_file(self):
        """Test VAR_FILE element initialization with real file."""
        barcode_file = os.path.join(TEST_DATA_DIR, "barcodes_small.tsv")
        element = MoleculeElement("Barcode", ElementType.VAR_FILE, barcode_file)

        assert element.element_name == "Barcode"
        assert element.element_type == ElementType.VAR_FILE
        assert len(element.element_value) == 8  # 8 barcodes in test file
        assert element.element_length == 16  # 16bp barcodes

    def test_init_var_any_separator(self):
        """Test VAR_ANY_SEPARATOR element initialization."""
        element = MoleculeElement("Sep", ElementType.VAR_ANY_SEPARATOR, "5")

        assert element.element_name == "Sep"
        assert element.element_type == ElementType.VAR_ANY_SEPARATOR
        assert element.element_value is None
        assert element.element_length == 5


class TestMoleculeStructure:
    """Test MoleculeStructure parsing."""

    def test_parse_simple_structure(self):
        """Test parsing simple molecule structure."""
        mdf_content = """Barcode:UMI:PolyT:cDNA
Barcode\tVAR_LIST\tAAAA,CCCC,GGGG
UMI\tVAR_ANY\t10
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        elements = list(structure)
        assert len(elements) == 4

        assert elements[0].element_name == "Barcode"
        assert elements[0].element_type == ElementType.VAR_LIST

        assert elements[1].element_name == "UMI"
        assert elements[1].element_type == ElementType.VAR_ANY

        assert elements[2].element_name == "PolyT"
        assert elements[2].element_type == ElementType.PolyT

        assert elements[3].element_name == "cDNA"
        assert elements[3].element_type == ElementType.cDNA

    def test_parse_with_const(self):
        """Test parsing structure with constant element."""
        mdf_content = """R1:Barcode:UMI:PolyT:cDNA
R1\tCONST\tCTACACGACGCTCTTCCGATCT
Barcode\tVAR_LIST\tAAAA,CCCC,GGGG
UMI\tVAR_ANY\t12
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        elements = list(structure)
        assert len(elements) == 5

        assert elements[0].element_name == "R1"
        assert elements[0].element_type == ElementType.CONST
        assert elements[0].element_value == "CTACACGACGCTCTTCCGATCT"

    def test_parse_from_file_simple(self):
        """Test parsing simple.mdf file."""
        mdf_path = os.path.join(TEST_DATA_DIR, "simple.mdf")

        with open(mdf_path) as f:
            structure = MoleculeStructure(f)

        elements = list(structure)

        # Simple.mdf: Barcode:UMI:PolyT:cDNA
        assert len(elements) == 4
        assert elements[0].element_name == "Barcode"
        assert elements[0].element_type == ElementType.VAR_LIST
        assert elements[1].element_name == "UMI"
        assert elements[1].element_type == ElementType.VAR_ANY
        assert elements[1].element_length == 10

    def test_parse_from_file_10x(self):
        """Test parsing 10x-like structure.

        Uses inline definition to avoid relative path issues with VAR_FILE.
        """
        # Use inline definition to avoid file path issues
        mdf_content = """R1:Barcode:UMI:PolyT:cDNA:TSO
R1\tCONST\tCTACACGACGCTCTTCCGATCT
Barcode\tVAR_LIST\tATCCTTAGTGTTTGTC,CTTAGGAGTGGTTAGC,TGAGGGCCAACGTGCT
UMI\tVAR_ANY\t12
TSO\tCONST\tCCCATGTACTCTGCGTTGATACCACTGCTT
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        elements = list(structure)

        # R1:Barcode:UMI:PolyT:cDNA:TSO
        assert len(elements) == 6
        assert elements[0].element_name == "R1"
        assert elements[0].element_type == ElementType.CONST
        assert elements[1].element_name == "Barcode"
        assert elements[1].element_type == ElementType.VAR_LIST

    def test_iteration(self):
        """Test iterating over structure elements."""
        mdf_content = """Barcode:UMI:cDNA
Barcode\tVAR_LIST\tAAAA,CCCC
UMI\tVAR_ANY\t8
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        names = [el.element_name for el in structure]
        assert names == ["Barcode", "UMI", "cDNA"]

    def test_from_element_list(self):
        """Test creating structure from element list."""
        elements = [
            MoleculeElement("Barcode", ElementType.VAR_LIST, "AAAA,CCCC"),
            MoleculeElement("UMI", ElementType.VAR_ANY, "10"),
            MoleculeElement("cDNA", ElementType.cDNA)
        ]

        structure = MoleculeStructure.from_element_list(elements)

        assert len(list(structure)) == 3
        assert list(structure)[0].element_name == "Barcode"

    def test_concat_delimiter(self):
        """Test concatenation delimiter constant."""
        assert MoleculeStructure.CONCAT_DELIM == "|"

    def test_dupl_delimiter(self):
        """Test duplication delimiter constant."""
        assert MoleculeStructure.DUPL_DELIM == "/"


class TestMoleculeStructureBarcodeUMI:
    """Test barcode/UMI element identification."""

    def test_barcode_elements_single(self):
        """Test single barcode element identification."""
        mdf_content = """Barcode:UMI:cDNA
Barcode\tVAR_LIST\tAAAA,CCCC
UMI\tVAR_ANY\t10
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert structure.barcode_elements == ["Barcode"]

    def test_umi_elements_single(self):
        """Test single UMI element identification."""
        mdf_content = """Barcode:UMI:cDNA
Barcode\tVAR_LIST\tAAAA,CCCC
UMI\tVAR_ANY\t10
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert structure.umi_elements == ["UMI"]

    def test_barcode_elements_multiple(self):
        """Test multiple barcode elements identification."""
        mdf_content = """barcode1:barcode2:UMI:cDNA
barcode1\tVAR_LIST\tAAAA,CCCC
barcode2\tVAR_LIST\tGGGG,TTTT
UMI\tVAR_ANY\t10
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert structure.barcode_elements == ["barcode1", "barcode2"]

    def test_umi_elements_multiple(self):
        """Test multiple UMI elements identification."""
        mdf_content = """Barcode:UMI1:UMI2:cDNA
Barcode\tVAR_LIST\tAAAA,CCCC
UMI1\tVAR_ANY\t6
UMI2\tVAR_ANY\t6
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert structure.umi_elements == ["UMI1", "UMI2"]

    def test_barcode_case_insensitive(self):
        """Test that barcode prefix is case-insensitive."""
        mdf_content = """BARCODE:Barcode2:barcode3:cDNA
BARCODE\tVAR_LIST\tAAAA,CCCC
Barcode2\tVAR_LIST\tGGGG,TTTT
barcode3\tVAR_LIST\tACGT,TGCA
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert len(structure.barcode_elements) == 3
        assert "BARCODE" in structure.barcode_elements
        assert "Barcode2" in structure.barcode_elements
        assert "barcode3" in structure.barcode_elements

    def test_umi_case_insensitive(self):
        """Test that UMI prefix is case-insensitive."""
        mdf_content = """Barcode:UMI:umi2:cDNA
Barcode\tVAR_LIST\tAAAA,CCCC
UMI\tVAR_ANY\t6
umi2\tVAR_ANY\t6
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert len(structure.umi_elements) == 2
        assert "UMI" in structure.umi_elements
        assert "umi2" in structure.umi_elements

    def test_no_barcode_elements(self):
        """Test structure with no barcode elements."""
        mdf_content = """R1:UMI:cDNA
R1\tCONST\tACGT
UMI\tVAR_ANY\t10
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert structure.barcode_elements == []

    def test_no_umi_elements(self):
        """Test structure with no UMI elements."""
        mdf_content = """Barcode:cDNA
Barcode\tVAR_LIST\tAAAA,CCCC
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert structure.umi_elements == []

    def test_from_element_list_identifies_barcodes(self):
        """Test that from_element_list also identifies barcode/UMI elements."""
        elements = [
            MoleculeElement("Barcode", ElementType.VAR_LIST, "AAAA,CCCC"),
            MoleculeElement("UMI", ElementType.VAR_ANY, "10"),
            MoleculeElement("cDNA", ElementType.cDNA)
        ]

        structure = MoleculeStructure.from_element_list(elements)

        assert structure.barcode_elements == ["Barcode"]
        assert structure.umi_elements == ["UMI"]


class TestMoleculeStructureEdgeCases:
    """Test edge cases and error handling."""

    def test_element_order_preserved(self):
        """Test that element order from header line is preserved."""
        mdf_content = """cDNA:UMI:Barcode:PolyT
Barcode\tVAR_LIST\tAAAA,CCCC
UMI\tVAR_ANY\t8
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        names = [el.element_name for el in structure]
        # Order should match header line, not definition order
        assert names == ["cDNA", "UMI", "Barcode", "PolyT"]

    def test_whitespace_handling(self):
        """Test that whitespace in MDF is handled correctly."""
        mdf_content = """  Barcode : UMI : cDNA
Barcode    VAR_LIST    AAAA,CCCC
UMI        VAR_ANY     8
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        elements = list(structure)
        assert len(elements) == 3

    def test_multiple_var_list_barcodes(self):
        """Test VAR_LIST with many barcodes."""
        barcodes = ",".join(["ACGT" * 4] * 100)  # 100 barcodes
        mdf_content = f"""Barcode:cDNA
Barcode\tVAR_LIST\t{barcodes}
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        elements = list(structure)
        assert len(elements[0].element_value) == 100


class TestConcatenatedElements:
    """Test concatenated element parsing and identification."""

    def test_parse_concatenated_elements(self):
        """Test parsing MDF with concatenated barcode parts."""
        mdf_content = """Barcode|1:Barcode|2:UMI:PolyT:cDNA
Barcode|1\tVAR_LIST\tAAAACCCC,GGGGTTTT\t4
Barcode|2\tVAR_LIST\tAAAACCCC,GGGGTTTT\t4
UMI\tVAR_ANY\t8
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert len(structure.ordered_elements) == 5
        assert "Barcode|1" in structure.elements_to_concatenate
        assert "Barcode|2" in structure.elements_to_concatenate
        assert structure.elements_to_concatenate["Barcode|1"] == ("Barcode", 1)
        assert structure.elements_to_concatenate["Barcode|2"] == ("Barcode", 2)

    def test_concatenated_elements_counts(self):
        """Test that concatenated_elements_counts has correct base name and count."""
        mdf_content = """Barcode|1:Barcode|2:cDNA
Barcode|1\tVAR_LIST\tAAAACCCC,GGGGTTTT\t4
Barcode|2\tVAR_LIST\tAAAACCCC,GGGGTTTT\t4
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert "Barcode" in structure.concatenated_elements_counts
        assert structure.concatenated_elements_counts["Barcode"] == 2

    def test_concatenated_three_parts(self):
        """Test three-part concatenated element."""
        mdf_content = """Barcode|1:Barcode|2:Barcode|3:cDNA
Barcode|1\tVAR_LIST\tAAAACCCCGGGG,TTTTAAAACCCC\t4
Barcode|2\tVAR_LIST\tAAAACCCCGGGG,TTTTAAAACCCC\t4
Barcode|3\tVAR_LIST\tAAAACCCCGGGG,TTTTAAAACCCC\t4
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert structure.concatenated_elements_counts["Barcode"] == 3

    def test_concatenated_barcode_umi_identification(self):
        """Test that concatenated barcode elements are identified by base name."""
        mdf_content = """Barcode|1:Barcode|2:UMI:cDNA
Barcode|1\tVAR_LIST\tAAAACCCC,GGGGTTTT\t4
Barcode|2\tVAR_LIST\tAAAACCCC,GGGGTTTT\t4
UMI\tVAR_ANY\t8
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        # Base name "Barcode" should be in barcode_elements, not "Barcode|1", "Barcode|2"
        assert structure.barcode_elements == ["Barcode"]
        assert "Barcode|1" not in structure.barcode_elements
        assert "Barcode|2" not in structure.barcode_elements
        assert structure.umi_elements == ["UMI"]

    def test_concatenated_umi_identification(self):
        """Test that concatenated UMI elements are identified by base name."""
        mdf_content = """Barcode:UMI|1:UMI|2:cDNA
Barcode\tVAR_LIST\tAAAA,CCCC
UMI|1\tVAR_ANY\t4
UMI|2\tVAR_ANY\t4
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        assert structure.barcode_elements == ["Barcode"]
        assert structure.umi_elements == ["UMI"]

    def test_concatenated_explicit_length(self):
        """Test that explicit length (4th field) is correctly parsed for concatenated parts."""
        mdf_content = """Barcode|1:Barcode|2:cDNA
Barcode|1\tVAR_LIST\tAAAACCCC,GGGGTTTT\t4
Barcode|2\tVAR_LIST\tAAAACCCC,GGGGTTTT\t4
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))

        # Each part should have length 4, not the full barcode length 8
        for el in structure:
            if el.element_name == "Barcode|1":
                assert el.element_length == 4
            elif el.element_name == "Barcode|2":
                assert el.element_length == 4

    def test_concatenated_elements_must_be_variable(self):
        """Test that non-variable concatenated elements are rejected."""
        mdf_content = """Adapter|1:Adapter|2:cDNA
Adapter|1\tCONST\tACGT
Adapter|2\tCONST\tTGCA
"""
        with pytest.raises(SystemExit):
            MoleculeStructure(iter(mdf_content.strip().split('\n')))

    def test_concatenated_nonconsecutive_indices_rejected(self):
        """Test that non-consecutive indices (e.g., 1 and 3) are rejected."""
        mdf_content = """Barcode|1:Barcode|3:cDNA
Barcode|1\tVAR_LIST\tAAAA,CCCC\t4
Barcode|3\tVAR_LIST\tAAAA,CCCC\t4
"""
        with pytest.raises(SystemExit):
            MoleculeStructure(iter(mdf_content.strip().split('\n')))


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
