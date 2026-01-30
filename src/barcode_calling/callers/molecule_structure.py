###########################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import sys
from enum import unique, Enum
from typing import Iterator, List
from collections import defaultdict

from ...error_codes import IsoQuantExitCode
from ..common import load_barcodes

logger = logging.getLogger('IsoQuant')


@unique
class ElementType(Enum):
    PolyT = 1
    cDNA = 2
    CONST = 10
    VAR_ANY = 20
    VAR_LIST = 21
    VAR_FILE = 22
    VAR_ANY_SEPARATOR = 31
    VAR_ANY_NON_T_SEPARATOR = 32

    def needs_correction(self):
        return self in (ElementType.VAR_FILE, ElementType.VAR_LIST)

    def needs_sequence_extraction(self):
        return self == ElementType.VAR_ANY

    def is_base_separator(self):
        return self in (ElementType.VAR_ANY_SEPARATOR, ElementType.VAR_ANY_NON_T_SEPARATOR)

    def is_variable(self):
        return self in (ElementType.VAR_ANY, ElementType.VAR_LIST, ElementType.VAR_FILE)

    def is_constant(self):
        return self == ElementType.CONST

    def needs_only_coordinates(self):
        """Return True if element only needs start/end coordinates (no sequence extraction)."""
        return self.is_constant()


class MoleculeElement:
    def __init__(self, element_name: str, element_type: ElementType, element_value_1 = None, element_value_2 = None):
        self.element_name = element_name
        self.element_type = element_type

        if self.element_type in [ElementType.PolyT, ElementType.cDNA]:
            self.element_value = None
            self.element_length = -1
        elif self.element_type == ElementType.CONST:
            self.element_value = element_value_1
            self.element_length = len(self.element_value)
        elif self.element_type == ElementType.VAR_FILE:
            self.element_value = load_barcodes(element_value_1, needs_iterator=False)
            if element_value_2 is None:
                self.element_length = len(self.element_value[0])
            else:
                try:
                    self.element_length = int(element_value_1)
                except ValueError:
                    logger.critical("Incorrectly specified length for element %s: %s" % (self.element_name, element_value_1))
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
        elif self.element_type == ElementType.VAR_LIST:
            self.element_value = element_value_1.split(',')
            if element_value_2 is None:
                self.element_length = len(self.element_value[0])
            else:
                try:
                    self.element_length = int(element_value_1)
                except ValueError:
                    logger.critical(
                        "Incorrectly specified length for element %s: %s" % (self.element_name, element_value_1))
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
        elif self.element_type in (ElementType.VAR_ANY, ElementType.VAR_ANY_SEPARATOR, ElementType.VAR_ANY_NON_T_SEPARATOR):
            self.element_value = None
            try:
                self.element_length = int(element_value_1)
            except ValueError:
                logger.critical("Incorrectly specified length for element %s: %s" % (self.element_name, element_value_1))
                sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
        else:
            logger.critical("Wrong element type %s" % element_type)
            sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)


class MoleculeStructure:
    """
    Defines the structure of a molecule for barcode/UMI extraction.

    Parses molecule definition files (MDF) and identifies which elements
    are barcodes and UMIs based on naming convention:
    - Elements with prefix "barcode" (case-insensitive) are barcodes
    - Elements with prefix "umi" (case-insensitive) are UMIs
    """

    CONCAT_DELIM = "|"
    DUPL_DELIM = "/"

    def __init__(self, str_iterator):
        self.ordered_elements: List[MoleculeElement] = []
        self.elements_to_concatenate = {}
        self.concatenated_elements_counts = {}
        self.duplicated_elements = {}
        self.duplicated_elements_counts = {}
        # Barcode/UMI element names identified by naming convention
        self._barcode_elements: List[str] = []
        self._umi_elements: List[str] = []

        l = next(str_iterator)
        elements = list(map(lambda x: x.strip(), l.strip().split(':')))
        element_properties = {}
        for l in str_iterator:
            v = l.strip().split()
            element_name = v[0]
            if element_name in element_properties:
                logger.critical("Duplicated element name %s" % element_name)
                sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)

            if len(v) == 3:
                element_properties[element_name] = (v[1], v[2], None)
            elif len(v) == 4:
                element_properties[element_name] = (v[1], v[2], v[3])
            else:
                logger.critical("Incorrect number of properties in line %s" % l.strip())
                sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)

        element_name: str
        for element_name in elements:
            if element_name not in element_properties:
                if element_name not in ElementType.__dict__:
                    logger.critical("Molecule element %s was not described in the format file" % element_name)
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
                element_type = ElementType[element_name]
                self.ordered_elements.append(MoleculeElement(element_name, element_type))
            else:
                element_type, element_val1, element_val2 = element_properties[element_name]
                if element_type not in ElementType.__dict__:
                    logger.critical("Molecule element type %s is not among the possible types" % element_type)
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
                element_type = ElementType[element_type]
                self.ordered_elements.append(MoleculeElement(element_name, element_type, element_val1, element_val2))

            if self.CONCAT_DELIM in element_name:
                if not (element_type.needs_sequence_extraction() or element_type.needs_correction()):
                    logger.critical("Concatenated elements must be variable (fix element %s)" % element_name)
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
                v = element_name.split(self.CONCAT_DELIM)
                if len(v) != 2:
                    logger.critical("Incorrect concatenated element %s" % element_name)
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
                try:
                    self.elements_to_concatenate[element_name] = (v[0], int(v[1]))
                except ValueError:
                    logger.critical("Incorrectly specified index for concatenated element %s: %s" % (element_name, v[1]))
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)

            elif self.DUPL_DELIM in element_name:
                if not (element_type.needs_sequence_extraction() or element_type.needs_correction()):
                    logger.critical("Concatenated elements must be variable (fix element %s)" % element_name)
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
                v = element_name.split(self.DUPL_DELIM)
                if len(v) != 2:
                    logger.critical("Incorrect duplicated element %s" % element_name)
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
                try:
                    self.duplicated_elements[element_name] = (v[0], int(v[1]))
                except ValueError:
                    logger.critical("Incorrectly specified index for duplicated element %s: %s" % (element_name, v[1]))
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)

        self._check_linked_elements(self.elements_to_concatenate)
        self.concatenated_elements_counts = self._check_indices(self.elements_to_concatenate)
        self._check_linked_elements(self.duplicated_elements)
        self.duplicated_elements_counts = self._check_indices(self.duplicated_elements)
        self._identify_barcode_umi_elements()

    def _check_linked_elements(self, element_dict):
        # check that duplicated or concatenated elements have the same values (e.g. barcode whitelist)
        # base element name -> List[MoleculeElement]
        linked_elements = defaultdict(list)
        for el in self.ordered_elements:
            if el.element_name in element_dict:
                base_element_name = element_dict[el][0]
                linked_elements[base_element_name].append(el)

        for element_list in linked_elements.values():
            if len(set(el.element_value for el in element_list)) != 1:
                logger.critical("Linked elements must have the same type (fix elements: %s)" % ", ".join(el.element_name for el in element_list))
                sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
            if len(set(el.element_value for el in element_list)) != 1:
                logger.critical("Linked elements must have the same values (fix elements: %s)" % ", ".join(el.element_name for el in element_list))
                sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)

    def _check_indices(self, element_dict):
        # check for consecutive indices
        index_dict = defaultdict(list)
        for el in element_dict:
            base_element_name = element_dict[el][0]
            index_dict[base_element_name].append(element_dict[el][1])

        for el in index_dict:
            if sorted(index_dict[el]) != list(range(1, len(index_dict[el]) + 1)):
                logger.critical("Concatenated elements must be consecutive from 1 to %d (fix element %s, indices: %s)" % (len(index_dict[el]), el, str(index_dict[el])))
                sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)

        return {el : len(val) for el, val in element_dict.items()}

    def _identify_barcode_umi_elements(self) -> None:
        """
        Identify barcode/UMI elements by naming convention.

        Elements with prefix "barcode" (case-insensitive) are barcodes.
        Elements with prefix "umi" (case-insensitive) are UMIs.
        """
        self._barcode_elements = []
        self._umi_elements = []
        for el in self.ordered_elements:
            name_lower = el.element_name.lower()
            if name_lower.startswith("barcode"):
                self._barcode_elements.append(el.element_name)
            elif name_lower.startswith("umi"):
                self._umi_elements.append(el.element_name)

    @property
    def barcode_elements(self) -> List[str]:
        """Get list of barcode element names (identified by 'barcode' prefix)."""
        return self._barcode_elements

    @property
    def umi_elements(self) -> List[str]:
        """Get list of UMI element names (identified by 'umi' prefix)."""
        return self._umi_elements

    @classmethod
    def from_element_list(cls, element_list):
        """Create MoleculeStructure from a list of MoleculeElement objects."""
        ms = cls.__new__(cls)
        ms.ordered_elements = element_list
        ms.elements_to_concatenate = {}
        ms.duplicated_elements = {}
        ms._barcode_elements = []
        ms._umi_elements = []
        ms._identify_barcode_umi_elements()
        return ms

    def __iter__(self) -> Iterator[MoleculeElement]:
        for e in self.ordered_elements:
            yield e

