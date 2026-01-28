###########################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import sys
from enum import unique, Enum
from typing import Iterator

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
    CONCAT_DELIM = "|"
    DUPL_DELIM = "/"

    def __init__(self, str_iterator):
        self.ordered_elements: list[MoleculeElement] = []
        self.elements_to_concatenate = {}
        self.duplicated_elements = {}

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

        for el in elements:
            if el not in element_properties:
                if el not in ElementType.__dict__:
                    logger.critical("Molecule element %s was not described in the format file" % el)
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
                element_type = ElementType[el]
                self.ordered_elements.append(MoleculeElement(el, element_type))
            else:
                element_type, element_val1, element_val2 = element_properties[el]
                if element_type not in ElementType.__dict__:
                    logger.critical("Molecule element type %s is not among the possible types" % element_type)
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
                element_type = ElementType[element_type]
                self.ordered_elements.append(MoleculeElement(el, element_type, element_val1, element_val2))

            if self.CONCAT_DELIM in el:
                v = el.split(self.CONCAT_DELIM)
                if len(v) != 2:
                    logger.critical("Incorrect concatenated element %s" % el)
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
                self.elements_to_concatenate[el] = v[0]
            elif self.DUPL_DELIM in el:
                v = el.split(self.DUPL_DELIM)
                if len(v) != 2:
                    logger.critical("Incorrect duplicated element %s" % el)
                    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)
                self.duplicated_elements[el] = v[0]

    @classmethod
    def from_element_list(cls, element_list):
        ms = cls.__new__(cls)
        ms.ordered_elements = element_list
        return ms

    def __iter__(self) -> Iterator[MoleculeElement]:
        for e in self.ordered_elements:
            yield e

    def header(self):
        printed_elements = set()
        header = "#read_id\tstrand"
        for el in self.ordered_elements:
            if el.element_type == ElementType.cDNA: continue
            if el.element_type == ElementType.PolyT:
                header += "\tpolyT_start\tpolyT_end"
            elif el.element_type.is_constant():
                header += "\t%s_start\t%s_end" % (el.element_name, el.element_name)
            elif el.element_type.needs_sequence_extraction() or el.element_type.needs_correction():
                if el.element_name in self.elements_to_concatenate.keys():
                    if self.elements_to_concatenate[el.element_name] in printed_elements:
                        continue
                    printed_elements.add(self.elements_to_concatenate[el.element_name])
                    header += "\t%s_start\t%s_end\t%s_sequence" % (el.element_name, el.element_name,
                                                                   self.elements_to_concatenate[el.element_name])
                else:
                    header += "\t%s_start\t%s_end\t%s_sequence" % (el.element_name, el.element_name, el.element_name)
        return header
