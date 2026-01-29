###########################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import sys

from ...error_codes import IsoQuantExitCode
from .extraction_result import DetectedElement, ExtractionResult
from ..indexers.base import KmerIndexer
from ..common import find_polyt, reverese_complement, detect_exact_positions, find_candidate_with_max_score_ssw
from .molecule_structure import ElementType, MoleculeStructure, MoleculeElement

logger = logging.getLogger('IsoQuant')


class UniversalSingleMoleculeExtractor:
    MIN_SCORE_COEFF = 0.75
    MIN_SCORE_COEFF_TERMMINAL = 0.5
    TERMINAL_MATCH_DELTA = 3
    MAX_LEN_DIFF = 0.25

    def __init__(self, molecule_structure: MoleculeStructure):
        self.correct_sequences = True
        self.molecule_structure = molecule_structure
        self.index_dict = {}
        self.min_scores = {}
        self.has_polyt = False
        self.has_cdna = False
        self.constant_elements_to_detect = set()
        self.elements_to_extract = set()
        self.elements_to_correct = set()

        for el in self.molecule_structure:
            if el.element_type == ElementType.PolyT:
                if self.has_polyt:
                    logger.critical(
                        "Current version only supports a single polyT, extraction results may be suboptimal")
                    sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
                self.has_polyt = True
            elif el.element_type == ElementType.cDNA:
                if self.has_cdna:
                    logger.critical(
                        "Current version only supports a single cDNA, reads will not be split into molecules")
                    sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
                self.has_cdna = True
            elif el.element_type.is_constant():
                self.constant_elements_to_detect.add(el.element_name)
                self.index_dict[el.element_name] = KmerIndexer([el.element_value],
                                                               max(6, int(el.element_length / 2) - 2))
            elif el.element_type.needs_sequence_extraction():
                self.elements_to_extract.add(el.element_name)
            elif el.element_type.needs_correction():
                if self.correct_sequences:
                    self.index_dict[el.element_name] = self.prepare_barcode_index(el)
                    # FIXME: min score should be based on the length of the barcode and their sparcity
                    # e.g. find a way to estimate min_score so that error rate does not exceed 1%
                    self.min_scores[el.element_name] = el.element_length - 2
                    self.elements_to_correct.add(el.element_name)
                else:
                    self.elements_to_extract.add(el.element_name)

        if not self.has_cdna:
            logger.critical("Molecule must include a cDNA")
            sys.exit(IsoQuantExitCode.INVALID_PARAMETER)

    def prepare_barcode_index(self, element: MoleculeElement):
        if element.element_type not in (ElementType.VAR_LIST, ElementType.VAR_FILE):
            return None
        barcode_list = element.element_value

        # TODO: use optimal indices for long barocde lists
        return KmerIndexer(barcode_list, max(6, int(element.element_length / 2) - 2))

    def result_type(self):
        return ExtractionResult

    def header(self):
        """Get header matching ExtractionResult output format."""
        # Create a dummy result to get the correct header format
        dummy_result = ExtractionResult(self.molecule_structure, "", ".", {})
        return dummy_result.header()

    def find_barcode_umi(self, read_id, sequence):
        detected_elements_fwd = self._find_patterns_fwd(read_id, sequence)
        rev_seq = reverese_complement(sequence)
        detected_elements_rev = self._find_patterns_fwd(read_id, rev_seq)

        if self.has_polyt:
            fwd_has_polyt = ElementType.PolyT.name in detected_elements_fwd
            rev_has_polyt = ElementType.PolyT.name in detected_elements_rev
            if fwd_has_polyt and not rev_has_polyt:
                return ExtractionResult(self.molecule_structure, read_id, '+', detected_elements_fwd)
            if not fwd_has_polyt and rev_has_polyt:
                return ExtractionResult(self.molecule_structure, read_id, '-', detected_elements_rev)

        if len(detected_elements_fwd) >= len(detected_elements_rev):
            return ExtractionResult(self.molecule_structure, read_id, '+', detected_elements_fwd)
        return ExtractionResult(self.molecule_structure, read_id, '-', detected_elements_rev)

    def _find_patterns_fwd(self, read_id, sequence):
        logger.debug("== read id %s ==" % read_id)
        # TODO: make a list of DetectedElement rather that dict of str->DetectedElement for speed
        detected_elements = {}
        polyt_start, polyt_end = -1, -1
        if self.has_polyt:
            polyt_start, polyt_end = find_polyt(sequence)
            detected_elements[ElementType.PolyT.name] = DetectedElement(polyt_start, polyt_end, 0)
        logger.debug("PolyT %d" % polyt_start)

        self.detect_const_elements(sequence, polyt_start, polyt_end, detected_elements)
        self.extract_variable_elements5(detected_elements, sequence)
        self.extract_variable_elements3(detected_elements, sequence)
        logger.debug("== end read id %s ==" % read_id)
        return detected_elements

    def concatenate_elements(self, detected_results):
        pass

    def extract_variable_elements5(self, detected_elements, sequence):
        first_detected_const_element = None
        for i, el in enumerate(self.molecule_structure):
            if el.element_type == ElementType.cDNA:
                break
            if (el.element_type == ElementType.PolyT or el.element_type.is_constant()) and el.element_name in detected_elements:
                first_detected_const_element = i
                break

        if first_detected_const_element is None: return

        # extracting elements preceding first detected const element
        current_pos = detected_elements[
                          self.molecule_structure.ordered_elements[first_detected_const_element].element_name].start - 1
        for i in range(first_detected_const_element - 1, -1, -1):
            el: MoleculeElement = self.molecule_structure.ordered_elements[i]

            if el.element_type.is_base_separator():
                current_pos -= el.element_length
                continue
            if not el.element_type.is_variable():
                # we do not skip undetected const elements on the open side (e.g. here we go to the left of the first detected const element)
                break

            potential_end = current_pos
            potential_start = potential_end - el.element_length + 1
            if potential_start < 0: break

            current_pos = potential_start - 1
            if el.element_name in self.elements_to_extract:
                detected_elements[el.element_name] = DetectedElement(potential_start, potential_end, 0,
                                                                     seq=sequence[potential_start:potential_end + 1])

            elif el.element_name in self.elements_to_correct:
                potential_seq = sequence[potential_start:potential_end + 1]
                matching_sequences = self.index_dict[el.element_name].get_occurrences(potential_seq)
                corrected_seq, seq_score, seq_start, seq_end = \
                    find_candidate_with_max_score_ssw(matching_sequences, potential_seq, min_score=self.min_scores[el.element_name])

                if corrected_seq is not None:
                    read_seq_start = potential_start + seq_start
                    read_seq_end = potential_start + seq_end - 1
                    current_pos = read_seq_start - 1
                    detected_elements[el.element_name] = DetectedElement(read_seq_start, read_seq_end, seq_score,
                                                                         corrected_seq)
                else:
                    detected_elements[el.element_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)

        current_pos = detected_elements[
                          self.molecule_structure.ordered_elements[first_detected_const_element].element_name].end + 1
        for i in range(first_detected_const_element + 1, len(self.molecule_structure.ordered_elements)):
            if current_pos >= len(sequence):
                break

            el: MoleculeElement = self.molecule_structure.ordered_elements[i]
            if el.element_type in [ElementType.cDNA, ElementType.PolyT]:
                break
            elif el.element_type.is_base_separator():
                current_pos += el.element_length
                continue
            elif el.element_type.is_constant():
                if el.element_name in detected_elements:
                    current_pos = detected_elements[el.element_name].end + 1
                else:
                    current_pos += el.element_length
                continue

            potential_start = current_pos
            if i + 1 == len(self.molecule_structure.ordered_elements):
                potential_end = potential_start + el.element_length - 1
            else:
                next_el = self.molecule_structure.ordered_elements[i + 1]
                if next_el.element_name in detected_elements:
                    potential_end = detected_elements[next_el.element_name].start - 1
                else:
                    potential_end = potential_start + el.element_length - 1
            if potential_end >= len(sequence):
                potential_end = len(sequence) - 1

            potential_len = potential_end - potential_start + 1
            current_pos = potential_end + 1
            if abs(potential_len - el.element_length) <= self.MAX_LEN_DIFF * el.element_length:
                if el.element_name in self.elements_to_extract:
                    detected_elements[el.element_name] = DetectedElement(potential_start, potential_end, 0,
                                                                         seq=sequence[
                                                                             potential_start:potential_end + 1])

                elif el.element_name in self.elements_to_correct:
                    potential_seq = sequence[potential_start:potential_end + 1]
                    matching_sequences = self.index_dict[el.element_name].get_occurrences(potential_seq)
                    corrected_seq, seq_score, seq_start, seq_end = \
                        find_candidate_with_max_score_ssw(matching_sequences, potential_seq,
                                                          min_score=self.min_scores[el.element_name])

                    if corrected_seq is not None:
                        read_seq_start = potential_start + seq_start
                        read_seq_end = potential_start + seq_end - 1
                        current_pos = read_seq_end + 1
                        detected_elements[el.element_name] = DetectedElement(read_seq_start, read_seq_end, seq_score,
                                                                             corrected_seq)
                    else:
                        detected_elements[el.element_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
            else:
                detected_elements[el.element_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)


    def extract_variable_elements3(self, detected_elements, sequence):
        last_detected_const_element = None
        for i in range(len(self.molecule_structure.ordered_elements) - 1, -1, -1):
            el = self.molecule_structure.ordered_elements[i]
            if el.element_type == ElementType.cDNA:
                break
            if (el.element_type == ElementType.PolyT or el.element_type.is_constant()) and el.element_name in detected_elements:
                last_detected_const_element = i
                break

        if last_detected_const_element is None: return

        # extracting elements following the last detected const element
        current_pos = detected_elements[
                          self.molecule_structure.ordered_elements[last_detected_const_element].element_name].end + 1
        for i in range(last_detected_const_element + 1, len(self.molecule_structure.ordered_elements)):
            el: MoleculeElement = self.molecule_structure.ordered_elements[i]

            if el.element_type.is_base_separator():
                current_pos += el.element_length
                continue
            if not el.element_type.is_variable():
                # we do not skip undetected const elements on the open side (e.g. here we go to the left of the first detected const element)
                break

            potential_start = current_pos
            potential_end = potential_start + el.element_length - 1
            if potential_end >= len(sequence): break
            current_pos = potential_end + 1

            if el.element_name in self.elements_to_extract:
                detected_elements[el.element_name] = DetectedElement(potential_start, potential_end, 0,
                                                                     seq=sequence[potential_start:potential_end + 1])

            elif el.element_name in self.elements_to_correct:
                potential_seq = sequence[potential_start:potential_end + 1]
                matching_sequences = self.index_dict[el.element_name].get_occurrences(potential_seq)
                corrected_seq, seq_score, seq_start, seq_end = \
                    find_candidate_with_max_score_ssw(matching_sequences, potential_seq, min_score=self.min_scores[el.element_name])

                if corrected_seq is not None:
                    read_seq_start = potential_start + seq_start
                    read_seq_end = potential_start + seq_end - 1
                    current_pos = read_seq_end + 1
                    detected_elements[el.element_name] = DetectedElement(read_seq_start, read_seq_end, seq_score,
                                                                         corrected_seq)
                else:
                    detected_elements[el.element_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)

        current_pos = detected_elements[
                          self.molecule_structure.ordered_elements[last_detected_const_element].element_name].start - 1
        for i in range(last_detected_const_element - 1, -1, -1):
            if current_pos <= 0:
                break
            el: MoleculeElement = self.molecule_structure.ordered_elements[i]

            if el.element_type in [ElementType.cDNA, ElementType.PolyT]:
                break
            elif el.element_type.is_base_separator():
                current_pos -= el.element_length
                continue
            elif el.element_type.is_constant():
                if el.element_name in detected_elements:
                    current_pos = detected_elements[el.element_name].start - 1
                else:
                    current_pos -= el.element_length
                continue

            potential_end = current_pos
            if i == 0:
                potential_start = potential_end - el.element_length + 1
            else:
                prev_el = self.molecule_structure.ordered_elements[i - 1]
                if prev_el.element_name in detected_elements:
                    potential_start = detected_elements[prev_el].end + 1
                else:
                    potential_start = potential_end - el.element_length + 1
            if potential_start < 0:
                potential_start = 0

            potential_len = potential_end - potential_start + 1
            current_pos = potential_start - 1
            if abs(potential_len - el.element_length) <= self.MAX_LEN_DIFF * el.element_length:
                if el.element_name in self.elements_to_extract:
                    detected_elements[el.element_name] = DetectedElement(potential_start, potential_end, 0,
                                                                         seq=sequence[
                                                                             potential_start:potential_end + 1])

                elif el.element_name in self.elements_to_correct:
                    potential_seq = sequence[potential_start:potential_end + 1]
                    matching_sequences = self.index_dict[el.element_name].get_occurrences(potential_seq)
                    corrected_seq, seq_score, seq_start, seq_end = \
                        find_candidate_with_max_score_ssw(matching_sequences, potential_seq,
                                                          min_score=self.min_scores[el.element_name])

                    if corrected_seq is not None:
                        read_seq_start = potential_start + seq_start
                        read_seq_end = potential_start + seq_end - 1
                        current_pos = read_seq_start - 1
                        detected_elements[el.element_name] = DetectedElement(read_seq_start, read_seq_end, seq_score,
                                                                             corrected_seq)
                    else:
                        detected_elements[el.element_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
            else:
                detected_elements[el.element_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)


    def detect_const_elements(self, sequence, polyt_start, polyt_end, detected_elements):
        # searching left of cDNA
        first_element_detected = False
        current_search_start = 0
        for el in self.molecule_structure:
            if el.element_type in [ElementType.cDNA, ElementType.PolyT]:
                break
            if el.element_name not in self.constant_elements_to_detect:
                continue

            search_start = current_search_start
            search_end = len(sequence) if polyt_start == -1 else polyt_start + 1
            element_occurrences = self.index_dict[el.element_name].get_occurrences_substr(sequence, search_start,
                                                                                          search_end)
            if not first_element_detected:
                min_score = int(el.element_length * self.MIN_SCORE_COEFF_TERMMINAL)
                start_delta = -1
            else:
                min_score = int(el.element_length * self.MIN_SCORE_COEFF)
                start_delta = self.TERMINAL_MATCH_DELTA

            element_start, element_end = detect_exact_positions(sequence,
                                                                search_start, search_end,
                                                                self.index_dict[el.element_name].k,
                                                                el.element_value,
                                                                element_occurrences,
                                                                min_score=min_score,
                                                                start_delta=start_delta,
                                                                end_delta=self.TERMINAL_MATCH_DELTA)

            if element_start is not None:
                current_search_start = element_end
                if not first_element_detected:
                    first_element_detected = True

                detected_elements[el.element_name] = DetectedElement(element_start, element_end)

        last_element_detected = False
        current_search_end = len(sequence)
        for i in range(len(self.molecule_structure.ordered_elements) - 1, -1, -1):
            el = self.molecule_structure.ordered_elements[i]
            if el.element_type in [ElementType.cDNA, ElementType.PolyT]:
                break
            if el.element_name not in self.constant_elements_to_detect:
                continue

            search_start = 0 if polyt_start == -1 else polyt_end + 1
            search_end = current_search_end
            element_occurrences = self.index_dict[el.element_name].get_occurrences_substr(sequence, search_start,
                                                                                          search_end)

            if not last_element_detected:
                min_score = int(el.element_length * self.MIN_SCORE_COEFF_TERMMINAL)
                end_delta = -1
            else:
                min_score = int(el.element_length * self.MIN_SCORE_COEFF)
                end_delta = self.TERMINAL_MATCH_DELTA

            element_start, element_end = detect_exact_positions(sequence,
                                                                search_start, search_end,
                                                                self.index_dict[el.element_name].k,
                                                                el.element_value,
                                                                element_occurrences,
                                                                min_score=min_score,
                                                                start_delta=self.TERMINAL_MATCH_DELTA,
                                                                end_delta=end_delta)
            if element_start is not None:
                current_search_end = element_end
                if not last_element_detected:
                    last_element_detected = True

                detected_elements[el.element_name] = DetectedElement(element_start, element_end)
