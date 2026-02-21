###########################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import math
import sys

from ..indexers import ArrayKmerIndexer, Array2BitKmerIndexer, Dict2BitKmerIndexer, KmerIndexer
from ..common import batch_str_to_2bit
from ...error_codes import IsoQuantExitCode
from .extraction_result import DetectedElement, ExtractionResult
from ..common import find_polyt, reverese_complement, detect_exact_positions, find_candidate_with_max_score_ssw
from .molecule_structure import ElementType, MoleculeStructure, MoleculeElement

logger = logging.getLogger('IsoQuant')


class UniversalSingleMoleculeExtractor:
    MIN_SCORE_COEFF = 0.75
    MIN_SCORE_COEFF_TERMMINAL = 0.5
    TERMINAL_MATCH_DELTA = 3
    MAX_LEN_DIFF = 0.25
    DEFAULT_ERROR_RATE = 0.03 # rought estimate for modern nanopore, doesn't matter that much
    MAX_HITS = 10
    SCORE_DIFF = 1

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
                if el.element_name in self.molecule_structure.elements_to_concatenate:
                    # Concatenated parts: skip individual index, handled after loop
                    pass
                elif el.element_name in self.molecule_structure.duplicated_elements:
                    # Duplicated parts: skip individual index, handled after loop
                    pass
                elif self.correct_sequences:
                    self.prepare_barcode_index_for_element(el)
                    self.elements_to_correct.add(el.element_name)
                else:
                    self.elements_to_extract.add(el.element_name)

        # Build one index per concatenated base element using the full whitelist
        if self.correct_sequences:
            for base_name, count in self.molecule_structure.concatenated_elements_counts.items():
                first_part_el = None
                total_length = 0
                for el in self.molecule_structure:
                    if el.element_name in self.molecule_structure.elements_to_concatenate:
                        if self.molecule_structure.elements_to_concatenate[el.element_name][0] == base_name:
                            if first_part_el is None:
                                first_part_el = el
                            total_length += el.element_length
                if first_part_el is not None and first_part_el.element_value is not None:
                    barcode_list = first_part_el.element_value
                    self.prepare_barcode_index(base_name, barcode_list, total_length)

        # Build one index per duplicated base element using the full whitelist
        if self.correct_sequences:
            for base_name, count in self.molecule_structure.duplicated_elements_counts.items():
                first_part_el = None
                for el in self.molecule_structure:
                    if el.element_name in self.molecule_structure.duplicated_elements:
                        if self.molecule_structure.duplicated_elements[el.element_name][0] == base_name:
                            first_part_el = el
                            break
                if first_part_el is not None and first_part_el.element_value is not None:
                    barcode_list = first_part_el.element_value
                    self.prepare_barcode_index(base_name, barcode_list, first_part_el.element_length)

        if not self.has_cdna:
            logger.critical("Molecule must include a cDNA")
            sys.exit(IsoQuantExitCode.INVALID_PARAMETER)

    def prepare_barcode_index(self, base_name, barcode_list, barcode_length):
        # Check whether all barcodes have the same length
        barcode_lengths = set(len(b) for b in barcode_list)
        variable_length = len(barcode_lengths) > 1
        if variable_length:
            logger.warning("Barcodes for element %s have variable lengths (%s), "
                           "performance may be suboptimal" %
                           (base_name, ", ".join(str(l) for l in sorted(barcode_lengths))))
            barcode_length = min(barcode_lengths)

        barcode_count = len(barcode_list)
        error_rate = self.DEFAULT_ERROR_RATE
        barcode_sparsity = math.pow(4, barcode_length) / barcode_count
        filling_edit_distance = math.ceil(math.log(barcode_sparsity, 8 * barcode_length))
        error_probability = math.pow(1 - error_rate, barcode_length - filling_edit_distance) * math.pow(error_rate, filling_edit_distance) * math.comb(barcode_length, filling_edit_distance)
        density = barcode_count * math.pow(8 * barcode_length, filling_edit_distance) / math.pow(4, barcode_length)
        if error_probability > 0.01 or density > 1:
            filling_edit_distance -= 1

        if filling_edit_distance == 0:
            self.min_scores[base_name] = barcode_length - 1
        self.min_scores[base_name] = barcode_length - filling_edit_distance
        if variable_length:
            self.min_scores[base_name] = min(self.min_scores[base_name] + 1, barcode_length)
        logger.info("Minimal score for element %s is set to %d" % (base_name, self.min_scores[base_name]))

        if barcode_count > 1000000:
            logger.warning("The number of barcodes for element %s is too large: %d, barcode calling may take substantial amount of time and RAM", (base_name, barcode_count))

        kmer_size = min(max(6, int(barcode_length / 2) - 1), 15)
        if kmer_size <= 8:
            if barcode_count < 100000 or variable_length:
                self.index_dict[base_name] = ArrayKmerIndexer(barcode_list, kmer_size)
            else:
                kmer_size += 1
                self.index_dict[base_name] = Array2BitKmerIndexer(batch_str_to_2bit(barcode_list, barcode_length),
                                                                  kmer_size, seq_len=barcode_length)
        else:
            if barcode_count < 100000 or variable_length:
                self.index_dict[base_name] = KmerIndexer(barcode_list, kmer_size)
            else:
                self.index_dict[base_name] = Dict2BitKmerIndexer(batch_str_to_2bit(barcode_list, barcode_length),
                                                                 kmer_size, seq_len=barcode_length)
        logger.info("Indexed %d barcodes for element %s" % (barcode_count, base_name))

    def prepare_barcode_index_for_element(self, element: MoleculeElement):
        if element.element_type not in (ElementType.VAR_LIST, ElementType.VAR_FILE):
            return
        barcode_list = element.element_value
        self.prepare_barcode_index(element.element_name, barcode_list, element.element_length)

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
        concat_element_storage = {}
        dupl_element_storage = {}
        self.extract_variable_elements5(sequence, detected_elements, concat_element_storage, dupl_element_storage)
        self.extract_variable_elements3(sequence, detected_elements, concat_element_storage, dupl_element_storage)
        self.process_concatenated_elements(concat_element_storage, detected_elements, sequence)
        self.process_duplicated_elements(dupl_element_storage, detected_elements, sequence)
        logger.debug("== end read id %s ==" % read_id)
        return detected_elements

    def process_concatenated_elements(self, concatenated_element_storage, detected_results, sequence):
        for base_name, parts in concatenated_element_storage.items():
            # Check all parts were detected
            if any(start == -1 for start, end in parts):
                detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
                continue

            # Concatenate sequences from all parts
            concatenated_seq = ""
            for start, end in parts:
                concatenated_seq += sequence[start:end + 1]

            # If correction is enabled and an index exists, correct against the whitelist
            if self.correct_sequences and base_name in self.index_dict:
                matching_sequences = self.index_dict[base_name].get_occurrences(concatenated_seq, max_hits=self.MAX_HITS)
                corrected_seq, seq_score, _, _ = find_candidate_with_max_score_ssw(
                    matching_sequences, concatenated_seq,
                    min_score=self.min_scores[base_name], score_diff=self.SCORE_DIFF)
                if corrected_seq is not None:
                    detected_results[base_name] = DetectedElement(
                        parts[0][0], parts[-1][1], seq_score, corrected_seq)
                else:
                    detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
            else:
                # No correction: just store the raw concatenated sequence
                detected_results[base_name] = DetectedElement(
                    parts[0][0], parts[-1][1], 0, concatenated_seq)

    def process_duplicated_elements(self, duplicated_element_storage, detected_results, sequence):
        for base_name, parts in duplicated_element_storage.items():
            # Extract sequences for detected parts
            part_seqs = []
            part_coords = []
            for start, end in parts:
                if start != -1:
                    part_seqs.append(sequence[start:end + 1])
                    part_coords.append((start, end))
                else:
                    part_seqs.append(None)
                    part_coords.append((-1, -1))

            detected_seqs = [s for s in part_seqs if s is not None]
            if not detected_seqs:
                detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
                continue

            # Use first detected copy's coordinates for reporting
            first_detected_idx = next(i for i, s in enumerate(part_seqs) if s is not None)

            if self.correct_sequences and base_name in self.index_dict:
                # Correct each copy independently
                corrected = []
                for seq in part_seqs:
                    if seq is None:
                        corrected.append((None, -1))
                        continue
                    matching = self.index_dict[base_name].get_occurrences(seq, max_hits=self.MAX_HITS)
                    corrected_seq, score, _, _ = find_candidate_with_max_score_ssw(
                        matching, seq, min_score=self.min_scores[base_name], score_diff=self.SCORE_DIFF)
                    corrected.append((corrected_seq, score))

                successful = [(seq, score) for seq, score in corrected if seq is not None]

                if not successful:
                    # Neither copy corrected
                    detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
                elif len(successful) == 1:
                    # One succeeded
                    seq, score = successful[0]
                    detected_results[base_name] = DetectedElement(
                        part_coords[first_detected_idx][0],
                        part_coords[first_detected_idx][1],
                        score, seq)
                else:
                    # Multiple succeeded - check agreement
                    unique_seqs = set(seq for seq, _ in successful)
                    if len(unique_seqs) == 1:
                        # All agree
                        best_score = max(score for _, score in successful)
                        detected_results[base_name] = DetectedElement(
                            part_coords[first_detected_idx][0],
                            part_coords[first_detected_idx][1],
                            best_score, successful[0][0])
                    else:
                        # Copies disagree after correction
                        detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
            else:
                # No correction: compare raw sequences
                if len(set(detected_seqs)) == 1:
                    # All match
                    detected_results[base_name] = DetectedElement(
                        part_coords[first_detected_idx][0],
                        part_coords[first_detected_idx][1],
                        0, detected_seqs[0])
                else:
                    # Differ: report all comma-separated
                    comma_sep = ",".join(detected_seqs)
                    detected_results[base_name] = DetectedElement(
                        part_coords[first_detected_idx][0],
                        part_coords[first_detected_idx][1],
                        0, comma_sep)

    def correct_element(self, element: MoleculeElement, sequence: str, potential_start: int, potential_end: int) -> DetectedElement:
        potential_seq = sequence[potential_start:potential_end + 1]
        matching_sequences = self.index_dict[element.element_name].get_occurrences(potential_seq, max_hits=self.MAX_HITS)
        corrected_seq, seq_score, seq_start, seq_end = \
            find_candidate_with_max_score_ssw(matching_sequences, potential_seq,
                                              min_score=self.min_scores[element.element_name], score_diff=self.SCORE_DIFF)
        if corrected_seq is not None:
            read_seq_start = potential_start + seq_start
            read_seq_end = potential_start + seq_end - 1
            return DetectedElement(read_seq_start, read_seq_end, seq_score, corrected_seq)
        else:
            return DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)

    def _backup_concatenated_element(self, element, potential_start, potential_end, element_storage):
        base_name, part_index = self.molecule_structure.elements_to_concatenate[element.element_name]
        if base_name not in element_storage:
            element_storage[base_name] = [(-1, -1) for _ in range(
                self.molecule_structure.concatenated_elements_counts[base_name])]
        # part_index is 1-based in MDF, convert to 0-based
        element_storage[base_name][part_index - 1] = (potential_start, potential_end)

    def _backup_duplicated_element(self, element, potential_start, potential_end, element_storage):
        base_name, part_index = self.molecule_structure.duplicated_elements[element.element_name]
        if base_name not in element_storage:
            element_storage[base_name] = [(-1, -1) for _ in range(
                self.molecule_structure.duplicated_elements_counts[base_name])]
        # part_index is 1-based in MDF, convert to 0-based
        element_storage[base_name][part_index - 1] = (potential_start, potential_end)

    def _detect_variable_element(self, el, sequence, potential_start, potential_end,
                                 detected_elements, concat_element_storage, dupl_element_storage,
                                 len_diff_threshold=None):
        potential_len = potential_end - potential_start + 1
        if len_diff_threshold is None or abs(potential_len - el.element_length) <= len_diff_threshold * el.element_length:
            if el.element_name in self.molecule_structure.elements_to_concatenate:
                self._backup_concatenated_element(el, potential_start, potential_end, concat_element_storage)
                return potential_start, potential_end

            if el.element_name in self.molecule_structure.duplicated_elements:
                self._backup_duplicated_element(el, potential_start, potential_end, dupl_element_storage)
                return potential_start, potential_end

            if el.element_name in self.elements_to_extract:
                detected_elements[el.element_name] = DetectedElement(potential_start, potential_end, 0,
                                                                     seq=sequence[potential_start:potential_end + 1])

            elif el.element_name in self.elements_to_correct:
                detected_element = self.correct_element(el, sequence, potential_start, potential_end)
                if detected_element.start != -1:
                    detected_elements[el.element_name] = detected_element
                    return detected_element.start, detected_element.end
        else:
            detected_elements[el.element_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)

        return potential_start, potential_end

    def _detect_potential_start(self, element, element_index, potential_end, detected_elements):
        if element_index == 0:
            potential_start = potential_end - element.element_length + 1
        else:
            prev_el = self.molecule_structure.ordered_elements[element_index - 1]
            if prev_el.element_name in detected_elements:
                potential_start = detected_elements[prev_el].end + 1
            else:
                potential_start = potential_end - element.element_length + 1
        if potential_start < 0:
            potential_start = 0

        return potential_start

    def _detect_potential_end(self, sequence, element, element_index, potential_start, detected_elements):
        if element_index + 1 == len(self.molecule_structure.ordered_elements):
            potential_end = potential_start + element.element_length - 1
        else:
            next_el = self.molecule_structure.ordered_elements[element_index + 1]
            if next_el.element_name in detected_elements:
                potential_end = detected_elements[next_el.element_name].start - 1
            else:
                potential_end = potential_start + element.element_length - 1
        if potential_end >= len(sequence):
            potential_end = len(sequence) - 1
        return potential_end

    def extract_variable_elements5(self, sequence, detected_elements, concat_element_storage, dupl_element_storage):
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

            el_start, _ = self._detect_variable_element(el, sequence, potential_start, potential_end,
                                                        detected_elements, concat_element_storage, dupl_element_storage)
            current_pos = el_start - 1

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
            potential_end = self._detect_potential_end(sequence, el, i, potential_start, detected_elements)

            _, el_end = self._detect_variable_element(el, sequence, potential_start, potential_end,
                                                      detected_elements, concat_element_storage, dupl_element_storage,
                                                      len_diff_threshold=self.MAX_LEN_DIFF)
            current_pos = el_end + 1


    def extract_variable_elements3(self, sequence, detected_elements, concat_element_storage, dupl_element_storage):
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

            _, el_end = self._detect_variable_element(el, sequence, potential_start, potential_end,
                                                      detected_elements, concat_element_storage, dupl_element_storage)
            current_pos = el_end + 1

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
            potential_start = self._detect_potential_start(el, i, potential_end, detected_elements)

            el_start, _ = self._detect_variable_element(el, sequence, potential_start, potential_end,
                                                        detected_elements, concat_element_storage, dupl_element_storage,
                                                        len_diff_threshold=self.MAX_LEN_DIFF)
            current_pos = el_start - 1

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
            element_occurrences = self.index_dict[el.element_name].get_occurrences_substr(sequence,
                                                                                          search_start,
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
