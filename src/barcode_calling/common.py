############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from ssw import AlignmentMgr


def find_polyt_start(seq, window_size = 16, polya_fraction = 0.75):
    polyA_count = int(window_size * polya_fraction)

    if len(seq) < window_size:
        return -1
    i = 0
    a_count = seq[0:window_size].count('T')
    while i < len(seq) - window_size:
        if a_count >= polyA_count:
            break
        first_base_a = seq[i] == 'T'
        new_base_a = i + window_size < len(seq) and seq[i + window_size] == 'T'
        if first_base_a and not new_base_a:
            a_count -= 1
        elif not first_base_a and new_base_a:
            a_count += 1
        i += 1

    if i >= len(seq) - window_size:
        return -1

    return i + max(0, seq[i:].find('TTT'))


base_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', " ": " "}


def reverese_complement(my_seq):  ## obtain reverse complement of a sequence
    lms = list(map(lambda x: base_comp[x], my_seq))[::-1]
    return ''.join(lms)


def align_pattern_ssw(sequence, start, end, pattern, min_score=0):
    seq = sequence[start:end]
    align_mgr = AlignmentMgr(match_score=1, mismatch_penalty=1)
    align_mgr.set_read(pattern)
    align_mgr.set_reference(seq)
    alignment = align_mgr.align(gap_open=1, gap_extension=1)
    if alignment.optimal_score < min_score:
        return None, None, None, None, None
    return start + alignment.reference_start, start + alignment.reference_end, \
        alignment.read_start, alignment.read_end, alignment.optimal_score


def find_candidate_with_max_score_ssw(barcode_matches: list, read_sequence, min_score=10, score_diff=0, sufficient_score=0):
    best_match = [0, 0, 0]
    best_barcode = None
    second_best_score = 0

    align_mgr = AlignmentMgr(match_score=1, mismatch_penalty=1)
    align_mgr.set_reference(read_sequence)
    for barcode_match in barcode_matches:
        barcode = barcode_match[0]
        align_mgr.set_read(barcode)
        alignment = align_mgr.align(gap_open=1, gap_extension=1)
        if alignment.optimal_score < min_score:
            continue

        if alignment.optimal_score > best_match[0]:
            best_barcode = barcode
            second_best_score = best_match[0]
            best_match[0] = alignment.optimal_score
            best_match[1] = alignment.reference_start
            best_match[2] = alignment.reference_end
        elif alignment.optimal_score  == best_match[0] and alignment.reference_start < best_match[1]:
            best_barcode = barcode
            second_best_score = best_match[0]
            best_match[1] = alignment.reference_start
            best_match[2] = alignment.reference_end

        if alignment.optimal_score > sufficient_score > 0:
            # dirty hack to select first "sufficiently good" alignment
            break

    if best_match[0] - second_best_score < score_diff:
        return None, 0, 0, 0

    return best_barcode, best_match[0], best_match[1], best_match[2]


def detect_exact_positions(sequence, start, end, kmer_size, pattern, pattern_occurrences: list,
                           min_score=0, start_delta=-1, end_delta=-1):
    pattern_index = None
    for i, p in enumerate(pattern_occurrences):
        if p[0] == pattern:
            pattern_index = i
            break
    if pattern_index is None:
        return None, None

    start_pos, end_pos, pattern_start, pattern_end, score  = None, None, None, None, 0
    last_potential_pos = -2*len(pattern)
    for match_position in pattern_occurrences[pattern_index][2]:
        if match_position - last_potential_pos < len(pattern):
            continue

        potential_start = start + match_position - len(pattern) + kmer_size
        potential_start = max(start, potential_start)
        potential_end = start + match_position + len(pattern) + 1
        potential_end = min(end, potential_end)
        alignment = \
            align_pattern_ssw(sequence, potential_start, potential_end, pattern, min_score)
        if alignment[4] is not None and alignment[4] > score:
            start_pos, end_pos, pattern_start, pattern_end, score = alignment

    if start_pos is None:
        return None, None

    if start_delta >= 0 and pattern_start > start_delta:
        return None, None
    if end_delta >= 0 and len(pattern) - pattern_end - 1 > end_delta:
        return None, None
    leftover_bases = len(pattern) - pattern_end - 1
    skipped_bases = pattern_start
    return start_pos - skipped_bases, end_pos + leftover_bases


NUCL2BIN = {'A': 0, 'C': 1, 'G': 3, 'T': 2, 'a': 0, 'c': 1, 'g': 3, 't': 2}
BIN2NUCL = ["A", "C", "T", "G"]


def str_to_2bit(seq):
    kmer_idx = 0
    seq_len = len(seq)
    for i in range(seq_len):
        kmer_idx |= ((ord(seq[i]) & 6) >> 1) << ((seq_len - i - 1) * 2)
    return kmer_idx


def bit_to_str(seq, seq_len):
    str_seq = ""
    for i in range(seq_len):
        str_seq += BIN2NUCL[(seq >> ((seq_len - i - 1) * 2)) & 3]
    return str_seq