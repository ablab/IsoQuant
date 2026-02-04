############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import copy
import os
import sys
import logging
import gzip
import numpy as np
from ssw import AlignmentMgr

from ..error_codes import IsoQuantExitCode


logger = logging.getLogger('IsoQuant')

# Try to import numba for JIT compilation of barcode conversion
try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except (ImportError, SystemError):
    NUMBA_AVAILABLE = False


if NUMBA_AVAILABLE:
    @njit(cache=True, parallel=True)
    def _compute_2bit_numba(encoded, shifts, n_barcodes, seq_len):
        """Numba JIT-compiled 2-bit encoding computation with parallel execution."""
        result = np.zeros(n_barcodes, dtype=np.uint64)
        for i in prange(n_barcodes):
            val = np.uint64(0)
            for pos in range(seq_len):
                val |= np.uint64(encoded[i, pos]) << shifts[pos]
            result[i] = val
        return result


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

    return i + max(0, seq[i:].find('TTTT'))


def find_polyt(seq, window_size = 16, polya_fraction = 0.75):
    polyA_count = int(window_size * polya_fraction)

    if len(seq) < window_size:
        return -1, -1
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
        return -1, -1
    polyt_start = i + max(0, seq[i:].find('TTTT'))

    while i < len(seq) - window_size:
        if a_count < polyA_count:
            break
        first_base_a = seq[i] == 'T'
        new_base_a = i + window_size < len(seq) and seq[i + window_size] == 'T'
        if first_base_a and not new_base_a:
            a_count -= 1
        elif not first_base_a and new_base_a:
            a_count += 1
        i += 1

    if i >= len(seq) - window_size:
        polyt_end = len(seq) - 1
    else:
        last_t_pos = seq[i:i + window_size].rfind('TTTT')
        polyt_end = i + window_size - 1 if last_t_pos == -1 else i + last_t_pos + 4

    return polyt_start, polyt_end


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
            best_match[1] = alignment.reference_start - alignment.read_start
            best_match[2] = alignment.reference_end + (len(barcode) - alignment.read_end)
        elif alignment.optimal_score  == best_match[0] and alignment.reference_start < best_match[1]:
            best_barcode = barcode
            second_best_score = best_match[0]
            best_match[1] = alignment.reference_start - alignment.read_start
            best_match[2] = alignment.reference_end + (len(barcode) - alignment.read_end)

        if alignment.optimal_score > sufficient_score > 0:
            # dirty hack to select first "sufficiently good" alignment
            break

    if best_match[0] - second_best_score < score_diff:
        return None, 0, 0, 0

    return best_barcode, best_match[0], best_match[1], best_match[2]


def find_candidate_with_max_score_ssw_var_len(barcode_matches: list, read_sequence, min_score=14, score_diff=1):
    best_match = [0, 0, 0, 0]
    second_best_match = [0, 0, 0, 0]
    best_barcode = None

    align_mgr = AlignmentMgr(match_score=1, mismatch_penalty=1)
    align_mgr.set_reference(read_sequence)
    for barcode_match in barcode_matches:
        barcode = barcode_match[0]
        align_mgr.set_read(barcode)
        alignment = align_mgr.align(gap_open=1, gap_extension=1)
        if alignment.optimal_score < min_score:
            continue

        ed = len(barcode) - alignment.optimal_score
        if alignment.optimal_score > best_match[0]:
            second_best_match = copy.copy(best_match)
            best_barcode = barcode
            best_match[0] = alignment.optimal_score
            best_match[1] = alignment.reference_start - alignment.read_start
            best_match[2] = alignment.reference_end + (len(barcode) - alignment.read_end)
            best_match[3] = ed
        elif alignment.optimal_score == best_match[0] and ed <= best_match[3]:
            second_best_match = copy.copy(best_match)
            best_barcode = barcode
            best_match[1] = alignment.reference_start - alignment.read_start
            best_match[2] = alignment.reference_end + (len(barcode) - alignment.read_end)
            best_match[3] = ed

    if best_barcode and best_match[0] < len(best_barcode) and best_match[0] - second_best_match[0] < score_diff:
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
    return max(0, start_pos - skipped_bases), min(end_pos + leftover_bases, len(sequence) - 1)


def detect_first_exact_positions(sequence, start, end, kmer_size, pattern, pattern_occurrences: list,
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
        if alignment[4] is not None:
            start_pos, end_pos, pattern_start, pattern_end, score = alignment
            break

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


def batch_str_to_2bit(barcodes, seq_len):
    """
    Convert a list/iterator of barcode strings to 2-bit encoded integers.

    Uses NumPy vectorization for ~50-100x speedup over individual str_to_2bit calls.
    For 500M barcodes, this reduces conversion time from hours to minutes.

    Args:
        barcodes: Iterable of barcode strings (all must be same length)
        seq_len: Length of each barcode (default 25 for Stereo-seq)

    Returns:
        numpy.ndarray of uint64 containing 2-bit encoded barcodes
    """
    # Convert iterator to list if needed, reading in chunks to show progress
    if not isinstance(barcodes, (list, np.ndarray)):
        logger.info("Loading barcodes from iterator...")
        barcode_list = []
        chunk_size = 10_000_000
        for bc in barcodes:
            barcode_list.append(bc)
            if len(barcode_list) % chunk_size == 0:
                logger.info("  Loaded %d barcodes..." % len(barcode_list))
        barcodes = barcode_list
        logger.info("  Loaded %d barcodes total" % len(barcodes))

    n_barcodes = len(barcodes)
    if n_barcodes == 0:
        return np.array([], dtype=np.uint64)

    logger.info("Converting %d barcodes to 2-bit encoding..." % n_barcodes)

    # Fast path: join all barcodes into single string, convert to bytes at once
    # This avoids Python loop overhead for encode()
    all_bytes = ''.join(barcodes).encode('ascii')
    byte_array = np.frombuffer(all_bytes, dtype=np.uint8).reshape(n_barcodes, seq_len)

    # Apply the encoding: (ord(char) & 6) >> 1
    # A=65 -> 0, C=67 -> 1, T=84 -> 2, G=71 -> 3
    encoded = (byte_array & 6) >> 1

    # Compute bit shifts for each position (highest position = leftmost)
    shifts = np.arange(seq_len - 1, -1, -1, dtype=np.uint64) * 2

    # Compute 2-bit encoding - use Numba if available for parallel execution
    if NUMBA_AVAILABLE:
        result = _compute_2bit_numba(encoded, shifts, n_barcodes, seq_len)
    else:
        # Fallback: vectorized NumPy loop
        result = np.zeros(n_barcodes, dtype=np.uint64)
        for pos in range(seq_len):
            result |= encoded[:, pos].astype(np.uint64) << shifts[pos]

    logger.info("Barcode conversion complete")
    return result


def batch_str_to_2bit_chunked(barcode_iterator, seq_len, chunk_size=50_000_000):
    """
    Convert barcodes to 2-bit encoding in memory-efficient chunks.

    Processes barcodes in chunks to avoid loading all strings into memory at once.
    Returns a single concatenated numpy array.

    Args:
        barcode_iterator: Iterator yielding barcode strings
        seq_len: Length of each barcode
        chunk_size: Number of barcodes to process per chunk

    Returns:
        numpy.ndarray of uint64 containing 2-bit encoded barcodes
    """
    all_results = []
    current_chunk = []
    total_processed = 0

    logger.info("Converting barcodes to 2-bit encoding (chunked mode)...")

    for bc in barcode_iterator:
        current_chunk.append(bc)
        if len(current_chunk) >= chunk_size:
            # Process this chunk
            chunk_result = _convert_chunk_to_2bit(current_chunk, seq_len)
            all_results.append(chunk_result)
            total_processed += len(current_chunk)
            logger.info("  Processed %d barcodes..." % total_processed)
            current_chunk = []

    # Process remaining barcodes
    if current_chunk:
        chunk_result = _convert_chunk_to_2bit(current_chunk, seq_len)
        all_results.append(chunk_result)
        total_processed += len(current_chunk)

    logger.info("Barcode conversion complete: %d barcodes" % total_processed)

    if not all_results:
        return np.array([], dtype=np.uint64)
    return np.concatenate(all_results)


def _convert_chunk_to_2bit(barcodes, seq_len):
    """Convert a chunk of barcodes to 2-bit encoding (internal helper)."""
    n_barcodes = len(barcodes)
    if n_barcodes == 0:
        return np.array([], dtype=np.uint64)

    # Fast path: join all barcodes into single string, convert to bytes at once
    all_bytes = ''.join(barcodes).encode('ascii')
    byte_array = np.frombuffer(all_bytes, dtype=np.uint8).reshape(n_barcodes, seq_len)

    # Apply encoding: (ord(char) & 6) >> 1
    encoded = (byte_array & 6) >> 1

    # Compute shifts
    shifts = np.arange(seq_len - 1, -1, -1, dtype=np.uint64) * 2

    # Compute result - use Numba if available for parallel execution
    if NUMBA_AVAILABLE:
        result = _compute_2bit_numba(encoded, shifts, n_barcodes, seq_len)
    else:
        result = np.zeros(n_barcodes, dtype=np.uint64)
        for pos in range(seq_len):
            result |= encoded[:, pos].astype(np.uint64) << shifts[pos]

    return result


def bit_to_str(seq, seq_len):
    # Convert numpy scalar to Python int for bitwise operations
    seq = int(seq)
    str_seq = ""
    for i in range(seq_len):
        str_seq += BIN2NUCL[(seq >> ((seq_len - i - 1) * 2)) & 3]
    return str_seq



def load_h5_barcodes_bit(h5_file_path, dataset_name='bpMatrix_1'):
    raise NotImplementedError()
    import h5py
    barcode_list = []
    with h5py.File(h5_file_path, 'r') as h5_file:
        dataset = numpy.array(h5_file[dataset_name])
        for row in dataset:
            for col in row:
                barcode_list.append(bit_to_str(int(col[0])))
    return barcode_list


def load_barcodes(inf, needs_iterator=False):
    if not os.path.isfile(inf):
        logger.critical("Barcode file '%s' does not exist.", inf)
        sys.exit(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)

    if inf.endswith("h5") or inf.endswith("hdf5"):
        return load_h5_barcodes_bit(inf)

    if inf.endswith("gz") or inf.endswith("gzip"):
        handle = gzip.open(inf, "rt")
    else:
        handle = open(inf, "r")

    barcode_iterator = iter(l.strip().split()[0] for l in handle)
    if needs_iterator:
        return barcode_iterator

    return list(barcode_iterator)
