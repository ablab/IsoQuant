############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import os
import re
import subprocess

logger = logging.getLogger('IsoQuant')

# key, value
def get_first_best_from_sorted(sorted_list_of_pairs):
    best_value = sorted_list_of_pairs[0][1]
    result = []
    for x in sorted_list_of_pairs:
        if x[1] > best_value:
            break
        result.append(x[0])
    return result


def list_to_str(element_list, element_delim=','):
    if len(element_list) == 0:
        return "."
    return element_delim.join(list(map(str, element_list)))


def rindex(l, el):
    for i in range(len(l) - 1, -1, -1):
        if l[i] == el:
            return i
    raise ValueError(f"{el} is not in list")


def proper_plural_form(name, count):
    return str(count) + " " + name + "" if count == 1 else "s"

# check whether genes overlap and should be processed together
def genes_overlap(gene_db1, gene_db2):
    if (gene_db1.seqid != gene_db2.seqid):
        return False
    return overlaps((gene_db1.start, gene_db1.end), (gene_db2.start, gene_db2.end))

def argmin(l):
    if len(l) == 0:
        return -1
    min_v = l[0]
    min_i = 0
    for i in range(len(l)):
        if l[i] < min_v:
            min_v = l[i]
            min_i = i
    return  min_i

# == range operations ==
def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def overlaps_at_least(range1, range2, delta = 0):
    cutoff = min([delta, range1[1] - range1[0] + 1, range2[1] - range2[0] + 1])
    overlap = min(range1[1], range2[1]) + 1 - max(range1[0], range2[0])
    return overlap >= cutoff


def left_of(range1, range2):
    return range1[1] < range2[0]


def equal_ranges(range1, range2, delta=0):
    return abs(range1[0] - range2[0]) <= delta and abs(range1[1] - range2[1]) <= delta


def covers_end(bigger_range, smaller_range):
    return bigger_range[1] <= smaller_range[1] and bigger_range[0] <= smaller_range[0]


def covers_start(bigger_range, smaller_range):
    return bigger_range[0] >= smaller_range[0] and bigger_range[1] >= smaller_range[1]


def contains(bigger_range, smaller_range):
    return bigger_range[1] >= smaller_range[1] and bigger_range[0] <= smaller_range[0]


def contains_approx(bigger_range, smaller_range, delta = 1):
    return bigger_range[1] + delta >= smaller_range[1] and bigger_range[0] - delta <= smaller_range[0]


def range_list_to_str(range_list, element_delim=',', coord_delim='-'):
    return element_delim.join(list(map(lambda x: str(x[0]) + coord_delim + str(x[1]), range_list)))


# == working with alignment blocks ==
def junctions_from_blocks(sorted_blocks):
    junctions = []
    if len(sorted_blocks) >= 2:
        for i in range(0, len(sorted_blocks) - 1):
            if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
                junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


def concat_gapless_blocks(blocks, cigar_tuples):
    cigar_index = 0
    block_index = 0
    resulting_blocks = []

    current_block = None
    deletions_before_block = 0

    while cigar_index < len(cigar_tuples) and block_index < len(blocks):
        # init new block
        if current_block is None:
            # init new block from match
            if cigar_tuples[cigar_index][0] == 0:
                current_block = (blocks[block_index][0] - deletions_before_block, blocks[block_index][1])
                deletions_before_block = 0
                block_index += 1
            # keep track of deletions before matched block
            elif cigar_tuples[cigar_index][0] == 2:
                deletions_before_block = cigar_tuples[cigar_index][1]
        # found intron, add current block
        elif cigar_tuples[cigar_index][0] == 3:
            resulting_blocks.append(current_block)
            current_block = None
        # add deletion to block
        elif cigar_tuples[cigar_index][0] == 2:
            current_block = (current_block[0], current_block[1] + cigar_tuples[cigar_index][1])
        # found match - merge blocks
        elif cigar_tuples[cigar_index][0] == 0:
            if abs(current_block[1] - blocks[block_index][0]) > 1:
                logger.debug("Distant blocks")
                logger.debug(current_block, blocks[block_index])
            current_block = (current_block[0], blocks[block_index][1])

            block_index += 1
        cigar_index += 1

    if current_block is not None:
        resulting_blocks.append(current_block)

    return resulting_blocks


def correct_bam_coords(blocks):
    return list(map(lambda x: (x[0] + 1, x[1]), blocks))


# == profile functions ==
def is_subprofile(short_isoform_profile, long_isoform_profile):
    if len(short_isoform_profile) != len(long_isoform_profile):
        return None

    if all(el == -1 for el in short_isoform_profile):
        return None

    short_range_start = None
    short_range_end = None
    for i in range(len(short_isoform_profile)):
        if short_isoform_profile[i] == 1:
            if short_range_start is None:
                short_range_start = i
            short_range_end = i

    for i in range(short_range_start, short_range_end + 1):
        if short_isoform_profile[i] != long_isoform_profile[i]:
            return False
    return True


def difference_in_present_features(profile1, profile2):
    """ computes Hamming distance for two profiles, returns -1 if profiles have different length

    Parameters
    ----------
    profile1: list of int
    profile2: list of int

    Returns
    -------
    result: int
        number of different items in lists or -1 if profiles have different length

    """
    if len(profile1) != len(profile2):
        return -1
    d = 0
    for i in range(len(profile1)):
        if profile1[i] == 0 or profile2[i] == 0:
            continue
        if profile1[i] != profile2[i]:
            d += 1
    return d


def find_matching_positions(profile1, profile2):
    if len(profile1) != len(profile2):
        return -1
    matches = [0 for i in range(len(profile1))]
    for i in range(len(profile1)):
        if profile1[i] == profile2[i]:
            matches[i] = 1
    return matches


def mask_profile(read_profile, true_profile):
    if len(read_profile) != len(true_profile):
        return None
    masked_profile = []
    for i in range(len(true_profile)):
        if true_profile[i] == 1:
            masked_profile.append(read_profile[i])
    return masked_profile


def get_blocks_from_profile(features, profile):
    if len(features) != len(profile):
        return None
    profile_features = []
    for i in range(len(profile)):
        if profile[i] == 1:
            profile_features.append(features[i])
    return profile_features


# assumes there are no contradictions
def left_truncated(read_profile, isoform_profile):
    if 1 not in read_profile or 1 not in isoform_profile:
        return True
    return read_profile.index(1) > isoform_profile.index(1)


# assumes there are no contradictions
def right_truncated(read_profile, isoform_profile):
    if 1 not in read_profile or 1 not in isoform_profile:
        return True
    return rindex(read_profile, 1) < rindex(isoform_profile, 1)


def get_path_to_program(program, dirpath=None, min_version=None):
    """
    returns the path to an executable or None if it can't be found
    """
    def is_exe(fpath):
        if os.path.isfile(fpath) and os.access(fpath, os.X_OK):
            if not min_version or check_version(fpath, min_version):
                return True

    def check_version(fpath, min_version):
        p = subprocess.run([fpath, '--version'], capture_output=True)
        version_pattern = re.compile(r'(?P<major_version>\d+)\.(?P<minor_version>\d+)')
        v = version_pattern.search(str(p.stdout))
        if not v.group('major_version') or not v.group('minor_version'):
            return False
        version, minor_version = map(int, min_version.split('.'))
        if int(v.group('major_version')) == version and int(v.group('minor_version')) >= minor_version:
            return True

    if dirpath:
        exe_file = os.path.join(dirpath, program)
        if is_exe(exe_file):
            return exe_file

    for path in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
            return exe_file
    return None


