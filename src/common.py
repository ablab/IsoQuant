############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import os
import re
import subprocess
import threading
import math
from collections import defaultdict

logger = logging.getLogger('IsoQuant')


class AtomicCounter(object):
    def __init__(self):
        self.value = 0
        self._lock = threading.Lock()

    def increment(self):
        with self._lock:
            self.value += 1
            return self.value


# key, value
def get_first_best_from_sorted(sorted_list_of_pairs):
    if not sorted_list_of_pairs:
        return []
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
    raise ValueError(str(el) + " is not in list")


def get_top_count(dict):
    if not dict:
        return None
    return max(dict.items(), key=lambda x: x[1])[0]


def find_closest(value, value_list):
    if not value_list:
        return None, math.inf
    best_diff = math.inf
    best_el = None
    for v in value_list:
        diff = abs(v - value)
        if diff < best_diff:
            best_diff = diff
            best_el = v
    return best_el, best_diff


def rreplace(s, old, new):
    return new.join(s.rsplit(old, 1))


# dict -> counts dict
def get_best_from_count_dicts(dict):
    res = {}
    for k in dict.keys():
        res[k] = get_top_count(dict[k])
    return res


def get_collective_property(value_collection, property_dict):
    counts = defaultdict(int)
    for v in value_collection:
        if v not in property_dict:
            continue
        counts[property_dict[v]] += 1
    return get_top_count(counts)


def proper_plural_form(name, count):
    return str(count) + " " + name + ("" if count == 1 else "s")


# check whether genes overlap and should be processed together
def genes_overlap(gene_db1, gene_db2):
    if gene_db1.seqid != gene_db2.seqid:
        return False
    return overlaps((gene_db1.start, gene_db1.end), (gene_db2.start, gene_db2.end))


def genes_contain(gene_db1, gene_db2):
    if gene_db1.seqid != gene_db2.seqid:
        return False
    return contains((gene_db1.start, gene_db1.end), (gene_db2.start, gene_db2.end))


def argmin(l):
    if len(l) == 0:
        return -1
    min_v = l[0]
    min_i = 0
    for i in range(len(l)):
        if l[i] < min_v:
            min_v = l[i]
            min_i = i
    return min_i


def cmp(x,y):
    if x<y:
        return -1
    elif x>y:
        return 1
    return 0


# == range operations ==
def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def overlaps_at_least(range1, range2, delta=0):
    cutoff = min([delta, range1[1] - range1[0] + 1, range2[1] - range2[0] + 1])
    overlap = min(range1[1], range2[1]) + 1 - max(range1[0], range2[0])
    return overlap >= cutoff


def left_of(range1, range2):
    return range1[1] < range2[0]


def equal_ranges(range1, range2, delta=0):
    return abs(range1[0] - range2[0]) <= delta and abs(range1[1] - range2[1]) <= delta


def covers_end(bigger_range, smaller_range):
    return bigger_range[0] <= smaller_range[0] <= bigger_range[1] <= smaller_range[1]


def covers_start(bigger_range, smaller_range):
    return smaller_range[0] <= bigger_range[0] <= smaller_range[1] <= bigger_range[1]


def contains(bigger_range, smaller_range):
    return bigger_range[1] >= smaller_range[1] and bigger_range[0] <= smaller_range[0]


def contains_well_inside(bigger_range, smaller_range, delta=1):
    return bigger_range[1] >= smaller_range[1]+delta and bigger_range[0] <= smaller_range[0]-delta


def contains_approx(bigger_range, smaller_range, delta = 1):
    return bigger_range[1] + delta >= smaller_range[1] and bigger_range[0] - delta <= smaller_range[0]


def range_list_to_str(range_list, element_delim=',', coord_delim='-'):
    return element_delim.join(list(map(lambda x: str(x[0]) + coord_delim + str(x[1]), range_list)))


def max_range(range1, range2):
    return (min(range1[0], range2[0]), max(range1[1], range2[1]))


def interval_len(interval):
    return interval[1] - interval[0] + 1


def intervals_total_length(sorted_range_list):
    total_len = 0
    for r in sorted_range_list:
        total_len += interval_len(r)
    return total_len


# check that all exons are sorted and have correct coordinates
def validate_exons(self, exons):
    return exons == sorted(exons) and all(0 < x[0] <= x[1] for x in exons)


# sum up all inervals until pos, pos not included
def sum_intervals_to_point(sorted_range_list, pos):
    if pos <= sorted_range_list[0][0]:
        return 0
    elif pos > sorted_range_list[-1][1]:
        return intervals_total_length(sorted_range_list)

    i = 0
    total_len = 0
    while i < len(sorted_range_list) and sorted_range_list[i][0] < pos:
        if sorted_range_list[i][0] <= pos <= sorted_range_list[i][1]:
            total_len += pos - sorted_range_list[i][0]
        else:
            total_len += sorted_range_list[i][1] - sorted_range_list[i][0] + 1
        i += 1
    return total_len


# sum up all intervals from pos, pos not included
def sum_intervals_from_point(sorted_range_list, pos):
    if pos < sorted_range_list[0][0]:
        return intervals_total_length(sorted_range_list)
    elif pos > sorted_range_list[-1][1]:
        return 0

    i = len(sorted_range_list) - 1
    total_len = 0
    while i >= 0 and sorted_range_list[i][1] > pos:
        if sorted_range_list[i][0] <= pos <= sorted_range_list[i][1]:
            total_len += sorted_range_list[i][1] - pos
        else:
            total_len += sorted_range_list[i][1] - sorted_range_list[i][0] + 1
        i -= 1
    return total_len


def jaccard_similarity(sorted_range_list1, sorted_range_list2):
    union = 0
    intersection = 0
    pos1 = 0
    pos2 = 0
    included1 = [0 for i in range(len(sorted_range_list1))]
    included2 = [0 for i in range(len(sorted_range_list2))]

    while pos1 < len(sorted_range_list1) and pos2 < len(sorted_range_list2):
        block1 = sorted_range_list1[pos1]
        block2 = sorted_range_list2[pos2]
        if overlaps(block1, block2):
            assert (included2[pos2] == 0 or included1[pos1] == 0)

            intersection += min(block1[1], block2[1]) - max(block1[0], block2[0]) + 1
            if included2[pos2] == 0 and included1[pos1] == 0:
                # both blocks were not counted
                union += max(block1[1], block2[1]) - min(block1[0], block2[0]) + 1
            elif included2[pos2] == 1:
                # second block was included, take extra bite from block1
                union += max(0, block1[1] - block2[1])
            else:
                # first block was included, take extra bite from block2
                union += max(0, block2[1] - block1[1])

            included1[pos1] = 1
            included2[pos2] = 1
            if block2[1] < block1[1]:
                pos2 += 1
            else:
                pos1 += 1
        elif left_of(block2, block1):
            if included2[pos2] == 0:
                union += block2[1] - block2[0] + 1
                included2[pos2] = 1
            pos2 += 1
        else:
            if included1[pos1] == 0:
                union += block1[1] - block1[0] + 1
                included1[pos1] = 1
            pos1 += 1

    while pos1 < len(sorted_range_list1):
        if included1[pos1] == 0:
            union += sorted_range_list1[pos1][1] - sorted_range_list1[pos1][0] + 1
        pos1 += 1

    while pos2 < len(sorted_range_list2):
        if included2[pos2] == 0:
            union += sorted_range_list2[pos2][1] - sorted_range_list2[pos2][0] + 1
        pos2 += 1

    assert (union != 0)
    return float(intersection) / float(union)


def read_coverage_fraction(read_range_list, isoform_range_list):
    intersection = 0
    pos1 = 0
    pos2 = 0

    while pos1 < len(read_range_list) and pos2 < len(isoform_range_list):
        block1 = read_range_list[pos1]
        block2 = isoform_range_list[pos2]
        if overlaps(block1, block2):
            intersection += min(block1[1], block2[1]) - max(block1[0], block2[0]) + 1
            if block2[1] < block1[1]:
                pos2 += 1
            else:
                pos1 += 1
        elif left_of(block2, block1):
            pos2 += 1
        else:
            pos1 += 1

    read_length = intervals_total_length(read_range_list)
    return float(intersection) / float(read_length)


def extra_exon_percentage(isoform_region, read_exons):
    total_read_len = 0
    outside_read_len = 0
    for e in read_exons:
        total_read_len += e[1] - e[0] + 1
        if e[0] < isoform_region[0]:
            outside_read_len += min(e[1], isoform_region[0] - 1) - e[0] + 1
        if e[1] > isoform_region[1]:
            outside_read_len += e[1] - max(e[0], isoform_region[1] + 1) + 1

    return float(outside_read_len) / float(total_read_len)


# == working with alignment blocks ==
def junctions_from_blocks(sorted_blocks):
    junctions = []
    if len(sorted_blocks) >= 2:
        for i in range(0, len(sorted_blocks) - 1):
            if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
                junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


def get_following_exon_from_junctions(region, introns, intron_position):
    if intron_position == len(introns) - 1 or intron_position == -1:
        # intron precedes the last exon
        following_exon_end = region[1]
    else:
        following_exon_end = introns[intron_position + 1][0] - 1
    return (introns[intron_position][1] + 1, following_exon_end)


def get_preceding_exon_from_junctions(region, introns, intron_position):
    assert intron_position <= len(introns)
    if intron_position == 0:
        # intron precedes the first exon
        preceding_exon_start = region[0]
    else:
        preceding_exon_start = introns[intron_position - 1][1] + 1

    if intron_position == len(introns):
        return (preceding_exon_start, region[1])
    return (preceding_exon_start, introns[intron_position][0] - 1)


def get_exon(read_region, read_junctions, exon_position):
    assert exon_position <= len(read_junctions)
    if exon_position < 0:
        exon_position = len(read_junctions) + exon_position + 1

    if exon_position == 0:
        return (read_region[0], read_junctions[0][0] - 1)
    elif exon_position == len(read_junctions):
        return (read_junctions[-1][1] + 1, read_region[1])
    else:
        return (read_junctions[exon_position - 1][1] + 1, read_junctions[exon_position][0] - 1)


def get_exons(read_region, read_introns):
    return junctions_from_blocks([(-math.inf, read_region[0] - 1)] + read_introns + [(read_region[1] + 1, math.inf)])


def concat_gapless_blocks(blocks, cigar_tuples):
    cigar_index = 0
    block_index = 0
    resulting_blocks = []

    current_block = None
    deletions_before_block = 0

    while cigar_index < len(cigar_tuples) and block_index < len(blocks):
        # init new block
        cigar_event = cigar_tuples[cigar_index][0]
        if current_block is None:
            # init new block from match
            if cigar_event in [0, 7, 8]:
                current_block = (blocks[block_index][0] - deletions_before_block, blocks[block_index][1])
                deletions_before_block = 0
                block_index += 1
            # keep track of deletions before matched block
            elif cigar_event == 2:
                deletions_before_block = cigar_tuples[cigar_index][1]
        # found intron, add current block
        elif cigar_event == 3:
            resulting_blocks.append(current_block)
            current_block = None
        # add deletion to block
        elif cigar_event == 2:
            current_block = (current_block[0], current_block[1] + cigar_tuples[cigar_index][1])
        # found match - merge blocks
        elif cigar_event in [0, 7, 8]:
            # if abs(current_block[1] - blocks[block_index][0]) > 1:
            #    logger.debug("Distant blocks")
            #    logger.debug(current_block, blocks[block_index])
            current_block = (current_block[0], blocks[block_index][1])

            block_index += 1
        cigar_index += 1

    if current_block is not None:
        resulting_blocks.append(current_block)

    return resulting_blocks


def get_read_blocks(ref_start, cigar_tuples):
    read_pos = 0
    ref_pos = ref_start + 1
    cigar_index = 0
    current_ref_block_start = None
    current_read_block_start = None
    current_cigar_block_start = None
    has_match=False
    ref_blocks = []
    cigar_blocks = []
    read_blocks = []

    while cigar_index < len(cigar_tuples):
        cigar_event = cigar_tuples[cigar_index][0]
        event_len = cigar_tuples[cigar_index][1]

        if current_ref_block_start is None and cigar_event in [0, 1, 2, 7, 8]:
            # init new block from match
            current_ref_block_start = ref_pos
            current_read_block_start = read_pos
            current_cigar_block_start = cigar_index
            if cigar_event == 1:
                read_pos += event_len
            elif cigar_event == 2:
                ref_pos += event_len
            else:
                read_pos += event_len
                ref_pos += event_len
                has_match = True
        # found intron, add current block
        elif cigar_event in [0, 7, 8]:
            read_pos += event_len
            ref_pos += event_len
            has_match = True
        elif cigar_event == 1:
            read_pos += event_len
        elif cigar_event == 2:
            ref_pos += event_len
        elif cigar_event == 3:
            if current_ref_block_start:
                if has_match:
                    ref_blocks.append((current_ref_block_start, ref_pos - 1))
                    read_blocks.append((current_read_block_start, read_pos - 1))
                    cigar_blocks.append((current_cigar_block_start, cigar_index - 1))
                has_match = False
                current_ref_block_start = None
                current_read_block_start = None
                current_cigar_block_start = None
            ref_pos += event_len
        elif cigar_event == 4:
            if current_ref_block_start:
                if has_match:
                    ref_blocks.append((current_ref_block_start, ref_pos - 1))
                    read_blocks.append((current_read_block_start, read_pos - 1))
                    cigar_blocks.append((current_cigar_block_start, cigar_index - 1))
                has_match = False
                current_ref_block_start = None
                current_read_block_start = None
                current_cigar_block_start = None
            read_pos += event_len

        cigar_index += 1

    if current_ref_block_start and has_match:
        ref_blocks.append((current_ref_block_start, ref_pos - 1))
        read_blocks.append((current_read_block_start, read_pos - 1))
        cigar_blocks.append((current_cigar_block_start, cigar_index - 1))

    return ref_blocks, read_blocks, cigar_blocks


def correct_bam_coords(blocks):
    return list(map(lambda x: (x[0] + 1, x[1]), blocks))


def truncate_read_to_polya(read_exons, polya_pos, polyt_pos):
    end_index = len(read_exons) - 1
    end_position = read_exons[-1][1]
    if polya_pos != -1:
        end_position = polya_pos
        while end_index >= 0:
            if read_exons[end_index][0] < polya_pos:
                break
            end_index -= 1

    start_index = 0
    start_position = read_exons[0][0]
    if polyt_pos != -1:
        start_position = polyt_pos
        while start_index <= end_index:
            if read_exons[start_index][1] > polyt_pos:
                break
            start_index += 1

    if start_position == read_exons[0][0] and end_position == read_exons[-1][1]:
        return read_exons

    if start_index == end_index:
        return [(start_position, end_position)]

    new_exons = [(start_position, read_exons[start_index][1])]
    new_exons += read_exons[start_index+1:end_index]
    new_exons.append((read_exons[end_index][0], end_position))
    return new_exons


# == profile functions ==
def is_subprofile(short_isoform_profile, long_isoform_profile):
    assert len(short_isoform_profile) == len(long_isoform_profile)
    assert 0 not in short_isoform_profile

    short_range_start = None
    short_range_end = None
    for i in range(len(short_isoform_profile)):
        if short_isoform_profile[i] in {-1, 1}:
            if short_range_start is None:
                short_range_start = i
            short_range_end = i

    assert short_range_end is not None

    for i in range(short_range_start, short_range_end + 1):
        if short_isoform_profile[i] != long_isoform_profile[i]:
            return False
    return True


def difference_in_present_features(profile1, profile2, diff_limit=-1):
    """ computes Hamming distance for two profiles

    Parameters
    ----------
    profile1: list(int)
    profile2: list(int)

    Returns
    -------
    result: int
        number of different items in lists

    """
    assert len(profile1) == len(profile2)
    d = 0
    if diff_limit == -1:
        diff_limit = len(profile1) + 1
    for i in range(len(profile1)):
        if profile2[i] == 0 or profile1[i] == 0:
            continue
        if profile1[i] != profile2[i]:
            d += 1
        if d > diff_limit:
            break
    return d


def count_both_present_features(profile1, profile2):
    assert len(profile1) == len(profile2)
    d = 0
    for i in range(len(profile1)):
        if profile1[i] == profile2[i] == 1:
            d += 1
    return d


def all_features_present(isoform_profile, read_profile):
    assert len(isoform_profile) == len(read_profile)
    for i in range(len(isoform_profile)):
        if isoform_profile[i] == 1 and read_profile[i] != 1:
            return False
    return True


def find_matching_positions(profile1, profile2):
    assert len(profile1) == len(profile2)
    matches = [0 for i in range(len(profile1))]
    for i in range(len(profile1)):
        if profile1[i] == profile2[i]:
            matches[i] = 1
    return matches


def has_overlapping_features(profile1, profile2):
    assert len(profile1) == len(profile2)
    for i in range(len(profile1)):
        if profile1[i] == profile2[i] == 1:
            return True
    return False


def has_inconsistent_features(read_profile, gene_profile):
    assert len(read_profile) == len(gene_profile)
    for i in range(len(read_profile)):
        if read_profile[i] != gene_profile[i] and read_profile[i] != 0:
            return True
    return False


def mask_profile(read_profile, true_profile):
    assert len(read_profile) == len(true_profile)

    masked_profile = []
    for i in range(len(true_profile)):
        if true_profile[i] == 1:
            masked_profile.append(read_profile[i])
        else:
            masked_profile.append(0)
    return masked_profile


def get_blocks_from_profile(features, profile):
    assert len(features) == len(profile)
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


def count_noncanonincal(introns, reference_region, strand, ref_region_start=0):
    count = 0
    for intron in introns:
        intron_left_pos = intron[0] - ref_region_start
        intron_right_pos = intron[1] - ref_region_start
        left_site = reference_region[intron_left_pos:intron_left_pos + 2]
        right_site = reference_region[intron_right_pos - 1:intron_right_pos + 1]
        if strand == '+':
            count += 0 if (left_site, right_site) == ("GT", "AG") else 1
        else:
            count += 0 if (left_site, right_site) == ("CT", "AC") else 1
    return count


def get_intron_strand(intron, reference_region, ref_region_start=1):
    intron_left_pos = intron[0] - ref_region_start
    intron_right_pos = intron[1] - ref_region_start
    left_site = str(reference_region[intron_left_pos:intron_left_pos + 2].seq)
    right_site = str(reference_region[intron_right_pos - 1:intron_right_pos + 1].seq)
    is_fwd = (left_site, right_site) == ("GT", "AG")
    is_rev = (left_site, right_site) == ("CT", "AC")

    if is_fwd == is_rev:
        return '.'
    return '+' if is_fwd else '-'


def get_strand(introns, reference_region, ref_region_start=1):
    if len(introns) == 0:
        return '.'

    count_fwd = 0
    count_rev = 0
    for intron in introns:
        intron_left_pos = intron[0] - ref_region_start
        intron_right_pos = intron[1] - ref_region_start
        left_site = str(reference_region[intron_left_pos:intron_left_pos + 2].seq)
        right_site = str(reference_region[intron_right_pos - 1:intron_right_pos + 1].seq)
        count_fwd += 1 if (left_site, right_site) == ("GT", "AG") else 0
        count_rev += 1 if (left_site, right_site) == ("CT", "AC") else 0

    if count_fwd == count_rev:
        return '.'
    return '+' if count_rev < count_fwd else '-'
