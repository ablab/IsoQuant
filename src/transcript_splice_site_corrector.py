import logging
logger = logging.getLogger('IsoQuant')

def extract_location_from_cigar_string(cigartuples: list,
                                    read_start: int,
                                    read_end: int,
                                    splice_site_location: int):
    """
    Extract location from cigar string.

    Args:
        cigar_tuples (list): list of cigar tuples (cigar code, aligned position).
        See pysam documentation for more information
        read_start (int): the start location for the read (base-1)
        read_end (int): the end location for the read (base-1)
        splice_site_location (int): location of interest (base-1)

    Returns:
        _type_: _description_
    """
    relative_position = splice_site_location - read_start
    alignment_position = 0
    ref_position = 0

    for cigar_code in cigartuples:

        if cigar_code[0] in [0, 2, 3, 7, 8]:
            ref_position += cigar_code[1]
        if ref_position <= relative_position and not \
                read_start + ref_position == read_end:
            alignment_position += cigar_code[1]
        else:
            return alignment_position + (cigar_code[1] - (ref_position - relative_position))

    return -1


def count_deletions_from_cigar_codes_in_given_window(cigartuples: list,
                                                aligned_location: int,
                                                location_is_end: bool,
                                                splice_site_data: dict,
                                                window_size: int):
    """
    Get cigar codes in a given window.

    Args:
        cigar_tuples (list): list of cigar tuples (cigar code, aligned position). See 
        pysam documentation for more information
        aligned_location (int): aligned location
        loc_type (str): type of location (start or end)
    """

    count_of_deletions = 0
    

    cigar_code_list = []
    location = 0

    if location_is_end:
        aligned_location = aligned_location - window_size + 1

    for cigar_code in cigartuples:
        if window_size == len(cigar_code_list):
            break
        if location + cigar_code[1] > aligned_location:
            overlap = location + \
                cigar_code[1] - (aligned_location + len(cigar_code_list))
            cigar_code_list.extend(
                [cigar_code[0] for _ in range(min(window_size -
                                                len(cigar_code_list), overlap))])
        location += cigar_code[1]

    for i in range(window_size):
        if i >= len(cigar_code_list):
            break
        if cigar_code_list[i] == 2:
            count_of_deletions += 1
            splice_site_data["del_pos_distr"][i] += 1
    
    if count_of_deletions not in splice_site_data["deletions"]:
        splice_site_data["deletions"][count_of_deletions] = 0
    
    splice_site_data["deletions"][count_of_deletions] += 1


def extract_splice_site_locations_within_aligned_read(read_start: int, read_end: int, exons:list):
    matching_locations = []
    for exon_start, exon_end in exons:
        if read_start <= exon_start <= read_end:
            location_is_end = False
            matching_locations.append((exon_start, location_is_end))
        if read_start <= exon_end <= read_end:
            location_is_end = True
            matching_locations.append((exon_end, location_is_end))
        if read_end <= exon_end:
            break
    return matching_locations


def count_deletions_for_splice_site_locations(
        read_start: int, 
        read_end: int, 
        cigartuples: list, 
        exons: list, 
        splice_site_cases: dict):
    """

    Args:
        assigned_read (ReadAssignment): read assignment
        exons (list): tuple of exons (start, end)
        splice_site_cases (dict): a dictionary for storing splice site cases
    """

    
    # Constant window size for counting deletions
    WINDOW_SIZE = 8
    
    # Extract splice site locations within aligned read
    matching_locations = extract_splice_site_locations_within_aligned_read(read_start, read_end, exons)
    
    logger.debug(f"Matching locations: {matching_locations}")
    # Count deletions for each splice site location
    for splice_site_location, location_type in matching_locations:
        if splice_site_location not in splice_site_cases:
            splice_site_cases[splice_site_location] = {
                'location_is_end': location_type,  
                'deletions': {},
                'del_pos_distr': [0 for _ in range(WINDOW_SIZE)],
                'most_common_del': -1,
                'del_location_has_canonical_nucleotides': False
            }
        
        # Processing cigartuples
        # 1. Find the aligned location
        aligned_location = extract_location_from_cigar_string(cigartuples, read_start, read_end, splice_site_location)
        # 2. Count deletions in a predefined window
        count_deletions_from_cigar_codes_in_given_window(
            cigartuples, 
            aligned_location, 
            location_type, 
            splice_site_cases[splice_site_location], 
            WINDOW_SIZE)



def compute_most_common_case_of_deletions(deletions: dict, location_is_end: bool):
    del_most_common_case = [k for k, v in deletions.items(
    ) if v == max(deletions.values())]
    if len(del_most_common_case) == 1:
        if location_is_end:
            return -del_most_common_case[0]
        return del_most_common_case[0]
    return -1    


def extract_nucleotides_from_most_common_del_location(
        location: int, 
        splice_site_data: dict, 
        chr_record, 
        strand: str):
    most_common_del = splice_site_data["most_common_del"]
    idx_correction = -1
    extraction_start = location + most_common_del + idx_correction
    extraction_end = location + most_common_del + 2 + idx_correction
    try:
        extracted_canonicals =  chr_record[extraction_start:extraction_end]
    except KeyError:
        extracted_canonicals = 'XX'
    
    
    canonical_pairs = {
        '+': {
            'start': ['AG', 'AC'],
            'end': ['GT', 'GC', 'AT']
        },
        '-': {
            'start': ['AC', 'GC', 'AC'],
            'end': ['CT', 'GT']
        }
    }
    if splice_site_data["location_is_end"]:
        possible_canonicals = canonical_pairs[strand]['end']
    else:
        possible_canonicals = canonical_pairs[strand]['start']
    if extracted_canonicals in possible_canonicals:
        splice_site_data["del_location_has_canonical_nucleotides"] = True

def compute_most_common_del_and_verify_nucleotides(
        splice_site_location: int, 
        splice_site_data: dict, 
        chr_record,
        ACCEPTED_DEL_CASES: list,
        strand: str,):
    
    
    # Compute most common case of deletions
    splice_site_data["most_common_del"] = compute_most_common_case_of_deletions(
        splice_site_data["deletions"],
        splice_site_data["location_is_end"])
    
    # Extract nucleotides from most common deletion location if it is an accepted case
    if abs(splice_site_data["most_common_del"]) in ACCEPTED_DEL_CASES:
        extract_nucleotides_from_most_common_del_location(
            splice_site_location, 
            splice_site_data, 
            chr_record, 
            strand)



def threshold_for_del_cases_exceeded(
        del_pos_distr: list,
        deletions: dict,
        most_common_del: int,
        THRESHOLD_CASES_AT_LOCATION):
    total_cases = sum(deletions.values())
    nucleotides_exceeding_treshold = 0
    for value in del_pos_distr:
        if value / total_cases > THRESHOLD_CASES_AT_LOCATION:
            nucleotides_exceeding_treshold += 1
    return bool(nucleotides_exceeding_treshold >= abs(most_common_del))

def sublist_largest_values_exists(lst, n):
    """
        Verifies that there is a sublist of size n that contains the largest values in the list.
        Not currently in use, but may be included in the error prediction strategy for stricter prediction. 
    Args:
        lst (int): list of deletion distribution
        n (int): most common case of deletions

    Returns:
        _type_: _description_
    """
    largest_values = set(sorted(lst, reverse=True)[:n])
    count = 0

    for num in lst:
        if num in largest_values:
            count += 1
            if count >= n:
                return True
        else:
            count = 0

    return False


def correct_splice_site_errors(
        splice_site_cases: dict,
        MIN_N_OF_ALIGNED_READS: int,
        ACCEPTED_DEL_CASES: list,
        THRESHOLD_CASES_AT_LOCATION: float,
        MORE_CONSERVATIVE_STRATEGY: bool,
        strand: str,
        chr_record):
    """ 1. Count most common deletion at each splice site location
        2. For interesting cases count nucleotides at deletion positions
        3. If canonical nucleotides are found, correct splice site

    Args:
        splice_site_cases (dict): collected splice site cases
        MIN_N_OF_ALIGNED_READS (int): constant for minimum number of aligned reads
        ACCEPTED_DEL_CASES (list): constant for accepted cases of deletions
        MORE_CONSERVATIVE_STRATEGY (bool): constant for more conservative strategy
        strand (str): transcript strand (extracted from first ReadAssignment-object in read_assignments list)
        chr_record (Fasta): FASTA recored, i.e. a single chromosome from a reference
    """

    
    
    locations_with_errors = []
    for splice_site_location, splice_site_data in splice_site_cases.items():
        
        reads = sum(splice_site_data["deletions"].values())
        if reads < MIN_N_OF_ALIGNED_READS:
            continue
        
        compute_most_common_del_and_verify_nucleotides(
            splice_site_location, 
            splice_site_data, 
            chr_record,
            ACCEPTED_DEL_CASES,
            strand
            )
        if MORE_CONSERVATIVE_STRATEGY:
            if not sublist_largest_values_exists(
                splice_site_data["del_pos_distr"],
                abs(splice_site_data["most_common_del"])):
                continue
            if not threshold_for_del_cases_exceeded(
                splice_site_data["del_pos_distr"],
                splice_site_data["deletions"],
                splice_site_data["most_common_del"],
                THRESHOLD_CASES_AT_LOCATION):
                continue

        if splice_site_data["del_location_has_canonical_nucleotides"]:
            locations_with_errors.append(splice_site_location)
    
    return locations_with_errors

def generate_updated_exon_list(
        splice_site_cases: dict,
        locations_with_errors: list,
        exons: list):
    updated_exons = []
    for exon in exons:
            updated_exon = exon
            if exon[0] in locations_with_errors:
                corrected_location = exon[0] + splice_site_cases[exon[0]]["most_common_del"]
                updated_exon = (corrected_location, exon[1])
            if exon[1] in locations_with_errors:
                corrected_location = exon[1] + splice_site_cases[exon[1]]["most_common_del"]
                updated_exon = (exon[0], corrected_location)
            updated_exons.append(updated_exon)
    return updated_exons