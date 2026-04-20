############################################################################
# Copyright (c) 2024-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
from enum import IntEnum


class IsoQuantExitCode(IntEnum):
    """Exit codes for IsoQuant."""

    # Success
    SUCCESS = 0

    # Input/Output Errors (1-19)
    INPUT_FILE_NOT_FOUND = 1
    OUTPUT_ALREADY_EXISTS = 2
    BAM_NOT_INDEXED = 3
    GENE_DB_NOT_FOUND = 4
    REFERENCE_NOT_FOUND = 5
    NO_INPUT_DATA = 6
    FILE_FORMAT_ERROR = 7
    YAML_PARSING_ERROR = 8
    INVALID_FILE_FORMAT = 9

    # Configuration Errors (20-39)
    INVALID_PARAMETER = 20
    MISSING_REQUIRED_OPTION = 21
    INCOMPATIBLE_OPTIONS = 22
    RESUME_CONFIG_NOT_FOUND = 23
    CORRUPTED_GTF = 24
    INVALID_CONVERSION_MODE = 25

    # Data Validation Errors (40-59)
    LABEL_COUNT_MISMATCH = 40
    DUPLICATE_FILES = 41
    CHROMOSOME_MISMATCH = 42
    REPLICA_COUNT_MISMATCH = 43

    # Alignment Errors (60-79)
    ALIGNER_NOT_FOUND = 60
    INDEXING_FAILED = 61
    ALIGNMENT_FAILED = 62
    SAMTOOLS_FAILED = 63
    SUBPROCESS_FAILED = 64

    # Barcode Calling Errors (80-89)
    BARCODE_FILE_COUNT_MISMATCH = 80
    BARCODE_WHITELIST_MISSING = 81

    # Test/Runtime Errors (90-99)
    TEST_FAILED = 90
    UNCAUGHT_EXCEPTION = 99

    # Standard Error (for utilities)
    STANDARD_ERROR = 1


def exit_with_code(code: IsoQuantExitCode, message: str = None):
    """
    Exit with a specific error code and optional message.

    Parameters
    ----------
    code : IsoQuantExitCode
        The exit code to use
    message : str, optional
        Additional error message to log
    """
    if message:
        import logging
        logger = logging.getLogger('IsoQuant')
        if code == IsoQuantExitCode.SUCCESS:
            logger.info(message)
        else:
            logger.error(message)
    sys.exit(code)
