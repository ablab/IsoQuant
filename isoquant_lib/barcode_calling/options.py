############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Resolving and validating the barcode calling command line options.

Shared by isoquant.py and isoquant_detect_barcodes.py so that both entry points agree on
what --mode, --split_molecules, --n_cells and --barcode_correction mean.
"""

import gzip
import logging
import sys

import pysam

from ..modes import (AUTO_BARCODES, BarcodeCorrectionMethod, DEPRECATED_MODE_ALIASES,
                     IsoQuantMode, LARGE_WHITELIST_SIZE, SPLIT_MOLECULES_AUTO,
                     SPLIT_MOLECULES_FALSE, SPLIT_MOLECULES_TRUE)
from ..utils.error_codes import IsoQuantExitCode
from ..utils.file_utils import check_file_exists
from .detect_barcodes import get_umi_length

logger = logging.getLogger('IsoQuant')


def count_whitelist_barcodes(whitelist_files):
    total = 0
    for file_name in whitelist_files:
        handle = gzip.open(file_name, "rt") if file_name.endswith(("gz", "gzip")) else open(file_name)
        with handle:
            total += sum(1 for _ in handle)
    return total


def resolve_n_cells(args):
    """Normalise --n_cells to an int, the string AUTO_BARCODES, or None."""
    if args.n_cells is None:
        return None
    if args.n_cells == AUTO_BARCODES:
        return AUTO_BARCODES
    try:
        n_cells = int(args.n_cells)
    except ValueError:
        logger.critical('--n_cells must be a positive integer or "%s"' % AUTO_BARCODES)
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
    if n_cells <= 0:
        logger.critical("--n_cells must be positive")
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
    return n_cells


def barcode_whitelist_of(args):
    """The whitelist, under whichever option name this entry point uses for it."""
    if hasattr(args, "barcode_whitelist"):
        return args.barcode_whitelist
    return args.barcodes


def resolve_barcode_correction(args, whitelist_option="--barcode_whitelist"):
    """Decide whether cell barcodes are supplied or detected from the data.

    --n_cells decides what the whitelist means: unset it is the cell barcodes themselves,
    set it is a pool to select them from, which costs an extra pass over the reads.
    Sets args.detect_cell_barcodes and normalises args.n_cells.
    """
    args.detect_cell_barcodes = False
    args.n_cells = resolve_n_cells(args)
    whitelist = barcode_whitelist_of(args)
    whitelist_is_auto = whitelist == [AUTO_BARCODES]

    # nothing to detect when barcodes come ready-made
    ready_made = getattr(args, "barcoded_reads", None) or getattr(args, "barcoded_bam", False)
    if ready_made:
        if args.n_cells is not None:
            logger.warning("--n_cells is ignored: barcodes are taken from %s, so there is nothing "
                           "to detect" % ("--barcoded_bam" if args.barcoded_bam else "--barcoded_reads"))
        return

    requested = BarcodeCorrectionMethod[args.barcode_correction]
    if requested == BarcodeCorrectionMethod.whitelist:
        if whitelist_is_auto:
            logger.critical("--barcode_correction whitelist cannot be used with %s %s"
                            % (whitelist_option, AUTO_BARCODES))
            sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)
        if args.n_cells is not None:
            logger.warning("--n_cells is ignored: --barcode_correction %s matches reads against the "
                           "whitelist as given" % BarcodeCorrectionMethod.whitelist.name)
        return

    # a detected whitelist needs a cell count; without one, estimate it
    if whitelist_is_auto and args.n_cells is None:
        args.n_cells = AUTO_BARCODES

    detect = whitelist_is_auto or args.n_cells is not None or requested == BarcodeCorrectionMethod.detect
    if not detect:
        if whitelist:
            warn_on_large_whitelist(args, whitelist_option)
        return

    if not args.mode.supports_cell_barcode_detection():
        logger.critical("Detecting cell barcodes from the data is not supported for mode %s" % args.mode.name)
        sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)

    if requested == BarcodeCorrectionMethod.detect and args.n_cells is None:
        args.n_cells = AUTO_BARCODES
    args.detect_cell_barcodes = True


def warn_on_large_whitelist(args, whitelist_option="--barcode_whitelist"):
    """A big whitelist taken as the cell list makes per-read matching demand exact matches."""
    barcode_count = count_whitelist_barcodes(barcode_whitelist_of(args))
    if barcode_count <= LARGE_WHITELIST_SIZE:
        return
    logger.warning("%s contains %d barcodes and is treated as the list of cell barcodes. "
                   "Matching every read against a list this large effectively requires an exact match, "
                   "so reads carrying a sequencing error in the barcode will be lost."
                   % (whitelist_option, barcode_count))
    logger.warning('Set --n_cells (or --n_cells %s) to select cell barcodes from the whitelist instead.'
                   % AUTO_BARCODES)


def resolve_deprecated_mode(args):
    """Translate a superseded mode name into a chemistry plus a --split_molecules default."""
    alias = DEPRECATED_MODE_ALIASES.get(args.mode)
    if not alias:
        return
    mode_name, splits = alias
    replacement = "--mode %s --split_molecules %s" % (
        mode_name, SPLIT_MOLECULES_TRUE if splits else SPLIT_MOLECULES_FALSE)
    logger.warning("Mode %s is deprecated, use `%s` instead" % (args.mode, replacement))
    args.mode = mode_name
    # only a default: an explicit --split_molecules wins
    if args.split_molecules is None:
        args.split_molecules = SPLIT_MOLECULES_TRUE if splits else SPLIT_MOLECULES_FALSE


def resolve_split_molecules(args):
    """Turn --split_molecules into a bool, refusing to silently ignore an impossible request."""
    if isinstance(args.split_molecules, bool):
        # already resolved, e.g. args restored from a previous run by --resume
        return
    requested = args.split_molecules or SPLIT_MOLECULES_AUTO
    supported = args.mode.supports_molecule_splitting()

    if requested == SPLIT_MOLECULES_FALSE:
        args.split_molecules = False
        return
    if requested == SPLIT_MOLECULES_TRUE and not supported:
        logger.critical("Mode %s cannot split reads into separate molecules. Drop "
                        "--split_molecules %s, or use --split_molecules %s to let IsoQuant "
                        "decide per protocol." % (args.mode.name, SPLIT_MOLECULES_TRUE, SPLIT_MOLECULES_AUTO))
        sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)

    args.split_molecules = supported
    if args.split_molecules:
        logger.info("Reads will be split into separate cDNA molecules")


def validate_barcode_whitelist(args):
    """The whitelist is either a set of files or the literal AUTO_BARCODES, never both."""
    whitelist = barcode_whitelist_of(args)
    if not whitelist or AUTO_BARCODES not in whitelist:
        return
    if len(whitelist) > 1:
        logger.critical('--barcode_whitelist %s cannot be combined with whitelist files' % AUTO_BARCODES)
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
    if not args.mode.supports_cell_barcode_detection():
        logger.critical('--barcode_whitelist %s is not supported for mode %s, please provide a whitelist file'
                        % (AUTO_BARCODES, args.mode.name))
        sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)


def validate_barcode_calling(args):
    args.umi_length = 0
    args.detect_cell_barcodes = False
    if not args.mode.needs_barcode_calling():
        return
    validate_barcode_whitelist(args)
    barcode_sources = sum([bool(args.barcode_whitelist), bool(args.barcoded_reads), bool(args.barcoded_bam)])
    if barcode_sources > 1:
        logger.critical("Options --barcode_whitelist, --barcoded_reads, and --barcoded_bam are mutually exclusive")
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
    if args.mode == IsoQuantMode.custom_sc:
        if not any([args.molecule, args.barcoded_reads, args.barcoded_bam]):
            logger.critical("custom_sc mode requires --molecule, --barcoded_reads, or --barcoded_bam")
            sys.exit(IsoQuantExitCode.BARCODE_WHITELIST_MISSING)
    elif not any([args.barcode_whitelist, args.barcoded_reads, args.barcoded_bam]):
        logger.critical("You have chosen single-cell/spatial mode %s, please specify barcode whitelist, "
                        "file with barcoded reads, or --barcoded_bam" % args.mode.name)
        sys.exit(IsoQuantExitCode.BARCODE_WHITELIST_MISSING)
    if args.barcoded_bam:
        args.umi_length = detect_umi_length_from_bam(args.input_data.samples[0].file_list[0][0], args.umi_tag)
    else:
        args.umi_length = get_umi_length(args.mode)

    resolve_barcode_correction(args)


def detect_umi_length_from_bam(bam_path: str, umi_tag: str) -> int:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.has_tag(umi_tag):
                return len(read.get_tag(umi_tag))
    return 0


def check_barcode_input_files(args):
    # Check barcoded reads files (from args, not sample - sample.barcoded_reads is set later)
    if hasattr(args, 'barcoded_reads') and args.barcoded_reads:
        if isinstance(args.barcoded_reads, list):
            for bc_file in args.barcoded_reads:
                check_file_exists(bc_file, "Barcoded reads file")
        else:
            check_file_exists(args.barcoded_reads, "Barcoded reads file")

    # Check molecule definition file
    if hasattr(args, 'molecule') and args.molecule:
        check_file_exists(args.molecule, "Molecule definition file")

    # Check barcode whitelist files; AUTO_BARCODES asks for detection, it is not a path
    if hasattr(args, 'barcode_whitelist') and args.barcode_whitelist:
        for wl_file in args.barcode_whitelist:
            if wl_file == AUTO_BARCODES:
                continue
            check_file_exists(wl_file, "Barcode whitelist file")


def check_barcode_mapping_files(args):
    from isoquant_lib.assignment.read_groups import parse_barcode2spot_spec

    # Check barcode2spot file (parse spec to extract filename)
    if hasattr(args, 'barcode2spot') and args.barcode2spot:
        bc2spot_file, _, _ = parse_barcode2spot_spec(args.barcode2spot)
        check_file_exists(bc2spot_file, "Barcode to spot mapping file")

    # Check barcode2barcode file (parse spec to extract filename)
    if hasattr(args, 'barcode2barcode') and args.barcode2barcode:
        bc2bc_file, _, _ = parse_barcode2spot_spec(args.barcode2barcode)
        check_file_exists(bc2bc_file, "Barcode to barcode mapping file")
