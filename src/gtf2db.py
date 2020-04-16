############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import gffutils
import argparse
from traceback import print_exc


def gtf2db(inf, outf):
    gffutils.create_db(inf, outf, force=True, keep_order=True, merge_strategy='merge',
                       sort_attribute_values=True, disable_infer_transcripts=True, disable_infer_genes=True)


def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: " + sys.argv[0] + " <GFF/GTF file> <output db file> \n")
        exit(0)
    gtf2db(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
