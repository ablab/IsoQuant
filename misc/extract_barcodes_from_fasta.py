############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
from Bio import SeqIO

def process_fasta(input_fasta, out_fasta, out_tsv):
    with open(out_fasta, "w") as output_handle, open(out_tsv, "w") as output_tsv_handle:
        new_contigs = []
        for seq_record in SeqIO.parse(input_fasta, "fasta"):
            delim_pos = seq_record.id.find('_barcodeIDs_')
            if delim_pos == -1:
                print("Malformed contig id: " + seq_record.id)
            contig_id = seq_record.id[:delim_pos]
            suffix_pos = seq_record.id.find('_f_')
            if delim_pos == -1:
                barcodes = seq_record.id[delim_pos+len('_barcodeIDs_'):]
            else:
                barcodes = seq_record.id[delim_pos + len('_barcodeIDs_'):suffix_pos]
            seq_record.id = contig_id
            seq_record.description = ""
            new_contigs.append(seq_record)

            output_tsv_handle.write(contig_id + "\t" + barcodes + "\n")

        SeqIO.write(new_contigs, output_handle, "fasta")


def main(args):
    if len(args) != 2:
        print("Usage: " + sys.argv[0] + " <rnaSPAdes barcoded contigs>")
        sys.exit(0)

    input_fasta = sys.argv[1]
    base, ext = os.path.splitext(input_fasta)
    out_fasta = base + ".short_names.fasta"
    out_tsv = base + ".barcodes.tsv"
    process_fasta(input_fasta, out_fasta, out_tsv)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main(sys.argv)
    except SystemExit:
        raise
    except:
        sys.exit(-1)