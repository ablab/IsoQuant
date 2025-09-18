############################################################################
# Copyright (c) 2025 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os


def convert_chr_id_to_file_name_str(chr_id: str):
    return chr_id.replace('/', '_')


def reads_collected_lock_file_name(sample_out_raw, chr_id):
    return "{}_{}_collected".format(sample_out_raw, convert_chr_id_to_file_name_str(chr_id))


def reads_processed_lock_file_name(dump_filename, chr_id):
    chr_dump_file = dump_filename + "_" + convert_chr_id_to_file_name_str(chr_id)
    return "{}_processed".format(chr_dump_file)


def read_group_lock_filename(sample):
    return sample.read_group_file + "_lock"


def split_barcodes_lock_filename(sample):
    return sample.barcodes_split_reads + "_lock"


def clean_locks(chr_ids, base_name, fname_function):
    for chr_id in chr_ids:
        fname = fname_function(base_name, chr_id)
        if os.path.exists(fname):
            os.remove(fname)


def saves_file_name(out_raw_file: str, chr_id: str):
    return out_raw_file + "_" + convert_chr_id_to_file_name_str(chr_id)


def multimappers_file_name(out_raw_file: str, chr_id: str):
    return out_raw_file + "_multimappers_" + convert_chr_id_to_file_name_str(chr_id)


def filtered_reads_file_name(out_raw_file: str, chr_id: str):
    return out_raw_file + "_filtered_" + chr_id