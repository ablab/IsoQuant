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


def umi_filtered_reads_file_name(out_umi_filtered_tmp: str, chr_id:str, edit_distance:int):
    return out_umi_filtered_tmp + ("_%s_ED%d" % (chr_id, edit_distance))


def umi_filtered_lock_file_name(out_umi_filtered_done: str, chr_id: str, edit_distance: int):
    return out_umi_filtered_done + ("_%s_ED%d" % (chr_id, edit_distance))


def allinfo_file_name(out_umi_filtered_tmp: str, chr_id: str, edit_distance: int):
    return umi_filtered_reads_file_name(out_umi_filtered_tmp, chr_id, edit_distance) + ".allinfo"


def allinfo_stats_file_name(out_umi_filtered_tmp: str, chr_id: str, edit_distance: int):
    return umi_filtered_reads_file_name(out_umi_filtered_tmp, chr_id, edit_distance) + ".stats.tsv"


def umi_filtered_global_lock_file_name(out_umi_filtered_done: str):
    return out_umi_filtered_done + ".lock"


def dynamic_pools_file_name(out_raw_file: str, chr_id: str):
    """File name for per-chromosome dynamic string pools (for read groups from BAM tags/read IDs)."""
    return out_raw_file + "_dynamic_pools_" + convert_chr_id_to_file_name_str(chr_id)


# Read collection auxiliary files

def read_groups_file_name(save_file: str):
    """File name for storing read group assignments during collection."""
    return save_file + "_groups"


def bamstat_file_name(save_file: str):
    """File name for BAM alignment statistics."""
    return save_file + "_bamstat"


def info_file_name(out_raw_file: str):
    """File name for sample collection info (resume functionality)."""
    return out_raw_file + "_info"


def collection_lock_file_name(out_raw_file: str):
    """File name for sample collection lock (resume functionality)."""
    return out_raw_file + "_lock"


# Model construction auxiliary files

def read_stat_file_name(saves_file: str):
    """File name for read assignment statistics."""
    return saves_file + "_read_stat"


def transcript_stat_file_name(saves_file: str):
    """File name for transcript model statistics."""
    return saves_file + "_transcript_stat"


def read_group_file_from_saves(saves_prefix: str):
    """Convert saves file prefix to read_group file prefix."""
    if saves_prefix.endswith('.save'):
        return saves_prefix.replace('.save', '.read_group')
    return None


# UMI filtering output files

def umi_output_prefix(out_umi_filtered: str, edit_distance: int):
    """Get UMI filtered output prefix with edit distance encoding."""
    if edit_distance < 0:
        return out_umi_filtered + ".ALL"
    return out_umi_filtered + ".ED%d" % edit_distance


# Counter output files

def counts_prefix(output_prefix: str):
    """Get counts file prefix from output prefix."""
    return output_prefix + "_counts"


def tpm_prefix(output_prefix: str):
    """Get TPM file prefix from output prefix."""
    return output_prefix + "_tpm"


def counts_file_name(counts_prefix: str, linear: bool = False):
    """Get counts file name with optional linear format."""
    return counts_prefix + (".linear.tsv" if linear else ".tsv")


def tpm_file_name(tpm_prefix: str):
    """Get TPM file name."""
    return tpm_prefix + ".tsv"


def counts_stats_file_name(counts_file: str):
    """Get stats file name from counts file."""
    return counts_file + ".stats"


def counts_usable_file_name(counts_file: str):
    """Get usable reads file name from counts file."""
    return counts_file + ".usable"


# Matrix Market format files (for grouped counts)

def mtx_matrix_file(output_prefix: str):
    """Get Matrix Market file name."""
    return output_prefix + ".matrix.mtx"


def mtx_features_file(output_prefix: str):
    """Get features file name for Matrix Market format."""
    return output_prefix + ".features.tsv"


def mtx_barcodes_file(output_prefix: str):
    """Get barcodes file name for Matrix Market format."""
    return output_prefix + ".barcodes.tsv"