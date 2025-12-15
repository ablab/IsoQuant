############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import pysam
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import logging

from .common import proper_plural_form


logger = logging.getLogger('IsoQuant')

CHUNK_SIZE = 500000

def load_barcodes_chunked(in_file, chunk_size=CHUNK_SIZE, use_untrusted_umis=False,
                          barcode_column=1, umi_column=2):
    """Generator that yields chunks of barcodes lazily"""
    in_files = in_file if isinstance(in_file, list) else [in_file]

    stats = {
        'read_count': 0,
        'barcoded': 0,
        'hq_barcoded': 0,
        'trusted_umi': 0
    }

    barcode_chunk = {}

    for f in in_files:
        logger.info(f"Loading barcodes from {f}")
        max_columns = max(barcode_column, umi_column)

        with open(f) as file:
            for line in file:
                if line.startswith("#"):
                    continue

                stats['read_count'] += 1
                v = line.strip().split("\t")

                if len(v) < max_columns:
                    continue

                barcode = v[barcode_column]
                if barcode == "*":
                    continue


                stats['barcoded'] += 1
                umi = v[umi_column]
                barcode_chunk[v[0]] = (barcode, umi)

                # Yield chunk when it reaches the specified size
                if len(barcode_chunk) >= chunk_size:
                    yield barcode_chunk
                    barcode_chunk = {}

    # Yield remaining barcodes
    if barcode_chunk:
        yield barcode_chunk

    logger.info(f"Total reads: {stats['read_count']}")
    logger.info(f"Barcoded: {stats['barcoded']}")


def load_table_chunked(in_file, chunk_size=CHUNK_SIZE, read_column=0, group_columns=(1,)):
    """Generator that yields chunks of barcodes lazily"""
    in_files = in_file if isinstance(in_file, list) else [in_file]

    read_chunk = {}
    malformed_lines = 0
    min_columns = max(read_column, max(group_columns))

    for f in in_files:
        logger.info(f"Loading barcodes from {f}")

        with open(f) as file:
            for line in file:
                if line.startswith("#"):
                    continue

                v = line.strip().split("\t")

                if len(v) < min_columns:
                    malformed_lines += 1
                    continue

                read_chunk[v[read_column]] = tuple(v[c] for c in group_columns)

                # Yield chunk when it reaches the specified size
                if len(read_chunk) >= chunk_size:
                    yield read_chunk
                    read_chunk = {}

    # Yield remaining barcodes
    if read_chunk:
        yield read_chunk

    if malformed_lines > 0:
        logger.warning(f"Read group file {in_files} has {malformed_lines} malformed lines")


def collect_chromosome_reads(chr_id, bam_files):
    """
    Collect all read IDs for a specific chromosome.
    This is done once per chromosome at the start.
    """
    read_ids = set()
    for bam_file in bam_files:
        try:
            bam = pysam.AlignmentFile(bam_file, "rb")
            try:
                for read in bam.fetch(chr_id):
                    read_ids.add(read.query_name)
            finally:
                bam.close()
        except:
            pass

    return read_ids


def process_barcode_chunk_with_read_cache(chr_id, read_chunk, read_ids_cache, output_file):
    """
    Process a chunk of barcodes using pre-cached read IDs.
    Much faster than reading BAM files each time.
    """
    with open(output_file, "a") as out_f:
        for read_id, vals in read_chunk.items():
            if read_id in read_ids_cache:
                out_f.write("%s\t%s\n" % (read_id, "\t".join(vals) + "\n"))

    return chr_id, len(read_chunk)


def combine_generators(chromosome_list, read_chung_generator):
    while True:
        try:
            read_chunk = next(read_chung_generator)
            for chr_id in chromosome_list:
                yield chr_id, read_chunk
        except StopIteration:
            return


def split_read_table_parallel(sample,
                              input_tsvs,
                              split_reads_file_names,
                              num_threads,
                              load_func=load_table_chunked):
    """
    Parallel version using ProcessPoolExecutor with lazy chunk processing.

    Args:
        sample: Sample object with barcoded_reads and file_list
        split_reads_file_names: Dict mapping chr_id to output file path
        num_threads: Number of parallel workers
        load_func: function to load read chunks taking one argument, a list of file paths
        chunk_size: Number of barcodes to process in each chunk
    """
    logger.info(f"Loading barcodes from " + str(input_tsvs))
    bam_files = list(map(lambda x: x[0], sample.file_list))

    # Determine which chromosomes are present in BAM files
    chromosomes_to_process = []
    for bam_file in bam_files:
        bam = pysam.AlignmentFile(bam_file, "rb")
        for chr_id in bam.references:
            if chr_id in split_reads_file_names and chr_id not in chromosomes_to_process:
                chromosomes_to_process.append(chr_id)
        bam.close()

    logger.info(f"Splitting barcodes")

    # Initialize output files (truncate if they exist)
    for chr_id in chromosomes_to_process:
        open(split_reads_file_names[chr_id], "w").close()

    # First, cache all read IDs for each chromosome (one-time cost)
    logger.info("Caching read IDs for each chromosome...")
    read_id_cache = {}
    for chr_id in chromosomes_to_process:
        read_id_cache[chr_id] = collect_chromosome_reads(chr_id, bam_files)

    # Create lazy barcode chunk generator
    barcode_chunk_gen = combine_generators(chromosomes_to_process, load_func(input_tsvs))

    # Process barcode chunks with the lazy pattern
    chunks_left = True
    with ProcessPoolExecutor(max_workers=num_threads) as proc:
        future_results = []
        chunk_counter = 0
        total_reads_processed = 0

        # Submit the initial batch of jobs (one per chromosome)
        while chunk_counter < num_threads:
            try:
                chr_id, read_chunk = next(barcode_chunk_gen)
                future_results.append(proc.submit(
                    process_barcode_chunk_with_read_cache,
                    chr_id,
                    read_chunk,
                    read_id_cache[chr_id],
                    split_reads_file_names[chr_id],
                ))
                chunk_counter += 1
            except StopIteration:
                chunks_left = False
                break

        # Process chunks as they complete
        while future_results:
            completed_futures, _ = wait(future_results, return_when=FIRST_COMPLETED)

            for future in completed_futures:
                if future.exception() is not None:
                    raise future.exception()

                chr_id, chunk_size_processed = future.result()
                total_reads_processed += chunk_size_processed
                logger.debug(f"Chromosome {chr_id}: processed chunk")
                future_results.remove(future)

                # Submit next chunk if available
                if chunks_left:
                    try:
                        next_chr, read_chunk = next(barcode_chunk_gen)

                        future_results.append(proc.submit(
                            process_barcode_chunk_with_read_cache,
                            next_chr,
                            read_chunk,
                            read_id_cache[next_chr],
                            split_reads_file_names[next_chr],
                        ))

                        chunk_counter += 1
                    except StopIteration:
                        chunks_left = False

        logger.info(f"Completed processing {chunk_counter} " + proper_plural_form("chunk", chunk_counter))
