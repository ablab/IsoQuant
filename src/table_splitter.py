############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import pysam
from concurrent.futures import ProcessPoolExecutor
import logging
import gzip
import os


logger = logging.getLogger('IsoQuant')


def get_chromosome_read_counts(bam_files, chromosomes):
    """
    Quickly get read counts per chromosome from BAM indices.
    Uses pysam.get_index_statistics() which reads from .bai without scanning BAM.

    Returns: dict mapping chr_id -> total_read_count
    """
    chr_counts = {chr_id: 0 for chr_id in chromosomes}

    for bam_file in bam_files:
        try:
            bam = pysam.AlignmentFile(bam_file, "rb")
            try:
                # Get statistics from index
                stats = bam.get_index_statistics()
                for stat in stats:
                    chr_id = bam.get_reference_name(stat.contig)
                    if chr_id in chr_counts:
                        # Use mapped + unmapped counts
                        chr_counts[chr_id] += stat.mapped + stat.unmapped
            finally:
                bam.close()
        except Exception as e:
            logger.warning(f"Could not read index statistics from {bam_file}: {e}")
            # Fallback: estimate from BAM header
            try:
                bam = pysam.AlignmentFile(bam_file, "rb")
                for chr_id in chromosomes:
                    if chr_id in bam.references:
                        # Rough estimate: count via fetch (slower but works without index stats)
                        try:
                            count = bam.count(chr_id)
                            chr_counts[chr_id] += count
                        except:
                            pass
                bam.close()
            except:
                pass

    return chr_counts


def distribute_chromosomes_weighted(chromosomes, num_workers, chr_read_counts):
    """
    Distribute chromosomes to workers using greedy bin packing.
    Balances load by assigning large chromosomes first.

    Args:
        chromosomes: list of chromosome IDs
        num_workers: number of worker processes
        chr_read_counts: dict mapping chr_id -> read_count

    Returns:
        List of lists, where workers[i] contains chromosome IDs for worker i
    """
    if not chromosomes:
        return [[] for _ in range(num_workers)]

    # Sort chromosomes by read count (descending)
    sorted_chrs = sorted(chromosomes, key=lambda c: chr_read_counts.get(c, 0), reverse=True)

    # Initialize workers
    workers = [[] for _ in range(num_workers)]
    worker_loads = [0] * num_workers

    # Greedy assignment: assign each chromosome to least-loaded worker
    for chr_id in sorted_chrs:
        min_worker_idx = worker_loads.index(min(worker_loads))
        workers[min_worker_idx].append(chr_id)
        worker_loads[min_worker_idx] += chr_read_counts.get(chr_id, 0)

    # Log distribution for debugging
    for i, (chrs, load) in enumerate(zip(workers, worker_loads)):
        if chrs:
            logger.debug(f"Worker {i}: {len(chrs)} chromosomes, {load:,} reads")

    return workers


def collect_chromosome_reads(chr_id, bam_files):
    """
    Collect all read IDs for a specific chromosome by scanning BAM files.
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
        except Exception as e:
            logger.warning(f"Could not fetch reads from {bam_file} for chromosome {chr_id}: {e}")

    return read_ids


def process_table_for_chromosomes(worker_id, input_tsvs, my_chromosomes, bam_files,
                                  output_files, read_column, group_columns, delim):
    """
    Worker function: processes entire table for assigned chromosomes.

    Each worker:
    1. Builds read ID cache for its assigned chromosomes
    2. Streams through table line-by-line (memory efficient)
    3. Writes matching reads to chromosome-specific output files

    Args:
        worker_id: Worker identifier for logging
        input_tsvs: TSV file(s) to read from
        my_chromosomes: List of chromosome IDs assigned to this worker
        bam_files: List of BAM files to scan for read IDs
        output_files: Dict mapping chr_id -> output_file_path
        read_column: Column index for read ID
        group_columns: Tuple of column indices for group values
        delim: Column delimiter
    """
    if not my_chromosomes:
        return 0, 0

    logger.debug(f"Worker {worker_id}: processing {len(my_chromosomes)} chromosomes: {', '.join(my_chromosomes)}")

    # Step 1: Build read ID cache for assigned chromosomes
    logger.debug(f"Worker {worker_id}: building read ID cache...")
    read_cache = {}
    for chr_id in my_chromosomes:
        read_cache[chr_id] = collect_chromosome_reads(chr_id, bam_files)
        logger.debug(f"Worker {worker_id}: cached {len(read_cache[chr_id]):,} reads for {chr_id}")

    # Step 2: Open output files
    out_handles = {}
    for chr_id in my_chromosomes:
        out_handles[chr_id] = open(output_files[chr_id], 'w')

    # Step 3: Stream through table line-by-line
    total_reads_processed = 0
    total_reads_written = 0
    min_columns = max(read_column, max(group_columns)) + 1

    try:
        input_files = input_tsvs if isinstance(input_tsvs, list) else [input_tsvs]

        for input_file in input_files:
            # Handle gzipped files
            _, ext = os.path.splitext(input_file)
            if ext.lower() in ['.gz', '.gzip']:
                file_handle = gzip.open(input_file, 'rt')
            else:
                file_handle = open(input_file)

            try:
                for line in file_handle:
                    if line.startswith("#") or not line.strip():
                        continue

                    total_reads_processed += 1
                    columns = line.rstrip('\n').split(delim)

                    if len(columns) < min_columns:
                        continue

                    read_id = columns[read_column]

                    # Check if this read belongs to any of my chromosomes
                    for chr_id in my_chromosomes:
                        if read_id in read_cache[chr_id]:
                            # Extract group values
                            group_vals = delim.join(columns[c] for c in group_columns)
                            out_handles[chr_id].write(f"{read_id}\t{group_vals}\n")
                            total_reads_written += 1
                            break  # Read can only be on one chromosome
            finally:
                file_handle.close()
    finally:
        # Close output files
        for f in out_handles.values():
            f.close()

    logger.debug(f"Worker {worker_id}: processed {total_reads_processed:,} table entries, "
                 f"wrote {total_reads_written:,} reads")

    return total_reads_processed, total_reads_written


def split_read_table_parallel(sample, input_tsvs, split_reads_file_names, num_threads,
                              read_column=0, group_columns=(1,), delim='\t'):
    """
    Improved parallel table splitting algorithm.

    Strategy: Assign chromosomes to workers, each worker streams table line-by-line for its chromosomes.

    Advantages over old algorithm:
    - No redundant work (each line processed once by each worker)
    - No chunking overhead (line-by-line streaming)
    - Minimal memory usage (no intermediate dicts)
    - Better memory distribution (each worker caches only its chromosomes)
    - True parallelism (no generator bottleneck)

    Args:
        sample: Sample object with file_list
        input_tsvs: TSV file(s) containing read groups
        split_reads_file_names: Dict mapping chr_id -> output_file_path
        num_threads: Number of worker processes
        read_column: Column index for read ID (default: 0)
        group_columns: Tuple of column indices for group values (default: (1,))
        delim: Column delimiter (default: '\t')
    """
    logger.info(f"Splitting table {input_tsvs} across {num_threads} workers")

    bam_files = list(map(lambda x: x[0], sample.file_list))
    chromosomes = list(split_reads_file_names.keys())

    if not chromosomes:
        logger.warning("No chromosomes to process")
        return

    # Step 1: Get chromosome read counts from BAM indices for load balancing
    logger.info("Analyzing chromosome sizes from BAM indices...")
    chr_read_counts = get_chromosome_read_counts(bam_files, chromosomes)

    total_reads = sum(chr_read_counts.values())
    logger.info(f"Total reads across {len(chromosomes)} chromosomes: {total_reads:,}")

    # Step 2: Distribute chromosomes to workers (weighted by read count)
    num_workers = min(num_threads, len(chromosomes))  # No point having more workers than chromosomes
    chr_assignments = distribute_chromosomes_weighted(chromosomes, num_workers, chr_read_counts)

    # Filter out empty assignments
    chr_assignments = [chrs for chrs in chr_assignments if chrs]
    num_workers = len(chr_assignments)

    logger.info(f"Using {num_workers} workers to process {len(chromosomes)} chromosomes")

    # Step 3: Initialize output files (truncate)
    for chr_id in chromosomes:
        open(split_reads_file_names[chr_id], 'w').close()

    # Step 4: Launch workers
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for worker_id, my_chrs in enumerate(chr_assignments):
            future = executor.submit(
                process_table_for_chromosomes,
                worker_id,
                input_tsvs,
                my_chrs,
                bam_files,
                split_reads_file_names,
                read_column,
                group_columns,
                delim
            )
            futures.append(future)

        # Wait for all workers and collect results
        total_processed = 0
        total_written = 0
        for future in futures:
            processed, written = future.result()
            total_processed += processed
            total_written += written

    logger.info(f"Table splitting complete: {total_written:,} reads written to {len(chromosomes)} chromosome files")
