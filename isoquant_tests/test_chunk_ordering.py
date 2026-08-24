############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Chunk results reach the caller in completion order, carrying the index to reorder by."""

import time
from types import SimpleNamespace

from isoquant_lib.barcode_calling.detect_barcodes import run_chunks_in_parallel

CHUNKS = 6
THREADS = 3


def finish_in_reverse(index: int) -> int:
    """Make later chunks finish first, so completion order is not submission order."""
    time.sleep(0.05 * (CHUNKS - index))
    return index


class TestRunChunksInParallel:
    def test_every_chunk_is_handled_once_with_its_index(self):
        completion_order = []
        collected = {}

        def submit(pool, chunk, num):
            return pool.submit(finish_in_reverse, num)

        def handle_result(result, chunk_index):
            completion_order.append(chunk_index)
            collected[chunk_index] = result
            return 1

        run_chunks_in_parallel(iter(range(CHUNKS)), SimpleNamespace(threads=THREADS),
                               None, submit, handle_result)

        assert sorted(collected) == list(range(CHUNKS))
        # the index really identifies its own chunk, so a caller can key on it
        assert all(index == result for index, result in collected.items())
        # and it was needed: the chunks did not come back in the order they went out
        assert completion_order != sorted(completion_order)

    def test_submission_order_is_recoverable(self):
        """What _process_single_file_in_parallel does to keep its merge deterministic."""
        outputs = {}

        def submit(pool, chunk, num):
            return pool.submit(finish_in_reverse, num)

        def handle_result(result, chunk_index):
            outputs[chunk_index] = "chunk_%d" % result
            return 1

        run_chunks_in_parallel(iter(range(CHUNKS)), SimpleNamespace(threads=THREADS),
                               None, submit, handle_result)

        assert [outputs[i] for i in sorted(outputs)] == ["chunk_%d" % i for i in range(CHUNKS)]
