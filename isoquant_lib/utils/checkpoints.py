############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Resume checkpoints.

A run is a fixed sequence of stages. Each finished stage leaves a marker file in a
checkpoint directory, and a resumed run (--resume) skips every stage whose marker
exists. Two directories are used: <output>/checkpoints/ for run-level stages and
<sample>/aux/checkpoints/ for the stages of one experiment.

Rules (see .claude/RESUME_CHECKPOINTS.md):
- a fresh run wipes the stores (reset), so on resume a missing marker means "not done";
- a stage writes all of its outputs and closes them before its marker is written;
- files consumed only by one stage are removed by its cleanup, after the marker;
- state later code needs is returned as a small JSON payload and re-applied by restore.
"""

from __future__ import annotations

import datetime
import json
import logging
import os
import shutil
from enum import Enum
from typing import Any, Callable, Dict, Iterable, Optional, Type

from isoquant_lib.utils.file_naming import convert_chr_id_to_file_name_str

logger = logging.getLogger('IsoQuant')

Payload = Optional[Dict[str, Any]]

MARKER_SUFFIX = ".done"
DEBUG_FAIL_ENV = "ISOQUANT_DEBUG_FAIL"
RUN_CHECKPOINT_DIR = "checkpoints"


class DebugFailure(RuntimeError):
    """Raised by the ISOQUANT_DEBUG_FAIL hook to emulate a crash at a stage boundary."""


def marker_name(*parts: str) -> str:
    """Build a marker name from parts; every part is made safe for a file name."""
    return "/".join(convert_chr_id_to_file_name_str(str(p)) for p in parts)


def _isoquant_version() -> str:
    try:
        from importlib.metadata import version
        return version("isoquant")
    except Exception:
        return "unknown"


class CheckpointStore:
    """A directory of marker files.

    Holds only a path and the resume flag, so it is cheap to pickle and can be
    passed to ProcessPoolExecutor workers; every marker is a separate file written
    atomically, so workers can mark their chromosomes concurrently.
    """

    def __init__(self, root: str, resume: bool) -> None:
        self.root = root
        self.resume = resume

    def path(self, name: str) -> str:
        return os.path.join(self.root, name + MARKER_SUFFIX)

    def reset(self) -> None:
        """Remove all markers; called once by the driver at the start of a fresh run."""
        shutil.rmtree(self.root, ignore_errors=True)
        os.makedirs(self.root, exist_ok=True)

    def ensure_dir(self, name: str) -> None:
        os.makedirs(os.path.join(self.root, name), exist_ok=True)

    def is_done(self, name: str) -> bool:
        return self.resume and os.path.exists(self.path(name))

    def load(self, name: str) -> Payload:
        if not self.is_done(name):
            return None
        with open(self.path(name)) as f:
            return json.load(f).get("payload")

    def mark_done(self, name: str, payload: Payload = None) -> None:
        marker = self.path(name)
        os.makedirs(os.path.dirname(marker), exist_ok=True)
        tmp = "%s.tmp.%d" % (marker, os.getpid())
        with open(tmp, "w") as f:
            json.dump({"isoquant_version": _isoquant_version(),
                       "finished_at": datetime.datetime.now().isoformat(timespec="seconds"),
                       "payload": payload}, f)
        os.replace(tmp, marker)


def check_debug_failure(name: str, position: str) -> None:
    """Raise when ISOQUANT_DEBUG_FAIL=<name>[:before|:after] matches; default position is after."""
    value = os.environ.get(DEBUG_FAIL_ENV)
    if not value:
        return
    target, sep, pos = value.rpartition(":")
    if not sep or pos not in ("before", "after"):
        # chromosome names may contain ':', the suffix only counts if it is a position
        target, pos = value, "after"
    if target == name and pos == position:
        raise DebugFailure("%s=%s: emulated crash %s marking %s" % (DEBUG_FAIL_ENV, value, position, name))


def mark_stage_done(store: CheckpointStore, name: str, payload: Payload = None) -> None:
    """mark_done wrapped by the debug hook; per-chromosome workers use it directly."""
    check_debug_failure(name, "before")
    store.mark_done(name, payload)
    check_debug_failure(name, "after")


def run_stage(store: CheckpointStore, name: str, run: Callable[[], Payload],
              restore: Callable[[Payload], None] | None = None,
              cleanup: Callable[[], None] | None = None,
              enabled: bool = True, keep_tmp: bool = False) -> None:
    """Run a stage unless the previous run finished it.

    cleanup removes files consumed only by this stage; it runs after the marker and
    again when the stage is skipped, so it must be idempotent.
    """
    if not enabled:
        return
    if store.is_done(name):
        logger.info("%s: done in the previous run, skipping" % name)
        if restore is not None:
            restore(store.load(name))
    else:
        payload = run()
        mark_stage_done(store, name, payload)
    if cleanup is not None and not keep_tmp:
        cleanup()


def remove_files(file_names: Iterable[str]) -> None:
    for f in file_names:
        if os.path.exists(f):
            os.remove(f)


def enum_stats_to_payload(stats_dict: Dict[Any, int]) -> Dict[str, int]:
    return {key.name: int(value) for key, value in stats_dict.items()}


def payload_to_enum_stats(payload: Dict[str, int], enum_class: Type[Enum]) -> Dict[Enum, int]:
    return {enum_class[name]: int(value) for name, value in payload.items()}


def run_checkpoint_dir(output_dir: str) -> str:
    """The run-level store; sample-level stores live in <sample>/aux/checkpoints."""
    return os.path.join(output_dir, RUN_CHECKPOINT_DIR)
