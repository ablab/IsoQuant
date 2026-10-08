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
- a fresh run wipes the stores (reset) and writes a run record (run.json: checkpoint format,
  run id, IsoQuant version); the run id is also saved in .params, and --resume refuses a store
  whose record is missing, of another format, or of another run;
- on resume a missing or unreadable marker means "not done";
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
import uuid
from enum import Enum
from typing import Any, Callable, Dict, Optional, Type

from isoquant_lib.utils.file_naming import convert_chr_id_to_file_name_str

logger = logging.getLogger('IsoQuant')

Payload = Optional[Dict[str, Any]]

MARKER_SUFFIX = ".done"
DEBUG_FAIL_ENV = "ISOQUANT_DEBUG_FAIL"
RUN_CHECKPOINT_DIR = "checkpoints"
RUN_RECORD = "run.json"
# bump whenever stage names, marker contents or the files a stage leaves behind change:
# --resume then refuses stores written by other formats instead of trusting their markers
CHECKPOINT_FORMAT = 1


class DebugFailure(RuntimeError):
    """Raised by the ISOQUANT_DEBUG_FAIL hook to emulate a crash at a stage boundary."""


def marker_name(*parts: str) -> str:
    """Build a marker name from parts; every part is made safe for a file name."""
    return "/".join(convert_chr_id_to_file_name_str(str(p)) for p in parts)


def _write_json_atomically(path: str, content: Dict[str, Any]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = "%s.tmp.%d" % (path, os.getpid())
    with open(tmp, "w") as f:
        json.dump(content, f)
    os.replace(tmp, path)


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

    def _read(self, name: str) -> Optional[Dict[str, Any]]:
        """The marker's content, or None when it is missing or unreadable.

        Markers are replaced atomically, so a killed process never leaves a partial one;
        a damaged file (e.g. after a power loss, nothing is fsynced) only reruns the stage.
        """
        path = self.path(name)
        if not os.path.exists(path):
            return None
        try:
            with open(path) as f:
                content = json.load(f)
        except (OSError, ValueError) as e:
            logger.warning("Checkpoint marker %s is unreadable (%s), the stage will run again" % (path, e))
            return None
        return content if isinstance(content, dict) else None

    def is_done(self, name: str) -> bool:
        return self.resume and self._read(name) is not None

    def load(self, name: str) -> Payload:
        if not self.resume:
            return None
        content = self._read(name)
        return content.get("payload") if content else None

    def mark_done(self, name: str, payload: Payload = None) -> None:
        _write_json_atomically(self.path(name),
                               {"finished_at": datetime.datetime.now().isoformat(timespec="seconds"),
                                "payload": payload})


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


def enum_stats_to_payload(stats_dict: Dict[Any, int]) -> Dict[str, int]:
    return {key.name: int(value) for key, value in stats_dict.items()}


def payload_to_enum_stats(payload: Dict[str, int], enum_class: Type[Enum]) -> Dict[Enum, int]:
    return {enum_class[name]: int(value) for name, value in payload.items()}


def run_checkpoint_dir(output_dir: str) -> str:
    """The run-level store; sample-level stores live in <sample>/aux/checkpoints."""
    return os.path.join(output_dir, RUN_CHECKPOINT_DIR)


def write_run_record(run_dir: str, isoquant_version: str) -> str:
    """Record which run owns a freshly reset run store; returns the run id to save in .params."""
    run_id = uuid.uuid4().hex
    _write_json_atomically(os.path.join(run_dir, RUN_RECORD),
                           {"checkpoint_format": CHECKPOINT_FORMAT, "run_id": run_id,
                            "isoquant_version": isoquant_version,
                            "started_at": datetime.datetime.now().isoformat(timespec="seconds")})
    return run_id


def check_resumable(run_dir: str, run_id: Optional[str], isoquant_version: str) -> Optional[str]:
    """Why the checkpoints in run_dir cannot be trusted for this run, or None if they can.

    run_id comes from .params. It is missing when another IsoQuant version wrote .params
    (e.g. an older one ran in this folder after a newer one), and differs when the store
    belongs to another run; either way the markers describe outputs this run did not make.
    """
    record_path = os.path.join(run_dir, RUN_RECORD)
    try:
        with open(record_path) as f:
            record = json.load(f)
    except (OSError, ValueError):
        return "it has no readable checkpoint record (%s); it was probably started by an older IsoQuant version" \
            % record_path
    if not isinstance(record, dict) or record.get("checkpoint_format") != CHECKPOINT_FORMAT:
        return "its checkpoints have format %s, this IsoQuant version uses format %d" \
            % (record.get("checkpoint_format") if isinstance(record, dict) else "unknown", CHECKPOINT_FORMAT)
    if not run_id or record.get("run_id") != run_id:
        return "its checkpoints do not belong to the run described by .params " \
               "(the folder was probably reused by another IsoQuant version)"
    if record.get("isoquant_version") != isoquant_version:
        logger.warning("The run was started by IsoQuant %s and is resumed by IsoQuant %s; stages finished "
                       "before will not be recomputed" % (record.get("isoquant_version"), isoquant_version))
    return None
