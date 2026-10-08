############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from __future__ import annotations

import json
import os
import pickle
from enum import Enum, unique

import pytest

from isoquant_lib.utils.checkpoints import (
    CHECKPOINT_FORMAT,
    CheckpointStore,
    DebugFailure,
    DEBUG_FAIL_ENV,
    Payload,
    RUN_RECORD,
    check_debug_failure,
    check_resumable,
    enum_stats_to_payload,
    marker_name,
    payload_to_enum_stats,
    run_stage,
    write_run_record,
)


@unique
class Colour(Enum):
    red = 1
    blue = 2


def fresh_and_resumed(root: str) -> tuple[CheckpointStore, CheckpointStore]:
    fresh = CheckpointStore(root, resume=False)
    fresh.reset()
    return fresh, CheckpointStore(root, resume=True)


class TestCheckpointStore:
    def test_payload_alias_is_runtime_safe(self):
        assert Payload is not None

    def test_fresh_store_never_reports_done(self, tmp_path):
        fresh, resumed = fresh_and_resumed(str(tmp_path / "cp"))
        fresh.mark_done("stage")
        assert not fresh.is_done("stage")
        assert resumed.is_done("stage")

    def test_payload_round_trip(self, tmp_path):
        fresh, resumed = fresh_and_resumed(str(tmp_path / "cp"))
        fresh.mark_done("collect", {"a": 1, "b": {"x": 2}})
        assert resumed.load("collect") == {"a": 1, "b": {"x": 2}}
        assert resumed.load("missing") is None

    def test_nested_names_and_no_tmp_left(self, tmp_path):
        fresh, resumed = fresh_and_resumed(str(tmp_path / "cp"))
        fresh.mark_done(marker_name("umi_bc2bc", "col0", "ED3", "chr1"))
        assert resumed.is_done("umi_bc2bc/col0/ED3/chr1")
        leftovers = [f for _, _, files in os.walk(fresh.root) for f in files if ".tmp." in f]
        assert leftovers == []

    def test_reset_wipes_markers(self, tmp_path):
        fresh, resumed = fresh_and_resumed(str(tmp_path / "cp"))
        fresh.mark_done("a/b")
        fresh.reset()
        assert not resumed.is_done("a/b")
        assert os.path.isdir(fresh.root)

    def test_picklable(self, tmp_path):
        store = CheckpointStore(str(tmp_path), resume=True)
        clone = pickle.loads(pickle.dumps(store))
        assert clone.root == store.root and clone.resume

    def test_unreadable_marker_is_not_done(self, tmp_path):
        fresh, resumed = fresh_and_resumed(str(tmp_path / "cp"))
        fresh.mark_done("collect", {"n": 1})
        with open(fresh.path("collect"), "w") as f:
            f.write('{"payl')  # e.g. truncated by a power loss
        assert not resumed.is_done("collect")
        assert resumed.load("collect") is None
        open(fresh.path("collect"), "w").close()
        assert not resumed.is_done("collect")

    def test_marker_name_sanitises_every_part(self):
        assert marker_name("collect", "chrUn/random") == "collect/chrUn_random"
        assert marker_name("sample", "a/b", "c") == "sample/a_b/c"


class TestRunStage:
    def test_runs_marks_and_cleans(self, tmp_path):
        fresh, resumed = fresh_and_resumed(str(tmp_path / "cp"))
        calls = []
        run_stage(fresh, "s", run=lambda: calls.append("run") or {"n": 3},
                  cleanup=lambda: calls.append("cleanup"))
        assert calls == ["run", "cleanup"]
        assert resumed.load("s") == {"n": 3}

    def test_skip_restores_and_cleans_again(self, tmp_path):
        fresh, resumed = fresh_and_resumed(str(tmp_path / "cp"))
        fresh.mark_done("s", {"n": 3})
        calls = []
        run_stage(resumed, "s", run=lambda: calls.append("run"),
                  restore=lambda p: calls.append(("restore", p)),
                  cleanup=lambda: calls.append("cleanup"))
        assert calls == [("restore", {"n": 3}), "cleanup"]

    def test_disabled_does_nothing(self, tmp_path):
        fresh, resumed = fresh_and_resumed(str(tmp_path / "cp"))
        calls = []
        run_stage(fresh, "s", run=lambda: calls.append("run"), enabled=False)
        assert calls == [] and not resumed.is_done("s")

    def test_keep_tmp_skips_cleanup(self, tmp_path):
        fresh, _ = fresh_and_resumed(str(tmp_path / "cp"))
        calls = []
        run_stage(fresh, "s", run=lambda: None, cleanup=lambda: calls.append("cleanup"), keep_tmp=True)
        assert calls == []

    def test_failure_leaves_no_marker_and_no_cleanup(self, tmp_path):
        fresh, resumed = fresh_and_resumed(str(tmp_path / "cp"))
        calls = []

        def boom():
            raise RuntimeError("crash")

        with pytest.raises(RuntimeError):
            run_stage(fresh, "s", run=boom, cleanup=lambda: calls.append("cleanup"))
        assert calls == [] and not resumed.is_done("s")


class TestDebugHook:
    def test_after_is_default(self, tmp_path, monkeypatch):
        fresh, resumed = fresh_and_resumed(str(tmp_path / "cp"))
        monkeypatch.setenv(DEBUG_FAIL_ENV, "s")
        with pytest.raises(DebugFailure):
            run_stage(fresh, "s", run=lambda: None)
        assert resumed.is_done("s")

    def test_before_leaves_no_marker(self, tmp_path, monkeypatch):
        fresh, resumed = fresh_and_resumed(str(tmp_path / "cp"))
        monkeypatch.setenv(DEBUG_FAIL_ENV, "s:before")
        with pytest.raises(DebugFailure):
            run_stage(fresh, "s", run=lambda: None)
        assert not resumed.is_done("s")

    def test_colon_in_chromosome_name(self, monkeypatch):
        monkeypatch.setenv(DEBUG_FAIL_ENV, "collect/HLA-A*01:01")
        with pytest.raises(DebugFailure):
            check_debug_failure("collect/HLA-A*01:01", "after")
        check_debug_failure("collect/HLA-A*01:01", "before")
        monkeypatch.setenv(DEBUG_FAIL_ENV, "collect/HLA-A*01:01:before")
        with pytest.raises(DebugFailure):
            check_debug_failure("collect/HLA-A*01:01", "before")

    def test_other_stage_unaffected(self, monkeypatch):
        monkeypatch.setenv(DEBUG_FAIL_ENV, "other")
        check_debug_failure("s", "after")


def test_enum_payload_round_trip():
    stats = {Colour.red: 5, Colour.blue: 0}
    assert payload_to_enum_stats(enum_stats_to_payload(stats), Colour) == stats


class TestRunRecord:
    def test_own_record_is_resumable(self, tmp_path):
        run_id = write_run_record(str(tmp_path), "4.1.0")
        assert check_resumable(str(tmp_path), run_id, "4.1.0") is None

    def test_other_version_only_warns(self, tmp_path, caplog):
        run_id = write_run_record(str(tmp_path), "4.1.0")
        assert check_resumable(str(tmp_path), run_id, "4.2.0") is None

    def test_missing_record(self, tmp_path):
        assert "no readable checkpoint record" in check_resumable(str(tmp_path), "abc", "4.1.0")

    def test_params_of_another_run(self, tmp_path):
        write_run_record(str(tmp_path), "4.1.0")
        assert "do not belong" in check_resumable(str(tmp_path), "another-run", "4.1.0")
        # .params written by a version without checkpoints has no run id at all
        assert "do not belong" in check_resumable(str(tmp_path), None, "4.1.0")

    def test_other_checkpoint_format(self, tmp_path):
        run_id = write_run_record(str(tmp_path), "4.1.0")
        path = tmp_path / RUN_RECORD
        record = json.loads(path.read_text())
        record["checkpoint_format"] = CHECKPOINT_FORMAT + 1
        path.write_text(json.dumps(record))
        assert "format" in check_resumable(str(tmp_path), run_id, "4.1.0")

    def test_reset_store_keeps_no_record(self, tmp_path):
        run_dir = str(tmp_path / "checkpoints")
        write_run_record(run_dir, "4.1.0")
        CheckpointStore(run_dir, resume=False).reset()
        assert not os.path.exists(os.path.join(run_dir, RUN_RECORD))
