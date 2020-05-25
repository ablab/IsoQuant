import os
import subprocess
from collections import namedtuple

import pytest

from src.common import get_path_to_program


class TestGetPathToProgram:
    def file_exist(self, path):
        return True

    def file_doesnot_exist(self, path):
        return False

    def is_access_true(self, path, mode):
        return True

    def test_in_dir(self, monkeypatch):
        def correct_versions(*args, **kwargs):
            out = namedtuple("CompletedProcess", ("stdout"))
            return out("1.38")

        monkeypatch.setattr(subprocess, "run", correct_versions)
        monkeypatch.setattr(os.path, "isfile", self.file_exist)
        monkeypatch.setattr(os, "access", self.is_access_true)
        program_name = "my_prog"
        dir = "/folder"
        assert os.path.join(dir, program_name) == get_path_to_program(program_name, dir, "1.38")

    @pytest.mark.parametrize("min_version", ["1.39", "2.0"])
    def test_wrong_version(self, monkeypatch, min_version):
        def correct_versions(*args, **kwargs):
            out = namedtuple("CompletedProcess", ("stdout"))
            return out("1.38")

        monkeypatch.setattr(subprocess, "run", correct_versions)
        monkeypatch.setattr(os.path, "isfile", self.file_exist)
        monkeypatch.setattr(os, "access", self.is_access_true)
        program_name = "my_prog"
        dir = "/folder"
        assert None is get_path_to_program(program_name, dir, "1.39")

    @pytest.mark.parametrize("min_version", ["1.39", "2.0"])
    def test_not_in_path(self, monkeypatch, min_version):
        monkeypatch.setenv("PATH", os.pathsep.join(["path1", "path2", "path3"]))
        monkeypatch.setattr(os.path, "isfile", self.file_doesnot_exist)
        monkeypatch.setattr(os, "access", self.is_access_true)
        program_name = "my_prog"
        assert None is get_path_to_program(program_name, "1.39")

    # TODO: add more tests
