import subprocess

import pytest


def test_run_without_parameters():
    result = subprocess.run(["python", "isoquant.py"], capture_output=True)
    assert result.returncode == 2
    assert b"usage" in result.stderr
    assert b"error: the following arguments are required: --genedb/-g, --data_type/-d" in result.stderr


@pytest.mark.parametrize("option", ["-h", "--help", "--full_help"])
def test__help(option):
    result = subprocess.run(["python", "isoquant.py", option], capture_output=True)
    assert result.returncode == 0
    assert b"usage" in result.stdout
    assert b"optional arguments:" in result.stdout
