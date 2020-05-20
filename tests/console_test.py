import subprocess

import pytest


def test_run_without_parameters():
    result = subprocess.run(["python", "isoquant.py"], capture_output=True)
    assert result.returncode == 2
    assert b"usage" in result.stderr
    assert b"isoquant.py: error: the following arguments are required: --data_type/-d, --genedb/-g" in result.stderr


@pytest.mark.parametrize("option", ["-h", "--help", "--full-help"])
def test__help(option):
    result = subprocess.run(["python", "isoquant.py", option], capture_output=True)
    assert result.returncode == 0
    assert b"usage" in result.stdout
    assert b"optional arguments:" in result.stdout
