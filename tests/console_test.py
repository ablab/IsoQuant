import subprocess

import pytest


def test_run_without_parameters():
    result = subprocess.run(["python", "IsoQuant/isoquant.py"], capture_output=True)
    assert result.returncode == 255
    assert result.stdout == b'Output folder was not specified\n'
