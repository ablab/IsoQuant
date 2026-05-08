############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Integration tests for tenX_v3_split and tenX_v2_split barcode detection pipeline.

Tests the full detect_barcodes pipeline (isoquant_detect_barcodes.py) with
split mode on small synthetic concatenated 10x reads, verifying:
- Split reads FASTA is produced
- Barcode TSV has expected columns
- Multi-molecule reads yield multiple output rows
- Quality assessment runs without error
"""

import os
import subprocess
import tempfile

import pytest


ISOQUANT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DETECT_BARCODES = os.path.join(ISOQUANT_DIR, "isoquant_detect_barcodes.py")
ASSESS_QUALITY = os.path.join(ISOQUANT_DIR, "misc/assess_barcode_quality.py")

# 10x v3 constants
R1 = "CTACACGACGCTCTTCCGATCT"           # 22bp
TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT"  # 30bp RC of TSO oligo

# 5 known barcodes (16bp each) and 12bp UMIs
BARCODES = [
    "AAACCCAAGAATCGAT",
    "AAACCCAAGCTCCACG",
    "AAACCCAAGGGTTTCT",
    "AAACGAAAGAAACCAT",
    "AAACGAAAGAAACCCA",
]
UMI = "AACCGGTTAACC"  # 12bp


def make_10x_template(barcode: str, cdna: str, umi: str = UMI) -> str:
    """Build a single-molecule 10x template: R1 + BC + UMI + polyT + RC(cDNA) + TSO."""
    polya = "T" * 30
    rc_map = str.maketrans("ACGT", "TGCA")
    rc_cdna = cdna[::-1].translate(rc_map)
    return R1 + barcode + umi + polya + rc_cdna + TSO


def make_concat_read(molecules: list[tuple[str, str]]) -> str:
    """Concatenate multiple molecules into a single read."""
    return "".join(make_10x_template(bc, cdna) for bc, cdna in molecules)


# Minimal synthetic cDNA sequences (~200bp each)
CDNA1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" \
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" \
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
CDNA2 = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA" \
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA" \
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"


@pytest.fixture()
def whitelist(tmp_path):
    """Write barcode whitelist to a temp file."""
    wl = tmp_path / "whitelist.txt"
    wl.write_text("\n".join(BARCODES) + "\n")
    return str(wl)


@pytest.fixture()
def single_molecule_fasta(tmp_path):
    """FASTA with one single-molecule 10x read per barcode."""
    fasta = tmp_path / "single.fasta"
    lines = []
    for i, bc in enumerate(BARCODES):
        seq = make_10x_template(bc, CDNA1)
        lines.append(f">READ_{i}_{bc}_{UMI}_+")
        lines.append(seq)
    fasta.write_text("\n".join(lines) + "\n")
    return str(fasta)


@pytest.fixture()
def concat_fasta(tmp_path):
    """FASTA with 2-molecule concatenated reads."""
    fasta = tmp_path / "concat.fasta"
    lines = []
    for i in range(3):
        bc1, bc2 = BARCODES[i], BARCODES[i + 1]
        seq = make_concat_read([(bc1, CDNA1), (bc2, CDNA2)])
        lines.append(f">CONCAT_{i}_{bc1}_{UMI}_+_{bc2}_{UMI}_+")
        lines.append(seq)
    fasta.write_text("\n".join(lines) + "\n")
    return str(fasta)


class TestDetectBarcodesScript:
    """Tests for isoquant_detect_barcodes.py CLI with tenX_v3_split mode."""

    def _run_detect(self, fasta, whitelist, output_prefix, mode="tenX_v3_split"):
        cmd = [
            "python3", DETECT_BARCODES,
            "--input", fasta,
            "--mode", mode,
            "--barcodes", whitelist,
            "-t", "1",
            "-o", output_prefix,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result

    def test_single_molecule_runs_successfully(self, single_molecule_fasta, whitelist, tmp_path):
        out = str(tmp_path / "out_single")
        result = self._run_detect(single_molecule_fasta, whitelist, out)
        assert result.returncode == 0, result.stderr

    def test_single_molecule_produces_tsv(self, single_molecule_fasta, whitelist, tmp_path):
        out = str(tmp_path / "out_single")
        self._run_detect(single_molecule_fasta, whitelist, out)
        tsv = out + ".barcoded_reads.tsv"
        assert os.path.exists(tsv), f"Expected {tsv}"

    def test_split_mode_produces_fasta(self, single_molecule_fasta, whitelist, tmp_path):
        """tenX_v3_split mode should produce split_reads.fasta (produces_new_fasta=True)."""
        out = str(tmp_path / "out_single")
        self._run_detect(single_molecule_fasta, whitelist, out)
        split_fasta = out + ".split_reads.fasta"
        assert os.path.exists(split_fasta), f"Expected split FASTA: {split_fasta}"

    def test_tsv_has_tso_start_column(self, single_molecule_fasta, whitelist, tmp_path):
        """TSV header should include tso_start column for tenX_v3_split."""
        out = str(tmp_path / "out_single")
        self._run_detect(single_molecule_fasta, whitelist, out)
        tsv = out + ".barcoded_reads.tsv"
        with open(tsv) as f:
            header = f.readline()
        assert "tso_start" in header

    def test_concat_reads_yield_multiple_rows(self, concat_fasta, whitelist, tmp_path):
        """Concatenated reads should produce more TSV rows than input reads."""
        out = str(tmp_path / "out_concat")
        result = self._run_detect(concat_fasta, whitelist, out)
        assert result.returncode == 0, result.stderr
        tsv = out + ".barcoded_reads.tsv"
        with open(tsv) as f:
            rows = [l for l in f if not l.startswith("read_id")]
        # 3 input reads × 2 molecules each → expect >3 rows
        assert len(rows) > 3, f"Expected >3 rows for concatenated reads, got {len(rows)}"

    def test_barcodes_detected_in_whitelist(self, single_molecule_fasta, whitelist, tmp_path):
        """All detected barcodes should be from the whitelist."""
        out = str(tmp_path / "out_single")
        self._run_detect(single_molecule_fasta, whitelist, out)
        tsv = out + ".barcoded_reads.tsv"
        detected = set()
        with open(tsv) as f:
            for line in f:
                if line.startswith("read_id"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) > 1 and cols[1] != "*":
                    detected.add(cols[1])
        assert detected.issubset(set(BARCODES)), f"Unexpected barcodes: {detected - set(BARCODES)}"

    def test_v3_split_mode_recognized(self, single_molecule_fasta, whitelist, tmp_path):
        """tenX_v3_split is a valid mode (no error about unknown mode)."""
        out = str(tmp_path / "out_v3")
        result = self._run_detect(single_molecule_fasta, whitelist, out, mode="tenX_v3_split")
        assert result.returncode == 0

    def test_v2_split_mode_recognized(self, single_molecule_fasta, whitelist, tmp_path):
        """tenX_v2_split is a valid mode."""
        out = str(tmp_path / "out_v2")
        result = self._run_detect(single_molecule_fasta, whitelist, out, mode="tenX_v2_split")
        assert result.returncode == 0


class TestAssessBarcode:
    """Tests for assess_barcode_quality.py with tenX_v3_split mode."""

    def _run_assess(self, tsv, mode, output=None):
        cmd = ["python3", ASSESS_QUALITY, "--mode", mode, "--input", tsv]
        if output:
            cmd += ["--output", output]
        return subprocess.run(cmd, capture_output=True, text=True)

    def test_assess_tenx_split_runs(self, single_molecule_fasta, whitelist, tmp_path):
        # First generate TSV
        out = str(tmp_path / "out")
        subprocess.run([
            "python3", DETECT_BARCODES,
            "--input", single_molecule_fasta,
            "--mode", "tenX_v3_split",
            "--barcodes", whitelist,
            "-t", "1", "-o", out,
        ], check=True)
        tsv = out + ".barcoded_reads.tsv"
        result = self._run_assess(tsv, "tenX_v3_split")
        assert result.returncode == 0, result.stderr

    def test_assess_produces_metrics(self, single_molecule_fasta, whitelist, tmp_path):
        out = str(tmp_path / "out")
        subprocess.run([
            "python3", DETECT_BARCODES,
            "--input", single_molecule_fasta,
            "--mode", "tenX_v3_split",
            "--barcodes", whitelist,
            "-t", "1", "-o", out,
        ], check=True)
        tsv = out + ".barcoded_reads.tsv"
        report = str(tmp_path / "report.tsv")
        self._run_assess(tsv, "tenX_v3_split", output=report)
        assert os.path.exists(report)
        with open(report) as f:
            content = f.read()
        assert "total_reads" in content
        assert "assignment_precision" in content
