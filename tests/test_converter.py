"""Tests for vcf_to_23andme.converter."""

from __future__ import annotations

import gzip
import os
import tempfile
from pathlib import Path

import pytest

from vcf_to_23andme.converter import convert_vcf_to_23andme, open_vcf

FIXTURES = Path(__file__).parent / "fixtures"


# --- open_vcf ---

def test_open_vcf_plain(tmp_path: Path) -> None:
    vcf = tmp_path / "test.vcf"
    vcf.write_text("hello\n", encoding="utf-8")
    with open_vcf(str(vcf)) as f:
        assert f.read() == "hello\n"


def test_open_vcf_gzip(tmp_path: Path) -> None:
    vcf = tmp_path / "test.vcf.gz"
    with gzip.open(str(vcf), "wt", encoding="utf-8") as f:
        f.write("hello\n")
    with open_vcf(str(vcf)) as f:
        assert f.read() == "hello\n"


# --- convert basic ---

def test_convert_basic(tmp_path: Path) -> None:
    output = tmp_path / "out.txt"
    count = convert_vcf_to_23andme(str(FIXTURES / "sample.vcf"), str(output))
    assert count == 4

    expected = (FIXTURES / "sample_23andme.txt").read_text(encoding="utf-8")
    actual = output.read_text(encoding="utf-8")
    assert actual == expected


def test_convert_returns_count(tmp_path: Path) -> None:
    output = tmp_path / "out.txt"
    count = convert_vcf_to_23andme(str(FIXTURES / "sample.vcf"), str(output))
    assert isinstance(count, int)
    assert count == 4


# --- skips ---

def test_skip_indels(tmp_path: Path) -> None:
    """Indels (AT->A) should not appear in output."""
    output = tmp_path / "out.txt"
    convert_vcf_to_23andme(str(FIXTURES / "sample.vcf"), str(output))
    content = output.read_text(encoding="utf-8")
    assert "rs_indel" not in content


def test_skip_missing_genotype(tmp_path: Path) -> None:
    """Variants with ./. genotype should not appear."""
    output = tmp_path / "out.txt"
    convert_vcf_to_23andme(str(FIXTURES / "sample.vcf"), str(output))
    content = output.read_text(encoding="utf-8")
    assert "rs_missing" not in content


def test_skip_multi_allelic(tmp_path: Path) -> None:
    """Multi-allelic sites (comma in ALT) should be skipped."""
    output = tmp_path / "out.txt"
    convert_vcf_to_23andme(str(FIXTURES / "sample.vcf"), str(output))
    content = output.read_text(encoding="utf-8")
    assert "rs_multi" not in content


def test_skip_no_gt_field(tmp_path: Path) -> None:
    """Variants without a GT field in FORMAT should be skipped."""
    output = tmp_path / "out.txt"
    convert_vcf_to_23andme(str(FIXTURES / "sample.vcf"), str(output))
    content = output.read_text(encoding="utf-8")
    assert "rs_no_gt" not in content


# --- multi-sample ---

def test_multi_sample_selection(tmp_path: Path) -> None:
    output = tmp_path / "out.txt"
    count = convert_vcf_to_23andme(
        str(FIXTURES / "sample.vcf"), str(output), sample_name="SampleB",
    )
    content = output.read_text(encoding="utf-8")
    # SampleB has rs429358 = 1/1 -> TT
    assert "rs429358\t1\t12345\tTT" in content
    assert count > 0


def test_sample_not_found(tmp_path: Path) -> None:
    output = tmp_path / "out.txt"
    with pytest.raises(ValueError, match="not found"):
        convert_vcf_to_23andme(
            str(FIXTURES / "sample.vcf"), str(output), sample_name="NoSuch",
        )


# --- error cases ---

def test_missing_chrom_header(tmp_path: Path) -> None:
    vcf = tmp_path / "bad.vcf"
    vcf.write_text("##fileformat=VCFv4.1\n1\t100\trs1\tA\tT\t.\t.\t.\tGT\t0/1\n")
    output = tmp_path / "out.txt"
    with pytest.raises(ValueError, match="#CHROM"):
        convert_vcf_to_23andme(str(vcf), str(output))


def test_file_not_found(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        convert_vcf_to_23andme(str(tmp_path / "nope.vcf"), str(tmp_path / "out.txt"))


# --- gzip round-trip ---

def test_convert_gzip(tmp_path: Path) -> None:
    gz_path = tmp_path / "sample.vcf.gz"
    plain = (FIXTURES / "sample.vcf").read_text(encoding="utf-8")
    with gzip.open(str(gz_path), "wt", encoding="utf-8") as f:
        f.write(plain)

    output = tmp_path / "out.txt"
    count = convert_vcf_to_23andme(str(gz_path), str(output))
    assert count == 4
