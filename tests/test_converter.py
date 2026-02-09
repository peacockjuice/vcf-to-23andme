"""Tests for vcf_to_23andme.converter."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from vcf_to_23andme.converter import convert_vcf_to_23andme, open_vcf

FIXTURES = Path(__file__).parent / "fixtures"


# --- open_vcf (parameterized) ---

@pytest.mark.parametrize("use_gzip", [False, True], ids=["plain", "gzip"])
def test_open_vcf(tmp_path: Path, use_gzip: bool) -> None:
    if use_gzip:
        vcf = tmp_path / "test.vcf.gz"
        with gzip.open(str(vcf), "wt", encoding="utf-8") as f:
            f.write("hello\n")
    else:
        vcf = tmp_path / "test.vcf"
        vcf.write_text("hello\n", encoding="utf-8")
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


# --- skips (parameterized) ---

@pytest.mark.parametrize(
    "rsid",
    ["rs_indel", "rs_missing", "rs_multi", "rs_no_gt"],
    ids=["indels", "missing_genotype", "multi_allelic", "no_gt_field"],
)
def test_skip_variants(tmp_path: Path, rsid: str) -> None:
    """Non-SNP and missing-genotype variants should not appear in output."""
    output = tmp_path / "out.txt"
    convert_vcf_to_23andme(str(FIXTURES / "sample.vcf"), str(output))
    content = output.read_text(encoding="utf-8")
    assert rsid not in content


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


# --- NEW: allele index > 1 maps to "N" ---

def test_allele_index_beyond_one_maps_to_N(tmp_path: Path) -> None:
    """When an allele index is not 0 or 1, it should be mapped to 'N'."""
    vcf = tmp_path / "allele_n.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
        "1\t100\trs_test\tA\tG\t50\tPASS\tDB\tGT\t2/0\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    convert_vcf_to_23andme(str(vcf), str(output))
    content = output.read_text(encoding="utf-8")
    assert "rs_test\t1\t100\tNA" in content


# --- NEW: phased genotype ---

def test_phased_genotype_explicit(tmp_path: Path) -> None:
    """Phased genotype separator '|' should be handled like '/'."""
    vcf = tmp_path / "phased.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
        "1\t200\trs_phased_test\tA\tG\t50\tPASS\tDB\tGT\t1|0\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    convert_vcf_to_23andme(str(vcf), str(output))
    content = output.read_text(encoding="utf-8")
    assert "rs_phased_test\t1\t200\tGA" in content


# --- NEW: homozygous REF and ALT ---

def test_homozygous_ref_and_alt(tmp_path: Path) -> None:
    """Homozygous genotypes 0/0 and 1/1 should produce doubled REF/ALT."""
    vcf = tmp_path / "homo.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
        "1\t100\trs_homo_ref\tC\tT\t50\tPASS\tDB\tGT\t0/0\n"
        "1\t200\trs_homo_alt\tC\tT\t50\tPASS\tDB\tGT\t1/1\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    convert_vcf_to_23andme(str(vcf), str(output))
    content = output.read_text(encoding="utf-8")
    assert "rs_homo_ref\t1\t100\tCC" in content
    assert "rs_homo_alt\t1\t200\tTT" in content


# --- NEW: header with too few columns ---

def test_header_too_few_columns(tmp_path: Path) -> None:
    """A #CHROM header with fewer than 10 columns should raise ValueError."""
    vcf = tmp_path / "short_header.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    with pytest.raises(ValueError, match="no sample data"):
        convert_vcf_to_23andme(str(vcf), str(output))


# --- NEW: VCF with header only, no data lines ---

def test_vcf_header_only_no_data(tmp_path: Path) -> None:
    """A VCF with a valid header but no data lines should return count 0."""
    vcf = tmp_path / "header_only.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    count = convert_vcf_to_23andme(str(vcf), str(output))
    assert count == 0
