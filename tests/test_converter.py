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
    assert count == 5

    expected = (FIXTURES / "sample_23andme.txt").read_text(encoding="utf-8")
    actual = output.read_text(encoding="utf-8")
    assert actual == expected


# --- skips (parameterized) ---

@pytest.mark.parametrize(
    "rsid",
    ["rs_indel", "rs_missing", "rs_no_gt"],
    ids=["indels", "missing_genotype", "no_gt_field"],
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
    assert count == 5


# --- NEW: unsupported genotypes are skipped ---

@pytest.mark.parametrize("gt", ["2/0", "0/.", "./1", "bad"])
def test_unsupported_or_partial_genotype_is_skipped(tmp_path: Path, gt: str) -> None:
    """Unsupported or partial GT values should not be emitted as fake data."""
    vcf = tmp_path / "allele_n.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
        f"1\t100\trs_test\tA\tG\t50\tPASS\tDB\tGT\t{gt}\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    count = convert_vcf_to_23andme(str(vcf), str(output))
    content = output.read_text(encoding="utf-8")
    assert count == 0
    assert "rs_test" not in content


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


def test_chromosome_names_are_normalized(tmp_path: Path) -> None:
    """Common chr-prefixed chromosome names should match 23andMe-like output."""
    vcf = tmp_path / "chroms.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
        "chr1\t100\trs_chr1\tA\tG\t50\tPASS\tDB\tGT\t0/1\n"
        "chrX\t200\trs_chrx\tC\tT\t50\tPASS\tDB\tGT\t1/1\n"
        "chrY\t300\trs_chry\tG\tA\t50\tPASS\tDB\tGT\t0/0\n"
        "chrM\t400\trs_chrm\tT\tC\t50\tPASS\tDB\tGT\t0/1\n"
        "M\t500\trs_m\tA\tC\t50\tPASS\tDB\tGT\t1/1\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    convert_vcf_to_23andme(str(vcf), str(output))
    content = output.read_text(encoding="utf-8")
    assert "rs_chr1\t1\t100\tAG" in content
    assert "rs_chrx\tX\t200\tTT" in content
    assert "rs_chry\tY\t300\tGG" in content
    assert "rs_chrm\tMT\t400\tTC" in content
    assert "rs_m\tMT\t500\tCC" in content


def test_marker_id_filtering(tmp_path: Path) -> None:
    """Missing IDs are skipped, while rsIDs and internal i IDs are retained."""
    vcf = tmp_path / "ids.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
        "1\t100\t.\tA\tG\t50\tPASS\tDB\tGT\t0/1\n"
        "1\t200\trs_ok\tA\tG\t50\tPASS\tDB\tGT\t0/1\n"
        "1\t300\ti12345\tA\tG\t50\tPASS\tDB\tGT\t0/1\n"
        "1\t400\tcustom_marker\tA\tG\t50\tPASS\tDB\tGT\t0/1\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    count = convert_vcf_to_23andme(str(vcf), str(output))
    content = output.read_text(encoding="utf-8")
    assert count == 2
    assert "\t.\t" not in content
    assert "rs_ok\t1\t200\tAG" in content
    assert "i12345\t1\t300\tAG" in content
    assert "custom_marker" not in content


def test_custom_marker_ids_can_be_included(tmp_path: Path) -> None:
    """Custom marker IDs should be opt-in."""
    vcf = tmp_path / "custom_ids.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
        "1\t100\tcustom_marker\tA\tG\t50\tPASS\tDB\tGT\t0/1\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    count = convert_vcf_to_23andme(
        str(vcf), str(output), include_custom_ids=True,
    )
    content = output.read_text(encoding="utf-8")
    assert count == 1
    assert "custom_marker\t1\t100\tAG" in content


def test_multi_allelic_snp_selected_alleles_convert(tmp_path: Path) -> None:
    """Multi-allelic SNPs can be converted when selected alleles are simple bases."""
    vcf = tmp_path / "multi.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
        "1\t100\trs_multi_test\tA\tC,G\t50\tPASS\tDB\tGT\t2/0\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    count = convert_vcf_to_23andme(str(vcf), str(output))
    content = output.read_text(encoding="utf-8")
    assert count == 1
    assert "rs_multi_test\t1\t100\tGA" in content


def test_multi_allelic_indel_selected_allele_is_skipped(tmp_path: Path) -> None:
    """Multi-allelic records are skipped when the selected allele is not a SNP."""
    vcf = tmp_path / "multi_indel.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"
        "1\t100\trs_multi_indel\tA\tC,GA\t50\tPASS\tDB\tGT\t2/0\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    count = convert_vcf_to_23andme(str(vcf), str(output))
    content = output.read_text(encoding="utf-8")
    assert count == 0
    assert "rs_multi_indel" not in content


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
