"""Tests for vcf_to_23andme.analyzer."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from vcf_to_23andme.analyzer import (
    analyze_dna_file,
    generate_report,
    load_known_markers,
)
from vcf_to_23andme.markers import KNOWN_MARKERS

FIXTURES = Path(__file__).parent / "fixtures"


# --- markers ---

def test_load_builtin_markers() -> None:
    markers = load_known_markers()
    assert markers == KNOWN_MARKERS


def test_load_custom_markers_json(tmp_path: Path) -> None:
    custom = {"rs999": {"gene": "TEST", "description": "test marker"}}
    json_file = tmp_path / "markers.json"
    json_file.write_text(json.dumps(custom), encoding="utf-8")
    markers = load_known_markers(str(json_file))
    assert markers == custom


# --- analyze ---

def test_analyze_file(tmp_path: Path) -> None:
    data_file = tmp_path / "data.txt"
    data_file.write_text(
        "# comment\n"
        "rs429358\t1\t12345\tCT\n"
        "rs7412\t1\t23456\tCC\n"
        "rsUnknown\t3\t99999\tAA\n",
        encoding="utf-8",
    )
    results = analyze_dna_file(str(data_file), KNOWN_MARKERS)
    assert "rs429358" in results
    assert "rs7412" in results
    assert "rsUnknown" not in results
    assert results["rs429358"]["genotype"] == "CT"
    assert results["rs429358"]["gene"] == "APOE"


# --- analyze: skip handling (parameterized) ---

@pytest.mark.parametrize(
    "content, expected_count",
    [
        ("# header\n# comment\nrs429358\t1\t12345\tCT\n", 1),
        ("rs429358\t1\n", 0),
    ],
    ids=["skips_comments", "skips_short_lines"],
)
def test_analyze_skips_invalid_lines(
    tmp_path: Path, content: str, expected_count: int,
) -> None:
    data_file = tmp_path / "data.txt"
    data_file.write_text(content, encoding="utf-8")
    results = analyze_dna_file(str(data_file), KNOWN_MARKERS)
    assert len(results) == expected_count


# --- NEW: empty data file ---

def test_analyze_empty_data_file(tmp_path: Path) -> None:
    """An empty file or file with only comments should return empty dict."""
    data_file = tmp_path / "empty.txt"
    data_file.write_text("# comment only\n", encoding="utf-8")
    results = analyze_dna_file(str(data_file), KNOWN_MARKERS)
    assert results == {}


# --- report ---

def test_report_format() -> None:
    results = {
        "rs429358": {
            "chromosome": "1",
            "position": "12345",
            "genotype": "CT",
            "gene": "APOE",
            "description": "Associated with Alzheimer's disease risk",
        }
    }
    report = generate_report(results)
    assert "Genetic Marker Analysis Report" in report
    assert "rs429358" in report
    assert "APOE" in report
    assert "CT" in report
    assert "Chromosome: 1" in report


def test_report_empty() -> None:
    report = generate_report({})
    assert "No known markers found" in report


def test_analyze_with_sample_23andme() -> None:
    """Integration: analyze the sample fixture file."""
    results = analyze_dna_file(str(FIXTURES / "sample_23andme.txt"), KNOWN_MARKERS)
    assert "rs429358" in results
    assert "rs7412" in results
