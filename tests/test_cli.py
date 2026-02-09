"""Integration tests for CLI entry points."""

from __future__ import annotations

from pathlib import Path

import pytest

from vcf_to_23andme.cli import analyze_main, convert_main

FIXTURES = Path(__file__).parent / "fixtures"


def test_convert_cli(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    output = tmp_path / "out.txt"
    convert_main([str(FIXTURES / "sample.vcf"), str(output)])
    captured = capsys.readouterr()
    assert "Converted 4 variants" in captured.out
    assert output.exists()


def test_convert_cli_sample(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    output = tmp_path / "out.txt"
    convert_main([str(FIXTURES / "sample.vcf"), str(output), "--sample", "SampleB"])
    captured = capsys.readouterr()
    assert "Converted" in captured.out
    # Verify SampleB's genotype was actually used
    content = output.read_text(encoding="utf-8")
    assert "rs429358\t1\t12345\tTT" in content


def test_convert_cli_error(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as exc_info:
        convert_main([str(tmp_path / "nope.vcf"), str(tmp_path / "out.txt")])
    assert exc_info.value.code == 1


def test_analyze_cli(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    report = tmp_path / "report.txt"
    analyze_main([
        str(FIXTURES / "sample_23andme.txt"),
        "--output", str(report),
    ])
    captured = capsys.readouterr()
    assert "Report saved to" in captured.out
    assert report.exists()
    content = report.read_text(encoding="utf-8")
    assert "Genetic Marker Analysis Report" in content


def test_analyze_cli_error(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as exc_info:
        analyze_main([str(tmp_path / "nope.txt")])
    assert exc_info.value.code == 1
