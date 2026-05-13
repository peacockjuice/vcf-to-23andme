"""Command-line entry points for vcf-to-23andme."""

from __future__ import annotations

import argparse
import logging
import sys


def convert_main(argv: list[str] | None = None) -> None:
    """CLI entry point: convert a VCF file to 23andMe v5 format."""
    parser = argparse.ArgumentParser(
        prog="vcf-to-23andme",
        description="Convert a VCF file to 23andMe v5 import format",
    )
    parser.add_argument("input", help="Path to input VCF file (.vcf or .vcf.gz)")
    parser.add_argument("output", help="Path for the output file")
    parser.add_argument(
        "--sample",
        help="Sample name to extract (defaults to first sample)",
    )
    parser.add_argument(
        "--include-custom-ids",
        action="store_true",
        help="Include marker IDs that are not rsIDs or 23andMe-style internal IDs",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable debug logging",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )

    try:
        from vcf_to_23andme.converter import convert_vcf_to_23andme

        count = convert_vcf_to_23andme(
            args.input,
            args.output,
            sample_name=args.sample,
            include_custom_ids=args.include_custom_ids,
        )
        print(f"Converted {count} variants to {args.output}")
    except (ValueError, FileNotFoundError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


def analyze_main(argv: list[str] | None = None) -> None:
    """CLI entry point: analyze a 23andMe-format file for known markers."""
    parser = argparse.ArgumentParser(
        prog="vcf-to-23andme-analyze",
        description="Analyze a genetic data file for known markers",
    )
    parser.add_argument(
        "input",
        help="Path to input genetic data file (23andMe format)",
    )
    parser.add_argument(
        "--markers",
        help="Path to a JSON file with custom markers (optional)",
        default=None,
    )
    parser.add_argument(
        "--output",
        help="Path for the analysis report (default: analysis_report.txt)",
        default="analysis_report.txt",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable debug logging",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )

    try:
        from vcf_to_23andme.analyzer import (
            analyze_dna_file,
            generate_report,
            load_known_markers,
        )

        known_markers = load_known_markers(args.markers)
        results = analyze_dna_file(args.input, known_markers)
        report = generate_report(results)

        with open(args.output, "w", encoding="utf-8") as f:
            f.write(report)

        print(f"Report saved to {args.output}")
    except (ValueError, FileNotFoundError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
