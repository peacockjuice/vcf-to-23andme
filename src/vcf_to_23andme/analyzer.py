"""Analyze genetic marker data from 23andMe-format files."""

from __future__ import annotations

import json
import logging

from vcf_to_23andme.markers import KNOWN_MARKERS

log = logging.getLogger(__name__)


def load_known_markers(markers_file: str | None = None) -> dict[str, dict[str, str]]:
    """Load known genetic markers.

    If *markers_file* is given, load markers from a JSON file.  Otherwise
    return the built-in marker database.
    """
    if markers_file:
        with open(markers_file, "r", encoding="utf-8") as f:
            return json.load(f)
    return dict(KNOWN_MARKERS)


def analyze_dna_file(
    input_file: str,
    known_markers: dict[str, dict[str, str]],
) -> dict[str, dict[str, str]]:
    """Analyze a genetic data file and extract information for known markers.

    The input file is expected to be tab-separated with columns:
    rsid, chromosome, position, genotype.

    Returns:
        Dict mapping rsid to marker details (chromosome, position, genotype,
        gene, description).
    """
    results: dict[str, dict[str, str]] = {}
    with open(input_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            rsid, chrom, pos, genotype = parts[0], parts[1], parts[2], parts[3]
            if rsid in known_markers:
                results[rsid] = {
                    "chromosome": chrom,
                    "position": pos,
                    "genotype": genotype,
                    "gene": known_markers[rsid]["gene"],
                    "description": known_markers[rsid]["description"],
                }
    return results


def generate_report(analysis_results: dict[str, dict[str, str]]) -> str:
    """Generate a plain-text report for matched markers."""
    lines: list[str] = [
        "Genetic Marker Analysis Report",
        "=" * 40,
    ]
    if not analysis_results:
        lines.append("No known markers found in the input file.")
    else:
        for rsid, data in analysis_results.items():
            lines.append(
                f"{rsid} ({data['gene']}): {data['genotype']} "
                f"[Chromosome: {data['chromosome']}, Position: {data['position']}]"
            )
            lines.append(f"  Description: {data['description']}")
            lines.append("")
    return "\n".join(lines)
