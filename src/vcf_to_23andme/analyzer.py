"""Analyze genetic marker data from 23andMe-format files."""

from __future__ import annotations

import json

KNOWN_MARKERS: dict[str, dict[str, str]] = {
    "rs429358": {
        "gene": "APOE",
        "description": "Associated with Alzheimer's disease risk",
    },
    "rs7412": {
        "gene": "APOE",
        "description": "APOE variant affecting lipid profile",
    },
    "rs1801133": {
        "gene": "MTHFR",
        "description": "Associated with folate metabolism",
    },
    "rs1799971": {
        "gene": "OPRM1",
        "description": "May affect pain sensitivity",
    },
    "rs671": {
        "gene": "ALDH2",
        "description": "Associated with alcohol metabolism",
    },
    "rs6265": {
        "gene": "BDNF",
        "description": "Affects neurotrophic activity",
    },
    "rs4680": {
        "gene": "COMT",
        "description": "Involved in dopamine metabolism",
    },
    "rs9939609": {
        "gene": "FTO",
        "description": "Associated with obesity predisposition",
    },
    "rs1805007": {
        "gene": "MC1R",
        "description": "Associated with pigmentation and hair color",
    },
    "rs12913832": {
        "gene": "HERC2",
        "description": "Associated with eye color",
    },
    "rs1229984": {
        "gene": "ADH1B",
        "description": "Affects alcohol metabolism",
    },
    "rs2228145": {
        "gene": "IL6R",
        "description": "Associated with inflammatory processes",
    },
    "rs7903146": {
        "gene": "TCF7L2",
        "description": "Associated with type 2 diabetes risk",
    },
    "rs1801131": {
        "gene": "MTHFR",
        "description": "Alternative variant associated with folate metabolism",
    },
    "rs662": {
        "gene": "PON1",
        "description": "Involved in antioxidant defense",
    },
    "rs1042713": {
        "gene": "ADRB2",
        "description": "Associated with stress response and receptor function",
    },
    "rs1800562": {
        "gene": "HFE",
        "description": "Associated with hemochromatosis",
    },
    "rs1800795": {
        "gene": "IL6",
        "description": "Associated with inflammatory processes",
    },
}


def load_known_markers(markers_file: str | None = None) -> dict[str, dict[str, str]]:
    """Load known genetic markers from a JSON file or return the built-in database."""
    if markers_file:
        with open(markers_file, "r", encoding="utf-8") as f:
            return json.load(f)
    return dict(KNOWN_MARKERS)


def analyze_dna_file(
    input_file: str,
    known_markers: dict[str, dict[str, str]],
) -> dict[str, dict[str, str]]:
    """Analyze a genetic data file and extract information for known markers."""
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
