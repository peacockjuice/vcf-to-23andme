"""Convert VCF files to 23andMe v5 import format."""

from __future__ import annotations

import gzip
import logging
from typing import IO

log = logging.getLogger(__name__)


def open_vcf(filename: str) -> IO[str]:
    """Open a VCF file. Supports both plain .vcf and compressed .vcf.gz."""
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt", encoding="utf-8")
    return open(filename, "r", encoding="utf-8")


def convert_vcf_to_23andme(
    input_file: str,
    output_file: str,
    sample_name: str | None = None,
) -> int:
    """Convert a VCF file to 23andMe v5 import format.

    Only SNPs where REF and ALT are single characters are processed.

    Args:
        input_file: Path to the input VCF file (.vcf or .vcf.gz).
        output_file: Path for the output file.
        sample_name: Optional sample name when the VCF contains multiple
            samples.  If not given, the first sample (column index 9) is used.

    Returns:
        Number of variants written.

    Raises:
        ValueError: If the ``#CHROM`` header line is missing or the
            requested sample is not found.
    """
    count = 0
    header_found = False

    with open_vcf(input_file) as fin, open(output_file, "w", encoding="utf-8") as fout:
        # Write 23andMe v5 header
        fout.write("# This file was generated for 23andMe v5 import\n")
        fout.write("# rsid\tchromosome\tposition\tgenotype\n")

        sample_index: int | None = None

        for line in fin:
            line = line.strip()

            # Skip VCF metadata lines
            if line.startswith("##"):
                continue

            # Parse the column header line
            if line.startswith("#CHROM"):
                header_found = True
                columns = line.split("\t")
                if len(columns) < 10:
                    raise ValueError(
                        "VCF file has no sample data (need at least 10 columns)."
                    )

                if sample_name:
                    if sample_name in columns:
                        sample_index = columns.index(sample_name)
                    else:
                        raise ValueError(
                            f"Sample '{sample_name}' not found in VCF header."
                        )
                else:
                    sample_index = 9
                continue

            if not header_found:
                continue

            # Split variant line into columns
            parts = line.split("\t")
            if len(parts) < 10:
                log.debug("Skipping malformed line: %s", line[:80])
                continue

            chrom = parts[0]
            pos = parts[1]
            rsid = parts[2]
            ref = parts[3]
            alt = parts[4]

            # Only process SNPs: skip indels and multi-allelic sites
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                log.debug("Skipping non-SNP variant: %s", rsid)
                continue

            # Find GT index in the FORMAT column
            format_fields = parts[8].split(":")
            try:
                gt_index = format_fields.index("GT")
            except ValueError:
                log.debug("No GT field for variant %s, skipping", rsid)
                continue

            # Extract the sample genotype
            if sample_index is None:
                raise ValueError("sample_index was not set; this should not happen")
            sample_data = parts[sample_index].split(":")
            if len(sample_data) <= gt_index:
                continue
            gt_field = sample_data[gt_index]

            # Skip missing genotypes
            if gt_field in (".", "./.", ".|."):
                log.debug("Missing genotype for %s, skipping", rsid)
                continue

            # Decode genotype: split on '/' or '|', substitute REF/ALT
            alleles: list[str] = []
            for allele in gt_field.replace("|", "/").split("/"):
                if allele == "0":
                    alleles.append(ref)
                elif allele == "1":
                    alleles.append(alt)
                else:
                    alleles.append("N")
            genotype = "".join(alleles)

            fout.write(f"{rsid}\t{chrom}\t{pos}\t{genotype}\n")
            count += 1

    if not header_found:
        raise ValueError("VCF file is missing the #CHROM header line.")

    return count
