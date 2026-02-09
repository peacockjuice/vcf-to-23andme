"""Convert VCF files to 23andMe v5 import format."""

from __future__ import annotations

import gzip
from typing import IO


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
    """Convert a VCF file to 23andMe v5 import format."""
    count = 0
    header_found = False

    with open_vcf(input_file) as fin, open(output_file, "w", encoding="utf-8") as fout:
        fout.write("# This file was generated for 23andMe v5 import\n")
        fout.write("# rsid\tchromosome\tposition\tgenotype\n")

        sample_index: int | None = None

        for line in fin:
            line = line.strip()

            if line.startswith("##"):
                continue

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

            parts = line.split("\t")
            if len(parts) < 10:
                continue

            chrom = parts[0]
            pos = parts[1]
            rsid = parts[2]
            ref = parts[3]
            alt = parts[4]

            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue

            format_fields = parts[8].split(":")
            try:
                gt_index = format_fields.index("GT")
            except ValueError:
                continue

            sample_data = parts[sample_index].split(":")
            if len(sample_data) <= gt_index:
                continue
            gt_field = sample_data[gt_index]

            if gt_field in (".", "./.", ".|."):
                continue

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
