# vcf-to-23andme

VCF → 23andMe-like raw genotype text converter + genetic marker analyzer.
Python 3.9+, no dependencies.

The converter produces a 4-column text file shaped like 23andMe raw data
(`rsid`, `chromosome`, `position`, `genotype`). It is not an official 23andMe
export and does not add 23andMe-specific validation metadata.

The analyzer is a marker lookup helper, not medical interpretation. It does not
perform clinical risk modeling, strand/build validation, or genotype-specific
medical guidance.

## Quick start

```bash
pip install git+https://github.com/peacockjuice/vcf-to-23andme.git

# convert
vcf-to-23andme input.vcf output.txt

# analyze against 18 built-in markers (APOE, MTHFR, COMT, FTO, etc.)
vcf-to-23andme-analyze output.txt --output report.txt
```

## Options

```bash
# multi-sample VCF — pick a specific sample
vcf-to-23andme input.vcf output.txt --sample SampleName

# include custom marker IDs that are not rsIDs or 23andMe-style internal IDs
vcf-to-23andme input.vcf output.txt --include-custom-ids

# custom marker database
vcf-to-23andme-analyze output.txt --markers my_markers.json --output report.txt
```

Custom markers JSON:
```json
{
  "rs429358": { "gene": "APOE", "description": "Alzheimer's risk" }
}
```

## Development

```bash
pip install -e . && pip install pytest
pytest tests/ -v
```

## License

MIT
