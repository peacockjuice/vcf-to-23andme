# vcf-to-23andme

VCF → 23andMe v5 converter + genetic marker analyzer. Python 3.9+, no dependencies.

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
