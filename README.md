# vcf-to-23andme

Convert VCF genomic files to 23andMe v5 format and analyze genetic markers.

## Features

- Convert VCF files (plain or gzip-compressed) to 23andMe v5 import format
- Analyze converted files against 18 built-in genetic markers
- Support for multi-sample VCF files
- Custom marker databases via JSON
- Zero external dependencies (Python 3.9+ stdlib only)

## Installation

```bash
pip install git+https://github.com/peacockjuice/vcf-to-23andme.git
```

Or install locally for development:

```bash
git clone https://github.com/peacockjuice/vcf-to-23andme.git
cd vcf-to-23andme
pip install -e .
```

## Usage

### Convert VCF to 23andMe format

```bash
vcf-to-23andme input.vcf output.txt
```

With a specific sample from a multi-sample VCF:

```bash
vcf-to-23andme input.vcf output.txt --sample SampleName
```

### Analyze genetic markers

```bash
vcf-to-23andme-analyze output.txt --output report.txt
```

With a custom marker database:

```bash
vcf-to-23andme-analyze output.txt --markers my_markers.json --output report.txt
```

Both commands support `-v` / `--verbose` for debug logging.

### Custom markers JSON format

```json
{
  "rs429358": {
    "gene": "APOE",
    "description": "Associated with Alzheimer's disease risk"
  },
  "rs7412": {
    "gene": "APOE",
    "description": "APOE variant affecting lipid profile"
  }
}
```

## Built-in markers

18 markers are included by default (APOE, MTHFR, COMT, FTO, and others). See [`src/vcf_to_23andme/markers.py`](src/vcf_to_23andme/markers.py) for the full list.

## Development

```bash
pip install -e .
pip install pytest
python -m pytest tests/ -v
```

## License

MIT
