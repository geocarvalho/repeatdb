# repeatdb

Analysis of **repeat expansions** from VCFs produced by [ExpansionHunter](https://github.com/Illumina/ExpansionHunter): batch ingestion, per-locus outlier detection, and reports (TSV, JSON, or HTML).

## Requirements

- **Python 3.10 or newer**
- VCFs in the usual ExpansionHunter layout (`FORMAT` fields such as `REPCN`, `INFO` with `REPID`, etc.)

On many systems (macOS with Homebrew, some Linux distros) the system Python is **externally managed (PEP 668)** and does not allow global `pip install`. Always use a **virtual environment**:

```bash
python3 -m venv .venv
source .venv/bin/activate   # Linux / macOS
# .venv\Scripts\activate    # Windows
```

## Installation

From the repository root (with the venv activated):

```bash
pip install -e .
```

This installs the `repeatdb` package and the **`repeatdb`** command on the venv’s PATH.

To run tests with **pytest**:

```bash
pip install -e ".[dev]"
pytest
```

## Quick start (CLI)

Show help:

```bash
repeatdb --help
repeatdb ingest --help
repeatdb detect --help
repeatdb query --help
```

### 1. Ingest VCFs → Parquet

Create a text file with **one absolute or relative VCF path per line** (avoid blank lines in the middle; they are skipped):

```text
/path/to/sample1.vcf
/path/to/sample2.vcf
```

Then:

```bash
repeatdb ingest --vcf-list vcfs.txt --output calls.parquet --workers 4
```

The result is a Parquet file with columns `sample_id`, `locus_id`, `allele1`, `allele2`, `filter`, `genotype`.

### 2. Detect outliers → report

On the calls Parquet:

```bash
repeatdb detect --input calls.parquet --output outliers.tsv --method mad --threshold 3.5 --format tsv
```

- **`--method`**: `mad` (default), `iqr`, or `zscore`. For `iqr`, a typical threshold is `1.5` (Tukey fences); for `zscore`, values around `2.5`–`3.0` are common, depending on the data.
- **`--format`**: `tsv` (default), `json`, or `html` (interactive table with DataTables via CDN).

### 3. Query calls by locus

The **`--gene`** filter matches the Parquet **`locus_id`** column (REPID-style identifier, e.g. `HTT`).

```bash
repeatdb query --input calls.parquet --gene HTT
repeatdb query --input calls.parquet --gene HTT --sample SAMPLE1
```

Output is a formatted table in the terminal (Rich).

## Library usage

You can import modules (`repeatdb.vcf_parser`, `repeatdb.store`, `repeatdb.aggregator`, `repeatdb.outlier`, `repeatdb.reporter`) in your own pipelines; the CLI above is the most direct end-to-end path.

## Package layout

| Module | Role |
|--------|------|
| `schema.py` | Pydantic models (`Locus`, `RepeatCall`, `Sample`) |
| `vcf_parser.py` | VCF parsing (cyvcf2) |
| `store.py` | Parallel ingestion and Parquet I/O |
| `aggregator.py` | Per-locus aggregation (max allele) |
| `outlier.py` | Outlier detection |
| `reporter.py` | TSV / JSON / HTML |
| `cli.py` | Typer CLI (`repeatdb`) |
