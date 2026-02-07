# compbiologist

Minimal starter for a virtual computational biologist with:
- Public data ingestion (UniProt API)
- Public GEO dataset search (NCBI E-utilities)
- Local data processing (pandas)
- Web UI for exploration (Flask)

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m app.main
```

Open [http://127.0.0.1:8000](http://127.0.0.1:8000).

## Web UI workflows

- `Refresh Protein Data`: fetches UniProt records and re-runs local protein processing.
- GEO search is at the top of the page and supports filters:
  - species text filter
  - experiment-type filter (`Single-cell RNA-seq`, `RNA-seq`, `Microarray`, etc.)
  - `Only analyzable datasets` pre-check (attempts to keep only GEO entries with parseable expression matrices)
- After search, select datasets and click `Load Selected Datasets`.
- GEO detail tables/charts are shown only for loaded datasets.
- Run `State Comparison Analysis` to generate:
  - volcano plot
  - top up/down genes
  - pathway enrichment (up and down)
  - real GEO matrix differential expression (when a selected GSE has usable series-matrix + sample-group labels)
- if automatic grouping fails, the UI now requires manual sample-group selection before real DE can run
  - current implementation analyzes one loaded GSE at a time (manual group selection supported when auto-detection fails)

## Pipeline only (CLI)

```bash
python -m app.main --refresh
```

With optional GEO fetch in the same run:

```bash
python -m app.main --refresh --geo-query "single cell glioblastoma" --geo-retmax 20
```

## Developer commands

```bash
make venv
make install
make check
make test
make run
```

## Local artifacts

- Protein raw TSV: `data/raw/uniprot_kinase.tsv`
- Protein processed CSV: `data/processed/processed_proteins.csv`
- Protein summary JSON: `data/processed/summary.json`
- GEO raw JSON: `data/raw/geo_last_search.json`
- GEO processed CSV: `data/processed/geo_datasets.csv`
- GEO summary JSON: `data/processed/geo_summary.json`

## Notes

- Primary public sources: UniProt REST API and NCBI E-utilities (GEO).
- If either API is unavailable, the app uses bundled sample data so the UI still works.
