# compbiologist

Minimal starter for a virtual computational biologist with:
- Public data ingestion (UniProt API)
- Local GEO metadata search (GEOmetadb SQLite)
- Local data processing (pandas)
- Web UI for exploration (Flask)

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
# Place GEOmetadb.sqlite at project root for local metadata search.
# Optional but recommended for faster/higher-rate NCBI calls:
export NCBI_API_KEY="your_ncbi_api_key"
export NCBI_EMAIL="your_email@example.com"
python -m app.main
```

Open [http://127.0.0.1:8000](http://127.0.0.1:8000).

## Run with Docker

1. Ensure `GEOmetadb.sqlite` is present at project root.
2. Optional: set NCBI credentials for higher API limits.
3. Build and run:

```bash
docker compose up -d --build
```

Then open [http://127.0.0.1:8000](http://127.0.0.1:8000).

The compose setup mounts:
- `./data` -> `/app/data` (persistent app cache/results)
- `./GEOmetadb.sqlite` -> `/app/GEOmetadb.sqlite` (read-only metadata database)

Stop with:

```bash
docker compose down
```

## Deploy on a server (shared access)

On a Linux VM/server with Docker + Compose:

```bash
git clone <your-repo-url>
cd compbiologist

# Put GEOmetadb.sqlite at repository root.
# Optional credentials for higher NCBI throughput:
export NCBI_API_KEY="your_ncbi_api_key"
export NCBI_EMAIL="your_email@example.com"

docker compose up -d --build
```

If you run a firewall:
- allow inbound TCP `8000`, or
- put a reverse proxy (Caddy/Nginx) in front and expose only `80/443`.

Minimum practical disk recommendation:
- `>= 40 GB` (the GEO SQLite file alone is large, plus cache and matrix downloads).

## Web UI workflows

- `Refresh Protein Data`: fetches UniProt records and re-runs local protein processing.
- GEO search is at the top of the page and supports filters:
  - species text filter
  - experiment-type filter (`Single-cell RNA-seq`, `RNA-seq`, `Microarray`, etc.)
  - `Only analyzable datasets` pre-check (attempts to keep only GEO entries with parseable expression matrices)
  - non-analyzable rows are not selectable for loading
- After search, select datasets and click `Load Selected Datasets`.
- GEO detail tables/charts are shown only for loaded datasets.
- Run `State Comparison Analysis` to generate:
  - volcano plot
  - top up/down genes
  - pathway enrichment (up and down)
  - real GEO matrix differential expression using DESeq2 (`pydeseq2`, counts-only)
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

- Primary public sources: UniProt REST API and GEOmetadb SQLite (metadata), plus NCBI GEO FTP for matrix files.
- GEO search responses and analyzability checks are cached locally in `data/processed/geo_cache.json`.
- If local SQLite metadata search fails, the UI shows `Source: geo_sqlite_error`.
