# compbiologist

Virtual computational biology web app with:
- GEO metadata search from local `GEOmetadb.sqlite`
- analyzability checks and selective dataset loading
- real differential expression (DESeq2 via `pydeseq2`) on GEO count-like matrices
- interactive volcano plot and pathway enrichment
- gene/pathway external links (Ensembl, GO/Reactome/WikiPathways/KEGG)
- Docker-ready deployment

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
On first load, protein summary artifacts are initialized automatically if missing.

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

Workflow is explicit and step-based:

1. `Step 1: Run Search`
- Filters: query, species (top species + Other), experiment type, state filter, count.
- Optional `Only analyzable datasets` to pre-filter to datasets likely usable for real DE.
- GEO IDs are linked to GEO and PubMed when available.

2. `Step 2: Run Load Selected`
- Select analyzable rows from Step 1 table.
- Loaded list is revalidated and restricted to analyzable datasets.

3. `Step 3: Run Analysis`
- Controls: comparison state, `padj` cutoff, `|log2FC|` cutoff, enrichment mode (`auto`, `gprofiler`, `local`).
- Outputs: interactive volcano plot, top up/down genes, enrichment up/down tables.
- Top genes include external Ensembl links and best-effort gene annotations.

Manual grouping mode:
- If automatic sample grouping fails, UI switches to manual group assignment.
- In this mode, the main Step 3 run button is suppressed and a single run button appears in the manual grouping card.
- Real analysis currently runs on one selected GSE at a time.

Workflow reset:
- `Start From Scratch` clears cached GEO search/load/analysis artifacts.

## Pipeline only (CLI)

```bash
python -m app.main --refresh
```

With optional GEO fetch in the same run:

```bash
python -m app.main --refresh --geo-query "single cell glioblastoma" --geo-retmax 20
```

Run web server on a different interface/port:

```bash
python -m app.main --host 0.0.0.0 --port 8000
```

## Developer commands

```bash
make venv
make install
make check
make test
make run
make docker-up
make docker-down
make docker-logs
```

## Local artifacts

- Protein raw TSV: `data/raw/uniprot_kinase.tsv`
- Protein processed CSV: `data/processed/processed_proteins.csv`
- Protein summary JSON: `data/processed/summary.json`
- GEO raw JSON: `data/raw/geo_last_search.json`
- GEO processed CSV: `data/processed/geo_datasets.csv`
- GEO summary JSON: `data/processed/geo_summary.json`
- Loaded GEO selection JSON/CSV: `data/processed/geo_loaded.json`, `data/processed/geo_loaded.csv`
- Analysis result JSON/CSV: `data/processed/analysis_result.json`, `data/processed/analysis_dge.csv`
- Gene annotation cache: `data/processed/gene_info_cache.json`

## Notes

- Primary public sources: GEOmetadb SQLite metadata, NCBI GEO FTP matrices, UniProt REST API, optional g:Profiler, optional Ensembl REST lookup.
- GEO search responses and analyzability checks are cached in `data/processed/geo_cache.json`.
- If local SQLite metadata search fails, source label can show `geo_sqlite_error` or fallback source labels.
- In restricted/offline environments, Ensembl gene names may be unavailable; IDs and external links remain available.
