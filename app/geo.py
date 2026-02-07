from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import requests

from app.config import (
    GEO_DEFAULT_FETCH_SIZE,
    GEO_DEFAULT_QUERY,
    GEO_ESEARCH_ENDPOINT,
    GEO_ESUMMARY_ENDPOINT,
    GEO_SAMPLE_PATH,
)

EXPERIMENT_TYPE_OPTIONS = [
    "All",
    "Single-cell RNA-seq",
    "RNA-seq",
    "Microarray",
    "ChIP-seq",
    "ATAC-seq",
    "Other",
]


@dataclass
class GEOSearchResult:
    source: str
    query: str
    total_found: int
    items: list[dict]


def _build_geo_link(accession: str) -> str:
    if not accession:
        return ""
    return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"


def _build_pubmed_link(pubmed_id: str) -> str:
    if not pubmed_id:
        return ""
    return f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"


def _infer_experiment_type(title: str, summary: str, gds_type: str) -> str:
    text = f"{title} {summary} {gds_type}".lower()
    if any(token in text for token in ["single cell", "single-cell", "scrna", "snrna"]):
        return "Single-cell RNA-seq"
    if any(token in text for token in ["rna-seq", "rnaseq", "high throughput sequencing", "transcriptome sequencing"]):
        return "RNA-seq"
    if any(token in text for token in ["microarray", "expression profiling by array"]):
        return "Microarray"
    if "chip-seq" in text:
        return "ChIP-seq"
    if "atac-seq" in text:
        return "ATAC-seq"
    return "Other"


def _parse_geo_summary_payload(payload: dict, id_list: list[str]) -> list[dict]:
    result = payload.get("result", {})
    items: list[dict] = []

    for uid in id_list:
        row = result.get(uid, {})
        title = row.get("title", "")
        summary = row.get("summary", "")
        gds_type = row.get("gdstype", "")
        items.append(
            {
                "uid": uid,
                "accession": row.get("accession", ""),
                "title": title,
                "summary": summary,
                "organism": row.get("taxon", ""),
                "n_samples": row.get("n_samples", ""),
                "gse": row.get("gse", ""),
                "pubmed_id": row.get("pubmed_id", ""),
                "pdat": row.get("pdat", ""),
                "gds_type": gds_type,
                "experiment_type": _infer_experiment_type(title, summary, gds_type),
            }
        )

    return items


def enrich_geo_items(items: list[dict]) -> list[dict]:
    enriched: list[dict] = []
    for row in items:
        row_copy = dict(row)
        accession = str(row_copy.get("accession", "")).strip()
        pubmed_id = str(row_copy.get("pubmed_id", "")).strip()
        row_copy["geo_link"] = _build_geo_link(accession)
        row_copy["pubmed_link"] = _build_pubmed_link(pubmed_id)
        if not row_copy.get("experiment_type"):
            row_copy["experiment_type"] = _infer_experiment_type(
                str(row_copy.get("title", "")),
                str(row_copy.get("summary", "")),
                str(row_copy.get("gds_type", "")),
            )
        enriched.append(row_copy)
    return enriched


def filter_geo_items(items: list[dict], species_filter: str = "", experiment_filter: str = "All") -> list[dict]:
    species_filter = species_filter.strip().lower()
    experiment_filter = experiment_filter.strip()

    filtered = []
    for row in items:
        organism = str(row.get("organism", "")).lower()
        experiment_type = str(row.get("experiment_type", "Other"))

        species_ok = not species_filter or species_filter in organism
        exp_ok = experiment_filter in ("", "All") or experiment_type == experiment_filter

        if species_ok and exp_ok:
            filtered.append(row)

    return filtered


def build_geo_insights(items: list[dict]) -> dict:
    if not items:
        return {
            "organism_distribution": {},
            "sample_distribution": {},
            "experiment_distribution": {},
        }

    df = pd.DataFrame(items)
    if "organism" not in df.columns:
        df["organism"] = ""
    if "n_samples" not in df.columns:
        df["n_samples"] = ""
    if "experiment_type" not in df.columns:
        df["experiment_type"] = "Other"

    organism_distribution = (
        df["organism"]
        .fillna("Unknown")
        .replace("", "Unknown")
        .value_counts()
        .head(8)
        .to_dict()
    )

    experiment_distribution = (
        df["experiment_type"]
        .fillna("Other")
        .replace("", "Other")
        .value_counts()
        .to_dict()
    )

    sample_bins = [0, 10, 50, 100, 500, 1000, 5000, 100000]
    sample_labels = ["1-10", "11-50", "51-100", "101-500", "501-1000", "1001-5000", "5001+"]
    n_samples_numeric = pd.to_numeric(df["n_samples"], errors="coerce")
    binned = pd.cut(n_samples_numeric, bins=sample_bins, labels=sample_labels, include_lowest=True)
    sample_distribution = (
        binned.value_counts(dropna=False)
        .reindex(sample_labels, fill_value=0)
        .astype(int)
        .to_dict()
    )

    return {
        "organism_distribution": organism_distribution,
        "sample_distribution": sample_distribution,
        "experiment_distribution": experiment_distribution,
    }


def search_geo_datasets(query: str = GEO_DEFAULT_QUERY, retmax: int = GEO_DEFAULT_FETCH_SIZE) -> GEOSearchResult:
    """Search GEO datasets (GDS) via NCBI E-utilities and return normalized records."""
    esearch_params = {
        "db": "gds",
        "term": query,
        "retmode": "json",
        "retmax": retmax,
        "sort": "relevance",
    }

    try:
        esearch_resp = requests.get(GEO_ESEARCH_ENDPOINT, params=esearch_params, timeout=20)
        esearch_resp.raise_for_status()
        esearch_data = esearch_resp.json()

        search_result = esearch_data.get("esearchresult", {})
        id_list = search_result.get("idlist", [])
        total_found = int(search_result.get("count", 0))

        if not id_list:
            return GEOSearchResult(source="geo", query=query, total_found=total_found, items=[])

        esummary_params = {
            "db": "gds",
            "id": ",".join(id_list),
            "retmode": "json",
        }
        esummary_resp = requests.get(GEO_ESUMMARY_ENDPOINT, params=esummary_params, timeout=20)
        esummary_resp.raise_for_status()
        esummary_data = esummary_resp.json()

        items = _parse_geo_summary_payload(esummary_data, id_list)
        return GEOSearchResult(source="geo", query=query, total_found=total_found, items=items)
    except Exception:
        sample = json.loads(GEO_SAMPLE_PATH.read_text(encoding="utf-8"))
        items = enrich_geo_items(sample.get("items", []))
        return GEOSearchResult(
            source="sample_fallback",
            query=query,
            total_found=sample.get("total_found", len(items)),
            items=items,
        )


def write_geo_artifacts(raw_path: Path, processed_path: Path, summary_path: Path, result: GEOSearchResult, species_filter: str = "", experiment_filter: str = "All") -> dict:
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    processed_path.parent.mkdir(parents=True, exist_ok=True)

    raw_payload = {
        "query": result.query,
        "source": result.source,
        "total_found": result.total_found,
        "species_filter": species_filter,
        "experiment_filter": experiment_filter,
        "items": result.items,
    }
    raw_path.write_text(json.dumps(raw_payload, indent=2), encoding="utf-8")

    df = pd.DataFrame(result.items)
    if not df.empty:
        columns = [
            "accession",
            "title",
            "organism",
            "experiment_type",
            "n_samples",
            "gse",
            "pubmed_id",
            "pdat",
            "summary",
        ]
        for col in columns:
            if col not in df.columns:
                df[col] = ""
        df = df[columns]
    df.to_csv(processed_path, index=False)

    summary = {
        "query": result.query,
        "source": result.source,
        "returned": len(result.items),
        "total_found": result.total_found,
        "species_filter": species_filter,
        "experiment_filter": experiment_filter,
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    return summary


def save_loaded_geo_selection(loaded_json_path: Path, loaded_csv_path: Path, source_payload: dict, selected_ids: list[str]) -> dict:
    loaded_json_path.parent.mkdir(parents=True, exist_ok=True)
    loaded_csv_path.parent.mkdir(parents=True, exist_ok=True)

    selected = set(selected_ids)
    items = source_payload.get("items", [])

    loaded_items = []
    for row in items:
        identifier = str(row.get("accession") or row.get("uid") or "")
        if identifier and identifier in selected:
            loaded_items.append(row)

    loaded_payload = {
        "query": source_payload.get("query", ""),
        "source": source_payload.get("source", "not_fetched"),
        "species_filter": source_payload.get("species_filter", ""),
        "experiment_filter": source_payload.get("experiment_filter", "All"),
        "returned": len(loaded_items),
        "selected_ids": sorted(selected),
        "items": loaded_items,
    }
    loaded_json_path.write_text(json.dumps(loaded_payload, indent=2), encoding="utf-8")

    df = pd.DataFrame(loaded_items)
    if not df.empty:
        columns = [
            "accession",
            "title",
            "organism",
            "experiment_type",
            "n_samples",
            "gse",
            "pubmed_id",
            "pdat",
            "summary",
        ]
        for col in columns:
            if col not in df.columns:
                df[col] = ""
        df = df[columns]
    df.to_csv(loaded_csv_path, index=False)

    return loaded_payload


def load_cached_geo(raw_path: Path) -> dict | None:
    if not raw_path.exists():
        return None
    return json.loads(raw_path.read_text(encoding="utf-8"))


def load_cached_loaded_geo(loaded_json_path: Path) -> dict | None:
    if not loaded_json_path.exists():
        return None
    return json.loads(loaded_json_path.read_text(encoding="utf-8"))
