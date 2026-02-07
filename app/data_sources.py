from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import requests

from app.config import DEFAULT_FETCH_SIZE, DEFAULT_QUERY, SAMPLE_DATA_PATH, UNIPROT_ENDPOINT


@dataclass
class FetchResult:
    source: str
    bytes_written: int


def fetch_uniprot_data(output_path: Path, query: str = DEFAULT_QUERY, size: int = DEFAULT_FETCH_SIZE) -> FetchResult:
    """Fetch protein records from a public UniProt endpoint and save TSV locally."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    params = {
        "query": query,
        "format": "tsv",
        "fields": "accession,id,protein_name,gene_names,length,annotation_score",
        "size": size,
    }

    try:
        response = requests.get(UNIPROT_ENDPOINT, params=params, timeout=20)
        response.raise_for_status()
        text = response.text
        if "Entry" not in text:
            raise ValueError("Unexpected UniProt response format")
        output_path.write_text(text, encoding="utf-8")
        return FetchResult(source="uniprot", bytes_written=len(text.encode("utf-8")))
    except Exception:
        fallback_text = SAMPLE_DATA_PATH.read_text(encoding="utf-8")
        output_path.write_text(fallback_text, encoding="utf-8")
        return FetchResult(source="sample_fallback", bytes_written=len(fallback_text.encode("utf-8")))
