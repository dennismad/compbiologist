from __future__ import annotations

import re
from urllib.parse import quote_plus
from typing import Any

import pandas as pd

GP_DEFAULT_SOURCES = ["GO:BP", "REAC", "WP", "KEGG"]


def infer_gprofiler_organism(loaded_items: list[dict]) -> str:
    text = " ".join(str(item.get("organism", "")) for item in loaded_items).lower()
    if "homo sapiens" in text or "human" in text:
        return "hsapiens"
    if "mus musculus" in text or "mouse" in text:
        return "mmusculus"
    if "rattus norvegicus" in text or "rat" in text:
        return "rnorvegicus"
    return "hsapiens"


def _normalize_identifier(raw: str) -> str:
    token = str(raw).strip()
    if not token:
        return ""

    if re.match(r"^ENS[A-Z0-9]+\.\d+$", token, flags=re.IGNORECASE):
        token = token.split(".", 1)[0]

    if re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]-\d+$", token, flags=re.IGNORECASE):
        token = token.split("-", 1)[0]

    return token.strip()


def _prepare_gene_list(gene_ids: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for raw in gene_ids:
        norm = _normalize_identifier(raw)
        if norm and norm not in seen:
            seen.add(norm)
            out.append(norm)
    return out


def _load_gprofiler() -> tuple[Any | None, str | None]:
    try:
        from gprofiler import GProfiler  # type: ignore
    except Exception as exc:
        return None, f"gprofiler-official unavailable: {exc}"
    return GProfiler, None


def _build_pathway_url(source: str, native: str, name: str) -> str:
    src = str(source).strip().upper()
    term = str(native).strip()
    if src.startswith("GO") and term.startswith("GO:"):
        return f"https://www.ebi.ac.uk/QuickGO/term/{term}"
    if src == "REAC" and term:
        return f"https://reactome.org/content/detail/{term}"
    if src == "KEGG" and term:
        return f"https://www.kegg.jp/entry/{term}"
    if src == "WP" and term:
        return f"https://www.wikipathways.org/pathways/{term}.html"
    query = quote_plus(f"{source} {native} {name}".strip())
    return f"https://reactome.org/content/query?q={query}"


def run_gprofiler_enrichment(
    gene_ids: list[str],
    organism: str,
    max_terms: int = 10,
    user_threshold: float = 0.05,
) -> tuple[list[dict], str | None]:
    query_genes = _prepare_gene_list(gene_ids)
    if len(query_genes) < 2:
        return [], None

    GProfiler, import_error = _load_gprofiler()
    if GProfiler is None:
        return [], import_error

    try:
        gp = GProfiler(return_dataframe=True)
        res = gp.profile(
            organism=organism,
            query=query_genes,
            sources=GP_DEFAULT_SOURCES,
            user_threshold=min(max(float(user_threshold), 1e-6), 1.0),
            no_iea=False,
        )
    except Exception as exc:
        return [], f"g:Profiler enrichment failed: {exc}"

    if res is None:
        return [], None
    if not isinstance(res, pd.DataFrame):
        res = pd.DataFrame(res)
    if res.empty:
        return [], None

    if "p_value" in res.columns:
        res = res.sort_values("p_value", ascending=True)

    rows: list[dict] = []
    for _, row in res.head(max_terms).iterrows():
        source = str(row.get("source", "")).strip()
        native = str(row.get("native", "")).strip()
        name = str(row.get("name", "")).strip()
        pathway = f"{source}:{native} {name}".strip()
        overlap = int(pd.to_numeric(row.get("intersection_size", 0), errors="coerce") or 0)
        set_size = int(pd.to_numeric(row.get("term_size", 0), errors="coerce") or 0)
        padj = float(pd.to_numeric(row.get("p_value", 1.0), errors="coerce") or 1.0)
        rows.append(
            {
                "pathway": pathway,
                "overlap": overlap,
                "set_size": set_size,
                "padj": padj,
                "pathway_url": _build_pathway_url(source, native, name),
            }
        )

    return rows, None
