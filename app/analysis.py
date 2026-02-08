from __future__ import annotations

import hashlib
import json
import math
import os
import random
from pathlib import Path
from urllib.parse import quote_plus
import re

import pandas as pd
import requests

from app.config import GENE_INFO_CACHE_PATH, GEO_MATRIX_CACHE_DIR
from app.enrichment import infer_gprofiler_organism, run_gprofiler_enrichment
from app.geo_expression import build_group_choice_context, run_real_geo_de

ENRICHMENT_MODE_OPTIONS = ["auto", "gprofiler", "local"]
DEFAULT_ENRICHMENT_MODE = "auto"

GENES = [
    "IL6", "TNF", "CXCL8", "STAT3", "NFKB1", "NFKBIA", "IFNG", "JAK2", "SOCS3", "CCL2",
    "PPARG", "SLC2A4", "INSR", "AKT1", "MTOR", "MAPK1", "MAPK3", "HIF1A", "VEGFA", "KDR",
    "MKI67", "PCNA", "CDK1", "CDK2", "CCND1", "BCL2", "BAX", "CASP3", "TP53", "MYC",
    "EPCAM", "KRT8", "KRT18", "VIM", "CDH1", "CDH2", "MMP2", "MMP9", "COL1A1", "COL3A1",
    "PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "GZMB", "PRF1", "FOXP3", "CD3D", "CD8A",
    "HLA-A", "HLA-B", "HLA-C", "B2M", "IRF1", "IRF7", "ISG15", "MX1", "OAS1", "IFI44",
    "EGFR", "ERBB2", "KRAS", "NRAS", "BRAF", "PIK3CA", "PTEN", "RB1", "SMAD4", "APC",
    "GAPDH", "ACTB", "RPLP0", "RPS6", "EEF1A1", "TUBA1B", "TUBB", "LDHA", "PKM", "ENO1",
]

PATHWAYS = {
    "Inflammatory_Response": {"IL6", "TNF", "CXCL8", "CCL2", "NFKB1", "NFKBIA", "STAT3", "IFNG"},
    "Interferon_Signaling": {"IFNG", "IRF1", "IRF7", "ISG15", "MX1", "OAS1", "IFI44"},
    "Cell_Cycle_Proliferation": {"MKI67", "PCNA", "CDK1", "CDK2", "CCND1", "MYC"},
    "Apoptosis": {"BCL2", "BAX", "CASP3", "TP53"},
    "PI3K_AKT_MTOR": {"PIK3CA", "PTEN", "AKT1", "MTOR", "INSR", "SLC2A4"},
    "EMT_Extracellular_Matrix": {"VIM", "CDH2", "CDH1", "MMP2", "MMP9", "COL1A1", "COL3A1"},
    "Immune_Checkpoint": {"PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "CD3D", "CD8A"},
    "Hypoxia_Angiogenesis": {"HIF1A", "VEGFA", "KDR", "LDHA", "ENO1"},
}

LOCAL_PATHWAY_URLS = {
    "Inflammatory_Response": "https://www.ebi.ac.uk/QuickGO/term/GO:0006954",
    "Interferon_Signaling": "https://reactome.org/content/query?q=Interferon%20Signaling",
    "Cell_Cycle_Proliferation": "https://reactome.org/content/query?q=Cell%20Cycle",
    "Apoptosis": "https://reactome.org/content/query?q=Apoptosis",
    "PI3K_AKT_MTOR": "https://reactome.org/content/query?q=PI3K%20AKT%20mTOR%20signaling",
    "EMT_Extracellular_Matrix": "https://reactome.org/content/query?q=Extracellular%20matrix%20organization",
    "Immune_Checkpoint": "https://reactome.org/content/query?q=PD-1%20signaling",
    "Hypoxia_Angiogenesis": "https://reactome.org/content/query?q=HIF-1%20signaling",
}

STATE_SIGNAL_GENES = {
    "Disease vs Healthy": ["IL6", "TNF", "NFKB1", "STAT3", "MKI67", "MYC", "HIF1A", "VEGFA"],
    "Disease only": ["IL6", "TNF", "CCL2", "NFKBIA", "MKI67", "PCNA"],
    "Healthy only": ["SLC2A4", "PPARG", "INSR", "CDH1"],
    "Treated vs Untreated": ["IFNG", "IRF1", "ISG15", "MX1", "PDCD1", "CD274"],
}


def _seed_from_context(state_profile: str, dataset_ids: list[str]) -> int:
    key = f"{state_profile}|{'|'.join(sorted(dataset_ids))}"
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
    return int(digest[:8], 16)


def _bh_adjust(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    qvals = [1.0] * n
    min_coeff = 1.0
    for rank in range(n, 0, -1):
        idx, p = indexed[rank - 1]
        coeff = min(min_coeff, p * n / rank)
        min_coeff = coeff
        qvals[idx] = min(1.0, coeff)
    return qvals


def _sanitize_enrichment_mode(mode: str | None) -> str:
    value = str(mode or DEFAULT_ENRICHMENT_MODE).strip().lower()
    if value not in ENRICHMENT_MODE_OPTIONS:
        return DEFAULT_ENRICHMENT_MODE
    return value


def _ensembl_gene_url(gene_id: str) -> str:
    gene = str(gene_id).strip()
    if not gene:
        return ""
    if re.match(r"^ENS[A-Z0-9]+(?:\.\d+)?$", gene, flags=re.IGNORECASE):
        gene = gene.split(".", 1)[0]
        return f"https://www.ensembl.org/id/{quote_plus(gene)}"
    return f"https://www.ensembl.org/Multi/Search/Results?q={quote_plus(gene)};site=ensembl"


def _normalize_ensembl_gene_id(gene_id: str) -> str:
    gene = str(gene_id).strip()
    if not gene:
        return ""
    if re.match(r"^ENS[A-Z0-9]+(?:\.\d+)?$", gene, flags=re.IGNORECASE):
        return gene.split(".", 1)[0].upper()
    return ""


def _clean_gene_description(raw: str) -> str:
    text = str(raw or "").strip()
    if not text:
        return ""
    return re.sub(r"\s*\[Source:.*?\]\s*$", "", text).strip()


def _load_gene_info_cache(cache_path: Path | None = None) -> dict[str, dict]:
    path = cache_path or GENE_INFO_CACHE_PATH
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    entries = payload.get("entries", {})
    if isinstance(entries, dict):
        return entries
    return {}


def _write_gene_info_cache(entries: dict[str, dict], cache_path: Path | None = None) -> None:
    path = cache_path or GENE_INFO_CACHE_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"entries": entries}
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _fetch_ensembl_gene_annotations(ensembl_ids: list[str], timeout_s: float = 10.0) -> dict[str, dict]:
    if not ensembl_ids:
        return {}
    endpoint = "https://rest.ensembl.org/lookup/id"
    try:
        response = requests.post(
            endpoint,
            headers={"Accept": "application/json", "Content-Type": "application/json"},
            json={"ids": ensembl_ids, "expand": False},
            timeout=timeout_s,
        )
        response.raise_for_status()
        payload = response.json()
    except Exception:
        return {}

    if not isinstance(payload, dict):
        return {}

    out: dict[str, dict] = {}
    for raw_id in ensembl_ids:
        rec = payload.get(raw_id)
        if not isinstance(rec, dict):
            continue
        resolved_id = _normalize_ensembl_gene_id(str(rec.get("id", "")).strip()) or raw_id
        symbol = str(rec.get("display_name", "")).strip() or resolved_id
        out[raw_id] = {
            "ensembl_id": resolved_id,
            "gene_symbol": symbol,
            "gene_name": _clean_gene_description(rec.get("description", "")),
            "gene_biotype": str(rec.get("biotype", "")).strip(),
        }
    return out


def _annotate_top_gene_rows(rows: list[dict], allow_remote_lookup: bool) -> tuple[list[dict], str | None]:
    out_rows = [dict(row) for row in rows]
    ens_ids = sorted(
        {
            _normalize_ensembl_gene_id(str(row.get("gene", "")))
            for row in out_rows
            if _normalize_ensembl_gene_id(str(row.get("gene", "")))
        }
    )

    cache_entries = _load_gene_info_cache()
    missing = [gene_id for gene_id in ens_ids if gene_id not in cache_entries]

    allow_remote_in_tests = os.getenv("COMPBIO_FORCE_REMOTE_GENE_LOOKUP", "").strip() == "1"
    if allow_remote_lookup and missing and (allow_remote_in_tests or not os.getenv("PYTEST_CURRENT_TEST")):
        fetched = _fetch_ensembl_gene_annotations(missing)
        for key in missing:
            cache_entries[key] = fetched.get(
                key,
                {
                    "ensembl_id": key,
                    "gene_symbol": key,
                    "gene_name": "",
                    "gene_biotype": "",
                },
            )
        _write_gene_info_cache(cache_entries)

    missing_annotations = 0
    for row in out_rows:
        gene_raw = str(row.get("gene", "")).strip()
        ens_id = _normalize_ensembl_gene_id(gene_raw)
        cached = cache_entries.get(ens_id, {}) if ens_id else {}
        gene_symbol = str(cached.get("gene_symbol", "")).strip() if cached else ""
        gene_name = str(cached.get("gene_name", "")).strip() if cached else ""
        gene_biotype = str(cached.get("gene_biotype", "")).strip() if cached else ""

        if ens_id and not cached:
            missing_annotations += 1

        row["gene_symbol"] = gene_symbol or gene_raw
        row["gene_name"] = gene_name
        row["gene_biotype"] = gene_biotype
        if not row.get("gene_url"):
            row["gene_url"] = _ensembl_gene_url(cached.get("ensembl_id", gene_raw) if ens_id else gene_raw)

    note: str | None = None
    if allow_remote_lookup and ens_ids and missing_annotations:
        note = f"Gene annotation lookup unavailable for {missing_annotations} Ensembl IDs."
    return out_rows, note


def _normalize_pathway_identifier(pathway: str) -> str:
    text = str(pathway).strip()
    if not text:
        return ""
    parts = text.split(" ", 1)[0]
    return parts


def _pathway_search_url(pathway_name: str) -> str:
    pathway = str(pathway_name).strip()
    if not pathway:
        return ""

    if pathway in LOCAL_PATHWAY_URLS:
        return LOCAL_PATHWAY_URLS[pathway]

    go_match = re.search(r"\bGO:\d{7}\b", pathway, flags=re.IGNORECASE)
    if go_match:
        return f"https://www.ebi.ac.uk/QuickGO/term/{go_match.group(0).upper()}"

    reactome_match = re.search(r"\bR-[A-Z]{3}-\d+\b", pathway, flags=re.IGNORECASE)
    if reactome_match:
        return f"https://reactome.org/content/detail/{reactome_match.group(0).upper()}"

    wp_match = re.search(r"\bWP\d+\b", pathway, flags=re.IGNORECASE)
    if wp_match:
        return f"https://www.wikipathways.org/pathways/{wp_match.group(0).upper()}.html"

    token = _normalize_pathway_identifier(pathway)
    if token.upper().startswith("KEGG:"):
        kegg_id = token.split(":", 1)[1].strip()
        if kegg_id:
            return f"https://www.kegg.jp/entry/{quote_plus(kegg_id)}"

    query = quote_plus(pathway)
    return f"https://reactome.org/content/query?q={query}"


def _hydrate_top_gene_links(rows: list[dict]) -> list[dict]:
    out: list[dict] = []
    for row in rows:
        item = dict(row)
        if not item.get("gene_url"):
            item["gene_url"] = _ensembl_gene_url(str(item.get("gene", "")))
        out.append(item)
    return out


def _hydrate_pathway_links(rows: list[dict]) -> list[dict]:
    out: list[dict] = []
    for row in rows:
        item = dict(row)
        if not item.get("pathway_url"):
            item["pathway_url"] = _pathway_search_url(str(item.get("pathway", "")))
        out.append(item)
    return out


def hydrate_analysis_links(payload: dict | None) -> dict | None:
    if payload is None:
        return None
    out = dict(payload)
    top_up = _hydrate_top_gene_links(list(out.get("top_up", [])))
    top_down = _hydrate_top_gene_links(list(out.get("top_down", [])))
    top_up, _ = _annotate_top_gene_rows(top_up, allow_remote_lookup=True)
    top_down, _ = _annotate_top_gene_rows(top_down, allow_remote_lookup=True)
    out["top_up"] = top_up
    out["top_down"] = top_down
    out["enrichment_up"] = _hydrate_pathway_links(list(out.get("enrichment_up", [])))
    out["enrichment_down"] = _hydrate_pathway_links(list(out.get("enrichment_down", [])))
    return out


def generate_differential_expression(state_profile: str, dataset_ids: list[str]) -> pd.DataFrame:
    rng = random.Random(_seed_from_context(state_profile, dataset_ids))
    target = set(STATE_SIGNAL_GENES.get(state_profile, []))

    rows = []
    for gene in GENES:
        baseline = rng.gauss(0.0, 0.8)
        if gene in target:
            baseline += rng.choice([1.2, 1.5, 1.8])
        if state_profile == "Healthy only" and gene in {"IL6", "TNF", "MKI67"}:
            baseline -= rng.choice([1.0, 1.4])

        if gene in target:
            p = max(1e-12, min(1.0, math.exp(-abs(baseline) * rng.uniform(2.8, 4.2)) * rng.uniform(1e-5, 0.01)))
        else:
            p = max(1e-8, min(1.0, math.exp(-abs(baseline) * rng.uniform(1.0, 2.2)) * rng.uniform(0.01, 0.6)))
        rows.append({"gene": gene, "log2fc": baseline, "pvalue": p})

    df = pd.DataFrame(rows)
    df["padj"] = _bh_adjust(df["pvalue"].tolist())
    df["neg_log10_padj"] = -df["padj"].clip(lower=1e-300).map(math.log10)
    return df


def _hypergeom_sf(k_minus_1: int, m: int, n: int, N: int) -> float:
    max_k = min(m, N)
    denom = math.comb(m + n, N)
    prob = 0.0
    for k in range(k_minus_1 + 1, max_k + 1):
        prob += (math.comb(m, k) * math.comb(n, N - k)) / denom
    return min(1.0, max(0.0, prob))


def run_enrichment(
    df: pd.DataFrame,
    direction: str,
    pathways: dict[str, set[str]],
    padj_cutoff: float,
    log2fc_cutoff: float,
    organism: str,
    enrichment_mode: str,
) -> tuple[list[dict], str | None]:
    if direction == "up":
        sig = df[(df["padj"] <= padj_cutoff) & (df["log2fc"] >= log2fc_cutoff)]
    else:
        sig = df[(df["padj"] <= padj_cutoff) & (df["log2fc"] <= -log2fc_cutoff)]

    sig_genes = set(sig["gene"].tolist())
    if not sig_genes:
        return [], None

    mode = _sanitize_enrichment_mode(enrichment_mode)
    module_note: str | None = None
    if mode in {"auto", "gprofiler"}:
        module_rows, module_note = run_gprofiler_enrichment(
            gene_ids=list(sig_genes),
            organism=organism,
            max_terms=10,
            user_threshold=padj_cutoff,
        )
        if module_rows:
            for row in module_rows:
                if not row.get("pathway_url"):
                    row["pathway_url"] = _pathway_search_url(str(row.get("pathway", "")))
            return module_rows, None
        if mode == "gprofiler":
            return [], module_note

    universe = set(df["gene"].tolist())
    M = len(universe)
    N = len(sig_genes)

    rows = []
    for pathway, genes in pathways.items():
        gset = genes & universe
        if not gset:
            continue
        overlap = sig_genes & gset
        if not overlap:
            continue
        m = len(gset)
        n = M - m
        pval = _hypergeom_sf(len(overlap) - 1, m, n, N)
        rows.append(
            {
                "pathway": pathway,
                "overlap": len(overlap),
                "set_size": m,
                "gene_ratio": round(len(overlap) / m, 3),
                "pvalue": pval,
                "pathway_url": _pathway_search_url(pathway),
            }
        )

    if not rows:
        return [], module_note

    pvals = [r["pvalue"] for r in rows]
    qvals = _bh_adjust(pvals)
    for i, q in enumerate(qvals):
        rows[i]["padj"] = q

    rows = sorted(rows, key=lambda x: x["padj"])
    return rows[:10], module_note


def build_volcano_points(df: pd.DataFrame, padj_cutoff: float, log2fc_cutoff: float) -> tuple[list[dict], float, float]:
    if df.empty:
        return [], 1.5, 2.0

    max_abs_fc = max(1.5, float(df["log2fc"].abs().max()))
    max_neglog = max(2.0, float(df["neg_log10_padj"].max()))

    points = []
    for row in df.to_dict(orient="records"):
        x = float(row["log2fc"])
        y = float(row["neg_log10_padj"])
        is_sig = row["padj"] <= padj_cutoff and abs(row["log2fc"]) >= log2fc_cutoff
        if not is_sig:
            color = "#8c9a91"
            group = "Not Significant"
        elif row["log2fc"] > 0:
            color = "#c0392b"
            group = "Upregulated"
        else:
            color = "#1f78b4"
            group = "Downregulated"

        px = 420 + (x / max_abs_fc) * 360
        py = 300 - (y / max_neglog) * 250
        points.append(
            {
                "gene": row["gene"],
                "x": x,
                "y": y,
                "px": round(px, 2),
                "py": round(py, 2),
                "color": color,
                "group": group,
            }
        )

    return points, max_abs_fc, max_neglog


def _empty_analysis_payload(
    state_profile: str,
    dataset_ids: list[str],
    padj_cutoff: float,
    log2fc_cutoff: float,
    enrichment_mode: str,
    mode: str,
    notes: list[str],
    group_choice: dict | None = None,
) -> dict:
    return {
        "metadata": {
            "state_profile": state_profile,
            "padj_cutoff": padj_cutoff,
            "log2fc_cutoff": log2fc_cutoff,
            "n_datasets": len(dataset_ids),
            "dataset_ids": dataset_ids,
            "n_genes": 0,
            "n_up": 0,
            "n_down": 0,
            "enrichment_mode": enrichment_mode,
            "mode": mode,
            "real_dataset_used": "",
            "group_a_name": "",
            "group_b_name": "",
            "group_a_n": 0,
            "group_b_n": 0,
            "notes": notes,
        },
        "group_choice": group_choice,
        "volcano": {"points": [], "max_abs_fc": 1.5, "max_neglog": 2.0},
        "top_up": [],
        "top_down": [],
        "enrichment_up": [],
        "enrichment_down": [],
        "dge_rows": [],
    }


def _dedupe_messages(messages: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for msg in messages:
        text = str(msg).strip()
        if text and text not in seen:
            seen.add(text)
            out.append(text)
    return out


def _has_missing_dependency_error(messages: list[str]) -> bool:
    return any("pydeseq2 is required" in str(msg).lower() for msg in messages)


def _merge_manual_choice_into_group_context(choice_context: dict | None, manual_choice: dict | None) -> dict | None:
    if not choice_context or not manual_choice:
        return choice_context

    manual_gse = str(manual_choice.get("gse", "")).strip()
    context_gse = str(choice_context.get("gse", "")).strip()
    if manual_gse and context_gse and manual_gse != context_gse:
        return choice_context

    sample_options = choice_context.get("sample_options", [])
    valid_ids = {str(x.get("id", "")).strip() for x in sample_options}

    group_a = [str(x).strip() for x in manual_choice.get("group_a_samples", []) if str(x).strip() in valid_ids]
    group_b = [str(x).strip() for x in manual_choice.get("group_b_samples", []) if str(x).strip() in valid_ids and str(x).strip() not in set(group_a)]

    merged = dict(choice_context)
    merged["suggested_group_a"] = group_a
    merged["suggested_group_b"] = group_b
    if str(manual_choice.get("group_a_name", "")).strip():
        merged["group_a_name"] = str(manual_choice.get("group_a_name")).strip()
    if str(manual_choice.get("group_b_name", "")).strip():
        merged["group_b_name"] = str(manual_choice.get("group_b_name")).strip()
    return merged


def run_state_comparison_analysis(
    state_profile: str,
    loaded_items: list[dict],
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    enrichment_mode: str = DEFAULT_ENRICHMENT_MODE,
    use_real_geo: bool = True,
    manual_choice: dict | None = None,
) -> dict:
    dataset_ids = [str(x.get("accession") or x.get("gse") or x.get("uid") or "") for x in loaded_items]
    dataset_ids = [x for x in dataset_ids if x]
    enrichment_mode = _sanitize_enrichment_mode(enrichment_mode)

    if not use_real_geo:
        df = generate_differential_expression(state_profile=state_profile, dataset_ids=dataset_ids)
        mode = "prototype_signature_bootstrap"
        notes: list[str] = []
        real_gse = ""
        group_a_name = "Case"
        group_b_name = "Control"
        group_a_n = 0
        group_b_n = 0
    else:
        real_result, errors = run_real_geo_de(
            loaded_items=loaded_items,
            state_profile=state_profile,
            cache_dir=GEO_MATRIX_CACHE_DIR,
            manual_choice=manual_choice,
        )
        if real_result is None:
            choice_context, choice_errors = build_group_choice_context(
                loaded_items=loaded_items,
                state_profile=state_profile,
                cache_dir=GEO_MATRIX_CACHE_DIR,
            )
            choice_context = _merge_manual_choice_into_group_context(choice_context, manual_choice)
            notes = _dedupe_messages(list(errors) + list(choice_errors))
            if _has_missing_dependency_error(notes):
                notes.append("Install pydeseq2 and rerun analysis to enable real DE.")
                mode = "missing_dependency"
                choice_context = None
            elif choice_context:
                notes.append("Automatic grouping failed. Choose sample groups manually and rerun real analysis.")
                mode = "needs_group_selection"
            else:
                notes.append("No usable expression matrix/sample metadata found for manual grouping in selected datasets.")
                mode = "no_usable_matrix_data"
            return _empty_analysis_payload(
                state_profile=state_profile,
                dataset_ids=dataset_ids,
                padj_cutoff=padj_cutoff,
                log2fc_cutoff=log2fc_cutoff,
                enrichment_mode=enrichment_mode,
                mode=mode,
                notes=notes,
                group_choice=choice_context,
            )

        df = real_result.dataframe.copy()
        mode = "geo_series_matrix_manual" if manual_choice else "geo_series_matrix"
        notes = list(errors)
        real_gse = real_result.gse
        group_a_name = real_result.group_a_name
        group_b_name = real_result.group_b_name
        group_a_n = len(real_result.group_a_samples)
        group_b_n = len(real_result.group_b_samples)

    df["direction"] = "neutral"
    df.loc[(df["padj"] <= padj_cutoff) & (df["log2fc"] >= log2fc_cutoff), "direction"] = "up"
    df.loc[(df["padj"] <= padj_cutoff) & (df["log2fc"] <= -log2fc_cutoff), "direction"] = "down"
    df["gene_url"] = df["gene"].map(_ensembl_gene_url)

    n_up = int((df["direction"] == "up").sum())
    n_down = int((df["direction"] == "down").sum())

    top_up = (
        df[df["direction"] == "up"]
        .sort_values(["padj", "log2fc"], ascending=[True, False])
        .head(10)
        [["gene", "gene_url", "log2fc", "padj"]]
        .to_dict(orient="records")
    )
    top_down = (
        df[df["direction"] == "down"]
        .sort_values(["padj", "log2fc"], ascending=[True, True])
        .head(10)
        [["gene", "gene_url", "log2fc", "padj"]]
        .to_dict(orient="records")
    )
    top_up, gene_info_note_up = _annotate_top_gene_rows(top_up, allow_remote_lookup=True)
    top_down, gene_info_note_down = _annotate_top_gene_rows(top_down, allow_remote_lookup=True)
    if gene_info_note_up:
        notes.append(gene_info_note_up)
    if gene_info_note_down:
        notes.append(gene_info_note_down)

    organism = infer_gprofiler_organism(loaded_items)
    enrichment_up, enrich_note_up = run_enrichment(
        df,
        direction="up",
        pathways=PATHWAYS,
        padj_cutoff=padj_cutoff,
        log2fc_cutoff=log2fc_cutoff,
        organism=organism,
        enrichment_mode=enrichment_mode,
    )
    enrichment_down, enrich_note_down = run_enrichment(
        df,
        direction="down",
        pathways=PATHWAYS,
        padj_cutoff=padj_cutoff,
        log2fc_cutoff=log2fc_cutoff,
        organism=organism,
        enrichment_mode=enrichment_mode,
    )
    if enrich_note_up:
        notes.append(enrich_note_up)
    if enrich_note_down:
        notes.append(enrich_note_down)
    notes = _dedupe_messages(notes)

    volcano_points, max_abs_fc, max_neglog = build_volcano_points(df, padj_cutoff=padj_cutoff, log2fc_cutoff=log2fc_cutoff)

    return {
        "metadata": {
            "state_profile": state_profile,
            "padj_cutoff": padj_cutoff,
            "log2fc_cutoff": log2fc_cutoff,
            "n_datasets": len(dataset_ids),
            "dataset_ids": dataset_ids,
            "n_genes": int(df.shape[0]),
            "n_up": n_up,
            "n_down": n_down,
            "enrichment_mode": enrichment_mode,
            "mode": mode,
            "real_dataset_used": real_gse,
            "group_a_name": group_a_name,
            "group_b_name": group_b_name,
            "group_a_n": group_a_n,
            "group_b_n": group_b_n,
            "notes": notes,
        },
        "group_choice": None,
        "volcano": {
            "points": volcano_points,
            "max_abs_fc": max_abs_fc,
            "max_neglog": max_neglog,
        },
        "top_up": top_up,
        "top_down": top_down,
        "enrichment_up": enrichment_up,
        "enrichment_down": enrichment_down,
        "dge_rows": df[["gene", "gene_url", "log2fc", "pvalue", "padj", "neg_log10_padj", "direction"]].to_dict(orient="records"),
    }


def save_analysis_outputs(result_path: Path, dge_path: Path, payload: dict) -> None:
    result_path.parent.mkdir(parents=True, exist_ok=True)
    dge_path.parent.mkdir(parents=True, exist_ok=True)

    result_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    df = pd.DataFrame(payload.get("dge_rows", []))
    df.to_csv(dge_path, index=False)


def load_analysis_outputs(result_path: Path) -> dict | None:
    if not result_path.exists():
        return None
    payload = json.loads(result_path.read_text(encoding="utf-8"))
    return hydrate_analysis_links(payload)
