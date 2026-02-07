from __future__ import annotations

import json
from datetime import datetime, timezone

from app.analysis import load_analysis_outputs, run_state_comparison_analysis, save_analysis_outputs
from app.config import (
    ANALYSIS_DGE_PATH,
    ANALYSIS_RESULT_PATH,
    GEO_DEFAULT_FETCH_SIZE,
    GEO_MATRIX_CACHE_DIR,
    GEO_LOADED_CSV_PATH,
    GEO_LOADED_JSON_PATH,
    GEO_PROCESSED_PATH,
    GEO_RAW_JSON_PATH,
    GEO_SUMMARY_PATH,
    PROCESSED_DATA_PATH,
    RAW_DATA_PATH,
    SUMMARY_PATH,
)
from app.geo_expression import annotate_items_analyzable
from app.data_sources import fetch_uniprot_data
from app.geo import (
    GEOSearchResult,
    build_geo_insights,
    enrich_geo_items,
    load_cached_geo,
    load_cached_loaded_geo,
    save_loaded_geo_selection,
    search_geo_datasets,
    write_geo_artifacts,
)
from app.processing import process_protein_table


def _dedupe_geo_items(items: list[dict]) -> list[dict]:
    out: list[dict] = []
    seen: set[str] = set()
    for row in items:
        key = str(row.get("accession") or row.get("uid") or "").strip()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def _dataset_id(row: dict) -> str:
    return str(row.get("accession") or row.get("uid") or "").strip()


def ensure_data_dirs() -> None:
    RAW_DATA_PATH.parent.mkdir(parents=True, exist_ok=True)
    PROCESSED_DATA_PATH.parent.mkdir(parents=True, exist_ok=True)


def run_pipeline() -> dict:
    ensure_data_dirs()

    fetch_result = fetch_uniprot_data(RAW_DATA_PATH)
    processed = process_protein_table(str(RAW_DATA_PATH))

    processed.dataframe.to_csv(PROCESSED_DATA_PATH, index=False)

    metadata = {
        "last_updated_utc": datetime.now(timezone.utc).isoformat(),
        "data_source": fetch_result.source,
        "raw_file": str(RAW_DATA_PATH),
        "processed_file": str(PROCESSED_DATA_PATH),
        "bytes_written": fetch_result.bytes_written,
    }

    payload = {
        "summary": processed.summary,
        "metadata": metadata,
    }

    SUMMARY_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def run_geo_pipeline(
    query: str,
    retmax: int,
    species_filter: str = "",
    experiment_filter: str = "All",
    state_filter: str = "All",
    only_analyzable: bool = True,
) -> dict:
    target = max(1, int(retmax))
    page_size = min(200, max(target * 3, 50))
    retstart = 0
    total_found = 0
    source = "geo"
    accepted: list[dict] = []
    accepted_before_analyzable = 0
    max_pages = 12

    for _ in range(max_pages):
        page = search_geo_datasets(
            query=query,
            retmax=page_size,
            retstart=retstart,
            species_filter=species_filter,
            experiment_filter=experiment_filter,
            state_filter=state_filter,
        )
        source = page.source
        total_found = page.total_found
        if not page.items:
            break

        filtered_items = enrich_geo_items(page.items)
        check_limit = len(filtered_items) if only_analyzable else 40
        annotated_page = annotate_items_analyzable(
            filtered_items,
            cache_dir=GEO_MATRIX_CACHE_DIR,
            check_limit=check_limit,
            strict_validation=False,
        )
        annotated_page = _dedupe_geo_items(annotated_page)
        accepted_before_analyzable += len(annotated_page)

        if only_analyzable:
            annotated_page = [x for x in annotated_page if x.get("analyzable")]
        accepted.extend(annotated_page)
        accepted = _dedupe_geo_items(accepted)
        if len(accepted) >= target:
            break

        retstart += page_size
        if retstart >= total_found:
            break
        if page.source != "geo":
            break

    accepted = accepted[:target]
    filtered_result = GEOSearchResult(
        source=source,
        query=query,
        total_found=total_found,
        items=accepted,
    )

    summary = write_geo_artifacts(
        raw_path=GEO_RAW_JSON_PATH,
        processed_path=GEO_PROCESSED_PATH,
        summary_path=GEO_SUMMARY_PATH,
        result=filtered_result,
        species_filter=species_filter,
        experiment_filter=experiment_filter,
        state_filter=state_filter,
        only_analyzable=only_analyzable,
        requested_retmax=target,
        returned_before_analyzable_filter=accepted_before_analyzable,
    )
    return {
        "summary": summary,
        "items": accepted,
    }


def run_geo_load_selection(selected_ids: list[str]) -> dict:
    source_payload = load_cached_geo_payload()
    selected_set = {str(x).strip() for x in selected_ids if str(x).strip()}
    source_items = list(source_payload.get("items", []))
    selected_rows = [row for row in source_items if _dataset_id(row) in selected_set]

    if selected_rows:
        validated = annotate_items_analyzable(
            selected_rows,
            cache_dir=GEO_MATRIX_CACHE_DIR,
            check_limit=len(selected_rows),
            strict_validation=True,
        )
        validated_map = {_dataset_id(row): row for row in validated}
        source_payload["items"] = [validated_map.get(_dataset_id(row), row) for row in source_items]

    return save_loaded_geo_selection(
        loaded_json_path=GEO_LOADED_JSON_PATH,
        loaded_csv_path=GEO_LOADED_CSV_PATH,
        source_payload=source_payload,
        selected_ids=selected_ids,
    )


def run_comparison_analysis(
    state_profile: str,
    padj_cutoff: float,
    log2fc_cutoff: float,
    manual_choice: dict | None = None,
) -> dict:
    loaded = load_cached_loaded_geo_payload()
    payload = run_state_comparison_analysis(
        state_profile=state_profile,
        loaded_items=loaded.get("items", []),
        padj_cutoff=padj_cutoff,
        log2fc_cutoff=log2fc_cutoff,
        use_real_geo=True,
        manual_choice=manual_choice,
    )
    save_analysis_outputs(ANALYSIS_RESULT_PATH, ANALYSIS_DGE_PATH, payload)
    return payload


def load_cached_summary() -> dict | None:
    if not SUMMARY_PATH.exists():
        return None
    return json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))


def load_processed_dataframe() -> list[dict]:
    if not PROCESSED_DATA_PATH.exists():
        return []

    import pandas as pd

    df = pd.read_csv(PROCESSED_DATA_PATH)
    return df.head(25).fillna("").to_dict(orient="records")


def load_cached_geo_payload() -> dict:
    payload = load_cached_geo(GEO_RAW_JSON_PATH)
    if payload is None:
        return {
            "query": "",
            "source": "not_fetched",
            "total_found": 0,
            "species_filter": "",
            "experiment_filter": "All",
            "state_filter": "All",
            "only_analyzable": True,
            "requested_retmax": GEO_DEFAULT_FETCH_SIZE,
            "items": [],
            "insights": {
                "organism_distribution": {},
                "sample_distribution": {},
                "experiment_distribution": {},
            },
        }
    enriched_items = enrich_geo_items(payload.get("items", []))
    for row in enriched_items:
        if "analyzable" not in row:
            row["analyzable"] = False
        if "analyzable_detail" not in row:
            row["analyzable_detail"] = "unknown"
    payload["items"] = enriched_items
    payload["insights"] = build_geo_insights(enriched_items)
    return payload


def load_cached_loaded_geo_payload() -> dict:
    payload = load_cached_loaded_geo(GEO_LOADED_JSON_PATH)
    if payload is None:
        return {
            "query": "",
            "source": "not_loaded",
            "species_filter": "",
            "experiment_filter": "All",
            "state_filter": "All",
            "returned": 0,
            "selected_ids": [],
            "items": [],
            "insights": {
                "organism_distribution": {},
                "sample_distribution": {},
                "experiment_distribution": {},
            },
        }

    enriched_items = enrich_geo_items(payload.get("items", []))
    annotated_items = annotate_items_analyzable(enriched_items, cache_dir=GEO_MATRIX_CACHE_DIR, check_limit=len(enriched_items))
    analyzable_items = [x for x in annotated_items if x.get("analyzable")]
    payload["items"] = analyzable_items
    payload["insights"] = build_geo_insights(analyzable_items)
    payload["returned"] = len(analyzable_items)
    return payload


def load_cached_analysis_payload() -> dict | None:
    return load_analysis_outputs(ANALYSIS_RESULT_PATH)
