from __future__ import annotations

import argparse
from pathlib import Path

from flask import Flask, redirect, render_template, request, url_for

from app.analysis import DEFAULT_ENRICHMENT_MODE, ENRICHMENT_MODE_OPTIONS
from app.config import GEO_DEFAULT_FETCH_SIZE, GEO_DEFAULT_QUERY
from app.geo import EXPERIMENT_TYPE_OPTIONS, STATE_FILTER_OPTIONS, get_common_species_options
from app.pipeline import (
    load_cached_analysis_payload,
    load_cached_geo_payload,
    load_cached_loaded_geo_payload,
    load_cached_summary,
    reset_workflow_state,
    run_comparison_analysis,
    run_geo_load_selection,
    run_geo_pipeline,
    run_pipeline,
)


def create_app() -> Flask:
    base_dir = Path(__file__).resolve().parent.parent
    app = Flask(
        __name__,
        template_folder=str(base_dir / "templates"),
        static_folder=str(base_dir / "static"),
    )

    @app.get("/")
    def index():
        payload = load_cached_summary()
        if payload is None:
            payload = run_pipeline()

        geo_search_payload = load_cached_geo_payload()
        geo_loaded_payload = load_cached_loaded_geo_payload()
        analysis_payload = load_cached_analysis_payload()
        species_options = get_common_species_options(limit=5)
        current_species_filter = str(geo_search_payload.get("species_filter", "")).strip()
        selected_species_option = current_species_filter if current_species_filter in species_options else ("__OTHER__" if current_species_filter else "")
        other_species_value = current_species_filter if selected_species_option == "__OTHER__" else ""

        return render_template(
            "index.html",
            geo_search_summary={
                "query": geo_search_payload.get("query", ""),
                "source": geo_search_payload.get("source", "not_fetched"),
                "error": geo_search_payload.get("error", ""),
                "total_found": geo_search_payload.get("total_found", 0),
                "returned": len(geo_search_payload.get("items", [])),
                "requested_retmax": geo_search_payload.get("requested_retmax", GEO_DEFAULT_FETCH_SIZE),
                "returned_before_analyzable_filter": geo_search_payload.get("returned_before_analyzable_filter", len(geo_search_payload.get("items", []))),
                "species_filter": geo_search_payload.get("species_filter", ""),
                "experiment_filter": geo_search_payload.get("experiment_filter", "All"),
                "state_filter": geo_search_payload.get("state_filter", "All"),
                "only_analyzable": geo_search_payload.get("only_analyzable", True),
            },
            geo_search_items=geo_search_payload.get("items", [])[:80],
            geo_loaded_summary={
                "query": geo_loaded_payload.get("query", ""),
                "source": geo_loaded_payload.get("source", "not_loaded"),
                "returned": geo_loaded_payload.get("returned", 0),
                "species_filter": geo_loaded_payload.get("species_filter", ""),
                "experiment_filter": geo_loaded_payload.get("experiment_filter", "All"),
                "state_filter": geo_loaded_payload.get("state_filter", "All"),
            },
            geo_loaded_items=geo_loaded_payload.get("items", [])[:80],
            geo_default_query=GEO_DEFAULT_QUERY,
            geo_default_fetch_size=GEO_DEFAULT_FETCH_SIZE,
            species_options=species_options,
            selected_species_option=selected_species_option,
            other_species_value=other_species_value,
            experiment_type_options=EXPERIMENT_TYPE_OPTIONS,
            state_filter_options=STATE_FILTER_OPTIONS,
            enrichment_mode_options=ENRICHMENT_MODE_OPTIONS,
            analysis=analysis_payload,
            protein_metadata=payload.get("metadata", {}),
        )

    @app.post("/refresh")
    def refresh():
        run_pipeline()
        return redirect(url_for("index"))

    @app.post("/workflow/reset")
    def workflow_reset():
        reset_workflow_state()
        return redirect(url_for("index"))

    @app.post("/geo/search")
    def geo_search():
        query = request.form.get("query", "").strip() or GEO_DEFAULT_QUERY
        species_choice = request.form.get("species_choice", "").strip()
        species_other = request.form.get("species_other", "").strip()
        legacy_species = request.form.get("species", "").strip()
        if species_choice == "__OTHER__":
            species_filter = species_other
        elif species_choice:
            species_filter = species_choice
        else:
            species_filter = legacy_species
        experiment_filter = request.form.get("experiment_type", "All").strip() or "All"
        state_filter = request.form.get("state_filter", "All").strip() or "All"
        only_analyzable = request.form.get("only_analyzable", "") == "on"

        retmax_raw = request.form.get("retmax", str(GEO_DEFAULT_FETCH_SIZE)).strip()
        try:
            retmax = int(retmax_raw)
        except ValueError:
            retmax = GEO_DEFAULT_FETCH_SIZE
        retmax = max(1, min(retmax, 100))

        run_geo_pipeline(
            query=query,
            retmax=retmax,
            species_filter=species_filter,
            experiment_filter=experiment_filter,
            state_filter=state_filter,
            only_analyzable=only_analyzable,
        )
        return redirect(url_for("index"))

    @app.post("/geo/load-selected")
    def geo_load_selected():
        selected_ids = request.form.getlist("selected_ids")
        run_geo_load_selection(selected_ids)
        return redirect(url_for("index"))

    @app.post("/analysis/run")
    def run_analysis():
        state_profile = request.form.get("analysis_state", "Disease vs Healthy").strip() or "Disease vs Healthy"
        enrichment_mode = request.form.get("enrichment_mode", DEFAULT_ENRICHMENT_MODE).strip().lower()
        padj_raw = request.form.get("padj_cutoff", "0.05").strip()
        lfc_raw = request.form.get("log2fc_cutoff", "1.0").strip()

        try:
            padj_cutoff = float(padj_raw)
        except ValueError:
            padj_cutoff = 0.05
        try:
            log2fc_cutoff = float(lfc_raw)
        except ValueError:
            log2fc_cutoff = 1.0

        padj_cutoff = min(max(padj_cutoff, 1e-6), 1.0)
        log2fc_cutoff = min(max(log2fc_cutoff, 0.1), 5.0)

        run_comparison_analysis(
            state_profile=state_profile,
            padj_cutoff=padj_cutoff,
            log2fc_cutoff=log2fc_cutoff,
            enrichment_mode=enrichment_mode,
            manual_choice=None,
        )
        return redirect(url_for("index"))

    @app.post("/analysis/run-manual")
    def run_analysis_manual():
        state_profile = request.form.get("analysis_state", "Disease vs Healthy").strip() or "Disease vs Healthy"
        enrichment_mode = request.form.get("enrichment_mode", DEFAULT_ENRICHMENT_MODE).strip().lower()
        padj_raw = request.form.get("padj_cutoff", "0.05").strip()
        lfc_raw = request.form.get("log2fc_cutoff", "1.0").strip()

        try:
            padj_cutoff = float(padj_raw)
        except ValueError:
            padj_cutoff = 0.05
        try:
            log2fc_cutoff = float(lfc_raw)
        except ValueError:
            log2fc_cutoff = 1.0

        padj_cutoff = min(max(padj_cutoff, 1e-6), 1.0)
        log2fc_cutoff = min(max(log2fc_cutoff, 0.1), 5.0)

        sample_ids = request.form.getlist("sample_ids")
        sample_groups = request.form.getlist("sample_groups")
        group_a_samples: list[str] = []
        group_b_samples: list[str] = []
        if sample_ids and len(sample_ids) == len(sample_groups):
            for sid, group_code in zip(sample_ids, sample_groups):
                sid_val = sid.strip()
                if not sid_val:
                    continue
                if group_code == "A":
                    group_a_samples.append(sid_val)
                elif group_code == "B":
                    group_b_samples.append(sid_val)
        else:
            # Backward-compatible fallback if older form fields are posted.
            group_a_samples = request.form.getlist("group_a_samples")
            group_b_samples = request.form.getlist("group_b_samples")

        manual_choice = {
            "gse": request.form.get("manual_gse", "").strip(),
            "group_a_name": request.form.get("group_a_name", "Group A").strip() or "Group A",
            "group_b_name": request.form.get("group_b_name", "Group B").strip() or "Group B",
            "group_a_samples": group_a_samples,
            "group_b_samples": group_b_samples,
        }

        run_comparison_analysis(
            state_profile=state_profile,
            padj_cutoff=padj_cutoff,
            log2fc_cutoff=log2fc_cutoff,
            enrichment_mode=enrichment_mode,
            manual_choice=manual_choice,
        )
        return redirect(url_for("index"))

    return app


def main() -> None:
    parser = argparse.ArgumentParser(description="Virtual computational biologist web app")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--refresh", action="store_true", help="Run protein data pipeline once and exit")
    parser.add_argument("--geo-query", default="", help="Optional GEO query term for CLI fetch")
    parser.add_argument("--geo-retmax", type=int, default=GEO_DEFAULT_FETCH_SIZE)
    args = parser.parse_args()

    if args.refresh:
        run_pipeline()
        if args.geo_query:
            run_geo_pipeline(query=args.geo_query, retmax=max(1, min(args.geo_retmax, 100)))
        return

    app = create_app()
    app.run(host=args.host, port=args.port, debug=True)


if __name__ == "__main__":
    main()
