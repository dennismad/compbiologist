from __future__ import annotations

import argparse
from pathlib import Path

from flask import Flask, redirect, render_template, request, url_for

from app.config import GEO_DEFAULT_FETCH_SIZE, GEO_DEFAULT_QUERY
from app.geo import EXPERIMENT_TYPE_OPTIONS
from app.pipeline import (
    load_cached_geo_payload,
    load_cached_loaded_geo_payload,
    load_cached_summary,
    load_processed_dataframe,
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

        summary = payload["summary"]
        metadata = payload["metadata"]
        proteins = load_processed_dataframe()
        max_bin = max(summary["length_distribution"].values() or [1])

        return render_template(
            "index.html",
            summary=summary,
            metadata=metadata,
            proteins=proteins,
            max_bin=max_bin,
            geo_search_summary={
                "query": geo_search_payload.get("query", ""),
                "source": geo_search_payload.get("source", "not_fetched"),
                "total_found": geo_search_payload.get("total_found", 0),
                "returned": len(geo_search_payload.get("items", [])),
                "species_filter": geo_search_payload.get("species_filter", ""),
                "experiment_filter": geo_search_payload.get("experiment_filter", "All"),
            },
            geo_search_items=geo_search_payload.get("items", [])[:50],
            geo_loaded_summary={
                "query": geo_loaded_payload.get("query", ""),
                "source": geo_loaded_payload.get("source", "not_loaded"),
                "returned": geo_loaded_payload.get("returned", 0),
                "species_filter": geo_loaded_payload.get("species_filter", ""),
                "experiment_filter": geo_loaded_payload.get("experiment_filter", "All"),
                "organism_distribution": geo_loaded_payload.get("insights", {}).get("organism_distribution", {}),
                "sample_distribution": geo_loaded_payload.get("insights", {}).get("sample_distribution", {}),
                "experiment_distribution": geo_loaded_payload.get("insights", {}).get("experiment_distribution", {}),
            },
            geo_loaded_items=geo_loaded_payload.get("items", [])[:50],
            geo_organism_max=max(geo_loaded_payload.get("insights", {}).get("organism_distribution", {}).values() or [1]),
            geo_sample_max=max(geo_loaded_payload.get("insights", {}).get("sample_distribution", {}).values() or [1]),
            geo_experiment_max=max(geo_loaded_payload.get("insights", {}).get("experiment_distribution", {}).values() or [1]),
            geo_default_query=GEO_DEFAULT_QUERY,
            geo_default_fetch_size=GEO_DEFAULT_FETCH_SIZE,
            experiment_type_options=EXPERIMENT_TYPE_OPTIONS,
        )

    @app.post("/refresh")
    def refresh():
        run_pipeline()
        return redirect(url_for("index"))

    @app.post("/geo/search")
    def geo_search():
        query = request.form.get("query", "").strip() or GEO_DEFAULT_QUERY
        species_filter = request.form.get("species", "").strip()
        experiment_filter = request.form.get("experiment_type", "All").strip() or "All"

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
        )
        return redirect(url_for("index"))

    @app.post("/geo/load-selected")
    def geo_load_selected():
        selected_ids = request.form.getlist("selected_ids")
        run_geo_load_selection(selected_ids)
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
