from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import app.pipeline as pipeline


SAMPLE_TSV = "\n".join(
    [
        "Entry\tEntry Name\tProtein names\tGene Names\tLength\tAnnotation",
        "P1\tA_HUMAN\tAlpha\tGENEA\t100\t5",
        "P2\tB_HUMAN\tBeta\tGENEB\t450\t4",
    ]
)


def test_run_pipeline_writes_summary_and_processed(monkeypatch, tmp_path):
    raw_path = tmp_path / "raw.tsv"
    processed_path = tmp_path / "processed.csv"
    summary_path = tmp_path / "summary.json"

    def fake_fetch(output_path):
        output_path.write_text(SAMPLE_TSV, encoding="utf-8")
        return SimpleNamespace(source="test", bytes_written=len(SAMPLE_TSV.encode("utf-8")))

    monkeypatch.setattr(pipeline, "RAW_DATA_PATH", raw_path)
    monkeypatch.setattr(pipeline, "PROCESSED_DATA_PATH", processed_path)
    monkeypatch.setattr(pipeline, "SUMMARY_PATH", summary_path)
    monkeypatch.setattr(pipeline, "fetch_uniprot_data", fake_fetch)

    payload = pipeline.run_pipeline()

    assert payload["metadata"]["data_source"] == "test"
    assert payload["summary"]["n_proteins"] == 2
    assert processed_path.exists()
    assert summary_path.exists()


def test_load_cached_geo_payload_enriches_items(monkeypatch):
    payload = {
        "query": "breast cancer",
        "source": "geo",
        "total_found": 1,
        "items": [{"accession": "GSE1", "pubmed_id": "100"}],
    }

    monkeypatch.setattr(pipeline, "load_cached_geo", lambda _: payload)

    out = pipeline.load_cached_geo_payload()

    assert out["items"][0]["geo_link"].endswith("acc=GSE1")
    assert out["items"][0]["pubmed_link"].endswith("/100/")
    assert "organism_distribution" in out["insights"]
    assert "sample_distribution" in out["insights"]


def test_run_geo_load_selection_saves_subset(monkeypatch, tmp_path):
    selected_payload = {
        "query": "brain",
        "source": "geo",
        "species_filter": "Homo sapiens",
        "experiment_filter": "RNA-seq",
        "items": [
            {"accession": "GSE1", "title": "A", "analyzable": False},
            {"accession": "GSE2", "title": "B", "analyzable": True},
        ],
    }

    monkeypatch.setattr(pipeline, "GEO_LOADED_JSON_PATH", tmp_path / "loaded.json")
    monkeypatch.setattr(pipeline, "GEO_LOADED_CSV_PATH", tmp_path / "loaded.csv")
    monkeypatch.setattr(pipeline, "load_cached_geo_payload", lambda: selected_payload)
    monkeypatch.setattr(
        pipeline,
        "annotate_items_analyzable",
        lambda items, cache_dir, check_limit, strict_validation=True: items,
    )

    out = pipeline.run_geo_load_selection(["GSE2"])

    assert out["returned"] == 1
    assert out["items"][0]["accession"] == "GSE2"
    assert (tmp_path / "loaded.json").exists()
    assert (tmp_path / "loaded.csv").exists()


def test_run_geo_pipeline_overfetches_for_analyzable(monkeypatch):
    pages = [
        {
            "source": "geo",
            "query": "diabetes",
            "total_found": 200,
            "items": [{"accession": "GSE1"}, {"accession": "GSE2"}],
        },
        {
            "source": "geo",
            "query": "diabetes",
            "total_found": 200,
            "items": [{"accession": "GSE3"}, {"accession": "GSE4"}],
        },
    ]
    calls = []

    def fake_search(query, retmax, retstart, species_filter, experiment_filter, state_filter):
        calls.append(retstart)
        idx = 0 if retstart == 0 else 1
        payload = pages[idx]
        return pipeline.GEOSearchResult(
            source=payload["source"],
            query=payload["query"],
            total_found=payload["total_found"],
            items=payload["items"],
        )

    def fake_annotate(items, cache_dir, check_limit, strict_validation=True):
        out = []
        for row in items:
            copy = dict(row)
            copy["analyzable"] = copy["accession"] in {"GSE2", "GSE3", "GSE4"}
            copy["analyzable_detail"] = "series_matrix" if copy["analyzable"] else "missing_matrix"
            out.append(copy)
        return out

    monkeypatch.setattr(pipeline, "search_geo_datasets", fake_search)
    monkeypatch.setattr(pipeline, "annotate_items_analyzable", fake_annotate)
    monkeypatch.setattr(pipeline, "GEO_RAW_JSON_PATH", Path("/tmp/geo_raw_test.json"))
    monkeypatch.setattr(pipeline, "GEO_PROCESSED_PATH", Path("/tmp/geo_processed_test.csv"))
    monkeypatch.setattr(pipeline, "GEO_SUMMARY_PATH", Path("/tmp/geo_summary_test.json"))

    out = pipeline.run_geo_pipeline(
        query="diabetes",
        retmax=3,
        species_filter="Homo sapiens",
        experiment_filter="RNA-seq",
        state_filter="Disease vs Healthy",
        only_analyzable=True,
    )

    assert len(out["items"]) == 3
    assert [x["accession"] for x in out["items"]] == ["GSE2", "GSE3", "GSE4"]
    assert calls == [0, 50]
