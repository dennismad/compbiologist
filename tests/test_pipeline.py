from __future__ import annotations

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
            {"accession": "GSE1", "title": "A"},
            {"accession": "GSE2", "title": "B"},
        ],
    }

    monkeypatch.setattr(pipeline, "GEO_LOADED_JSON_PATH", tmp_path / "loaded.json")
    monkeypatch.setattr(pipeline, "GEO_LOADED_CSV_PATH", tmp_path / "loaded.csv")
    monkeypatch.setattr(pipeline, "load_cached_geo_payload", lambda: selected_payload)

    out = pipeline.run_geo_load_selection(["GSE2"])

    assert out["returned"] == 1
    assert out["items"][0]["accession"] == "GSE2"
    assert (tmp_path / "loaded.json").exists()
    assert (tmp_path / "loaded.csv").exists()
