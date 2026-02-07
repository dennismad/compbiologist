from __future__ import annotations

import sqlite3
from pathlib import Path

import app.geo as geo
from app.geo import build_geo_insights, enrich_geo_items, filter_geo_items


def test_enrich_geo_items_adds_links():
    items = [
        {"accession": "GSE123", "pubmed_id": "12345"},
        {"accession": "", "pubmed_id": ""},
    ]

    enriched = enrich_geo_items(items)

    assert enriched[0]["geo_link"].endswith("acc=GSE123")
    assert enriched[0]["pubmed_link"].endswith("/12345/")
    assert enriched[1]["geo_link"] == ""
    assert enriched[1]["pubmed_link"] == ""


def test_filter_geo_items_species_and_experiment():
    items = enrich_geo_items(
        [
            {
                "accession": "GSE1",
                "title": "single-cell RNA seq breast",
                "summary": "",
                "organism": "Homo sapiens",
                "n_samples": 10,
            },
            {
                "accession": "GSE2",
                "title": "RNA-seq mouse atlas",
                "summary": "",
                "organism": "Mus musculus",
                "n_samples": 20,
            },
        ]
    )

    filtered_species = filter_geo_items(items, species_filter="homo", experiment_filter="All")
    filtered_experiment = filter_geo_items(items, species_filter="", experiment_filter="RNA-seq")

    assert len(filtered_species) == 1
    assert filtered_species[0]["accession"] == "GSE1"
    assert len(filtered_experiment) == 1
    assert filtered_experiment[0]["accession"] == "GSE2"


def test_filter_geo_items_state():
    items = enrich_geo_items(
        [
            {
                "accession": "GSE10",
                "title": "Diabetes patient and healthy control cohort",
                "summary": "",
                "organism": "Homo sapiens",
            },
            {
                "accession": "GSE11",
                "title": "Cancer patient biopsy only",
                "summary": "",
                "organism": "Homo sapiens",
            },
        ]
    )
    filtered = filter_geo_items(items, state_filter="Disease vs Healthy")
    assert len(filtered) == 1
    assert filtered[0]["accession"] == "GSE10"


def test_infer_state_profile_prioritizes_treated_over_disease_healthy():
    state = geo._infer_state_profile(
        "Diabetes patients treated with metformin and healthy control baseline",
        "",
    )
    assert state == "Treated vs Untreated"


def test_build_geo_insights_distributions():
    items = [
        {"organism": "Homo sapiens", "n_samples": 8, "experiment_type": "Single-cell RNA-seq"},
        {"organism": "Homo sapiens", "n_samples": 40, "experiment_type": "RNA-seq"},
        {"organism": "Mus musculus", "n_samples": 520, "experiment_type": "RNA-seq"},
        {"organism": "", "n_samples": "bad", "experiment_type": "Other"},
    ]

    insights = build_geo_insights(items)

    assert insights["organism_distribution"]["Homo sapiens"] == 2
    assert insights["organism_distribution"]["Mus musculus"] == 1
    assert insights["organism_distribution"]["Unknown"] == 1

    assert insights["sample_distribution"]["1-10"] == 1
    assert insights["sample_distribution"]["11-50"] == 1
    assert insights["sample_distribution"]["501-1000"] == 1
    assert insights["experiment_distribution"]["RNA-seq"] == 2


def test_search_geo_datasets_returns_error_source_when_remote_fails(monkeypatch):
    monkeypatch.setattr(geo, "GEO_SQLITE_PATH", Path("/tmp/nonexistent-geometadb.sqlite"))
    monkeypatch.setattr(geo, "get_cached_search", lambda **kwargs: None)
    monkeypatch.setattr(
        geo.requests,
        "get",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("network down")),
    )

    result = geo.search_geo_datasets(query="diabetes", retmax=5)

    assert result.source == "geo_error"
    assert result.total_found == 0
    assert result.items == []


def test_search_geo_datasets_uses_local_sqlite(monkeypatch, tmp_path):
    db_path = tmp_path / "GEOmetadb.sqlite"
    with sqlite3.connect(str(db_path)) as conn:
        conn.executescript(
            """
            CREATE TABLE gse (
              ID REAL,
              title TEXT,
              gse TEXT,
              status TEXT,
              submission_date TEXT,
              last_update_date TEXT,
              pubmed_id INTEGER,
              summary TEXT,
              type TEXT
            );
            CREATE TABLE gse_gsm (gse TEXT, gsm TEXT);
            CREATE TABLE gsm (ID REAL, title TEXT, gsm TEXT, organism_ch1 TEXT);
            INSERT INTO gse (ID, title, gse, last_update_date, pubmed_id, summary, type)
              VALUES (1, 'Diabetes treated cohort', 'GSEX1', '2024-01-01', 12345, 'treated and untreated study', 'Expression profiling by high throughput sequencing');
            INSERT INTO gse (ID, title, gse, last_update_date, pubmed_id, summary, type)
              VALUES (2, 'Unrelated cancer cohort', 'GSEX2', '2024-01-02', 0, 'tumor study', 'Expression profiling by array');
            INSERT INTO gse_gsm (gse, gsm) VALUES ('GSEX1', 'GSM1'), ('GSEX1', 'GSM2'), ('GSEX1', 'GSM3'), ('GSEX1', 'GSM4');
            INSERT INTO gsm (ID, title, gsm, organism_ch1) VALUES
              (1, 's1', 'GSM1', 'Homo sapiens'),
              (2, 's2', 'GSM2', 'Homo sapiens'),
              (3, 's3', 'GSM3', 'Homo sapiens'),
              (4, 's4', 'GSM4', 'Homo sapiens');
            """
        )
        conn.commit()

    monkeypatch.setattr(geo, "GEO_SQLITE_PATH", db_path)
    monkeypatch.setattr(geo, "get_cached_search", lambda **kwargs: None)

    result = geo.search_geo_datasets(
        query="diabetes",
        retmax=10,
        species_filter="Homo sapiens",
        experiment_filter="RNA-seq",
        state_filter="Treated vs Untreated",
    )

    assert result.source == "geo_sqlite"
    assert result.total_found == 1
    assert len(result.items) == 1
    assert result.items[0]["accession"] == "GSEX1"


def test_sqlite_state_filter_returns_only_matching_inferred_state(monkeypatch, tmp_path):
    db_path = tmp_path / "GEOmetadb.sqlite"
    with sqlite3.connect(str(db_path)) as conn:
        conn.executescript(
            """
            CREATE TABLE gse (
              ID REAL,
              title TEXT,
              gse TEXT,
              status TEXT,
              submission_date TEXT,
              last_update_date TEXT,
              pubmed_id INTEGER,
              summary TEXT,
              type TEXT
            );
            CREATE TABLE gse_gsm (gse TEXT, gsm TEXT);
            CREATE TABLE gsm (ID REAL, title TEXT, gsm TEXT, organism_ch1 TEXT);

            INSERT INTO gse (ID, title, gse, last_update_date, pubmed_id, summary, type)
              VALUES (1, 'Diabetes patient and healthy control cohort', 'GSEX_DH', '2024-01-03', 11111, 'disease vs healthy', 'Expression profiling by high throughput sequencing');
            INSERT INTO gse (ID, title, gse, last_update_date, pubmed_id, summary, type)
              VALUES (2, 'Diabetes treated cohort with healthy controls', 'GSEX_TU', '2024-01-02', 22222, 'treated baseline intervention', 'Expression profiling by high throughput sequencing');
            """
        )
        conn.commit()

    monkeypatch.setattr(geo, "GEO_SQLITE_PATH", db_path)
    monkeypatch.setattr(geo, "get_cached_search", lambda **kwargs: None)

    result = geo.search_geo_datasets(
        query="diabetes",
        retmax=20,
        state_filter="Disease vs Healthy",
        experiment_filter="All",
    )

    assert result.source == "geo_sqlite"
    assert len(result.items) >= 1
    assert all(item.get("state_profile") == "Disease vs Healthy" for item in result.items)
