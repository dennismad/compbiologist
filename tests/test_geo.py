from __future__ import annotations

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
