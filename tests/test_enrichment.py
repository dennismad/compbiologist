from __future__ import annotations

import pandas as pd

from app.enrichment import infer_gprofiler_organism, run_gprofiler_enrichment


def test_infer_gprofiler_organism():
    assert infer_gprofiler_organism([{"organism": "Homo sapiens"}]) == "hsapiens"
    assert infer_gprofiler_organism([{"organism": "Mus musculus"}]) == "mmusculus"
    assert infer_gprofiler_organism([{"organism": "Rattus norvegicus"}]) == "rnorvegicus"
    assert infer_gprofiler_organism([{"organism": "Unknown"}]) == "hsapiens"


def test_run_gprofiler_enrichment_normalizes_ensembl_and_uniprot(monkeypatch):
    captured: dict = {}

    class FakeGP:
        def __init__(self, return_dataframe: bool = True):
            assert return_dataframe is True

        def profile(self, **kwargs):
            captured.update(kwargs)
            return pd.DataFrame(
                [
                    {
                        "source": "GO:BP",
                        "native": "GO:0006954",
                        "name": "inflammatory response",
                        "p_value": 1e-4,
                        "intersection_size": 3,
                        "term_size": 120,
                    }
                ]
            )

    monkeypatch.setattr("app.enrichment._load_gprofiler", lambda: (FakeGP, None))

    rows, note = run_gprofiler_enrichment(
        gene_ids=["ENSG00000141510.18", "P04637-2", "P04637-2"],
        organism="hsapiens",
    )

    assert note is None
    assert rows[0]["pathway"].startswith("GO:BP:GO:0006954")
    assert rows[0]["overlap"] == 3
    assert rows[0]["set_size"] == 120
    assert rows[0]["padj"] == 1e-4
    assert captured["query"] == ["ENSG00000141510", "P04637"]


def test_run_gprofiler_enrichment_reports_missing_dependency(monkeypatch):
    monkeypatch.setattr("app.enrichment._load_gprofiler", lambda: (None, "gprofiler-official unavailable"))
    rows, note = run_gprofiler_enrichment(
        gene_ids=["ENSG00000141510", "ENSG00000157764"],
        organism="hsapiens",
    )
    assert rows == []
    assert "gprofiler-official unavailable" in str(note)
