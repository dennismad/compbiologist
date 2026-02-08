from __future__ import annotations

from app.analysis import (
    _annotate_top_gene_rows,
    _ensembl_gene_url,
    _pathway_search_url,
    hydrate_analysis_links,
    run_state_comparison_analysis,
)


def test_run_state_comparison_analysis_returns_expected_sections_prototype_mode():
    loaded = [
        {"accession": "GSE1"},
        {"accession": "GSE2"},
    ]

    result = run_state_comparison_analysis(
        state_profile="Disease vs Healthy",
        loaded_items=loaded,
        padj_cutoff=0.05,
        log2fc_cutoff=1.0,
        use_real_geo=False,
    )

    assert "metadata" in result
    assert "volcano" in result
    assert "top_up" in result
    assert "top_down" in result
    assert "enrichment_up" in result
    assert "enrichment_down" in result
    assert "dge_rows" in result

    assert result["metadata"]["n_datasets"] == 2
    assert result["metadata"]["n_genes"] > 0
    assert len(result["volcano"]["points"]) == result["metadata"]["n_genes"]
    assert result["metadata"]["mode"] == "prototype_signature_bootstrap"


def test_analysis_requires_manual_groups_when_real_fails(monkeypatch):
    loaded = [{"accession": "GSE290559"}]

    monkeypatch.setattr(
        "app.analysis.run_real_geo_de",
        lambda loaded_items, state_profile, cache_dir, manual_choice=None: (None, ["GSE290559: no auto groups"]),
    )
    monkeypatch.setattr(
        "app.analysis.build_group_choice_context",
        lambda loaded_items, state_profile, cache_dir: (
            {
                "gse": "GSE290559",
                "source": "supplementary",
                "sample_options": [{"id": "S1", "label": "s1"}, {"id": "S2", "label": "s2"}, {"id": "S3", "label": "s3"}, {"id": "S4", "label": "s4"}],
                "suggested_group_a": [],
                "suggested_group_b": [],
                "group_a_name": "Treated",
                "group_b_name": "Untreated",
            },
            [],
        ),
    )

    result = run_state_comparison_analysis(
        state_profile="Treated vs Untreated",
        loaded_items=loaded,
        padj_cutoff=0.05,
        log2fc_cutoff=1.0,
        use_real_geo=True,
    )

    assert result["metadata"]["mode"] == "needs_group_selection"
    assert result["metadata"]["n_genes"] == 0
    assert result["group_choice"]["gse"] == "GSE290559"


def test_analysis_reports_missing_pydeseq2_dependency(monkeypatch):
    loaded = [{"accession": "GSE1"}]

    monkeypatch.setattr(
        "app.analysis.run_real_geo_de",
        lambda loaded_items, state_profile, cache_dir, manual_choice=None: (
            None,
            ["GSE1: pydeseq2 is required for real DE analysis. Install dependencies and rerun."],
        ),
    )
    monkeypatch.setattr(
        "app.analysis.build_group_choice_context",
        lambda loaded_items, state_profile, cache_dir: (
            {
                "gse": "GSE1",
                "source": "supplementary",
                "sample_options": [{"id": "S1", "label": "s1"}, {"id": "S2", "label": "s2"}],
                "suggested_group_a": [],
                "suggested_group_b": [],
                "group_a_name": "A",
                "group_b_name": "B",
            },
            [],
        ),
    )

    result = run_state_comparison_analysis(
        state_profile="Treated vs Untreated",
        loaded_items=loaded,
        padj_cutoff=0.05,
        log2fc_cutoff=1.0,
        use_real_geo=True,
    )

    assert result["metadata"]["mode"] == "missing_dependency"
    assert result["group_choice"] is None


def test_manual_group_choices_are_preserved_when_validation_fails(monkeypatch):
    loaded = [{"accession": "GSE232764"}]

    monkeypatch.setattr(
        "app.analysis.run_real_geo_de",
        lambda loaded_items, state_profile, cache_dir, manual_choice=None: (
            None,
            ["GSE232764: Manual groups need >=2 samples per group; got 1 vs 1"],
        ),
    )
    monkeypatch.setattr(
        "app.analysis.build_group_choice_context",
        lambda loaded_items, state_profile, cache_dir: (
            {
                "gse": "GSE232764",
                "source": "supplementary",
                "sample_options": [
                    {"id": "S1", "label": "sample 1 long annotation"},
                    {"id": "S2", "label": "sample 2 long annotation"},
                    {"id": "S3", "label": "sample 3 long annotation"},
                    {"id": "S4", "label": "sample 4 long annotation"},
                ],
                "suggested_group_a": ["S1", "S2"],
                "suggested_group_b": ["S3", "S4"],
                "group_a_name": "Case",
                "group_b_name": "Control",
            },
            [],
        ),
    )

    manual_choice = {
        "gse": "GSE232764",
        "group_a_name": "Disease",
        "group_b_name": "Healthy",
        "group_a_samples": ["S1"],
        "group_b_samples": ["S3"],
    }

    result = run_state_comparison_analysis(
        state_profile="Disease vs Healthy",
        loaded_items=loaded,
        padj_cutoff=0.05,
        log2fc_cutoff=1.0,
        use_real_geo=True,
        manual_choice=manual_choice,
    )

    assert result["metadata"]["mode"] == "needs_group_selection"
    assert result["group_choice"]["group_a_name"] == "Disease"
    assert result["group_choice"]["group_b_name"] == "Healthy"
    assert result["group_choice"]["suggested_group_a"] == ["S1"]
    assert result["group_choice"]["suggested_group_b"] == ["S3"]


def test_analysis_uses_gprofiler_results_when_available(monkeypatch):
    monkeypatch.setattr(
        "app.analysis.run_gprofiler_enrichment",
        lambda gene_ids, organism, max_terms=10, user_threshold=0.05: (
            [
                {
                    "pathway": "GO:BP:GO:0006954 inflammatory response",
                    "overlap": 3,
                    "set_size": 100,
                    "padj": 0.001,
                }
            ],
            None,
        ),
    )

    result = run_state_comparison_analysis(
        state_profile="Disease vs Healthy",
        loaded_items=[{"accession": "GSE1", "organism": "Homo sapiens"}],
        padj_cutoff=0.05,
        log2fc_cutoff=1.0,
        use_real_geo=False,
    )

    assert result["enrichment_up"]
    assert result["enrichment_up"][0]["pathway"].startswith("GO:BP:")


def test_analysis_local_enrichment_mode_skips_gprofiler(monkeypatch):
    def _fail_if_called(*args, **kwargs):
        raise AssertionError("gprofiler should not be called in local mode")

    monkeypatch.setattr("app.analysis.run_gprofiler_enrichment", _fail_if_called)

    result = run_state_comparison_analysis(
        state_profile="Disease vs Healthy",
        loaded_items=[{"accession": "GSE1"}],
        padj_cutoff=0.05,
        log2fc_cutoff=1.0,
        enrichment_mode="local",
        use_real_geo=False,
    )

    assert result["metadata"]["enrichment_mode"] == "local"
    assert isinstance(result["enrichment_up"], list)


def test_analysis_gprofiler_only_mode_does_not_fallback(monkeypatch):
    monkeypatch.setattr(
        "app.analysis.run_gprofiler_enrichment",
        lambda gene_ids, organism, max_terms=10, user_threshold=0.05: ([], "gprofiler unavailable for test"),
    )

    result = run_state_comparison_analysis(
        state_profile="Disease vs Healthy",
        loaded_items=[{"accession": "GSE1"}],
        padj_cutoff=0.05,
        log2fc_cutoff=1.0,
        enrichment_mode="gprofiler",
        use_real_geo=False,
    )

    assert result["metadata"]["enrichment_mode"] == "gprofiler"
    assert result["enrichment_up"] == []
    assert any("gprofiler unavailable for test" in note for note in result["metadata"]["notes"])


def test_ensembl_gene_url_strips_ensembl_version_suffix():
    url = _ensembl_gene_url("ENSG00000141510.18")
    assert url == "https://www.ensembl.org/id/ENSG00000141510"


def test_ensembl_gene_url_searches_non_ensembl_ids():
    url = _ensembl_gene_url("TP53")
    assert "https://www.ensembl.org/Multi/Search/Results" in url
    assert "TP53" in url


def test_pathway_search_url_supports_go_reactome_wp_and_kegg():
    assert _pathway_search_url("GO:BP:GO:0006954 inflammatory response") == "https://www.ebi.ac.uk/QuickGO/term/GO:0006954"
    assert _pathway_search_url("REAC:R-HSA-168256 Immune System") == "https://reactome.org/content/detail/R-HSA-168256"
    assert _pathway_search_url("REAC:REAC:R-HSA-9031628 NGF-stimulated transcription") == "https://reactome.org/content/detail/R-HSA-9031628"
    assert _pathway_search_url("WP:WP554 DNA damage response") == "https://www.wikipathways.org/pathways/WP554.html"
    assert _pathway_search_url("KEGG:hsa04110 Cell cycle") == "https://www.kegg.jp/entry/hsa04110"


def test_pathway_search_url_uses_curated_local_pathway_links():
    assert _pathway_search_url("Inflammatory_Response") == "https://www.ebi.ac.uk/QuickGO/term/GO:0006954"


def test_hydrate_analysis_links_backfills_missing_gene_and_pathway_urls():
    payload = {
        "top_up": [{"gene": "ENSG00000141510", "log2fc": 1.2, "padj": 0.01}],
        "top_down": [{"gene": "TP53", "log2fc": -1.0, "padj": 0.02}],
        "enrichment_up": [{"pathway": "GO:BP:GO:0006954 inflammatory response", "overlap": 3, "set_size": 50, "padj": 1e-4}],
        "enrichment_down": [{"pathway": "REAC:REAC:R-HSA-9031628 NGF-stimulated transcription", "overlap": 2, "set_size": 40, "padj": 2e-3}],
    }

    out = hydrate_analysis_links(payload)

    assert out is not None
    assert out["top_up"][0]["gene_url"] == "https://www.ensembl.org/id/ENSG00000141510"
    assert "https://www.ensembl.org/Multi/Search/Results" in out["top_down"][0]["gene_url"]
    assert out["top_up"][0]["gene_symbol"] == "ENSG00000141510"
    assert out["top_up"][0]["gene_name"] == ""
    assert out["top_down"][0]["gene_symbol"] == "TP53"
    assert out["enrichment_up"][0]["pathway_url"] == "https://www.ebi.ac.uk/QuickGO/term/GO:0006954"
    assert out["enrichment_down"][0]["pathway_url"] == "https://reactome.org/content/detail/R-HSA-9031628"


def test_annotate_top_gene_rows_applies_fetched_ensembl_metadata(monkeypatch):
    monkeypatch.setenv("COMPBIO_FORCE_REMOTE_GENE_LOOKUP", "1")
    monkeypatch.setattr("app.analysis._load_gene_info_cache", lambda: {})

    written: dict = {}

    def _capture(entries, cache_path=None):
        written.update(entries)

    monkeypatch.setattr("app.analysis._write_gene_info_cache", _capture)
    monkeypatch.setattr(
        "app.analysis._fetch_ensembl_gene_annotations",
        lambda ensembl_ids, timeout_s=10.0: {
            "ENSG00000141510": {
                "ensembl_id": "ENSG00000141510",
                "gene_symbol": "TP53",
                "gene_name": "tumor protein p53",
                "gene_biotype": "protein_coding",
            }
        },
    )

    rows, note = _annotate_top_gene_rows(
        [{"gene": "ENSG00000141510.18", "log2fc": 1.2, "padj": 0.01}],
        allow_remote_lookup=True,
    )

    assert note is None
    assert rows[0]["gene_symbol"] == "TP53"
    assert rows[0]["gene_name"] == "tumor protein p53"
    assert rows[0]["gene_biotype"] == "protein_coding"
    assert rows[0]["gene_url"] == "https://www.ensembl.org/id/ENSG00000141510"
    assert "ENSG00000141510" in written
