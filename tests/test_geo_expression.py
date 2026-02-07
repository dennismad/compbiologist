from __future__ import annotations

from pathlib import Path
import sys
import types

import pandas as pd
import pytest

from app.geo_expression import _compute_de, _iter_expression_sources, _parse_delimited_expression


def test_iter_expression_sources_falls_through_to_supplementary(monkeypatch):
    monkeypatch.setattr(
        "app.geo_expression._download_series_matrix_files",
        lambda gse, cache_dir: [Path("/tmp/fake_matrix.txt.gz")],
    )
    monkeypatch.setattr(
        "app.geo_expression._parse_series_matrix",
        lambda path: (_ for _ in ()).throw(ValueError("Series matrix table is empty")),
    )
    monkeypatch.setattr(
        "app.geo_expression._download_supplementary_files",
        lambda gse, cache_dir: [Path("/tmp/fake_suppl.tsv")],
    )
    monkeypatch.setattr(
        "app.geo_expression._read_text",
        lambda path: "gene\ts1_count\ts2_count\ts3_count\ts4_count\nA\t1\t2\t3\t4\nB\t3\t2\t1\t5\n",
    )

    rows = list(_iter_expression_sources("GSE_FAKE", Path("/tmp")))

    assert len(rows) == 1
    expr, sample_cols, _sample_text, source = rows[0]
    assert source == "supplementary"
    assert expr.shape[0] == 2
    assert len(sample_cols) == 4


def test_parse_delimited_expression_prefers_count_columns():
    text = (
        "gene\tS1_FPKM\tS2_FPKM\tS3_FPKM\tS4_FPKM\tS1_count\tS2_count\tS3_count\tS4_count\n"
        "A\t1.1\t2.2\t3.3\t4.4\t10\t20\t30\t40\n"
        "B\t2.1\t3.2\t4.3\t5.4\t11\t21\t31\t41\n"
    )

    expr, cols = _parse_delimited_expression(text)

    assert cols == ["S1_count", "S2_count", "S3_count", "S4_count"]
    assert isinstance(expr, pd.DataFrame)
    assert expr.shape == (2, 4)


def test_compute_de_requires_count_columns():
    expr = pd.DataFrame(
        {
            "A1_fpkm": [1.0, 2.0],
            "A2_fpkm": [1.2, 2.1],
            "B1_fpkm": [0.5, 1.5],
            "B2_fpkm": [0.6, 1.7],
        },
        index=["GENE1", "GENE2"],
    )

    with pytest.raises(ValueError, match="count matrix"):
        _compute_de(expr, ["A1_fpkm", "A2_fpkm"], ["B1_fpkm", "B2_fpkm"], "A", "B")


def test_compute_de_uses_pydeseq2_when_available(monkeypatch):
    class FakeDDS:
        def __init__(self, counts, metadata, design_factors, refit_cooks, **kwargs):
            self.counts = counts
            self.metadata = metadata
            self.design_factors = design_factors
            self.refit_cooks = refit_cooks
            self.kwargs = kwargs

        def deseq2(self):
            return None

    class FakeStats:
        def __init__(self, dds, contrast, **kwargs):
            self.dds = dds
            self.contrast = contrast
            self.kwargs = kwargs
            self.results_df = pd.DataFrame(
                {
                    "log2FoldChange": [1.2, -1.4],
                    "pvalue": [0.001, 0.02],
                    "padj": [0.01, 0.04],
                },
                index=["GENE1", "GENE2"],
            )

        def summary(self):
            return None

    fake_dds_mod = types.ModuleType("pydeseq2.dds")
    fake_ds_mod = types.ModuleType("pydeseq2.ds")
    fake_inf_mod = types.ModuleType("pydeseq2.default_inference")
    fake_pkg_mod = types.ModuleType("pydeseq2")
    fake_dds_mod.DeseqDataSet = FakeDDS
    fake_ds_mod.DeseqStats = FakeStats
    fake_inf_mod.DefaultInference = lambda **kwargs: object()

    monkeypatch.setitem(sys.modules, "pydeseq2", fake_pkg_mod)
    monkeypatch.setitem(sys.modules, "pydeseq2.dds", fake_dds_mod)
    monkeypatch.setitem(sys.modules, "pydeseq2.ds", fake_ds_mod)
    monkeypatch.setitem(sys.modules, "pydeseq2.default_inference", fake_inf_mod)

    expr = pd.DataFrame(
        {
            "A1_count": [10, 50],
            "A2_count": [11, 52],
            "B1_count": [4, 10],
            "B2_count": [5, 12],
        },
        index=["GENE1", "GENE2"],
    )

    out = _compute_de(expr, ["A1_count", "A2_count"], ["B1_count", "B2_count"], "A", "B")

    assert list(out.columns) == ["gene", "log2fc", "pvalue", "padj", "neg_log10_padj"]
    assert out.shape[0] == 2
    assert set(out["gene"]) == {"GENE1", "GENE2"}


def test_compute_de_accepts_count_like_values_without_count_suffix(monkeypatch):
    class FakeDDS:
        def __init__(self, counts, metadata, design_factors, refit_cooks, **kwargs):
            self.counts = counts
            self.metadata = metadata
            self.design_factors = design_factors
            self.refit_cooks = refit_cooks
            self.kwargs = kwargs

        def deseq2(self):
            return None

    class FakeStats:
        def __init__(self, dds, contrast, **kwargs):
            self.kwargs = kwargs
            self.results_df = pd.DataFrame(
                {
                    "log2FoldChange": [1.0],
                    "pvalue": [0.01],
                    "padj": [0.02],
                },
                index=["GENE1"],
            )

        def summary(self):
            return None

    fake_dds_mod = types.ModuleType("pydeseq2.dds")
    fake_ds_mod = types.ModuleType("pydeseq2.ds")
    fake_inf_mod = types.ModuleType("pydeseq2.default_inference")
    fake_pkg_mod = types.ModuleType("pydeseq2")
    fake_dds_mod.DeseqDataSet = FakeDDS
    fake_ds_mod.DeseqStats = FakeStats
    fake_inf_mod.DefaultInference = lambda **kwargs: object()

    monkeypatch.setitem(sys.modules, "pydeseq2", fake_pkg_mod)
    monkeypatch.setitem(sys.modules, "pydeseq2.dds", fake_dds_mod)
    monkeypatch.setitem(sys.modules, "pydeseq2.ds", fake_ds_mod)
    monkeypatch.setitem(sys.modules, "pydeseq2.default_inference", fake_inf_mod)

    expr = pd.DataFrame(
        {
            "Baseline_1": [100, 200],
            "Baseline_2": [110, 210],
            "Treat_1": [80, 150],
            "Treat_2": [85, 145],
        },
        index=["GENE1", "GENE2"],
    )

    out = _compute_de(expr, ["Treat_1", "Treat_2"], ["Baseline_1", "Baseline_2"], "Treat", "Baseline")
    assert out.shape[0] == 1
    assert out.iloc[0]["gene"] == "GENE1"
