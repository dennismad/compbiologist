from __future__ import annotations

from app.processing import process_protein_table


def test_process_protein_table_summary_and_columns(tmp_path):
    tsv = "\n".join(
        [
            "Entry\tEntry Name\tProtein names\tGene Names\tLength\tAnnotation",
            "P1\tA_HUMAN\tAlpha\tGENEA\t100\t5",
            "P2\tB_HUMAN\tBeta\tGENEB\t450\t4",
            "P3\tC_HUMAN\tGamma\tGENEA GENEA2\t800\t3",
            "P4\tD_HUMAN\tDelta\t\tbad\t2",
        ]
    )
    raw_path = tmp_path / "proteins.tsv"
    raw_path.write_text(tsv, encoding="utf-8")

    result = process_protein_table(str(raw_path))

    assert list(result.dataframe.columns) == [
        "accession",
        "protein_name",
        "gene_names",
        "length",
        "annotation_score",
    ]
    assert result.summary["n_proteins"] == 3
    assert result.summary["avg_length"] == 450.0
    assert result.summary["median_length"] == 450
    assert result.summary["top_genes"]["GENEA"] == 2
    assert result.summary["length_distribution"]["1-300"] == 1
    assert result.summary["length_distribution"]["301-600"] == 1
    assert result.summary["length_distribution"]["601-900"] == 1
