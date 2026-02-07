from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass
class ProcessedResult:
    dataframe: pd.DataFrame
    summary: dict


def process_protein_table(raw_tsv_path: str) -> ProcessedResult:
    df = pd.read_csv(raw_tsv_path, sep="\t")

    rename_map = {
        "Entry": "accession",
        "Entry Name": "entry_name",
        "Protein names": "protein_name",
        "Gene Names": "gene_names",
        "Length": "length",
        "Annotation": "annotation_score",
    }
    df = df.rename(columns=rename_map)

    required_cols = ["accession", "protein_name", "gene_names", "length", "annotation_score"]
    for col in required_cols:
        if col not in df.columns:
            df[col] = ""

    df["length"] = pd.to_numeric(df["length"], errors="coerce")
    df["annotation_score"] = pd.to_numeric(df["annotation_score"], errors="coerce")
    df = df.dropna(subset=["length"]).copy()
    df["length"] = df["length"].astype(int)

    bins = [0, 300, 600, 900, 1200, 2000]
    labels = ["1-300", "301-600", "601-900", "901-1200", "1201-2000"]
    df["length_bin"] = pd.cut(df["length"], bins=bins, labels=labels, include_lowest=True)

    distribution = (
        df["length_bin"]
        .value_counts(dropna=False)
        .reindex(labels, fill_value=0)
        .astype(int)
        .to_dict()
    )

    top_genes = (
        df["gene_names"]
        .fillna("")
        .astype(str)
        .str.split()
        .str[0]
        .replace("", pd.NA)
        .dropna()
        .value_counts()
        .head(10)
        .to_dict()
    )

    summary = {
        "n_proteins": int(df.shape[0]),
        "avg_length": round(float(df["length"].mean()), 1) if not df.empty else 0,
        "median_length": int(df["length"].median()) if not df.empty else 0,
        "avg_annotation_score": round(float(df["annotation_score"].mean()), 2) if not df.empty else 0,
        "length_distribution": distribution,
        "top_genes": top_genes,
    }

    display_cols = ["accession", "protein_name", "gene_names", "length", "annotation_score"]
    return ProcessedResult(dataframe=df[display_cols], summary=summary)
