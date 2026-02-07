from __future__ import annotations

import csv
import gzip
import io
import os
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests

from app.config import DATA_DIR
from app.geo_cache import get_cached_analyzable, set_cached_analyzable


@dataclass
class RealDEResult:
    gse: str
    group_a_name: str
    group_b_name: str
    group_a_samples: list[str]
    group_b_samples: list[str]
    dataframe: pd.DataFrame


def _extract_item_gse(item: dict) -> str:
    for key in ("accession", "gse"):
        value = str(item.get(key, "")).strip()
        if value.startswith("GSE"):
            return value
    return ""


def _series_prefix(gse: str) -> str:
    num = int(gse.replace("GSE", ""))
    return f"GSE{num // 1000}nnn"


def _extract_gse_candidates(loaded_items: list[dict]) -> list[str]:
    out: list[str] = []
    for row in loaded_items:
        value = _extract_item_gse(row)
        if value and value not in out:
            out.append(value)
    return out


def _split_geo_line(line: str) -> list[str]:
    reader = csv.reader([line], delimiter="\t", quotechar='"')
    return next(reader)


def _resolve_listing_urls(listing_url: str, pattern: str) -> list[str]:
    resp = requests.get(listing_url, timeout=25)
    resp.raise_for_status()
    names = re.findall(pattern, resp.text)
    urls = []
    for name in names:
        if name.startswith("http"):
            urls.append(name)
        else:
            urls.append(listing_url + name)
    return urls


def _download_url(url: str, out_path: Path, timeout: int = 120) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
    resp = requests.get(url, timeout=timeout)
    resp.raise_for_status()
    out_path.write_bytes(resp.content)
    return out_path


def _resolve_series_matrix_urls(gse: str) -> list[str]:
    prefix = _series_prefix(gse)
    listing_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/"
    return _resolve_listing_urls(listing_url, r'href="([^\"]+_series_matrix\.txt\.gz)"')


def _download_series_matrix_files(gse: str, cache_dir: Path) -> list[Path]:
    local = sorted((cache_dir / gse).glob(f"{gse}_series_matrix_*.txt.gz"))
    local = [p for p in local if p.is_file() and p.stat().st_size > 0]
    if local:
        return local

    urls = _resolve_series_matrix_urls(gse)
    paths = []
    for i, url in enumerate(urls):
        name = f"{gse}_series_matrix_{i}.txt.gz"
        paths.append(_download_url(url, cache_dir / gse / name, timeout=120))
    return paths


def _resolve_supplementary_urls(gse: str) -> list[str]:
    prefix = _series_prefix(gse)
    listing_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/suppl/"
    urls = _resolve_listing_urls(listing_url, r'href="([^\"]+)"')

    keep = []
    for url in urls:
        low = url.lower()
        if any(x in low for x in ["matrix", "count", "counts", "expr", "expression", "fpkm", "tpm"]) and any(
            low.endswith(ext) for ext in [".txt", ".txt.gz", ".tsv", ".tsv.gz", ".csv", ".csv.gz"]
        ):
            keep.append(url)
    return keep


def _download_supplementary_files(gse: str, cache_dir: Path) -> list[Path]:
    local = sorted((cache_dir / gse).glob(f"{gse}_suppl_*"))
    local = [p for p in local if p.is_file() and p.stat().st_size > 0]
    if local:
        return local

    paths: list[Path] = []
    for i, url in enumerate(_resolve_supplementary_urls(gse)):
        suffix = ".txt"
        for ext in [".txt.gz", ".tsv.gz", ".csv.gz", ".txt", ".tsv", ".csv"]:
            if url.lower().endswith(ext):
                suffix = ext
                break
        name = f"{gse}_suppl_{i}{suffix}"
        try:
            paths.append(_download_url(url, cache_dir / gse / name, timeout=180))
        except Exception:
            continue
    return paths


def _read_text(path: Path) -> str:
    if path.suffix == ".gz" or str(path).endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read()
    return path.read_text(encoding="utf-8", errors="replace")


def _parse_delimited_expression(text: str) -> tuple[pd.DataFrame, list[str]]:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 3:
        raise ValueError("File has too few lines")

    sample = "\n".join(lines[:200])
    sep = "\t" if sample.count("\t") >= sample.count(",") else ","

    df = pd.read_csv(io.StringIO("\n".join(lines)), sep=sep, engine="python", on_bad_lines="skip")
    if df.empty or df.shape[1] < 4:
        raise ValueError("Expression candidate did not parse into a valid table")

    id_col = df.columns[0]
    candidate_cols = [c for c in df.columns if c != id_col]

    numeric_cols = []
    for col in candidate_cols:
        numeric = pd.to_numeric(df[col], errors="coerce")
        frac = float(numeric.notna().mean())
        if frac >= 0.6:
            numeric_cols.append(col)

    if len(numeric_cols) < 4:
        raise ValueError("Not enough numeric sample columns")

    # Prefer one coherent assay family instead of mixing count + FPKM/TPM columns.
    lower_cols = {c: str(c).lower() for c in numeric_cols}
    assay_priority = [
        ("count", lambda s: "count" in s),
        ("tpm", lambda s: "tpm" in s),
        ("fpkm", lambda s: "fpkm" in s),
        ("rpkm", lambda s: "rpkm" in s),
        ("cpm", lambda s: "cpm" in s),
    ]
    for _, matcher in assay_priority:
        subset = [c for c in numeric_cols if matcher(lower_cols[c])]
        if len(subset) >= 4:
            numeric_cols = subset
            break

    out = df[[id_col] + numeric_cols].copy()
    out[id_col] = out[id_col].astype(str)
    for col in numeric_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    out = out.dropna(subset=numeric_cols, how="all")
    if out.empty:
        raise ValueError("Expression table became empty after numeric filtering")

    return out.set_index(id_col), numeric_cols


def _parse_series_matrix(path: Path) -> tuple[pd.DataFrame, list[str], dict[str, list[str]]]:
    sample_meta: dict[str, list[str]] = {}
    table_lines: list[str] = []
    in_table = False

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break

            if in_table:
                table_lines.append(line)
            elif line.startswith("!Sample_"):
                parts = _split_geo_line(line)
                key = parts[0].replace("!", "")
                sample_meta[key] = [x.strip().strip('"') for x in parts[1:]]

    if len(table_lines) < 3:
        raise ValueError("Series matrix table is empty")

    expr, sample_cols = _parse_delimited_expression("\n".join(table_lines))
    return expr, sample_cols, sample_meta


def _sample_annotation_map(sample_ids: list[str], sample_meta: dict[str, list[str]]) -> dict[str, str]:
    meta_ids = sample_meta.get("Sample_geo_accession", sample_ids)
    idx_map = {sid: i for i, sid in enumerate(meta_ids)}

    annotation_fields = [k for k in sample_meta if k.startswith("Sample_characteristics")]
    out: dict[str, str] = {}
    for sid in sample_ids:
        i = idx_map.get(sid)
        parts: list[str] = [sid]

        for field in ["Sample_title", "Sample_source_name_ch1"] + annotation_fields:
            values = sample_meta.get(field, [])
            if i is not None and i < len(values):
                val = str(values[i]).strip()
                if val:
                    parts.append(val)
        out[sid] = " | ".join(parts).lower()

    return out


def _assign_groups(sample_text: dict[str, str], state_profile: str) -> tuple[list[str], list[str], str, str]:
    disease_tokens = ["disease", "patient", "tumor", "cancer", "diabetes", "case"]
    healthy_tokens = ["healthy", "normal", "control", "non-disease"]
    treated_tokens = ["treated", "treatment", "drug", "therapy", "stimulated"]
    untreated_tokens = ["untreated", "vehicle", "placebo", "mock", "baseline", "control"]

    disease = [s for s, t in sample_text.items() if any(tok in t for tok in disease_tokens)]
    healthy = [s for s, t in sample_text.items() if any(tok in t for tok in healthy_tokens)]
    treated = [s for s, t in sample_text.items() if any(tok in t for tok in treated_tokens)]
    untreated = [s for s, t in sample_text.items() if any(tok in t for tok in untreated_tokens)]

    if state_profile == "Disease vs Healthy":
        a, b, an, bn = disease, healthy, "Disease", "Healthy"
    elif state_profile == "Treated vs Untreated":
        a, b, an, bn = treated, untreated, "Treated", "Untreated"
    elif state_profile == "Disease only":
        rest = [s for s in sample_text if s not in set(disease)]
        a, b, an, bn = disease, rest, "Disease", "Other"
    elif state_profile == "Healthy only":
        rest = [s for s in sample_text if s not in set(healthy)]
        a, b, an, bn = healthy, rest, "Healthy", "Other"
    else:
        a, b, an, bn = disease, healthy, "Disease", "Healthy"

    a_set = list(dict.fromkeys(a))
    b_set = [x for x in dict.fromkeys(b) if x not in set(a_set)]
    return a_set, b_set, an, bn


def _normalize_gene_id(raw: str) -> str:
    token = str(raw).strip()
    for sep in [" /// ", "//", ";", "|", ",", " "]:
        if sep in token:
            token = token.split(sep)[0]
            break
    token = token.strip().upper()
    token = re.sub(r"[^A-Z0-9_.-]", "", token)
    return token or "UNKNOWN"


def _bh_adjust(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    qvals = [1.0] * n
    min_coeff = 1.0
    for rank in range(n, 0, -1):
        idx, p = indexed[rank - 1]
        coeff = min(min_coeff, p * n / rank)
        min_coeff = coeff
        qvals[idx] = min(1.0, coeff)
    return qvals


def _is_count_like_column(name: str) -> bool:
    return "count" in str(name).strip().lower()


def _is_count_like_matrix(expr: pd.DataFrame, sample_cols: list[str]) -> bool:
    if not sample_cols:
        return False

    sub = expr[sample_cols].copy()
    for col in sample_cols:
        sub[col] = pd.to_numeric(sub[col], errors="coerce")

    values = sub.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return False

    frac_nonneg = float(np.mean(finite >= 0))
    frac_int = float(np.mean(np.abs(finite - np.round(finite)) < 1e-6))
    return frac_nonneg >= 0.99 and frac_int >= 0.95


def _coerce_gene_count_matrix(expr: pd.DataFrame, sample_cols: list[str]) -> pd.DataFrame:
    out = expr[sample_cols].copy()
    for col in sample_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.fillna(0.0).clip(lower=0.0)
    out.index = [_normalize_gene_id(x) for x in out.index]
    out = out.groupby(level=0).sum()
    out = out[(out >= 10).sum(axis=1) >= 2]
    if out.empty:
        raise ValueError("No genes passed minimal count filter (>=10 reads in at least 2 samples)")
    return out


def _load_pydeseq2() -> tuple[Any, Any, Any]:
    mpl_dir = DATA_DIR / "processed" / ".mplconfig"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
    try:
        from pydeseq2.dds import DeseqDataSet  # type: ignore
        from pydeseq2.ds import DeseqStats  # type: ignore
        from pydeseq2.default_inference import DefaultInference  # type: ignore
    except Exception as exc:
        raise ImportError(
            "pydeseq2 is required for real DE analysis. Install dependencies and rerun."
        ) from exc

    # pydeseq2 passes `inner_max_num_threads` to joblib.parallel_backend for all
    # backends. With "threading", some joblib versions reject that kwarg.
    # Patch once to drop it only for threading backend.
    infer_mod = sys.modules.get(getattr(DefaultInference, "__module__", ""))
    if infer_mod and not getattr(infer_mod, "_compbio_safe_backend_patch", False):
        original_parallel_backend = getattr(infer_mod, "parallel_backend", None)
        if callable(original_parallel_backend):
            def _safe_parallel_backend(backend: str, *args, **kwargs):
                if backend == "threading":
                    kwargs.pop("inner_max_num_threads", None)
                return original_parallel_backend(backend, *args, **kwargs)

            setattr(infer_mod, "parallel_backend", _safe_parallel_backend)
            setattr(infer_mod, "_compbio_safe_backend_patch", True)
    return DeseqDataSet, DeseqStats, DefaultInference


def _compute_de(
    expr: pd.DataFrame,
    group_a: list[str],
    group_b: list[str],
    group_a_name: str,
    group_b_name: str,
) -> pd.DataFrame:
    sample_cols = list(group_a) + list(group_b)
    has_count_names = all(_is_count_like_column(col) for col in sample_cols)
    if not has_count_names and not _is_count_like_matrix(expr, sample_cols):
        raise ValueError("Real DE requires raw count matrix columns or integer count-like values.")

    gene_counts = _coerce_gene_count_matrix(expr, sample_cols=sample_cols)
    counts_df = gene_counts[sample_cols].T.round().astype(int)
    counts_df = counts_df.loc[:, counts_df.sum(axis=0) > 0]
    if counts_df.shape[1] == 0:
        raise ValueError("No non-zero genes remain after filtering.")

    metadata = pd.DataFrame(
        {
            "condition": [group_a_name] * len(group_a) + [group_b_name] * len(group_b),
        },
        index=sample_cols,
    )

    DeseqDataSet, DeseqStats, DefaultInference = _load_pydeseq2()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        inference = DefaultInference(n_cpus=1, backend="threading")
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design_factors="condition",
            refit_cooks=True,
            n_cpus=1,
            inference=inference,
            quiet=True,
            low_memory=True,
        )
        dds.deseq2()
        stat_res = DeseqStats(
            dds,
            contrast=["condition", group_a_name, group_b_name],
            inference=inference,
            n_cpus=1,
        )
        stat_res.summary()

    result_df = stat_res.results_df.copy()
    if result_df.empty:
        raise ValueError("DESeq2 returned no results.")

    result_df = result_df.reset_index().rename(
        columns={
            result_df.index.name or "index": "gene",
            "log2FoldChange": "log2fc",
        }
    )
    if "gene" not in result_df.columns:
        result_df = result_df.rename(columns={"index": "gene"})
    for col in ("log2fc", "pvalue", "padj"):
        if col not in result_df.columns:
            result_df[col] = np.nan

    result_df = result_df[["gene", "log2fc", "pvalue", "padj"]].copy()
    result_df["gene"] = result_df["gene"].astype(str).map(_normalize_gene_id)
    result_df["log2fc"] = pd.to_numeric(result_df["log2fc"], errors="coerce").fillna(0.0)
    result_df["pvalue"] = pd.to_numeric(result_df["pvalue"], errors="coerce").fillna(1.0).clip(lower=1e-300, upper=1.0)
    result_df["padj"] = pd.to_numeric(result_df["padj"], errors="coerce").fillna(1.0).clip(lower=1e-300, upper=1.0)
    result_df["neg_log10_padj"] = -np.log10(result_df["padj"])

    if result_df.empty:
        raise ValueError("No DE rows could be computed.")
    return result_df


def _iter_expression_sources(gse: str, cache_dir: Path):
    try:
        matrix_paths = _download_series_matrix_files(gse, cache_dir=cache_dir)
    except Exception:
        matrix_paths = []
    for matrix_path in matrix_paths:
        try:
            expr, sample_ids, sample_meta = _parse_series_matrix(matrix_path)
            sample_text = _sample_annotation_map(sample_ids, sample_meta)
            yield expr, list(expr.columns), sample_text, "series_matrix"
        except Exception:
            continue

    try:
        supplementary_paths = _download_supplementary_files(gse, cache_dir=cache_dir)
    except Exception:
        supplementary_paths = []
    for path in supplementary_paths:
        try:
            text = _read_text(path)
            expr, sample_cols = _parse_delimited_expression(text)
            sample_text = {sid: sid.lower() for sid in sample_cols}
            yield expr, sample_cols, sample_text, "supplementary"
        except Exception:
            continue


def _validate_manual_groups(sample_cols: list[str], group_a_samples: list[str], group_b_samples: list[str]) -> tuple[list[str], list[str]]:
    group_a = [x for x in group_a_samples if x in sample_cols]
    group_b = [x for x in group_b_samples if x in sample_cols and x not in set(group_a)]
    if len(group_a) < 2 or len(group_b) < 2:
        raise ValueError(f"Manual groups need >=2 samples per group; got {len(group_a)} vs {len(group_b)}")
    return group_a, group_b


def build_group_choice_context(loaded_items: list[dict], state_profile: str, cache_dir: Path) -> tuple[dict | None, list[str]]:
    errors: list[str] = []

    for gse in _extract_gse_candidates(loaded_items):
        try:
            for expr, sample_cols, sample_text, source in _iter_expression_sources(gse, cache_dir=cache_dir):
                suggested_a, suggested_b, group_a_name, group_b_name = _assign_groups(sample_text, state_profile)
                suggested_a = [x for x in suggested_a if x in sample_cols]
                suggested_b = [x for x in suggested_b if x in sample_cols and x not in set(suggested_a)]

                options = [{"id": sid, "label": sample_text.get(sid, sid)} for sid in sample_cols]
                return {
                    "gse": gse,
                    "source": source,
                    "sample_options": options,
                    "suggested_group_a": suggested_a,
                    "suggested_group_b": suggested_b,
                    "group_a_name": group_a_name,
                    "group_b_name": group_b_name,
                }, errors
        except Exception as exc:
            errors.append(f"{gse}: {exc}")

    return None, errors


def run_real_geo_de(
    loaded_items: list[dict],
    state_profile: str,
    cache_dir: Path,
    manual_choice: dict | None = None,
) -> tuple[RealDEResult | None, list[str]]:
    errors: list[str] = []

    candidates = _extract_gse_candidates(loaded_items)
    if manual_choice and manual_choice.get("gse"):
        selected = str(manual_choice["gse"])
        candidates = [selected] + [x for x in candidates if x != selected]

    for gse in candidates:
        gse_errors: list[str] = []
        try:
            for expr, sample_cols, sample_text, source in _iter_expression_sources(gse, cache_dir=cache_dir):
                if manual_choice and manual_choice.get("gse") == gse:
                    group_a, group_b = _validate_manual_groups(
                        sample_cols=sample_cols,
                        group_a_samples=manual_choice.get("group_a_samples", []),
                        group_b_samples=manual_choice.get("group_b_samples", []),
                    )
                    group_a_name = str(manual_choice.get("group_a_name", "Group A"))
                    group_b_name = str(manual_choice.get("group_b_name", "Group B"))
                else:
                    group_a, group_b, group_a_name, group_b_name = _assign_groups(sample_text, state_profile)
                    group_a = [x for x in group_a if x in sample_cols]
                    group_b = [x for x in group_b if x in sample_cols and x not in set(group_a)]
                    if len(group_a) < 2 or len(group_b) < 2:
                        raise ValueError(f"Insufficient samples for {state_profile}: {len(group_a)} vs {len(group_b)}")

                df = _compute_de(
                    expr=expr,
                    group_a=group_a,
                    group_b=group_b,
                    group_a_name=group_a_name,
                    group_b_name=group_b_name,
                )
                return RealDEResult(
                    gse=gse,
                    group_a_name=group_a_name,
                    group_b_name=group_b_name,
                    group_a_samples=group_a,
                    group_b_samples=group_b,
                    dataframe=df,
                ), errors
        except Exception as exc:
            gse_errors.append(str(exc))

        errors.append(f"{gse}: {' | '.join(gse_errors) or 'no usable expression matrix found'}")

    return None, errors


def assess_gse_analyzable(gse: str, cache_dir: Path) -> tuple[bool, str]:
    cached = get_cached_analyzable(gse)
    if cached is not None:
        return cached

    errors: list[str] = []

    try:
        for matrix_path in _download_series_matrix_files(gse, cache_dir=cache_dir):
            try:
                expr, sample_ids, sample_meta = _parse_series_matrix(matrix_path)
                sample_text = _sample_annotation_map(sample_ids, sample_meta)
                if expr.shape[0] > 0 and len(sample_text) >= 4:
                    set_cached_analyzable(gse, True, "series_matrix")
                    return True, "series_matrix"
            except Exception as exc:
                errors.append(f"series_matrix: {exc}")
    except Exception as exc:
        errors.append(f"series_matrix: {exc}")

    try:
        for path in _download_supplementary_files(gse, cache_dir=cache_dir):
            try:
                text = _read_text(path)
                expr, sample_cols = _parse_delimited_expression(text)
                if expr.shape[0] > 0 and len(sample_cols) >= 4:
                    set_cached_analyzable(gse, True, "supplementary")
                    return True, "supplementary"
            except Exception as exc:
                errors.append(f"supplementary: {exc}")
    except Exception as exc:
        errors.append(f"supplementary: {exc}")

    reason = " | ".join(dict.fromkeys(errors)) if errors else "No parseable matrix found"
    set_cached_analyzable(gse, False, reason)
    return False, reason


def assess_gse_analyzable_fast(gse: str, cache_dir: Path, row: dict | None = None) -> tuple[bool, str]:
    cached = get_cached_analyzable(gse)
    if cached is not None:
        return cached

    local_matrix = sorted((cache_dir / gse).glob(f"{gse}_series_matrix_*.txt.gz"))
    if any(p.is_file() and p.stat().st_size > 0 for p in local_matrix):
        return True, "series_matrix_cached"

    local_suppl = sorted((cache_dir / gse).glob(f"{gse}_suppl_*"))
    if any(p.is_file() and p.stat().st_size > 0 for p in local_suppl):
        return True, "supplementary_cached"

    # Search-time check is intentionally metadata-only to keep latency low.
    experiment = str((row or {}).get("experiment_type", "")).strip()
    n_samples = int(pd.to_numeric((row or {}).get("n_samples", 0), errors="coerce") or 0)
    if experiment in {"RNA-seq", "Single-cell RNA-seq"} and n_samples >= 4:
        return True, "fast_metadata_predicted"
    return False, "fast_metadata_not_analyzable"


def annotate_items_analyzable(
    items: list[dict],
    cache_dir: Path,
    check_limit: int = 40,
    strict_validation: bool = True,
) -> list[dict]:
    out: list[dict] = []
    cache: dict[str, tuple[bool, str]] = {}

    for idx, item in enumerate(items):
        row = dict(item)
        gse = _extract_item_gse(row)
        if not gse:
            row["analyzable"] = False
            row["analyzable_detail"] = "No GSE accession"
            out.append(row)
            continue

        if gse not in cache:
            if idx < check_limit:
                if strict_validation:
                    cache[gse] = assess_gse_analyzable(gse, cache_dir=cache_dir)
                else:
                    cache[gse] = assess_gse_analyzable_fast(gse, cache_dir=cache_dir, row=row)
            else:
                cache[gse] = (False, "unchecked_not_validated")

        ok, detail = cache[gse]
        row["analyzable"] = bool(ok)
        row["analyzable_detail"] = detail
        out.append(row)

    return out
