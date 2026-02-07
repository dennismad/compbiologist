from __future__ import annotations

import json
from hashlib import sha256
from datetime import datetime, timedelta, timezone
from pathlib import Path

from app.config import GEO_CACHE_DB_PATH


def _utcnow() -> datetime:
    return datetime.now(timezone.utc)


def _cache_path(path: Path | None = None) -> Path:
    db_path = path or GEO_CACHE_DB_PATH
    db_path.parent.mkdir(parents=True, exist_ok=True)
    return db_path


def _load_cache(path: Path | None = None) -> dict:
    cache_path = _cache_path(path)
    if not cache_path.exists():
        return {}
    try:
        return json.loads(cache_path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _write_cache(payload: dict, path: Path | None = None) -> None:
    cache_path = _cache_path(path)
    cache_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _search_key(
    query: str,
    species_filter: str,
    experiment_filter: str,
    state_filter: str,
    retmax: int,
    retstart: int,
) -> str:
    cache_version = "v10"
    raw = "|".join(
        [
            cache_version,
            query.strip().lower(),
            species_filter.strip().lower(),
            experiment_filter.strip().lower(),
            state_filter.strip().lower(),
            str(int(retmax)),
            str(int(retstart)),
        ]
    )
    return sha256(raw.encode("utf-8")).hexdigest()


def get_cached_search(
    query: str,
    species_filter: str,
    experiment_filter: str,
    state_filter: str,
    retmax: int,
    retstart: int,
    ttl_minutes: int = 120,
) -> dict | None:
    cache = _load_cache()
    records = cache.get("search", {})
    key = _search_key(
        query=query,
        species_filter=species_filter,
        experiment_filter=experiment_filter,
        state_filter=state_filter,
        retmax=retmax,
        retstart=retstart,
    )
    rec = records.get(key)
    if not rec:
        return None

    updated_raw = rec.get("updated_at")
    if not updated_raw:
        return None

    try:
        updated = datetime.fromisoformat(updated_raw)
        if updated.tzinfo is None:
            updated = updated.replace(tzinfo=timezone.utc)
    except Exception:
        return None

    if _utcnow() - updated > timedelta(minutes=ttl_minutes):
        return None

    payload = rec.get("payload")
    return payload if isinstance(payload, dict) else None


def set_cached_search(
    query: str,
    species_filter: str,
    experiment_filter: str,
    state_filter: str,
    retmax: int,
    retstart: int,
    payload: dict,
) -> None:
    cache = _load_cache()
    cache.setdefault("search", {})
    key = _search_key(
        query=query,
        species_filter=species_filter,
        experiment_filter=experiment_filter,
        state_filter=state_filter,
        retmax=retmax,
        retstart=retstart,
    )
    cache["search"][key] = {
        "updated_at": _utcnow().isoformat(),
        "payload": payload,
    }
    _write_cache(cache)


def get_cached_analyzable(gse: str, ttl_days: int = 30) -> tuple[bool, str] | None:
    cache = _load_cache()
    records = cache.get("analyzable", {})
    rec = records.get(gse)
    if not rec:
        return None

    updated_raw = rec.get("updated_at")
    if not updated_raw:
        return None

    try:
        updated = datetime.fromisoformat(updated_raw)
        if updated.tzinfo is None:
            updated = updated.replace(tzinfo=timezone.utc)
    except Exception:
        return None

    if _utcnow() - updated > timedelta(days=ttl_days):
        return None

    return bool(rec.get("ok", False)), str(rec.get("detail", ""))


def set_cached_analyzable(gse: str, ok: bool, detail: str) -> None:
    cache = _load_cache()
    cache.setdefault("analyzable", {})
    cache["analyzable"][gse] = {
        "ok": bool(ok),
        "detail": str(detail)[:2000],
        "updated_at": _utcnow().isoformat(),
    }
    _write_cache(cache)
