"""ChEMBL access -- ports A8-A12, with a single canonical cache-key function.

The R original built cache filenames inline (A17D:570-576) and, across
revisions, used two different zero-padding conventions. The result was
``activities_CHEMBL941_0.json`` and ``activities_CHEMBL941_00000.json`` sitting
side by side with *different response schemas* -- 47 fields versus 35, with
``standard_upper_value`` present in one and absent from the other. See
CODE_REVIEW.md S2-1.

Here :func:`cache_key` is the only place a cache filename is constructed.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Iterator

import pandas as pd

from .config import CHEMBL

__all__ = [
    "cache_key",
    "ChemblClient",
    "activities_to_frame",
    "audit_cache_dir",
]


def cache_key(molecule_chembl_id: str, offset: int, kind: str = "activities") -> str:
    """Build a cache filename. The single source of truth for this format.

    >>> cache_key("CHEMBL941", 0)
    'activities_CHEMBL941_00000.json'
    >>> cache_key("CHEMBL941", 2000)
    'activities_CHEMBL941_02000.json'
    >>> cache_key("CHEMBL25", 0, kind="molecule")
    'molecule_CHEMBL25_00000.json'
    """
    width = CHEMBL.cache_offset_width
    return f"{kind}_{molecule_chembl_id}_{offset:0{width}d}.json"


def audit_cache_dir(cache_dir: Path | str) -> pd.DataFrame:
    """Report cache files whose names do not match the canonical convention.

    Use this to find orphaned files from an older padding convention before
    trusting a cache directory.

    Columns:

    ``paged``
        Whether the filename ends in a numeric offset token at all. Files that
        do not (e.g. ``molecule_alprazolam.json``) are not offset-paged and
        their padding is not meaningful -- ``canonical_name`` is NA for these.
    ``canonical_name``
        For paged files only: whether the offset is padded to the configured
        width.
    ``n_fields_union`` / ``n_fields_first`` / ``n_fields_common``
        Key counts across the payload's records: the union, the first record's,
        and the intersection. Divergence between union and common means records
        within one file carry different key sets -- typically because the
        serialiser omitted null-valued fields rather than emitting them
        explicitly. Comparing these across conventions is how a schema
        difference becomes visible.
    """
    cache_dir = Path(cache_dir)
    width = CHEMBL.cache_offset_width
    rows: list[dict[str, Any]] = []
    for path in sorted(cache_dir.glob("*.json")):
        parts = path.stem.rsplit("_", 1)
        offset_str = parts[1] if len(parts) == 2 else ""
        paged = offset_str.isdigit()
        record_key_sets: list[set[str]] = []
        try:
            payload = json.loads(path.read_text())
            records = payload.get("activities") or []
            record_key_sets = [set(r) for r in records if isinstance(r, dict)]
        except (json.JSONDecodeError, OSError):
            records = []
        rows.append(
            {
                "file": path.name,
                "paged": paged,
                "canonical_name": (len(offset_str) == width) if paged else pd.NA,
                "offset": int(offset_str) if paged else pd.NA,
                "n_records": len(record_key_sets),
                "n_fields_union": (
                    len(set().union(*record_key_sets)) if record_key_sets else 0
                ),
                "n_fields_first": len(record_key_sets[0]) if record_key_sets else 0,
                "n_fields_common": (
                    len(set.intersection(*record_key_sets)) if record_key_sets else 0
                ),
                "bytes": path.stat().st_size,
            }
        )
    return pd.DataFrame(rows)


class ChemblClient:
    """Paged ChEMBL client with an on-disk JSON cache.

    Caching writes the raw response bytes, matching the R original's
    ``writeLines(body, ..., useBytes = TRUE)`` -- that part of A17D was correct
    and is preserved deliberately.
    """

    def __init__(self, cache_dir: Path | str, *, offline: bool = False) -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.offline = offline

    def _get_json(self, url: str, cache_name: str) -> dict[str, Any]:
        path = self.cache_dir / cache_name
        if path.exists() and path.stat().st_size > 0:
            try:
                return json.loads(path.read_text())
            except json.JSONDecodeError:
                path.unlink()
        if self.offline:
            raise FileNotFoundError(
                f"offline=True and no cached response at {path}. "
                "Run once with offline=False to populate the cache."
            )
        import requests  # imported lazily so offline use needs no network stack

        last_error: Exception | str = "no attempt made"
        for attempt in range(CHEMBL.max_attempts):
            try:
                response = requests.get(
                    url,
                    headers={"User-Agent": CHEMBL.user_agent},
                    timeout=CHEMBL.timeout_s,
                )
                if 200 <= response.status_code < 300:
                    path.write_bytes(response.content)
                    time.sleep(CHEMBL.request_pause_s)
                    return json.loads(response.text)
                last_error = f"HTTP {response.status_code}"
            except Exception as exc:  # noqa: BLE001 - retried below
                last_error = exc
            if attempt < CHEMBL.max_attempts - 1:
                time.sleep(min(2**attempt, 16))
        raise RuntimeError(f"ChEMBL request failed for {url}. Last error: {last_error}")

    def activities(self, molecule_chembl_id: str) -> Iterator[dict[str, Any]]:
        """Yield every activity record for a molecule, following pagination."""
        offset = 0
        while True:
            url = (
                f"{CHEMBL.base_url}/activity.json"
                f"?molecule_chembl_id={molecule_chembl_id}"
                f"&limit={CHEMBL.page_limit}&offset={offset}"
            )
            payload = self._get_json(
                url, cache_key(molecule_chembl_id, offset, "activities")
            )
            records = payload.get("activities") or []
            if not records:
                return
            yield from records
            if len(records) < CHEMBL.page_limit:
                return
            offset += CHEMBL.page_limit


#: Fields carried through to the analysis frame. ``standard_relation`` is
#: included explicitly -- it is the field the R ranking dropped.
ACTIVITY_FIELDS: tuple[str, ...] = (
    "activity_id",
    "molecule_chembl_id",
    "target_chembl_id",
    "target_pref_name",
    "target_organism",
    "assay_chembl_id",
    "assay_description",
    "assay_type",
    "document_chembl_id",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "standard_upper_value",
    "pchembl_value",
    "action_type",
    "data_validity_comment",
    "potential_duplicate",
)


def activities_to_frame(records: list[dict[str, Any]]) -> pd.DataFrame:
    """Flatten activity records to a frame with a stable column set."""
    frame = pd.DataFrame(records)
    for field in ACTIVITY_FIELDS:
        if field not in frame.columns:
            frame[field] = None
    return frame[list(ACTIVITY_FIELDS)]
