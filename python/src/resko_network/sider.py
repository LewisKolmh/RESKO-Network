"""SIDER access and seed coverage -- ports A2 / A3 / A17A / A17B.

The R original matched seed drug names against SIDER with an unanchored
``str_detect`` (A2:35, A3:20). No false positive fires on the current seed list,
but the pattern is unsafe: ``"gam"`` matches both ``gamma-aminobutyric`` and
``sugammadex``. See CODE_REVIEW.md S3-1.

:func:`match_drug_names` matches on a normalised key with explicit boundaries,
and :func:`coverage_report` reports the resulting coverage with unmatchable
category labels excluded from the denominator.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

__all__ = [
    "normalise_name",
    "load_drug_names",
    "load_side_effects",
    "match_drug_names",
    "coverage_report",
]

#: Strings that mark a seed entry as a category label rather than a compound.
#: Such an entry can never match a database row, so counting it in a coverage
#: denominator understates coverage.
_CATEGORY_MARKERS = ("analogue", "analog", "derivative", "synthetic ")


def normalise_name(name: str) -> str:
    """Lowercase, collapse whitespace, and strip punctuation used inconsistently.

    >>> normalise_name("Ternatin-4")
    'ternatin 4'
    >>> normalise_name("  Didemnin  B ")
    'didemnin b'
    """
    text = str(name).lower().strip()
    text = re.sub(r"[-_/(),.]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def is_category_label(name: str) -> bool:
    """True when a seed entry names a compound class, not a compound.

    >>> is_category_label("Didemnin analogues")
    True
    >>> is_category_label("didemnin B")
    False
    """
    return any(marker in str(name).lower() for marker in _CATEGORY_MARKERS)


def load_drug_names(path: Path | str) -> pd.DataFrame:
    """Load SIDER ``drug_names.tsv`` (no header: STITCH flat id, name)."""
    return pd.read_csv(
        path, sep="\t", names=["stitch_flat_id", "drug_name"], dtype=str
    )


def load_side_effects(path: Path | str) -> pd.DataFrame:
    """Load SIDER ``meddra_all_se.tsv[.gz]``.

    Returns only preferred-term ("PT") rows, which is the granularity the
    RESKO side-effect signature comparison operates on.
    """
    path = Path(path)
    columns = [
        "stitch_flat_id",
        "stitch_stereo_id",
        "umls_cui_label",
        "meddra_concept_type",
        "umls_cui_meddra",
        "side_effect_name",
    ]
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        frame = pd.read_csv(handle, sep="\t", names=columns, dtype=str)
    return frame[frame["meddra_concept_type"] == "PT"].reset_index(drop=True)


def match_drug_names(
    seeds: list[str], drug_names: pd.DataFrame, *, name_column: str = "drug_name"
) -> pd.DataFrame:
    """Match seed names against a SIDER name table with word-boundary anchoring.

    A seed matches a SIDER row when its normalised form appears in the row's
    normalised name as a whole-token sequence -- so ``"gam"`` will not match
    ``"sugammadex"``, but ``"didemnin b"`` still matches
    ``"Didemnin B (NSC-325319)"``.
    """
    normalised = drug_names[name_column].map(normalise_name)
    rows: list[dict[str, object]] = []
    for seed in seeds:
        key = normalise_name(seed)
        pattern = rf"(?:^|\s){re.escape(key)}(?:\s|$)"
        hits = drug_names[normalised.str.contains(pattern, regex=True, na=False)]
        if hits.empty:
            rows.append(
                {
                    "seed": seed,
                    "normalised_seed": key,
                    "is_category_label": is_category_label(seed),
                    "matched": False,
                    "sider_drug_name": None,
                    "stitch_flat_id": None,
                }
            )
            continue
        for _, hit in hits.iterrows():
            rows.append(
                {
                    "seed": seed,
                    "normalised_seed": key,
                    "is_category_label": is_category_label(seed),
                    "matched": True,
                    "sider_drug_name": hit[name_column],
                    "stitch_flat_id": hit.get("stitch_flat_id"),
                }
            )
    return pd.DataFrame(rows)


def coverage_report(matches: pd.DataFrame) -> dict[str, object]:
    """Summarise seed coverage, excluding category labels from the denominator.

    The R original reported 0/15 for the expanded seed set. Three of those 15
    were category labels, so the honest figure is 0/12 -- the conclusion is
    unchanged, but the denominator is defensible.
    """
    per_seed = matches.groupby("seed").agg(
        matched=("matched", "any"),
        is_category_label=("is_category_label", "first"),
    )
    compounds = per_seed[~per_seed["is_category_label"]]
    return {
        "n_seeds_supplied": int(len(per_seed)),
        "n_category_labels_excluded": int(per_seed["is_category_label"].sum()),
        "n_compounds_searched": int(len(compounds)),
        "n_matched": int(compounds["matched"].sum()),
        "coverage_fraction": (
            float(compounds["matched"].mean()) if len(compounds) else 0.0
        ),
        "matched_names": sorted(
            matches.loc[matches["matched"], "sider_drug_name"].dropna().unique().tolist()
        ),
    }
