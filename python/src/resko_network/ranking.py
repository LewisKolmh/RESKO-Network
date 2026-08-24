"""Candidate ranking -- ports ``A15_rank_candidates.R`` with relation-aware scoring.

The R original scored potency from the numeric activity value alone, discarding
ChEMBL's ``standard_relation`` operator. A record of ``IC50 >= 10000 nM`` --
meaning no inhibition was observed up to 10 uM -- therefore scored the same as a
genuine 9 uM measurement. See CODE_REVIEW.md S1-1.

This module keeps the R scoring ladder identical but makes two additions:

1. ``potency_score`` inspects the relation and returns 0 for censored records.
2. Records whose assay description indicates a proteomics/pulldown format do not
   qualify for an inhibition-strength bonus.

Both additions are individually switchable so the original behaviour can be
reproduced for comparison -- see :func:`rank_candidates`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .config import (
    CENSORED_RELATIONS,
    EEF1_COMPLEX,
    EEF1A_PARALOGUES,
    EVIDENCE,
    POTENCY,
    PROTEOMICS_ASSAY_PATTERNS,
)

__all__ = [
    "ScoringMode",
    "is_censored",
    "is_proteomics_assay",
    "potency_score",
    "evidence_score",
    "development_score",
    "score_evidence_table",
    "rank_candidates",
]


@dataclass(frozen=True)
class ScoringMode:
    """Toggles for the two corrections, so before/after is reproducible.

    ``ScoringMode.legacy()`` reproduces the R original's behaviour exactly;
    ``ScoringMode()`` (the default) applies both corrections.
    """

    relation_aware: bool = True
    demote_proteomics: bool = True

    @classmethod
    def legacy(cls) -> ScoringMode:
        """The R original: relation ignored, assay format ignored."""
        return cls(relation_aware=False, demote_proteomics=False)


# --------------------------------------------------------------------------
# Record-level predicates
# --------------------------------------------------------------------------


def is_censored(relation: object) -> bool:
    """True when a ChEMBL ``standard_relation`` marks a censored measurement.

    Censored means the assay found no activity up to the tested ceiling, so the
    numeric value is a bound rather than a measurement.

    >>> is_censored(">=")
    True
    >>> is_censored("=")
    False
    >>> is_censored(None)
    False
    """
    if relation is None or (isinstance(relation, float) and np.isnan(relation)):
        return False
    return str(relation).strip() in CENSORED_RELATIONS


def is_proteomics_assay(description: object) -> bool:
    """True when an assay description indicates a pulldown/proteomics format.

    >>> is_proteomics_assay("Inhibition of EEF1A1 by colloidal coomassie LC-MS/MS")
    True
    >>> is_proteomics_assay("Inhibition of recombinant EEF1A1 GTPase activity")
    False
    """
    if description is None or (
        isinstance(description, float) and np.isnan(description)
    ):
        return False
    text = str(description).lower()
    return any(pattern in text for pattern in PROTEOMICS_ASSAY_PATTERNS)


# --------------------------------------------------------------------------
# Scoring functions (ported from A15:92-120)
# --------------------------------------------------------------------------


def potency_score(
    value_nM: object,
    relation: object = None,
    *,
    relation_aware: bool = True,
) -> int:
    """Score an activity value on the A15 potency ladder.

    Ports ``potency_score`` (A15:101-110), adding the relation check. With
    ``relation_aware=False`` the behaviour is byte-for-byte the R original's.

    >>> potency_score(3.74, "=")
    4
    >>> potency_score(240.0, ">=")            # censored -> no credit
    0
    >>> potency_score(240.0, ">=", relation_aware=False)   # the R behaviour
    2
    """
    if relation_aware and is_censored(relation):
        return 0
    if value_nM is None or (isinstance(value_nM, float) and np.isnan(value_nM)):
        return 0
    try:
        value = float(value_nM)
    except (TypeError, ValueError):
        return 0
    for threshold, score in POTENCY.breakpoints:
        if value <= threshold:
            return score
    return 0


def evidence_score(
    relationship: object,
    assay_description: object = None,
    *,
    demote_proteomics: bool = True,
) -> int:
    """Score a relationship class.

    Ports ``evidence_score`` (A15:92-99). When ``demote_proteomics`` is set, an
    inhibition claim resting on a pulldown/proteomics assay is demoted to the
    ``BINDS_TO`` weight, since co-precipitation evidences association at most.

    >>> evidence_score("HAS_INHIBITORY_ACTIVITY_AGAINST")
    3
    >>> evidence_score("HAS_INHIBITORY_ACTIVITY_AGAINST", "colloidal coomassie LC-MS/MS")
    2
    >>> evidence_score("SOMETHING_ELSE")
    0
    """
    base = EVIDENCE.weights.get(str(relationship).strip(), 0)
    if (
        demote_proteomics
        and str(relationship).strip() == "HAS_INHIBITORY_ACTIVITY_AGAINST"
        and is_proteomics_assay(assay_description)
    ):
        return EVIDENCE.weights["BINDS_TO"]
    return base


def development_score(max_phase: object, status: object = None) -> int:
    """Score clinical development stage. Ports A15:112-120."""
    try:
        phase = float(max_phase)
    except (TypeError, ValueError):
        phase = float("nan")
    if not np.isnan(phase):
        if phase >= 2:
            return 2
        if phase == 1:
            return 1
    text = "" if status is None else str(status).lower()
    if any(k in text for k in ("phase 2", "phase 3", "approved", "launched")):
        return 2
    if "phase 1" in text:
        return 1
    return 0


# --------------------------------------------------------------------------
# Table-level scoring
# --------------------------------------------------------------------------

#: Column aliases, so the loader tolerates both the R pipeline's output schema
#: and a freshly-fetched ChEMBL frame.
_COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "molecule_chembl_id": ("molecule_chembl_id",),
    "display_name": ("display_name", "pref_name", "molecule_chembl_id"),
    "protein": ("queried_protein", "target_pref_name", "protein"),
    "relationship": ("corrected_relationship", "relationship"),
    "value_nM": (
        "standard_value_numeric",
        "minimum_activity_nM",
        "standard_value",
    ),
    "relation": ("standard_relations", "standard_relation"),
    "assay_description": ("assay_description",),
    "max_phase": ("max_phase",),
    "development_status": ("development_status",),
    "assay_id": ("assay_chembl_id",),
    "document_id": ("final_document_chembl_id", "document_chembl_id"),
}


def _resolve(frame: pd.DataFrame, key: str) -> pd.Series:
    """Return the first aliased column present, or an all-NA series."""
    for name in _COLUMN_ALIASES[key]:
        if name in frame.columns:
            return frame[name]
    return pd.Series([None] * len(frame), index=frame.index, dtype="object")


def score_evidence_table(
    evidence: pd.DataFrame, mode: ScoringMode | None = None
) -> pd.DataFrame:
    """Attach per-record scores and provenance flags to an evidence table.

    Accepts the R pipeline's ``A12B_compound_protein_evidence_corrected.csv``
    schema directly. Returns a new frame; the input is not mutated.
    """
    mode = mode or ScoringMode()
    out = pd.DataFrame(index=evidence.index)
    for key in _COLUMN_ALIASES:
        out[key] = _resolve(evidence, key).values

    out["is_censored"] = out["relation"].map(is_censored)
    out["is_proteomics_assay"] = out["assay_description"].map(is_proteomics_assay)

    out["potency"] = [
        potency_score(v, r, relation_aware=mode.relation_aware)
        for v, r in zip(out["value_nM"], out["relation"])
    ]
    out["evidence"] = [
        evidence_score(rel, desc, demote_proteomics=mode.demote_proteomics)
        for rel, desc in zip(out["relationship"], out["assay_description"])
    ]
    out["development"] = [
        development_score(p, s)
        for p, s in zip(out["max_phase"], out["development_status"])
    ]

    protein = out["protein"].astype(str).str.upper()
    out["is_eef1a"] = protein.isin(EEF1A_PARALOGUES)
    out["is_eef1_complex"] = protein.isin(EEF1_COMPLEX)
    out["is_ribosomal"] = protein.str.match(r"^RP[LS]\d")

    #: Target relevance, ported from the A15 relevance tiers.
    out["relevance"] = np.select(
        [out["is_eef1a"], out["is_eef1_complex"], out["is_ribosomal"]],
        [4, 3, 2],
        default=1,
    )
    out["record_score"] = (
        out["relevance"] + out["evidence"] + out["potency"] + out["development"]
    )
    return out


def rank_candidates(
    evidence: pd.DataFrame, mode: ScoringMode | None = None
) -> pd.DataFrame:
    """Aggregate scored records to one ranked row per compound.

    The ranking key is the compound's best single record score, tie-broken by
    total evidence breadth (distinct proteins) then by summed potency -- the
    same ordering intent as A15's ``network_total_score``.
    """
    mode = mode or ScoringMode()
    scored = score_evidence_table(evidence, mode)

    grouped = scored.groupby("molecule_chembl_id", dropna=False)
    summary = pd.DataFrame(
        {
            "display_name": grouped["display_name"].first(),
            "n_records": grouped.size(),
            "n_proteins": grouped["protein"].nunique(),
            "best_record_score": grouped["record_score"].max(),
            "total_potency": grouped["potency"].sum(),
            "max_evidence": grouped["evidence"].max(),
            "hits_eef1a": grouped["is_eef1a"].any(),
            "n_censored": grouped["is_censored"].sum(),
            "n_proteomics": grouped["is_proteomics_assay"].sum(),
            "best_protein": grouped.apply(
                lambda g: g.loc[g["record_score"].idxmax(), "protein"],
                include_groups=False,
            ),
        }
    ).reset_index()

    #: A candidate whose every record is censored has no positive evidence at
    #: all; surface that rather than letting it rank on relevance alone.
    summary["all_evidence_censored"] = summary["n_censored"] == summary["n_records"]

    summary = summary.sort_values(
        ["best_record_score", "n_proteins", "total_potency"],
        ascending=False,
        kind="mergesort",
    ).reset_index(drop=True)
    summary.insert(0, "rank", np.arange(1, len(summary) + 1))
    return summary
