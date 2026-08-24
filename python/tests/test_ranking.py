"""Tests for relation-aware scoring.

The first test class is the regression test for CODE_REVIEW.md S1-1: it encodes
the exact molibresib records that produced the spurious rank #1, and asserts
that the corrected scoring rejects them while the legacy mode reproduces the
bug. If someone later "simplifies" potency_score back to a value-only ladder,
these fail.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from resko_network.ranking import (
    ScoringMode,
    development_score,
    evidence_score,
    is_censored,
    is_proteomics_assay,
    potency_score,
    rank_candidates,
    score_evidence_table,
)

# The two real censored records from
# results/A12B_compound_protein_evidence_corrected.csv, plus a genuine
# measurement from a different compound for contrast.
MOLIBRESIB_ASSAY = (
    "Inhibition of EEF1A1 in human MV4-11 cells assessed as protein "
    "abundance by colloidal coomassie staining based LC-MS/MS analysis"
)

EVIDENCE_FIXTURE = pd.DataFrame(
    [
        {
            "molecule_chembl_id": "CHEMBL1232461",
            "display_name": "MOLIBRESIB",
            "queried_protein": "EEF1A1",
            "corrected_relationship": "HAS_INHIBITORY_ACTIVITY_AGAINST",
            "minimum_activity_nM": 10000.0,
            "standard_relations": ">=",
            "assay_description": MOLIBRESIB_ASSAY,
            "max_phase": 2.0,
            "development_status": None,
        },
        {
            "molecule_chembl_id": "CHEMBL1232461",
            "display_name": "MOLIBRESIB",
            "queried_protein": "RPS3",
            "corrected_relationship": "HAS_INHIBITORY_ACTIVITY_AGAINST",
            "minimum_activity_nM": 240.0,
            "standard_relations": ">=",
            "assay_description": MOLIBRESIB_ASSAY,
            "max_phase": 2.0,
            "development_status": None,
        },
        {
            "molecule_chembl_id": "CHEMBL1802814",
            "display_name": "CHEMBL1802814",
            "queried_protein": "EEF1A1",
            "corrected_relationship": "BINDS_TO",
            "minimum_activity_nM": 3.74,
            "standard_relations": "=",
            "assay_description": "Binding affinity to human EEF1A1 by SPR",
            "max_phase": np.nan,
            "development_status": None,
        },
    ]
)


class TestCensoredRelations:
    """S1-1 regression: censored records must not earn potency credit."""

    @pytest.mark.parametrize("relation", [">", ">=", ">>"])
    def test_censored_relations_score_zero(self, relation):
        assert potency_score(240.0, relation) == 0

    @pytest.mark.parametrize("relation", ["=", "~", "<", "<="])
    def test_measured_relations_score_normally(self, relation):
        # 240 nM sits in the <=1000 bucket -> 2, matching A15's ladder.
        assert potency_score(240.0, relation) == 2

    def test_legacy_mode_reproduces_the_bug(self):
        """The R original's behaviour, kept switchable for comparison."""
        assert potency_score(240.0, ">=", relation_aware=False) == 2
        assert potency_score(10000.0, ">=", relation_aware=False) == 1

    def test_the_specific_molibresib_records(self):
        """Both real censored records score 0 corrected, >0 under legacy.

        The legacy scores 1 and 2 are exactly the ``direct_potency_score`` and
        ``network_potency_score`` recorded for molibresib in
        ``A15_candidate_score_components.csv``, which is what confirms this
        fixture reproduces the real bug rather than an invented one.
        """
        assert potency_score(10000.0, ">=") == 0
        assert potency_score(240.0, ">=") == 0
        assert potency_score(10000.0, ">=", relation_aware=False) == 1
        assert potency_score(240.0, ">=", relation_aware=False) == 2

    def test_missing_relation_is_not_treated_as_censored(self):
        """Absent relation means unknown, not censored -- keep the ladder."""
        assert potency_score(3.74, None) == 4
        assert potency_score(3.74, np.nan) == 4

    def test_is_censored_predicate(self):
        assert is_censored(">=") is True
        assert is_censored(" > ") is True
        assert is_censored("=") is False
        assert is_censored(None) is False


class TestPotencyLadder:
    """The ladder itself must stay identical to A15:101-110."""

    @pytest.mark.parametrize(
        ("value", "expected"),
        [(1.0, 4), (10.0, 4), (10.001, 3), (100.0, 3), (1000.0, 2),
         (10000.0, 1), (10001.0, 0), (1e9, 0)],
    )
    def test_breakpoints(self, value, expected):
        assert potency_score(value, "=") == expected

    def test_missing_value_scores_zero(self):
        assert potency_score(None, "=") == 0
        assert potency_score(np.nan, "=") == 0
        assert potency_score("not a number", "=") == 0


class TestProteomicsDemotion:
    """S1-1, second half: a pulldown does not evidence inhibition."""

    def test_proteomics_assay_detected(self):
        assert is_proteomics_assay(MOLIBRESIB_ASSAY) is True

    def test_functional_assay_not_flagged(self):
        assert is_proteomics_assay("Inhibition of EEF1A1 GTPase activity") is False

    def test_inhibition_claim_demoted_to_binding(self):
        assert evidence_score("HAS_INHIBITORY_ACTIVITY_AGAINST") == 3
        assert evidence_score("HAS_INHIBITORY_ACTIVITY_AGAINST", MOLIBRESIB_ASSAY) == 2

    def test_demotion_can_be_disabled(self):
        assert (
            evidence_score(
                "HAS_INHIBITORY_ACTIVITY_AGAINST",
                MOLIBRESIB_ASSAY,
                demote_proteomics=False,
            )
            == 3
        )

    def test_unknown_relationship_scores_zero(self):
        assert evidence_score("SOMETHING_ELSE") == 0


class TestDevelopmentScore:
    def test_phase_from_number(self):
        assert development_score(3.0) == 2
        assert development_score(1.0) == 1
        assert development_score(0.0) == 0

    def test_phase_from_status_text(self):
        assert development_score(None, "Approved") == 2
        assert development_score(np.nan, "Phase 1 clinical") == 1
        assert development_score(None, "Preclinical") == 0


class TestRankingIntegration:
    """End-to-end: the correction must actually change the ranking."""

    def test_corrected_ranking_demotes_molibresib(self):
        ranked = rank_candidates(EVIDENCE_FIXTURE, ScoringMode())
        assert ranked.loc[0, "molecule_chembl_id"] == "CHEMBL1802814"
        molibresib = ranked[ranked.molecule_chembl_id == "CHEMBL1232461"].iloc[0]
        assert molibresib["rank"] == 2
        assert bool(molibresib["all_evidence_censored"]) is True
        assert int(molibresib["n_censored"]) == 2

    def test_legacy_ranking_reproduces_molibresib_first(self):
        ranked = rank_candidates(EVIDENCE_FIXTURE, ScoringMode.legacy())
        assert ranked.loc[0, "molecule_chembl_id"] == "CHEMBL1232461"

    def test_censored_flag_is_surfaced_per_record(self):
        scored = score_evidence_table(EVIDENCE_FIXTURE)
        assert scored["is_censored"].tolist() == [True, True, False]
        assert scored["potency"].tolist() == [0, 0, 4]

    def test_target_relevance_tiers(self):
        scored = score_evidence_table(EVIDENCE_FIXTURE)
        assert scored.loc[0, "relevance"] == 4  # EEF1A1
        assert scored.loc[1, "relevance"] == 2  # RPS3, ribosomal

    def test_alias_resolution_accepts_raw_chembl_schema(self):
        """A freshly-fetched frame uses different column names; both must work."""
        raw = pd.DataFrame(
            [
                {
                    "molecule_chembl_id": "CHEMBL1232461",
                    "target_pref_name": "EEF1A1",
                    "relationship": "HAS_INHIBITORY_ACTIVITY_AGAINST",
                    "standard_value": 240.0,
                    "standard_relation": ">=",
                    "assay_description": MOLIBRESIB_ASSAY,
                }
            ]
        )
        scored = score_evidence_table(raw)
        assert scored.loc[0, "is_censored"] is True or scored.loc[0, "is_censored"]
        assert scored.loc[0, "potency"] == 0
