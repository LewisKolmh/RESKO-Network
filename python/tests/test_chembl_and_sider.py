"""Tests for cache-key canonicalisation (S2-1) and seed matching (S3-1)."""

from __future__ import annotations

import json

import pandas as pd
import pytest

from resko_network.chembl import ACTIVITY_FIELDS, activities_to_frame, audit_cache_dir, cache_key
from resko_network.sider import (
    coverage_report,
    is_category_label,
    match_drug_names,
    normalise_name,
)


class TestCacheKey:
    """S2-1: exactly one filename convention, enforced by test."""

    @pytest.mark.parametrize(
        ("offset", "expected"),
        [(0, "activities_CHEMBL941_00000.json"),
         (1000, "activities_CHEMBL941_01000.json"),
         (2000, "activities_CHEMBL941_02000.json"),
         (99000, "activities_CHEMBL941_99000.json")],
    )
    def test_offsets_are_zero_padded_to_five(self, offset, expected):
        assert cache_key("CHEMBL941", offset) == expected

    def test_unpadded_form_is_never_produced(self):
        """The orphaned-file bug was an unpadded name; assert it can't recur."""
        assert cache_key("CHEMBL941", 0) != "activities_CHEMBL941_0.json"

    def test_kind_prefix_is_configurable(self):
        assert cache_key("CHEMBL25", 0, kind="molecule").startswith("molecule_")

    def test_audit_detects_non_canonical_names(self, tmp_path):
        payload = {"activities": [{"activity_id": 1, "standard_relation": ">="}]}
        (tmp_path / "activities_CHEMBL941_00000.json").write_text(json.dumps(payload))
        (tmp_path / "activities_CHEMBL941_0.json").write_text(json.dumps(payload))
        audit = audit_cache_dir(tmp_path).set_index("file")
        assert bool(audit.loc["activities_CHEMBL941_00000.json", "canonical_name"])
        assert not bool(audit.loc["activities_CHEMBL941_0.json", "canonical_name"])

    def test_non_paged_files_are_not_judged_on_padding(self, tmp_path):
        """``molecule_alprazolam.json`` has no offset; padding is meaningless."""
        (tmp_path / "molecule_alprazolam.json").write_text(json.dumps({"molecules": []}))
        audit = audit_cache_dir(tmp_path).set_index("file")
        assert not bool(audit.loc["molecule_alprazolam.json", "paged"])
        assert pd.isna(audit.loc["molecule_alprazolam.json", "canonical_name"])

    def test_audit_surfaces_field_count_divergence(self, tmp_path):
        """The real failure mode: same records, different serialised key sets."""
        rich = {"activities": [{"activity_id": 1, "standard_relation": ">=",
                               "standard_upper_value": 10000, "action_type": None}]}
        poor = {"activities": [{"activity_id": 1, "standard_relation": ">="}]}
        (tmp_path / "activities_X_00000.json").write_text(json.dumps(poor))
        (tmp_path / "activities_X_0.json").write_text(json.dumps(rich))
        audit = audit_cache_dir(tmp_path).set_index("file")
        assert audit.loc["activities_X_0.json", "n_fields_union"] == 4
        assert audit.loc["activities_X_00000.json", "n_fields_union"] == 2

    def test_audit_detects_within_file_key_variation(self, tmp_path):
        """union > common means the serialiser dropped null fields per record."""
        payload = {"activities": [
            {"activity_id": 1, "standard_relation": ">=", "standard_upper_value": 1},
            {"activity_id": 2, "standard_relation": "="},
        ]}
        (tmp_path / "activities_Y_00000.json").write_text(json.dumps(payload))
        row = audit_cache_dir(tmp_path).set_index("file").loc["activities_Y_00000.json"]
        assert row["n_fields_union"] == 3
        assert row["n_fields_common"] == 2


class TestActivityFrame:
    def test_standard_relation_is_always_present(self):
        """The field the R ranking dropped must survive flattening."""
        assert "standard_relation" in ACTIVITY_FIELDS
        assert "standard_upper_value" in ACTIVITY_FIELDS

    def test_missing_fields_are_filled_not_dropped(self):
        frame = activities_to_frame([{"activity_id": 1}])
        assert list(frame.columns) == list(ACTIVITY_FIELDS)
        assert frame.loc[0, "standard_relation"] is None


DRUG_NAMES = pd.DataFrame(
    {
        "stitch_flat_id": [f"CID{i:09d}" for i in range(6)],
        "drug_name": [
            "gamma-aminobutyric acid",
            "sugammadex",
            "alprazolam",
            "Didemnin B",
            "plitidepsin",
            "NSC",
        ],
    }
)


class TestSeedMatching:
    """S3-1: unanchored substring matching must not recur."""

    def test_short_probe_does_not_match_substrings(self):
        """'gam' matched 2 unrelated drugs under the R str_detect approach."""
        matched = match_drug_names(["gam"], DRUG_NAMES)
        assert not matched["matched"].any()

    def test_nsc_probe_matches_only_whole_token(self):
        matched = match_drug_names(["NSC"], DRUG_NAMES)
        assert matched.loc[matched.matched, "sider_drug_name"].tolist() == ["NSC"]

    def test_real_seed_still_matches(self):
        matched = match_drug_names(["Didemnin B"], DRUG_NAMES)
        assert matched.loc[matched.matched, "sider_drug_name"].tolist() == ["Didemnin B"]

    def test_hyphen_and_case_normalisation(self):
        assert normalise_name("Ternatin-4") == "ternatin 4"
        assert normalise_name("  Didemnin  B ") == "didemnin b"

    def test_unmatched_seed_yields_a_row(self):
        """Absence must be recorded, not silently dropped."""
        matched = match_drug_names(["ternatin-4"], DRUG_NAMES)
        assert len(matched) == 1
        assert not matched.loc[0, "matched"]


class TestCoverageDenominator:
    """S3-1, second half: category labels inflate the denominator."""

    def test_category_labels_detected(self):
        assert is_category_label("Didemnin analogues") is True
        assert is_category_label("Synthetic Ternatin analogues") is True
        assert is_category_label("didemnin B") is False

    def test_denominator_excludes_category_labels(self):
        seeds = ["plitidepsin", "ternatin-4", "Didemnin analogues"]
        report = coverage_report(match_drug_names(seeds, DRUG_NAMES))
        assert report["n_seeds_supplied"] == 3
        assert report["n_category_labels_excluded"] == 1
        assert report["n_compounds_searched"] == 2
        assert report["n_matched"] == 1  # plitidepsin only

    def test_zero_coverage_is_reported_as_zero_not_error(self):
        report = coverage_report(match_drug_names(["ternatin-4"], DRUG_NAMES))
        assert report["coverage_fraction"] == 0.0
        assert report["matched_names"] == []
