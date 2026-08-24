"""Tests for scaffold-aware similarity screening.

The core assertion is CODE_REVIEW.md S1-2: molibresib, alprazolam and triazolam
share an identical Bemis-Murcko scaffold, so their apparent similarity is a
structural echo. The screen must flag that automatically.
"""

from __future__ import annotations

import pandas as pd
import pytest
from rdkit import Chem

from resko_network.config import BANDS
from resko_network.similarity import (
    generic_scaffold_smiles,
    make_generator,
    murcko_smiles,
    prepare,
    screen,
)

SMILES = {
    "MOLIBRESIB": "CCNC(=O)C[C@@H]1N=C(c2ccc(Cl)cc2)c2cc(OC)ccc2-n2c(C)nnc21",
    "alprazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",
    "triazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2",
    "temazepam": "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21",
    "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
}


@pytest.fixture(scope="module")
def generator():
    return make_generator()


@pytest.fixture(scope="module")
def molecules(generator):
    frame = pd.DataFrame(
        {"name": list(SMILES), "smiles": [SMILES[k] for k in SMILES]}
    )
    return prepare(frame, "name", "smiles", generator)


class TestScaffoldIdentity:
    """The finding itself, asserted structurally rather than by eye."""

    def test_molibresib_shares_scaffold_with_alprazolam(self):
        a = Chem.MolFromSmiles(SMILES["MOLIBRESIB"])
        b = Chem.MolFromSmiles(SMILES["alprazolam"])
        assert murcko_smiles(a) == murcko_smiles(b)

    def test_molibresib_shares_scaffold_with_triazolam(self):
        a = Chem.MolFromSmiles(SMILES["MOLIBRESIB"])
        b = Chem.MolFromSmiles(SMILES["triazolam"])
        assert murcko_smiles(a) == murcko_smiles(b)

    def test_temazepam_scaffold_differs(self):
        """Temazepam is a benzodiazepine but lacks the fused triazole."""
        a = Chem.MolFromSmiles(SMILES["MOLIBRESIB"])
        b = Chem.MolFromSmiles(SMILES["temazepam"])
        assert murcko_smiles(a) != murcko_smiles(b)

    def test_kinase_inhibitor_scaffold_is_unrelated(self):
        a = Chem.MolFromSmiles(SMILES["MOLIBRESIB"])
        b = Chem.MolFromSmiles(SMILES["imatinib"])
        assert murcko_smiles(a) != murcko_smiles(b)
        assert generic_scaffold_smiles(a) != generic_scaffold_smiles(b)


class TestScreenFlagging:
    def test_artefact_flag_set_for_identical_scaffold(self, molecules):
        queries = [m for m in molecules if m.identifier == "MOLIBRESIB"]
        library = [m for m in molecules if m.identifier != "MOLIBRESIB"]
        result = screen(queries, library)
        flagged = set(result.loc[result.scaffold_artefact, "library_id"])
        assert flagged == {"alprazolam", "triazolam"}

    def test_scaffold_tanimoto_is_unity_for_flagged_pairs(self, molecules):
        queries = [m for m in molecules if m.identifier == "MOLIBRESIB"]
        library = [m for m in molecules if m.identifier in ("alprazolam", "triazolam")]
        result = screen(queries, library)
        assert (result["scaffold_tanimoto"] == 1.0).all()

    def test_full_molecule_similarity_stays_below_close_band(self, molecules):
        """The reported hits never reach the script's own 0.70 band."""
        queries = [m for m in molecules if m.identifier == "MOLIBRESIB"]
        library = [m for m in molecules if m.identifier != "MOLIBRESIB"]
        result = screen(queries, library)
        assert result["tanimoto"].max() < BANDS.close
        assert (result["similarity_band"] == "weak_or_partial_similarity").all()

    def test_threshold_is_respected(self, molecules):
        result = screen(molecules, molecules, threshold=0.99)
        # Only self-pairs survive a 0.99 threshold.
        assert (result.query_id == result.library_id).all()

    def test_empty_result_is_a_frame(self, molecules):
        result = screen(molecules, [], threshold=0.3)
        assert isinstance(result, pd.DataFrame)
        assert result.empty


class TestPrepare:
    def test_unparseable_smiles_are_skipped(self, generator):
        frame = pd.DataFrame(
            {"name": ["good", "bad", "empty"],
             "smiles": [SMILES["alprazolam"], "not_a_smiles((", ""]}
        )
        assert [m.identifier for m in prepare(frame, "name", "smiles", generator)] == [
            "good"
        ]

    def test_fingerprint_parameters_match_config(self, molecules):
        from resko_network.config import FINGERPRINT

        assert molecules[0].fingerprint.GetNumBits() == FINGERPRINT.n_bits
