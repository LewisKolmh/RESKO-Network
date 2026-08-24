"""Chemical similarity with scaffold-artefact detection -- ports A16 / A17C.

The R original computed Morgan/Tanimoto similarity between candidate compounds
and the SIDER drug library, reporting anything above 0.30 as a candidate. It had
no way to notice that a high-scoring pair shared an identical scaffold, and so
reported six "candidates" that were three benzodiazepines (matching molibresib's
triazolobenzodiazepine core) and three kinase inhibitors (matching a shared
benzamide/pyrimidine-amine motif). See CODE_REVIEW.md S1-2.

:func:`screen` returns the same similarity numbers plus the scaffold columns
needed to see that immediately.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import DataStructs
from rdkit.Chem import rdFingerprintGenerator as rfg
from rdkit.Chem.Scaffolds import MurckoScaffold

from .config import BANDS, FINGERPRINT

RDLogger.DisableLog("rdApp.warning")

__all__ = [
    "Molecule",
    "make_generator",
    "murcko_smiles",
    "generic_scaffold_smiles",
    "prepare",
    "screen",
]


@dataclass
class Molecule:
    """A parsed molecule with its scaffolds and fingerprint precomputed."""

    identifier: str
    smiles: str
    mol: Chem.Mol
    murcko: str
    generic: str
    fingerprint: DataStructs.ExplicitBitVect
    scaffold_fingerprint: DataStructs.ExplicitBitVect


def make_generator(cfg=FINGERPRINT):
    """Build a Morgan fingerprint generator from config.

    Kept as a function rather than a module-level constant so a caller can vary
    parameters without mutating shared state.
    """
    return rfg.GetMorganGenerator(
        radius=cfg.radius,
        fpSize=cfg.n_bits,
        includeChirality=cfg.include_chirality,
    )


def murcko_smiles(mol: Chem.Mol) -> str:
    """Canonical SMILES of the Bemis-Murcko scaffold."""
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))


def generic_scaffold_smiles(mol: Chem.Mol) -> str:
    """Murcko scaffold flattened to remove element and bond-order identity.

    Two molecules sharing a generic scaffold have the same ring topology even if
    heteroatoms differ -- a looser grouping than identical Murcko scaffolds.
    """
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(scaffold))


def prepare(
    frame: pd.DataFrame,
    id_column: str,
    smiles_column: str,
    generator=None,
) -> list[Molecule]:
    """Parse a frame of SMILES into :class:`Molecule` objects.

    Rows whose SMILES fail to parse are skipped silently -- call
    ``len(prepare(...))`` against ``len(frame)`` if you need to detect that.
    """
    generator = generator or make_generator()
    out: list[Molecule] = []
    for _, row in frame.iterrows():
        smiles = row[smiles_column]
        if not isinstance(smiles, str) or not smiles.strip():
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        out.append(
            Molecule(
                identifier=str(row[id_column]),
                smiles=smiles,
                mol=mol,
                murcko=Chem.MolToSmiles(scaffold),
                generic=generic_scaffold_smiles(mol),
                fingerprint=generator.GetFingerprint(mol),
                scaffold_fingerprint=generator.GetFingerprint(scaffold),
            )
        )
    return out


def screen(
    queries: list[Molecule],
    library: list[Molecule],
    *,
    threshold: float | None = None,
) -> pd.DataFrame:
    """Score every query x library pair, flagging scaffold-driven similarity.

    Returns one row per pair at or above ``threshold`` (default: the configured
    candidate threshold), with these columns beyond the raw similarity:

    ``scaffold_tanimoto``
        Similarity between the two Murcko scaffolds. When this approaches 1.0
        the full-molecule similarity is explained by the shared core.
    ``same_murcko`` / ``same_generic_scaffold``
        Exact scaffold identity at the two granularities.
    ``scaffold_artefact``
        True when the pair shares an identical Murcko scaffold -- i.e. the hit
        is structural echo, not evidence of shared target pharmacology.
    """
    threshold = BANDS.candidate if threshold is None else threshold
    rows: list[dict[str, object]] = []
    for query in queries:
        for member in library:
            similarity = DataStructs.TanimotoSimilarity(
                query.fingerprint, member.fingerprint
            )
            if similarity < threshold:
                continue
            scaffold_similarity = DataStructs.TanimotoSimilarity(
                query.scaffold_fingerprint, member.scaffold_fingerprint
            )
            same_murcko = query.murcko == member.murcko
            rows.append(
                {
                    "query_id": query.identifier,
                    "library_id": member.identifier,
                    "tanimoto": round(similarity, 4),
                    "scaffold_tanimoto": round(scaffold_similarity, 4),
                    "similarity_band": BANDS.band(similarity),
                    "same_murcko": same_murcko,
                    "same_generic_scaffold": query.generic == member.generic,
                    "scaffold_artefact": same_murcko,
                    "query_murcko": query.murcko,
                    "library_murcko": member.murcko,
                }
            )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return result.sort_values("tanimoto", ascending=False).reset_index(drop=True)
