"""Central configuration.

Every threshold, weight and magic number used anywhere in the pipeline lives
here. In the R original these were scattered across scripts (e.g. the 0.30
candidate threshold at A17C:1142, the potency ladder at A15:101-110), which made
it hard to see what the pipeline's parameters actually were. Changing a
threshold should mean editing exactly one line of exactly one file.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# --------------------------------------------------------------------------
# Seed definitions
# --------------------------------------------------------------------------

#: eEF1A paralogues -- the biological anchor of the whole project.
EEF1A_PARALOGUES: tuple[str, ...] = ("EEF1A1", "EEF1A2")

#: The eEF1 complex beyond the alpha subunits.
EEF1_COMPLEX: tuple[str, ...] = ("EEF1B2", "EEF1D", "EEF1G")

#: Validated eEF1A-binding chemotypes from the literature. These are the
#: compounds with actual eEF1A pharmacology -- the depsipeptide/macrolactam
#: families. Contrast with the ChEMBL bioactivity hits, which are not.
#:
#: Note: the R original's ``eef1a_master_seed_list.csv`` also contained three
#: category labels ("Synthetic Ternatin analogues", "Didemnin analogues",
#: "Plitidepsin analogues"). Those can never match a database row and inflated
#: the reported 0/15 coverage denominator, so they are excluded here. See
#: CODE_REVIEW.md S3-1.
VALIDATED_EEF1A_BINDERS: tuple[str, ...] = (
    "plitidepsin",
    "didemnin B",
    "ternatin-4",
    "ternatin",
    "nannocystin A",
    "nannocystin Ax",
    "cytotrienin A",
    "ansatrienin A",
    "ansatrienin B",
    "metarrestin",
    "gamendazole",
    "cordyheptapeptide A",
)


# --------------------------------------------------------------------------
# Chemical similarity (ports A16 / A17C)
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class FingerprintConfig:
    """Morgan fingerprint parameters.

    Defaults match the R original's reticulate/RDKit call so that ported
    results are numerically comparable. ``include_chirality=True`` matters:
    flipping it changes Tanimoto values in the third decimal.
    """

    radius: int = 2
    n_bits: int = 2048
    include_chirality: bool = True


@dataclass(frozen=True)
class SimilarityBands:
    """Tanimoto banding, ported from ``assign_similarity_band`` (A17C:194-199).

    ``candidate`` is the threshold above which a SIDER drug is reported as a
    candidate (A17C:1142). Note the gap between ``candidate`` (0.30) and
    ``close`` (0.70): in the observed data nothing exceeded 0.505, so every
    reported "candidate" sat in the weakest band. See CODE_REVIEW.md S1-2.
    """

    close: float = 0.70
    moderate: float = 0.50
    candidate: float = 0.30

    def band(self, x: float | None) -> str:
        if x is None:
            return "structurally_distant_under_this_fingerprint"
        if x >= self.close:
            return "close_structural_neighbour"
        if x >= self.moderate:
            return "moderate_similarity"
        if x >= self.candidate:
            return "weak_or_partial_similarity"
        return "structurally_distant_under_this_fingerprint"


#: Scaffold-artefact guard. When a query/hit pair shares an identical
#: Bemis-Murcko scaffold, similarity is explained by the shared core rather
#: than by anything target-specific, and the pair is flagged. This has no
#: counterpart in the R original -- it is the mechanism that would have caught
#: the molibresib/alprazolam/triazolam artefact automatically.
FLAG_IDENTICAL_MURCKO_SCAFFOLD: bool = True


# --------------------------------------------------------------------------
# Ranking (ports A15)
# --------------------------------------------------------------------------

#: ChEMBL ``standard_relation`` values denoting a *censored* measurement:
#: activity was not observed up to the stated concentration. These are evidence
#: AGAINST binding and must never contribute positive potency score.
#: The R original ignored this column entirely -- see CODE_REVIEW.md S1-1.
CENSORED_RELATIONS: frozenset[str] = frozenset({">", ">=", ">>"})

#: Relations denoting a real, bounded measurement.
MEASURED_RELATIONS: frozenset[str] = frozenset({"=", "~", "<", "<="})


@dataclass(frozen=True)
class PotencyLadder:
    """Potency scoring thresholds in nM, ported from A15:101-110."""

    breakpoints: tuple[tuple[float, int], ...] = (
        (10.0, 4),
        (100.0, 3),
        (1000.0, 2),
        (10000.0, 1),
    )


@dataclass(frozen=True)
class EvidenceWeights:
    """Relationship-class weights, ported from ``evidence_score`` (A15:92-99)."""

    weights: dict[str, int] = field(
        default_factory=lambda: {
            "HAS_INHIBITORY_ACTIVITY_AGAINST": 3,
            "BINDS_TO": 2,
            "HAS_ACTIVITY_AGAINST": 1,
        }
    )


#: Assay descriptions matching these patterns should not support an inhibition
#: claim. A protein identified in a pulldown/proteomics experiment is a
#: co-precipitant, not a demonstrated target. This is what made molibresib (a
#: BET inhibitor, via an LC-MS/MS pulldown in nature10509) look like an eEF1A1
#: inhibitor.
PROTEOMICS_ASSAY_PATTERNS: tuple[str, ...] = (
    "lc-ms/ms",
    "coomassie",
    "pulldown",
    "pull-down",
    "affinity capture",
    "chemoproteomic",
)


# --------------------------------------------------------------------------
# ChEMBL API
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ChemblConfig:
    base_url: str = "https://www.ebi.ac.uk/chembl/api/data"
    page_limit: int = 1000
    #: Zero-padding width for the offset component of a cache filename.
    #: The R original used two conventions across revisions ("%d" and "%05d"),
    #: leaving orphaned cache files with a different response schema. Exactly
    #: one convention is defined here. See CODE_REVIEW.md S2-1.
    cache_offset_width: int = 5
    max_attempts: int = 5
    timeout_s: int = 120
    request_pause_s: float = 0.25
    user_agent: str = "RESKO-Network-py/1.0"


# --------------------------------------------------------------------------
# Network construction (ports A4-A7)
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class StringConfig:
    """STRING interaction-network retrieval parameters."""

    species: int = 9606
    #: STRING combined_score floor (0-1000). 700 = "high confidence".
    required_score: int = 700
    base_url: str = "https://string-db.org/api"


FINGERPRINT = FingerprintConfig()
BANDS = SimilarityBands()
POTENCY = PotencyLadder()
EVIDENCE = EvidenceWeights()
CHEMBL = ChemblConfig()
STRING = StringConfig()
