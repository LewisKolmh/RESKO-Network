"""resko_network -- eEF1A allosteric-inhibitor screening pipeline.

A Python port of the R pipeline in ``RESKO-Network``, restructured so that the
stage boundaries are importable modules rather than numbered scripts. The
mapping from each R script to its Python home is in ``docs/R_TO_PYTHON_MAP.md``.

Three behavioural differences from the R original, each documented in
``CODE_REVIEW.md``:

* :mod:`resko_network.ranking` reads ChEMBL's ``standard_relation``, so censored
  measurements (``IC50 >= 10 uM``) no longer score as weak positives (S1-1).
* :mod:`resko_network.similarity` reports scaffold-level similarity alongside
  full-molecule similarity, so scaffold artefacts are visible (S1-2).
* :mod:`resko_network.chembl` has exactly one cache-key function, so a cache
  directory cannot accumulate same-offset files under two conventions (S2-1).
"""

from __future__ import annotations

__version__ = "1.0.0"

from . import chembl, config, ranking, sider, similarity

__all__ = ["chembl", "config", "ranking", "sider", "similarity", "__version__"]
