from __future__ import annotations

import hashlib
from typing import Any

import numpy as np
import pandas as pd


def subsample_individuals(
    adata: Any,
    indiv_frac: float,
    *,
    individual_col: str = "individual",
    seed_key: str = "",
    copy: bool = True,
):
    """
    Deterministically subsample individuals (keeping all their cells).

    Retain floor(indiv_frac * n_individuals) individuals (at least 1 if indiv_frac > 0),
    using a deterministic RNG seed derived from f"{seed_key}|{indiv_frac}".

    If `indiv_frac >= 1.0`, returns `adata` unchanged (or `adata.copy()` if copy=True).
    """
    indiv_frac = float(indiv_frac)

    if indiv_frac >= 1.0:
        return adata.copy() if copy else adata

    if indiv_frac < 0.0:
        raise ValueError(f"indiv_frac must be >= 0, got {indiv_frac}")

    if individual_col not in adata.obs.columns:
        raise KeyError(f"Missing obs column '{individual_col}'")

    indiv_col = adata.obs[individual_col].to_numpy()
    if indiv_col.size == 0:
        return adata.copy() if copy else adata

    # Sort for determinism (so the same seed yields the same individuals).
    indivs = sorted(pd.unique(indiv_col).tolist())
    n_indivs = len(indivs)
    if n_indivs == 0:
        return adata.copy() if copy else adata

    k = int(np.floor(indiv_frac * n_indivs))
    if indiv_frac > 0.0:
        k = max(1, k)
    k = min(n_indivs, k)

    seed_hex = hashlib.md5(f"{seed_key}|{indiv_frac}".encode("utf-8")).hexdigest()[:8]
    base_seed = int(seed_hex, 16)

    rng = np.random.default_rng(base_seed)
    chosen = rng.choice(np.asarray(indivs, dtype=object), size=k, replace=False).tolist()
    chosen_set = set(chosen)

    keep_idx = np.flatnonzero(np.isin(indiv_col, list(chosen_set)))
    view = adata[keep_idx, :]
    return view.copy() if copy else view

