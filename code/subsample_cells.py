
from __future__ import annotations

import hashlib
from typing import Any

import numpy as np
import pandas as pd


def subsample_cells_within_individuals(
    adata: Any,
    cell_frac: float,
    *,
    individual_col: str = "individual",
    seed_key: str = "",
    copy: bool = True,
):
    """
    Deterministically subsample cells within each individual.

    For each value in `adata.obs[individual_col]`, retain
    floor(cell_frac * n_indiv_cells) cells (at least 1 if cell_frac > 0), using a
    deterministic RNG seed derived from f"{seed_key}|{cell_frac}" and the
    individual's position in sorted unique individuals.

    If `cell_frac >= 1.0`, returns `adata` unchanged (or `adata.copy()` if copy=True).
    """
    cell_frac = float(cell_frac)

    if cell_frac >= 1.0:
        return adata.copy() if copy else adata

    if cell_frac < 0.0:
        raise ValueError(f"cell_frac must be >= 0, got {cell_frac}")

    if individual_col not in adata.obs.columns:
        raise KeyError(f"Missing obs column '{individual_col}'")

    indiv_col = adata.obs[individual_col].to_numpy()
    if indiv_col.size == 0:
        return adata.copy() if copy else adata

    seed_hex = hashlib.md5(f"{seed_key}|{cell_frac}".encode("utf-8")).hexdigest()[:8]
    base_seed = int(seed_hex, 16)

    # Sort for determinism (so the same seed yields the same cells).
    indivs = sorted(pd.unique(indiv_col).tolist())
    keep_idx_chunks: list[np.ndarray] = []

    for i, indiv in enumerate(indivs):
        indiv_idx = np.flatnonzero(indiv_col == indiv)
        n = indiv_idx.size
        if n == 0:
            continue

        k = int(np.floor(cell_frac * n))
        if cell_frac > 0.0:
            k = max(1, k)
        k = min(n, k)

        rng = np.random.default_rng(base_seed + i)
        chosen_idx = rng.choice(indiv_idx, size=k, replace=False)
        keep_idx_chunks.append(chosen_idx)

    if keep_idx_chunks:
        keep_idx = np.sort(np.concatenate(keep_idx_chunks))
    else:
        keep_idx = np.array([], dtype=int)

    view = adata[keep_idx, :]
    return view.copy() if copy else view

