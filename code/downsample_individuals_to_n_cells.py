from __future__ import annotations

import hashlib
from typing import Any, Optional

import numpy as np
import pandas as pd


def _parse_target(value) -> Optional[int]:
    """Return None for unset targets (None or -1)."""
    if value is None:
        return None
    value = int(value)
    if value < 0:
        return None
    return value


def downsample_individuals_to_n_cells(
    adata: Any,
    n_cells_target=None,
    *,
    individual_col: str = "individual",
    seed_key: str = "",
    copy: bool = True,
):
    """
    Keep individuals with at least `n_cells_target` cells, then downsample each
    to exactly `n_cells_target` cells.

    A target of None or -1 leaves `adata` unchanged (or returns a copy if
    copy=True). Sampling is deterministic given `seed_key` and `n_cells_target`.
    """
    n_cells_target = _parse_target(n_cells_target)

    if n_cells_target is None:
        return adata.copy() if copy else adata

    if n_cells_target == 0:
        raise ValueError("n_cells_target must be > 0 when set")

    if individual_col not in adata.obs.columns:
        raise KeyError(f"Missing obs column '{individual_col}'")

    indiv_col = adata.obs[individual_col].to_numpy()
    if indiv_col.size == 0:
        return adata.copy() if copy else adata

    counts = pd.Series(indiv_col).value_counts()
    keep_indivs = set(counts.index[counts >= n_cells_target].tolist())
    if not keep_indivs:
        keep_idx = np.array([], dtype=int)
        view = adata[keep_idx, :]
        return view.copy() if copy else view

    seed_hex = hashlib.md5(
        f"{seed_key}|n_cells_target={n_cells_target}".encode("utf-8")
    ).hexdigest()[:8]
    base_seed = int(seed_hex, 16)

    indivs = sorted(keep_indivs)
    keep_idx_chunks: list[np.ndarray] = []

    for i, indiv in enumerate(indivs):
        indiv_idx = np.flatnonzero(indiv_col == indiv)
        n = indiv_idx.size
        if n < n_cells_target:
            continue

        rng = np.random.default_rng(base_seed + i)
        chosen_idx = rng.choice(indiv_idx, size=n_cells_target, replace=False)
        keep_idx_chunks.append(chosen_idx)

    if keep_idx_chunks:
        keep_idx = np.sort(np.concatenate(keep_idx_chunks))
    else:
        keep_idx = np.array([], dtype=int)

    view = adata[keep_idx, :]
    return view.copy() if copy else view
