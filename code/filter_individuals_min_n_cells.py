"""Keep donors with at least `min_n_cells` cells; drop the rest."""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
import pandas as pd


def filter_individuals_min_n_cells(
    adata: Any,
    min_n_cells: Optional[int] = None,
    *,
    individual_col: str = "individual",
    copy: bool = True,
):
    """
    Keep individuals with at least `min_n_cells` cells in the current `adata`.

    Intended to run *after* cell/individual fraction subsampling so the
    threshold is on post-subsample counts. None or a negative value leaves
    `adata` unchanged (or returns a copy if copy=True).
    """
    if min_n_cells is None:
        return adata.copy() if copy else adata

    min_n_cells = int(min_n_cells)
    if min_n_cells < 0:
        return adata.copy() if copy else adata

    if min_n_cells == 0:
        raise ValueError("min_n_cells must be > 0 when set")

    if individual_col not in adata.obs.columns:
        raise KeyError(f"Missing obs column '{individual_col}'")

    indiv_col = adata.obs[individual_col].to_numpy()
    if indiv_col.size == 0:
        return adata.copy() if copy else adata

    counts = pd.Series(indiv_col).value_counts()
    keep_indivs = set(counts.index[counts >= min_n_cells].tolist())
    if not keep_indivs:
        keep_idx = np.array([], dtype=int)
    else:
        keep_idx = np.flatnonzero(np.isin(indiv_col, list(keep_indivs)))

    view = adata[keep_idx, :]
    return view.copy() if copy else view
