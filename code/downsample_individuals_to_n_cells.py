from __future__ import annotations

import hashlib
from typing import Any, Optional

import numpy as np
import pandas as pd


# Targets at or above this value keep donors with at least
# (target - offset) cells and do not downsample those donors' cells.
# Smaller positive targets also downsample each donor to that exact count.
DONOR_ONLY_OFFSET = 1_000_000


def _parse_target(value) -> tuple[Optional[int], bool]:
    """Return (min_cells, equalize).

    None or a negative value leaves the data unchanged.
    A positive value below DONOR_ONLY_OFFSET keeps donors with at least that
    many cells and downsamples each of them to exactly that count.
    A value of DONOR_ONLY_OFFSET + n keeps donors with at least n cells and
    keeps every cell from those donors.
    """
    if value is None:
        return None, False
    value = int(value)
    if value < 0:
        return None, False
    if value >= DONOR_ONLY_OFFSET:
        min_cells = value - DONOR_ONLY_OFFSET
        if min_cells <= 0:
            raise ValueError(
                "donor-only n_cells_target must be DONOR_ONLY_OFFSET + n with n > 0"
            )
        return min_cells, False
    return value, True


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
    copy=True). A target of DONOR_ONLY_OFFSET + n keeps the same individuals
    but all of their cells. Sampling is deterministic given `seed_key` and
    `n_cells_target`.
    """
    n_cells_target, equalize = _parse_target(n_cells_target)

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

        if equalize:
            rng = np.random.default_rng(base_seed + i)
            chosen_idx = rng.choice(indiv_idx, size=n_cells_target, replace=False)
        else:
            chosen_idx = indiv_idx
        keep_idx_chunks.append(chosen_idx)

    if keep_idx_chunks:
        keep_idx = np.sort(np.concatenate(keep_idx_chunks))
    else:
        keep_idx = np.array([], dtype=int)

    view = adata[keep_idx, :]
    return view.copy() if copy else view
