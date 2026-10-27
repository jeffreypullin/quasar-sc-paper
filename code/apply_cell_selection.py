"""Apply indiv/cell subsampling, then the min_n_cells donor filter.

Order matters: cell_frac first, then min_n_cells on the *post*-subsample
counts so donors that fall below the threshold after thinning are dropped.
"""

from __future__ import annotations

from typing import Any

from filter_individuals_min_n_cells import filter_individuals_min_n_cells
from subsample_cells import subsample_cells_within_individuals
from subsample_individuals import subsample_individuals


def apply_cell_selection(
    adata: Any,
    *,
    cell_frac: float = 1.0,
    indiv_frac: float = 1.0,
    min_n_cells=None,
    seed_key: str = "",
    individual_col: str = "individual",
    copy: bool = True,
):
    if indiv_frac != 1.0:
        adata = subsample_individuals(
            adata,
            indiv_frac,
            seed_key=seed_key,
            individual_col=individual_col,
            copy=copy,
        )
        copy = True

    if cell_frac != 1.0:
        adata = subsample_cells_within_individuals(
            adata,
            cell_frac,
            seed_key=seed_key,
            individual_col=individual_col,
            copy=copy,
        )
        copy = True

    return filter_individuals_min_n_cells(
        adata,
        min_n_cells,
        individual_col=individual_col,
        copy=copy,
    )
