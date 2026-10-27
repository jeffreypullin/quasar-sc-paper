#!/usr/bin/env python3

"""Per-gene sparsity under donor min-cells filters for OneK1K cell types."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from cell_label_subset import cell_label_mask  # noqa: E402
from filter_individuals_min_n_cells import filter_individuals_min_n_cells  # noqa: E402

CELL_TYPES = [
    "CD4_NC",
    "CD4_ET",
    "CD4_SOX4",
    "CD8_ET",
    "CD8_NC",
    "CD8_S100B",
    "NK",
    "NK_R",
    "Plasma",
    "B_Mem",
    "B_IN",
    "Mono_C",
    "Mono_NC",
    "DC",
]

MIN_CELLS = [0, 5, 10, 20]


def nonzero_frac(X, n_obs: int) -> np.ndarray:
    if n_obs == 0:
        return np.zeros(X.shape[1], dtype=float)
    if sp.issparse(X):
        return np.asarray(X.tocsr().getnnz(axis=0)).ravel() / float(n_obs)
    return (np.asarray(X) != 0).mean(axis=0)


def main() -> None:
    adata_path = sys.argv[1]
    anno_path = sys.argv[2]

    adata = sc.read_h5ad(adata_path, backed="r")
    anno_ids = set(pd.read_csv(anno_path, sep="\t")["phenotype_id"])

    stats_rows: list[pd.DataFrame] = []
    donor_rows: list[dict] = []

    for cell_type in CELL_TYPES:
        subset = adata[
            cell_label_mask(adata.obs["cell_label"], cell_type).to_numpy()
        ].to_memory()

        for min_cells in MIN_CELLS:
            if min_cells == 0:
                filtered = subset
            else:
                filtered = filter_individuals_min_n_cells(
                    subset,
                    min_cells,
                    individual_col="individual",
                    copy=True,
                )

            n_donors = (
                int(filtered.obs["individual"].nunique()) if filtered.n_obs > 0 else 0
            )
            donor_rows.append(
                {
                    "cell_type": cell_type,
                    "min_cells": min_cells,
                    "n_donors": n_donors,
                    "n_cells": int(filtered.n_obs),
                }
            )

            nz = nonzero_frac(filtered.X, filtered.n_obs)
            feature_id = filtered.var_names.to_numpy()
            keep = np.isin(feature_id, list(anno_ids))
            stats_rows.append(
                pd.DataFrame(
                    {
                        "cell_type": cell_type,
                        "min_cells": min_cells,
                        "feature_id": feature_id[keep],
                        "sc_non_zero_frac": nz[keep],
                    }
                )
            )

    pd.concat(stats_rows, ignore_index=True).to_csv(
        "gene-filter-stats.tsv", sep="\t", index=False
    )
    pd.DataFrame(donor_rows).to_csv(
        "gene-filter-donor-summary.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
