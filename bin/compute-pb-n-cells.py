#! /usr/bin/env python3

import sys
from pathlib import Path
import anndata as ad
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from subsample_cells import subsample_cells_within_individuals  # noqa: E402
from subsample_individuals import subsample_individuals  # noqa: E402
from downsample_individuals_to_n_cells import downsample_individuals_to_n_cells  # noqa: E402
from cell_label_subset import cell_label_mask  # noqa: E402

cell_label = sys.argv[1]
cell_frac = float(sys.argv[2])
indiv_frac = float(sys.argv[3])
adata = sc.read_h5ad(sys.argv[4], backed="r")
n_cells_target = int(sys.argv[5])

obs = adata.obs.loc[cell_label_mask(adata.obs["cell_label"], cell_label).to_numpy()]
cells = ad.AnnData(obs=obs[["individual"]].copy())

# Same cell selection as compute-pb-logcounts.py.
cells = downsample_individuals_to_n_cells(
    cells, n_cells_target, individual_col="individual", seed_key=cell_label, copy=True
)
if indiv_frac != 1.0:
    cells = subsample_individuals(
        cells, indiv_frac, seed_key=cell_label, individual_col="individual", copy=True
    )
if cell_frac != 1.0:
    cells = subsample_cells_within_individuals(
        cells, cell_frac, seed_key=cell_label, individual_col="individual", copy=True
    )

n_cells = (
    cells.obs["individual"]
    .astype(str)
    .value_counts(sort=False)
    .rename_axis("sample_id")
    .reset_index(name="n_cells")
)
n_cells = n_cells[n_cells["n_cells"] > 0]
n_cells.to_csv(f"{cell_label}-pb-n-cells.tsv", sep="\t", index=False)
