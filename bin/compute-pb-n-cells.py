#! /usr/bin/env python3

import sys
from pathlib import Path
import anndata as ad
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from apply_cell_selection import apply_cell_selection  # noqa: E402
from cell_label_subset import cell_label_mask  # noqa: E402

cell_label = sys.argv[1]
cell_frac = float(sys.argv[2])
indiv_frac = float(sys.argv[3])
adata = sc.read_h5ad(sys.argv[4], backed="r")
min_n_cells = int(sys.argv[5])

obs = adata.obs.loc[cell_label_mask(adata.obs["cell_label"], cell_label).to_numpy()]
cells = ad.AnnData(obs=obs[["individual"]].copy())

# Same cell selection as compute-pb-logcounts.py.
cells = apply_cell_selection(
    cells,
    cell_frac=cell_frac,
    indiv_frac=indiv_frac,
    min_n_cells=min_n_cells,
    seed_key=cell_label,
    individual_col="individual",
    copy=True,
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
