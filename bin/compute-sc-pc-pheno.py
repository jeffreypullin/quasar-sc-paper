#! /usr/bin/env python3

import sys
from pathlib import Path

import scanpy as sc
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from subsample_cells import subsample_cells_within_individuals  # noqa: E402
from subsample_individuals import subsample_individuals  # noqa: E402
from cell_label_subset import cell_label_mask  # noqa: E402

cell_type = sys.argv[1]
cell_frac = float(sys.argv[2])
indiv_frac = float(sys.argv[3])
adata = sc.read_h5ad(sys.argv[4], backed="r")

subset = adata[
    cell_label_mask(adata.obs["cell_label"], cell_type).to_numpy()
].to_memory()

if indiv_frac != 1.0:
    subset = subsample_individuals(
        subset,
        indiv_frac,
        seed_key=cell_type,
        individual_col="individual",
        copy=True,
    )

if cell_frac != 1.0:
    subset = subsample_cells_within_individuals(
        subset,
        cell_frac,
        seed_key=cell_type,
        individual_col="individual",
        copy=True,
    )

sc.pp.filter_genes(subset, min_cells=3)
sc.pp.normalize_total(subset)
sc.pp.log1p(subset)

sc.tl.pca(subset)

pc_df = pd.DataFrame(
    subset.obsm["X_pca"][:, :10],
    index=subset.obs_names,
    columns=[f"scPC_{i}" for i in range(1, 11)],
)
pc_df["individual"] = subset.obs["individual"].values

donor_pc = (
    pc_df.groupby("individual", observed=True)
    .mean(numeric_only=True)
    .T
)
donor_pc.index.name = "phenotype_id"

donor_pc.to_csv(f"{cell_type}-sc-pc-pheno.tsv", sep="\t")
