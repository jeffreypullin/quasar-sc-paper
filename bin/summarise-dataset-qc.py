#!/usr/bin/env python3

from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

dataset = sys.argv[1]
adata = sc.read_h5ad(sys.argv[2])

X = adata.X
n_counts = np.asarray(X.sum(axis=1)).ravel()
if sparse.issparse(X):
    nnz = np.diff(X.tocsr().indptr)
else:
    nnz = np.count_nonzero(np.asarray(X), axis=1)
zero_frac = 1 - nnz / adata.n_vars
cell_label = adata.obs["cell_label"].astype(str).to_numpy()

pd.DataFrame({
    "dataset": dataset,
    "cell_label": cell_label,
    "n_counts": n_counts,
}).to_csv("cell-counts.tsv", sep="\t", index=False)
pd.DataFrame({
    "dataset": dataset,
    "cell_label": cell_label,
    "zero_frac": zero_frac,
}).to_csv("zero-frac.tsv", sep="\t", index=False)

cells_per_indiv = (
    adata.obs[["individual"]]
    .assign(individual=lambda df: df["individual"].astype(str))
    .groupby("individual", sort=True)
    .size()
    .reset_index(name="n_cells")
)
cells_per_indiv.insert(0, "dataset", dataset)
cells_per_indiv.to_csv("cells-per-indiv.tsv", sep="\t", index=False)
