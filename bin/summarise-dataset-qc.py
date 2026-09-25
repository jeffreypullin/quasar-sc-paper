#!/usr/bin/env python3

from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import scanpy as sc

dataset = sys.argv[1]
adata = sc.read_h5ad(sys.argv[2])

n_counts = np.asarray(adata.X.sum(axis=1)).ravel()

pd.DataFrame({"dataset": dataset, "n_counts": n_counts}).to_csv(
    "cell-counts.tsv", sep="\t", index=False
)

cells_per_indiv = (
    adata.obs[["individual"]]
    .assign(individual=lambda df: df["individual"].astype(str))
    .groupby("individual", sort=True)
    .size()
    .reset_index(name="n_cells")
)
cells_per_indiv.insert(0, "dataset", dataset)
cells_per_indiv.to_csv("cells-per-indiv.tsv", sep="\t", index=False)
