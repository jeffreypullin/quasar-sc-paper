#!/usr/bin/env python3

from __future__ import annotations

import sys
import pandas as pd
import scanpy as sc

adata = sc.read_h5ad(sys.argv[1])
df = adata.obs[["cell_label", "individual"]].copy()
res = (
    df.groupby("cell_label", sort=True)
    .agg(n_cells=("cell_label", "size"), n_indiv=("individual", pd.Series.nunique))
    .reset_index()
 )

res.to_csv("cluster-sizes.tsv", sep="\t", index=False)
