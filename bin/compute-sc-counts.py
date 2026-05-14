#! /usr/bin/env python3

import sys
import os
from pathlib import Path
import numpy as np
import scanpy as sc
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from subsample_cells import subsample_cells_within_individuals  # noqa: E402
from subsample_individuals import subsample_individuals  # noqa: E402
from cell_label_subset import cell_label_mask  # noqa: E402

cell_type = sys.argv[1]
cell_frac = float(sys.argv[2])
indiv_frac = float(sys.argv[3])
adata = sc.read_h5ad(sys.argv[4])
gene_prop_data = pd.read_csv(sys.argv[5], sep="\t")

keep_genes = gene_prop_data.loc[gene_prop_data['sc_non_zero_frac'] > 0.005, 'feature_id'].values
row_idx = np.flatnonzero(
    cell_label_mask(adata.obs["cell_label"], cell_type).to_numpy()
)
adata_rows = adata[row_idx, :].copy()
cell_type_subset = adata_rows[:, list(keep_genes)].copy()

if indiv_frac != 1.0:
    cell_type_subset = subsample_individuals(
        cell_type_subset,
        indiv_frac,
        seed_key=cell_type,
        individual_col="individual",
        copy=True,
    )

if cell_frac != 1.0:
    cell_type_subset = subsample_cells_within_individuals(
        cell_type_subset,
        cell_frac,
        seed_key=cell_type,
        individual_col="individual",
        copy=True,
    )

Xf = cell_type_subset.X
if sp.issparse(Xf):
    Xf = Xf.toarray()

counts_df = pd.DataFrame(
    Xf,
    index=cell_type_subset.obs['individual'],
    columns=cell_type_subset.var_names
)

counts_df.insert(0, 'cell_id', cell_type_subset.obs.index.to_numpy())
counts_df.sort_index(inplace=True)
counts_df.index.name = 'sample_id'

counts_df.to_csv(f"{cell_type}-sc-pheno.tsv", sep='\t')
