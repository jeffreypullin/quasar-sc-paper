#! /usr/bin/env python3

import sys
from pathlib import Path
import numpy as np
import scanpy as sc
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from subsample_cells import subsample_cells_within_individuals  # noqa: E402
from subsample_individuals import subsample_individuals  # noqa: E402
from downsample_individuals_to_n_cells import downsample_individuals_to_n_cells  # noqa: E402
from binomial_downsample_counts import binomial_downsample_counts  # noqa: E402
from cell_label_subset import cell_label_mask  # noqa: E402

cell_type = sys.argv[1]
cell_frac = float(sys.argv[2])
indiv_frac = float(sys.argv[3])
adata = sc.read_h5ad(sys.argv[4])
gene_prop_data = pd.read_csv(sys.argv[5], sep="\t")
n_cells_target = int(sys.argv[6])
count_frac = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0

keep_genes = gene_prop_data.loc[gene_prop_data['sc_non_zero_frac'] > 0.01, 'feature_id'].values
row_idx = np.flatnonzero(
    cell_label_mask(adata.obs["cell_label"], cell_type).to_numpy()
)
adata_rows = adata[row_idx, :].copy()
cell_type_subset = adata_rows[:, list(keep_genes)].copy()

cell_type_subset = downsample_individuals_to_n_cells(
    cell_type_subset,
    n_cells_target,
    individual_col="individual",
    seed_key=cell_type,
    copy=True,
)

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

cell_type_subset = binomial_downsample_counts(
    cell_type_subset,
    count_frac,
    seed_key=cell_type,
    copy=True,
)

sc.pp.normalize_total(cell_type_subset)
sc.pp.log1p(cell_type_subset)

Xf = cell_type_subset.X
if sp.issparse(Xf):
    Xf = Xf.toarray()

logcounts_df = pd.DataFrame(
    Xf,
    index=cell_type_subset.obs['individual'],
    columns=cell_type_subset.var_names
)

logcounts_df.insert(0, 'cell_id', cell_type_subset.obs.index.to_numpy())
logcounts_df.sort_index(inplace=True)
logcounts_df.index.name = 'sample_id'

logcounts_df.to_csv(f"{cell_type}-sc-logcounts.tsv", sep='\t')
