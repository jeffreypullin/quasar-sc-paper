#! /usr/bin/env python3

import sys
from pathlib import Path

import scanpy as sc
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from cell_label_subset import cell_label_mask  # noqa: E402

def flatten(xss):
    return [x for xs in xss for x in xs]

cell_label = sys.argv[1]
adata = sc.read_h5ad(sys.argv[2], backed="r")
adata = adata[
    cell_label_mask(adata.obs["cell_label"], cell_label).to_numpy()
].to_memory()

sc.pp.filter_genes(adata, min_cells=3)

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

pbs = []
for indiv in adata.obs.individual.unique():
    indiv_cell_subset = adata[adata.obs['individual'] == indiv]
    
    rep_adata = sc.AnnData(X = indiv_cell_subset.X.mean(axis = 0),
                           var = indiv_cell_subset.var[[]])
    rep_adata.layers['counts'] = indiv_cell_subset.layers['counts'].sum(axis = 0)
    rep_adata.obs_names = [indiv]
    rep_adata.obs['individual'] = indiv
    rep_adata.obs['sex'] = indiv_cell_subset.obs['sex'].iloc[0]
    rep_adata.obs['age'] = indiv_cell_subset.obs['age'].iloc[0]
    pbs.append(rep_adata)

pb = sc.concat(pbs)

# Compute PCA.
sc.tl.pca(pb)

# We only save the first 20 PCs.
pc_data = pb.obsm['X_pca'][0:, 0:20]

pc_df = pd.DataFrame(
  data=pc_data[0:,0:],
  index=pb.obs_names,
  columns=[f"PC_{i}" for i in range(1, 21)]
)

# Write out covariates.
pb_cov = pb.obs[['sex', 'age']].copy()
pb_cov['sex'] = pb_cov['sex'] - 1
pb_cov = pd.concat([pb_cov, pc_df], axis=1)
pb_cov.index.name = 'sample_id'
pb_cov.to_csv(f"{cell_label}-expr-covs.tsv", sep = "\t")
