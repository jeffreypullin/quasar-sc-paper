#! /usr/bin/env python3

import sys
import os
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from subsample_cells import subsample_cells_within_individuals  # noqa: E402
from subsample_individuals import subsample_individuals  # noqa: E402
from cell_label_subset import cell_label_mask  # noqa: E402

def flatten(xss):
    return [x for xs in xss for x in xs]

cell_label = sys.argv[1]
cell_frac = float(sys.argv[2])
indiv_frac = float(sys.argv[3])
adata = sc.read_h5ad(sys.argv[4])
gene_prop_data = pd.read_csv(sys.argv[5], sep="\t")

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

keep_genes = gene_prop_data.loc[gene_prop_data['sc_non_zero_frac'] > 0.005, 'feature_id'].values
row_idx = np.flatnonzero(
    cell_label_mask(adata.obs["cell_label"], cell_label).to_numpy()
)
adata_rows = adata[row_idx, :].copy()
cell_type_subset = adata_rows[:, list(keep_genes)].copy()

if indiv_frac != 1.0:
    cell_type_subset = subsample_individuals(
        cell_type_subset,
        indiv_frac,
        seed_key=cell_label,
        individual_col="individual",
        copy=True,
    )

if cell_frac != 1.0:
    cell_type_subset = subsample_cells_within_individuals(
        cell_type_subset,
        cell_frac,
        seed_key=cell_label,
        individual_col="individual",
        copy=True,
    )

pbs = []
for indiv in cell_type_subset.obs.individual.unique():
    indiv_cell_subset = cell_type_subset[cell_type_subset.obs['individual'] == indiv]
    
    rep_adata = sc.AnnData(X = indiv_cell_subset.X.mean(axis = 0),
                           var = indiv_cell_subset.var[[]])
    rep_adata.obs_names = [indiv]
    rep_adata.obs['individual'] = indiv
    pbs.append(rep_adata)

pb = sc.concat(pbs)

gene_ids = pb.var.index.tolist()
indiv_ids = pb.obs['individual'].tolist()

first_line = indiv_ids
first_line.insert(0, 'feature_id')
file_name = f"{cell_label}-pb-pheno.tsv"

with open(file_name, 'w') as fsum:
    fsum.write('\t'.join(first_line) + '\n') 
    for id in gene_ids:
        gene_res_list = flatten(pb[:, id].X.toarray().tolist())
        line = gene_res_list
        line.insert(0, id)
        fsum.write('\t'.join(map(str, line)) + '\n') 
