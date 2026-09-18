#! /usr/bin/env python3

import sys
from pathlib import Path
import scanpy as sc
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from subsample_cells import subsample_cells_within_individuals  # noqa: E402
from subsample_individuals import subsample_individuals  # noqa: E402
from downsample_individuals_to_n_cells import downsample_individuals_to_n_cells  # noqa: E402
from binomial_downsample_counts import binomial_downsample_counts  # noqa: E402
from cell_label_subset import cell_label_mask  # noqa: E402

def flatten(xss):
    return [x for xs in xss for x in xs]

cell_label = sys.argv[1]
cell_frac = float(sys.argv[2])
indiv_frac = float(sys.argv[3])
adata = sc.read_h5ad(sys.argv[4], backed="r")
gene_prop_data = pd.read_csv(sys.argv[5], sep="\t")
n_cells_target = int(sys.argv[6])
count_frac = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0

keep_genes = gene_prop_data.loc[gene_prop_data['sc_non_zero_frac'] > 0.01, 'feature_id'].values
adata_rows = adata[
    cell_label_mask(adata.obs["cell_label"], cell_label).to_numpy()
].to_memory()
cell_type_subset = adata_rows[:, list(keep_genes)].copy()

cell_type_subset = downsample_individuals_to_n_cells(
    cell_type_subset,
    n_cells_target,
    individual_col="individual",
    seed_key=cell_label,
    copy=True,
)

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

cell_type_subset = binomial_downsample_counts(
    cell_type_subset,
    count_frac,
    seed_key=cell_label,
    copy=True,
)

sc.pp.normalize_total(cell_type_subset)
sc.pp.log1p(cell_type_subset)

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
