#! /usr/bin/env python3

import sys

import scanpy as sc
import pandas as pd

cell_type = sys.argv[1]
groups_file = sys.argv[2]
adata = sc.read_h5ad(sys.argv[3])

groups = pd.read_csv(groups_file, sep="\t")

individual = adata.obs["individual"].reindex(groups["cell_id"]).to_numpy()
groups = groups.assign(individual=individual)

counts = (
    groups.groupby(["group", "individual"], observed=True)
    .size()
    .unstack(fill_value=0)
    .sort_index()
)
counts.index.name = "phenotype_id"

counts.to_csv(f"{cell_type}-csaqtl-pheno.tsv", sep="\t")
