#!/usr/bin/env python3
# Build Castonguay-style TensorQTL covariates for Tremor PB INT phenotypes:
# sex, age, condition, seq_batch, geno PC1-3, expression PC1-30.

import sys

import numpy as np
import pandas as pd
import scanpy as sc

pheno_path, geno_pcs_path, h5ad_path, out_path = sys.argv[1:5]

pheno = pd.read_csv(pheno_path, sep="\t", index_col=0)
donors = pheno.columns.astype(str).tolist()

Y = pheno.to_numpy(dtype=float).T
Y = Y - Y.mean(0, keepdims=True)
U, S, _ = np.linalg.svd(Y, full_matrices=False)
n_pc = min(30, U.shape[1], max(1, U.shape[0] - 1))
expr_pcs = U[:, :n_pc] * S[:n_pc]

geno = pd.read_csv(geno_pcs_path, sep="\t", dtype={"sample_id": str}).set_index(
    "sample_id"
)
gcols = [c for c in geno.columns if c.startswith("geno_pc")][:3]

adata = sc.read_h5ad(h5ad_path, backed="r")
meta = (
    adata.obs[
        [
            "donor_id",
            "age",
            "sex_ontology_term_id",
            "disease_ontology_term_id",
            "seq_batch",
        ]
    ]
    .drop_duplicates("donor_id")
    .assign(donor_id=lambda d: d.donor_id.astype(str))
    .set_index("donor_id")
    .reindex(donors)
)

# Authors: sex F=1 M=2; condition case=1 control=2; seq_batch WT=1 else 2
sex = np.where(meta["sex_ontology_term_id"].astype(str) == "PATO:0000383", 1.0, 2.0)
condition = np.where(
    meta["disease_ontology_term_id"].astype(str) == "MONDO_0003233", 1.0, 2.0
)
seq_batch = np.where(meta["seq_batch"].astype(str) == "WT", 1.0, 2.0)
age = meta["age"].astype(float).to_numpy()
age = (age - np.nanmean(age)) / np.nanstd(age)

cov = pd.DataFrame(
    {
        "sample_id": donors,
        "condition": condition,
        "sex": sex,
        "age": age,
        "seq_batch": seq_batch,
    }
)
for c in gcols:
    cov[c] = geno.reindex(donors)[c].to_numpy()
for i in range(n_pc):
    cov[f"PC_{i + 1}"] = expr_pcs[:, i]
cov.to_csv(out_path, sep="\t", index=False)
