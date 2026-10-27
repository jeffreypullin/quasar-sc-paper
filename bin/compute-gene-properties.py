#! /usr/bin/env python3

import sys
from pathlib import Path

import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from apply_cell_selection import apply_cell_selection  # noqa: E402
from binomial_downsample_counts import binomial_downsample_counts  # noqa: E402
from cell_label_subset import cell_label_mask  # noqa: E402

def compute_gene_stats(X):
    n_obs, n_vars = X.shape

    if n_obs == 0:
        mean = np.zeros(n_vars, dtype=float)
        var = np.zeros(n_vars, dtype=float)
        cv = np.full(n_vars, np.nan, dtype=float)
        nonzero_frac = np.zeros(n_vars, dtype=float)
        return mean, var, cv, nonzero_frac

    if sp.issparse(X):
        X = X.tocsr()
        mean = np.asarray(X.mean(axis=0)).ravel()
        ex2 = np.asarray(X.multiply(X).mean(axis=0)).ravel()
        var = ex2 - mean**2
        var[var < 0] = 0.0
        nonzero_frac = np.asarray(X.getnnz(axis=0)).ravel() / float(n_obs)
    else:
        X = np.asarray(X)
        mean = X.mean(axis=0)
        var = X.var(axis=0)
        nonzero_frac = (X != 0).mean(axis=0)

    sd = np.sqrt(var)
    cv = np.full(mean.shape, np.nan, dtype=float)
    nonzero_mean = mean != 0
    cv[nonzero_mean] = sd[nonzero_mean] / mean[nonzero_mean]

    return mean, var, cv, nonzero_frac

cell_label = sys.argv[1]
cell_frac = float(sys.argv[2])
indiv_frac = float(sys.argv[3])
adata = sc.read_h5ad(sys.argv[4], backed="r")
anno_df = pd.read_csv(sys.argv[5], sep="\t")
min_n_cells = int(sys.argv[6])
count_frac = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0

subset = adata[
    cell_label_mask(adata.obs["cell_label"], cell_label).to_numpy()
].to_memory()

subset = apply_cell_selection(
    subset,
    cell_frac=cell_frac,
    indiv_frac=indiv_frac,
    min_n_cells=min_n_cells,
    seed_key=cell_label,
    individual_col="individual",
    copy=True,
)

if subset.n_obs == 0:
    raise ValueError(
        "No cells found for cell_label="
        f"{cell_label} after filtering/subsampling "
        f"cell_frac={cell_frac}, indiv_frac={indiv_frac}, "
        f"min_n_cells={min_n_cells}, count_frac={count_frac}"
    )

# Expression before count thinning — used to fix the gene set across count_frac.
_, _, _, sc_non_zero_frac_pre_count = compute_gene_stats(subset.X)

subset = binomial_downsample_counts(
    subset,
    count_frac,
    seed_key=cell_label,
    copy=True,
)

sc_mean, sc_var, sc_cv, sc_non_zero_frac = compute_gene_stats(subset.X)

pbs = []
for indiv in subset.obs["individual"].unique():
    indiv_cells = subset[subset.obs["individual"] == indiv, :]
    X_rep = indiv_cells.X.sum(axis=0)

    rep_adata = sc.AnnData(X=X_rep, var=subset.var[[]].copy())
    rep_adata.obs_names = [str(indiv)]
    rep_adata.obs["individual"] = str(indiv)
    pbs.append(rep_adata)

pb = sc.concat(pbs)
pb_mean, pb_var, pb_cv, pb_non_zero_frac = compute_gene_stats(pb.X)

df = pd.DataFrame(
    {
        "feature_id": subset.var_names.to_numpy(),
        "sc_mean": sc_mean,
        "sc_var": sc_var,
        "sc_cv": sc_cv,
        "sc_non_zero_frac": sc_non_zero_frac,
        "sc_non_zero_frac_pre_count": sc_non_zero_frac_pre_count,
        "pb_mean": pb_mean,
        "pb_var": pb_var,
        "pb_cv": pb_cv,
        "pb_non_zero_frac": pb_non_zero_frac,
    }
    )

# Remove non-autosomal genes from all downstream analysis.
df = df[df["feature_id"].isin(anno_df["phenotype_id"])]

df.to_csv(f"{cell_label}-gene-properties.tsv", sep="\t", index=False)
