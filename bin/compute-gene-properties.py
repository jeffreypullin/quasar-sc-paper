#! /usr/bin/env python3

import sys
from pathlib import Path

import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from subsample_cells import subsample_cells_within_individuals  # noqa: E402
from subsample_individuals import subsample_individuals  # noqa: E402
from downsample_individuals_to_n_cells import downsample_individuals_to_n_cells  # noqa: E402
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
adata = sc.read_h5ad(sys.argv[4])
anno_df = pd.read_csv(sys.argv[5], sep="\t")
n_cells_target = int(sys.argv[6])
count_frac = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0

subset = adata[
    cell_label_mask(adata.obs["cell_label"], cell_label), :
].copy()

subset = downsample_individuals_to_n_cells(
    subset,
    n_cells_target,
    individual_col="individual",
    seed_key=cell_label,
    copy=True,
)

if indiv_frac != 1.0:
    subset = subsample_individuals(
        subset, indiv_frac, seed_key=cell_label, individual_col="individual", copy=True
    )

if cell_frac != 1.0:
    subset = subsample_cells_within_individuals(
        subset, cell_frac, seed_key=cell_label, individual_col="individual", copy=True
    )

subset = binomial_downsample_counts(
    subset,
    count_frac,
    seed_key=cell_label,
    copy=True,
)

if subset.n_obs == 0:
    raise ValueError(
        "No cells found for cell_label="
        f"{cell_label} after filtering/subsampling "
        f"cell_frac={cell_frac}, indiv_frac={indiv_frac}, "
        f"n_cells_target={n_cells_target}, count_frac={count_frac}"
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
        "pb_mean": pb_mean,
        "pb_var": pb_var,
        "pb_cv": pb_cv,
        "pb_non_zero_frac": pb_non_zero_frac,
    }
    )

# Remove non-autosomal genes from all downstream analysis.
df = df[df["feature_id"].isin(anno_df["phenotype_id"])]

df.to_csv(f"{cell_label}-gene-properties.tsv", sep="\t", index=False)
