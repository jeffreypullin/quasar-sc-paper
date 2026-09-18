#!/usr/bin/env python3
# Pseudobulk raw counts → TMM → filter mean CPM > 6 → inverse-normal transform.
# Matches Castonguay et al. / Bryois phenotype processing used for TensorQTL.

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from binomial_downsample_counts import binomial_downsample_counts  # noqa: E402
from cell_label_subset import cell_label_mask  # noqa: E402
from downsample_individuals_to_n_cells import downsample_individuals_to_n_cells  # noqa: E402
from subsample_cells import subsample_cells_within_individuals  # noqa: E402
from subsample_individuals import subsample_individuals  # noqa: E402

cell_label = sys.argv[1]
cell_frac = float(sys.argv[2])
indiv_frac = float(sys.argv[3])
adata = sc.read_h5ad(sys.argv[4], backed="r")
gene_prop_data = pd.read_csv(sys.argv[5], sep="\t")
n_cells_target = int(sys.argv[6])
count_frac = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0

keep_genes = gene_prop_data.loc[
    gene_prop_data["sc_non_zero_frac"] > 0.01, "feature_id"
].values
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

indivs = cell_type_subset.obs["individual"].astype(str).unique().tolist()
X = cell_type_subset.X
if sparse.issparse(X):
    X = X.tocsr()
n_genes = X.shape[1]
counts = np.zeros((len(indivs), n_genes), dtype=np.float64)
indiv_labels = cell_type_subset.obs["individual"].astype(str).to_numpy()
for i, indiv in enumerate(indivs):
    rows = np.where(indiv_labels == indiv)[0]
    counts[i] = np.asarray(X[rows].sum(axis=0)).ravel()

lib = counts.sum(1, keepdims=True)
lib[lib == 0] = 1
cpm0 = counts / lib * 1e6
keep = cpm0.mean(0) > 6
counts = counts[:, keep]
gene_ids = np.asarray(cell_type_subset.var_names)[keep]


def tmm_norm_factors(counts_mat):
    """edgeR-like TMM factors (Robinson & Oshlack)."""
    libsize = counts_mat.sum(1)
    ref = int(np.argsort(libsize)[len(libsize) // 2])
    factors = []
    for i in range(counts_mat.shape[0]):
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.log2(
                (counts_mat[i] / libsize[i]) / (counts_mat[ref] / libsize[ref])
            )
            a = 0.5 * np.log2(
                (counts_mat[i] / libsize[i]) * (counts_mat[ref] / libsize[ref])
            )
        ok = (
            np.isfinite(r)
            & np.isfinite(a)
            & (counts_mat[i] > 0)
            & (counts_mat[ref] > 0)
        )
        r, a = r[ok], a[ok]
        if len(r) < 10:
            factors.append(1.0)
            continue
        keep_a = (a > np.quantile(a, 0.05)) & (a < np.quantile(a, 0.95))
        r = r[keep_a]
        lo, hi = np.quantile(r, [0.3, 0.7])
        r = r[(r >= lo) & (r <= hi)]
        factors.append(float(2 ** np.mean(r)) if len(r) else 1.0)
    factors = np.asarray(factors, float)
    factors = factors / np.exp(np.mean(np.log(factors)))
    return factors


nf = tmm_norm_factors(counts)
cpm = counts / (counts.sum(1, keepdims=True) * nf[:, None]) * 1e6
n = cpm.shape[0]
pheno = np.empty_like(cpm)
for j in range(cpm.shape[1]):
    ranks = pd.Series(cpm[:, j]).rank(method="average").to_numpy()
    pheno[:, j] = norm.ppf(ranks / (n + 1.0))

out = pd.DataFrame(pheno.T, index=gene_ids, columns=indivs)
out.index.name = "feature_id"
out.to_csv(f"{cell_label}-pb-pheno.tsv", sep="\t")
