#!/usr/bin/env python3

import sys
import types


def _use_cli_tqdm() -> None:
    """SEACells imports tqdm.notebook; redirect to std tqdm for batch runs."""
    import tqdm.std as s

    m = types.ModuleType("tqdm.notebook")
    m.tqdm = m.trange = m.tqdm_notebook = s.tqdm
    sys.modules["tqdm.notebook"] = m


def _subsample_indices(individuals, n_cells: int, rng):
    """Stratified subsample by individual, returning up to n_cells indices."""
    import numpy as np

    inds = np.arange(len(individuals))
    if len(inds) <= n_cells:
        return inds

    frac = n_cells / len(inds)
    keep = []
    for ind in individuals.unique():
        idx = inds[individuals == ind]
        n_keep = max(1, round(len(idx) * frac))
        n_keep = min(n_keep, len(idx))
        keep.extend(rng.choice(idx, size=n_keep, replace=False))

    keep = np.array(keep)
    if len(keep) > n_cells:
        keep = rng.choice(keep, size=n_cells, replace=False)
    return np.sort(keep)


def _assign_to_nearest_metacell(pca_all, pca_fit, metacell_labels):
    """Assign each cell to the nearest metacell centroid in PCA space."""
    from scipy.spatial.distance import cdist

    groups = sorted(metacell_labels.unique())
    centroids = np.vstack([pca_fit[metacell_labels == g].mean(axis=0) for g in groups])
    nearest = cdist(pca_all, centroids).argmin(axis=1)
    return np.array(groups)[nearest]


_use_cli_tqdm()

import numpy as np
import pandas as pd
import scanpy as sc
import SEACells

MAX_CELLS = 30_000
RNG = np.random.default_rng(42)

h5ad_file = sys.argv[1]
cell_type = sys.argv[2]

adata = sc.read_h5ad(h5ad_file)

# Exclude ribosomal and mitochondrial genes from HVG candidates.
# Including them risks circular confounding: strong eQTLs at e.g. RPS26 shift
# cells in PCA space by genotype, causing spurious csaQTL signals.
exclude_genes = set(
    g for g in adata.var_names
    if g.startswith(("RPS", "RPL", "MRPS", "MRPL", "MT-"))
)
sc.pp.highly_variable_genes(adata, n_top_genes=2500, subset=False)
adata.var["highly_variable"] = (
    adata.var["highly_variable"] & ~adata.var_names.isin(exclude_genes)
)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

if adata.n_obs > MAX_CELLS:
    fit_idx = _subsample_indices(adata.obs["individual"], MAX_CELLS, RNG)
    adata_fit = adata[fit_idx].copy()
else:
    adata_fit = adata

n_seacells = max(1, round(adata_fit.n_obs / 1000))

model = SEACells.core.SEACells(
    adata_fit,
    build_kernel_on="X_pca",
    n_SEACells=n_seacells,
    n_waypoint_eigs=min(10, n_seacells - 1),
    convergence_epsilon=1e-3,
)
model.construct_kernel_matrix()
model.initialize_archetypes()
model.fit(min_iter=10, max_iter=100)

if adata_fit.n_obs < adata.n_obs:
    adata.obs["SEACell"] = _assign_to_nearest_metacell(
        adata.obsm["X_pca"],
        adata_fit.obsm["X_pca"],
        adata_fit.obs["SEACell"],
    )
else:
    adata.obs["SEACell"] = adata_fit.obs["SEACell"].values

pd.DataFrame({
    "group": adata.obs["SEACell"].values,
    "cell_id": adata.obs_names,
}).to_csv(f"{cell_type}-seacells-cell-groups.tsv", sep="\t", index=False)

# Compute UMAP on the full adata to obtain per-metacell coordinates.
sc.pp.neighbors(adata, use_rep="X_pca")
sc.tl.umap(adata)

umap_coords = pd.DataFrame(
    adata.obsm["X_umap"],
    index=adata.obs_names,
    columns=["umap_1", "umap_2"],
).assign(metacell=adata.obs["SEACell"].values)

metacell_umap = (
    umap_coords
    .groupby("metacell")[["umap_1", "umap_2"]]
    .mean()
)

n_cells = (
    adata.obs["SEACell"]
    .value_counts()
    .rename("n_cells")
)

(
    metacell_umap
    .join(n_cells)
    .sort_index()
    .rename_axis("metacell")
    .reset_index()
    .to_csv(f"{cell_type}-seacells-info.tsv", sep="\t", index=False)
)
