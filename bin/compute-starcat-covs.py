#! /usr/bin/env python3
# pyright: reportMissingImports=false

import sys
from pathlib import Path

import matplotlib
import pandas as pd
import scanpy as sc
from starcat import starCAT

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from cell_label_subset import cell_label_mask  # noqa: E402

score_cols = [
    "Cytotoxic",
    "Th22",
    "MAIT",
    "TEMRA",
    "CD4-CM",
    "CD8-EM",
    "CD4-Naive",
    "Th2-Activated",
    "Th2-Resting",
    "Th1-Like",
    "CD8-Trm",
    "Th17-Activated",
    "Tfh-2",
    "Tph",
    "Exhaustion",
    "Tfh-1",
    "CellCycle-S",
    "CellCycle-G2M",
]
cov_names = {col: "starcat_" + col.replace("-", "_") for col in score_cols}

cell_label = sys.argv[1]
adata = sc.read_h5ad(sys.argv[2])
starcat_cache = sys.argv[3] if len(sys.argv) > 3 else "cache"

subset = adata[cell_label_mask(adata.obs["cell_label"], cell_label).to_numpy(), :].copy()

gene_symbols = subset.var["GeneSymbol"].astype(str).str.strip()
keep_genes = (gene_symbols != "") & (gene_symbols.str.lower() != "nan")
subset = subset[:, keep_genes.to_numpy()].copy()
subset.var_names = gene_symbols[keep_genes].to_numpy()
subset.var_names_make_unique()

usage, _ = starCAT(reference="TCAT.V1", cachedir=starcat_cache).fit_transform(subset)

covs = usage[score_cols].copy()
covs = covs.apply(pd.to_numeric, errors="coerce")
covs.rename(columns=cov_names, inplace=True)
covs.insert(0, "cell_id", covs.index)
covs.to_csv(f"{cell_label}-starcat-covs.tsv", sep="\t", index=False)

n_plot = min(50000, subset.n_obs)
plot_cells = (
    pd.Series(subset.obs_names)
    .sample(n=n_plot, random_state=1)
    .to_numpy()
)
plot_subset = subset[plot_cells, :].copy()
plot_labels = plot_subset.obs["cell_label"].astype(str)
label_order = plot_labels.value_counts().index.tolist()
sc.pp.normalize_total(plot_subset, target_sum=1e4)
sc.pp.log1p(plot_subset)
sc.tl.pca(plot_subset)
coords = plot_subset.obsm["X_pca"][:, :2]
plot_usage = usage.reindex(plot_cells)[score_cols]

n_panels = 1 + len(score_cols)
n_cols = 3
n_rows = (n_panels + n_cols - 1) // n_cols
fig, axes = plt.subplots(
    n_rows,
    n_cols,
    figsize=(7 * n_cols, 6 * n_rows),
    constrained_layout=True,
)
axes_flat = axes.flatten()

label_ax = axes_flat[0]
label_cmap = plt.get_cmap("tab20", max(len(label_order), 1))
for idx, label in enumerate(label_order):
    mask = (plot_labels == label).to_numpy()
    label_ax.scatter(
        coords[mask, 0],
        coords[mask, 1],
        s=2,
        color=label_cmap(idx),
        label=label,
        linewidths=0,
    )
label_ax.set_title("Cell label")
label_ax.set_xlabel("PC 1")
label_ax.set_ylabel("PC 2")
label_ax.legend(markerscale=4, fontsize="small", frameon=False)

for i, col in enumerate(score_cols):
    ax = axes_flat[i + 1]
    vals = pd.to_numeric(plot_usage[col], errors="coerce").to_numpy()
    pcm = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=vals,
        s=2,
        cmap="viridis",
        linewidths=0,
    )
    ax.set_title(col)
    ax.set_xlabel("PC 1")
    ax.set_ylabel("PC 2")
    fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.04)

for j in range(n_panels, len(axes_flat)):
    axes_flat[j].set_visible(False)

fig.savefig(f"{cell_label}-starcat-score-plots.pdf")
