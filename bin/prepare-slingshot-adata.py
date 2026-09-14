#!/usr/bin/env python3

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from cell_label_subset import cell_label_mask  # noqa: E402

cell_type = sys.argv[1]
adata = sc.read_h5ad(sys.argv[2])

subset = adata[
    cell_label_mask(adata.obs["cell_label"], cell_type), :
].copy()

sc.pp.filter_genes(subset, min_cells=3)
sc.pp.normalize_total(subset)
sc.pp.log1p(subset)

sc.tl.pca(subset)

# Export PCA + labels as TSV. anndataR cannot reliably read scanpy h5ad
# (nullable-string-array index + sparse X → SummarizedExperiment errors).
n_pcs = subset.obsm["X_pca"].shape[1]
pca = pd.DataFrame(
    np.asarray(subset.obsm["X_pca"]),
    index=subset.obs_names.astype(str),
    columns=[f"PC{i}" for i in range(1, n_pcs + 1)],
)
out = pca.copy()
out.insert(0, "cell_label", subset.obs["cell_label"].astype(str).to_numpy())
out.insert(0, "cell_id", out.index)
out.to_csv(f"{cell_type}-slingshot-input.tsv", sep="\t", index=False)
