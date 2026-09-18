#!/usr/bin/env python3

import sys
from pathlib import Path

import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from cell_label_subset import cell_label_mask  # noqa: E402

cell_type = sys.argv[1]
adata = sc.read_h5ad(sys.argv[2], backed="r")

subset = adata[
    cell_label_mask(adata.obs["cell_label"], cell_type).to_numpy()
].to_memory()

sc.pp.filter_genes(subset, min_cells=3)
sc.pp.normalize_total(subset)
sc.pp.log1p(subset)

sc.tl.pca(subset)

subset.obs = subset.obs[["cell_label", "individual"]].copy()

subset.write_h5ad(f"{cell_type}-seacells-input.h5ad")
