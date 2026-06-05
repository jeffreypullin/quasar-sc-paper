#! /usr/bin/env python3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from cell_label_subset import cell_label_mask  # noqa: E402
from subsample_cells import subsample_cells_within_individuals  # noqa: E402
from subsample_individuals import subsample_individuals  # noqa: E402

GENES = [
    "ENSG00000197728",
    "ENSG00000196735",
    "ENSG00000143297",
    "ENSG00000179344",
    "ENSG00000215845",
    "ENSG00000171476",
    "ENSG00000196975",
    "ENSG00000153283",
    "ENSG00000138639",
    "ENSG00000132185",
    "ENSG00000127990",
    "ENSG00000184009",
    "ENSG00000114013",
    "ENSG00000224137",
    "ENSG00000033867",
    "ENSG00000185344",
    "ENSG00000183172",
    "ENSG00000204001",
    "ENSG00000196562",
    "ENSG00000132507",
    "ENSG00000161016",
]

cell_type = sys.argv[1]
cell_frac = float(sys.argv[2])
indiv_frac = float(sys.argv[3])
adata_path = sys.argv[4]
out_path = sys.argv[5] if len(sys.argv) > 5 else f"{cell_type}-sc-logcounts.tsv"

adata = sc.read_h5ad(adata_path)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

row_idx = np.flatnonzero(
    cell_label_mask(adata.obs["cell_label"], cell_type).to_numpy()
)
cell_type_subset = adata[row_idx, :].copy()

missing = [g for g in GENES if g not in cell_type_subset.var_names]
if missing:
    sys.exit(f"genes not in data: {', '.join(missing)}")

cell_type_subset = cell_type_subset[:, GENES].copy()

if indiv_frac != 1.0:
    cell_type_subset = subsample_individuals(
        cell_type_subset,
        indiv_frac,
        seed_key=cell_type,
        individual_col="individual",
        copy=True,
    )

if cell_frac != 1.0:
    cell_type_subset = subsample_cells_within_individuals(
        cell_type_subset,
        cell_frac,
        seed_key=cell_type,
        individual_col="individual",
        copy=True,
    )

Xf = cell_type_subset.X
if sp.issparse(Xf):
    Xf = Xf.toarray()

out_df = pd.DataFrame(
    Xf,
    index=cell_type_subset.obs["individual"],
    columns=GENES,
)
out_df.insert(0, "cell_id", cell_type_subset.obs.index.to_numpy())
out_df.sort_index(inplace=True)
out_df.index.name = "sample_id"
out_df.to_csv(out_path, sep="\t")
