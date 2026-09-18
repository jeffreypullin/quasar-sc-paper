#! /usr/bin/env python3
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from cell_label_subset import cell_label_mask  # noqa: E402
from subsample_cells import subsample_cells_within_individuals  # noqa: E402
from subsample_individuals import subsample_individuals  # noqa: E402

cell_type = sys.argv[1]
cell_frac = float(sys.argv[2])
indiv_frac = float(sys.argv[3])
adata_path = sys.argv[4]
genes_tsv = Path(sys.argv[5])
out_path = sys.argv[6] if len(sys.argv) > 6 else f"{cell_type}-sc-logcounts.tsv"

genes_df = pd.read_csv(genes_tsv, sep="\t")
if "feature_id" not in genes_df.columns:
    sys.exit(f"genes TSV must contain feature_id column; got: {', '.join(genes_df.columns)}")

if "cell_type" in genes_df.columns:
    genes_df = genes_df.loc[genes_df["cell_type"].astype(str) == cell_type]

GENES = list(dict.fromkeys(genes_df["feature_id"].dropna().astype(str).tolist()))
if not GENES:
    sys.exit(f"no feature_id values found in genes TSV for cell type {cell_type}")

adata = sc.read_h5ad(adata_path, backed="r")
cell_type_subset = adata[
    cell_label_mask(adata.obs["cell_label"], cell_type).to_numpy()
].to_memory()
sc.pp.normalize_total(cell_type_subset)
sc.pp.log1p(cell_type_subset)

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
