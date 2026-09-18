#!/usr/bin/env python3
"""Extract Ensembl feature_id → GeneSymbol from a Tremor h5ad."""

import sys

import anndata as ad
import pandas as pd

h5ad_path, out_path = sys.argv[1:3]

adata = ad.read_h5ad(h5ad_path, backed="r")
symbol_col = next(
    (c for c in ("GeneSymbol", "name", "feature_name", "gene_symbols") if c in adata.var.columns),
    None,
)
symbols = (
    adata.var[symbol_col].astype(str)
    if symbol_col is not None
    else adata.var_names.astype(str)
)

pd.DataFrame(
    {
        "feature_id": adata.var_names.astype(str),
        "gene_symbol": symbols,
    }
).to_csv(out_path, sep="\t", index=False)
