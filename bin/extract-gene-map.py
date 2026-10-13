#!/usr/bin/env python3
"""Extract Ensembl feature_id → gene symbol from an h5ad."""

import sys

import anndata as ad
import pandas as pd

h5ad_path, out_path = sys.argv[1:3]

adata = ad.read_h5ad(h5ad_path, backed="r")

pd.DataFrame(
    {
        "feature_id": adata.var_names.astype(str),
        "gene_symbol": adata.var["GeneSymbol"].astype(str).str.strip(),
    }
).to_csv(out_path, sep="\t", index=False)
