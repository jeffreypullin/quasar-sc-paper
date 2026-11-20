#!/usr/bin/env python3
"""Print Cytoimmgen h5ad and predictions CSV column names for config fill-in.

Usage on Sanger:
  python bin/inspect-cytoimmgen.py \\
    /lustre/scratch124/humgen/projects_v2/healthy_imm_expr/cytoimmgen/data/inhouse/scrna/processed/from_bucket/h5ad/cytoimmgen_qc.h5ad \\
    /lustre/scratch124/humgen/projects_v2/healthy_imm_expr/CellDivider/analysis/prediction_model/models/output/predictions_cytoimmgen_df.csv
"""

from __future__ import annotations

import sys

import pandas as pd
import scanpy as sc


def main() -> None:
    h5ad_path = sys.argv[1]
    preds_path = sys.argv[2] if len(sys.argv) > 2 else None

    adata = sc.read_h5ad(h5ad_path, backed="r")
    print("=== h5ad obs columns ===")
    for col in adata.obs.columns:
        nuniq = adata.obs[col].nunique(dropna=False)
        sample = adata.obs[col].dropna().astype(str).unique()[:8].tolist()
        print(f"  {col}\tnuniq={nuniq}\texamples={sample}")
    print("\n=== h5ad var columns ===")
    for col in adata.var.columns:
        print(f"  {col}")
    print(f"\nvar_names example: {list(adata.var_names[:5])}")
    print(f"obs_names example: {list(adata.obs_names[:5])}")
    print(f"n_obs={adata.n_obs} n_vars={adata.n_vars}")

    if preds_path is None:
        return
    preds = pd.read_csv(preds_path, nrows=5)
    print("\n=== predictions CSV columns ===")
    print(list(preds.columns))
    print(preds.head().to_string())


if __name__ == "__main__":
    main()
