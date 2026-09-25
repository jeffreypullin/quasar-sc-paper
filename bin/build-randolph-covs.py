#!/usr/bin/env python3
"""Build Randolph paper-matched covariates for QuASAR log-count mapping.

Covariates (paper Methods / MatrixEQTL scripts):
  - scaled age
  - condition/cell-type bimodality proportion (geneProp)
  - genotype PC1–2
  - condition/cell-type-specific expression PCs (counts from submit_*_mapping.sh)

Sex is omitted (male-only cohort). A constant sex placeholder would be collinear
with the intercept QuASAR adds, so every score-test p-value becomes NaN.
The ctc_metadata argument is accepted for API stability but unused.

Usage:
  build-randolph-covs.py \\
      <pb_pheno.tsv> <geno_pcs.tsv> <metadata.csv> <ctc_metadata.csv> \\
      <cell_type> <condition> <out.tsv>
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

# Expression PCs to retain: (flu, NI) from authors' submit_*_mapping.sh.
EXPR_PCS = {
    "B": (3, 6),
    "CD4_T": (2, 4),
    "CD8_T": (4, 6),
    "monocytes_combined": (7, 10),
    "NK_combined": (2, 2),
}

GENEPROP_COL = {
    "B": "B_geneProp",
    "CD4_T": "CD4_T_geneProp",
    "CD8_T": "CD8_T_geneProp",
    "monocytes_combined": "monocytes_combined_geneProp",
    "NK_combined": "NK_combined_geneProp",
}


def main() -> None:
    if len(sys.argv) < 8:
        raise SystemExit(
            "Usage: build-randolph-covs.py <pb_pheno> <geno_pcs> <meta> "
            "<ctc_meta> <cell_type> <condition> <out>"
        )
    pheno_path, geno_path, meta_path, ctc_path, cell_type, condition, out_path = (
        sys.argv[1:8]
    )

    if cell_type not in EXPR_PCS:
        raise SystemExit(f"Unknown cell_type for Randolph PC map: {cell_type}")
    if condition not in ("NI", "flu"):
        raise SystemExit(f"condition must be NI or flu, got {condition}")

    pheno = pd.read_csv(pheno_path, sep="\t", index_col=0)
    donors = pheno.columns.astype(str).tolist()

    n_flu, n_ni = EXPR_PCS[cell_type]
    n_pc = n_flu if condition == "flu" else n_ni

    Y = pheno.to_numpy(dtype=float).T
    Y = Y - Y.mean(0, keepdims=True)
    # PCA via SVD of samples × genes.
    u, s, _ = np.linalg.svd(Y, full_matrices=False)
    n_pc = min(n_pc, u.shape[1], max(1, u.shape[0] - 1))
    expr_pcs = u[:, :n_pc] * s[:n_pc]

    geno = pd.read_csv(geno_path, sep="\t", dtype={"sample_id": str}).set_index(
        "sample_id"
    )
    gcols = [c for c in geno.columns if c.startswith("geno_pc")][:2]

    meta = pd.read_csv(meta_path)
    # One row per infection_ID; restrict to this condition then key by donor.
    meta = meta.loc[meta["infection_status"] == condition].copy()
    meta["indiv_ID"] = meta["indiv_ID"].astype(str)
    meta = meta.drop_duplicates("indiv_ID").set_index("indiv_ID")

    age = meta.reindex(donors)["age_Scale"].astype(float).to_numpy()
    gp_col = GENEPROP_COL[cell_type]
    if gp_col is not None and gp_col in meta.columns:
        gene_prop = meta.reindex(donors)[gp_col].astype(float).to_numpy()
    else:
        gene_prop = np.zeros(len(donors), dtype=float)

    cov = pd.DataFrame({"sample_id": donors, "age": age, "geneProp": gene_prop})
    for c in gcols:
        cov[c] = geno.reindex(donors)[c].to_numpy()
    for i in range(n_pc):
        cov[f"PC_{i + 1}"] = expr_pcs[:, i]

    cov.to_csv(out_path, sep="\t", index=False)
    print(
        f"Wrote {out_path}: {len(donors)} samples, "
        f"{n_pc} expr PCs, condition={condition}, cell_type={cell_type}"
    )


if __name__ == "__main__":
    main()
