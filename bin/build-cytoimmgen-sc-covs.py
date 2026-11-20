#!/usr/bin/env python3
"""Build Cytoimmgen single-cell covariate files for QuASAR.

Covariates: sex, age, culture-batch one-hots (drop first level), 15 pseudobulk
expression PCs, genotype PCs 1-6. Interaction specs additionally include the
tested term as a main-effect covariate.

Usage:
  build-cytoimmgen-sc-covs.py \\
      <sc_pheno.tsv> <pb_expr_covs.tsv> <geno_pcs.txt> <harmonised.h5ad> \\
      <cell_type> <cov_spec> <out.tsv> [timepoints_csv]

cov_spec is one of: cyto, cyto+timepoint, cyto+pseudotime, cyto+division_state
"""

from __future__ import annotations

import sys

import pandas as pd
import scanpy as sc


def one_hot(series: pd.Series, prefix: str, drop_first: bool = True) -> pd.DataFrame:
    levels = sorted(series.astype(str).unique().tolist())
    if drop_first and levels:
        levels = levels[1:]
    out = pd.DataFrame(0, index=series.index, columns=[f"{prefix}_{lv}" for lv in levels])
    for lv in levels:
        out[f"{prefix}_{lv}"] = (series.astype(str) == lv).astype(int)
    return out


def main() -> None:
    (
        sc_pheno_path,
        pb_cov_path,
        geno_pcs_path,
        h5ad_path,
        cell_type,
        cov_spec,
        out_path,
    ) = sys.argv[1:8]
    timepoints = (
        sys.argv[8].split(",") if len(sys.argv) > 8 and sys.argv[8] else None
    )

    cells = pd.read_csv(sc_pheno_path, sep="\t", usecols=[0, 1], dtype=str)
    cells.columns = ["sample_id", "cell_id"]

    pb = pd.read_csv(pb_cov_path, sep="\t", dtype={"sample_id": str})
    geno = pd.read_csv(geno_pcs_path, sep="\t", dtype={"sample_id": str})

    pc_cols = [f"PC_{i}" for i in range(1, 16)]
    missing_pcs = [c for c in pc_cols if c not in pb.columns]
    if missing_pcs:
        raise SystemExit(f"pb expr covs missing {missing_pcs}; have {list(pb.columns)}")
    geno_cols = [f"geno_pc{i}" for i in range(1, 7)]
    missing_geno = [c for c in geno_cols if c not in geno.columns]
    if missing_geno:
        raise SystemExit(
            f"geno PCs missing {missing_geno}; have {list(geno.columns)}"
        )

    donor = pb[["sample_id", "sex", "age", *pc_cols]].merge(
        geno[["sample_id", *geno_cols]], on="sample_id", how="inner"
    )

    adata = sc.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.copy()
    obs["cell_id"] = obs.index.astype(str)
    need = ["individual", "batch", "timepoint", "pseudotime", "division_state"]
    missing_obs = [c for c in need if c not in obs.columns]
    if missing_obs:
        raise SystemExit(f"harmonised h5ad missing {missing_obs}")

    cell_meta = obs.loc[
        obs["cell_id"].isin(cells["cell_id"]),
        ["cell_id", "individual", "batch", "timepoint", "pseudotime", "division_state"],
    ].copy()
    cell_meta["individual"] = cell_meta["individual"].astype(str)

    out = cells.merge(cell_meta, on="cell_id", how="inner")
    out = out.merge(donor, left_on="individual", right_on="sample_id", how="inner")
    # sample_id must be the donor id for QuASAR.
    out["sample_id"] = out["individual"]

    batch_oh = one_hot(out["batch"], "batch", drop_first=True)
    out = pd.concat([out.reset_index(drop=True), batch_oh.reset_index(drop=True)], axis=1)

    base_cols = (
        ["sample_id", "cell_id", "sex", "age"]
        + list(batch_oh.columns)
        + pc_cols
        + geno_cols
    )

    if cov_spec == "cyto":
        cols = base_cols
    elif cov_spec == "cyto+timepoint":
        # Ordered timepoints from config when provided; else sorted unique values.
        levels = timepoints if timepoints else sorted(out["timepoint"].astype(str).unique())
        # Drop first level as reference; remaining levels become GxE columns.
        dummies = levels[1:]
        for lv in dummies:
            out[f"timepoint_{lv}"] = (out["timepoint"].astype(str) == lv).astype(int)
        cols = base_cols + [f"timepoint_{lv}" for lv in dummies]
    elif cov_spec == "cyto+pseudotime":
        out["pseudotime"] = pd.to_numeric(out["pseudotime"], errors="coerce")
        cols = base_cols + ["pseudotime"]
    elif cov_spec == "cyto+division_state":
        out["division_state"] = pd.to_numeric(out["division_state"], errors="coerce")
        cols = base_cols + ["division_state"]
    else:
        raise SystemExit(f"Unknown cov_spec: {cov_spec}")

    out = out.dropna(subset=[c for c in cols if c not in ("sample_id", "cell_id")])
    out = out.sort_values(["sample_id", "cell_id"])
    out[cols].to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {out_path} ({out.shape[0]} cells, cov_spec={cov_spec})")


if __name__ == "__main__":
    main()
