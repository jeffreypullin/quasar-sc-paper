#!/usr/bin/env python3
"""Join Tremor TensorQTL donor-level covariates onto single-cell rows.

Usage:
  expand-tremor-sc-covs.py <tqtl_covs.tsv> <sc_pheno.tsv> <out.tsv>
"""

from __future__ import annotations

import sys

import pandas as pd


def main() -> None:
    cov_path, sc_path, out_path = sys.argv[1:4]
    cov = pd.read_csv(cov_path, sep="\t", dtype={"sample_id": str})
    sc = pd.read_csv(sc_path, sep="\t", usecols=[0, 1], dtype=str)
    sc.columns = ["sample_id", "cell_id"]
    out = sc.merge(cov, on="sample_id", how="inner")
    cols = ["sample_id", "cell_id"] + [c for c in cov.columns if c != "sample_id"]
    out[cols].to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
