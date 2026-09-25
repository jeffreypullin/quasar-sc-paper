#!/usr/bin/env python3
"""Join Randolph donor-level covariates onto single-cell phenotype rows.

Usage:
  expand-randolph-sc-covs.py <pb_covs.tsv> <sc_pheno.tsv> <out.tsv>
"""

from __future__ import annotations

import sys

import pandas as pd


def main() -> None:
    pb_path, sc_path, out_path = sys.argv[1:4]
    pb = pd.read_csv(pb_path, sep="\t", dtype={"sample_id": str})
    pb = pb.drop(columns=["sex"], errors="ignore")
    sc = pd.read_csv(sc_path, sep="\t", usecols=[0, 1], dtype=str)
    sc.columns = ["sample_id", "cell_id"]
    out = sc.merge(pb, on="sample_id", how="inner")
    cols = ["sample_id", "cell_id"] + [c for c in pb.columns if c != "sample_id"]
    out[cols].to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
