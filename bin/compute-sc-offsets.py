#! /usr/bin/env python3

"""Build the three log-scale offset files used by the offset-precision experiment.

percell   : offset_ic = log(n_ic), the per-cell library size quasar computes internally.
donorflat : offset_ic = log(N_i / C_i), constant within a donor. Total exposure per
            donor, sum_c exp(offset_ic) = N_i, is preserved exactly, so this is the
            pseudobulk-equivalent offset spread uniformly over that donor's cells.
constant  : offset_ic = log(N / C), the same value for every cell, so exposure is
            proportional to cell count alone.
"""

import sys

import numpy as np
import pandas as pd

CHUNK_ROWS = 1000

cell_type = sys.argv[1]
sc_pheno_file = sys.argv[2]

sample_id_chunks = []
cell_id_chunks = []
total_chunks = []

for chunk in pd.read_csv(sc_pheno_file, sep="\t", chunksize=CHUNK_ROWS):
    gene_cols = chunk.columns.drop(["sample_id", "cell_id"])
    sample_id_chunks.append(chunk["sample_id"].to_numpy(dtype=str))
    cell_id_chunks.append(chunk["cell_id"].to_numpy(dtype=str))
    total_chunks.append(chunk[gene_cols].to_numpy(dtype=np.float64).sum(axis=1))

sample_id = np.concatenate(sample_id_chunks)
cell_id = np.concatenate(cell_id_chunks)
cell_total = np.concatenate(total_chunks)

n_zero_cells = int((cell_total <= 0).sum())
if n_zero_cells > 0:
    print(
        f"Warning: {n_zero_cells} cells have zero total counts. Assigning them "
        "offset 0 in the per-cell file, matching quasar's internal behaviour. "
        "Donor total exposure is therefore not preserved exactly for the "
        "affected donors.",
        file=sys.stderr,
    )

offsets = pd.DataFrame({"sample_id": sample_id, "cell_id": cell_id})

with np.errstate(divide="ignore", invalid="ignore"):
    offsets["percell"] = np.where(cell_total > 0, np.log(cell_total), 0.0)

    donor_total = pd.Series(cell_total).groupby(sample_id).transform("sum").to_numpy()
    donor_n_cells = pd.Series(cell_total).groupby(sample_id).transform("size").to_numpy()
    donor_mean = donor_total / donor_n_cells
    offsets["donorflat"] = np.where(donor_mean > 0, np.log(donor_mean), 0.0)

grand_mean = cell_total.sum() / len(cell_total)
offsets["constant"] = np.log(grand_mean) if grand_mean > 0 else 0.0

for spec in ["percell", "donorflat", "constant"]:
    offsets[["sample_id", "cell_id", spec]].rename(columns={spec: "offset"}).to_csv(
        f"{cell_type}-offset-{spec}.tsv",
        sep="\t",
        index=False,
        float_format="%.17g",
    )

# Within-donor depth variability is the quantity that would drive an offset effect
# if one existed, so record it alongside the offsets.
depth = pd.DataFrame({"sample_id": sample_id, "total": cell_total})
per_donor = depth.groupby("sample_id")["total"].agg(["size", "mean", "std", "sum"])
per_donor["cv"] = per_donor["std"] / per_donor["mean"]

summary = pd.DataFrame(
    {
        "cell_type": [cell_type],
        "n_cells": [len(cell_total)],
        "n_donors": [per_donor.shape[0]],
        "median_cells_per_donor": [per_donor["size"].median()],
        "median_cell_total": [float(np.median(cell_total))],
        "median_within_donor_depth_cv": [per_donor["cv"].median()],
        "min_within_donor_depth_cv": [per_donor["cv"].min()],
        "max_within_donor_depth_cv": [per_donor["cv"].max()],
        "between_donor_mean_depth_cv": [per_donor["mean"].std() / per_donor["mean"].mean()],
        "n_zero_cells": [n_zero_cells],
    }
)
summary.to_csv(f"{cell_type}-offset-summary.tsv", sep="\t", index=False)
