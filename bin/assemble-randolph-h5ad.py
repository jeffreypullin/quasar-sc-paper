#!/usr/bin/env python3
"""Combine Randolph cluster MTX bundles into a condition-specific H5AD.

Usage:
  assemble-randolph-h5ad.py <condition> <out_h5ad> <cluster_mtx_dir>...
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def read_tsv(path, **kwargs):
    return pd.read_csv(path, sep="\t", keep_default_na=False, **kwargs)


def read_bundle(mtx_dir):
    counts = sp.csr_matrix(sio.mmread(mtx_dir / "matrix.mtx.gz", spmatrix=True)).T.tocsr()
    barcodes = read_tsv(
        mtx_dir / "barcodes.tsv.gz",
        header=None,
        names=["cell_id"],
        dtype=str,
    )
    features = read_tsv(
        mtx_dir / "features.tsv.gz",
        header=None,
        names=["gene_id", "GeneSymbol"],
        dtype=str,
    )
    obs = read_tsv(mtx_dir / "obs.tsv.gz", dtype={"cell_id": str})

    if not np.array_equal(barcodes["cell_id"].to_numpy(), obs["cell_id"].to_numpy()):
        obs = obs.set_index("cell_id").reindex(barcodes["cell_id"]).reset_index()

    return counts, barcodes, features, obs


def union_genes(bundles):
    gene_ids = []
    symbols = {}
    seen = set()
    for _, _, features, _ in bundles:
        for gene_id, symbol in zip(features["gene_id"], features["GeneSymbol"]):
            if gene_id in seen:
                continue
            seen.add(gene_id)
            gene_ids.append(gene_id)
            symbols[gene_id] = symbol
    return gene_ids, symbols


def align_counts(counts, features, gene_ids):
    union_pos = {gene_id: i for i, gene_id in enumerate(gene_ids)}
    local_to_union = np.array(
        [union_pos[gene_id] for gene_id in features["gene_id"]], dtype=np.int32
    )
    coo = counts.tocoo()
    return sp.csr_matrix(
        (coo.data, (coo.row, local_to_union[coo.col])),
        shape=(counts.shape[0], len(gene_ids)),
    )


def main():
    if len(sys.argv) < 4:
        raise SystemExit(
            "Usage: assemble-randolph-h5ad.py "
            "<condition> <out_h5ad> <cluster_mtx_dir>..."
        )

    condition = sys.argv[1]
    out_h5ad = Path(sys.argv[2])
    condition_dirs = [
        Path(cluster_dir) / condition
        for cluster_dir in sys.argv[3:]
        if (Path(cluster_dir) / condition).is_dir()
    ]
    if not condition_dirs:
        raise ValueError(f"No Randolph cluster bundles found for {condition}")
    bundles = [read_bundle(condition_dir) for condition_dir in condition_dirs]
    gene_ids, symbols = union_genes(bundles)

    counts = sp.vstack(
        [align_counts(counts, features, gene_ids) for counts, _, features, _ in bundles],
        format="csr",
    )
    barcodes = pd.concat([bundle[1] for bundle in bundles], ignore_index=True)
    obs = pd.concat([bundle[3] for bundle in bundles], ignore_index=True)
    if barcodes["cell_id"].duplicated().any():
        raise ValueError("Duplicate cell IDs across Randolph cluster bundles")

    adata = ad.AnnData(
        X=counts.astype(np.float32),
        obs=obs.set_index("cell_id"),
        var=pd.DataFrame(
            {"GeneSymbol": [symbols[gene_id] for gene_id in gene_ids]},
            index=gene_ids,
        ),
    )
    adata.var.index = adata.var.index.astype(str)
    adata.var["GeneSymbol"] = adata.var["GeneSymbol"].astype(str)
    adata.obs["individual"] = adata.obs["individual"].astype(str)
    adata.obs["cell_label"] = adata.obs["cell_label"].astype(str)
    adata.obs["sex"] = adata.obs["sex"].astype(int)
    adata.obs["age"] = adata.obs["age"].astype(float)

    if sp.issparse(adata.X):
        adata.X = adata.X.tocsr()

    adata.write_h5ad(out_h5ad, compression="gzip")
    print(
        f"Wrote {out_h5ad}: {adata.n_obs} cells × {adata.n_vars} genes, "
        f"{adata.obs['individual'].nunique()} donors"
    )


if __name__ == "__main__":
    main()
